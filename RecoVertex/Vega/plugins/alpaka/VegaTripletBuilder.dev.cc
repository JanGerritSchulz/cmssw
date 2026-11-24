#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoVertex/Vega/interface/alpaka/IndexPairSoACollection.h"
#include "RecoVertex/Vega/plugins/alpaka/VegaGeometry.h"
#include "RecoVertex/Vega/plugins/VegaConstants.h"

#include "VegaTripletBuilder.h"

#define VEGA_TRIPLETS_DEBUG 0
#define VEGA_TRIPLETS_WARNINGS 0
#define VEGA_TRIPLETS_PRINTALL 0

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  // -------------------------------------------------------------
  // TRIPLET BUILDING KERNELS
  // -------------------------------------------------------------
  //
  // Vocabulary used below
  //   pair    : two tracks fitted to a common vertex (index into `pairs`)
  //   triplet : three tracks = two pairs that share exactly one track
  //   "path"  : (pairA, pairB) sharing track `sharedTrk`; the other tracks are `outerTrkA` (in pairA)
  //             and `outerTrkB` (in pairB). A triplet {s,a,b} is reached by 1 path if open, 3 paths if closed.
  //   closed  : the third pair {outerTrkA, outerTrkB} also exists (all three edges of the triangle present)

  // Is `sharedTrk`-centred path (outerTrkA, outerTrkB) the unique representative of its triplet?
  // Requires outerTrkA < outerTrkB (removes the A<->B mirror of every path) and, for closed
  // triplets (three paths), additionally sharedTrk < outerTrkA, i.e. the path centred on the lowest index.
  ALPAKA_FN_ACC inline bool isCanonicalPath(uint32_t sharedTrk, uint32_t outerTrkA, uint32_t outerTrkB, bool closed) {
    if (outerTrkB < outerTrkA)
      return false;
    return !closed || sharedTrk < outerTrkA;
  }

  template <KernelMode Mode>
  class KernelBuildTripletCandidates {
  public:
    using TripletsArg = std::conditional_t<Mode == KernelMode::Form, TripletCandSoAView, NoInput>;

    ALPAKA_FN_ACC void operator()(Acc2D const& acc,
                                  TripletParams const& params,
                                  VegaPairSoAConstView pairs,
                                  TrkToTrk const* __restrict__ trk1ToTrk2,
                                  TrkToPair const* __restrict__ trkToPair,
                                  uint32_t* nTriplets,  // Count writes this, Form writes this as well
                                  TripletsArg triplets  // Form writes this
    ) const {
      // one work item per (pairA, side): `side` selects which of pairA's two tracks is the shared one
      for (auto item : cms::alpakatools::uniform_elements_y(acc, 2u * pairs.nPairs())) {
        auto const pairIdxA = item / 2u;
        bool const sharedIsTrk1 = (item % 2 == 0);
        auto const pairA = pairs[pairIdxA];
        auto const sharedTrk = sharedIsTrk1 ? pairA.trk1() : pairA.trk2();
        auto const outerTrkA = sharedIsTrk1 ? pairA.trk2() : pairA.trk1();

        auto const* __restrict__ sharingPairs = trkToPair->begin(sharedTrk);
        auto const nSharingPairs = trkToPair->size(sharedTrk);
        if (nSharingPairs <= 1)
          continue;

        for (auto j : cms::alpakatools::uniform_elements_x(acc, nSharingPairs)) {
          auto const pairIdxB = sharingPairs[j];
          if (pairIdxB == pairIdxA)
            continue;
 
          auto const pairB = pairs[pairIdxB];
          auto const outerTrkB = (pairB.trk1() != sharedTrk) ? pairB.trk1() : pairB.trk2();
 
          // cheap early exit: this ordering condition is necessary for any canonical path
          if (outerTrkB < outerTrkA)
            continue;
 
          // closed loop <=> third pair {outerTrkA, outerTrkB} exists (look up in the lower-index-keyed map;
          // outerTrkA < outerTrkB here, so outerTrkA is the key)
          bool closed = false;
          for (auto it = trk1ToTrk2->begin(outerTrkA); it != trk1ToTrk2->end(outerTrkA); ++it) {
            if (*it == outerTrkB) {
              closed = true;
              break;
            }
          }
 
          if (!isCanonicalPath(sharedTrk, outerTrkA, outerTrkB, closed))
            continue;
          if (params.requireClosedLoop && !closed)
            continue;

          // cut on the significance of the separation of the two pair vertices:
          //   sep = vtxA - vtxB,  sepSig = |sep| / sigma_along_sep = |sep|^2 / sqrt(sep^T (CA+CB) sep)
          // evaluated squared and without divisions: reject if sepSig > maxSepSig
          // NOTE: for closed triplets this uses only the path centred on the lowest track index, so the
          //       decision still depends on track labelling. Consider max/min over the three paths instead.
          auto const sepX = pairA.x() - pairB.x();
          auto const sepY = pairA.y() - pairB.y();
          auto const sepZ = pairA.z() - pairB.z();
          auto const sep2 = sepX * sepX + sepY * sepY + sepZ * sepZ;
          auto const covSum = pairA.covariance() + pairB.covariance();
          using Cov = ::reco::VertexSoACovIdx;
          auto const sepCovSep = sepX * sepX * covSum(Cov::kVarX) + sepY * sepY * covSum(Cov::kVarY) +
                                 sepZ * sepZ * covSum(Cov::kVarZ) + 2 * sepX * sepY * covSum(Cov::kCovXY) +
                                 2 * sepX * sepZ * covSum(Cov::kCovXZ) + 2 * sepY * sepZ * covSum(Cov::kCovYZ);
 
          if (params.maxSepSig * params.maxSepSig * sepCovSep < sep2 * sep2)
            continue;

          // count or fill the triplet candidate
          if constexpr (Mode == KernelMode::Count) {
            alpaka::atomicAdd(acc, nTriplets, 1u, alpaka::hierarchy::Blocks{});
          } else if constexpr (Mode == KernelMode::Form) {
            auto tripletIdx = alpaka::atomicAdd(acc, nTriplets, 1u, alpaka::hierarchy::Blocks{});

            if (tripletIdx >= static_cast<uint32_t>(triplets.metadata().size())) {
#if VEGA_TRIPLETS_WARNINGS
              printf(
                  "vega::KernelBuildTripletCandidates: Warning!!!! Too many triplet candidates (maxNumOfTriplets = "
                  "%d)!\n",
                  triplets.metadata().size());
#endif
              alpaka::atomicSub(acc, nTriplets, 1u, alpaka::hierarchy::Blocks{});
              continue;
            }

            auto triplet = triplets[tripletIdx];
            triplet.idx1() = pairIdxA;
            triplet.idx2() = pairIdxB;
          }
        }
      }
    }
  };

  template <KernelMode Mode>
  class KernelBuildTriplets {
  public:
    using TripletsArg = std::conditional_t<Mode == KernelMode::Form, VegaTripletSoAView, NoInput>;
    using MapArg = std::conditional_t<Mode == KernelMode::Form, PairToTriplet*, NoInput>;

    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  TripletParams const& params,
                                  VegaPairSoAView pairs,
                                  uint32_t const* nTripletCands,
                                  TripletCandSoAConstView tripletCands,
                                  uint32_t* nTriplets,   // Count writes this, Form writes this as well
                                  TripletsArg triplets,  // Form writes this
                                  MapArg pair1ToTriplet  // Form writes this
    ) const {
      for (auto candIdx : cms::alpakatools::uniform_elements(acc, *nTripletCands)) {
        auto cand = tripletCands[candIdx];
        auto pair1 = pairs[cand.idx1()];
        auto pair2 = pairs[cand.idx2()];
        auto trkIdx1 = pair1.trk1();
        auto trkIdx2 = pair1.trk2();
        auto trkIdx3 = pair2.trk2();

        // FIXME: add triplet cuts

        // count or fill the triplet candidate
        if constexpr (Mode == KernelMode::Count) {
          alpaka::atomicAdd(acc, nTriplets, 1u, alpaka::hierarchy::Blocks{});
        } else if constexpr (Mode == KernelMode::Form) {
          auto tripletIdx = alpaka::atomicAdd(acc, nTriplets, 1u, alpaka::hierarchy::Blocks{});

          if (tripletIdx >= static_cast<uint32_t>(triplets.metadata().size())) {
#if VEGA_TRIPLETS_WARNINGS
            printf(
                "vega::KernelBuildTriplets: Warning!!!! Too many triplets (maxNumOfTriplets = "
                "%d)!\n",
                triplets.metadata().size());
#endif
            alpaka::atomicSub(acc, nTriplets, 1u, alpaka::hierarchy::Blocks{});
            continue;
          }
          pair1.used() = true;
          pair2.used() = true;
          pair2.start() = false;

          auto triplet = triplets[tripletIdx];
          triplet.pair1() = cand.idx1();
          triplet.pair2() = cand.idx2();
          triplet.trk1() = trkIdx1;
          triplet.trk2() = trkIdx2;
          triplet.trk3() = trkIdx3;
          triplet.x() = pair1.x();
          triplet.y() = pair1.y();
          triplet.z() = pair1.z();
          triplet.chi2() = pair1.chi2();
          triplet.covariance() = pair1.covariance();
          triplet.used() = false;

          pair1ToTriplet->count(acc, cand.idx1());
        }
      }
    }
  };

  // -------------------------------------------------------------
  // MAP FILLING KERNELS
  // -------------------------------------------------------------

  class KernelFillPair1ToTriplet {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  uint32_t* nTriplets,
                                  VegaTripletSoAConstView triplets,
                                  PairToTriplet* pair1ToTriplet) const {
      for (auto tripletIdx : cms::alpakatools::uniform_elements(acc, *nTriplets)) {
#if VEGA_TRIPLETS_PRINTALL
        printf("vega::KernelFillPair1ToTriplet: fill %d (pairIdx) -> %d (tripletIdx)\n",
               triplets[tripletIdx].pair1(),
               tripletIdx);
#endif
        pair1ToTriplet->fill(acc, triplets[tripletIdx].pair1(), tripletIdx);
      }
    }
  };

  // -------------------------------------------------------------
  // MAIN FUNCTION
  // -------------------------------------------------------------

  TripletProducts TripletBuilder::build(TrkSoAConstView trks,
                                        TrkExtraSoAConstView trksExtra,
                                        uint32_t const nPairs,
                                        VegaPairSoAView pairs,
                                        TrkToTrk const* trk1ToTrk2,
                                        TrkToPair const* trkToPair) const {
    // Find triplets of tracks compatible with vertexing
    auto nTripletCands = cms::alpakatools::make_device_buffer<uint32_t>(queue_);
    alpaka::memset(queue_, nTripletCands, 0u);
    auto nTriplets = cms::alpakatools::make_device_buffer<uint32_t>(queue_);
    alpaka::memset(queue_, nTriplets, 0u);

    // 2D work division: each block handles a stride-sized set of first-pairs,
    // with the x-threads scanning the candidate second-pairs.
    constexpr uint32_t nThreadsPerBlock = 64u;
    constexpr uint32_t stride = 4u;
    constexpr uint32_t nPairsPerBlock = nThreadsPerBlock / stride;

    using Vec2D = alpaka::Vec<alpaka::DimInt<2u>, uint32_t>;
    Vec2D blocksPerGrid{cms::alpakatools::divide_up_by(nPairs, nPairsPerBlock), 1u};  // index 0 = outer/y
    Vec2D threadsPerBlock{nThreadsPerBlock, stride};
    const auto workDivBuildTripletCands = cms::alpakatools::make_workdiv<Acc2D>(blocksPerGrid, threadsPerBlock);

#if VEGA_TRIPLETS_DEBUG
    printf("TripletBuilder::build: made work division\n");
    alpaka::wait(queue_);
#endif

    // Size the triplet collection: exact count from the Count kernel (requires a
    // device-to-host copy + sync), or the heuristic upper bound.
    uint32_t tripletCandCapacity = nTracks_;
    if constexpr (::vega::useExactNTripletsOnHost) {
      alpaka::exec<Acc2D>(queue_,
                          workDivBuildTripletCands,
                          KernelBuildTripletCandidates<KernelMode::Count>{},
                          params_,
                          pairs,
                          trk1ToTrk2,
                          trkToPair,
                          nTripletCands.data(),
                          NoInput{});  // Count mode, no triplets are filled

      auto nTripletCandsHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
      alpaka::memcpy(queue_, nTripletCandsHost, nTripletCands);
      alpaka::wait(queue_);  // host sync: triplet count is now valid
      tripletCandCapacity = *nTripletCandsHost.data();

#if VEGA_TRIPLETS_DEBUG
      printf("TripletBuilder::build: counted %d triplet candidates\n", *nTripletCandsHost.data());
      alpaka::wait(queue_);
#endif
    }

    auto tripletCands = IndexPairSoACollection(queue_, tripletCandCapacity);

    // reset nTripletCands counter for the filling
    alpaka::memset(queue_, nTripletCands, 0u);

    alpaka::exec<Acc2D>(queue_,
                        workDivBuildTripletCands,
                        KernelBuildTripletCandidates<KernelMode::Form>{},
                        params_,
                        pairs,
                        trk1ToTrk2,
                        trkToPair,
                        nTripletCands.data(),
                        tripletCands.view());

#if VEGA_TRIPLETS_DEBUG
    printf("TripletBuilder::build: formed triplet candidates\n");
    alpaka::wait(queue_);
#endif

    const auto workDivBuildTriplets =
        cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(tripletCandCapacity, 64u), 64u);

    // Size the triplet collection: exact count from the Count kernel (requires a
    // device-to-host copy + sync), or the heuristic upper bound.
    uint32_t tripletCapacity = tripletCandCapacity;
    if constexpr (false) {  // currently disabled since expensive, could be ::vega::useExactNTripletsOnHost again
      alpaka::exec<Acc1D>(queue_,
                          workDivBuildTriplets,
                          KernelBuildTriplets<KernelMode::Count>{},
                          params_,
                          pairs,
                          nTripletCands.data(),
                          tripletCands.view(),
                          nTriplets.data(),
                          NoInput{},
                          NoInput{});  // Count mode, no triplets are filled

      auto nTripletHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
      alpaka::memcpy(queue_, nTripletHost, nTriplets);
      alpaka::wait(queue_);  // host sync: triplet count is now valid
      tripletCapacity = *nTripletHost.data();

#if VEGA_TRIPLETS_DEBUG
      printf("TripletBuilder::build: counted %d triplets\n", *nTripletHost.data());
      alpaka::wait(queue_);
#endif
    }

    // triplets: collection of track triplets
    // pair1ToTriplet: map from track to pairs in which this track fills idx1
    auto triplets = VegaTripletSoACollection(queue_, tripletCapacity);
    auto pair1ToTriplet = PortableOneToManyAssoc<PairToTriplet>(queue_, nPairs, tripletCapacity);

    // reset nTriplets counter for the filling
    alpaka::memset(queue_, nTriplets, 0u);

    alpaka::exec<Acc1D>(queue_,
                        workDivBuildTriplets,
                        KernelBuildTriplets<KernelMode::Form>{},
                        params_,
                        pairs,
                        nTripletCands.data(),
                        tripletCands.view(),
                        nTriplets.data(),
                        triplets.view(),
                        pair1ToTriplet.data());

#if VEGA_TRIPLETS_DEBUG
    auto nTripletHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
    alpaka::memcpy(queue_, nTripletHost, nTriplets);
    alpaka::wait(queue_);  // host sync: triplet count is now valid
    tripletCapacity = *nTripletHost.data();

    printf("TripletBuilder::build: formed %d triplets\n", *nTripletHost.data());
    alpaka::wait(queue_);
#endif

    // publish the actual triplet count into the SoA (device-to-device copy)
    auto nTripletsView = cms::alpakatools::make_device_view(alpaka::getDev(queue_), triplets.view().nTriplets());
    alpaka::memcpy(queue_, nTripletsView, nTriplets);

    // fill the pair1ToTriplet map
    const auto workDivFillPair1ToTriplet =
        cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(tripletCapacity, 128u), 128u);
    PairToTriplet::template launchFinalize<Acc1D>(pair1ToTriplet.view(), queue_);
    alpaka::exec<Acc1D>(queue_,
                        workDivFillPair1ToTriplet,
                        KernelFillPair1ToTriplet{},
                        nTriplets.data(),
                        triplets.view(),
                        pair1ToTriplet.data());

#if VEGA_TRIPLETS_DEBUG
    printf("TripletBuilder::build: done\n");
    alpaka::wait(queue_);
#endif

    return TripletProducts{tripletCapacity, std::move(triplets), std::move(pair1ToTriplet)};
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
