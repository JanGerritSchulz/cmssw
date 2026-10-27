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

  template <KernelMode Mode>
  class KernelBuildTripletCandidates {
  public:
    using TripletsArg = std::conditional_t<Mode == KernelMode::Form, TripletCandSoAView, NoInput>;

    ALPAKA_FN_ACC void operator()(Acc2D const& acc,
                                  TripletParams const& params,
                                  VegaPairSoAConstView pairs,
                                  TrkToPair const* __restrict__ trk1ToPair,
                                  uint32_t* nTriplets,  // Count writes this, Form writes this as well
                                  TripletsArg triplets  // Form writes this
    ) const {
      for (auto pairIdx1 : cms::alpakatools::uniform_elements_y(acc, pairs.nPairs())) {
        auto pair1 = pairs[pairIdx1];
        auto trk2Idx = pair1.idx2();  // second trk in first pair == first trk in second pair
        auto const* __restrict__ potentialPairs2 = trk1ToPair->begin(trk2Idx);
        auto const nPotentialPairs2 = trk1ToPair->size(trk2Idx);

        if (nPotentialPairs2 == 0)
          continue;

        for (auto j : cms::alpakatools::uniform_elements_x(acc, nPotentialPairs2)) {
          auto pairIdx2 = potentialPairs2[j];
          auto pair2 = pairs[pairIdx2];

          // cut on the distance significance between two pairs
          auto dx = pair1.x() - pair2.x();
          auto dy = pair1.y() - pair2.y();
          auto dz = pair1.z() - pair2.z();
          auto dist2 = dx * dx + dy * dy + dz * dz;
          auto covTot = pair1.covariance() + pair2.covariance();
          using C = ::reco::VertexSoACovIdx;
          auto distErr2TimesDist2 = dx * dx * covTot(C::kVarX) + dy * dy * covTot(C::kVarY) +
                                    dz * dz * covTot(C::kVarZ) + 2 * dx * dy * covTot(C::kCovXY) +
                                    2 * dx * dz * covTot(C::kCovXZ) + 2 * dy * dz * covTot(C::kCovYZ);

          if (params.maxDistanceSig * params.maxDistanceSig * distErr2TimesDist2 < dist2 * dist2)
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
            triplet.idx1() = pairIdx1;
            triplet.idx2() = pairIdx2;
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
        auto trkIdx1 = pair1.idx1();
        auto trkIdx2 = pair1.idx2();
        auto trkIdx3 = pair2.idx2();

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
          triplet.pairIdx1() = cand.idx1();
          triplet.pairIdx2() = cand.idx2();
          triplet.idx1() = trkIdx1;
          triplet.idx2() = trkIdx2;
          triplet.idx3() = trkIdx3;
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
               triplets[tripletIdx].pairIdx1(),
               tripletIdx);
#endif
        pair1ToTriplet->fill(acc, triplets[tripletIdx].pairIdx1(), tripletIdx);
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
                                        TrkToPair const* trk1ToPair) const {
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
                          trk1ToPair,
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
                        trk1ToPair,
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
