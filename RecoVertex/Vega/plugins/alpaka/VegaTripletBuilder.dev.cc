#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

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
                                  VegaPairSoAConstView const pairs,
                                  TrackToPair const* __restrict__ trk1ToPair,
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
          //   auto pair2 = pairs[pairIdx2];

          // FIXME: add triplet cuts

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

    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  TripletParams const& params,
                                  VegaPairSoAView const pairs,
                                  uint32_t const* nTripletCands,
                                  TripletCandSoAConstView tripletCands,
                                  uint32_t* nTriplets,  // Count writes this, Form writes this as well
                                  TripletsArg triplets  // Form writes this
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

          auto triplet = triplets[tripletIdx];
          triplet.idx1() = trkIdx1;
          triplet.idx2() = trkIdx2;
          triplet.idx3() = trkIdx3;
          triplet.x() = pair1.x();
          triplet.y() = pair1.y();
          triplet.z() = pair1.z();
          triplet.chi2() = pair1.chi2();
          triplet.covariance() = pair1.covariance();
        }
      }
    }
  };

  // -------------------------------------------------------------
  // MAIN FUNCTION
  // -------------------------------------------------------------

  VegaTripletSoACollection TripletBuilder::build(TrkSoAConstView trks,
                                                 TrkExtraSoAConstView trksExtra,
                                                 uint32_t const nPairs,
                                                 VegaPairSoAView pairs,
                                                 TrackToPair const* trk1ToPair) const {
    // Find triplets of tracks compatible with vertexing
    auto nTripletCands = cms::alpakatools::make_device_buffer<uint32_t>(queue_);
    alpaka::memset(queue_, nTripletCands, 0u);
    auto nTriplets = cms::alpakatools::make_device_buffer<uint32_t>(queue_);
    alpaka::memset(queue_, nTriplets, 0u);

    const uint32_t nThreadsPerBlock = 64u;
    const uint32_t stride = 4u;
    const auto nThreadsPerStride = nThreadsPerBlock / stride;
    const uint32_t nBlocks = cms::alpakatools::divide_up_by(nPairs + nThreadsPerStride - 1u, nThreadsPerStride);

    using Vec2D = alpaka::Vec<alpaka::DimInt<2u>, uint32_t>;
    Vec2D blocksPerGrid{nBlocks, 1u};  // (row tiles, col tiles) -- Alpaka convention: index 0 = outer/y
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

    const uint32_t blockSize = 64u;
    const uint32_t numberOfBlocks = cms::alpakatools::divide_up_by(tripletCandCapacity + blockSize - 1u, blockSize);
    const auto workDivBuildTriplets = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    // Size the triplet collection: exact count from the Count kernel (requires a
    // device-to-host copy + sync), or the heuristic upper bound.
    uint32_t tripletCapacity = tripletCandCapacity;
    if constexpr (false) {  // currently disabled since expensive, could be ::vega::useExactNTripletsOnHost again
      // const auto workDivFindTripletCands = cms::alpakatools::make_workdiv<Acc2D>(1u, 256u);
      alpaka::exec<Acc1D>(queue_,
                          workDivBuildTriplets,
                          KernelBuildTriplets<KernelMode::Count>{},
                          params_,
                          pairs,
                          nTripletCands.data(),
                          tripletCands.view(),
                          nTriplets.data(),
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
    auto triplets = VegaTripletSoACollection(queue_, tripletCapacity);

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
                        triplets.view());

#if VEGA_TRIPLETS_DEBUG
    auto nPairHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
    alpaka::memcpy(queue_, nPairHost, nPairs);
    alpaka::wait(queue_);  // host sync: triplet count is now valid
    tripletCapacity = *nPairHost.data();

    printf("TripletBuilder::build: formed %d triplets\n", *nPairHost.data());
    alpaka::wait(queue_);
#endif

    // publish the actual triplet count into the SoA (device-to-device copy)
    auto nTripletsView = cms::alpakatools::make_device_view(alpaka::getDev(queue_), triplets.view().nTriplets());
    alpaka::memcpy(queue_, nTripletsView, nTriplets);

#if VEGA_TRIPLETS_DEBUG
    printf("TripletBuilder::build: done\n");
    alpaka::wait(queue_);
#endif

    return triplets;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
