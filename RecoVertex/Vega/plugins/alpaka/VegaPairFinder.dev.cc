#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "RecoVertex/Vega/interface/alpaka/IndexPairSoACollection.h"
#include "RecoVertex/Vega/plugins/VegaConstants.h"

#include "VegaPairFinder.h"

#define VEGA_PAIRS_DEBUG 1
#define VEGA_PAIRS_WARNINGS 1

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  // -------------------------------------------------------------
  // TRACK SORTING KERNELS
  // -------------------------------------------------------------

  class KernelSortTracksByDz {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, TrkIdx* sortInd, TrkIdx const nTracks) const {
      for (auto it : cms::alpakatools::uniform_elements(acc, nTracks)) {
        sortInd[it] = it;
      }

      auto trkStateRow0 = trks[0].state();
      // FIXME: index 4 is hardcoded for dz, should be replaced with a proper accessor
      float const* dz = trkStateRow0.data() + 4 * trkStateRow0.stride();

      auto& sortWorkSpace = alpaka::declareSharedVar<TrkIdx[::vega::maxTracksForVertexing], __COUNTER__>(acc);
      // sort using only 24 bits (meaning it might happen that very close tracks in dz are ordered incorrectly,
      // but this is not a problem for the pair finding, also rare-ish: ~10 in 1k tracks)
      cms::alpakatools::radixSort<Acc1D, float, 3>(acc, dz, sortInd, sortWorkSpace, nTracks);
    }
  };

  class KernelPrintSortedTracks {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, TrkIdx* sortInd, TrkIdx const nTracks) const {
      float dzPrev = -9999.f;
      int counterWrongOrder = 0;
      for (int i{0}; i < nTracks; ++i) {
        TrkIdx it = sortInd[i];
        float dz = ::reco::zip(trks, it);
        bool wrongOrder = (dz < dzPrev);
        if (wrongOrder) {
          ++counterWrongOrder;
        }
        printf("KernelPrintSortedTracks: (< than prev)=%d, i=%d, it=%d, dz=%f\n", wrongOrder, i, it, dz);
        dzPrev = dz;
      }
      printf("KernelPrintSortedTracks: Total wrong orders found: %d out of %d\n", counterWrongOrder, nTracks);
    }
  };

  // -------------------------------------------------------------
  // PAIR FINDING KERNELS
  // -------------------------------------------------------------

  template <KernelMode Mode>
  class KernelFindPairsLoose {
  public:
    using CountArg = std::conditional_t<Mode == KernelMode::Form, int const*, int*>;
    using PairsArg = std::conditional_t<Mode == KernelMode::Form, ::reco::IndexPairSoAView, NoInput>;
    static constexpr TrkIdx T = 32u;  // tile size

    ALPAKA_FN_ACC void operator()(Acc2D const& acc,
                                  PairParams const& params,
                                  TrkSoAConstView trks,
                                  TrkIdx const* sortInd,
                                  TrkIdx const nTracks,
                                  int* nPairs,    // Count writes this, Form reads this
                                  PairsArg pairs  // Form writes this
    ) const {
      // if (cms::alpakatools::once_per_grid(acc)) {
      //   *nPairs = 0;
      // }
      // Explore the 2D grid of nTracks x nTracks in tiles of size (tileSize x tileSize)
      // auto blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc);  // Vec2D: (bi, bj)

      // // get block indices
      // TrkIdx bi = blockIdx[0u];  // row tile
      // TrkIdx bj = blockIdx[1u];  // col tile

      // return early if the block is in the lower triangle of the 2D grid
      // (we only need to explore one half of the grid and choose the upper triangle)
      // if (bi > bj)
      //   return;

      // const TrkIdx iMin = bi * T;
      // const TrkIdx jMin = bj * T;

      // // return early if the tile does not contain any tracks
      // if (iMin >= nTracks || jMin >= nTracks)
      //   return;

      // const TrkIdx iMax = alpaka::math::min(acc, static_cast<TrkIdx>(iMin + T), nTracks) - 1u;

      // const float iMaxTileDZ = ::reco::zip(trks, sortInd[iMax]);
      // const float jMinTileDZ = ::reco::zip(trks, sortInd[jMin]);

      // // tiles are sorted by z, so tile boundaries directly bound the tile's range -- no reduction needed
      // if (jMinTileDZ - iMaxTileDZ > params.maxDZ)
      //   return;

      // // define shared memory for the tile's track parameters
      // auto& idz = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      // auto& jdz = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      // auto& iphi = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      // auto& jphi = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      // auto& ieta = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      // auto& jeta = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);

      // --- cooperative tile load: flat loop over T elements, portable ---
      // for (uint32_t t : cms::alpakatools::uniform_elements_x(acc, T)) {
      //   const TrkIdx i = iMin + t;
      //   const TrkIdx j = jMin + t;

      //   // fill shared memory
      //   if (i < nTracks) {
      //     TrkIdx s = sortInd[i];
      //     idz[t] = ::reco::zip(trks, s);
      //     iphi[t] = ::reco::phi(trks, s);
      //     ieta[t] = trks[s].eta();
      //   }
      //   if (j < nTracks) {
      //     TrkIdx s = sortInd[j];
      //     jdz[t] = ::reco::zip(trks, s);
      //     jphi[t] = ::reco::phi(trks, s);
      //     jeta[t] = trks[s].eta();
      //   }
      // }
      // alpaka::syncBlockThreads(acc);

      for (TrkIdx i : cms::alpakatools::uniform_elements_y(acc, nTracks)) {
        for (TrkIdx j : cms::alpakatools::uniform_elements_x(acc, nTracks)) {
          // find pairs by applying compatibility cuts
          const TrkIdx si = sortInd[i], sj = sortInd[j];
          const float idz = ::reco::zip(trks, si);
          const float iphi = ::reco::phi(trks, si);
          const float ieta = trks[si].eta();

          const float jdz = ::reco::zip(trks, sj);
          const float jphi = ::reco::phi(trks, sj);
          const float jeta = trks[sj].eta();

          // only upper triangle
          if (i >= j)
            continue;

          // dz cut
          if (alpaka::math::abs(acc, idz - jdz) > params.maxDZ)
            continue;

          // phi cut
          if (alpaka::math::abs(acc, iphi - jphi) > params.maxDPhi)
            continue;

          // eta cut
          if (alpaka::math::abs(acc, ieta - jeta) > params.maxDEta)
            continue;

          // count or fill the pair
          if constexpr (Mode == KernelMode::Count) {
            alpaka::atomicAdd(acc, nPairs, 1, alpaka::hierarchy::Blocks{});
          } else {
            // fill the pair
            auto idx = alpaka::atomicAdd(acc, nPairs, 1, alpaka::hierarchy::Blocks{});
            if (idx >= pairs.metadata().size()) {
#if VEGA_PAIRS_WARNINGS
              printf("Warning!!!! Too many pairs (maxNumOfPairs = %d)!\n", pairs.metadata().size());
#endif
              alpaka::atomicSub(acc, nPairs, 1, alpaka::hierarchy::Blocks{});
              continue;
            }

            pairs[idx].idx1() = alpaka::math::min(acc, si, sj);
            pairs[idx].idx2() = alpaka::math::max(acc, si, sj);

#if VEGA_PAIRS_DEBUG
            printf("vega::KernelFindPairsLoose: new track pair %d: (%d, %d)\n",
                   idx,
                   alpaka::math::min(acc, si, sj),
                   alpaka::math::max(acc, si, sj));
#endif
          }
        }
      }
    }
  };

  void PairFinder::find(TrkSoAConstView trks) const {
    // Get sorted track indices by dz (z of closest approach)
    auto sortedTrackIndices = cms::alpakatools::make_device_buffer<TrkIdx[]>(queue_, nTracks_);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: created sorted track indices\n");
    alpaka::wait(queue_);
#endif

    const auto workDivSortTracksByDz = cms::alpakatools::make_workdiv<Acc1D>(1u, 256u);
    alpaka::exec<Acc1D>(
        queue_, workDivSortTracksByDz, KernelSortTracksByDz{}, trks, sortedTrackIndices.data(), nTracks_);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: sorted tracks by dz\n");
    const auto workDivPrintSortedTracks = cms::alpakatools::make_workdiv<Acc1D>(1u, 1u);
    alpaka::exec<Acc1D>(
        queue_, workDivPrintSortedTracks, KernelPrintSortedTracks{}, trks, sortedTrackIndices.data(), nTracks_);
#endif

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: done sorting and printing\n");
    alpaka::wait(queue_);
#endif

    // Find pairs of tracks compatible with vertexing
    auto nLoosePairs = cms::alpakatools::make_device_buffer<int>(queue_);
    alpaka::memset(queue_, nLoosePairs, 0);

    constexpr TrkIdx T = KernelFindPairsLoose<KernelMode::Count>::T;  // tile size

    const uint32_t nTiles = cms::alpakatools::divide_up_by(nTracks_ + T - 1u, T);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: making work division\n");
    alpaka::wait(queue_);
#endif

    using Vec2D = alpaka::Vec<alpaka::DimInt<2u>, uint32_t>;
    Vec2D blocksPerGrid{nTiles, nTiles};  // (row tiles, col tiles) -- Alpaka convention: index 0 = outer/y
    Vec2D threadsPerBlock{T, T};

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: still making work division\n");
    alpaka::wait(queue_);
#endif

    const auto workDivFindPairsLoose = cms::alpakatools::make_workdiv<Acc2D>(blocksPerGrid, threadsPerBlock);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: made work division\n");
    alpaka::wait(queue_);
#endif

    // const auto workDivFindPairsLoose = cms::alpakatools::make_workdiv<Acc2D>(1u, 256u);
    alpaka::exec<Acc2D>(queue_,
                        workDivFindPairsLoose,
                        KernelFindPairsLoose<KernelMode::Count>{},
                        params_,
                        trks,
                        sortedTrackIndices.data(),
                        nTracks_,
                        nLoosePairs.data(),
                        NoInput{});  // Count mode, no pairs are filled
#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: ran counting kernel\n");
    alpaka::wait(queue_);
#endif

    // FIXME: decide on better to allocate the exact size after counting or just use a preset size
    auto loosePairs = reco::IndexPairSoACollection(queue_, nTracks_ * 10u);

    // reset nLoosePairs counter for the filling
    alpaka::memset(queue_, nLoosePairs, 0);

    alpaka::exec<Acc2D>(queue_,
                        workDivFindPairsLoose,
                        KernelFindPairsLoose<KernelMode::Form>{},
                        params_,
                        trks,
                        sortedTrackIndices.data(),
                        nTracks_,
                        nLoosePairs.data(),  // nPairs is not used in Form mode
                        loosePairs.view());

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: done\n");
    alpaka::wait(queue_);
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
