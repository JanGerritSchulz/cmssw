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
      auto blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc);     // Vec2D: (bi, bj)
      auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc);  // Vec2D: (ti, tj)

      // get block indices
      TrkIdx bi = blockIdx[0u];  // row tile
      TrkIdx bj = blockIdx[1u];  // col tile

      // return early if the block is in the lower triangle of the 2D grid
      // (we only need to explore one half of the grid and choose the upper triangle)
      if (bi > bj)
        return;

      TrkIdx iMin = bi * T;
      TrkIdx iMax = alpaka::math::min(acc, static_cast<TrkIdx>(iMin + T), nTracks) - 1u;
      TrkIdx jMin = bj * T;

      // return early if the tile does not contain any tracks
      if (iMin >= nTracks || jMin >= nTracks)
        return;

      float iMaxTileDZ = ::reco::zip(trks, sortInd[iMax]);
      float jMinTileDZ = ::reco::zip(trks, sortInd[jMin]);

      // tiles are sorted by z, so tile boundaries directly bound the tile's range -- no reduction needed
      if (jMinTileDZ - iMaxTileDZ > params.maxDZ)
        return;

      // get thread indices
      TrkIdx ti = threadIdx[0u];  // row thread
      TrkIdx tj = threadIdx[1u];  // col thread

      // get the index of the tracks of the thread
      TrkIdx i = iMin + ti;
      TrkIdx j = jMin + tj;

      // define shared memory for the tile's track parameters
      auto& idz = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      auto& jdz = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      auto& iphi = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      auto& jphi = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      auto& ieta = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);
      auto& jeta = alpaka::declareSharedVar<float[T], __COUNTER__>(acc);

      // fill shared memory
      if (i < nTracks && tj == 0u) {  // load row tile once per row
        TrkIdx srcIdx = sortInd[i];
        idz[ti] = ::reco::zip(trks, srcIdx);
        iphi[ti] = ::reco::phi(trks, srcIdx);
        ieta[ti] = trks[srcIdx].eta();
      }
      if (j < nTracks && ti == 0u) {  // load col tile once per column
        TrkIdx srcIdx = sortInd[j];
        jdz[tj] = ::reco::zip(trks, srcIdx);
        jphi[tj] = ::reco::phi(trks, srcIdx);
        jeta[tj] = trks[srcIdx].eta();
      }
      alpaka::syncBlockThreads(acc);

      // find pairs by applying compatibility cuts
      // both tracks valid
      if (i >= nTracks || j >= nTracks)
        return;

      // only upper triangle
      if (bi == bj && i >= j)
        return;

      // dz cut
      if (alpaka::math::abs(acc, idz[ti] - jdz[tj]) > params.maxDZ)
        return;

      // phi cut
      if (alpaka::math::abs(acc, iphi[ti] - jphi[tj]) > params.maxDPhi)
        return;

      // eta cut
      if (alpaka::math::abs(acc, ieta[ti] - jeta[tj]) > params.maxDEta)
        return;

      // count or fill the pair
      if constexpr (Mode == KernelMode::Count) {
        alpaka::atomicAdd(acc, nPairs, 1, alpaka::hierarchy::Blocks{});
      } else if constexpr (Mode == KernelMode::Form) {
        // fill the pair
        auto ind = alpaka::atomicAdd(acc, nPairs, 1, alpaka::hierarchy::Blocks{});
        if (ind >= pairs.metadata().size()) {
#if VEGA_PAIRS_WARNINGS
          printf("Warning!!!! Too many pairs (maxNumOfPairs = %d)!\n", pairs.metadata().size());
#endif
          alpaka::atomicSub(acc, nPairs, 1, alpaka::hierarchy::Blocks{});
          return;
        }

        TrkIdx srcIdxI = sortInd[i];
        TrkIdx srcIdxJ = sortInd[j];
        TrkIdx idx1 = alpaka::math::min(acc, srcIdxI, srcIdxJ);
        TrkIdx idx2 = alpaka::math::max(acc, srcIdxI, srcIdxJ);
        pairs[ind].idx1() = idx1;
        pairs[ind].idx2() = idx2;

#if VEGA_PAIRS_DEBUG
        printf("vega::KernelFindPairsLoose: new track pair (%d, %d)\n", idx1, idx2);
#endif
      }
    }
  };

  void PairFinder::find(TrkSoAConstView trks) const {
    // Get sorted track indices by dz (z of closest approach)
    auto sortedTrackIndices = cms::alpakatools::make_device_buffer<TrkIdx[]>(queue, nTracks);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: created sorted track indices\n");
    alpaka::wait(queue);
#endif

    const auto workDivSortTracksByDz = cms::alpakatools::make_workdiv<Acc1D>(1u, 256u);
    alpaka::exec<Acc1D>(queue, workDivSortTracksByDz, KernelSortTracksByDz{}, trks, sortedTrackIndices.data(), nTracks);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: sorted tracks by dz\n");
    const auto workDivPrintSortedTracks = cms::alpakatools::make_workdiv<Acc1D>(1u, 1u);
    alpaka::exec<Acc1D>(
        queue, workDivPrintSortedTracks, KernelPrintSortedTracks{}, trks, sortedTrackIndices.data(), nTracks);
#endif

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: done sorting and printing\n");
    alpaka::wait(queue);
#endif

    // Find pairs of tracks compatible with vertexing
    auto nLoosePairs = cms::alpakatools::make_device_buffer<int>(queue);
    alpaka::memset(queue, nLoosePairs, 0);

    constexpr TrkIdx T = KernelFindPairsLoose<KernelMode::Count>::T;  // tile size

    const uint32_t nTiles = cms::alpakatools::divide_up_by(nTracks + T - 1u, T);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: making work division\n");
    alpaka::wait(queue);
#endif

    using Vec2D = alpaka::Vec<alpaka::DimInt<2u>, uint32_t>;
    Vec2D blocksPerGrid{nTiles, nTiles};  // (row tiles, col tiles) -- Alpaka convention: index 0 = outer/y
    Vec2D threadsPerBlock{T, T};

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: still making work division\n");
    alpaka::wait(queue);
#endif

    const auto workDivFindPairsLoose = cms::alpakatools::make_workdiv<Acc2D>(blocksPerGrid, threadsPerBlock);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: made work division\n");
    alpaka::wait(queue);
#endif

    // const auto workDivFindPairsLoose = cms::alpakatools::make_workdiv<Acc2D>(1u, 256u);
    alpaka::exec<Acc2D>(queue,
                        workDivFindPairsLoose,
                        KernelFindPairsLoose<KernelMode::Count>{},
                        params_d,
                        trks,
                        sortedTrackIndices.data(),
                        nTracks,
                        nLoosePairs.data(),
                        NoInput{});  // Count mode, no pairs are filled
#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: ran counting kernel\n");
    alpaka::wait(queue);
#endif

    // FIXME: decide on better to allocate the exact size after counting or just use a preset size
    auto loosePairs = reco::IndexPairSoACollection(queue, nTracks * 10u);

    // reset nLoosePairs counter for the filling
    alpaka::memset(queue, nLoosePairs, 0);

    alpaka::exec<Acc2D>(queue,
                        workDivFindPairsLoose,
                        KernelFindPairsLoose<KernelMode::Form>{},
                        params_d,
                        trks,
                        sortedTrackIndices.data(),
                        nTracks,
                        nLoosePairs.data(),  // nPairs is not used in Form mode
                        loosePairs.view());

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: done\n");
    alpaka::wait(queue);
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
