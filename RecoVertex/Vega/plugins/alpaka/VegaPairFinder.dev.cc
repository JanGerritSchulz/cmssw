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
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, TrkIdx* sortIdx, TrkIdx const nTracks) const {
      for (auto it : cms::alpakatools::uniform_elements(acc, nTracks)) {
        sortIdx[it] = it;
      }

      auto trkStateRow0 = trks[0].state();
      // FIXME: index 4 is hardcoded for dz, should be replaced with a proper accessor
      float const* dz = trkStateRow0.data() + 4 * trkStateRow0.stride();

      auto& sortWorkSpace = alpaka::declareSharedVar<TrkIdx[::vega::maxTracksForVertexing], __COUNTER__>(acc);
      // sort using only 24 bits (meaning it might happen that very close tracks in dz are ordered incorrectly,
      // but this is not a problem for the pair finding, also rare-ish: ~10 in 1k tracks)
      cms::alpakatools::radixSort<Acc1D, float, 3>(acc, dz, sortIdx, sortWorkSpace, nTracks);
    }
  };

  class KernelPrintSortedTracks {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, TrkIdx* sortIdx, TrkIdx const nTracks) const {
      float dzPrev = -9999.f;
      int counterWrongOrder = 0;
      for (int i{0}; i < nTracks; ++i) {
        TrkIdx it = sortIdx[i];
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
                                  TrkIdx const* sortIdx,
                                  TrkIdx const nTracks,
                                  int* nPairs,    // Count writes this, Form writes this as well
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

      // const float iMaxTileDZ = ::reco::zip(trks, sortIdx[iMax]);
      // const float jMinTileDZ = ::reco::zip(trks, sortIdx[jMin]);

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
      //     TrkIdx s = sortIdx[i];
      //     idz[t] = ::reco::zip(trks, s);
      //     iphi[t] = ::reco::phi(trks, s);
      //     ieta[t] = trks[s].eta();
      //   }
      //   if (j < nTracks) {
      //     TrkIdx s = sortIdx[j];
      //     jdz[t] = ::reco::zip(trks, s);
      //     jphi[t] = ::reco::phi(trks, s);
      //     jeta[t] = trks[s].eta();
      //   }
      // }
      // alpaka::syncBlockThreads(acc);

      for (TrkIdx i : cms::alpakatools::uniform_elements_y(acc, nTracks)) {
        const TrkIdx si = sortIdx[i];
        const float idxy = ::reco::tip(trks, si);
        const float iz = ::reco::zip(trks, si);
        const float iphi = ::reco::phi(trks, si);
        const float ieta = trks[si].eta();
        const float icottheta = trks[si].state()(3);
        const float icosphi = alpaka::math::cos(acc, iphi);
        const float isinphi = alpaka::math::sin(acc, iphi);
        const float ix = -idxy * isinphi;
        const float iy = idxy * icosphi;

        for (TrkIdx j : cms::alpakatools::uniform_elements_x(acc, nTracks)) {
          // only upper triangle
          if (i >= j)
            continue;

          // find pairs by applying compatibility cuts
          const TrkIdx sj = sortIdx[j];

          // dz cut
          const float jz = ::reco::zip(trks, sj);
          if (alpaka::math::abs(acc, iz - jz) > params.maxDZ)
            continue;

          // phi cut
          const float jphi = ::reco::phi(trks, sj);
          if (alpaka::math::abs(acc, iphi - jphi) > params.maxDPhi)
            continue;

          // eta cut
          const float jeta = trks[sj].eta();
          if (alpaka::math::abs(acc, ieta - jeta) > params.maxDEta)
            continue;

          // linear distance cut
          const float jdxy = ::reco::tip(trks, sj);
          const float jcottheta = trks[sj].state()(3);
          const float jcosphi = alpaka::math::cos(acc, jphi);
          const float jsinphi = alpaka::math::sin(acc, jphi);
          const float jx = -jdxy * jsinphi;
          const float jy = jdxy * jcosphi;

          const float nx = isinphi * jcottheta - jsinphi * icottheta;
          const float ny = icottheta * jcosphi - jcottheta * icosphi;
          const float nz = icosphi * jsinphi - jcosphi * isinphi;
          const float dx = ix - jx;
          const float dy = iy - jy;
          const float dz = iz - jz;

          const float maxLinDistance2 = (params.maxLinDistance * params.maxLinDistance) * (nx * nx + ny * ny + nz * nz);
          const float linearApproxDistance = dx * nx + dy * ny + dz * nz;
          if (linearApproxDistance * linearApproxDistance > maxLinDistance2)
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

    // Size the pair collection: exact count from the Count kernel (requires a
    // device-to-host copy + sync), or the heuristic upper bound.
    uint32_t pairCapacity = nTracks_ * 10u;
    if constexpr (::vega::useExactNPairsOnHost) {
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

      auto nLoosePairsHost = cms::alpakatools::make_host_buffer<int, Platform>();
      alpaka::memcpy(queue_, nLoosePairsHost, nLoosePairs);
      alpaka::wait(queue_);  // host sync: pair count is now valid
      pairCapacity = static_cast<uint32_t>(*nLoosePairsHost.data());

#if VEGA_PAIRS_DEBUG
      printf("PairFinder::find: counted %d loose pairs\n", *nLoosePairsHost.data());
      alpaka::wait(queue_);
#endif
    }

    auto loosePairs = reco::IndexPairSoACollection(queue_, pairCapacity);

    // reset nLoosePairs counter for the filling
    alpaka::memset(queue_, nLoosePairs, 0);

    alpaka::exec<Acc2D>(queue_,
                        workDivFindPairsLoose,
                        KernelFindPairsLoose<KernelMode::Form>{},
                        params_,
                        trks,
                        sortedTrackIndices.data(),
                        nTracks_,
                        nLoosePairs.data(),
                        loosePairs.view());

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: done\n");
    alpaka::wait(queue_);
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
