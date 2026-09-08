#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "RecoVertex/Vega/plugins/VegaConstants.h"

#include "VegaPairFinder.h"

#define VEGA_PAIRS_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  class KernelSortTracksByDz {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, uint16_t* sortInd, int const nTracks) const {
      for (auto it : cms::alpakatools::uniform_elements(acc, nTracks)) {
        sortInd[it] = it;
      }

      auto trkStateRow0 = trks[0].state();
      // FIXME: index 4 is hardcoded for dz, should be replaced with a proper accessor
      float const* dz = trkStateRow0.data() + 4 * trkStateRow0.stride();

      auto& sortWorkSpace = alpaka::declareSharedVar<uint16_t[::vega::maxTracksForVertexing], __COUNTER__>(acc);
      // sort using only 24 bits
      cms::alpakatools::radixSort<Acc1D, float, 3>(acc, dz, sortInd, sortWorkSpace, nTracks);
    }
  };

  class KernelPrintSortedTracks {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, uint16_t* sortInd, int const nTracks) const {
      float dzPrev = -9999.f;
      int counterWrongOrder = 0;
      for (int i{0}; i < nTracks; ++i) {
        int it = sortInd[i];
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

  void PairFinder::find(TrkSoAConstView trks) const {
    // Get sorted track indices by dz (z of closest approach)
    auto sortedTrackIndices = cms::alpakatools::make_device_buffer<uint16_t[]>(queue, nTracks);

#ifdef VEGA_PAIRS_DEBUG
    printf("PairFinder::find: created sorted track indices\n");
#endif

    const auto workDivSortTracksByDz = cms::alpakatools::make_workdiv<Acc1D>(1u, 256u);
    alpaka::exec<Acc1D>(queue, workDivSortTracksByDz, KernelSortTracksByDz{}, trks, sortedTrackIndices.data(), nTracks);

#ifdef VEGA_PAIRS_DEBUG
    printf("PairFinder::find: sorted tracks by dz\n");
    const auto workDivPrintSortedTracks = cms::alpakatools::make_workdiv<Acc1D>(1u, 1u);
    alpaka::exec<Acc1D>(
        queue, workDivPrintSortedTracks, KernelPrintSortedTracks{}, trks, sortedTrackIndices.data(), nTracks);
#endif

    // Find pairs of tracks compatible with vertexing
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
