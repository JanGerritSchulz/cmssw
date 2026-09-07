#ifndef RecoVertex_Vega_plugins_alpaka_PairFinder_h
#define RecoVertex_Vega_plugins_alpaka_PairFinder_h

#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "VegaStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace ::vega::structs;

  class PairFinder {
  public:
    PairFinder(Params const& params, int const nTracks, Queue& queue)
        : params_d(params.pair), nTracks(nTracks), queue(queue) {}

    ~PairFinder() = default;

    void find(TrkSoAConstView trks) const;

  private:
    // parameters
    PairParams const& params_d;
    int const nTracks;

    // alpaka queue
    Queue& queue;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif  // RecoVertex_Vega_plugins_alpaka_PairFinder_h
