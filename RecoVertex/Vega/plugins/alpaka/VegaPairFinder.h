#ifndef RecoVertex_Vega_plugins_alpaka_VegaPairFinder_h
#define RecoVertex_Vega_plugins_alpaka_VegaPairFinder_h

#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoVertex/Vega/interface/alpaka/PortableOneToManyAssoc.h"
#include "RecoVertex/Vega/interface/alpaka/VegaPairSoACollection.h"

#include "VegaStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace ::vega::structs;

  struct PairProducts {
    uint32_t nPairs;                                 // number of pairs (host side)
    VegaPairSoACollection pairs;                     // the track pairs
    PortableOneToManyAssoc<TrkToPair> trk1ToPair;  // track -> pairs map
  };

  class PairFinder {
  public:
    PairFinder(Params const& params, TrkIdx const nTracks, Queue& queue)
        : params_(params.pair), nTracks_(nTracks), queue_(queue) {}

    ~PairFinder() = default;

    // Non-copyable due to Queue& member
    PairFinder(const PairFinder&) = delete;
    PairFinder& operator=(const PairFinder&) = delete;

    PairProducts find(TrkSoAConstView trks, TrkExtraSoAConstView trksExtra) const;

  private:
    // parameters (stored by value to avoid dangling reference)
    PairParams const params_;
    TrkIdx const nTracks_;

    // alpaka queue
    Queue& queue_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif  // RecoVertex_Vega_plugins_alpaka_VegaPairFinder_h
