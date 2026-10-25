#ifndef RecoVertex_Vega_plugins_alpaka_VegaTripletBuilder_h
#define RecoVertex_Vega_plugins_alpaka_VegaTripletBuilder_h

#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoVertex/Vega/interface/alpaka/VegaTripletSoACollection.h"

#include "VegaStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace ::vega::structs;

  class TripletBuilder {
  public:
    TripletBuilder(Params const& params, TrkIdx const nTracks, Queue& queue)
        : params_(params.triplet), nTracks_(nTracks), queue_(queue) {}

    ~TripletBuilder() = default;

    // Non-copyable due to Queue& member
    TripletBuilder(const TripletBuilder&) = delete;
    TripletBuilder& operator=(const TripletBuilder&) = delete;

    VegaTripletSoACollection build(TrkSoAConstView trks,
                                   TrkExtraSoAConstView trksExtra,
                                   uint32_t const nPairs,
                                   VegaPairSoAView pairs,
                                   TrackToPair const* trk1ToPair) const;

  private:
    // parameters (stored by value to avoid dangling reference)
    TripletParams const params_;
    TrkIdx const nTracks_;

    // alpaka queue
    Queue& queue_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif  // RecoVertex_Vega_plugins_alpaka_VegaTripletBuilder_h
