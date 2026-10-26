#ifndef RecoVertex_Vega_plugins_VegaParams_h
#define RecoVertex_Vega_plugins_VegaParams_h

#include <type_traits>
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"

namespace vega {

  struct PairParams {
    const pixelTrack::Quality minTrkQuality;
    float minPt;
    float maxDPhi;
    float maxDEta;
    float maxDZ;
    float maxLinDistance;
    float maxChi2;
    float minDXYSig;
    float maxDXYSig;
  };

  struct TripletParams {
    float maxDistanceSig;
  };

  struct VertexParams {
    float maxChi2;
    float minDXY;
    float maxDXY;
    float minDXYSig;
    float maxDXYSig;
  };

  struct VegaParams {
    PairParams pair;
    TripletParams triplet;
    VertexParams vertex;
  };

  static_assert(std::is_trivially_copyable_v<VegaParams>);

}  // namespace vega

#endif
