#ifndef RecoVertex_Vega_interface_VegaParams_h
#define RecoVertex_Vega_interface_VegaParams_h

#include <type_traits>

namespace vega {

  struct PairParams {
    float minPt;
    float maxDPhi;
    float maxDEta;
    float maxDZ;
    float maxLinDistance;
    float max3DDistance;
  };

  struct TripletParams {
    float max3DDistance;
  };

  struct VegaParams {
    PairParams pair;
    TripletParams triplet;
  };

  static_assert(std::is_trivially_copyable_v<VegaParams>);

}  // namespace vega

#endif
