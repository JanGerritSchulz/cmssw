#ifndef RecoVertex_Vega_interface_TrackExtraSoA_h
#define RecoVertex_Vega_interface_TrackExtraSoA_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace vega {

  GENERATE_SOA_LAYOUT(TrackExtraSoALayout,
                      SOA_COLUMN(float, cx),        // circle center in x-y plane (x-coordinate)
                      SOA_COLUMN(float, cy),        // circle center in x-y plane (y-coordinate)
                      SOA_COLUMN(float, r),         // circle radius
                      SOA_COLUMN(float, q),         // charge/sign
                      SOA_COLUMN(float, refAngle))  // reference angle as seen from circle center

  using TrackExtraSoA = TrackExtraSoALayout<>;
  using TrackExtraSoAView = TrackExtraSoA::View;
  using TrackExtraSoAConstView = TrackExtraSoA::ConstView;

}  // namespace vega

#endif  // RecoVertex_Vega_interface_TrackExtraSoA_h
