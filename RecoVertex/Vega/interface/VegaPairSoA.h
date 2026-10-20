#ifndef RecoVertex_Vega_interface_VegaPairsSoA_h
#define RecoVertex_Vega_interface_VegaPairsSoA_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace vega {

  // pair of tracks (two vertex-compatible tracks)
  GENERATE_SOA_LAYOUT(VegaPairSoALayout,
                      SOA_COLUMN(uint32_t, idx1),    // index of track 1
                      SOA_COLUMN(uint32_t, idx2),    // index of track 2
                      SOA_COLUMN(float, arcLen1),    // arc length of track 1 at closest approach
                      SOA_COLUMN(float, arcLen2),    // arc length of track 2 at closest approach
                      SOA_SCALAR(uint32_t, nPairs))  // number of pairs

  using VegaPairSoA = VegaPairSoALayout<>;
  using VegaPairSoAView = VegaPairSoA::View;
  using VegaPairSoAConstView = VegaPairSoA::ConstView;

}  // namespace vega

#endif  // RecoVertex_Vega_interface_VegaPairsSoA_h
