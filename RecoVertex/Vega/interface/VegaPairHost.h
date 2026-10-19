#ifndef RecoVertex_Vega_interface_VegaPairHost_h
#define RecoVertex_Vega_interface_VegaPairHost_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "VegaPairSoA.h"

namespace vega {
  using VegaPairHost = PortableHostCollection<VegaPairSoA>;
}
#endif  // RecoVertex_Vega_interface_VegaPairHost_h