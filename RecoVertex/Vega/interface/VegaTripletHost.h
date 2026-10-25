#ifndef RecoVertex_Vega_interface_VegaTripletHost_h
#define RecoVertex_Vega_interface_VegaTripletHost_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "VegaTripletSoA.h"

namespace vega {
  using VegaTripletHost = PortableHostCollection<VegaTripletSoA>;
}
#endif  // RecoVertex_Vega_interface_VegaTripletHost_h