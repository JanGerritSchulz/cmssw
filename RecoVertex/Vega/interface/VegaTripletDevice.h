#ifndef RecoVertex_Vega_interface_VegaTripletDevice_h
#define RecoVertex_Vega_interface_VegaTripletDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "VegaTripletSoA.h"

namespace vega {
  template <typename TDev>
  using VegaTripletDevice = PortableDeviceCollection<TDev, VegaTripletSoA>;
}
#endif  // RecoVertex_Vega_interface_VegaTripletDevice_h