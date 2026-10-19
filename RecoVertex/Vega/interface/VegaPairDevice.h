#ifndef RecoVertex_Vega_interface_VegaPairDevice_h
#define RecoVertex_Vega_interface_VegaPairDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "VegaPairSoA.h"

namespace vega {
  template <typename TDev>
  using VegaPairDevice = PortableDeviceCollection<TDev, VegaPairSoA>;
}
#endif  // RecoVertex_Vega_interface_VegaPairDevice_h