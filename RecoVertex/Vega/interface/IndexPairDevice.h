#ifndef RecoVertex_Vega_interface_IndexPairDevice_h
#define RecoVertex_Vega_interface_IndexPairDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "IndexPairSoA.h"

namespace reco {
  template <typename TDev>
  using IndexPairDevice = PortableDeviceCollection<TDev, IndexPairSoA>;
}
#endif  // RecoVertex_Vega_interface_IndexPairDevice_h