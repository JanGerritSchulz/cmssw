#ifndef RecoVertex_Vega_interface_TrackExtraDevice_h
#define RecoVertex_Vega_interface_TrackExtraDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "TrackExtraSoA.h"

namespace vega {
  template <typename TDev>
  using TrackExtraDevice = PortableDeviceCollection<TDev, TrackExtraSoA>;
}
#endif  // RecoVertex_Vega_interface_TrackExtraDevice_h