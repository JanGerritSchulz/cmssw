#ifndef RecoVertex_Vega_interface_TrackExtraHost_h
#define RecoVertex_Vega_interface_TrackExtraHost_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "TrackExtraSoA.h"

namespace vega {
  using TrackExtraHost = PortableHostCollection<TrackExtraSoA>;
}
#endif  // RecoVertex_Vega_interface_TrackExtraHost_h