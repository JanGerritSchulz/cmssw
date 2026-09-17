#ifndef RecoVertex_Vega_interface_IndexPairHost_h
#define RecoVertex_Vega_interface_IndexPairHost_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "IndexPairSoA.h"

namespace reco {
  using IndexPairHost = PortableHostCollection<IndexPairSoA>;
}
#endif  // RecoVertex_Vega_interface_IndexPairHost_h