#ifndef DataFormats_VertexSoA_interface_VertexDevice_h
#define DataFormats_VertexSoA_interface_VertexDevice_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "DataFormats/VertexSoA/interface/VertexHost.h"
#include "DataFormats/Portable/interface/PortableDeviceCollection.h"

namespace reco {
  template <typename TDev>
  using VertexDevice = PortableDeviceCollection<TDev, reco::VertexSoABlocks>;
}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_VertexDevice_h
