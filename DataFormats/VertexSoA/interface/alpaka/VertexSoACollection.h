#ifndef DataFormats_VertexSoA_interface_VertexSoACollection_h
#define DataFormats_VertexSoA_interface_VertexSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/VertexSoA/interface/VertexDevice.h"
#include "DataFormats/VertexSoA/interface/VertexHost.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::reco {

  using ::reco::VertexDevice;
  using ::reco::VertexHost;
  using VertexSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, VertexHost, VertexDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::reco

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(reco::VertexSoACollection, reco::VertexHost);

#endif  // DataFormats_VertexSoA_interface_VertexSoACollection_h
