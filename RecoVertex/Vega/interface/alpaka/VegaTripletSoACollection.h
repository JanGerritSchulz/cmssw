#ifndef RecoVertex_Vega_interface_alpaka_VegaTripletSoACollection_h
#define RecoVertex_Vega_interface_alpaka_VegaTripletSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "RecoVertex/Vega/interface/VegaTripletDevice.h"
#include "RecoVertex/Vega/interface/VegaTripletHost.h"
#include "RecoVertex/Vega/interface/VegaTripletSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  using ::vega::VegaTripletDevice;
  using ::vega::VegaTripletHost;
  using VegaTripletSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, VegaTripletHost, VegaTripletDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(vega::VegaTripletSoACollection, ::vega::VegaTripletHost);

#endif  // RecoVertex_Vega_interface_alpaka_VegaTripletSoACollection_h
