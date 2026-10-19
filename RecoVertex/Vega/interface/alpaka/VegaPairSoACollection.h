#ifndef RecoVertex_Vega_interface_alpaka_VegaPairSoACollection_h
#define RecoVertex_Vega_interface_alpaka_VegaPairSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "RecoVertex/Vega/interface/VegaPairDevice.h"
#include "RecoVertex/Vega/interface/VegaPairHost.h"
#include "RecoVertex/Vega/interface/VegaPairSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  using ::vega::VegaPairDevice;
  using ::vega::VegaPairHost;
  using VegaPairSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, VegaPairHost, VegaPairDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(vega::VegaPairSoACollection,
                                      ::vega::VegaPairHost);

#endif  // RecoVertex_Vega_interface_alpaka_VegaPairSoACollection_h
