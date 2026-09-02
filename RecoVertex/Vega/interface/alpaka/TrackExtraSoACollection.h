#ifndef RecoVertex_Vega_interface_alpaka_TrackExtraSoACollection_h
#define RecoVertex_Vega_interface_alpaka_TrackExtraSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "RecoVertex/Vega/interface/TrackExtraDevice.h"
#include "RecoVertex/Vega/interface/TrackExtraHost.h"
#include "RecoVertex/Vega/interface/TrackExtraSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  using ::vega::TrackExtraDevice;
  using ::vega::TrackExtraHost;
  using TrackExtraSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, TrackExtraHost, TrackExtraDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(vega::TrackExtraSoACollection,
                                      ::vega::TrackExtraHost);

#endif  // RecoVertex_Vega_interface_alpaka_TrackExtraSoACollection_h
