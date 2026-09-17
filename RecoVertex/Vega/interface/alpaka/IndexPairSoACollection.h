#ifndef RecoVertex_Vega_interface_alpaka_IndexPairSoACollection_h
#define RecoVertex_Vega_interface_alpaka_IndexPairSoACollection_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "RecoVertex/Vega/interface/IndexPairDevice.h"
#include "RecoVertex/Vega/interface/IndexPairHost.h"
#include "RecoVertex/Vega/interface/IndexPairSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::reco {

  using ::reco::IndexPairDevice;
  using ::reco::IndexPairHost;
  using IndexPairSoACollection =
      std::conditional_t<std::is_same_v<Device, alpaka::DevCpu>, IndexPairHost, IndexPairDevice<Device>>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::reco

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(reco::IndexPairSoACollection,
                                      ::reco::IndexPairHost);

#endif  // RecoVertex_Vega_interface_alpaka_IndexPairSoACollection_h
