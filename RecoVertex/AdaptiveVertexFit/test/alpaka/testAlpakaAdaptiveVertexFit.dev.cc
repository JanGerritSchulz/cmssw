#include <iostream>

#include <alpaka/alpaka.hpp>

#include "FWCore/Utilities/interface/stringize.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/warpsize.h"

#include "../../interface/alpaka/VertexFitter.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;



int main() {
  // get the list of devices on the current platform
  auto const& devices = cms::alpakatools::devices<Platform>();

  if (devices.empty()) {
    std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
      "the test will be skipped.\n";
    exit(EXIT_FAILURE);
  }

  for (auto const& device : devices) {
    std::cout << "Test prefix scan on " << alpaka::getName(device) << '\n';
    auto queue = Queue(device);
    const auto warpSize = alpaka::getPreferredWarpSize(device);
  }

  return 0;
}
