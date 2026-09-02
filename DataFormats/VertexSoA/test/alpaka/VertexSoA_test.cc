/**
   Simple test for the reco::VertexSoA data structure
   which inherits from Portable{Host}Collection.

   Creates an instance of the class (automatically allocates
   memory on device), passes the view of the SoA data to
   the kernels which:
   - Fill the SoA with data.
   - Verify that the data written is correct.

   Then, the SoA data are copied back to Host, where
   a temporary host-side view (tmp_view) is created using
   the same Layout to access the data on host and print it.
 */

#include <cstdlib>
#include <unistd.h>

#include <alpaka/alpaka.hpp>

#include "DataFormats/VertexSoA/interface/VertexHost.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "FWCore/Utilities/interface/stringize.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/devices.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "VertexSoA_test.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;
using namespace ALPAKA_ACCELERATOR_NAMESPACE::reco;

// Run 3 values, used for testing
constexpr uint32_t maxTracks = 32 * 1024;
constexpr uint32_t maxVertices = 1024;

int main() {
  // Get the list of devices on the current platform
  auto const& devices = cms::alpakatools::devices<Platform>();
  if (devices.empty()) {
    std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
      "the test will be skipped.\n";
    exit(EXIT_FAILURE);
  }

  // Run the test on each device
  for (const auto& device : devices) {
    Queue queue(device);

    // Inner scope to deallocate memory before destroying the stream
    {
      // Instantiate vertices on device. PortableCollection allocates
      // SoA on device automatically.
      VertexSoACollection vertex_d(queue, maxVertices, maxTracks);
      testVertexSoAT::runKernels(vertex_d.view(), queue);

      // If the device is actually the host, use the collection as-is.
      // Otherwise, copy the data from the device to the host.
#ifdef ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED
      VertexHost vertex_h = std::move(vertex_d);
#else
      VertexHost vertex_h = cms::alpakatools::CopyToHost<VertexSoACollection>::copyAsync(queue, vertex_d);
#endif
      alpaka::wait(queue);
      std::cout << vertex_h.view().vertex().metadata().size() << std::endl;

      // Print results
      std::cout << "id\t"
                << "x\t"
                << "y\t"
                << "z\t"
                << "t\t"
                << "chi2\t"
                << "ndof\t"
                << "nTracks\t"
                << "trackOffsets\t";

      auto vtx_v = vertex_h.view().vertex();
      auto trk_v = vertex_h.view().tracks();
      for (int i = 0; i < 10; ++i) {
        auto vi = vtx_v[i];
        auto ti = trk_v[i];
        std::cout << (int)ti.id() << "\t" << vi.x() << "\t" << vi.y() << "\t" << vi.z() << "\t" << vi.t() << "\t"
                  << vi.chi2() << "\t" << (int)vi.ndof() << "\t" << vi.nTracks() << "\t" << vi.trackOffsets()
                  << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
