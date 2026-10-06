#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include <alpaka/alpaka.hpp>

// TrackUtilities only included in order to compile SoALayout with Eigen columns
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoVertex/Vega/plugins/alpaka/VegaGeometry.h"

#include "VegaGeometry_test.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  using namespace vega::geometry;
  using TrackPair = std::array<HelixGeom, 2>;
  template <int N>
  using TrackPairArray = std::array<TrackPair, N>;

  namespace vega::geometry::test {

    // Track parameters used for testing
    constexpr int nTestPairs = 3;
    constexpr TrackPairArray<nTestPairs> testTrackPairs{{
        {{{.cx = 0.0f,
           .cy = 0.0f,
           .r = 90.0f,
           .q = 1.0f,
           .refAngle = 3.141592326f,
           .refZ = -0.2f,
           .cotTheta = 0.304520293f},
          {.cx = 150.0f,
           .cy = 0.0f,
           .r = 130.f,
           .q = 1.0f,
           .refAngle = 3.141592326f,
           .refZ = 0.3f,
           .cotTheta = -0.41075232f}}},
        {{{.cx = 0.0f,
           .cy = 0.0f,
           .r = 70.0f,
           .q = 1.0f,
           .refAngle = 3.141592326f,
           .refZ = 0.1f,
           .cotTheta = 0.521095305f},
          {.cx = 250.0f,
           .cy = 0.0f,
           .r = 90.0f,
           .q = 1.0f,
           .refAngle = 3.141592326f,
           .refZ = -0.1f,
           .cotTheta = 0.636653582f}}},
        {{{.cx = 0.0f,
           .cy = 0.0f,
           .r = 40.0f,
           .q = 1.0f,
           .refAngle = 3.141592326f,
           .refZ = 0.0f,
           .cotTheta = -0.2013360025f},
          {.cx = 30.0f,
           .cy = 0.0f,
           .r = 150.0f,
           .q = 1.0f,
           .refAngle = 3.141592326f,
           .refZ = 0.15f,
           .cotTheta = 0.7585837018f}}},
    }};

    constexpr std::array<VertexEstimate, nTestPairs> expectedResults{{
        {.x = 45.666667f, .y = 77.553566f, .z = 11.801519f, .distance3D = 91.268387f},
        {.x = 115.0f, .y = 0.0f, .z = 57.297394f, .distance3D = 145.869324f},
        {.x = -80.0f, .y = 0.0f, .z = 0.074975f, .distance3D = 80.000141f},
    }};

    class KernelVertexEstimator {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrackPair* trackPairs, VertexEstimate* results) const {
        for (auto it : uniform_elements(acc, nTestPairs)) {
          HelixGeom const& t1 = trackPairs[it][0];
          HelixGeom const& t2 = trackPairs[it][1];
          results[it] = estimatedVertex(acc, t1, t2);
        }
      }
    };

    // printer function for fitted vertices
    void printFitResult(std::array<VertexEstimate, nTestPairs> results, size_t i) {
      auto const& r = results[i];
      auto const& t = expectedResults[i];

      std::cout << "Vertex " << i << ":\n";
      std::cout << std::fixed << std::setprecision(6);

      std::cout << "  expected result (x, y, z) = (" << t.x << ", " << t.y << ", " << t.z << ")\n";

      std::cout << "  achieved result (x, y, z) = (" << r.x << ", " << r.y << ", " << r.z << ")\n";

      std::cout << "  expected distance = " << t.distance3D;
      std::cout << "  achieved distance = " << r.distance3D;
      std::cout << "\n";
    }

    // check function for fitted vertices
    void checkFitResult(std::array<VertexEstimate, nTestPairs> results, size_t i) {
      static constexpr float tolerance = 1e-5f;

      std::cout << "Checking vertex " << i << ":";

      auto const& r = results[i];
      auto const& t = expectedResults[i];

      ALPAKA_ASSERT_ACC(std::abs(r.x - t.x) <= (std::abs(r.x) + std::abs(t.x)) * tolerance);
      ALPAKA_ASSERT_ACC(std::abs(r.y - t.y) <= (std::abs(r.y) + std::abs(t.y)) * tolerance);
      ALPAKA_ASSERT_ACC(std::abs(r.z - t.z) <= (std::abs(r.z) + std::abs(t.z)) * tolerance);
      ALPAKA_ASSERT_ACC(std::abs(r.distance3D - t.distance3D) <
                        (std::abs(r.distance3D) + std::abs(t.distance3D)) * tolerance);

      std::cout << " passed the check.\n";
    }

    void runKernels(Queue& queue) {
      // copy Track parameters to device
      auto trackPairs_d = make_device_buffer<TrackPair[]>(queue, nTestPairs);
      auto trackPairs_h = make_host_view(testTrackPairs.data(), nTestPairs);
      alpaka::memcpy(queue, trackPairs_d, trackPairs_h);

      auto blockSize = nTestPairs;
      auto numberOfBlocks = 1;
      auto workDiv1D = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

      auto results_d = make_device_buffer<VertexEstimate[]>(queue, nTestPairs);

      alpaka::exec<Acc1D>(queue, workDiv1D, KernelVertexEstimator{}, trackPairs_d.data(), results_d.data());
      alpaka::wait(queue);

      std::array<VertexEstimate, nTestPairs> results_h{};
      alpaka::memcpy(queue, results_h, results_d);
      alpaka::wait(queue);

      for (size_t i{0}; i < nTestPairs; i++)
        printFitResult(results_h, i);

      alpaka::wait(queue);

      for (size_t i{0}; i < nTestPairs; i++)
        checkFitResult(results_h, i);
    }
  }  // namespace vega::geometry::test
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
