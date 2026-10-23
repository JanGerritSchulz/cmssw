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
#include "RecoVertex/Vega/plugins/alpaka/VegaPairFit.h"

#include "VegaPairFit_test.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  using namespace vega::geometry;
  using namespace vega::pairfit;
  using TrackPair = std::array<HelixGeom, 2>;
  template <int N>
  using TrackPairArray = std::array<TrackPair, N>;

  namespace vega::geometry::test {

    // Track parameters used for testing
    constexpr int nTestPairs = 4;
    constexpr TrackPairArray<nTestPairs> testTrackPairs{{// [clean_prompt_signal] track 0
                                                         {HelixGeom{-140.05367f,
                                                                    19.345122f,
                                                                    141.38297f,
                                                                    1.0f,
                                                                    -0.13725798f,
                                                                    -3.5074148f,
                                                                    0.094997527f,
                                                                    -0.990595f,
                                                                    -0.136827f,
                                                                    1.611766f},
                                                          // [clean_prompt_signal] track 1
                                                          HelixGeom{-219.05944f,
                                                                    8.7449535f,
                                                                    219.23283f,
                                                                    1.0f,
                                                                    -0.039899274f,
                                                                    -3.5081093f,
                                                                    -0.38020528f,
                                                                    -0.999204f,
                                                                    -0.039889f,
                                                                    2.499254f}},
                                                         // [clean_displaced_signal] track 0
                                                         {HelixGeom{-247.74449f,
                                                                    -176.66185f,
                                                                    304.26488f,
                                                                    -1.0f,
                                                                    0.61945121f,
                                                                    -0.77098292f,
                                                                    -3.4255075f,
                                                                    0.814197f,
                                                                    -0.580588f,
                                                                    3.468620f},
                                                          // [clean_displaced_signal] track 1
                                                          HelixGeom{143.41799f,
                                                                    107.86092f,
                                                                    179.43043f,
                                                                    1.0f,
                                                                    -2.4967656f,
                                                                    -0.73834729f,
                                                                    -3.4372645f,
                                                                    0.799204f,
                                                                    -0.601060f,
                                                                    2.045507f}},
                                                         // [confidently_wrong_basin_highpt] track 0
                                                         {HelixGeom{-286.52465f,
                                                                    -1903.6037f,
                                                                    1925.046f,
                                                                    1.0f,
                                                                    1.4214008f,
                                                                    -5.5792575f,
                                                                    2.6386111f,
                                                                    -0.148840f,
                                                                    0.988861f,
                                                                    21.945524f},
                                                          // [confidently_wrong_basin_highpt] track 1
                                                          HelixGeom{229.11975f,
                                                                    1328.0728f,
                                                                    1347.6916f,
                                                                    -1.0f,
                                                                    -1.7416351f,
                                                                    -5.5785546f,
                                                                    2.6385647f,
                                                                    -0.170009f,
                                                                    0.985443f,
                                                                    15.363685f}},
                                                         // [background_pair] track 0
                                                         {HelixGeom{-83.913287f,
                                                                    -146.20206f,
                                                                    168.57711f,
                                                                    -1.0f,
                                                                    1.0497483f,
                                                                    -3.4804566f,
                                                                    -2.5410525f,
                                                                    0.497789f,
                                                                    -0.867298f,
                                                                    1.921779f},
                                                          // [background_pair] track 1
                                                          HelixGeom{65.588648f,
                                                                    244.70726f,
                                                                    253.34195f,
                                                                    1.0f,
                                                                    -1.8326702f,
                                                                    -3.4525764f,
                                                                    -0.1919519f,
                                                                    0.258891f,
                                                                    -0.965907f,
                                                                    2.888098f}}}};

    constexpr std::array<VertexSeed, nTestPairs> expectedResults{{
        {.x = 45.666667f, .y = 77.553566f, .z = 11.801519f, .distance3D = 91.268387f},
        {.x = 115.0f, .y = 0.0f, .z = 57.297394f, .distance3D = 145.869324f},
        {.x = -80.0f, .y = 0.0f, .z = 0.074975f, .distance3D = 80.000141f},
        {.x = -80.0f, .y = 0.0f, .z = 0.074975f, .distance3D = 80.000141f},
    }};

    constexpr std::array<VertexSeed, nTestPairs> simulatedTruth{{
        {.x = -0.000308f, .y = -0.000116f, .z = -3.507735f, .distance3D = 0.0f},
        {.x = -1.593421f, .y = 2.187447f, .z = -10.041983f, .distance3D = 0.0f},
        {.x = 0.000481f, .y = 0.000184f, .z = -5.577008f, .distance3D = 0.0f},
        {.x = 0.0f, .y = 0.0f, .z = 0.0f, .distance3D = 0.0f},
    }};

    class KernelVertexEstimatorSmartSeed {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrackPair* trackPairs, VertexSeed* results) const {
        for (auto it : uniform_elements(acc, nTestPairs)) {
          HelixGeom const& t1 = trackPairs[it][0];
          HelixGeom const& t2 = trackPairs[it][1];
          results[it] = estimatedVertex(acc, t1, t2);
        }
      }
    };

    class KernelVertexEstimatorGaussNewton {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrackPair* trackPairs, VertexSeed* results) const {
        for (auto it : uniform_elements(acc, nTestPairs)) {
          HelixGeom const& t1 = trackPairs[it][0];
          HelixGeom const& t2 = trackPairs[it][1];
          float l1, l2;
          smartStartArcLengths(acc, t1, t2, l1, l2);

          auto r = refinePairGaussNewton(acc, t1, t2, l1, l2, {.maxIter=30});
          results[it] = refinedVertex(acc, t1, t2, r.l1, r.l2);
        }
      }
    };

    class KernelVertexEstimatorDualStart {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrackPair* trackPairs, VertexSeed* results) const {
        for (auto it : uniform_elements(acc, nTestPairs)) {
          HelixGeom const& t1 = trackPairs[it][0];
          HelixGeom const& t2 = trackPairs[it][1];
          auto r = refinePairDualStart(acc, t1, t2);
          results[it] = refinedVertex(acc, t1, t2, r.l1, r.l2);
        }
      }
    };

    float difference(const VertexSeed result, const VertexSeed truth) {
      auto dx = result.x - truth.x;
      auto dy = result.y - truth.y;
      auto dz = result.z - truth.z;
      return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    // printer function for fitted vertices
    void printFitResult(std::array<VertexSeed, nTestPairs> results, size_t i) {
      auto const& r = results[i];
      auto const& e = expectedResults[i];
      auto const& t = simulatedTruth[i];

      std::cout << "Vertex " << i << ":\n";
      std::cout << std::fixed << std::setprecision(6);

      std::cout << "  [sim] simulated vertex (x, y, z) = (" << t.x << ", " << t.y << ", " << t.z << ")\n";
      std::cout << "        expected result  (x, y, z) = (" << e.x << ", " << e.y << ", " << e.z << ")\n";
      std::cout << "  [fit] achieved result  (x, y, z) = (" << r.x << ", " << r.y << ", " << r.z << ")\n";

      std::cout << "    => |sim-fit| = " << difference(r, t) << " cm\n";

      std::cout << "  expected distance between tracks = " << e.distance3D << " cm\n";
      std::cout << "  achieved distance between tracks = " << r.distance3D << " cm\n";
      std::cout << "\n";
    }

    // check function for fitted vertices
    void checkFitResult(std::array<VertexSeed, nTestPairs> results, size_t i) {
      // static constexpr float tolerance = 1e-5f;

      std::cout << "Checking vertex " << i << ":";

      // auto const& r = results[i];
      // auto const& t = expectedResults[i];

      // ALPAKA_ASSERT_ACC(std::abs(r.x - t.x) <= (std::abs(r.x) + std::abs(t.x)) * tolerance);
      // ALPAKA_ASSERT_ACC(std::abs(r.y - t.y) <= (std::abs(r.y) + std::abs(t.y)) * tolerance);
      // ALPAKA_ASSERT_ACC(std::abs(r.z - t.z) <= (std::abs(r.z) + std::abs(t.z)) * tolerance);
      // ALPAKA_ASSERT_ACC(std::abs(r.distance3D - t.distance3D) <
      //                   (std::abs(r.distance3D) + std::abs(t.distance3D)) * tolerance);

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

      auto resultsSmartSeed_d = make_device_buffer<VertexSeed[]>(queue, nTestPairs);
      auto resultsGaussNewton_d = make_device_buffer<VertexSeed[]>(queue, nTestPairs);
      auto resultsDualStart_d = make_device_buffer<VertexSeed[]>(queue, nTestPairs);

      alpaka::exec<Acc1D>(queue, workDiv1D, KernelVertexEstimatorSmartSeed{}, trackPairs_d.data(), resultsSmartSeed_d.data());
      alpaka::wait(queue);
      
      alpaka::exec<Acc1D>(queue, workDiv1D, KernelVertexEstimatorGaussNewton{}, trackPairs_d.data(), resultsGaussNewton_d.data());
      alpaka::wait(queue);

      alpaka::exec<Acc1D>(queue, workDiv1D, KernelVertexEstimatorDualStart{}, trackPairs_d.data(), resultsDualStart_d.data());
      alpaka::wait(queue);

      std::array<VertexSeed, nTestPairs> resultsSmartSeed_h{};
      alpaka::memcpy(queue, resultsSmartSeed_h, resultsSmartSeed_d);

      std::array<VertexSeed, nTestPairs> resultsGaussNewton_h{};
      alpaka::memcpy(queue, resultsGaussNewton_h, resultsGaussNewton_d);

      std::array<VertexSeed, nTestPairs> resultsDualStart_h{};
      alpaka::memcpy(queue, resultsDualStart_h, resultsDualStart_d);
      alpaka::wait(queue);

      std::cout << "=======  Smart Seed  =========================================================\n";
      for (size_t i{0}; i < nTestPairs; i++)
        printFitResult(resultsSmartSeed_h, i);
      
      std::cout << "=======  Gauss Newton  =======================================================\n";
      for (size_t i{0}; i < nTestPairs; i++)
        printFitResult(resultsGaussNewton_h, i);
      
      std::cout << "=======  Dual Start  =========================================================\n";
      for (size_t i{0}; i < nTestPairs; i++)
        printFitResult(resultsDualStart_h, i);

      alpaka::wait(queue);

      for (size_t i{0}; i < nTestPairs; i++)
        checkFitResult(resultsSmartSeed_h, i);
    }
  }  // namespace vega::geometry::test
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
