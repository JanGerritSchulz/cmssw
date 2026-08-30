#include <iostream>

#include <alpaka/alpaka.hpp>

#include "FWCore/Utilities/interface/stringize.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/warpsize.h"

#include "../../interface/alpaka/VertexFitter.h"

constexpr int maxTracksPerVertex = 8;

using namespace ALPAKA_ACCELERATOR_NAMESPACE;
using namespace cms::alpakatools;

using TrackParameters = VertexFitter<maxTracksPerVertex>::TrackParameters;
using TrackCovariances = VertexFitter<maxTracksPerVertex>::TrackCovariances;
using VertexPosition = VertexFitter<maxTracksPerVertex>::VertexPosition;

// Build a TrackParameters matrix (5 params x maxTracksPerVertex tracks) from
// human-readable "one row per track" data: {param0, param1, param2, param3, param4}.
template <typename Matrix, size_t nParam>
Matrix makeTrackParameters(std::array<std::array<float, nParam>, maxTracksPerVertex> const& tracksByRow) {
  Matrix m;
  for (std::size_t track = 0; track < maxTracksPerVertex; ++track) {
    m.phi[track] = tracksByRow[track][0];
    m.dxy[track] = tracksByRow[track][1];
    m.qOverPt[track] = tracksByRow[track][2];
    m.cotTheta[track] = tracksByRow[track][3];
    m.dz[track] = tracksByRow[track][4];
  }
  return m;
}

// bField = 3.8 T (geometry.DEFAULT_B_FIELD) - MUST match the C++ harness's own field constant
// N_TRACKS per test vertex = 8 (shared by every case, see module docstring)
// test case order:
//   [0] pure_wellConstrained_nearPV - pure, N=8, ptRange=(1.0, 3.0)
//   [1] pure_displacedSV_moderatePt - pure, N=8, ptRange=(1.0, 5.0)
//   [2] pure_displacedSV_heterogeneousPt_hard - pure, N=8, ptRange=(0.7, 25.0)
//   [3] mixed_classicMergedCandidate - MIXED: 5 tracks from VA=(-0.72, 0.227, -1.87) (reported truth) + 3 tracks from VB=(0.00108, -0.000183, -1.96) (not reported - only one truth slot)
//   [4] mixed_closeMerge_hardToSplit - MIXED: 4 tracks from VA=(-0.15, 0.1, -0.5) (reported truth) + 4 tracks from VB=(0.001, -0.001, -0.52) (not reported - only one truth slot)

// [pure_wellConstrained_nearPV] Python reference fit (fitVertexWeighted): V=(0.051695, -0.024150, 0.018764) cm, sigma=(0.004297, 0.005029, 0.009611) cm, chi2/ndof=0.7731
//   true vertex: (0.050000, -0.030000, 0.020000) cm
// [pure_displacedSV_moderatePt] Python reference fit (fitVertexWeighted): V=(0.796653, -0.500661, 1.204701) cm, sigma=(0.002987, 0.002836, 0.007942) cm, chi2/ndof=0.9181
//   true vertex: (0.800000, -0.500000, 1.200000) cm
// [pure_displacedSV_heterogeneousPt_hard] Python reference fit (fitVertexWeighted): V=(-0.598063, 0.898260, -2.005064) cm, sigma=(0.002272, 0.004117, 0.006565) cm, chi2/ndof=0.6065
//   true vertex: (-0.600000, 0.900000, -2.000000) cm
// [mixed_classicMergedCandidate] Python reference fit (fitVertexWeighted): V=(-0.324951, 0.048284, -2.048569) cm, sigma=(0.003613, 0.004241, 0.009118) cm, chi2/ndof=516.8090
//   true vertex: (-0.720000, 0.227000, -1.870000) cm
// [mixed_closeMerge_hardToSplit] Python reference fit (fitVertexWeighted): V=(-0.101394, 0.061800, -0.474222) cm, sigma=(0.004991, 0.003496, 0.008980) cm, chi2/ndof=21.7069
//   true vertex: (-0.150000, 0.100000, -0.500000) cm

constexpr float bField = 3.8f;
constexpr int nTestVertices = 1;
std::array<uint8_t, nTestVertices> nTracksPerVertex{8};

// testVertexTrackParameters: 5 test vertices x 40 floats each
std::array<std::array<std::array<float, 5>, maxTracksPerVertex>, nTestVertices> testTrackParametersByTrack = {{
    // V0
    {{{-0.35424806f, 0.014392954f, -0.54955008f, 1.9500064f, -0.11085489f},
      {2.266534f, 0.048376235f, 0.64028004f, -0.93060119f, -0.018712188f},
      {-1.3554805f, -0.054512728f, 0.43756635f, 1.0558566f, -0.034839579f},
      {-2.5618537f, -0.03251268f, 0.63244026f, 0.56689364f, 0.042290629f},
      {1.8452501f, 0.028539189f, -0.53646257f, -0.090866506f, 0.024347949f},
      {-0.48359425f, 0.010888745f, 0.55378452f, -1.9181496f, 0.1206391f},
      {-0.82108677f, -0.018609972f, 0.89783605f, 1.4773234f, -0.048106194f},
      {2.2920082f, 0.026734507f, -0.59739806f, 0.21975429f, 0.040328298f}}}
    // V1, V2, ... add more vertices here as {{ ...tracks... }}, if nTestVertices > 1
}};

// Transpose each vertex's rows into the Eigen matrix layout the kernel actually needs.
std::array<TrackParameters, nTestVertices> testTrackParameters = [] {
  std::array<TrackParameters, nTestVertices> result;
  for (std::size_t v = 0; v < nTestVertices; ++v) {
    result[v] = makeTrackParameters<TrackParameters>(testTrackParametersByTrack[v]);
  }
  return result;
}();

// testTrackCovariances: 5 test vertices x 120 floats each
std::array<TrackCovariances, nTestVertices> testTrackCovariances = {{
    // V0
    TrackCovariances{{{{{1.406688e-05f,
                         4.1178298e-05f,
                         0.00010569133f,
                         0.0f,
                         0.0f,
                         0.00012753452f,
                         0.00023455388f,
                         0.0f,
                         0.0f,
                         0.0020219088f,
                         0.0f,
                         0.0f,
                         9.5846589e-05f,
                         -0.00031277455f,
                         0.0010562802f}},
                       {{1.227701e-05f,
                         3.5938772e-05f,
                         9.2243167e-05f,
                         0.0f,
                         0.0f,
                         0.00011130703f,
                         0.00020470925f,
                         0.0f,
                         0.0f,
                         0.0017646412f,
                         0.0f,
                         0.0f,
                         8.3651071e-05f,
                         -0.00027297713f,
                         0.00092187911f}},
                       {{5.6197844e-06f,
                         1.6450923e-05f,
                         4.2224182e-05f,
                         0.0f,
                         0.0f,
                         5.0950641e-05f,
                         9.3705376e-05f,
                         0.0f,
                         0.0f,
                         0.00080776205f,
                         0.0f,
                         0.0f,
                         3.8291162e-05f,
                         -0.0001249549f,
                         0.00042198888f}},
                       {{8.9078948e-06f,
                         2.6076284e-05f,
                         6.692936e-05f,
                         0.0f,
                         0.0f,
                         8.076163e-05f,
                         0.00014853197f,
                         0.0f,
                         0.0f,
                         0.00128038f,
                         0.0f,
                         0.0f,
                         6.0695147e-05f,
                         -0.00019806545f,
                         0.00066889267f}},
                       {{1.3869153e-05f,
                         4.0599488e-05f,
                         0.00010420571f,
                         0.0f,
                         0.0f,
                         0.00012574188f,
                         0.00023125695f,
                         0.0f,
                         0.0f,
                         0.0019934885f,
                         0.0f,
                         0.0f,
                         9.4499352e-05f,
                         -0.00030837814f,
                         0.0010414329f}},
                       {{1.1558748e-05f,
                         3.3836186e-05f,
                         8.6846511e-05f,
                         0.0f,
                         0.0f,
                         0.00010479505f,
                         0.0001927328f,
                         0.0f,
                         0.0f,
                         0.0016614014f,
                         0.0f,
                         0.0f,
                         7.8757093e-05f,
                         -0.00025700669f,
                         0.00086794488f}},
                       {{6.7977157e-06f,
                         1.9899108e-05f,
                         5.1074555e-05f,
                         0.0f,
                         0.0f,
                         6.1630117e-05f,
                         0.00011334643f,
                         0.0f,
                         0.0f,
                         0.00097707249f,
                         0.0f,
                         0.0f,
                         4.6317156e-05f,
                         -0.00015114599f,
                         0.00051043959f}},
                       {{7.8329157e-06f,
                         2.2929473e-05f,
                         5.8852517e-05f,
                         0.0f,
                         0.0f,
                         7.1015549e-05f,
                         0.00013060756f,
                         0.0f,
                         0.0f,
                         0.0011258674f,
                         0.0f,
                         0.0f,
                         5.3370632e-05f,
                         -0.00017416348f,
                         0.00058817263f}}}}}
    // add more vertices here as {{ ...8 rows... }}, if nTestVertices > 1
}};

// testTrueVertex: 5 test vertices x 3 floats each
std::array<VertexPosition, nTestVertices> testVertexSeeds{{
    VertexPosition{{{0.0f, 0.0f, 0.0f}}}  // V0
}};

// testTrueVertex: 5 test vertices x 3 floats each
std::array<std::array<float, 3>, nTestVertices> testTrueVertices{{
    {0.05f, -0.03f, 0.02f}  // V0
}};

// printer function for fitted vertices
void printFitResult(std::array<vertexfit::VertexFitResult, nTestVertices> results, size_t i) {
  auto const& r = results[i];
  auto const& t = testTrueVertices[i];

  std::cout << "Vertex " << i << ":\n";
  std::cout << std::fixed << std::setprecision(6);

  std::cout << "  true position (x, y, z) = (" << t[0] << ", " << t[1] << ", " << t[2] << ")\n";

  std::cout << "  reco position (x, y, z) = (" << r.position.x() << ", " << r.position.y() << ", " << r.position.z()
            << ")\n";

  std::cout << "  covariance matrix:\n";
  for (int row = 0; row < 3; ++row) {
    std::cout << "    ";
    for (int col = 0; col < 3; ++col) {
      std::cout << std::setw(14) << r.covariances(row, col) << " ";
    }
    std::cout << "\n";
  }

  std::cout << "  chi2 / ndof = " << r.chi2 << " / " << r.ndof;
  if (r.ndof > 0) {
    std::cout << "  (chi2/ndof = " << (r.chi2 / static_cast<float>(r.ndof)) << ")";
  }
  std::cout << "\n";
}

int main() {
  // get the list of devices on the current platform
  auto const& devices = cms::alpakatools::devices<Platform>();

  if (devices.empty()) {
    std::cerr << "No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
      "the test will be skipped.\n";
    exit(EXIT_FAILURE);
  }

  //   // fill host containers of input vertex candidates
  //   std::array<std::array<TrackCovariances, maxTracksPerVertex>, nTestVertices> trackCov_h;
  //   for (size_t v{0}; v < nTestVertices; v++) {
  //     for (size_t t{0}; t < maxTracksPerVertex; t++) {
  //       auto& civ = testTrackCovariances[v][t];
  //       auto& cov = trackCov_h[v][t];
  //       cov(0, 0) = civ[0];
  //       cov(1, 0) = cov(0, 1) = civ[1];
  //       cov(2, 0) = cov(0, 2) = civ[2];
  //       cov(3, 0) = cov(0, 3) = civ[3];
  //       cov(4, 0) = cov(0, 4) = civ[4];
  //       cov(1, 1) = civ[5];
  //       cov(2, 1) = cov(1, 2) = civ[6];
  //       cov(3, 1) = cov(1, 3) = civ[7];
  //       cov(4, 1) = cov(1, 4) = civ[8];
  //       cov(2, 2) = civ[9];
  //       cov(3, 2) = cov(2, 3) = civ[10];
  //       cov(4, 2) = cov(2, 4) = civ[11];
  //       cov(3, 3) = civ[12];
  //       cov(4, 3) = cov(3, 4) = civ[13];
  //       cov(4, 4) = civ[14];
  //     }
  //   }

  for (auto const& device : devices) {
    std::cout << "Test prefix scan on " << alpaka::getName(device) << '\n';
    auto queue = Queue(device);
    // const auto warpSize = alpaka::getPreferredWarpSize(device);

    // copy Track numbers to device
    auto nTracks_d = make_device_buffer<uint8_t[]>(queue, nTestVertices);
    auto nTracks_h = make_host_view(nTracksPerVertex.data(), nTestVertices);
    alpaka::memcpy(queue, nTracks_d, nTracks_h);

    // copy Track parameters to device
    auto trackParams_d = make_device_buffer<TrackParameters[]>(queue, nTestVertices);
    auto trackParams_h = make_host_view(testTrackParameters.data(), nTestVertices);
    alpaka::memcpy(queue, trackParams_d, trackParams_h);

    // copy Track covariances to device
    auto trackCovs_d = make_device_buffer<TrackCovariances[]>(queue, nTestVertices);
    auto trackCovs_h = make_host_view(testTrackCovariances.data(), nTestVertices);
    alpaka::memcpy(queue, trackCovs_d, trackCovs_h);

    // copy Vertex seeds to device
    auto vertexSeed_d = make_device_buffer<VertexPosition[]>(queue, nTestVertices);
    auto vertexSeed_h = make_host_view(testVertexSeeds.data(), nTestVertices);
    alpaka::memcpy(queue, vertexSeed_d, vertexSeed_h);

    auto blockSize = 64;
    auto numberOfBlocks = 1;  //cms::alpakatools::divide_up_by(nTestVertices, blockSize);
    auto workDiv1D = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    auto results_d = make_device_buffer<vertexfit::VertexFitResult[]>(queue, nTestVertices);

    alpaka::exec<Acc1D>(queue,
                        workDiv1D,
                        VertexFitter<maxTracksPerVertex>{},
                        nTracks_d.data(),
                        trackParams_d.data(),
                        trackCovs_d.data(),
                        vertexSeed_d.data(),
                        bField,
                        nTestVertices,
                        results_d.data());

    std::array<vertexfit::VertexFitResult, nTestVertices> results_h{};
    alpaka::memcpy(queue, results_h, results_d);
    // alpaka::wait(queue);

    ALPAKA_ASSERT_ACC(8 == results_h[0].ndof);

    for (size_t i{0}; i < nTestVertices; i++)
      printFitResult(results_h, i);
  }

  return 0;
}
