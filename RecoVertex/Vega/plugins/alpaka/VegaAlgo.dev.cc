#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"
#include "RecoVertex/Vega/interface/alpaka/TrackExtraSoACollection.h"

#include "VegaAlgo.h"

#undef PIXVERTEX_DEBUG_PRODUCE
namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace vega {

    class Kernel_buildTrackExtra {
    public:
      static constexpr float SPEED_OF_LIGHT_FACTOR = 0.003f;

      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    TrkSoAConstView trks,
                                    TrkExtraSoAView trksExtra,
                                    const float bField) const {
        for (auto it : cms::alpakatools::uniform_elements(acc, trks.nTracks())) {
          const float phi = ::reco::phi(trks, it);
          const float dxy = ::reco::tip(trks, it);
          const float q = ::reco::charge(trks, it);
          const float r = 1 / (alpaka::math::abs(acc, q) * SPEED_OF_LIGHT_FACTOR * bField);
          trksExtra[it].q() = q;
          trksExtra[it].r() = r;
          const float sinPhi = alpaka::math::sin(acc, phi);
          const float cosPhi = alpaka::math::cos(acc, phi);
          trksExtra[it].cx() = (dxy + q * r) * sinPhi;
          trksExtra[it].cy() = (dxy - q * r) * cosPhi;
          trksExtra[it].refAngle() = reducePhiRange(acc, phi + q * std::numbers::pi_v<float> / 2.0f);
        }
      }
    };

    reco::VertexSoACollection VegaAlgo::makeAsync(Queue& queue,
                                                  TrkSoAConstView const& trks,
                                                  int maxVertices) const {
      const auto maxTracks = trks.metadata().size();
      const float bField = 3.8f;  // FIXME: get from EventSetup

      auto trksExtra = TrackExtraSoACollection(queue, maxTracks);

      // Compute and fill TracksExtras
      const uint32_t blockSize = 128;
      const uint32_t numberOfBlocks =
          cms::alpakatools::divide_up_by(trks.metadata().size() + blockSize - 1, blockSize);
      const auto buildTrackExtraWorkDiv = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);
      alpaka::exec<Acc1D>(
          queue, buildTrackExtraWorkDiv, Kernel_buildTrackExtra{}, trks, trksExtra.view(), bField);

      reco::VertexSoACollection vertexCollection(queue, maxVertices, maxTracks);
      vertexCollection.zeroInitialise(queue);
      auto vertices = vertexCollection.view().vertex();
      auto vtxTracks = vertexCollection.view().tracks();

      return vertexCollection;
    }
  }  // namespace vega
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
