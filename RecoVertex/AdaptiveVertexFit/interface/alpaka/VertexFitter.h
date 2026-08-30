#ifndef RecoVertex_ApdaptiveVertexFitter_interface_alpaka_VertexFitter
#define RecoVertex_ApdaptiveVertexFitter_interface_alpaka_VertexFitter

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"
#include "../VertexFitResult.h"
#include "../VertexFitUtils.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;

  // position at arc length
  template <alpaka::concepts::Acc TAcc, typename TrackParameters, typename TrackGeomParameters, typename VertexPosition>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void positionAtArcLength(const TAcc& acc,
                                                          const TrackParameters& trks,
                                                          const TrackGeomParameters& geos,
                                                          const uint8_t it,
                                                          const float l,
                                                          VertexPosition& position) {
    const float angleL = geos.refAngle[it] - geos.q[it] * l / geos.r[it];
    position(0) = geos.cx[it] + geos.r[it] * alpaka::math::cos(acc, angleL);
    position(1) = geos.cy[it] + geos.r[it] * alpaka::math::sin(acc, angleL);
    position(2) = trks.dz[it] + trks.cotTheta[it] * l;
  }

  // N: number of tracks to consider for the fit
  template <int N>
  class VertexFitter {
  public:
    using TrackParameters = ::vertexfit::TrackParams<N>;
    using TrackCovariances = ::vertexfit::TrackCov<N>;
    using TrackGeomParameters = ::vertexfit::TrackGeomParams<N>;
    using VertexPosition = Eigen::Vector3d;

    static constexpr float SPEED_OF_LIGHT_FACTOR = 0.003f;

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC ALPAKA_FN_INLINE void operator()(const TAcc& acc,
                                                   const uint8_t* nTracks,
                                                   const TrackParameters* tracks,
                                                   const TrackCovariances* trackCovs,
                                                   const VertexPosition* seedPos,
                                                   const float bField,
                                                   const size_t nVertices,
                                                   ::vertexfit::VertexFitResult* result) const {
      for (auto iv : cms::alpakatools::uniform_elements(acc, nVertices)) {
        const auto trks = tracks[iv];
        // const auto covs = trackCovs[iv];
        TrackGeomParameters geos;
        // get derrived geometric track parameters like radius R and center C=(cx,cy) in the x-y plane
        for (uint8_t it{0}; it < nTracks[iv]; it++) {
          geos.q[it] = trks.qOverPt[it] > 0.f ? 1.0f : -1.0f;
          geos.r[it] = 1 / (alpaka::math::abs(acc, trks.qOverPt[it]) * SPEED_OF_LIGHT_FACTOR * bField);
          const float sinPhi = alpaka::math::sin(acc, trks.phi[it]);
          const float cosPhi = alpaka::math::cos(acc, trks.phi[it]);
          geos.cx[it] = (trks.dxy[it] + geos.q[it] * geos.r[it]) * sinPhi;
          geos.cy[it] = (trks.dxy[it] - geos.q[it] * geos.r[it]) * cosPhi;
          geos.refAngle[it] = reducePhiRange(acc, trks.phi[it] + geos.q[it] * std::numbers::pi_v<float> / 2.0f);
        }

        //
        positionAtArcLength(acc, trks, geos, 0u, 2.0f, result[iv].position);
        // result[iv].position(0) = 0.0f;
        // result[iv].position(1) = 0.0f;
        // result[iv].position(2) = 10.0f;

        result[iv].covariances(0, 0) = 1.0f;
        result[iv].covariances(0, 1) = 2.0f;
        result[iv].covariances(0, 2) = 3.0f;
        result[iv].covariances(1, 0) = 4.0f;
        result[iv].covariances(1, 1) = 5.0f;
        result[iv].covariances(1, 2) = 6.0f;
        result[iv].covariances(2, 0) = 7.0f;
        result[iv].covariances(2, 1) = 8.0f;
        result[iv].covariances(2, 2) = 9.0f;

        result[iv].chi2 = 1.0f;
        result[iv].ndof = 8;
      }
    }
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
