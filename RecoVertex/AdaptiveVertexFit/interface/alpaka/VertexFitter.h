#ifndef RecoVertex_ApdaptiveVertexFitter_interface_alpaka_VertexFitter
#define RecoVertex_ApdaptiveVertexFitter_interface_alpaka_VertexFitter

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "../VertexFitResult.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexfit {
  using namespace cms::alpakatools;

  // N: number of tracks to consider for the fit
  template <int N>
  class VertexFit {
  public:
    using TrackParameters = Eigen::Matrix<float, 5, N>;
    using TrackCovariances = Eigen::Matrix<float, 15, N>;

    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC ALPAKA_FN_INLINE void operator()(const TAcc& acc,
                                                   const TrackParameters* tracks,
                                                   const TrackCovariances* trackCovs,
                                                   const float bField,
                                                   const size_t nVertices,
                                                   ::vertexFit::VertexFitResult* result) const {
      for (auto it : cms::alpakatools::uniform_elements(acc, nVertices)) {
        result[it].position(0) = 0.0f;
        result[it].position(1) = 0.0f;
        result[it].position(2) = 10.0f;

        result[it].covariances(0, 0) = 1.0f;
        result[it].covariances(0, 1) = 2.0f;
        result[it].covariances(0, 2) = 3.0f;
        result[it].covariances(1, 0) = 4.0f;
        result[it].covariances(1, 1) = 5.0f;
        result[it].covariances(1, 2) = 6.0f;
        result[it].covariances(2, 0) = 7.0f;
        result[it].covariances(2, 1) = 8.0f;
        result[it].covariances(2, 2) = 9.0f;

        result[it].chi2 = 1.0f;
        result[it].ndof = 8;
      }
    }
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexfit

#endif
