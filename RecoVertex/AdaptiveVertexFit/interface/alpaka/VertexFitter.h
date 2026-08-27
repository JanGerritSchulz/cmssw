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
                                                   ::vertexFit::VertexFitResult* result) const {}
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexfit

#endif
