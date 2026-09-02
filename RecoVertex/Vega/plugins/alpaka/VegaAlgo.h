#ifndef RecoVertex_Vega_plugins_alpaka_VegaAlgo_h
#define RecoVertex_Vega_plugins_alpaka_VegaAlgo_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoVertex/Vega/interface/TrackExtraSoA.h"

#include "VegaParams.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace cms::alpakatools;
  using VtxSoAView = ::reco::VertexSoAView;
  using VtxTrkSoAView = ::reco::VertexTracksSoAView;
  using TrkSoAConstView = ::reco::TrackSoAConstView;
  using TrkExtraSoAView = ::vega::TrackExtraSoAView;
  using VegaParams = ::vega::VegaParams;

  class VegaAlgo {
  public:
    VegaAlgo(const VegaParams& params) : params_(params) {}

    ~VegaAlgo() = default;

    reco::VertexSoACollection makeAsync(Queue &queue, TrkSoAConstView const &tracksView, int maxVertices) const;

  private:
    const VegaParams params_;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif
