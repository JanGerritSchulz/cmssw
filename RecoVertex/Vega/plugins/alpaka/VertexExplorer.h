#ifndef RecoVertex_Vega_plugins_alpaka_VertexExplorer
#define RecoVertex_Vega_plugins_alpaka_VertexExplorer

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace cms::alpakatools;
  using VtxSoAView = ::reco::VertexSoAView;
  using VtxTrkSoAView = ::reco::VertexTracksSoAView;
  using TrkSoAConstView = ::reco::TrackSoAConstView;

  class VertexExplorer {
  public:
    VertexExplorer() : test(false) {}

    ~VertexExplorer() = default;

    reco::VertexSoACollection makeAsync(Queue &queue, TrkSoAConstView const &tracksView, int maxVertices) const;

  private:
    bool test;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif
