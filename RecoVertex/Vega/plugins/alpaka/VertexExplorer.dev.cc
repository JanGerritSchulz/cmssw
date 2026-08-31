#include <alpaka/alpaka.hpp>

#include "VertexExplorer.h"

#undef PIXVERTEX_DEBUG_PRODUCE
namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace vega {

    reco::VertexSoACollection VertexExplorer::makeAsync(Queue& queue,
                                                        TrkSoAConstView const& tracksView,
                                                        int maxVertices) const {
      const auto maxTracks = tracksView.metadata().size();
      reco::VertexSoACollection vertexCollection(queue, maxVertices, maxTracks);
      vertexCollection.zeroInitialise(queue);
      auto vertices = vertexCollection.view().vertex();
      auto vtxTracks = vertexCollection.view().tracks();

      return vertexCollection;
    }
  }  // namespace vega
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
