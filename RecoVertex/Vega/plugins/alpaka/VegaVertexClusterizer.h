#ifndef RecoVertex_Vega_plugins_alpaka_VegaVertexClusterizer_h
#define RecoVertex_Vega_plugins_alpaka_VegaVertexClusterizer_h

#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "VegaStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace ::vega::structs;

  class VertexClusterizer {
  public:
    VertexClusterizer(Params const& params, TrkIdx const nTracks, TrkIdx const maxVertices, Queue& queue)
        : params_(params.vertex), nTracks_(nTracks), maxVertices_(maxVertices), queue_(queue) {}

    ~VertexClusterizer() = default;

    // Non-copyable due to Queue& member
    VertexClusterizer(const VertexClusterizer&) = delete;
    VertexClusterizer& operator=(const VertexClusterizer&) = delete;

    reco::VertexSoACollection cluster(TrkSoAConstView trks,
                                      TrkExtraSoAConstView trksExtra,
                                      VegaPairSoAConstView pairs) const;

  private:
    // parameters (stored by value to avoid dangling reference)
    VertexParams const params_;
    TrkIdx const nTracks_;
    TrkIdx const maxVertices_;

    // alpaka queue
    Queue& queue_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif  // RecoVertex_Vega_plugins_alpaka_VegaVertexClusterizer_h
