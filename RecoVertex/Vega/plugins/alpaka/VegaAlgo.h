#ifndef RecoVertex_Vega_plugins_alpaka_VegaAlgo_h
#define RecoVertex_Vega_plugins_alpaka_VegaAlgo_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "VegaStructures.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {
  using namespace cms::alpakatools;
  using namespace ::vega::structs;

  class VegaAlgo {
  public:
    VegaAlgo(const Params &params) : params_(params) {}

    ~VegaAlgo() = default;

    reco::VertexSoACollection makeVerticesAsync(Queue &queue,
                                                TrkSoAConstView const &tracksView,
                                                int const maxVertices,
                                                int const nTracksRaw) const;

  private:
    const Params params_;

    // VegaPairs (two vertex-compatible tracks)
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega

#endif
