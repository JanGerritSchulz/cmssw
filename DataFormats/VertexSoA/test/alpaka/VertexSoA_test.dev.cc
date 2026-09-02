#include <alpaka/alpaka.hpp>

#include "DataFormats/VertexSoA/interface/VertexDevice.h"
#include "DataFormats/VertexSoA/interface/VertexHost.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT {

  class TestFillKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, ::reco::VertexSoABlocksView view) const {
      ::reco::VertexSoAView vertex_view = view.vertex();
      ::reco::VertexTracksSoAView tracks_view = view.tracks();

      if (cms::alpakatools::once_per_grid(acc)) {
        vertex_view.nVertices() = 420;
      }

      for (int32_t j : cms::alpakatools::uniform_elements(acc, vertex_view.metadata().size())) {
        vertex_view[j].x() = (float)j;
        vertex_view[j].y() = (float)j;
        vertex_view[j].z() = (float)j;
        vertex_view[j].t() = (float)j;
        vertex_view[j].chi2() = (float)j;
        vertex_view[j].ndof() = (float)j;
        vertex_view[j].nTracks() = (uint16_t)j;
      }
      for (int32_t j : cms::alpakatools::uniform_elements(acc, tracks_view.metadata().size())) {
        tracks_view[j].id() = (uint32_t)j;
      }
    }
  };

  class TestVerifyKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, ::reco::VertexSoABlocksView view) const {
      ::reco::VertexSoAView vertex_view = view.vertex();
      ::reco::VertexTracksSoAView tracks_view = view.tracks();


      for (int32_t j : cms::alpakatools::uniform_elements(acc, vertex_view.nVertices())) {
        ALPAKA_ASSERT(vertex_view[j].x() - (float)j < 0.0001);
        ALPAKA_ASSERT(vertex_view[j].y() - (float)j < 0.0001);
        ALPAKA_ASSERT(vertex_view[j].z() - (float)j < 0.0001);
        ALPAKA_ASSERT(vertex_view[j].t() - (float)j < 0.0001);
        ALPAKA_ASSERT(vertex_view[j].chi2() - (float)j < 0.0001);
        ALPAKA_ASSERT(vertex_view[j].ndof() - (float)j < 0.0001);
        ALPAKA_ASSERT(vertex_view[j].nTracks() == uint32_t(j));
      }
      for (int32_t j : cms::alpakatools::uniform_elements(acc, tracks_view.metadata().size())) {
        ALPAKA_ASSERT(tracks_view[j].id() == (uint32_t)j);
      }
    }
  };

  void runKernels(::reco::VertexSoABlocksView view, Queue& queue) {
    uint32_t items = 64;
    uint32_t groups = cms::alpakatools::divide_up_by(view.vertex().metadata().size(), items);
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, TestFillKernel{}, view);
    alpaka::exec<Acc1D>(queue, workDiv, TestVerifyKernel{}, view);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT
