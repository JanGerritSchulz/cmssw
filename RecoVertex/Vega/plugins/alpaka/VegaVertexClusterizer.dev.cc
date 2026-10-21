#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoVertex/Vega/interface/alpaka/PortableOneToManyAssoc.h"
#include "RecoVertex/Vega/plugins/alpaka/VegaGeometry.h"
#include "RecoVertex/Vega/plugins/VegaConstants.h"

#include "VegaVertexClusterizer.h"

#define VEGA_VERTEX_DEBUG 0
#define VEGA_VERTEX_WARNINGS 1
#define VEGA_VERTEX_PRINTALL 0

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  // -------------------------------------------------------------
  // VERTEX CLUSTERING KERNELS
  // -------------------------------------------------------------

  class KernelClusterVertexCandidates {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VertexParams const& params,
                                  TrkSoAConstView trks,
                                  TrkExtraSoAConstView trksExtra,
                                  VegaPairSoAConstView pairs,
                                  VtxSoAView vertices,
                                  VtxToTrack* trkTuplets,
                                  cms::alpakatools::AtomicPairCounter* apc) const {
      for (auto ip : cms::alpakatools::uniform_elements(acc, pairs.nPairs())) {
        TrkIdx trkTuplet[::vega::maxTracksPerVertex];
        auto pair = pairs[ip];
        trkTuplet[0] = pair.idx1();
        trkTuplet[1] = pair.idx2();
        uint32_t nTracksInTuplet = 2;

        // put vertex candidate into the container
        auto iv = trkTuplets->bulkFill(acc, *apc, trkTuplet, nTracksInTuplet);

        if (iv != cms::alpakatools::kOverflow) {
          auto vtx = vertices[iv];
          vtx.x() = 0.0f;
          vtx.y() = 0.0f;
          vtx.z() = 0.0f;
          vtx.t() = 0.0f;
          // add trkToVertex association
        } else {
#if VEGA_VERTEX_WARNINGS
          printf("Warning!!!! Too many vertices (maxNumOfVertices = %d. ip = %d)!\n", vertices.metadata().size(), ip);
#endif
          break;
        }
      }
    }
  };

  // -------------------------------------------------------------
  // VERTEX FILLING KERNELS
  // -------------------------------------------------------------

  class KernelFillTracksIntoVertex {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VtxToTrack const* __restrict__ trkTuplets,
                                  cms::alpakatools::AtomicPairCounter* apc,
                                  VtxSoAView vertices,
                                  VtxTrkSoAView vtxTracks) const {
      // clamp the number of vertices to the capacity of the SoA
      auto nVertices = std::min<int>(apc->get().first, vertices.metadata().size() - 1);

      if (cms::alpakatools::once_per_grid(acc))
        vertices.nVertices() = nVertices;

      const auto maxTrkOffset = uint32_t(vtxTracks.metadata().size());
      for (auto iv : cms::alpakatools::uniform_elements(acc, nVertices)) {
        auto vtx = vertices[iv];
        vtx.nTracks() = static_cast<uint8_t>(trkTuplets->size(iv));
        vtx.tracksOffset() = alpaka::math::min(acc, trkTuplets->off[iv + 1], maxTrkOffset);
      }

      const auto maxTracks = alpaka::math::min(acc, maxTrkOffset, trkTuplets->size());
      for (auto it : cms::alpakatools::uniform_elements(acc, maxTracks)) {
        auto vtxTrks = vtxTracks[it];
        vtxTrks.id() = trkTuplets->content[it];
        vtxTrks.weight() = 1.0f;
      }
    }
  };

  class KernelPrintVertices {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, VtxSoAConstView vertices, VtxTrkSoAConstView vtxTracks) const {
      for (auto iv{0u}; iv < vertices.nVertices(); ++iv) {
        auto vtx = vertices[iv];
        auto beginTrk = (iv == 0u) ? 0u : vertices[iv - 1u].tracksOffset();
        auto endTrk = vtx.tracksOffset();
        printf("KernelPrintVertices: Vertex %d with nTrks=%d, x=%f, y=%f, z=%f, t=%f\n     Tracks (weights): ",
               iv,
               vtx.nTracks(),
               vtx.x(),
               vtx.y(),
               vtx.z(),
               vtx.t());
        for (auto it = beginTrk; it < endTrk; it++)
          printf("%d (%f), ", vtxTracks[it].id(), vtxTracks[it].weight());
        printf("\n");
      }
    }
  };

  // -------------------------------------------------------------
  // MAIN FUNCTION
  // -------------------------------------------------------------

  reco::VertexSoACollection VertexClusterizer::cluster(TrkSoAConstView trks,
                                                       TrkExtraSoAConstView trksExtra,
                                                       VegaPairSoAConstView pairs) const {
    reco::VertexSoACollection vertexCollection(queue_, maxVertices_, 4 * nTracks_);  // FIXME: arbitrary number
    vertexCollection.zeroInitialise(queue_);
    auto vertices = vertexCollection.view().vertex();
    auto vtxTracks = vertexCollection.view().tracks();

#if VEGA_VERTEX_DEBUG
    printf("Create PortableOneToManyAssoc with maxVertices=%d and nVertices*maxTracksPerVertex=%d\n",
           maxVertices_,
           (maxVertices_ * ::vega::maxTracksPerVertex));
#endif

    auto trkTuplets =
        PortableOneToManyAssoc<VtxToTrack>(queue_, maxVertices_, maxVertices_ * ::vega::maxTracksPerVertex);
    auto trkTupletsAPC = cms::alpakatools::make_device_buffer<cms::alpakatools::AtomicPairCounter>(queue_);
    alpaka::memset(queue_, trkTupletsAPC, 0);

    const uint32_t blockSize = 256u;
    const uint32_t numberOfBlocks = cms::alpakatools::divide_up_by(maxVertices_ + blockSize - 1u, blockSize);
    const auto workDivClusterCandidates = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    alpaka::exec<Acc1D>(queue_,
                        workDivClusterCandidates,
                        KernelClusterVertexCandidates{},
                        params_,
                        trks,
                        trksExtra,
                        pairs,
                        vertices,
                        trkTuplets.data(),
                        trkTupletsAPC.data());

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: KernelClusterVertexCandidates done\n");
    alpaka::wait(queue_);
#endif

    const auto workDivFillTracksIntoVertex = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    alpaka::exec<Acc1D>(queue_,
                        workDivFillTracksIntoVertex,
                        typename VtxToTrack::finalizeBulk{},
                        trkTupletsAPC.data(),
                        trkTuplets.data());

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: VtxToTrack::finalizeBulk done\n");
    alpaka::wait(queue_);
#endif

    alpaka::exec<Acc1D>(queue_,
                        workDivFillTracksIntoVertex,
                        KernelFillTracksIntoVertex{},
                        trkTuplets.data(),
                        trkTupletsAPC.data(),
                        vertices,
                        vtxTracks);

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: KernelFillTracksIntoVertex done\n");
    alpaka::wait(queue_);
#endif

#if VEGA_VERTEX_PRINTALL
    const auto workDivPrintVertices = cms::alpakatools::make_workdiv<Acc1D>(1u, 1u);
    alpaka::exec<Acc1D>(queue_, workDivPrintVertices, KernelPrintVertices{}, vertices, vtxTracks);

    printf("VertexClusterizer::cluster: done printing the vertices\n");
    alpaka::wait(queue_);
#endif

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: done\n");
    alpaka::wait(queue_);
#endif

    return vertexCollection;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
