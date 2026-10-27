#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoVertex/Vega/interface/alpaka/PortableOneToManyAssoc.h"
#include "RecoVertex/Vega/plugins/alpaka/VegaGeometry.h"
#include "RecoVertex/Vega/plugins/VegaConstants.h"

#include "VegaVertexClusterizer.h"

#define VEGA_VERTEX_DEBUG 0
#define VEGA_VERTEX_WARNINGS 0
#define VEGA_VERTEX_PRINTALL 0

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  // -------------------------------------------------------------
  // VERTEX CLUSTERING KERNELS
  // -------------------------------------------------------------

  // DEPTH represents the number of available slots for tracks in the given trkTuplet
  template <int DEPTH>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void collectTracks(Acc1D const& acc,
                                                    VegaPairSoAConstView pairs,
                                                    VegaTripletSoAConstView triplets,
                                                    PairToTriplet const* __restrict__ pair1ToTriplet,
                                                    VtxSoAView vertices,
                                                    VtxToTrk* trkTuplets,
                                                    cms::alpakatools::AtomicPairCounter* apc,
                                                    TrkIdx* trkTuplet,
                                                    uint32_t nTracksInTuplet,
                                                    uint32_t lastPairIdx) {
    // the building process for a vertex candidate ends if:
    //  - it has no connected triplet anymore
    //  - the maximum number of tracks per vertex is reached
    // the track tuplet is then saved
    if constexpr (DEPTH < 0) {
      printf("ERROR: KernelClusterVertexCandidates: collectTracks reached maxTracksPerVertex (DEPTH<=0)!\n");
      ALPAKA_ASSERT_ACC(false);
    } else {
      // if there are still slots available in the track tuplet
      // check if there are further connected triplets to be added
      if constexpr (DEPTH > 0) {
        auto const nNextTriplets = pair1ToTriplet->size(lastPairIdx);

        // if further triplets are connected, keep adding tracks
        if (nNextTriplets > 0) {
          auto const* __restrict__ nextTripletsIdxs = pair1ToTriplet->begin(lastPairIdx);

          for (auto i = 0u; i < nNextTriplets; i++) {
            auto nextTripletIdx = nextTripletsIdxs[i];
            auto triplet = triplets[nextTripletIdx];
            trkTuplet[nTracksInTuplet++] = triplet.idx3();
            auto lastPairIdx = triplet.pairIdx2();

            collectTracks<DEPTH - 1>(
                acc, pairs, triplets, pair1ToTriplet, vertices, trkTuplets, apc, trkTuplet, nTracksInTuplet, lastPairIdx);
          }
          // return after adding calling addions of next triplets, saving of the vertices is taken care of by them
          return;
        }
      }

      // no more slots for additional track OR no further triplets connected
      // => save the vertex candidate

      // put vertex candidate into the container
      auto iv = trkTuplets->bulkFill(acc, *apc, trkTuplet, nTracksInTuplet);

      if (iv != cms::alpakatools::kOverflow) {
        auto vtx = vertices[iv];
        // FIXME: these values should be fitted later, we fill them here just to fill something for the time being...
        auto pair = pairs[lastPairIdx];
        vtx.x() = pair.x();
        vtx.y() = pair.y();
        vtx.z() = pair.z();
        vtx.t() = 0.0f;
        vtx.ndof() = 1.0f;
        vtx.chi2() = pair.chi2();
        vtx.covariance() = pair.covariance();
      } else {
#if VEGA_VERTEX_WARNINGS
        printf("Warning!!!! Too many vertices (maxNumOfVertices = %d, iv = %d)!\n", vertices.metadata().size(), iv);
#endif
      }
    }
  }

  class KernelClusterVertexCandidates {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VegaPairSoAConstView pairs,
                                  VegaTripletSoAConstView triplets,
                                  PairToTriplet const* __restrict__ pair1ToTriplet,
                                  VtxSoAView vertices,
                                  VtxToTrk* trkTuplets,
                                  cms::alpakatools::AtomicPairCounter* apc) const {
      for (auto i : cms::alpakatools::uniform_elements(acc, triplets.nTriplets())) {
        TrkIdx trkTuplet[::vega::maxTracksPerVertex];
        auto triplet = triplets[i];

        // if the first pair of the triplet is not considered a valid start, don't start
        auto firstPair = pairs[triplet.pairIdx1()];
        if (!firstPair.start())
          continue;

        trkTuplet[0] = triplet.idx1();
        trkTuplet[1] = triplet.idx2();
        trkTuplet[2] = triplet.idx3();
        uint32_t nTracksInTuplet = 3;

        auto lastPairIdx = triplet.pairIdx2();

        // After adding the first 3 tracks to the vertex,
        // there are maxTracksPerVertex - 3 slots still available
        constexpr int DEPTH = ::vega::maxTracksPerVertex - 3;
        collectTracks<DEPTH>(
            acc, pairs, triplets, pair1ToTriplet, vertices, trkTuplets, apc, trkTuplet, nTracksInTuplet, lastPairIdx);
      }
    }
  };

  class KernelAddVertexCandidatesFromPairs {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  TrkSoAConstView trks,
                                  TrkExtraSoAConstView trksExtra,
                                  VegaPairSoAConstView pairs,
                                  VtxSoAView vertices,
                                  VtxToTrk* trkTuplets,
                                  cms::alpakatools::AtomicPairCounter* apc) const {
      for (auto ip : cms::alpakatools::uniform_elements(acc, pairs.nPairs())) {
        TrkIdx trkTuplet[::vega::maxTracksPerVertex];
        auto pair = pairs[ip];

        // skip pairs that are part of larger objects
        if (pair.used())
          continue;

        trkTuplet[0] = pair.idx1();
        trkTuplet[1] = pair.idx2();
        uint32_t nTracksInTuplet = 2;

        // put vertex candidate into the container
        auto iv = trkTuplets->bulkFill(acc, *apc, trkTuplet, nTracksInTuplet);

        if (iv != cms::alpakatools::kOverflow) {
          auto vtx = vertices[iv];
          vtx.x() = pair.x();
          vtx.y() = pair.y();
          vtx.z() = pair.z();
          vtx.t() = 0.0f;
          vtx.ndof() = 1.0f;
          vtx.chi2() = pair.chi2();
          vtx.covariance() = pair.covariance();
          // add trkToVertex association
        } else {
#if VEGA_VERTEX_WARNINGS
          printf("Warning!!!! Too many vertices (maxNumOfVertices = %d, iv = %d)!\n", vertices.metadata().size(), iv);
#endif
          break;
        }
      }
    }
  };

  class KernelAddVertexCandidatesFromTriplets {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VertexParams const& params,
                                  TrkSoAConstView trks,
                                  TrkExtraSoAConstView trksExtra,
                                  VegaTripletSoAConstView triplets,
                                  VtxSoAView vertices,
                                  VtxToTrk* trkTuplets,
                                  cms::alpakatools::AtomicPairCounter* apc) const {
      for (auto ip : cms::alpakatools::uniform_elements(acc, triplets.nTriplets())) {
        TrkIdx trkTuplet[::vega::maxTracksPerVertex];
        auto triplet = triplets[ip];

        // skip triplets that are part of larger objects
        if (triplet.used())
          continue;

        trkTuplet[0] = triplet.idx1();
        trkTuplet[1] = triplet.idx2();
        trkTuplet[2] = triplet.idx3();
        uint32_t nTracksInTuplet = 3;

        // put vertex candidate into the container
        auto iv = trkTuplets->bulkFill(acc, *apc, trkTuplet, nTracksInTuplet);

        if (iv != cms::alpakatools::kOverflow) {
          auto vtx = vertices[iv];
          vtx.x() = triplet.x();
          vtx.y() = triplet.y();
          vtx.z() = triplet.z();
          vtx.t() = 0.0f;
          vtx.ndof() = 1.0f;
          vtx.chi2() = triplet.chi2();
          vtx.covariance() = triplet.covariance();
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
                                  VtxToTrk const* __restrict__ trkTuplets,
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

  class KernelCountVerticesPerTrack {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VtxSoAConstView vertices,
                                  VtxToTrk const* __restrict__ trkTuplets,
                                  TrkToVtx* trackToVtx) const {
      for (auto idx : cms::alpakatools::uniform_elements(acc, vertices.nVertices())) {
        if (trkTuplets->size(idx) == 0)
          break;  // guard
        for (auto h = trkTuplets->begin(idx); h != trkTuplets->end(idx); ++h)
          trackToVtx->count(acc, *h);
      }
    }
  };

  class KernelFillVerticesPerTrack {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  VtxSoAConstView vertices,
                                  VtxToTrk const* __restrict__ trkTuplets,
                                  TrkToVtx* trackToVtx) const {
      for (auto idx : cms::alpakatools::uniform_elements(acc, vertices.nVertices())) {
        if (trkTuplets->size(idx) == 0)
          break;  // guard
        for (auto h = trkTuplets->begin(idx); h != trkTuplets->end(idx); ++h)
          trackToVtx->fill(acc, *h, idx);
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
                                                       uint32_t const nPairs,
                                                       VegaPairSoAConstView pairs,
                                                       uint32_t const nTriplets,
                                                       VegaTripletSoAConstView triplets,
                                                       PairToTriplet const* pair1ToTriplet) const {
    reco::VertexSoACollection vertexCollection(queue_, maxVertices_, maxVertices_ * avgTracksPerVertex_);
    vertexCollection.zeroInitialise(queue_);
    auto vertices = vertexCollection.view().vertex();
    auto vtxTracks = vertexCollection.view().tracks();

#if VEGA_VERTEX_DEBUG
    printf(
        "VertexClusterizer::cluster: Create PortableOneToManyAssoc with maxVertices=%d and "
        "maxVertices*avgTracksPerVertex=%d\n",
        maxVertices_,
        (maxVertices_ * avgTracksPerVertex_));
#endif

    auto trkTuplets = PortableOneToManyAssoc<VtxToTrk>(queue_, maxVertices_, maxVertices_ * avgTracksPerVertex_);
    auto trkTupletsAPC = cms::alpakatools::make_device_buffer<cms::alpakatools::AtomicPairCounter>(queue_);
    alpaka::memset(queue_, trkTupletsAPC, 0);

    auto trkToVtx = PortableOneToManyAssoc<TrkToVtx>(queue_, nTracks_, maxVertices_ * avgTracksPerVertex_);

    const uint32_t blockSize = 64u;
    const uint32_t numberOfBlocks = cms::alpakatools::divide_up_by(nTriplets + blockSize - 1u, blockSize);
    const auto workDivClusterCandidates = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    alpaka::exec<Acc1D>(queue_,
                        workDivClusterCandidates,
                        KernelClusterVertexCandidates{},
                        pairs,
                        triplets,
                        pair1ToTriplet,
                        vertices,
                        trkTuplets.data(),
                        trkTupletsAPC.data());

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: KernelAddVertexCandidatesFromTriplets done\n");
    alpaka::wait(queue_);
#endif

    const uint32_t blockSizeAddPairs = 64u;
    const uint32_t numberOfBlocksAddPairs =
        cms::alpakatools::divide_up_by(nPairs + blockSizeAddPairs - 1u, blockSizeAddPairs);
    const auto workDivAddPairs = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocksAddPairs, blockSizeAddPairs);

    alpaka::exec<Acc1D>(queue_,
                        workDivAddPairs,
                        KernelAddVertexCandidatesFromPairs{},
                        trks,
                        trksExtra,
                        pairs,
                        vertices,
                        trkTuplets.data(),
                        trkTupletsAPC.data());

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: KernelAddVertexCandidatesFromPairs done\n");
    alpaka::wait(queue_);
#endif

    const auto workDivFillTracksIntoVertex = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    alpaka::exec<Acc1D>(
        queue_, workDivFillTracksIntoVertex, typename VtxToTrk::finalizeBulk{}, trkTupletsAPC.data(), trkTuplets.data());

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: VtxToTrk::finalizeBulk done\n");
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

    const uint32_t blockSizeFill = 64u;
    const uint32_t numberOfBlocksFill =
        cms::alpakatools::divide_up_by(maxVertices_ + blockSizeFill - 1u, blockSizeFill);
    const auto workDivFill = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocksFill, blockSizeFill);

    alpaka::exec<Acc1D>(
        queue_, workDivFill, KernelCountVerticesPerTrack{}, vertices, trkTuplets.data(), trkToVtx.data());
    TrkToVtx::template launchFinalize<Acc1D>(trkToVtx.view(), queue_);
    alpaka::exec<Acc1D>(
        queue_, workDivFill, KernelFillVerticesPerTrack{}, vertices, trkTuplets.data(), trkToVtx.data());

#if VEGA_VERTEX_DEBUG
    printf("VertexClusterizer::cluster: KernelFillVerticesPerTrack done\n");
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
