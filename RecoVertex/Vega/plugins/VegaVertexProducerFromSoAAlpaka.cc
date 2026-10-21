#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackSoA/interface/alpaka/TrackUtilities.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexSoA/interface/VertexHost.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/TTHelpers.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#define VEGA_CONVERTER_DEBUG 0
#define VEGA_CONVERTER_WARNINGS 1

class VegaVertexProducerFromSoAAlpaka : public edm::global::EDProducer<> {
public:
  using IndToEdm = std::vector<uint32_t>;
  struct FitParamsAVF {
    double sigmacut;
    double tini;
    double ratio;
  };

  explicit VegaVertexProducerFromSoAAlpaka(const edm::ParameterSet &iConfig);
  ~VegaVertexProducerFromSoAAlpaka() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  FitParamsAVF getFitParams(const edm::ParameterSet &iConfig) const;
  void produce(edm::StreamID iStreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const override;

  edm::EDGetTokenT<reco::VertexHost> vertexSoAToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<IndToEdm> indToEdmToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> trackBuilderToken_;

  const float minTrackWeight_;  // minimum track weight for a track to be put in the vertex

  const bool refitVertices_;
  const FitParamsAVF fitParams_;
};

VegaVertexProducerFromSoAAlpaka::VegaVertexProducerFromSoAAlpaka(const edm::ParameterSet &iConfig)
    : vertexSoAToken_(consumes(iConfig.getParameter<edm::InputTag>("vertexSoACollection"))),
      tracksToken_(consumes(iConfig.getParameter<edm::InputTag>("trackCollection"))),
      indToEdmToken_(consumes(iConfig.getParameter<edm::InputTag>("trackCollection"))),
      trackBuilderToken_(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      minTrackWeight_(iConfig.getParameter<float>("minTrackWeight")),
      refitVertices_(iConfig.getParameter<bool>("refitVerticesWithAVF")),
      fitParams_(getFitParams(iConfig)) {
  produces<reco::VertexCollection>();
}

VegaVertexProducerFromSoAAlpaka::FitParamsAVF VegaVertexProducerFromSoAAlpaka::getFitParams(
    const edm::ParameterSet &iConfig) const {
  // configuration for track pair finding
  auto fitConfig = iConfig.getParameter<edm::ParameterSet>("paramsAVF");

  return FitParamsAVF{.sigmacut = fitConfig.getParameter<double>("sigmacut"),
                      .tini = fitConfig.getParameter<double>("tini"),
                      .ratio = fitConfig.getParameter<double>("ratio")};
}

void VegaVertexProducerFromSoAAlpaka::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("trackCollection", edm::InputTag("pixelTracks"));
  desc.add<edm::InputTag>("vertexSoACollection", edm::InputTag("vegaVerticesAlpaka"));

  desc.add<float>("minTrackWeight", 0.5f);

  // refitting options
  desc.add<bool>("refitVerticesWithAVF", false);
  edm::ParameterSetDescription fitParams;
  fitParams.add<double>("sigmacut", 3.0);
  fitParams.add<double>("tini", 256.0);
  fitParams.add<double>("ratio", 0.25);
  desc.add<edm::ParameterSetDescription>("paramsAVF", fitParams);

  descriptions.addWithDefaultLabel(desc);
}

void VegaVertexProducerFromSoAAlpaka::produce(edm::StreamID iStreamID,
                                              edm::Event &iEvent,
                                              const edm::EventSetup &iSetup) const {
  auto vertices = std::make_unique<reco::VertexCollection>();

  auto tracks = iEvent.getHandle(tracksToken_);
  auto const &indToEdm = iEvent.get(indToEdmToken_);
  auto const &vertexSoA = iEvent.get(vertexSoAToken_);
  auto const nVertices = vertexSoA.view().vertex().nVertices();
  auto const vtxTracksSoA = vertexSoA.view().tracks();

  std::vector<int32_t> vtxTrackIds;
  vtxTrackIds.reserve(64);  // avoid first relocations

  AdaptiveVertexFitter theAdaptiveFitter(GeometricAnnealing(fitParams_.sigmacut, fitParams_.tini, fitParams_.ratio),
                                         DefaultLinearizationPointFinder(),
                                         KalmanVertexUpdator<5>(),
                                         KalmanVertexTrackCompatibilityEstimator<5>(),
                                         KalmanVertexSmoother());

  edm::ESHandle<TransientTrackBuilder> trackBuilder = iSetup.getHandle(trackBuilderToken_);

  // eventually build transient tracks for refitting
  std::vector<reco::TransientTrack> transientTracks;
  if (refitVertices_) {
    for (auto track = tracks->begin(); track != tracks->end(); ++track) {
      transientTracks.push_back(tthelpers::buildTT(tracks, trackBuilder, track - tracks->begin()));
    }
  }

#if VEGA_CONVERTER_DEBUG
  printf("VegaVertexProducerFromSoAAlpaka: Everything is setup - start looping over vertices from SoA\n");
#endif

  for (auto iv = 0u; iv < nVertices; iv++) {
    vtxTrackIds.clear();

    assert(iv < nVertices);
    assert(vtxTrackIds.empty());

#if VEGA_CONVERTER_DEBUG
    printf("VegaVertexProducerFromSoAAlpaka: Start with Vertex %d\n", iv);
#endif

    auto vtx = vertexSoA.view().vertex()[iv];

    for (auto it = (iv == 0u) ? 0u : vertexSoA.view().vertex()[iv - 1u].tracksOffset(); it < vtx.tracksOffset(); it++) {
      if (vtxTracksSoA[it].weight() >= minTrackWeight_) {
        auto trackIdSoA = vtxTracksSoA[it].id();
        auto trackIdEdm = indToEdm[trackIdSoA];
        if (trackIdEdm != pixelTrack::skippedTrack)
          vtxTrackIds.push_back(trackIdEdm);

#if VEGA_CONVERTER_WARNINGS
        else
          printf(
              "VegaVertexProducerFromSoAAlpaka: Vertex %d has track %d (in SoA) which is not present in EDM - skip the "
              "track\n",
              iv,
              trackIdSoA);
#endif
      }
    }

#if VEGA_CONVERTER_DEBUG
    printf("VegaVertexProducerFromSoAAlpaka: Vertex %d - done with collecting tracks\n", iv);
#endif

    auto nt = vtxTrackIds.size();
    if (nt < 2) {
#if VEGA_CONVERTER_WARNINGS
      printf("VegaVertexProducerFromSoAAlpaka: Vertex %d with less than two effective tracks - skip it\n", iv);
#endif
      continue;
    }

    GlobalPoint position(vtx.x(), vtx.y(), vtx.z());

    TransientVertex vtxEdm;

#if VEGA_CONVERTER_DEBUG
    printf("VegaVertexProducerFromSoAAlpaka: Vertex %d - done with position\n", iv);
#endif

    if (refitVertices_) {
      // created TransientTrack vector of the vertex
      std::vector<reco::TransientTrack> vtxTracks;
      for (auto it : vtxTrackIds) {
        reco::TransientTrack &tt = (transientTracks)[it];
        vtxTracks.push_back(tt);
      }

      vtxEdm = theAdaptiveFitter.vertex(vtxTracks, position);

      if (!vtxEdm.isValid()) {
#if VEGA_CONVERTER_WARNINGS
        printf("VegaVertexProducerFromSoAAlpaka: AVF refit of vertex %d invalid - skip it\n", iv);
#endif
        continue;
      }
    }

#if VEGA_CONVERTER_DEBUG
    printf("VegaVertexProducerFromSoAAlpaka: Vertex %d - done with refitting\n", iv);
#endif

    vertices->push_back(reco::Vertex(vtxEdm));

#if VEGA_CONVERTER_DEBUG
    printf("VegaVertexProducerFromSoAAlpaka: Vertex %d - done with pushing it to the EDM collection\n", iv);
#endif
  }

  iEvent.put(std::move(vertices));

#if VEGA_CONVERTER_DEBUG
  printf("VegaVertexProducerFromSoAAlpaka: done\n");
#endif
}

DEFINE_FWK_MODULE(VegaVertexProducerFromSoAAlpaka);
