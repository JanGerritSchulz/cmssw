#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "DataFormats/TrackSoA/interface/TracksDevice.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "RecoVertex/Vega/plugins/VegaParams.h"

#include "VegaAlgo.h"

#define VEGA_PRODUCER_DEBUG 0

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  using namespace ::vega;

  class VegaVertexSoAProducer : public stream::SynchronizingEDProducer<> {
    using TrkSoADevice = reco::TracksSoACollection;
    using VtxSoADevice = reco::VertexSoACollection;
    using Algo = vega::VegaAlgo;

  public:
    explicit VegaVertexSoAProducer(const edm::ParameterSet& iConfig);
    ~VegaVertexSoAProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void acquire(const device::Event& iEvent, const device::EventSetup& iSetup) override;
    void produce(device::Event& iEvent, const device::EventSetup& iSetup) override;

  private:
    VegaParams getVegaParams(const edm::ParameterSet& iConfig) const;

    const Algo algo_;

    const int maxVertices_;
    const int avgTracksPerVertex_;
    cms::alpakatools::host_buffer<int> nTracks_;

    device::EDGetToken<TrkSoADevice> tracksDeviceToken_;
    device::EDPutToken<VtxSoADevice> verticesDeviceToken_;
  };

  VegaVertexSoAProducer::VegaVertexSoAProducer(const edm::ParameterSet& iConfig)
      : SynchronizingEDProducer(iConfig),
        algo_(getVegaParams(iConfig)),
        maxVertices_(iConfig.getParameter<int>("maxVertices")),
        avgTracksPerVertex_(iConfig.getParameter<int>("avgTracksPerVertex")),
        nTracks_{cms::alpakatools::make_host_buffer<int, Platform>()},
        tracksDeviceToken_(consumes(iConfig.getParameter<edm::InputTag>("trackSoACollection"))),
        verticesDeviceToken_(produces()) {}

  VegaParams VegaVertexSoAProducer::getVegaParams(const edm::ParameterSet& iConfig) const {
    // configuration for track pair finding
    auto pairConfig = iConfig.getParameter<edm::ParameterSet>("pairParams");
    ::vega::PairParams pairParams = {
        .minTrkQuality = pixelTrack::qualityByName(pairConfig.getParameter<std::string>("minTrackQuality")),
        .minPt = pairConfig.getParameter<float>("minPt"),
        .maxDPhi = pairConfig.getParameter<float>("maxDPhi"),
        .maxDEta = pairConfig.getParameter<float>("maxDEta"),
        .maxDZ = pairConfig.getParameter<float>("maxDZ"),
        .maxLinDistance = pairConfig.getParameter<float>("maxLinDistance"),
        .maxChi2 = pairConfig.getParameter<float>("maxChi2"),
        .minDXYSig = pairConfig.getParameter<float>("minDXYSig"),
        .maxDXYSig = pairConfig.getParameter<float>("maxDXYSig")};

    // configuration for triplet building
    auto tripletConfig = iConfig.getParameter<edm::ParameterSet>("tripletParams");
    ::vega::TripletParams tripletParams = {.max3DDistance = tripletConfig.getParameter<float>("max3DDistance")};

    // configuration for vertex clustering
    auto vertexConfig = iConfig.getParameter<edm::ParameterSet>("vertexParams");
    ::vega::VertexParams vertexParams = {.maxChi2 = vertexConfig.getParameter<float>("maxChi2")};

    return VegaParams{.pair = pairParams, .triplet = tripletParams, .vertex = vertexParams};
  }

  void VegaVertexSoAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("trackSoACollection", edm::InputTag("pixelTracksAlpaka"));
    desc.add<int>("maxVertices", 100);
    desc.add<int>("avgTracksPerVertex", 4);

    // ----------------------------------------------
    // Pair finding configuration
    // ----------------------------------------------
    edm::ParameterSetDescription pairParams;
    pairParams.add<std::string>("minTrackQuality", "tight")
        ->setComment("Minimum track quality for both tracks in the track pair");
    pairParams.add<float>("minPt", 0.9)->setComment("Minimum track pT for the leading track in the track pair");
    pairParams.add<float>("maxDPhi", 0.4)->setComment("Maximum difference in phi for the track pair");
    pairParams.add<float>("maxDEta", 0.4)->setComment("Maximum difference in eta for the track pair");
    pairParams.add<float>("maxDZ", 0.75)
        ->setComment("Maximum difference in z0 (reference point = PCA to beamspot) for the track pair");
    pairParams.add<float>("maxLinDistance", 0.1)
        ->setComment(
            "Maximum distance for the track pair using linear approx. around reference point (= PCA to beamspot)");
    pairParams.add<float>("maxChi2", 5.0)->setComment("Maximum chi2/ndof for the vertex fit on the track pair");
    pairParams.add<float>("minDXYSig", 3.0)
        ->setComment("Minimum dxy significance of the fitted vertex of the pair (-1=disabled)");
    pairParams.add<float>("maxDXYSig", -1.0)
        ->setComment("Maximum dxy significance of the fitted vertex of the pair (-1=disabled)");
    desc.add<edm::ParameterSetDescription>("pairParams", pairParams);

    // ----------------------------------------------
    // Triplet building configuration
    // ----------------------------------------------
    edm::ParameterSetDescription tripletParams;
    tripletParams.add<float>("max3DDistance", 0.02)
        ->setComment("Maximum 3D distance for the track <-> best vertex-estimate of this triplet (using helices)");
    desc.add<edm::ParameterSetDescription>("tripletParams", tripletParams);

    // ----------------------------------------------
    // Vertex clustering configuration
    // ----------------------------------------------
    edm::ParameterSetDescription vertexParams;
    vertexParams.add<float>("maxChi2", 5.0)->setComment("Maximum chi-squared for the vertex fit");
    desc.add<edm::ParameterSetDescription>("vertexParams", vertexParams);

    descriptions.addWithDefaultLabel(desc);
  }

  void VegaVertexSoAProducer::acquire(const device::Event& iEvent, const device::EventSetup& iSetup) {
    auto& tracksOnDevice = iEvent.get(tracksDeviceToken_);
    auto nTracksDevice =
        cms::alpakatools::make_device_view<int const>(iEvent.queue(), tracksOnDevice.view().tracks().nTracks());
    alpaka::memcpy(iEvent.queue(), nTracks_, nTracksDevice);
  }

  void VegaVertexSoAProducer::produce(device::Event& iEvent, const device::EventSetup& iSetup) {
#if VEGA_PRODUCER_DEBUG
    std::cout << "VegaVertexSoAProducer::produce: aquired nTracks = " << *nTracks_.data() << std::endl;
#endif

    auto const& tracksOnDevice = iEvent.get(tracksDeviceToken_);

    iEvent.emplace(
        verticesDeviceToken_,
        algo_.makeVerticesAsync(
            iEvent.queue(), tracksOnDevice.view().tracks(), maxVertices_, avgTracksPerVertex_, *nTracks_.data()));

#if VEGA_PRODUCER_DEBUG
    std::cout << "VegaVertexSoAProducer::produce: produced vertex collection with maxVertices = " << maxVertices_
              << std::endl;
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(VegaVertexSoAProducer);
