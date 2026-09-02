#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
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
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"

#include "VegaAlgo.h"
#include "VegaParams.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;
  using namespace ::vega;

  class VegaVertexSoAProducer : public global::EDProducer<> {
    using TrkSoADevice = reco::TracksSoACollection;
    using Algo = vega::VegaAlgo;

  public:
    explicit VegaVertexSoAProducer(const edm::ParameterSet& iConfig);
    ~VegaVertexSoAProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::StreamID streamID, device::Event& iEvent, const device::EventSetup& iSetup) const override;
    VegaParams getVegaParams(const edm::ParameterSet& iConfig) const;

    const Algo algo_;

    const int maxVertices_;

    device::EDGetToken<TrkSoADevice> tracks_deviceToken_;
    device::EDPutToken<reco::VertexSoACollection> vertex_deviceToken_;
  };

  VegaVertexSoAProducer::VegaVertexSoAProducer(const edm::ParameterSet& iConfig)
      : EDProducer(iConfig),
        algo_(getVegaParams(iConfig)),
        maxVertices_(iConfig.getParameter<int>("maxVertices")),
        tracks_deviceToken_(consumes(iConfig.getParameter<edm::InputTag>("trackSrc"))),
        vertex_deviceToken_(produces()) {}

  VegaParams VegaVertexSoAProducer::getVegaParams(const edm::ParameterSet& iConfig) const {
    // configuration for track pair finding
    auto pairConfig = iConfig.getParameter<edm::ParameterSet>("pairParams");
    ::vega::PairParams pairParams = {.minPt = pairConfig.getParameter<float>("minPt"),
                                     .maxDPhi = pairConfig.getParameter<float>("maxDPhi"),
                                     .maxDEta = pairConfig.getParameter<float>("maxDEta"),
                                     .maxDZ = pairConfig.getParameter<float>("maxDZ"),
                                     .maxLinDistance = pairConfig.getParameter<float>("maxLinDistance"),
                                     .max3DDistance = pairConfig.getParameter<float>("max3DDistance")};

    // configuration for triplet building
    auto tripletConfig = iConfig.getParameter<edm::ParameterSet>("tripletParams");
    ::vega::TripletParams tripletParams = {.max3DDistance = tripletConfig.getParameter<float>("max3DDistance")};

    return VegaParams{.pair = pairParams, .triplet = tripletParams};
  }

  void VegaVertexSoAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("trackSrc", edm::InputTag("pixelTracksAlpaka"));
    desc.add<int>("maxVertices", 100);

    // ----------------------------------------------
    // Pair finding configuration
    // ----------------------------------------------
    edm::ParameterSetDescription pairParams;
    pairParams.add<float>("minPt", 0.9)->setComment("Minimum track pT for the leading track in the track pair");
    pairParams.add<float>("maxDPhi", 0.4)->setComment("Maximum difference in phi for the track pair");
    pairParams.add<float>("maxDEta", 0.4)->setComment("Maximum difference in eta for the track pair");
    pairParams.add<float>("maxDZ", 0.75)
        ->setComment("Maximum difference in z0 (reference point = PCA to beamspot) for the track pair");
    pairParams.add<float>("maxLinDistance", 0.1)
        ->setComment(
            "Maximum distance for the track pair using linear approx. around reference point (= PCA to beamspot)");
    pairParams.add<float>("max3DDistance", 0.02)->setComment("Maximum 3D distance for the track pair (using helices)");
    desc.add<edm::ParameterSetDescription>("pairParams", pairParams);

    // ----------------------------------------------
    // Triplet building configuration
    // ----------------------------------------------
    edm::ParameterSetDescription tripletParams;
    tripletParams.add<float>("max3DDistance", 0.02)
        ->setComment("Maximum 3D distance for the track <-> best vertex-estimate of this triplet (using helices)");
    desc.add<edm::ParameterSetDescription>("tripletParams", tripletParams);

    descriptions.addWithDefaultLabel(desc);
  }

  void VegaVertexSoAProducer::produce(edm::StreamID streamID,
                                      device::Event& iEvent,
                                      const device::EventSetup& iSetup) const {
    auto const& hTracks = iEvent.get(tracks_deviceToken_);

    iEvent.emplace(vertex_deviceToken_, algo_.makeAsync(iEvent.queue(), hTracks.view().tracks(), maxVertices_));

    std::cout << "VegaVertexSoAProducer::produce: produced vertex collection with maxVertices = " << maxVertices_
              << std::endl;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(VegaVertexSoAProducer);
