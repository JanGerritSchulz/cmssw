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

#include "VertexExplorer.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using namespace cms::alpakatools;

  class VegaVertexSoAProducer : public global::EDProducer<> {
    using TrkSoADevice = reco::TracksSoACollection;
    using Algo = vega::VertexExplorer;

  public:
    explicit VegaVertexSoAProducer(const edm::ParameterSet& iConfig);
    ~VegaVertexSoAProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::StreamID streamID, device::Event& iEvent, const device::EventSetup& iSetup) const override;

    const Algo algo_;

    const int maxVertices_;

    device::EDGetToken<TrkSoADevice> tracks_deviceToken_;
    device::EDPutToken<reco::VertexSoACollection> vertex_deviceToken_;
  };

  VegaVertexSoAProducer::VegaVertexSoAProducer(const edm::ParameterSet& iConfig)
      : EDProducer(iConfig),
        algo_(),
        maxVertices_(iConfig.getParameter<int>("maxVertices")),
        tracks_deviceToken_(consumes(iConfig.getParameter<edm::InputTag>("trackSrc"))),
        vertex_deviceToken_(produces()) {}

  void VegaVertexSoAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("trackSrc", edm::InputTag("pixelTracksAlpaka"));
    desc.add<int>("maxVertices", 100);

    descriptions.addWithDefaultLabel(desc);
  }

  void VegaVertexSoAProducer::produce(edm::StreamID streamID,
                                      device::Event& iEvent,
                                      const device::EventSetup& iSetup) const {
    auto const& hTracks = iEvent.get(tracks_deviceToken_);

    iEvent.emplace(vertex_deviceToken_, algo_.makeAsync(iEvent.queue(), hTracks.view().tracks(), maxVertices_));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(VegaVertexSoAProducer);
