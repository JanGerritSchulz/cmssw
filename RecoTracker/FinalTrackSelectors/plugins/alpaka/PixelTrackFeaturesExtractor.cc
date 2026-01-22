#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/alpaka/BeamSpotDevice.h"
#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "RecoTracker/FinalTrackSelectors/interface/PixelTrackFeaturesSoA.h"
#include "RecoTracker/FinalTrackSelectors/interface/PixelRecHitFeaturesSoA.h"
#include "RecoTracker/FinalTrackSelectors/interface/alpaka/PixelTrackFeaturesDeviceCollection.h"
#include "RecoTracker/FinalTrackSelectors/interface/alpaka/PixelRecHiFeaturesDeviceCollection.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackFeaturesExtractorKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PixelTrackFeaturesExtractor : public stream::EDProducer<> {
    using TkSoADevice = reco::TracksSoACollection;
    using HitsOnDevice = reco::TrackingRecHitsSoACollection;

  public:
    explicit PixelTrackFeaturesExtractor(const edm::ParameterSet&);
    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    void produce(device::Event&, const device::EventSetup&) override;

    const device::EDGetToken<BeamSpotDevice> tBeamSpot_;
    const device::EDGetToken<HitsOnDevice> tokenHit_;
    const device::EDGetToken<TkSoADevice> tokenTrackIn_;

    const device::EDPutToken<PixelTrackFeaturesOnDevice> tokenTrackOut_;
    const device::EDPutToken<PixelRecHitFeaturesOnDevice> tokenHitOut_;
  };

  PixelTrackFeaturesExtractor::PixelTrackFeaturesExtractor(
      const edm::ParameterSet& iConfig)
      : EDProducer(iConfig),
        tBeamSpot_(consumes(iConfig.getParameter<edm::InputTag>("beamSpot"))),
        tokenHit_(consumes(iConfig.getParameter<edm::InputTag>("pixelRecHitSrc"))),
        tokenTrackIn_(consumes(iConfig.getParameter<edm::InputTag>("pixelTrackSrc"))),
        tokenTrackOut_(produces()),
        tokenHitOut_(produces()){}

  void PixelTrackFeaturesExtractor::produce(
      device::Event& iEvent,
      const device::EventSetup&) {

    const auto& bs     = iEvent.get(tBeamSpot_);
    const auto& hits   = iEvent.get(tokenHit_);
    const auto& tracks = iEvent.get(tokenTrackIn_);

    auto& queue = iEvent.queue();
      
    const uint32_t nTracks = tracks.view().nTracks();

    auto features =
        std::make_unique<PixelTrackFeaturesOnDevice>(nTracks, queue);

    launchTrackFeatureExtractorKernel(
      queue,
      tracks.view(),
      features->view()
    );

    iEvent.put(tokenTrackOut_, std::move(features));
  }

  void PixelTrackFeaturesExtractor::fillDescriptions(
      edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("beamSpot", {"hltPhase2OnlineBeamSpotDevice"});
    desc.add<edm::InputTag>("pixelRecHitSrc", {"hltPhase2SiPixelRecHitsSoA"});
    desc.add<edm::InputTag>("pixelTrackSrc", {"hltPhase2PixelTracksCAExtension"});
    descriptions.addWithDefaultLabel(desc);
  }
};

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PixelTrackFeaturesExtractor);
