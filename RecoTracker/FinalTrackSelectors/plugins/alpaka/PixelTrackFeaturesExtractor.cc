#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "RecoTracker/FinalTrackSelectors/interface/PixelTrackFeaturesSoA.h"
#include "RecoTracker/FinalTrackSelectors/interface/PixelRecHitFeaturesSoA.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackFeaturesDeviceCollection.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelRecHiFeaturesDeviceCollection.h"
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

    const device::EDGetToken<TkSoADevice> tokenTrackIn_;
    const int32_t maxTracks_;
    const int32_t maxHitsPerTrack_;
    const int32_t minNumberOfHits_;
    const pixelTrack::Quality minQuality_;

    const device::EDPutToken<TkSoADevice> tokenTrackOut_;
    //const device::EDPutToken<PixelTrackFeaturesOnDevice> tokenTrackOut_;
    //const device::EDPutToken<PixelRecHitFeaturesOnDevice> tokenHitOut_;
  };

  PixelTrackFeaturesExtractor::PixelTrackFeaturesExtractor(
      const edm::ParameterSet& iConfig)
      : EDProducer(iConfig),
        tokenTrackIn_(consumes(iConfig.getParameter<edm::InputTag>("pixelTrackSrc"))),
        maxTracks_(iConfig.getParameter<int>("maxTracks")),
        maxHitsPerTrack_(RecHitFeatures::MaxHitsPerTrack),
        minNumberOfHits_(iConfig.getParameter<int>("minNumberOfHits")),
        minQuality_(pixelTrack::qualityByName(iConfig.getParameter<std::string>("minQuality"))),
        tokenTrackOut_(produces())
  {
    if (minQuality_ == pixelTrack::Quality::notQuality) {
      throw cms::Exception("PixelTrackConfiguration")
        << iConfig.getParameter<std::string>("minQuality") + " is not a pixelTrack::Quality";
    }
    if (minQuality_ < pixelTrack::Quality::dup) {
      throw cms::Exception("PixelTrackConfiguration")
          << iConfig.getParameter<std::string>("minQuality") + " not supported";
    }
  }

  void PixelTrackFeaturesExtractor::produce(
      device::Event& iEvent,
      const device::EventSetup&) {
    std::cout << "PixelTrackFeaturesExtractor::produce" << std::endl;

    // Retrieve tokens
    auto& queue = iEvent.queue();
    
    const auto& hits   = iEvent.get(tokenHit_);
    const auto& tracks = iEvent.get(tokenTrackIn_);

    //Instanciate the features SoA
    std::cout << "PixelTrackFeaturesExtractor::Producing features collection" << std::endl;
    PixelTrackFeaturesOnDevice features(maxTracks_, queue);
    
    //Instanciate the variables
    auto d_nKeptTracks = cms::alpakatools::make_device_buffer<int>(queue);
    auto d_newIndex = cms::alpakatools::make_device_buffer<int[]>(queue, maxTracks_);
    auto d_nHitsPerKeptTrack = cms::alpakatools::make_device_buffer<int[]>(queue, maxTracks_);
    auto h_zero = cms::alpakatools::make_host_buffer<int>();
    *h_zero = 0;
    
    alpaka::memcpy(queue, d_nKeptTracks, h_zero);
    
    //Launch first kernel to look which tracks need to be filtered out
    // based on quality criteria from the CA
    std::cout << "PixelTrackFeaturesExtractor::Launching kernel" << std::endl;
    launchTrackFeatureExtractorKernel(
      queue,
      maxTracks_,
      minNumberOfHits_,
      minQuality_,
      tracks.view(),
      features.view(),
      alpaka::getPtrNative(d_nKeptTracks),
      alpaka::getPtrNative(d_newIndex),
      alpaka::getPtrNative(d_nHitsPerKeptTrack)
    );
    std::cout << "PixelTrackFeaturesExtractor::done" << std::endl;

    //Prepare to shrink SoA
    auto h_nKeptTracks = cms::alpakatools::make_host_buffer<int>();
    auto h_nHitsPerKeptTrack = cms::alpakatools::make_host_buffer<int[]>(queue, maxTracks_);

    alpaka::memcpy(queue, h_nKeptTracks, d_nKeptTracks);
    alpaka::memcpy(queue, h_nHitsPerKeptTrack, d_nHitsPerKeptTrack);
    alpaka::wait(queue);
    const int nKeptTracks = *h_nKeptTracks;

    std::cout << "Kept tracks: " << nKeptTracks << "\n";
    TkSoADevice outTracks(nKeptTracks, queue);
    auto h_hitOffsets = cms::alpakatools::make_host_buffer<int[]>(queue, nKeptTracks);
    auto d_hitOffsets = cms::alpakatools::make_host_buffer<int[]>(queue, nKeptTracks);

    //Pre compute the new hit offsets per track
    int sum = 0;
    for (int i = 0; i < nKeptTracks; ++i) {
      sum += h_nHitsPerKeptTrack[i];
      h_hitOffsets[i] = sum;
    }

    alpaka::memcpy(queue, d_hitOffsets, h_hitOffsets);
    
    //Launch the kernel to clean the SoA
    std::cout << "Prefiltering bad tracks" << std::endl;
    launchCompactKernel(
      queue, 
      maxTracks_, 
      tracks.view(), 
      outTracks.view(), 
      alpaka::getPtrNative(d_newIndex),
      alpaka::getPtrNative(d_hitOffsets)
    );

    //iEvent.emplace(tokenTrackOut_, std::move(outTracks));
  }

  void PixelTrackFeaturesExtractor::fillDescriptions(
      edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("pixelTrackSrc", {"hltPhase2PixelTracksSoA"});
    desc.add<int>("maxTracks", 100000);
    desc.add<int>("minNumberOfHits", 0);
    desc.add<std::string>("minQuality", "tight");
    descriptions.addWithDefaultLabel(desc);
  }
};

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PixelTrackFeaturesExtractor);
