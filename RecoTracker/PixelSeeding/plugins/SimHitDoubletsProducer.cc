#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "SimTracker/Common/interface/TrackingParticleSelector.h"           // include the selector for TrackingParticles

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"  // includes input data format: RecHit collection 
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"     // includes input data format: TrackingParticleCollection, TrackingParticle
#include "RecoTracker/TkHitPairs/interface/RecHitsSortedInPhi.h"            // includes output data format: HitDoublets, RecHitsSortedInPhi::Hit (BaseTrackerRecHit)

#include <iostream>
#include <vector>


class SimDoubletsProducer : public edm::stream::EDProducer<> {
public:
  using RecHit = RecHitsSortedInPhi::Hit;

  explicit SimDoubletsProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  TrackingParticleSelector trackingParticleSelector;

  const edm::EDGetTokenT<TrackingParticleCollection> trackingParticles_getToken_;
  const edm::EDGetTokenT<SiPixelRecHitCollection> pixelRecHits_getToken_;
  const edm::EDPutTokenT<SiPixelRecHitCollection> pixelRecHits_putToken_;
};





SimDoubletsProducer::SimDoubletsProducer(const edm::ParameterSet& pSet)
    : trackingParticles_getToken_(consumes(pSet.getParameter<edm::InputTag>("trackingParticleSrc"))),
      pixelRecHits_getToken_(consumes(pSet.getParameter<edm::InputTag>("pixelRecHitSrc"))),
      pixelRecHits_putToken_(produces<SiPixelRecHitCollection>()) {

  // initialize the selector for TrackingParticles used to create SimHitDoublets
  const edm::ParameterSet& pSetTPSel = pSet.getParameter<edm::ParameterSet>("TrackingParticleSelectionConfig");
  trackingParticleSelector = TrackingParticleSelector(pSetTPSel.getParameter<double>("ptMin"),
                                                      pSetTPSel.getParameter<double>("ptMax"),
                                                      pSetTPSel.getParameter<double>("minRapidity"),
                                                      pSetTPSel.getParameter<double>("maxRapidity"),
                                                      pSetTPSel.getParameter<double>("tip"),
                                                      pSetTPSel.getParameter<double>("lip"),
                                                      pSetTPSel.getParameter<int>("minHit"),
                                                      pSetTPSel.getParameter<bool>("signalOnly"),
                                                      pSetTPSel.getParameter<bool>("intimeOnly"),
                                                      pSetTPSel.getParameter<bool>("chargedOnly"),
                                                      pSetTPSel.getParameter<bool>("stableOnly"),
                                                      pSetTPSel.getParameter<std::vector<int>>("pdgId"),
                                                      pSetTPSel.getParameter<bool>("invertRapidityCut"),
                                                      pSetTPSel.getParameter<double>("minPhi"),
                                                      pSetTPSel.getParameter<double>("maxPhi"));
}



void SimDoubletsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  
  // sources for TrackingParticles and RecHits
  desc.add<edm::InputTag>("trackingParticleSrc", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("pixelRecHitSrc", edm::InputTag("hltSiPixelRecHits"));

  // FIXME : maybe move those settings to a separate file
  // parameter set for the selection of TrackingParticles that will be used for SimHitDoublets
  edm::ParameterSetDescription descTPSelector;
  descTPSelector.add<double>("ptMin", 0.005);
  descTPSelector.add<double>("ptMax", 1e100);
  descTPSelector.add<double>("minRapidity", -5.);
  descTPSelector.add<double>("maxRapidity", 5.);
  descTPSelector.add<double>("tip", 100.); // FIXME : for now it's just a random large value
  descTPSelector.add<double>("lip", 100.); // FIXME : for now it's just a random large value
  descTPSelector.add<int>("minHit", 0);
  descTPSelector.add<bool>("signalOnly", true);
  descTPSelector.add<bool>("intimeOnly", true);
  descTPSelector.add<bool>("chargedOnly", true);
  descTPSelector.add<bool>("stableOnly", false);
  descTPSelector.add<std::vector<int>>("pdgId", {});
  descTPSelector.add<bool>("invertRapidityCut", false);
  descTPSelector.add<double>("minPhi", -3.2);
  descTPSelector.add<double>("maxPhi", 3.2);
  desc.add<edm::ParameterSetDescription>("TrackingParticleSelectionConfig", descTPSelector);

  descriptions.addWithDefaultLabel(desc);
}


void SimDoubletsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup) {
  // get the pixel RecHit collection from the event
  SiPixelRecHitCollection const& hits = event.get(pixelRecHits_getToken_);

  // get trackingParticles from the event
  TrackingParticleCollection const& trackingParticles = event.get(trackingParticles_getToken_);

  // loop over trackingParticles
  size_t i = 0;
  for (TrackingParticle const& trackingParticle : trackingParticles) {
    if (trackingParticleSelector(trackingParticle)){
      i++;
    }
  }

  std::cout << "Size of TrackingParticleCollection : " << trackingParticles.size() << std::endl;
  std::cout << "i after looping through TrackingParticles : " << i << std::endl;
  
  

  // for now just write the RecHits to the event once more...
  SiPixelRecHitCollection output = hits;
  event.emplace(pixelRecHits_putToken_, std::move(output));
}

DEFINE_FWK_MODULE(SimDoubletsProducer);