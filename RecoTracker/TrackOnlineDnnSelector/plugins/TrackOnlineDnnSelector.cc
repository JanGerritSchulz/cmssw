// -*- C++ -*-
//
// Package:    RecoTracker/TrackOnlineDnnSelector
// Class:      TrackOnlineDnnSelector
//
/**\class TrackOnlineDnnSelector TrackOnlineDnnSelector.cc RecoTracker/TrackOnlineDnnSelector/plugins/TrackOnlineDnnSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emanuele Coradin
//         Created:  Wed, 26 Nov 2025 16:11:05 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include <fstream>

//
// class declaration
//

class TrackOnlineDnnSelector : public edm::stream::EDProducer<> {
public:
  explicit TrackOnlineDnnSelector(const edm::ParameterSet&);
  ~TrackOnlineDnnSelector() override;
  using MVACollection = std::vector<float>;
  using QualityMaskCollection = std::vector<unsigned char>;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  const uint32_t maxRecHits_;
  const double probThreshold_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

  const std::unique_ptr<cms::Ort::ONNXRuntime> onnxRuntimeInstance_;
  const cms::Ort::ONNXRuntime* onnxSession_;
  const std::vector<std::string> inputNames_;
  const std::vector<std::string> output_en_;
};

namespace {
  inline int nPixelHits(reco::Track const& tk) { return tk.hitPattern().numberOfValidPixelHits(); }
  inline int nTrkLays(reco::Track const& tk) { return tk.hitPattern().trackerLayersWithMeasurement(); }
}

//
// constructors and destructor
//
TrackOnlineDnnSelector::TrackOnlineDnnSelector(const edm::ParameterSet& params) 
    : maxRecHits_(params.getParameter<uint32_t>("maxRecHits")),
      probThreshold_(params.getParameter<double>("probThreshold")),
      tracks_(consumes<std::vector<reco::Track>>(params.getParameter<edm::InputTag>("tracksSrc"))),
      beamSpot_(consumes<reco::BeamSpot>(params.getParameter<edm::InputTag>("beamSpot"))),
      onnxRuntimeInstance_(std::make_unique<cms::Ort::ONNXRuntime>(
            params.getParameter<edm::FileInPath>("onnxModelPath").fullPath().c_str())),
      inputNames_(params.getParameter<std::vector<std::string>>("inputNames")),  // Define input names for inference
      output_en_(params.getParameter<std::vector<std::string>>("output_en"))    // Define output energy for inference
  
{
  onnxSession_ = onnxRuntimeInstance_.get();
  //TODO: these are placeholders taken from the MVA producer
  // decide what's the actual output of the producer
  produces<MVACollection>("MVAValues");
  produces<QualityMaskCollection>("QualityMasks");

  //now do what ever other initialization is needed
}

TrackOnlineDnnSelector::~TrackOnlineDnnSelector() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void TrackOnlineDnnSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::LogWarning("TrackOnlineDnnSelector") << "Hello world! Let's produce";
  
  //retrieve tokens
  auto tracksIn = iEvent.getHandle(tracks_);
  auto beamSpotIn = iEvent.getHandle(beamSpot_);

 
  //initialize to NaN
  static constexpr float default_value = std::numeric_limits<float>::quiet_NaN();
  
  const int64_t nTracks = tracksIn.isValid() ? (*tracksIn).size() : 0;
  constexpr int64_t nHitFeatures = 9;
  constexpr int64_t nTrackFeatures = 25;

  std::vector<int64_t> hits_shape{nTracks, maxRecHits_, nHitFeatures};
  std::vector<int64_t> tracks_shape{nTracks, nTrackFeatures};
  std::vector<float> hits_input(nTracks * maxRecHits_ * nHitFeatures, default_value);
  std::vector<float> tracks_input(nTracks * nTrackFeatures, default_value);

  std::vector<std::vector<float>> input_Data_;
  std::vector<std::vector<int64_t>> input_shapes_ = {hits_shape, tracks_shape}; 

  // Initialize output

  auto mvas  = std::make_unique<MVACollection>(nTracks, -99.f);
  auto quals = std::make_unique<QualityMaskCollection>(nTracks, 0);

  // Escape if there are no tracks:
  if(!tracksIn.isValid() || nTracks==0){
    iEvent.put(std::move(mvas), "MVAValues");
    iEvent.put(std::move(quals), "QualityMasks");
    return ;
  }
  
  
  // prepare input data 
  
  // retrieve the BS information
  math::XYZPoint point;
  if (beamSpotIn.isValid()) {
    const auto& beamSpot = *beamSpotIn;
    point = math::XYZPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());
  }
  else {
    edm::LogWarning("TrackOnlineDnnSelector") 
      << " Invalid handle for beam spot features in tracks input collection";
    point = math::XYZPoint(0, 0, 0);
  }

  
  const auto& tracks = *tracksIn;

  for (int64_t tkIndex = 0; tkIndex < nTracks; ++tkIndex) {
    const auto& track = tracks[tkIndex];

    // Retrieve recHit features
    for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
      auto hit = *it;
      auto const& globalPoint = hit->globalPosition();
      auto const& globalError = hit->globalPositionError();
      auto hitIndex = std::distance(track.recHitsBegin(), it);
      if (hitIndex >= maxRecHits_) {
        edm::LogWarning("TrackOnlineDnnSelector")
            << " Track " << tkIndex << " has more (" << track.recHitsSize() << ") than " << maxRecHits_
            << " recHits, skipping the rest.";
        break;
      }

      auto base = tkIndex * maxRecHits_ * nHitFeatures + hitIndex * nHitFeatures;
      hits_input[base + 0] = globalPoint.x();
      hits_input[base + 1] = globalPoint.y();
      hits_input[base + 2] = globalPoint.z();
      hits_input[base + 3] = globalError.cxx();
      hits_input[base + 4] = globalError.cyy();
      hits_input[base + 5] = globalError.czz();
      hits_input[base + 6] = globalPoint.perp();
      hits_input[base + 7] = globalPoint.eta();
      hits_input[base + 8] = globalPoint.phi();
    }

    // retrieve recoPixelTrack features 
    auto tbase = tkIndex * nTrackFeatures;
    tracks_input[tbase + 0] = nPixelHits(track);
    tracks_input[tbase + 1] = nTrkLays(track);
    tracks_input[tbase + 2] = track.charge();
    tracks_input[tbase + 3] = track.chi2();
    tracks_input[tbase + 4] = track.dxy();
    tracks_input[tbase + 5] = track.dz();
    tracks_input[tbase + 6] = track.dzError();
    tracks_input[tbase + 7] = track.dsz();
    tracks_input[tbase + 8] = track.dszError();
    tracks_input[tbase + 9] = track.dxyError();
    tracks_input[tbase + 10] = track.eta();
    tracks_input[tbase + 11] = track.etaError();
    tracks_input[tbase + 12] = track.lambdaError();
    tracks_input[tbase + 13] = track.ndof();
    tracks_input[tbase + 14] = track.phi();
    tracks_input[tbase + 15] = track.phiError();
    tracks_input[tbase + 16] = track.pt();
    tracks_input[tbase + 17] = track.ptError();
    tracks_input[tbase + 18] = track.qoverp();
    tracks_input[tbase + 19] = track.qoverpError();
    tracks_input[tbase + 20] = track.vx();
    tracks_input[tbase + 21] = track.vy();
    tracks_input[tbase + 22] = track.vz();

    // retrieve Beamspot related features
    if (beamSpotIn.isValid()) {
      tracks_input[tbase + 23] = track.dz(point);
      tracks_input[tbase + 24] = track.dxy(point);
    }
  } 

  input_Data_.push_back(hits_input);
  input_Data_.push_back(tracks_input);
  auto pidOutput = onnxSession_->run(inputNames_, input_Data_, input_shapes_, output_en_, nTracks);
  const auto& probTensor = pidOutput[0];
  const float* probabilities = probTensor.data();

  // Loop over tracks, fill probability and quality mask
  for (int i = 0; i < nTracks; i++) {

      float prob = probabilities[i];
      (*mvas)[i] = prob;

      unsigned char mask = 0;

      // Only highPurity is meaningful: 1 bit when prob > threshold
      if (prob > probThreshold_) {
          mask |= (1 << reco::TrackBase::highPurity);
      }

      // loose and tight intentionally left OFF

      (*quals)[i] = mask;
  }

  // put into event
  iEvent.put(std::move(mvas), "MVAValues");
  iEvent.put(std::move(quals), "QualityMasks");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TrackOnlineDnnSelector::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TrackOnlineDnnSelector::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
TrackOnlineDnnSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
TrackOnlineDnnSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TrackOnlineDnnSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TrackOnlineDnnSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackOnlineDnnSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::FileInPath>(
            "onnxModelPath",
            edm::FileInPath("RecoTracker/TrackOnlineDnnSelector/data/track_classifier.onnx"))
        ->setComment("Path to ONNX model");
   
  desc.add<std::vector<std::string>>("inputNames", {"hits", "tracks"});
  desc.add<std::vector<std::string>>("output_en", {"probabilities"});
  desc.add<uint32_t>("maxRecHits", 16);
  desc.add<double>("probThreshold", 0.095)->setComment("DNN probability threshold for highPurity");
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("hltInitialStepTracks"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("hltOnlineBeamSpot"));
  descriptions.add("trackOnlineDnnSelector", desc);  
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackOnlineDnnSelector);
