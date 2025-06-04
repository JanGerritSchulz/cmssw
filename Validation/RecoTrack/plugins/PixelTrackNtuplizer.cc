// system include files
#include <memory>

#include "TTree.h"
#include "TFile.h"

/* ################################################################################
To be used like this:
---------------------------------

1. Create your cmsconfig with cmsDriver.py, e.g.:
```
cmsDriver.py step2 -s L1P2GT,HLT:75e33,VALIDATION:@hltValidation --conditions auto:phase2_realistic_T33 \
--datatier GEN-SIM-DIGI-RAW,DQMIO -n 10 --eventcontent FEVTDEBUGHLT,DQMIO --geometry ExtendedRun4D110 \
--era Phase2C17I13M9 --procModifier alpaka --filein file:step2.root --fileout file:step3.root --no_exec \
--process HLTX
```

2. Modify the line with the validation step from this:
```
process.validation_step = cms.EndPath(process.hltvalidation)
```
to this:
```
process.TFileService = cms.Service('TFileService', fileName=cms.string('./pixelTrackNtuples.root'))
process.load('Validation.RecoTrack.pixelTrackNtuplizer_cfi')
process.validation_step = cms.EndPath(
    process.hltvalidation +
    process.pixelTrackNtuplizer)
```

3. Run the config as usual:
```
cmsRun step2_L1P2GT_HLT_VALIDATION.py
```

4. A file named "pixelTrackNtuples.root" is created.
*/

// user include files
// #include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/DynArray.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"

#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "SimTracker/Common/interface/TrackingParticleSelector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <cstddef>

namespace simdoublets {
  // function that determines the layerId from the detId for Phase 2
  unsigned int getLayerId(DetId const& detId, const TrackerTopology* trackerTopology) {
    // number of barrel layers
    constexpr unsigned int numBarrelLayers{4};
    // number of disks per endcap
    constexpr unsigned int numEndcapDisks = 12;

    // set default to 999 (invalid)
    unsigned int layerId{99};

    if (detId.subdetId() == PixelSubdetector::PixelBarrel) {
      // subtract 1 in the barrel to get, e.g. for Phase 2, from (1,4) to (0,3)
      layerId = trackerTopology->pxbLayer(detId) - 1;
    } else if (detId.subdetId() == PixelSubdetector::PixelEndcap) {
      if (trackerTopology->pxfSide(detId) == 1) {
        // add offset in the backward endcap to get, e.g. for Phase 2, from (1,12) to (16,27)
        layerId = trackerTopology->pxfDisk(detId) + numBarrelLayers + numEndcapDisks - 1;
      } else {
        // add offest in the forward endcap to get, e.g. for Phase 2, from (1,12) to (4,15)
        layerId = trackerTopology->pxfDisk(detId) + numBarrelLayers - 1;
      }
    }
    // return the determined Id
    return layerId;
  }
}  // namespace simdoublets

// -------------------------------------------------------------------------------------------------------------
// class declaration
// -------------------------------------------------------------------------------------------------------------

class PixelTrackNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PixelTrackNtuplizer(const edm::ParameterSet&);
  ~PixelTrackNtuplizer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void clearVectors();

  // ------------ member data ------------
  TrackingParticleSelector trackingParticleSelector;

  const TrackerTopology* trackerTopology_ = nullptr;
  const TrackerGeometry* trackerGeometry_ = nullptr;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topology_getToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geometry_getToken_;
  const edm::EDGetTokenT<SiPixelRecHitCollection> pixelRecHits_getToken_;
  const edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> otRecHits_getToken_;
  const edm::EDGetTokenT<edm::View<reco::Track>> tracks_getToken_;
  const edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> trackAssociatorToken_;
  const edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;
  const edm::EDGetTokenT<ClusterTPAssociation> clusterTPAssociation_getToken_;

  TTree* output_tree_;
  // general event information

  edm::EventNumber_t ev_event_;

  // TrackingParticles
  std::vector<float> trackingParticle_eta_;
  std::vector<float> trackingParticle_phi_;
  std::vector<float> trackingParticle_pt_;
  std::vector<int> trackingParticle_pid_;
  std::vector<bool> trackingParticle_isMissed_;
  std::vector<bool> trackingParticle_isSignal_;
  std::vector<bool> trackingParticle_isIntime_;
  std::vector<int> trackingParticle_nHits_;
  std::vector<int> trackingParticle_hitIds_;
  std::vector<std::pair<TrackingParticleRef::key_type, std::vector<int>>>
      trackingParticle_hitIds_Vector_;  // not stored directly

  // RecHits
  std::vector<float> recHit_x_;
  std::vector<float> recHit_y_;
  std::vector<float> recHit_z_;
  std::vector<uint64_t> recHit_detId_;
  std::vector<uint> recHit_layerId_;
  std::vector<int> recHit_trackingParticle_;
  std::vector<SiPixelRecHitRef::key_type> recHit_keysVector_;  // not stored directly

  // Tracks
  std::vector<float> track_eta_;
  std::vector<float> track_phi_;
  std::vector<float> track_pt_;
  std::vector<float> track_tip_;
  std::vector<float> track_zip_;
  std::vector<float> track_nchi2_;
  std::vector<float> track_ndof_;
  std::vector<int> track_nHits_;
  std::vector<int> track_hitIds_;
  std::vector<bool> track_isFake_;
  std::vector<bool> track_isSignal_;
};

void PixelTrackNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("hltPhase2PixelTracks"));
  desc.add<edm::InputTag>("pixelRecHitSrc", edm::InputTag("hltSiPixelRecHits"));
  desc.add<edm::InputTag>("otRecHitSrc", edm::InputTag("hltSiPhase2RecHits"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.addUntracked<edm::InputTag>("trackAssociator", edm::InputTag("hltTrackAssociatorByHits"));
  desc.add<edm::InputTag>("clusterTPAssociationSrc", edm::InputTag("hltTPClusterProducer"));

  // parameter set for the selection of TrackingParticles that will be used for SimHitDoublets
  edm::ParameterSetDescription descTPSelector;
  descTPSelector.add<double>("ptMin", 1e-3);
  descTPSelector.add<double>("ptMax", 1e100);
  descTPSelector.add<double>("minRapidity", -4.5);
  descTPSelector.add<double>("maxRapidity", 4.5);
  descTPSelector.add<double>("tip", 20.);  // NOTE: differs from HLT MultiTrackValidator
  descTPSelector.add<double>("lip", 300.);
  descTPSelector.add<int>("minHit", 0);
  descTPSelector.add<bool>("signalOnly", false);
  descTPSelector.add<bool>("intimeOnly", false);
  descTPSelector.add<bool>("chargedOnly", true);
  descTPSelector.add<bool>("stableOnly", false);
  descTPSelector.add<std::vector<int>>("pdgId", {});
  descTPSelector.add<bool>("invertRapidityCut", false);
  descTPSelector.add<double>("minPhi", -3.2);
  descTPSelector.add<double>("maxPhi", 3.2);
  desc.add<edm::ParameterSetDescription>("TrackingParticleSelectionConfig", descTPSelector);

  descriptions.addWithDefaultLabel(desc);
}

PixelTrackNtuplizer::PixelTrackNtuplizer(const edm::ParameterSet& pSet)
    : topology_getToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      geometry_getToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      pixelRecHits_getToken_(consumes(pSet.getParameter<edm::InputTag>("pixelRecHitSrc"))),
      otRecHits_getToken_(consumes(pSet.getParameter<edm::InputTag>("otRecHitSrc"))),
      tracks_getToken_(consumes<edm::View<reco::Track>>(pSet.getParameter<edm::InputTag>("trackSrc"))),
      trackAssociatorToken_(consumes<reco::TrackToTrackingParticleAssociator>(
          pSet.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
      trackingParticleToken_(
          consumes<TrackingParticleCollection>(pSet.getParameter<edm::InputTag>("trackingParticles"))),
      clusterTPAssociation_getToken_(
          consumes<ClusterTPAssociation>(pSet.getParameter<edm::InputTag>("clusterTPAssociationSrc"))) {
  // initialize the selector for TrackingParticles used to create SimDoublets
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

void PixelTrackNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // clear all vectors
  clearVectors();

  // get tracker topology and geometry
  trackerTopology_ = &iSetup.getData(topology_getToken_);
  trackerGeometry_ = &iSetup.getData(geometry_getToken_);

  // get tracks
  edm::Handle<View<reco::Track>> trackCollectionHandle;
  if (!iEvent.getByToken(tracks_getToken_, trackCollectionHandle)) {
    return;
  }
  const edm::View<reco::Track>& trackCollection = *trackCollectionHandle;

  // get cluster to TrackingParticle association
  ClusterTPAssociation const& clusterTPAssociation = iEvent.get(clusterTPAssociation_getToken_);

  // get track to TrackingParticle association
  auto const& associatorByHits = iEvent.get(trackAssociatorToken_);
  // get TrackingParticles
  auto TPCollectionH = iEvent.getHandle(trackingParticleToken_);

  // get the pixel RecHit collection from the event
  edm::Handle<SiPixelRecHitCollection> hits;
  iEvent.getByToken(pixelRecHits_getToken_, hits);
  if (!hits.isValid()) {
    return;
  }

  // get the OT RecHit collection from the event
  edm::Handle<Phase2TrackerRecHit1DCollectionNew> hitsOT;
  iEvent.getByToken(otRecHits_getToken_, hitsOT);
  if (!hitsOT.isValid()) {
    return;
  }

  // create track references
  edm::RefToBaseVector<reco::Track> trackRefs;
  for (edm::View<reco::Track>::size_type i = 0; i < trackCollection.size(); ++i) {
    trackRefs.push_back(trackCollection.refAt(i));
  }
  // select reasonable TrackingParticles
  TrackingParticleRefVector tpCollection;
  for (size_t i = 0, size = TPCollectionH->size(); i < size; ++i) {
    auto tp = TrackingParticleRef(TPCollectionH, i);
    if (trackingParticleSelector(*tp)) {
      tpCollection.push_back(tp);
    }
  }
  // associate Tracks and TrackingParticles
  reco::RecoToSimCollection recSimColl = associatorByHits.associateRecoToSim(trackRefs, tpCollection);
  reco::SimToRecoCollection simRecColl = associatorByHits.associateSimToReco(trackRefs, tpCollection);

  // loop over TrackingParticles
  for (const TrackingParticleRef& tp : tpCollection) {
    // add the TrackingParticle key with empty vector of hit Ids
    trackingParticle_hitIds_Vector_.emplace_back(
        std::make_pair<TrackingParticleRef::key_type, std::vector<int>>(tp.key(), {}));

    // add TrackingParticle properties to output branches
    trackingParticle_eta_.push_back(tp->eta());
    trackingParticle_phi_.push_back(tp->phi());
    trackingParticle_pt_.push_back(tp->pt());
    trackingParticle_pid_.push_back(tp->pdgId());

    // check if it is associated
    auto foundTrack = simRecColl.find(tp);
    if (foundTrack != simRecColl.end()) {
      trackingParticle_isMissed_.push_back(false);
    } else {
      trackingParticle_isMissed_.push_back(true);
    }

    // check if signal and intime
    trackingParticle_isSignal_.push_back((tp->eventId().bunchCrossing() == 0 && tp->eventId().event() == 0));
    trackingParticle_isIntime_.push_back(tp->eventId().bunchCrossing() == 0);
  }

  // initialize a map from detId to index of first RecHit
  std::unordered_map<unsigned int, int> detId2firstHitIdMap;
  std::unordered_map<unsigned int, int> detId2lastHitIdMap;
  // loop over RecHits
  int hitId = 0;
  // loop over pixel RecHit collections of the different pixel modules
  for (const auto& detSet : *hits) {
    // get detector Id
    uint detId = detSet.detId();

    // determine layer Id from detector Id
    uint layerId = simdoublets::getLayerId(detId, trackerTopology_);

    // put the current hitId in the map
    detId2firstHitIdMap.insert({detId, hitId});

    // loop over RecHits
    for (auto const& hit : detSet) {
      // fill hit properties
      GlobalPoint hitPosition = hit.globalPosition();
      recHit_x_.push_back(hitPosition.x());
      recHit_y_.push_back(hitPosition.y());
      recHit_z_.push_back(hitPosition.z());
      recHit_detId_.push_back(detId);
      recHit_layerId_.push_back(layerId);

      // find associated TrackingParticles
      auto range = clusterTPAssociation.equal_range(OmniClusterRef(hit.cluster()));

      // if the RecHit has associated TrackingParticles
      int tpId = -1;
      if (range.first != range.second) {
        for (auto assocTrackingParticleIter = range.first; assocTrackingParticleIter != range.second;
             assocTrackingParticleIter++) {
          const TrackingParticleRef assocTrackingParticle = (assocTrackingParticleIter->second);

          // loop over collection of SimDoublets and find the one of the associated TrackingParticle
          for (size_t j{0}; j < trackingParticle_hitIds_Vector_.size(); j++) {
            auto& tpHitsPair = trackingParticle_hitIds_Vector_.at(j);
            if (assocTrackingParticle.key() == tpHitsPair.first) {
              tpHitsPair.second.push_back(hitId);
              tpId = j;
              break;
            }
          }
        }
      }
      recHit_trackingParticle_.push_back(tpId);
      hitId++;
    }  // end loop over RecHits
    // put the current hitId in the map
    detId2lastHitIdMap.insert({detId, hitId});
  }  // end loop over pixel RecHit collections of the different pixel modules

  // loop over OT RecHit collections of the different OT modules
  for (const auto& detSet : *hitsOT) {
    // get detector Id
    uint detId = detSet.detId();
      DetId detIdObject(detId);

    // get layerId of the OT (barrel: 1,2,3,...; endcap: 101, 201, 102, 202)
    uint layerIdOT = trackerTopology_->getOTLayerNumber(detId);

    // only use the RecHits if the module is in the accepted range of layers and one of the Phase 2 PS, p-sensor
    if ((layerIdOT <= 3) &&
        (trackerGeometry_->getDetectorType(detIdObject) == TrackerGeometry::ModuleType::Ph2PSP)) {
      // layerId corrected by number of pixel layers
      uint layerId = layerIdOT + 27;
      // put the current hitId in the map
      detId2firstHitIdMap.insert({detId, hitId});

      // loop over RecHits
      for (auto const& hit : detSet) {
        // fill hit properties
        GlobalPoint hitPosition = hit.globalPosition();
        recHit_x_.push_back(hitPosition.x());
        recHit_y_.push_back(hitPosition.y());
        recHit_z_.push_back(hitPosition.z());
        recHit_detId_.push_back(detId);
        recHit_layerId_.push_back(layerId);

        // find associated TrackingParticles
        auto range = clusterTPAssociation.equal_range(OmniClusterRef(hit.cluster()));

        // if the RecHit has associated TrackingParticles
        int tpId = -1;
        if (range.first != range.second) {
          for (auto assocTrackingParticleIter = range.first; assocTrackingParticleIter != range.second;
               assocTrackingParticleIter++) {
            const TrackingParticleRef assocTrackingParticle = (assocTrackingParticleIter->second);

            // loop over collection of SimDoublets and find the one of the associated TrackingParticle
            for (size_t j{0}; j < trackingParticle_hitIds_Vector_.size(); j++) {
              auto& tpHitsPair = trackingParticle_hitIds_Vector_.at(j);
              if (assocTrackingParticle.key() == tpHitsPair.first) {
                tpHitsPair.second.push_back(hitId);
                tpId = j;
                break;
              }
            }
          }
        }
        recHit_trackingParticle_.push_back(tpId);
        hitId++;
      }  // end loop over RecHits
      // put the current hitId in the map
      detId2lastHitIdMap.insert({detId, hitId});
    }
  }  // end loop over OT RecHit collections of the different pixel modules

  // loop over Tracks
  for (auto const& track : trackRefs) {
    // check if the track is associated
    bool isFake{true};
    bool isSignal{false};
    auto foundTP = recSimColl.find(track);
    if (foundTP != recSimColl.end()) {
      const auto& tp = foundTP->val;
      if (!tp.empty()) {
        isFake = false;
        isSignal = (tp[0].first->eventId().bunchCrossing() == 0 && tp[0].first->eventId().event() == 0);
      }
    }
    track_isFake_.push_back(isFake);
    track_isSignal_.push_back(isSignal);
    track_nHits_.push_back(track->recHitsSize());
    track_eta_.push_back(track->eta());
    track_phi_.push_back(track->phi());
    track_pt_.push_back(track->pt());
    track_tip_.push_back(track->dxy());
    track_zip_.push_back(track->dz());
    track_nchi2_.push_back(track->normalizedChi2());
    track_ndof_.push_back(track->ndof());

    // loop over the RecHits of that Track
    for (size_t j{0}; j < track->recHitsSize(); j++) {
      uint detId = track->recHit(j)->geographicalId().rawId();
      auto globalPos = track->recHit(j)->globalPosition();
      // loop over the RecHits with that detId
      for (int k = detId2firstHitIdMap[detId]; k < detId2lastHitIdMap[detId]; k++) {
        if ((std::abs(recHit_x_.at(k) - globalPos.x()) < 1e-4) && (std::abs(recHit_y_.at(k) - globalPos.y()) < 1e-4) &&
            (std::abs(recHit_z_.at(k) - globalPos.z()) < 1e-4)) {
          track_hitIds_.push_back(k);
          break;
        }
      }
    }
  }

  // finalize the TrackingParticle hitIds
  for (auto const& tpHitsPair : trackingParticle_hitIds_Vector_) {
    trackingParticle_nHits_.push_back(tpHitsPair.second.size());
    for (auto const hitId_ : tpHitsPair.second) {
      trackingParticle_hitIds_.push_back(hitId_);
    }
  }

  // set general event information
  ev_event_ = iEvent.id().event();

  // fill the tree
  output_tree_->Fill();
}

void PixelTrackNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("PixelTracks", "Pixel Track Validation TTree");

  output_tree_->Branch("event", &ev_event_);
  output_tree_->Branch("trackingParticle_eta", &trackingParticle_eta_);
  output_tree_->Branch("trackingParticle_phi", &trackingParticle_phi_);
  output_tree_->Branch("trackingParticle_pt", &trackingParticle_pt_);
  output_tree_->Branch("trackingParticle_pid", &trackingParticle_pid_);
  output_tree_->Branch("trackingParticle_isMissed", &trackingParticle_isMissed_);
  output_tree_->Branch("trackingParticle_isSignal", &trackingParticle_isSignal_);
  output_tree_->Branch("trackingParticle_isIntime", &trackingParticle_isIntime_);
  output_tree_->Branch("trackingParticle_nHits", &trackingParticle_nHits_);
  output_tree_->Branch("trackingParticle_hitIds", &trackingParticle_hitIds_);
  output_tree_->Branch("recHit_x", &recHit_x_);
  output_tree_->Branch("recHit_y", &recHit_y_);
  output_tree_->Branch("recHit_z", &recHit_z_);
  output_tree_->Branch("recHit_detId", &recHit_detId_);
  output_tree_->Branch("recHit_layerId", &recHit_layerId_);
  output_tree_->Branch("recHit_trackingParticle", &recHit_trackingParticle_);
  output_tree_->Branch("track_eta", &track_eta_);
  output_tree_->Branch("track_phi", &track_phi_);
  output_tree_->Branch("track_pt", &track_pt_);
  output_tree_->Branch("track_tip", &track_tip_);
  output_tree_->Branch("track_zip", &track_zip_);
  output_tree_->Branch("track_nchi2", &track_nchi2_);
  output_tree_->Branch("track_ndof", &track_ndof_);
  output_tree_->Branch("track_nHits", &track_nHits_);
  output_tree_->Branch("track_hitIds", &track_hitIds_);
  output_tree_->Branch("track_isFake", &track_isFake_);
  output_tree_->Branch("track_isSignal", &track_isSignal_);
}

void PixelTrackNtuplizer::endJob() {}

void PixelTrackNtuplizer::clearVectors() {
  trackingParticle_eta_.clear();
  trackingParticle_phi_.clear();
  trackingParticle_pt_.clear();
  trackingParticle_pid_.clear();
  trackingParticle_isMissed_.clear();
  trackingParticle_isSignal_.clear();
  trackingParticle_isIntime_.clear();
  trackingParticle_nHits_.clear();
  trackingParticle_hitIds_.clear();
  recHit_x_.clear();
  recHit_y_.clear();
  recHit_z_.clear();
  recHit_detId_.clear();
  recHit_layerId_.clear();
  recHit_trackingParticle_.clear();
  track_eta_.clear();
  track_phi_.clear();
  track_pt_.clear();
  track_tip_.clear();
  track_zip_.clear();
  track_nchi2_.clear();
  track_ndof_.clear();
  track_nHits_.clear();
  track_hitIds_.clear();
  track_isFake_.clear();
  track_isSignal_.clear();

  trackingParticle_hitIds_Vector_.clear();
  recHit_keysVector_.clear();
}

DEFINE_FWK_MODULE(PixelTrackNtuplizer);
