// system include files
#include <cstdint>
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
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
#include <limits>

namespace {
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

  std::pair<int16_t, int16_t> getBestVertex(reco::Track const& trk, reco::VertexCollection const& vertices) {
    int16_t bestVertex{-1}, bestVertexWith3Tracks{-1};
    float dzmin = std::numeric_limits<float>::max();
    float dzminWith3Tracks = std::numeric_limits<float>::max();
    for (int16_t i{0}; auto const& vertex : vertices) {
      float dz = std::abs(trk.dz(vertex.position()));
      if (dz < dzmin) {
        bestVertex = i;
        dzmin = dz;
      }
      if ((dz < dzminWith3Tracks) && (vertex.tracksSize() >= 3)) {
        bestVertexWith3Tracks = i;
        dzminWith3Tracks = dz;
      }
      i++;
    }

    return {bestVertex, bestVertexWith3Tracks};
  }
}  // namespace

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
  const edm::EDGetTokenT<reco::VertexCollection> vertices_getToken_;
  const edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> trackAssociatorToken_;
  const edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;
  const edm::EDGetTokenT<ClusterTPAssociation> clusterTPAssociation_getToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_getToken_;

  TTree* output_tree_;

  // constants
  const static int nLayers_ = 31;

  // general event information
  edm::EventNumber_t ev_event_;
  float beamspot_x_;
  float beamspot_y_;
  float beamspot_z_;

  // TrackingParticles
  std::vector<float> trackingParticle_eta_;
  std::vector<float> trackingParticle_phi_;
  std::vector<float> trackingParticle_pt_;
  std::vector<int> trackingParticle_pid_;
  std::vector<bool> trackingParticle_isMissed_;
  std::vector<bool> trackingParticle_isSignal_;
  std::vector<bool> trackingParticle_isIntime_;
  std::vector<int> trackingParticle_nLayers_;
  std::vector<int> trackingParticle_nHits_;
  std::vector<int> trackingParticle_hitIds_;
  std::vector<std::pair<TrackingParticleRef::key_type, std::vector<int>>>
      trackingParticle_hitIds_Vector_;  // not stored directly

  // RecHits
  std::vector<float> recHit_x_;
  std::vector<float> recHit_y_;
  std::vector<float> recHit_z_;
  std::vector<uint64_t> recHit_detId_;
  std::vector<uint> recHit_moduleId_;
  std::vector<uint> recHit_layerId_;
  std::vector<int> recHit_trackingParticle_;
  std::vector<SiPixelRecHitRef::key_type> recHit_keysVector_;  // not stored directly

  // Tracks
  std::vector<float> track_eta_;
  std::vector<float> track_pt_;
  std::vector<float> track_tip_;
  std::vector<float> track_zip_;
  std::vector<float> track_qoverp_;
  std::vector<float> track_lambda_;
  std::vector<float> track_phi_;
  std::vector<float> track_dxy_;
  std::vector<float> track_dsz_;
  std::vector<float> track_cov_qoverp_;
  std::vector<float> track_cov_lambda_;
  std::vector<float> track_cov_phi_;
  std::vector<float> track_cov_dxy_;
  std::vector<float> track_cov_dsz_;
  std::vector<float> track_nchi2_;
  std::vector<float> track_ndof_;
  std::vector<int> track_nLayers_;
  std::vector<int> track_nHits_;
  std::vector<int> track_nPixelHits_;
  std::vector<int> track_nOTHits_;
  std::vector<int> track_nHitsInLayers_[nLayers_];
  std::vector<int> track_hitIds_;
  std::vector<bool> track_isFake_;
  std::vector<bool> track_isSignal_;
  std::vector<int> track_trackingParticle_;
  std::vector<int32_t> track_bestVertex_;
  std::vector<int32_t> track_bestVertexWith3Tracks_;

  // Vertices
  std::vector<float> vertex_x_;
  std::vector<float> vertex_y_;
  std::vector<float> vertex_z_;
  std::vector<int16_t> vertex_nTracks_;
};

void PixelTrackNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("hltPhase2PixelTracks"));
  desc.add<edm::InputTag>("verticesSrc", edm::InputTag("hltPhase2PixelVertices"));
  desc.add<edm::InputTag>("pixelRecHitSrc", edm::InputTag("hltSiPixelRecHits"));
  desc.add<edm::InputTag>("otRecHitSrc", edm::InputTag("hltSiPhase2RecHits"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.addUntracked<edm::InputTag>("trackAssociator", edm::InputTag("hltTrackAssociatorByHits"));
  desc.add<edm::InputTag>("clusterTPAssociationSrc", edm::InputTag("hltTPClusterProducer"));
  desc.add<edm::InputTag>("beamSpotSrc", edm::InputTag("hltOnlineBeamSpot"));

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
      vertices_getToken_(mayConsume<reco::VertexCollection>(pSet.getParameter<edm::InputTag>("verticesSrc"))),
      trackAssociatorToken_(consumes<reco::TrackToTrackingParticleAssociator>(
          pSet.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
      trackingParticleToken_(
          consumes<TrackingParticleCollection>(pSet.getParameter<edm::InputTag>("trackingParticles"))),
      clusterTPAssociation_getToken_(
          consumes<ClusterTPAssociation>(pSet.getParameter<edm::InputTag>("clusterTPAssociationSrc"))),
      beamSpot_getToken_(consumes(pSet.getParameter<edm::InputTag>("beamSpotSrc"))) {
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

  // get beamspot from the event
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpot_getToken_, beamSpot);
  if (!beamSpot.isValid())
    return;

  // get tracks
  edm::Handle<View<reco::Track>> trackCollectionHandle;
  if (!iEvent.getByToken(tracks_getToken_, trackCollectionHandle))
    return;
  const edm::View<reco::Track>& trackCollection = *trackCollectionHandle;

  // get vertices
  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByToken(vertices_getToken_, vertexCollection);

  // get cluster to TrackingParticle association
  ClusterTPAssociation const& clusterTPAssociation = iEvent.get(clusterTPAssociation_getToken_);

  // get track to TrackingParticle association
  auto const& associatorByHits = iEvent.get(trackAssociatorToken_);
  // get TrackingParticles
  auto TPCollectionH = iEvent.getHandle(trackingParticleToken_);

  // get the pixel RecHit collection from the event
  edm::Handle<SiPixelRecHitCollection> hits;
  iEvent.getByToken(pixelRecHits_getToken_, hits);
  if (!hits.isValid())
    return;

  // get the OT RecHit collection from the event
  edm::Handle<Phase2TrackerRecHit1DCollectionNew> hitsOT;
  iEvent.getByToken(otRecHits_getToken_, hitsOT);
  if (!hitsOT.isValid())
    return;

  // -------------------------------------------------------------------------------------
  //  Vertices
  // -------------------------------------------------------------------------------------
  if (vertexCollection.isValid()) {
    // loop over vertices
    for (auto const& vertex : *vertexCollection) {
      auto position = vertex.position();
      vertex_x_.push_back(position.x());
      vertex_y_.push_back(position.y());
      vertex_z_.push_back(position.z());
      vertex_nTracks_.push_back(vertex.tracksSize());
    }
  }

  // -------------------------------------------------------------------------------------
  //  TrackingParticles
  // -------------------------------------------------------------------------------------

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

  // -------------------------------------------------------------------------------------
  //  RecHits
  // -------------------------------------------------------------------------------------

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
    uint layerId = getLayerId(detId, trackerTopology_);

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
      auto moduleId = trackerGeometry_->idToDetUnit(detId)->index();
      recHit_moduleId_.push_back(moduleId);
      recHit_layerId_.push_back(layerId);

      // find associated TrackingParticles
      auto range = clusterTPAssociation.equal_range(OmniClusterRef(hit.cluster()));

      // if the RecHit has associated TrackingParticles
      int tpId = -1;
      if (range.first != range.second) {
        for (auto assocTrackingParticleIter = range.first; assocTrackingParticleIter != range.second;
             assocTrackingParticleIter++) {
          const TrackingParticleRef assocTrackingParticle = (assocTrackingParticleIter->second);

          // loop over collection of TrackingParticles and find the one of the associated TrackingParticle
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
    if ((layerIdOT <= 3) && (trackerGeometry_->getDetectorType(detIdObject) == TrackerGeometry::ModuleType::Ph2PSP)) {
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
        auto moduleId = trackerGeometry_->idToDetUnit(detId)->index();
        recHit_moduleId_.push_back(moduleId);

        // find associated TrackingParticles
        auto range = clusterTPAssociation.equal_range(OmniClusterRef(hit.cluster()));

        // if the RecHit has associated TrackingParticles
        int tpId = -1;
        if (range.first != range.second) {
          for (auto assocTrackingParticleIter = range.first; assocTrackingParticleIter != range.second;
               assocTrackingParticleIter++) {
            const TrackingParticleRef assocTrackingParticle = (assocTrackingParticleIter->second);

            // loop over collection and find the one of the associated TrackingParticle
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

  // -------------------------------------------------------------------------------------
  //  Tracks
  // -------------------------------------------------------------------------------------

  // loop over Tracks
  for (std::unordered_set<uint> layerSet{}; auto const& track : trackRefs) {
    // check if the track is associated
    int tpId = -1;
    bool isFake{true};
    bool isSignal{false};
    auto foundTP = recSimColl.find(track);
    if (foundTP != recSimColl.end()) {
      const auto& tp = foundTP->val;
      if (!tp.empty()) {
        isFake = false;
        isSignal = (tp[0].first->eventId().bunchCrossing() == 0 && tp[0].first->eventId().event() == 0);

        // find the correct TP
        for (size_t j{0}; auto const& tpHitsPair : trackingParticle_hitIds_Vector_) {
          if (tp[0].first.key() == tpHitsPair.first) {
            tpId = j;
            break;
          }
          j++;
        }
      }
    }
    track_isFake_.push_back(isFake);
    track_isSignal_.push_back(isSignal);
    track_trackingParticle_.push_back(tpId);
    track_nHits_.push_back(track->recHitsSize());
    track_eta_.push_back(track->eta());
    track_pt_.push_back(track->pt());
    track_tip_.push_back(track->dxy());
    track_zip_.push_back(track->dz());
    track_qoverp_.push_back(track->qoverp());
    track_lambda_.push_back(track->lambda());
    track_phi_.push_back(track->phi());
    track_dxy_.push_back(track->dxy());
    track_dsz_.push_back(track->dsz());
    track_cov_qoverp_.push_back(track->covariance(0, 0));
    track_cov_lambda_.push_back(track->covariance(1, 1));
    track_cov_phi_.push_back(track->covariance(2, 2));
    track_cov_dxy_.push_back(track->covariance(3, 3));
    track_cov_dsz_.push_back(track->covariance(4, 4));
    track_nchi2_.push_back(track->normalizedChi2());
    track_ndof_.push_back(track->ndof());

    // reset layerSet to be filled with RecHit layers of the current track
    layerSet.clear();
    int nPixelHits{0}, nOTHits{0};
    for (auto& nHitsInLayer : track_nHitsInLayers_)
      nHitsInLayer.push_back(0);

    // loop over the RecHits of that Track
    for (size_t j{0}; j < track->recHitsSize(); j++) {
      uint detId = track->recHit(j)->geographicalId().rawId();
      auto globalPos = track->recHit(j)->globalPosition();
      // loop over the RecHits with that detId
      for (int k = detId2firstHitIdMap[detId]; k < detId2lastHitIdMap[detId]; k++) {
        if ((std::abs(recHit_x_.at(k) - globalPos.x()) < 1e-4) && (std::abs(recHit_y_.at(k) - globalPos.y()) < 1e-4) &&
            (std::abs(recHit_z_.at(k) - globalPos.z()) < 1e-4)) {
          track_hitIds_.push_back(k);
          auto layerId = recHit_layerId_.at(k);
          layerSet.insert(layerId);
          track_nHitsInLayers_[layerId].back()++;
          (layerId > 27) ? nOTHits++ : nPixelHits++;
          break;
        }
      }
    }
    // fill the number of layers of the track
    track_nLayers_.push_back(layerSet.size());
    track_nPixelHits_.push_back(nPixelHits);
    track_nOTHits_.push_back(nOTHits);

    if (vertexCollection.isValid()) {
      // find best vertices
      auto bestVertices = getBestVertex(*track, *vertexCollection);
      track_bestVertex_.push_back(bestVertices.first);
      track_bestVertexWith3Tracks_.push_back(bestVertices.second);
    }
  }

  // finalize the TrackingParticle hitIds
  for (std::unordered_set<uint> layerSet{}; auto const& tpHitsPair : trackingParticle_hitIds_Vector_) {
    trackingParticle_nHits_.push_back(tpHitsPair.second.size());

    layerSet.clear();
    for (auto const hitId_ : tpHitsPair.second) {
      trackingParticle_hitIds_.push_back(hitId_);
      layerSet.insert(recHit_layerId_.at(hitId_));
    }
    // fill the number of layers of the trackingParticle
    trackingParticle_nLayers_.push_back(layerSet.size());
  }

  // -------------------------------------------------------------------------------------
  //  General event information
  // -------------------------------------------------------------------------------------

  // set general event information
  ev_event_ = iEvent.id().event();
  beamspot_x_ = beamSpot->x0();
  beamspot_y_ = beamSpot->y0();
  beamspot_z_ = beamSpot->z0();

  // fill the tree
  output_tree_->Fill();
}

void PixelTrackNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("PixelTracks", "Pixel Track Validation TTree");

  output_tree_->Branch("event", &ev_event_);
  output_tree_->Branch("beamspot_x", &beamspot_x_);
  output_tree_->Branch("beamspot_y", &beamspot_y_);
  output_tree_->Branch("beamspot_z", &beamspot_z_);
  output_tree_->Branch("trackingParticle_eta", &trackingParticle_eta_);
  output_tree_->Branch("trackingParticle_phi", &trackingParticle_phi_);
  output_tree_->Branch("trackingParticle_pt", &trackingParticle_pt_);
  output_tree_->Branch("trackingParticle_pid", &trackingParticle_pid_);
  output_tree_->Branch("trackingParticle_isMissed", &trackingParticle_isMissed_);
  output_tree_->Branch("trackingParticle_isSignal", &trackingParticle_isSignal_);
  output_tree_->Branch("trackingParticle_isIntime", &trackingParticle_isIntime_);
  output_tree_->Branch("trackingParticle_nLayers", &trackingParticle_nLayers_);
  output_tree_->Branch("trackingParticle_nHits", &trackingParticle_nHits_);
  output_tree_->Branch("trackingParticle_hitIds", &trackingParticle_hitIds_);
  output_tree_->Branch("recHit_x", &recHit_x_);
  output_tree_->Branch("recHit_y", &recHit_y_);
  output_tree_->Branch("recHit_z", &recHit_z_);
  output_tree_->Branch("recHit_detId", &recHit_detId_);
  output_tree_->Branch("recHit_moduleId", &recHit_moduleId_);
  output_tree_->Branch("recHit_layerId", &recHit_layerId_);
  output_tree_->Branch("recHit_trackingParticle", &recHit_trackingParticle_);
  output_tree_->Branch("track_eta", &track_eta_);
  output_tree_->Branch("track_pt", &track_pt_);
  output_tree_->Branch("track_tip", &track_tip_);
  output_tree_->Branch("track_zip", &track_zip_);
  output_tree_->Branch("track_qoverp", &track_qoverp_);
  output_tree_->Branch("track_lambda", &track_lambda_);
  output_tree_->Branch("track_phi", &track_phi_);
  output_tree_->Branch("track_dxy", &track_dxy_);
  output_tree_->Branch("track_dsz", &track_dsz_);
  output_tree_->Branch("track_cov_qoverp", &track_cov_qoverp_);
  output_tree_->Branch("track_cov_lambda", &track_cov_lambda_);
  output_tree_->Branch("track_cov_phi", &track_cov_phi_);
  output_tree_->Branch("track_cov_dxy", &track_cov_dxy_);
  output_tree_->Branch("track_cov_dsz", &track_cov_dsz_);
  output_tree_->Branch("track_nchi2", &track_nchi2_);
  output_tree_->Branch("track_ndof", &track_ndof_);
  output_tree_->Branch("track_nLayers", &track_nLayers_);
  output_tree_->Branch("track_nHits", &track_nHits_);
  output_tree_->Branch("track_nPixelHits", &track_nPixelHits_);
  output_tree_->Branch("track_nOTHits", &track_nOTHits_);
  for (int i{0}; i < nLayers_; i++)
    output_tree_->Branch(("track_nHitsInLayer" + std::to_string(i)).c_str(), &(track_nHitsInLayers_[i]));
  output_tree_->Branch("track_nOTHits", &track_nOTHits_);
  output_tree_->Branch("track_hitIds", &track_hitIds_);
  output_tree_->Branch("track_isFake", &track_isFake_);
  output_tree_->Branch("track_isSignal", &track_isSignal_);
  output_tree_->Branch("track_trackingParticle", &track_trackingParticle_);
  output_tree_->Branch("track_bestVertex", &track_bestVertex_);
  output_tree_->Branch("track_bestVertexWith3Tracks", &track_bestVertexWith3Tracks_);
  output_tree_->Branch("vertex_x", &vertex_x_);
  output_tree_->Branch("vertex_y", &vertex_y_);
  output_tree_->Branch("vertex_z", &vertex_z_);
  output_tree_->Branch("vertex_nTracks", &vertex_nTracks_);
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
  trackingParticle_nLayers_.clear();
  trackingParticle_nHits_.clear();
  trackingParticle_hitIds_.clear();
  recHit_x_.clear();
  recHit_y_.clear();
  recHit_z_.clear();
  recHit_detId_.clear();
  recHit_moduleId_.clear();
  recHit_layerId_.clear();
  recHit_trackingParticle_.clear();
  track_eta_.clear();
  track_pt_.clear();
  track_tip_.clear();
  track_zip_.clear();
  track_qoverp_.clear();
  track_lambda_.clear();
  track_phi_.clear();
  track_dxy_.clear();
  track_dsz_.clear();
  track_cov_qoverp_.clear();
  track_cov_lambda_.clear();
  track_cov_phi_.clear();
  track_cov_dxy_.clear();
  track_cov_dsz_.clear();
  track_nchi2_.clear();
  track_ndof_.clear();
  track_nLayers_.clear();
  track_nHits_.clear();
  track_nPixelHits_.clear();
  track_nOTHits_.clear();
  for (int i{0}; i < nLayers_; i++)
    track_nHitsInLayers_[i].clear();
  track_hitIds_.clear();
  track_isFake_.clear();
  track_isSignal_.clear();
  track_trackingParticle_.clear();
  track_bestVertex_.clear();
  track_bestVertexWith3Tracks_.clear();
  vertex_x_.clear();
  vertex_y_.clear();
  vertex_z_.clear();
  vertex_nTracks_.clear();

  trackingParticle_hitIds_Vector_.clear();
  recHit_keysVector_.clear();
}

DEFINE_FWK_MODULE(PixelTrackNtuplizer);