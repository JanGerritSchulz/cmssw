#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TLorentzVector.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include <vector>
#include <unordered_set>
#include <limits>
#include <tuple>
#include <cmath>

#include <set>
#include <iostream>
#include <iomanip>
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

// ── Remonte vers les ancêtres ──────────────────────────────────────────────
void printAncestors(const HepMC::GenVertex* vertex, int depth, std::set<int>& visited) {
  if (!vertex)
    return;
  if (visited.count(vertex->barcode()))
    return;
  visited.insert(vertex->barcode());

  for (auto mother = vertex->particles_in_const_begin(); mother != vertex->particles_in_const_end(); ++mother) {
    // Remonter d'abord (pour afficher du plus ancien au plus récent)
    const HepMC::GenVertex* prodVtx = (*mother)->production_vertex();
    printAncestors(prodVtx, depth - 1, visited);

    std::string indent(std::max(0, depth) * 4, ' ');
    std::cout << indent << "PDG=" << std::setw(6) << (*mother)->pdg_id() << "  status=" << (*mother)->status()
              << "  pT=" << std::setw(8) << (*mother)->momentum().perp();
    if (prodVtx)
      std::cout << "  position=(" << std::setw(8) << prodVtx->position().x() << ", " << std::setw(8)
                << prodVtx->position().y() << ", " << std::setw(8) << prodVtx->position().z() << ")";
    std::cout << std::endl;
  }
}

// ── Descend vers les descendants ───────────────────────────────────────────
void printDescendants(const HepMC::GenVertex* vertex, int depth, std::set<int>& visited) {
  if (!vertex)
    return;
  if (visited.count(vertex->barcode()))
    return;
  visited.insert(vertex->barcode());

  for (auto daughter = vertex->particles_out_const_begin(); daughter != vertex->particles_out_const_end(); ++daughter) {
    std::string indent(depth * 4, ' ');
    std::cout << indent << "|-- PDG=" << std::setw(6) << (*daughter)->pdg_id() << "  status=" << (*daughter)->status()
              << "  pT=" << std::setw(8) << (*daughter)->momentum().perp() << "  position=(" << std::setw(8)
              << vertex->position().x() << ", " << std::setw(8) << vertex->position().y() << ", " << std::setw(8)
              << vertex->position().z() << ")" << std::endl;

    const HepMC::GenVertex* endVtx = (*daughter)->end_vertex();
    if (endVtx) {
      std::cout << "  --> vtx " << endVtx->barcode() << std::endl;
      printDescendants(endVtx, depth + 1, visited);
    } else {
      std::cout << "  [stable]" << std::endl;
    }
  }
}

// ── Point d'entrée principal ───────────────────────────────────────────────
void printFullTree(const HepMC::GenVertex* vertex) {
  if (!vertex)
    return;

  // --- Ancêtres ---
  std::cout << "\n╔══ ANCESTORS ═══════════════════════════════╗\n";
  std::set<int> visitedUp;
  // depth=0 sera la particule mère directe, on remonte avec des valeurs négatives
  // qu'on recentre ensuite ; plus simple : on compte d'abord la profondeur max
  int ancestorDepth = 0;
  {
    // Passe rapide pour estimer la profondeur de la chaîne montante
    const HepMC::GenVertex* v = vertex;
    while (v && v->particles_in_const_begin() != v->particles_in_const_end()) {
      ++ancestorDepth;
      const HepMC::GenParticle* m = *v->particles_in_const_begin();
      v = m->production_vertex();
      if (ancestorDepth > 50)
        break;  // garde-fou
    }
  }
  printAncestors(vertex, ancestorDepth, visitedUp);

  // --- Vertex courant ---
  std::cout << "╠══ THIS VERTEX (barcode " << vertex->barcode() << ") ══════════╣\n";

  // --- Descendants ---
  std::cout << "╚══ DESCENDANTS ═════════════════════════════╝\n";
  std::set<int> visitedDown;
  printDescendants(vertex, 0, visitedDown);
  std::cout << std::endl;
}

//
// class declaration
//

class TrackingVertexTableProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackingVertexTableProducer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  int checkPDG(int abs_pdg) const;

  std::optional<std::tuple<float, float, float>> isAncestor(const reco::Candidate* mother,
                                                            const reco::Candidate* daughter) const;

  std::vector<std::vector<float>> computeDistanceMatrix(const std::vector<float>& SV_x,
                                                        const std::vector<float>& SV_y,
                                                        const std::vector<float>& SV_z,
                                                        const std::vector<float>& TrackingVertex_x,
                                                        const std::vector<float>& TrackingVertex_y,
                                                        const std::vector<float>& TrackingVertex_z);
  void printDistanceMatrix(const std::vector<std::vector<float>>& distances);
  std::pair<std::vector<int>, std::vector<float>> matchHadronsToSV(
      std::vector<std::vector<float>> distances,
      const std::vector<float>& SVtrk_pt,
      const std::vector<float>& SVtrk_eta,
      const std::vector<float>& SVtrk_phi,
      const std::vector<int>& SVtrk_SVidx,
      const std::vector<float>& Daughters_pt,   //genparticles
      const std::vector<float>& Daughters_eta,  //genparticles
      const std::vector<float>& Daughters_phi,  //genparticles
      const std::vector<int>& Daughters_GVidx,  // hadron index per daughter
      int n_Hadrons,
      int nRequiredCommonTracks,
      double dR_max,
      double relPt_max);

  const edm::EDGetTokenT<edm::HepMCProduct> hepmcToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  const edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  const edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  const edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  const edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  const edm::EDGetTokenT<reco::VertexToTrackingVertexAssociator<std::vector<reco::Vertex>>> vertexAssociatorToken_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
  const edm::EDGetTokenT<edm::View<reco::Candidate>> genToken_;
  const edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> svToken_;
  int nRequiredCommonTracks_;
  double dlenSigMin_;
  double dR_max_;
  double relPt_max_;
};

TrackingVertexTableProducer::TrackingVertexTableProducer(const edm::ParameterSet& iConfig)
    : hepmcToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("moduleLabelHepMC"))),
      simTrackToken_(
          consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simTrackCollection"))),
      trackingParticleCollectionToken_(consumes<TrackingParticleCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackingParticleCollection"))),
      trackingVertexCollectionToken_(
          consumes<TrackingVertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingVertexCollection"))),
      simToRecoAssociationToken_(
          consumes<reco::SimToRecoCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
      recoToSimAssociationToken_(
          consumes<reco::RecoToSimCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
      vertexAssociatorToken_(consumes<reco::VertexToTrackingVertexAssociator<std::vector<reco::Vertex>>>(
          iConfig.getUntrackedParameter<edm::InputTag>("vertexAssociator"))),
      pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
      genToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
      svToken_(consumes<std::vector<reco::VertexCompositePtrCandidate>>(
          iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
      nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
      dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
      dR_max_(iConfig.getParameter<double>("dR_max")),
      relPt_max_(iConfig.getParameter<double>("relPt_max")) {
  // produces<nanoaod::FlatTable>("GenPV");
  produces<nanoaod::FlatTable>("TV");
  produces<nanoaod::FlatTable>("TP");
  // produces<nanoaod::FlatTable>("GVDaughters");
  // produces<nanoaod::FlatTable>("SVDaughters");
}

// ------------ fill 'descriptions' with the allowed parameters for the module ------------
void TrackingVertexTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<double>("dlenMin", 0.0)->setComment("minimum value of dl to call a secondary vertex good");
  desc.add<double>("dlenSigMin", 3.0)->setComment("minimum value of dl significance to call a secondary vertex good");
  desc.add<double>("dR_max", 3.0);
  desc.add<double>("relPt_max", 3.0);
  desc.add<int>("nRequiredCommonTracks", 2);

  desc.add<edm::InputTag>("moduleLabelHepMC", edm::InputTag("generatorSmeared"))
      ->setComment("Input generated HepMC event after vtx smearing");
  desc.addUntracked<edm::InputTag>("simTrackCollection", edm::InputTag("g4SimHits"));
  desc.addUntracked<edm::InputTag>("trackingParticleCollection", edm::InputTag("mix", "MergedTrackTruth"));
  desc.addUntracked<edm::InputTag>("trackingVertexCollection", edm::InputTag("mix", "MergedTrackTruth"));
  desc.addUntracked<edm::InputTag>("trackAssociatorMap", edm::InputTag("tpToHLTGeneralTrackAssociation"));
  desc.addUntracked<edm::InputTag>("vertexAssociator",
                                   edm::InputTag("hltPVAssociatorByPositionAndTracks4GeneralTracks"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("mergedGenParticles"));
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("hltOfflinePrimaryVertices"))
      ->setComment("std::vector<reco::Vertex> and ValueMap<float> primary vertex input collections");
  desc.add<edm::InputTag>("pfSrc", edm::InputTag("hltParticleFlowTmp"))
      ->setComment("reco::PFCandidateCollection PF candidates input collections");
  desc.add<edm::InputTag>("secondaryVertices", edm::InputTag("hltDeepInclusiveMergedVerticesPF"))
      ->setComment("reco::VertexCompositePtrCandidate compatible secondary vertex input collection");

  descriptions.addWithDefaultLabel(desc);
}

void TrackingVertexTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  // const reco::RecoToSimCollection* r2s_;
  const reco::SimToRecoCollection* s2r_;

  edm::Handle<edm::HepMCProduct> mcEvtHandle;
  iEvent.getByToken(hepmcToken_, mcEvtHandle);
  const HepMC::GenEvent* evt = mcEvtHandle->GetEvent();

  edm::Handle<edm::View<reco::Candidate>> genHandle;
  iEvent.getByToken(genToken_, genHandle);
  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> svHandle;
  iEvent.getByToken(svToken_, svHandle);
  auto pvsIn = iEvent.getHandle(pvs_);

  // Get TrackingParticles + TrackingVertices and associators
  edm::Handle<edm::SimTrackContainer> simTracksH;
  iEvent.getByToken(simTrackToken_, simTracksH);
  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  if (!TPCollectionH.isValid())
    edm::LogWarning("PrimaryVertexAnalyzer4PUSlimmed") << "TPCollectionH is not valid";

  edm::Handle<TrackingVertexCollection> TVCollectionH;
  iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);
  if (!TVCollectionH.isValid())
    edm::LogWarning("PrimaryVertexAnalyzer4PUSlimmed") << "TVCollectionH is not valid";

  edm::Handle<reco::SimToRecoCollection> simToRecoH;
  iEvent.getByToken(simToRecoAssociationToken_, simToRecoH);
  if (simToRecoH.isValid())
    s2r_ = simToRecoH.product();
  else
    edm::LogWarning("PrimaryVertexAnalyzer4PUSlimmed") << "simToRecoH is not valid";

  // edm::Handle<reco::RecoToSimCollection> recoToSimH;
  // iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
  // if (recoToSimH.isValid())
  //   r2s_ = recoToSimH.product();
  // else
  //   edm::LogWarning("PrimaryVertexAnalyzer4PUSlimmed") << "recoToSimH is not valid";

  // Vertex associator
  edm::Handle<reco::VertexToTrackingVertexAssociator<std::vector<reco::Vertex>>> vertexAssociatorH;
  iEvent.getByToken(vertexAssociatorToken_, vertexAssociatorH);
  if (!vertexAssociatorH.isValid()) {
    edm::LogWarning("PrimaryVertexAnalyzer4PUSlimmed") << "vertexAssociatorH is not valid";
    return;
  }
  const reco::VertexToTrackingVertexAssociator<std::vector<reco::Vertex>>& vertexAssociator =
      *(vertexAssociatorH.product());

  // reco::VertexRecoToSimCollection vertex_r2s = vertexAssociator.associateRecoToSim(svHandle, TVCollectionH);
  // reco::VertexSimToRecoCollection vertex_s2r = vertexAssociator.associateSimToReco(svHandle, TVCollectionH);

  const auto& genParticles = genHandle;
  const auto& secondaryVertices = svHandle;

  // Output vectors
  std::vector<float> TrackingVertex_pt, TrackingVertex_eta, TrackingVertex_phi;
  std::vector<float> TrackingVertex_dlen, TrackingVertex_dxy;
  std::vector<float> SV_x, SV_y, SV_z;
  std::vector<float> TrackingVertex_x, TrackingVertex_y, TrackingVertex_z;
  std::vector<float> TrackingVertex_x_i, TrackingVertex_y_i, TrackingVertex_z_i;
  std::vector<int> TrackingVertex_pdgClass, TrackingVertex_isB, TrackingVertex_isD, TrackingVertex_isSig;
  std::vector<int> TrackingVertex_pdgId, TrackingVertex_nDaughters, TrackingVertex_nGoodDaughters,
      TrackingVertex_nGenVertices, TrackingVertex_nNearestGenVertices, TrackingVertex_nSimVertices,
      TrackingVertex_nTrackingParticles;
  std::vector<float> TrackingParticles_pt, TrackingParticles_eta, TrackingParticles_phi;
  std::vector<float> Daughters_pt, Daughters_eta, Daughters_phi;
  std::vector<int> Daughters_charge, Daughters_GVidx, TrackingParticles_pdgId;
  VertexDistance3D vdist;

  int ngv = 0;
  for (TrackingVertexCollection::const_iterator v = TVCollectionH->begin(); v != TVCollectionH->end(); ++v) {
    // require vertex to be in-time
    if (v->eventId().bunchCrossing() != 0)
      continue;

    // // could be a new vertex, check  all primaries found so far to avoid
    // // multiple entries
    // simPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z());
    // sv.eventId = v->eventId();
    // sv.sim_vertex = TrackingVertexRef(TVCollectionH, std::distance(TVCollectionH->begin(), v));

#ifdef TRACKINGVERTEX_DEBUG
    for (TrackingParticleRefVector::iterator iTrack = v->daughterTracks_begin(); iTrack != v->daughterTracks_end();
         ++iTrack) {
      // TODO(rovere) isn't it always the case? Is it really worth
      // checking this out?
      // sv.eventId = (**iTrack).eventId();
      assert((**iTrack).eventId().bunchCrossing() == 0);
    }
#endif

    // // TODO(rovere) maybe get rid of this old logic completely ... ?
    // simPrimaryVertex* vp = nullptr;  // will become non-NULL if a vertex
    //                                  // is found and then point to it
    // for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
    //   if ((sv.eventId == v0->eventId) && (fabs(sv.x - v0->x) < 1e-5) && (fabs(sv.y - v0->y) < 1e-5) &&
    //       (fabs(sv.z - v0->z) < 1e-5)) {
    //     vp = &(*v0);
    //     break;
    //   }
    // }
    // if (!vp) {
    //   // this is a new vertex, add it to the list of sim-vertices
    //   simpv.push_back(sv);
    //   vp = &simpv.back();
    //   if (verbose_) {
    //     std::cout << "this is a new vertex " << sv.eventId.event() << "   " << sv.x << " " << sv.y << " " << sv.z
    //               << std::endl;
    //   }
    // } else {
    //   if (verbose_) {
    //     std::cout << "this is not a new vertex" << sv.x << " " << sv.y << " " << sv.z << std::endl;
    //   }
    // }

    constexpr bool use_only_charged_tracks_ = true;
    constexpr bool use_reconstructable_simvertices_ = true;
    constexpr int reco_tracks_for_reconstructable_simvertices_ = 2;
    constexpr double minPt_ = 0.5;
    constexpr double maxEta_ = 4.5;

    int num_matched_reco_tracks{0}, nGenTrk{0};
    // Loop over daughter track(s) as Tracking Particles
    for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
      // auto momentum = (*(*iTP)).momentum();
      const reco::Track* matched_best_reco_track = nullptr;
      // double match_quality = -1;
      // if (use_only_charged_tracks_ && (**iTP).charge() == 0)
      //   continue;
      if (s2r_->find(*iTP) != s2r_->end()) {
        matched_best_reco_track = (*s2r_)[*iTP][0].first.get();
        // match_quality = (*s2r_)[*iTP][0].second;
      }
      // vp->ptot.setPx(vp->ptot.x() + momentum.x());
      // vp->ptot.setPy(vp->ptot.y() + momentum.y());
      // vp->ptot.setPz(vp->ptot.z() + momentum.z());
      // vp->ptot.setE(vp->ptot.e() + (**iTP).energy());
      // vp->ptsq += ((**iTP).pt() * (**iTP).pt());
      if (matched_best_reco_track) {
        num_matched_reco_tracks++;
      }
      // TODO(rovere) get rid of cuts on sim-tracks
      // TODO(rovere) be consistent between simulated tracks and
      // reconstructed tracks selection
      // count relevant particles
      if (((**iTP).pt() > minPt_) && (fabs((**iTP).eta()) < maxEta_) && (**iTP).charge() != 0) {
        nGenTrk++;
      }
    }  // End of for loop on daughters sim-particles

    // Remove the SimVertex if I cannot reconstruct it 'cause I miss at the very least reco_tracks_for_reconstructable_simvertices_ tracks
    if (use_reconstructable_simvertices_ && nGenTrk < reco_tracks_for_reconstructable_simvertices_)
      continue;

    bool isSig = (v->eventId().event() == 0);

    ngv++;

    TrackingVertex_x.push_back(v->position().x());
    TrackingVertex_y.push_back(v->position().y());
    TrackingVertex_z.push_back(v->position().z());
    TrackingVertex_nDaughters.push_back(nGenTrk);
    TrackingVertex_nGoodDaughters.push_back(num_matched_reco_tracks);
    TrackingVertex_nGenVertices.push_back(v->nGenVertices());
    TrackingVertex_nSimVertices.push_back(v->nG4Vertices());
    TrackingVertex_isSig.push_back((int)isSig);

    int pdgId = 0, nGVs = 0;
    if (v->nSourceTracks() > 0) {
      pdgId = 1;
    } else if ((isSig) && (v->nGenVertices() > 0)) {
      // Signal vertices without a valid parent TrackingParticle come directly from the generator instead. Get their pdgId by:
      // 1. Get from a daughter TrackingParticle a reference GenParticle of the decay
      // 2. Loop over the GenVertices associated to the TrackingVertex (can be many since association is fully distance based)
      //    to find the correct one (which has the reference particle as daughter)
      // 3. Get the pdgId of the mother GenParticle of the GenVertex found at step 2

      auto const daughterTP = v->daughterTracks_begin();
      auto const& daughterSimTrack = (*daughterTP)->g4Tracks()[0];
      HepMC::GenParticle* daughterGenParticle = evt->barcode_to_particle(daughterSimTrack.genpartIndex());
      HepMC::GenVertex* genVertex = daughterGenParticle->production_vertex();
      HepMC::GenVertex::particles_in_const_iterator motherGenParticle = genVertex->particles_in_const_begin();
      pdgId = (*motherGenParticle)->pdg_id();

      // if (daughterGenParticle != nullptr) {
      //   for (auto const& genVertex : v->genVertices()) {
      //     for (auto genParticle = genVertex->particles_out_const_begin();
      //          genParticle != genVertex->particles_out_const_end();
      //          ++genParticle) {
      //       if ((**genParticle) == (*daughterGenParticle)) {
      //         auto const motherParticle = genVertex->particles_in_const_begin();
      //         pdgId = (*motherParticle)->pdg_id();
      //         break;
      //       }
      //     }
      //     if (pdgId == -313)
      //       printFullTree(genVertex.get());
      //     if (pdgId != 4)
      //       break;
      //   }
      // }
      // double minDst = std::numeric_limits<double>::max();
      // for (auto const& genVertex : v->genVertices()) {
      //   //   pdgId = -5;
      //   //   if (genVertex->particles_in_size() > 0) {
      //   //     //     unsigned int parentId = genVertex.parentIndex();
      //   //     //     for (const auto& genParticle : *genParticles) {
      //   //     //       if (genParticle.genParticleIndex() == parentId) {
      //   //     //         pdgId = genParticle.pdgId();
      //   //     //         break;
      //   //     //       }
      //   //     pdgId = -6;
      //   //     for (auto genparticle = genVertex->particles_in_const_begin();
      //   //          genparticle != genVertex->particles_in_const_end();
      //   //          ++genparticle) {
      //   //       pdgId = (*genparticle)->pdg_id();
      //   //       if (pdgId == 111) {
      //   //         std::cout << "GenVertex for pdgId " << pdgId << " (" << genVertex->position().x() << ", "
      //   //                   << genVertex->position().y() << ", " << genVertex->position().z() << ")" << std::endl;
      //   //         std::cout << "TrackingVertex (" << v->position().x() << ", " << v->position().y() << ", "
      //   //                   << v->position().z() << ")" << std::endl;
      //   //         printFullTree(genVertex.get());
      //   //         break;
      //   //       }
      //   //       if (pdgId == 421)
      //   //         break;
      //   //     }
      //   //   }
      //   //   if (pdgId == 421)
      //   //     break;
      //   // loop and find the closest distance
      //   double dst = std::sqrt(std::pow(0.1 * genVertex->position().x() - v->position().x(), 2) +
      //                          std::pow(0.1 * genVertex->position().y() - v->position().y(), 2) +
      //                          std::pow(0.1 * genVertex->position().z() - v->position().z(), 2));
      //   if (dst < minDst) {
      //     minDst = dst;
      //   }
      // }
      // for (auto const& genVertex : v->genVertices()) {
      //   double dst = std::sqrt(std::pow(0.1 * genVertex->position().x() - v->position().x(), 2) +
      //                          std::pow(0.1 * genVertex->position().y() - v->position().y(), 2) +
      //                          std::pow(0.1 * genVertex->position().z() - v->position().z(), 2));
      //   if (dst == minDst) {
      //     nGVs++;
      //   }
      // }
    }

    TrackingVertex_pdgId.push_back(pdgId);
    TrackingVertex_nNearestGenVertices.push_back(nGVs);
    TrackingVertex_nTrackingParticles.push_back(v->nDaughterTracks());

    for (auto daughterTP = v->daughterTracks_begin(); daughterTP != v->daughterTracks_end(); ++daughterTP) {
      auto momentum = (*(*daughterTP)).p4();
      TrackingParticles_pt.push_back(momentum.Pt());
      TrackingParticles_eta.push_back(momentum.Eta());
      TrackingParticles_phi.push_back(momentum.Phi());
      TrackingParticles_pdgId.push_back((*(*daughterTP)).pdgId());
    }
  }  // End of for loop on tracking vertices

  //  Build FlatTables

  // auto genPVTable = std::make_unique<nanoaod::FlatTable>(1, "GenPV", true);
  // genPVTable->addColumn<float>("x", xGenPV_vec, "Gen PV x coordinate");
  // genPVTable->addColumn<float>("y", yGenPV_vec, "Gen PV y coordinate");
  // genPVTable->addColumn<float>("z", zGenPV_vec, "Gen PV z coordinate");

  auto gvTable = std::make_unique<nanoaod::FlatTable>(ngv, "TV", false);
  // gvTable->addColumn<float>("pt", TrackingVertex_pt, "Hadron pt");
  // gvTable->addColumn<float>("eta", TrackingVertex_eta, "Hadron eta");
  // gvTable->addColumn<float>("phi", TrackingVertex_phi, "Hadron phi");
  // gvTable->addColumn<float>("dlen", TrackingVertex_dlen, "GV decay length");
  // gvTable->addColumn<float>("dxy", TrackingVertex_dxy, "GV 2D decay length in xy");
  gvTable->addColumn<float>("x", TrackingVertex_x, "TV x");
  gvTable->addColumn<float>("y", TrackingVertex_y, "TV y");
  gvTable->addColumn<float>("z", TrackingVertex_z, "TV z");
  // gvTable->addColumn<float>("x_i", TrackingVertex_x_i, "Born x coordinate of GV ");
  // gvTable->addColumn<float>("y_i", TrackingVertex_y_i, "Born y coordinate of GV ");
  // gvTable->addColumn<float>("z_i", TrackingVertex_z_i, "Born z coordinate of GV ");
  // gvTable->addColumn<int>("SVIdx", TrackingVertex_SVIdx, "SVIdx");
  gvTable->addColumn<int>("pdgId", TrackingVertex_pdgId, "TrackingVertex pdgId");
  gvTable->addColumn<int>("nDaughters", TrackingVertex_nDaughters, "Number of daughters");
  gvTable->addColumn<int>(
      "nGoodDaughters", TrackingVertex_nGoodDaughters, "Number of good daughters (with reco track)");
  gvTable->addColumn<int>("nSimVertices", TrackingVertex_nSimVertices, "Number of associated sim vertices");
  gvTable->addColumn<int>("nGenVertices", TrackingVertex_nGenVertices, "Number of associated gen vertices");
  gvTable->addColumn<int>(
      "nNearestGenVertices", TrackingVertex_nNearestGenVertices, "Number of nearest associated gen vertices");
  gvTable->addColumn<int>(
      "nTrackingParticles", TrackingVertex_nTrackingParticles, "Number of associated tracking particles");
  gvTable->addColumn<int>("isSig", TrackingVertex_isSig, "Is signal vertex");
  // new class
  // gvTable->addColumn<int>("isB", TrackingVertex_isB, "isB");
  // gvTable->addColumn<int>("isD", TrackingVertex_isD, "isD");
  // gvTable->addColumn<int>("pdgClass", TrackingVertex_pdgId, "pdgClass");

  auto tpTable = std::make_unique<nanoaod::FlatTable>(TrackingParticles_pt.size(), "TP", false);
  tpTable->addColumn<float>("pt", TrackingParticles_pt, "TP pt");
  tpTable->addColumn<float>("eta", TrackingParticles_eta, "TP eta");
  tpTable->addColumn<float>("phi", TrackingParticles_phi, "TP phi");
  tpTable->addColumn<int>("pdgId", TrackingParticles_pdgId, "TP pdgId");

  // auto dauTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(), "GVDaughters", false);
  // dauTable->addColumn<float>("pt", Daughters_pt, "Daughter pt");
  // dauTable->addColumn<float>("eta", Daughters_eta, "Daughter eta");
  // dauTable->addColumn<float>("phi", Daughters_phi, "Daughter phi");
  // dauTable->addColumn<int>("charge", Daughters_charge, "Daughter charge");
  // dauTable->addColumn<int>("hadronIndex", Daughters_GVidx, "Hadron index");

  // auto svdauTable = std::make_unique<nanoaod::FlatTable>(SVtrk_pt.size(), "SVDaughters", false);
  // svdauTable->addColumn<float>("pt", SVtrk_pt, "Daughter pt");
  // svdauTable->addColumn<float>("eta", SVtrk_eta, "Daughter eta");
  // svdauTable->addColumn<float>("phi", SVtrk_phi, "Daughter phi");
  // //svdauTable->addColumn<int>("charge",SVtrk_charge,"Daughter charge");
  // svdauTable->addColumn<int>("SVIdx", SVtrk_SVidx, "Hadron index");

  //
  // iEvent.put(std::move(genPVTable), "GenPV");
  iEvent.put(std::move(gvTable), "TV");
  iEvent.put(std::move(tpTable), "TP");
  // iEvent.put(std::move(dauTable), "GVDaughters");
  // iEvent.put(std::move(svdauTable), "SVDaughters");
}

DEFINE_FWK_MODULE(TrackingVertexTableProducer);
