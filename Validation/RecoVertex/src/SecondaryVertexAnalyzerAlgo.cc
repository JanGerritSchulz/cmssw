#include "Validation/RecoVertex/interface/SecondaryVertexAnalyzerAlgo.h"

#include <algorithm>
#include <cmath>
#include <set>

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMBookingHelpers.h"
using namespace dqm::booking;

// =============================================================================
// Constructor
// =============================================================================

SecondaryVertexAnalyzerAlgo::SecondaryVertexAnalyzerAlgo(const Config &cfg) : cfg_(cfg) {}

// =============================================================================
// isReconstructable
// =============================================================================

bool SecondaryVertexAnalyzerAlgo::isReconstructable(const SimSecondaryVertex &sv, uint32_t skipCuts) const {
  if (!(skipCuts & kDecayLength) && sv.decayLength < cfg_.minDecayLength)
    return false;
  if (!(skipCuts & kNDaughters) && sv.nCharged < cfg_.minReconstructableDaughters)
    return false;
  if (!(skipCuts & kEta) && std::abs(std::atanh(sv.z / std::hypot(sv.r, sv.z))) > cfg_.absEtaMax)
    return false;
  if (!(skipCuts & kPdgId) && !cfg_.signalPdgIds.empty()) {
    const int absPdg = std::abs(sv.motherPdgId);
    const bool found = std::find(cfg_.signalPdgIds.begin(), cfg_.signalPdgIds.end(), absPdg) != cfg_.signalPdgIds.end();
    if (!found)
      return false;
  }
  return true;
}

// =============================================================================
// bookHistograms
// =============================================================================

void SecondaryVertexAnalyzerAlgo::bookHistograms(IBooker &ibook, const std::vector<std::string> &collectionLabels) {
  // Generic sim plots — booked once, collection-independent
  if (cfg_.doGenericSimPlots) {
    ibook.setCurrentFolder(cfg_.rootFolder);
    genericSimHistos_.h_decayLength =
        book1DLogX(ibook, "GenSV_decayLength", "All sim SVs;3D decay length L_{3D} [cm];Entries", 50, 1e-3, 100.);
    genericSimHistos_.h_r =
        book1DLogX(ibook, "GenSV_r", "All sim SVs;Transverse decay radius r_{T} [cm];Entries", 50, 1e-3, 50.);
    genericSimHistos_.h_nDaughters =
        ibook.book1D("GenSV_nDaughters", "All sim SVs;N charged daughters;Entries", 20, -0.5, 19.5);
    genericSimHistos_.h_motherPdgId =
        ibook.book1D("GenSV_motherPdgId", "All sim SVs;|mother PDG ID|;Entries", 600, 0., 6000.);
    genericSimHistos_.h_numAllSimSVs =
        ibook.book1D("GenSV_numAll", "N sim SVs per event (all);N sim SVs;Entries", 100, 0., 200.);
    genericSimHistos_.h_numSignalSimSVs =
        ibook.book1D("GenSV_numSignal", "N sim SVs per event (signal selected);N sim SVs;Entries", 100, 0., 200.);
  }

  // Per-collection histograms
  for (const auto &label : collectionLabels) {
    const std::string folder = cfg_.rootFolder + "/" + label;
    ibook.setCurrentFolder(folder);

    auto &ch = collectionHistos_[label];
    auto &me = mes_[label];

    // ----- Monitoring bundles -----
    // Decay length — suppress kDecayLength cut for sim-side fills
    ch.h_decayLength.bundle.book1DLogX(
        ibook, true, true, cfg_.doPerPdgPlots, "decayLength", "3D decay length L_{3D} [cm]", "Entries", 50, 1e-3, 100.);

    // Decay length significance
    ch.h_decayLengthSig.bundle.book1DLogX(
        ibook, true, true, cfg_.doPerPdgPlots, "decayLengthSig", "L_{3D}/#sigma_{L_{3D}}", "Entries", 50, 0.1, 1000.);

    // Transverse radius
    ch.h_r.bundle.book1DLogX(
        ibook, true, true, cfg_.doPerPdgPlots, "r", "Transverse decay radius r_{T} [cm]", "Entries", 50, 1e-3, 50.);

    // Track multiplicity — suppress kNDaughters cut
    ch.h_nTracks.bundle.book1D(
        ibook, true, true, cfg_.doPerPdgPlots, "nTracks", "N tracks at SV", "Entries", 20, -0.5, 19.5);

    // Eta — suppress kEta cut
    ch.h_eta.bundle.book1D(ibook,
                           true,
                           true,
                           cfg_.doPerPdgPlots,
                           "eta",
                           "#eta of SV position",
                           "Entries",
                           50,
                           -cfg_.absEtaMax * 1.2,
                           cfg_.absEtaMax * 1.2);

    // Normalised chi2
    ch.h_chi2ndof.bundle.book1D(ibook, true, true, false, "chi2ndof", "Normalised #chi^{2}", "Entries", 50, 0., 10.);

    // Invariant mass — CPC only; booked for all but only filled when available
    ch.h_mass.bundle.book1D(
        ibook, true, true, cfg_.doPerPdgPlots, "mass", "SV invariant mass [GeV]", "Entries", 100, 0., 10.);

    // ----- Resolution bundles -----
    // Bin axes: decay length [0,30 cm], r [0,10 cm], nTracks [0,20]
    ch.h_xRes.bookResolutions(ibook, 30, 0., 30., 20, 0., 10., 20, 0., 20., "x_", 100, -0.05, 0.05);
    ch.h_yRes.bookResolutions(ibook, 30, 0., 30., 20, 0., 10., 20, 0., 20., "y_", 100, -0.05, 0.05);
    ch.h_zRes.bookResolutions(ibook, 30, 0., 30., 20, 0., 10., 20, 0., 20., "z_", 100, -0.05, 0.05);
    ch.h_decayLengthRes.bookResolutions(ibook,
                                        20,
                                        0.,
                                        20.,
                                        50,
                                        -3.,
                                        3.,  // nTracks, eta, lRes
                                        50,
                                        -0.5,
                                        0.5,  // lRes axis
                                        50,
                                        -30.,
                                        30.);  // lSigRes axis
    ch.h_massRes.bookResolutions(ibook, 30, 0., 30., 20, 0., 10., 20, 0., 20., "mass_", 100, -1., 1.);

    // ----- Plain per-collection histograms -----
    me["numRecoSVs"] = ibook.book1D("numRecoSVs", "N reco SVs per event;N SVs;Entries", 100, 0., 200.);
    me["numSimSVsAll"] = ibook.book1D("numSimSVsAll", "N all sim SVs per event;N SVs;Entries", 100, 0., 200.);
    me["numSimSVsSignal"] = ibook.book1D("numSimSVsSignal", "N signal sim SVs per event;N SVs;Entries", 100, 0., 200.);
    me["sharedTrackFraction"] = ibook.book1D("sharedTrackFraction",
                                             "Shared track fraction (matched pairs);"
                                             "Shared fraction;Entries",
                                             50,
                                             0.,
                                             1.);
  }
}

// =============================================================================
// Sim vertex building
// =============================================================================

void SecondaryVertexAnalyzerAlgo::buildSimSVs(const TrackingVertexCollection &simVertices) {
  allSimSVs_ = buildAllSimSVs(simVertices);
  signalSimSVs_ = buildSignalSimSVs();

  if (cfg_.verbose) {
    LogDebug("SecondaryVertexAnalyzer") << "SimSecondaryVertex overview: " << allSimSVs_.size() << " all sim SVs, "
                                        << signalSimSVs_.size() << " signal sim SVs";
  }
}

void SecondaryVertexAnalyzerAlgo::resetSimSVs() {
  for (auto &sv : allSimSVs_) {
    sv.num_matched_reco_vertices = 0;
    sv.average_match_quality = 0.0;
    sv.matched_reco_shared_fractions.clear();
  }
}

void SecondaryVertexAnalyzerAlgo::clearSimSVs() {
  allSimSVs_.clear();
  signalSimSVs_.clear();
}

// TODO: fix this according to the logic I developed in the other package
int SecondaryVertexAnalyzerAlgo::motherPdgId(const TrackingVertex &tv) const {
  // Walk the first TrackingParticle's parent chain looking for a B or D hadron.
  // Fall back to the direct parent PDG ID if none found.
  if (tv.nSourceTracks() == 0)
    return 0;

  const TrackingParticleRef &tp = *tv.sourceTracks_begin();

  // Walk up through parent vertices
  int lastPdgId = 0;
  TrackingParticleRef current = tp;
  while (current.isNonnull()) {
    // Check parent vertices of this TP
    if (current->parentVertex().isNull())
      break;
    const TrackingVertex &parentVtx = *(current->parentVertex());
    if (parentVtx.nSourceTracks() == 0)
      break;

    const TrackingParticleRef &parent = *parentVtx.sourceTracks_begin();
    const int absPdg = std::abs(parent->pdgId());
    lastPdgId = parent->pdgId();

    // Stop at B hadrons (5xx, 5xxx) or D hadrons (4xx, 4xxx)
    const bool isB = (absPdg / 500 == 1) || (absPdg / 5000 == 1);
    const bool isC = (absPdg / 400 == 1) || (absPdg / 4000 == 1);
    if (isB || isC)
      return parent->pdgId();

    current = parent;
  }
  return lastPdgId;
}

double SecondaryVertexAnalyzerAlgo::decayLength(const TrackingVertex &tv, const TrackingVertex &pv) const {
  const auto &svPos = tv.position();
  const auto &pvPos = pv.position();
  const double dx = svPos.x() - pvPos.x();
  const double dy = svPos.y() - pvPos.y();
  const double dz = svPos.z() - pvPos.z();
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<SimSecondaryVertex> SecondaryVertexAnalyzerAlgo::buildAllSimSVs(
    const TrackingVertexCollection &simVertices) const {
  std::vector<SimSecondaryVertex> result;
  result.reserve(simVertices.size());

  // get the PV for the decay length calculation
  const TrackingVertex &pv = simVertices[0];

  // skip the PV in the SimSV filling
  for (size_t i = 1; i < simVertices.size(); ++i) {
    const TrackingVertex &tv = simVertices[i];

    // Skip the hard-scatter primary vertex (event 0, BX 0, first vertex)
    // and out-of-time pileup.
    if (tv.eventId().bunchCrossing() != 0)
      continue;
    // Skip vertices with no outgoing daughter tracks — these are not SVs
    if (tv.nDaughterTracks() == 0)
      continue;

    const auto &pos = tv.position();
    SimSecondaryVertex sv(pos.x(), pos.y(), pos.z());

    sv.decayLength = decayLength(tv, pv);
    sv.motherPdgId = motherPdgId(tv);
    sv.isFromPileup = (tv.eventId().event() != 0);
    sv.eventId = tv.eventId();

    // Count charged daughters
    sv.nCharged = 0;
    sv.nReconstructable = 0;
    for (auto iTP = tv.daughterTracks_begin(); iTP != tv.daughterTracks_end(); ++iTP) {
      if ((*iTP)->charge() != 0) {
        ++sv.nCharged;
        // A daughter is reconstructable if it has enough hits — use the
        // standard threshold of >=3 hits as a proxy; this can be made
        // configurable if needed.
        if ((*iTP)->numberOfTrackerHits() >= 3)
          ++sv.nReconstructable;
      }
    }

    // Store ref for later association lookup
    sv.simVertex = TrackingVertexRef(edm::Ref<TrackingVertexCollection>(&simVertices, i));

    result.push_back(std::move(sv));
  }
  return result;
}

std::vector<SimSecondaryVertex *> SecondaryVertexAnalyzerAlgo::buildSignalSimSVs() {
  std::vector<SimSecondaryVertex *> result;
  result.reserve(allSimSVs_.size());
  for (auto &sv : allSimSVs_) {
    // Since vertices from the signal interaction come first, break after the first PU vertex
    if (sv.isFromPileup)
      break;
    // Apply all cuts (no suppression — this is the global signal selection,
    // not the per-bundle variable-blind filling path).
    if (isReconstructable(sv, kNone))
      result.push_back(&sv);
  }
  return result;
}

// =============================================================================
// Reco vertex building
// =============================================================================

std::vector<RecoSecondaryVertex> SecondaryVertexAnalyzerAlgo::buildRecoSVs(
    const edm::View<reco::Vertex> &recoVertices) const {
  std::vector<RecoSecondaryVertex> result;
  result.reserve(recoVertices.size());
  for (size_t i = 0; i < recoVertices.size(); ++i) {
    const auto &vtx = recoVertices[i];
    if (vtx.isFake() || !vtx.isValid() || vtx.ndof() < 0.)
      continue;
    RecoSecondaryVertex rv(vtx.x(), vtx.y(), vtx.z());
    rv.chi2 = vtx.chi2();
    rv.ndof = vtx.ndof();
    rv.nTracks = static_cast<int>(vtx.tracksSize());
    rv.recoVertexRef = recoVertices.refAt(i);
    result.push_back(std::move(rv));
  }
  return result;
}

std::vector<RecoSecondaryVertex> SecondaryVertexAnalyzerAlgo::buildRecoSVs(
    const edm::View<reco::VertexCompositePtrCandidate> &recoVertices) const {
  std::vector<RecoSecondaryVertex> result;
  result.reserve(recoVertices.size());
  for (size_t i = 0; i < recoVertices.size(); ++i) {
    const auto &vtx = recoVertices[i];
    if (vtx.numberOfDaughters() == 0)
      continue;
    RecoSecondaryVertex rv(vtx.vx(), vtx.vy(), vtx.vz());
    rv.chi2 = vtx.vertexChi2();
    rv.ndof = vtx.vertexNdof();
    rv.nTracks = static_cast<int>(vtx.numberOfDaughters());
    rv.mass = vtx.mass();
    rv.recoVertexCPCRef = recoVertices.refAt(i);
    result.push_back(std::move(rv));
  }
  return result;
}

// =============================================================================
// Association and matching
// =============================================================================

template <typename AssociatorType>
void SecondaryVertexAnalyzerAlgo::matchSim2RecoVertices(const AssociatorType &simToReco) {
  for (auto &sv : allSimSVs_) {
    if (sv.simVertex.isNull())
      continue;
    auto it = simToReco.find(sv.simVertex);
    if (it == simToReco.end())
      continue;
    sv.num_matched_reco_vertices = static_cast<int>(it->val.size());
    sv.matched_reco_shared_fractions.clear();
    float qualSum = 0.f;
    for (const auto &recoAndQuality : it->val) {
      sv.matched_reco_shared_fractions.push_back(recoAndQuality.second);
      qualSum += recoAndQuality.second;
    }
    if (sv.num_matched_reco_vertices > 0)
      sv.average_match_quality = qualSum / static_cast<float>(sv.num_matched_reco_vertices);
  }
}

template <typename AssociatorType>
void SecondaryVertexAnalyzerAlgo::matchReco2SimVertices(std::vector<RecoSecondaryVertex> &recoSVs,
                                                        const AssociatorType &recoToSim) const {
  // key_type is the edm::RefToBase, and value_type either reco::Vertex or reco::VertexCompositePtrCandidate
  using VertexType = AssociatorType::key_type::value_type;

  // Build a lookup set of signal sim SV pointers for O(1) membership test.
  std::set<const SimSecondaryVertex *> signalSet(signalSimSVs_.begin(), signalSimSVs_.end());

  // Create lookup set for keeping track of reconstructed SimSVs and identify duplicates
  std::set<const SimSecondaryVertex *> reconstructedSimSVs;

  // Build a map from TrackingVertexRef key to SimSecondaryVertex* for
  // reverse lookup from the association map results.
  std::map<unsigned int, const SimSecondaryVertex *> keyToSimSV;
  for (const auto &sv : allSimSVs_) {
    if (sv.simVertex.isNonnull())
      keyToSimSV[sv.simVertex.key()] = &sv;
  }

  for (auto &rv : recoSVs) {
    // The association map is keyed by VertexBaseRef / VertexCompositePtrCandidateRef.
    // We stored the raw reco vertex pointer in rv.recVtxPtr; the association
    // map ref must be looked up by the plugin using the original handle index.
    // Here we iterate all map entries and match by pointer identity.
    // TODO: use the SV index in the collection for the association map
    auto it = recoToSim.find(rv.recoVertex<VertexType>());

    if (it == recoToSim.end()) {
      // No sim match found → fake
      rv.kind_of_vertex |= RecoSecondaryVertex::FAKE;
      continue;
    }

    rv.kind_of_vertex |= RecoSecondaryVertex::MATCHED;
    rv.num_matched_sim_vertices = static_cast<int>(it->val.size());

    bool anySignal = false;
    bool anyPileup = false;

    for (const auto &simAndQuality : it->val) {
      const unsigned int key = simAndQuality.first.key();
      auto simIt = keyToSimSV.find(key);
      if (simIt == keyToSimSV.end())
        continue;
      const SimSecondaryVertex *svPtr = simIt->second;
      rv.sim_vertices.push_back(svPtr);
      rv.sim_vertices_shared_fraction.push_back(simAndQuality.second);

      if (signalSet.count(svPtr))
        anySignal = true;
      if (svPtr->isFromPileup)
        anyPileup = true;
    }

    if (rv.sim_vertices.size() > 1)
      rv.kind_of_vertex |= RecoSecondaryVertex::MERGED;

    // A reco SV is flagged as pileup if ALL its matched sim SVs are pileup.
    // If at least one match is signal, it is not flagged as pileup.
    if (anyPileup && !anySignal)
      rv.isFromPileup = true;

    // Propagate mother PDG ID from the best-quality sim match.
    if (!rv.sim_vertices.empty()) {
      const SimSecondaryVertex *bestSim = rv.sim_vertices.front();
      rv.motherPdgId = bestSim->motherPdgId;

      // Duplicate check
      auto [it2, inserted] = reconstructedSimSVs.emplace(bestSim);
      if (!inserted) {
        // Another reco SV already claimed this sim SV — this one is a duplicate.
        rv.kind_of_vertex |= RecoSecondaryVertex::DUPLICATE;
      }
    }
  }
}

// =============================================================================
// Histogram filling
// =============================================================================

void SecondaryVertexAnalyzerAlgo::fillSimVertexHistograms(const std::string &label, const SimSecondaryVertex &sv) {
  const bool isMatched = (sv.num_matched_reco_vertices > 0);
  auto &ch = collectionHistos_.at(label);

  // Generic sim plots (collection-independent, filled once per all-sim pass)
  if (cfg_.doGenericSimPlots && genericSimHistos_.h_decayLength) {
    genericSimHistos_.h_decayLength->Fill(sv.decayLength);
    genericSimHistos_.h_r->Fill(sv.r);
    genericSimHistos_.h_nDaughters->Fill(sv.nCharged);
    genericSimHistos_.h_motherPdgId->Fill(std::abs(sv.motherPdgId));
  }

  // Helper lambda: fill a BundleWithCutMask, evaluating reconstructability
  // with the bundle's own skipCuts mask.
  auto fillBundle = [&](BundleWithCutMask &bwm, const double value) {
    const bool isReco = isReconstructable(sv, bwm.skipCuts);
    bwm.bundle.fillSimVertexHistos(isMatched, isReco, value);
    if (cfg_.doPerPdgPlots)
      bwm.bundle.fillSimVertexHistosByPdg(sv.motherPdgId, isMatched, value);
  };

  fillBundle(ch.h_decayLength, sv.decayLength);
  fillBundle(ch.h_r, sv.r);
  fillBundle(ch.h_nTracks, sv.nCharged);

  // Eta of SV: use pseudorapidity of position vector from origin
  const double eta = (sv.r > 0. || sv.z != 0.) ? std::atanh(sv.z / std::hypot(sv.r, sv.z)) : 0.;
  fillBundle(ch.h_eta, eta);
}

void SecondaryVertexAnalyzerAlgo::fillRecoVertexHistograms(const std::string &label, const RecoSecondaryVertex &rv) {
  const bool isMatched = (rv.kind_of_vertex & RecoSecondaryVertex::MATCHED) != 0;
  const bool isDuplicate = (rv.kind_of_vertex & RecoSecondaryVertex::DUPLICATE) != 0;
  const bool isFake = (rv.kind_of_vertex & RecoSecondaryVertex::FAKE) != 0;
  const bool isMerged = (rv.kind_of_vertex & RecoSecondaryVertex::MERGED) != 0;
  const bool isPileup = rv.isFromPileup;

  auto &ch = collectionHistos_.at(label);

  auto fillBundle = [&](BundleWithCutMask &bwm, const double value) {
    bwm.bundle.fillRecoVertexHistos(isMatched, isDuplicate, isFake, isMerged, isPileup, value);
  };

  fillBundle(ch.h_decayLength, rv.decayLength);
  fillBundle(ch.h_decayLengthSig, rv.decayLengthSignificance);
  fillBundle(ch.h_r, rv.r);
  fillBundle(ch.h_nTracks, rv.nTracks);

  const double eta = (rv.r > 0. || rv.z != 0.) ? std::atanh(rv.z / std::hypot(rv.r, rv.z)) : 0.;
  fillBundle(ch.h_eta, eta);
  fillBundle(ch.h_chi2ndof, rv.normalizedChi2());

  if (rv.mass.has_value())
    fillBundle(ch.h_mass, *rv.mass);
}

void SecondaryVertexAnalyzerAlgo::fillResolutionHistograms(const std::string &label,
                                                           const RecoSecondaryVertex &rv,
                                                           const SimSecondaryVertex &sv) {
  auto &ch = collectionHistos_.at(label);

  const double decayLen = sv.decayLength;
  const double r = sv.r;
  const double nTrk = static_cast<double>(sv.nCharged);
  const double eta = (sv.r > 0. || sv.z != 0.) ? std::atanh(sv.z / std::hypot(sv.r, sv.z)) : 0.;

  // Position residuals (reco - sim) and pulls
  // Pulls require vertex position errors — use covariance diagonal if available.
  // For now fill residuals; pulls are set to residual / 1 as placeholder
  // until covariance is wired through GeneralRecoVertex.
  // TODO: wire position errors from reco vertex covariance matrix.
  const double xRes = rv.x - sv.x;
  const double yRes = rv.y - sv.y;
  const double zRes = rv.z - sv.z;
  const double xPull = xRes;  // placeholder
  const double yPull = yRes;
  const double zPull = zRes;

  ch.h_xRes.fill(decayLen, r, nTrk, xRes, xPull);
  ch.h_yRes.fill(decayLen, r, nTrk, yRes, yPull);
  ch.h_zRes.fill(decayLen, r, nTrk, zRes, zPull);

  // Decay length residual and significance residual
  const double lRes = rv.decayLength - sv.decayLength;
  const double lPull = lRes;        // placeholder
  const double lSigRes = 0.;        //rv.decayLengthSignificance - sv.decayLengthSignificance;
  const double lSigPull = lSigRes;  // placeholder
  ch.h_decayLengthRes.fill(
      nTrk,
      eta,
      lRes,
      lPull,
      lSigRes,
      lSigPull);  // TODO: fix the decayLength residual histogram bundle, makes no sense to me especially the dlenSig...

  // Mass residual — only for CPC vertices
  if (rv.mass.has_value() && sv.decayLength > 0.) {
    // No sim-level mass available from TrackingVertex directly; the mass
    // residual requires knowing the true particle mass from the mother PDG ID.
    // TODO: look up PDG mass from motherPdgId and fill mass resolution.
  }
}

// =============================================================================
// analyzeImpl — shared logic for both vertex types
// =============================================================================

template <typename AssociatorType>
void SecondaryVertexAnalyzerAlgo::analyzeImpl(std::vector<RecoSecondaryVertex> recoSVs,
                                              const TrackingVertexCollection &simVertices,
                                              const AssociatorType &associator,
                                              const reco::RecoToSimCollection & /*trackRecoToSim*/,
                                              const reco::SimToRecoCollection & /*trackSimToReco*/,
                                              const std::string &collectionLabel) {
  // ------------------------------------------------------------------
  // 1. Run association (produces sim↔reco maps for this collection)
  // ------------------------------------------------------------------
  // Note: the associator here is the VertexToTrackingVertexAssociator
  // wrapper fetched from the event by the plugin. We need edm::Handle
  // views for the call — these are provided by the plugin; for now
  // the association maps produced upstream are passed directly.
  // TODO: thread the edm::Handle<edm::View<VertexType>> through to
  // call associator.associateRecoToSim / associateSimToReco directly
  // here, rather than relying on pre-produced maps.
  // For the current draft the plugin is expected to have already run
  // the association and stored results in the event; the vertex-level
  // maps are fetched below via the associator wrapper.
  //
  resetSimSVs();
  // Placeholder: direct calls will be wired in when the plugin passes
  // handles through analyzeImpl.
  (void)associator;
  // TODO: make sure the sorting according to the matching quality is performed

  // ------------------------------------------------------------------
  // 2. Build reco↔sim matching using pre-produced association maps
  // ------------------------------------------------------------------
  // matchSim2RecoVertices and matchReco2SimVertices are called after the
  // association maps are available. In the current architecture, the maps
  // are produced upstream and consumed by the plugin, which passes the
  // already-fetched map handles. The matching functions below receive the
  // collection-level SimToReco / RecoToSim maps.
  //
  // For now: construct empty maps as placeholders.
  // TODO: accept VertexSimToRecoCollection and VertexRecoToSimCollection
  // as arguments to analyzeImpl once the plugin wires them through.
  reco::VertexSimToRecoCollection simToReco;
  reco::VertexRecoToSimCollection recoToSim;

  matchSim2RecoVertices(simToReco);
  matchReco2SimVertices(recoSVs, recoToSim);

  // ------------------------------------------------------------------
  // 3. Fill generic sim histograms (collection-independent, first call only)
  // ------------------------------------------------------------------
  if (cfg_.doGenericSimPlots && genericSimHistos_.h_numAllSimSVs) {
    genericSimHistos_.h_numAllSimSVs->Fill(allSimSVs_.size());
    genericSimHistos_.h_numSignalSimSVs->Fill(signalSimSVs_.size());
  }

  auto &me = mes_.at(collectionLabel);
  me.at("numRecoSVs")->Fill(recoSVs.size());
  me.at("numSimSVsAll")->Fill(allSimSVs_.size());
  me.at("numSimSVsSignal")->Fill(signalSimSVs_.size());

  // ------------------------------------------------------------------
  // 4. Fill sim-side histograms
  //    Signal sim SVs → efficiency numerator/denominator.
  //    All sim SVs → generic sim plots already filled above.
  // ------------------------------------------------------------------
  for (const SimSecondaryVertex *sv : signalSimSVs_)
    fillSimVertexHistograms(collectionLabel, *sv);

  // ------------------------------------------------------------------
  // 5. Fill reco-side histograms and resolution plots
  // ------------------------------------------------------------------
  for (const auto &rv : recoSVs) {
    fillRecoVertexHistograms(collectionLabel, rv);

    // Fill shared track fraction histogram for matched vertices
    if (!rv.sim_vertices_shared_fraction.empty())
      me.at("sharedTrackFraction")->Fill(rv.sim_vertices_shared_fraction.front());

    // Resolution plots: use best-quality sim match (first in sorted list)
    if (!rv.sim_vertices.empty()) {
      fillResolutionHistograms(collectionLabel, rv, *rv.sim_vertices.front());
    }
  }
}

// =============================================================================
// Public analyze() overloads
// =============================================================================

void SecondaryVertexAnalyzerAlgo::analyze(
    const edm::View<reco::Vertex> &recoVertices,
    const TrackingVertexCollection &simVertices,
    const reco::VertexToTrackingVertexAssociator<std::vector<reco::Vertex>> &associator,
    const reco::RecoToSimCollection &trackRecoToSim,
    const reco::SimToRecoCollection &trackSimToReco,
    const std::string &collectionLabel) {
  analyzeImpl(buildRecoSVs(recoVertices), simVertices, associator, trackRecoToSim, trackSimToReco, collectionLabel);
}

void SecondaryVertexAnalyzerAlgo::analyze(
    const edm::View<reco::VertexCompositePtrCandidate> &recoVertices,
    const TrackingVertexCollection &simVertices,
    const reco::VertexToTrackingVertexAssociator<std::vector<reco::VertexCompositePtrCandidate>> &associator,
    const reco::RecoToSimCollection &trackRecoToSim,
    const reco::SimToRecoCollection &trackSimToReco,
    const std::string &collectionLabel) {
  analyzeImpl(buildRecoSVs(recoVertices), simVertices, associator, trackRecoToSim, trackSimToReco, collectionLabel);
}
