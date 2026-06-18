#include "Validation/RecoVertex/interface/SecondaryVertexAnalyzerAlgo.h"

#include <algorithm>
#include <cmath>
#include <set>

#include "DataFormats/Math/interface/deltaR.h"
#include "DQMServices/Core/interface/DQMBookingHelpers.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimTracker/TrackAssociation/interface/trackingVertexMotherPdgId.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

using namespace dqm::booking;

// =============================================================================
// Constructor
// =============================================================================

SecondaryVertexAnalyzerAlgo::SecondaryVertexAnalyzerAlgo(const Config &cfg) : cfg_(cfg) {}

// =============================================================================
// isEligibleForEff
// =============================================================================

bool SecondaryVertexAnalyzerAlgo::isEligibleForEff(const SimSecondaryVertex &sv, EfficiencyEligibility mask) const {
  // check if the required bit in the mask (if any) is set in the vertex
  return (sv.eligibility & mask) == mask;
}

// =============================================================================
// bookHistograms
// =============================================================================

void SecondaryVertexAnalyzerAlgo::bookHistograms(IBooker &ibook, const std::vector<std::string> &collectionLabels) {
  auto bins = SVResolutionBundle::BinConfig();
  // config = {nBins, min, max}
  bins.decayLength = {20, 0., 5.};
  bins.decayLength2D = bins.decayLength;
  bins.eta = {45, -4.5, 4.5};
  bins.pt = {40, 0.1, 1000.};
  bins.nTracks = {21, -0.5, 20.5};

  // Generic sim plots — booked once, collection-independent
  if (cfg_.doGenericSimPlots) {
    ibook.setCurrentFolder(cfg_.rootFolder);
    genericSimHistos_.h_decayLength =
        book1DLogX(ibook, "SimSV_decayLength", "All sim SVs;3D decay length L_{3D} [cm];Entries", 50, 1e-3, 100.);
    genericSimHistos_.h_r =
        book1DLogX(ibook, "SimSV_r", "All sim SVs;Transverse decay radius r_{T} [cm];Entries", 50, 1e-3, 50.);
    genericSimHistos_.h_nDaughters =
        ibook.book1D("SimSV_nDaughters", "All sim SVs;N charged daughters;Entries", 20, -0.5, 19.5);
    genericSimHistos_.h_motherPdgId =
        ibook.book1D("SimSV_motherPdgId", "All sim SVs;|mother PDG ID|;Entries", 600, 0., 6000.);
    genericSimHistos_.h_numAllSimSVs =
        ibook.book1D("SimSV_numAll", "N sim SVs per event (all);N sim SVs;Entries", 100, 0., 200.);
    genericSimHistos_.h_numSignalSimSVs =
        ibook.book1D("SimSV_numSignal", "N sim SVs per event (signal selected);N sim SVs;Entries", 100, 0., 200.);
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

    // pt of summed track 4-mometum vectors
    ch.h_pt.bundle.book1DLogX(ibook, true, true, cfg_.doPerPdgPlots, "pt", "SV p_{T} [GeV]", "Entries", 50, 0.1, 1000.);

    // Invariant mass — CPC only; booked for all but only filled when available
    ch.h_mass.bundle.book1DLogX(
        ibook, true, true, cfg_.doPerPdgPlots, "mass", "SV invariant mass [GeV]", "Entries", 50, 0.1, 1000.);

    // ----- Resolution bundles -----
    // Bin axes: decay length [0,30 cm], r [0,10 cm], nTracks [0,20]
    ch.h_xRes.bookResolutions(ibook, bins, "x", 100, -0.05, 0.05);
    ch.h_yRes.bookResolutions(ibook, bins, "y", 100, -0.05, 0.05);
    ch.h_zRes.bookResolutions(ibook, bins, "z", 100, -0.05, 0.05);
    ch.h_decayLengthRes.bookResolutions(ibook, bins, "decayLength", 50, -3., 3.);
    ch.h_ptRes.bookResolutions(ibook, bins, "pt", 100, -10., 10.);
    ch.h_etaRes.bookResolutions(ibook, bins, "eta", 100, -0.2, 0.2);
    ch.h_massRes.bookResolutions(ibook, bins, "mass", 100, -1., 1.);

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
// Set the primary vertex
// =============================================================================

void SecondaryVertexAnalyzerAlgo::setPrimaryVertex(const edm::Handle<reco::VertexCollection> &pvsHandle) {
  if (pvsHandle.isValid() && !pvsHandle->empty()) {
    pv_ = pvsHandle->front();
    return;
  }

  edm::LogWarning("SecondaryVertexAnalyzer") << "Primary vertex collection missing or empty — "
                                                "falling back to detector center for decay length calculations.";

  // Fallback: detector center, zero uncertainty.
  // reco::Vertex(position, error, chi2, ndof, size) — a minimal fake vertex.
  const reco::Vertex::Point origin(0., 0., 0.);
  const reco::Vertex::Error zeroError;  // default-constructed: all zeros
  pv_ = reco::Vertex(origin, zeroError, 0., 0., 0);
}

// =============================================================================
// Sim vertex building
// =============================================================================

void SecondaryVertexAnalyzerAlgo::prepareEventTruth(const edm::Handle<TrackingVertexCollection> &simVerticesH,
                                                    const HepMC::GenEvent *genEvent = nullptr) {
  allSimSVs_ = buildAllSimSVs(simVerticesH);
  signalSimSVs_ = buildSignalSimSVs(genEvent);

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

void SecondaryVertexAnalyzerAlgo::clearEventTruth() {
  allSimSVs_.clear();
  signalSimSVs_.clear();
}

double SecondaryVertexAnalyzerAlgo::decayLength(const TrackingVertex &tv, const TrackingVertex &pv) const {
  const auto &svPos = tv.position();
  const auto &pvPos = pv.position();
  const double dx = svPos.x() - pvPos.x();
  const double dy = svPos.y() - pvPos.y();
  const double dz = svPos.z() - pvPos.z();
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

EfficiencyPrecheck SecondaryVertexAnalyzerAlgo::precheckEligibility(const SimSecondaryVertex &sv) const {
  const int failsDecayLength = (sv.decayLength < cfg_.minDecayLength) || (sv.decayLength > cfg_.maxDecayLength);
  const int failsNDaughters = sv.nCharged < cfg_.minReconstructableDaughters;
  const int failsEta = std::abs(sv.eta()) > cfg_.absEtaMax;

  EfficiencyPrecheck result;
  result.nFailingCuts = failsDecayLength + failsNDaughters + failsEta;

  if (result.nFailingCuts - failsDecayLength <= 0)  // eligible for eff vs. decay length
    result.eligibility |= EfficiencyEligibility::kDecayLength;
  if (result.nFailingCuts - failsNDaughters <= 0)  // eligible for eff vs. nTracks
    result.eligibility |= EfficiencyEligibility::kNDaughters;
  if (result.nFailingCuts - failsEta <= 0)  // eligible for eff vs. eta
    result.eligibility |= EfficiencyEligibility::kEta;

  return result;
}

bool SecondaryVertexAnalyzerAlgo::finalizeEligibility(SimSecondaryVertex &sv,
                                                      const EfficiencyPrecheck &precheck) const {
  EfficiencyEligibility result = precheck.eligibility;

  bool passPdgIdCut = true;
  const auto pdgId = sv.motherPdgId;
  if (!(cfg_.bHadrons) && sim::isBHadron(pdgId))
    passPdgIdCut = false;
  if (!(cfg_.cHadrons) && sim::isCHadron(pdgId))
    passPdgIdCut = false;
  if (!(cfg_.otherParticles) && !(sim::isBHadron(pdgId) || sim::isCHadron(pdgId)))
    passPdgIdCut = false;
  if (!cfg_.signalPdgIds.empty()) {
    const int absPdg = std::abs(pdgId);
    passPdgIdCut = std::find(cfg_.signalPdgIds.begin(), cfg_.signalPdgIds.end(), absPdg) != cfg_.signalPdgIds.end();
  }

  // kPdgId bundle: eligible if all three cheap cuts pass (PDG cut itself
  // is suppressed for this bundle's own plot).
  if (precheck.nFailingCuts == 0)
    result |= EfficiencyEligibility::kPdgId;

  // The cheap-cut bundles additionally require the PDG cut to pass, since
  // they do NOT suppress it.
  if (!passPdgIdCut) {
    // Clear all cheap-cut bits if the PDG cut fails — those bundles do not
    // suppress the PDG cut, so failing it disqualifies them regardless of
    // the cheap-cut outcome.
    result = static_cast<EfficiencyEligibility>(static_cast<uint32_t>(result) &
                                                ~static_cast<uint32_t>(EfficiencyEligibility::kDecayLength |
                                                                       EfficiencyEligibility::kNDaughters |
                                                                       EfficiencyEligibility::kEta));
  }

  sv.eligibility = result;
  return result != EfficiencyEligibility::kNone;
}

std::vector<SimSecondaryVertex> SecondaryVertexAnalyzerAlgo::buildAllSimSVs(
    const edm::Handle<TrackingVertexCollection> &simVerticesH) const {
  const TrackingVertexCollection &simVertices = *simVerticesH;

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
    sv.motherPdgId = 0;  // mother pdgId is assigned (later) for signal used in efficiency only
    sv.isFromPileup = (tv.eventId().event() != 0);
    sv.eventId = tv.eventId();
    sv.eligibility = EfficiencyEligibility::kNone;

    // Count charged daughters
    sv.nCharged = 0;
    sv.nReconstructable = 0;
    for (auto iTP = tv.daughterTracks_begin(); iTP != tv.daughterTracks_end(); ++iTP) {
      if ((*iTP)->charge() != 0) {
        ++sv.nCharged;
        sv.chargedP4 += (*iTP)->p4();
        // A daughter is reconstructable if it has enough hits — use the
        // standard threshold of >=3 hits as a proxy; this can be made
        // configurable if needed.
        if (((*iTP)->numberOfTrackerHits() >= 3) && (*iTP)->pt() >= cfg_.minPtReconstructableDaughters)
          ++sv.nReconstructable;
      }
    }

    if (sv.nCharged == 0)
      continue;

    // Store ref for later association lookup
    sv.simVertex = TrackingVertexRef(simVerticesH, i);

    result.push_back(std::move(sv));
  }
  return result;
}

std::vector<SimSecondaryVertex *> SecondaryVertexAnalyzerAlgo::buildSignalSimSVs(
    const HepMC::GenEvent *genEvent = nullptr) {
  std::vector<SimSecondaryVertex *> result;
  result.reserve(allSimSVs_.size());
  for (auto &sv : allSimSVs_) {
    // Since vertices from the signal interaction come first, break after the first PU vertex
    if (sv.isFromPileup)
      break;

    // veto vertices with less than 2 charged daugters
    if (sv.nCharged < 2)
      continue;

    // Apply all cuts for checking the eligibility for efficiency calculation.
    // And determine the PDG ID for vertices worth checking (expensive HepMC tree search).
    auto preCheck = precheckEligibility(sv);

    if (!preCheck.potentiallyEligible())
      continue;

    sv.motherPdgId = sim::trackingVertexMotherPdgId(*(sv.simVertex), genEvent);

    if (finalizeEligibility(sv, preCheck))
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
  VertexDistance3D vdist;
  for (size_t i = 0; i < recoVertices.size(); ++i) {
    const auto &vtx = recoVertices[i];
    if (vtx.isFake() || !vtx.isValid() || vtx.ndof() < 0.)
      continue;
    RecoSecondaryVertex rv(vtx.x(), vtx.y(), vtx.z());
    Measurement1D dl =
        vdist.distance(pv_, VertexState(RecoVertex::convertPos(vtx.position()), RecoVertex::convertError(vtx.error())));
    rv.decayLength = dl.value();
    rv.decayLengthSignificance = dl.significance();
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
  VertexDistance3D vdist;
  for (size_t i = 0; i < recoVertices.size(); ++i) {
    const auto &vtx = recoVertices[i];
    if (vtx.numberOfDaughters() == 0)
      continue;
    RecoSecondaryVertex rv(vtx.vx(), vtx.vy(), vtx.vz());
    Measurement1D dl =
        vdist.distance(pv_, VertexState(RecoVertex::convertPos(vtx.position()), RecoVertex::convertError(vtx.error())));
    rv.decayLength = dl.value();
    rv.decayLengthSignificance = dl.significance();
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
  for (auto &svp : signalSimSVs_) {
    auto &sv = *svp;
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
  const bool isReconstructable = true;  // TODO: implement a proper reconstructable test
  auto &ch = collectionHistos_.at(label);

  // Generic sim plots (collection-independent, filled once per all-sim pass)
  if (cfg_.doGenericSimPlots && genericSimHistos_.h_decayLength) {
    genericSimHistos_.h_decayLength->Fill(sv.decayLength);
    genericSimHistos_.h_r->Fill(sv.r);
    genericSimHistos_.h_nDaughters->Fill(sv.nCharged);
    genericSimHistos_.h_motherPdgId->Fill(std::abs(sv.motherPdgId));
  }

  // Helper lambda: fill a BundleWithCutMask, evaluating reconstructability
  // with the bundle's own mask.
  auto fillBundle = [&](BundleWithCutMask &bwm, const double value) {
    if (!isEligibleForEff(sv, bwm.mask))
      return;  // not eligible for this set of plots
    bwm.bundle.fillSimVertexHistos(isMatched, isReconstructable, value);
    if (cfg_.doPerPdgPlots)
      bwm.bundle.fillSimVertexHistosByPdg(sv.motherPdgId, isMatched, value);
  };

  fillBundle(ch.h_decayLength, sv.decayLength);
  fillBundle(ch.h_r, sv.r);
  fillBundle(ch.h_nTracks, sv.nCharged);
  fillBundle(ch.h_eta, sv.eta());
  fillBundle(ch.h_mass, sv.mass());
  fillBundle(ch.h_pt, sv.pt());
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
  const double eta = sv.eta();
  const double pt = 0.;

  // Position residuals (reco - sim) and pulls
  // Pulls require vertex position errors — use covariance diagonal if available.
  // For now fill residuals; pulls are set to residual / 1 as placeholder
  // until covariance is wired through GeneralRecoVertex.
  // TODO: wire position errors from reco vertex covariance matrix.
  const double xRes = rv.x - sv.x;
  const double yRes = rv.y - sv.y;
  const double zRes = rv.z - sv.z;
  const double lRes = rv.decayLength - sv.decayLength;
  const double xPull = xRes;  // placeholder
  const double yPull = yRes;
  const double zPull = zRes;
  const double lPull = lRes;  // placeholder

  ch.h_xRes.fill(decayLen, r, eta, pt, nTrk, xRes, xPull);
  ch.h_yRes.fill(decayLen, r, eta, pt, nTrk, yRes, yPull);
  ch.h_zRes.fill(decayLen, r, eta, pt, nTrk, zRes, zPull);
  ch.h_decayLengthRes.fill(decayLen, r, eta, pt, nTrk, lRes, lPull);

  // Mass residual — only for CPC vertices
  if (rv.mass.has_value() && sv.decayLength > 0.) {
    const double mRes = rv.mass.value() - sv.mass();
    const double mPull = mRes;  // placeholder
    ch.h_massRes.fill(decayLen, r, eta, pt, nTrk, mRes, mPull);
  }
}

// =============================================================================
// analyzeImpl — shared logic for both vertex types
// =============================================================================

template <typename SimToRecoAssociationType, typename RecoToSimAssociationType>
void SecondaryVertexAnalyzerAlgo::analyzeImpl(std::vector<RecoSecondaryVertex> recoSVs,
                                              const RecoToSimAssociationType &recoToSim,
                                              const SimToRecoAssociationType &simToReco,
                                              const std::string &collectionLabel) {
  // ------------------------------------------------------------------
  // 1. Reset the reco dependent members of Sim SVs
  // ------------------------------------------------------------------
  resetSimSVs();

  // ------------------------------------------------------------------
  // 2. Build reco↔sim matching using pre-produced association maps
  // ------------------------------------------------------------------
  // matchSim2RecoVertices and matchReco2SimVertices are called after the
  // association maps are available. In the current architecture, the maps
  // are produced upstream and consumed by the plugin, which passes the
  // already-fetched map handles. The matching functions below receive the
  // collection-level SimToReco / RecoToSim maps.
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

void SecondaryVertexAnalyzerAlgo::analyze(const edm::View<reco::Vertex> &recoVertices,
                                          const RecoToSimCollectionVtx &recoToSim,
                                          const SimToRecoCollectionVtx &simToReco,
                                          const std::string &collectionLabel) {
  analyzeImpl(buildRecoSVs(recoVertices), recoToSim, simToReco, collectionLabel);
}

void SecondaryVertexAnalyzerAlgo::analyze(const edm::View<reco::VertexCompositePtrCandidate> &recoVertices,
                                          const RecoToSimCollectionCPC &recoToSim,
                                          const SimToRecoCollectionCPC &simToReco,
                                          const std::string &collectionLabel) {
  analyzeImpl(buildRecoSVs(recoVertices), recoToSim, simToReco, collectionLabel);
}
