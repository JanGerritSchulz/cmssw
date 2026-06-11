#ifndef Validation_RecoVertex_SecondaryVertexAnalyzerAlgo_h
#define Validation_RecoVertex_SecondaryVertexAnalyzerAlgo_h

// Package:    Validation/RecoVertex
// Class:      SecondaryVertexAnalyzerAlgo
//
/**\class SecondaryVertexAnalyzerAlgo
   Validation/RecoVertex/interface/SecondaryVertexAnalyzerAlgo.h

 Description: Algorithm class for secondary vertex validation. Owns all
              histogram booking and filling logic, completely decoupled from
              the EDM framework. Receives already-fetched collections and
              association maps from the plugin (SecondaryVertexAnalyzer.cc).

 Deliberately contains no EDM includes. All framework interaction is the
 responsibility of the plugin.

 Sim vertex handling philosophy
 ───────────────────────────────
 Two distinct sim vertex collections are maintained per event, mirroring the
 approach used in track validation:

   allSimSVs_     All non-PV TrackingVertices, with at least one charged 
                  daughter particle of a certain min pT but no physics selection
                  applied. Used as the truth reference for fake rate, duplicate
                  rate, and pileup rate estimates, so that these rates reflect
                  the full landscape of true secondary vertices in the event.

   signalSimSVs_  Subset of allSimSVs_ passing the configured signal selection
                  (minReconstructableDaughters, absEtaMax, PDG ID filter if
                  set). Used as the denominator for efficiency estimates.

 Variable-blind cut suppression
 ───────────────────────────────
 Efficiency plots whose x-axis is the same quantity as a selection cut must
 NOT apply that cut, otherwise the efficiency is trivially 1 or 0 at the
 boundaries. This is handled via the ReconstructabilityFlags bitmask: each
 monitoring bundle is associated with a set of flags that identifies which
 cut(s) to suppress when evaluating reconstructability for that bundle's
 sim-side fills. The isReconstructable() predicate accepts a SkipCuts value
 and bypasses the corresponding checks.

 Example:
   h_decayLength bundle  → SkipCuts::kDecayLength  (do not apply minDecayLength)
   h_nTracks bundle      → SkipCuts::kNDaughters   (do not apply minReconstructableDaughters)
   h_eta bundle          → SkipCuts::kEta          (do not apply absEtaMax)
   all other bundles     → SkipCuts::kNone         (apply all cuts)

 Original Author: Jan Schulz
*/

#include <cstdint>
#include <map>
#include <string>
#include <vector>

// DQM
#include "DQMServices/Core/interface/DQMStore.h"

// Vertex data formats
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// Sim truth
#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

// Internal types and bundles
#include "Validation/RecoVertex/interface/SVMonitoringBundle.h"
#include "Validation/RecoVertex/interface/SVResolutionBundle.h"
#include "Validation/RecoVertex/interface/SVValidationStructs.h"

class SecondaryVertexAnalyzerAlgo {
public:
  using IBooker = dqm::reco::DQMStore::IBooker;

  // =========================================================================
  // Configuration
  // =========================================================================

  /// Configuration struct — populated from the ParameterSet by the plugin
  /// and passed here so the algo has no PSet dependency.
  struct Config {
    std::string rootFolder;
    bool verbose;
    bool doGenericSimPlots;  // book/fill collection-independent sim plots
    bool doPerPdgPlots;      // book per-b/c/other efficiency breakdowns

    // Signal selection cuts applied to build signalSimSVs_.
    // Each cut is individually suppressed for the monitoring bundle
    // whose x-axis is that quantity (see SkipCuts below).
    double minDecayLength;            // minimum 3D decay length [cm]
    int minReconstructableDaughters;  // minimum charged daughters
    double absEtaMax;                 // maximum |eta| of SV position

    // Optional PDG ID filter: if non-empty, only sim SVs whose mother PDG ID
    // (absolute value) appears in this list are included in signalSimSVs_.
    // Empty means no PDG filter (accept all).
    std::vector<int> signalPdgIds;
  };

  // =========================================================================
  // Reconstructability cut suppression
  // =========================================================================

  /// Bitmask identifying which reconstructability cuts to suppress.
  /// Used to implement variable-blind efficiency plots: the bundle whose
  /// x-axis is quantity X suppresses the cut on X when evaluating
  /// reconstructability for its sim-side fills.
  enum SkipCuts : uint32_t {
    kNone = 0,
    kDecayLength = 1 << 0,  // suppress minDecayLength cut
    kNDaughters = 1 << 1,   // suppress minReconstructableDaughters cut
    kEta = 1 << 2,          // suppress absEtaMax cut
    kPdgId = 1 << 3,        // suppress PDG ID filter
  };

  /// Returns true if sv passes the reconstructability criteria, optionally
  /// suppressing individual cuts as indicated by skipCuts.
  bool isReconstructable(const simSecondaryVertex &sv, uint32_t skipCuts = kNone) const;

  // =========================================================================
  // Public interface
  // =========================================================================

  explicit SecondaryVertexAnalyzerAlgo(const Config &cfg);
  ~SecondaryVertexAnalyzerAlgo() = default;

  /// Called from DQMEDAnalyzer::bookHistograms.
  void bookHistograms(IBooker &ibook, const std::vector<std::string> &collectionLabels);

  /// Per-event entry point for reco::Vertex collections (track-based SVs).
  void analyze(const edm::View<reco::Vertex> &recoVertices,
               const TrackingVertexCollection &simVertices,
               const reco::VertexToTrackingVertexAssociator<std::vector<reco::Vertex>> &associator,
               const reco::RecoToSimCollection &trackRecoToSim,
               const reco::SimToRecoCollection &trackSimToReco,
               const std::string &collectionLabel);

  /// Per-event entry point for reco::VertexCompositePtrCandidate collections.
  void analyze(const edm::View<reco::VertexCompositePtrCandidate> &recoVertices,
               const TrackingVertexCollection &simVertices,
               const reco::VertexToTrackingVertexAssociator<std::vector<reco::VertexCompositePtrCandidate>> &associator,
               const reco::RecoToSimCollection &trackRecoToSim,
               const reco::SimToRecoCollection &trackSimToReco,
               const std::string &collectionLabel);

private:
  // =========================================================================
  // Sim vertex building
  // =========================================================================

  /// Build the full sim SV list: all non-PV TrackingVertices with decay
  /// length and mother PDG ID populated. No signal selection applied.
  /// Used as truth reference for fake/duplicate/pileup rate estimates.
  std::vector<simSecondaryVertex> buildAllSimSVs(const TrackingVertexCollection &simVertices) const;

  /// Apply signal selection to allSimSVs to produce the efficiency
  /// denominator. Applies minDecayLength, minReconstructableDaughters,
  /// absEtaMax, and signalPdgIds from Config.
  std::vector<simSecondaryVertex *> buildSignalSimSVs(std::vector<simSecondaryVertex> &allSimSVs) const;

  /// Walk the TrackingParticle parent chain to find the PDG ID of the
  /// decaying mother particle (stopping at B/D hadron level).
  int motherPdgId(const TrackingVertex &tv) const;

  /// Compute 3D decay length w.r.t. the hard-scatter primary vertex.
  /// Returns -1 if no PV can be identified in the sim collection.
  double decayLength(const TrackingVertex &tv, const TrackingVertexCollection &allSimVertices) const;

  // =========================================================================
  // Reco vertex building
  // =========================================================================

  std::vector<recoSecondaryVertex> buildRecoSVs(const edm::View<reco::Vertex> &recoVertices) const;

  std::vector<recoSecondaryVertex> buildRecoSVs(const edm::View<reco::VertexCompositePtrCandidate> &recoVertices) const;

  // =========================================================================
  // Association and matching
  // =========================================================================

  /// Sim→Reco direction: populates simSecondaryVertex::num_matched_reco_vertices
  /// and matched_reco_shared_fractions. Operates on allSimSVs so that all
  /// true SVs, including pileup, are considered for matching.
  void matchSim2RecoVertices(std::vector<simSecondaryVertex> &allSimSVs,
                             const reco::VertexSimToRecoCollection &simToReco) const;

  /// Reco→Sim direction: populates recoSecondaryVertex::kind_of_vertex,
  /// sim_vertices, and sim_vertices_shared_fraction.
  ///
  /// The full allSimSVs collection is used here so that:
  ///   - A reco SV matched only to a pileup sim SV is correctly flagged
  ///     as pileup rather than fake.
  ///   - A reco SV matched to multiple sim SVs (merged) is correctly
  ///     identified even if only one of those sim SVs passes signal selection.
  ///
  /// The signalSimSVs pointer set is used only to distinguish, for matched
  /// reco SVs, whether the matched sim SV is a signal vertex or not.
  void matchReco2SimVertices(std::vector<recoSecondaryVertex> &recoSVs,
                             const reco::VertexRecoToSimCollection &recoToSim,
                             const std::vector<simSecondaryVertex> &allSimSVs,
                             const std::vector<simSecondaryVertex *> &signalSimSVs) const;

  // =========================================================================
  // Histogram filling
  // =========================================================================

  /// Fill sim-side histograms for one simSecondaryVertex.
  /// Each bundle is filled with its associated SkipCuts mask applied to the
  /// reconstructability evaluation — this is the variable-blind mechanism.
  void fillSimVertexHistograms(const std::string &label, const simSecondaryVertex &sv);

  /// Fill reco-side histograms for one recoSecondaryVertex.
  void fillRecoVertexHistograms(const std::string &label, const recoSecondaryVertex &rv);

  /// Fill resolution/pull histograms for a matched reco-sim pair.
  void fillResolutionHistograms(const std::string &label, const recoSecondaryVertex &rv, const simSecondaryVertex &sv);

  // =========================================================================
  // Shared implementation
  // =========================================================================

  /// Internal template called by both public analyze() overloads after
  /// type-specific buildRecoSVs() has been called.
  template <typename VertexCollection, typename AssociatorType>
  void analyzeImpl(std::vector<recoSecondaryVertex> recoSVs,
                   const TrackingVertexCollection &simVertices,
                   const AssociatorType &associator,
                   const reco::RecoToSimCollection &trackRecoToSim,
                   const reco::SimToRecoCollection &trackSimToReco,
                   const std::string &collectionLabel);

  // =========================================================================
  // Histogram storage
  // =========================================================================

  Config cfg_;

  // Plain MonitorElement* histograms not belonging to a bundle,
  // keyed by [collectionLabel][histogramName].
  std::map<std::string, std::map<std::string, dqm::reco::MonitorElement *>> mes_;

  // Per-collection bundle structs. Each SVMonitoringBundle is associated
  // with a SkipCuts mask that is applied when evaluating reconstructability
  // during sim-side fills — this implements variable-blind cut suppression.
  struct BundleWithCutMask {
    SVMonitoringBundle bundle{};
    uint32_t skipCuts = kNone;  // which reconstructability cuts to suppress
  };

  struct CollectionHistograms {
    // Efficiency / fake rate monitoring bundles.
    // skipCuts encodes which cut is suppressed for each bundle's x-axis.
    BundleWithCutMask h_decayLength{.skipCuts = kDecayLength};
    BundleWithCutMask h_decayLengthSig{.skipCuts = kDecayLength};
    BundleWithCutMask h_r{.skipCuts = kNone};
    BundleWithCutMask h_nTracks{.skipCuts = kNDaughters};
    BundleWithCutMask h_eta{.skipCuts = kEta};
    BundleWithCutMask h_chi2ndof{.skipCuts = kNone};
    // mass only available for CPC — bundle present but unfilled for reco::Vertex
    BundleWithCutMask h_mass{.skipCuts = kNone};

    // Resolution bundles — filled only for matched reco-sim pairs,
    // no reconstructability evaluation needed here.
    SVResolutionBundle h_xRes;
    SVResolutionBundle h_yRes;
    SVResolutionBundle h_zRes;
    SVDecayLengthBundle h_decayLengthRes;
    SVResolutionBundle h_massRes;  // CPC only
  };

  std::map<std::string, CollectionHistograms> collectionHistos_;

  // Generic sim-side histograms booked once (collection-independent).
  // Only populated when cfg_.doGenericSimPlots is true.
  struct GenericSimHistograms {
    dqm::reco::MonitorElement *h_decayLength = nullptr;
    dqm::reco::MonitorElement *h_r = nullptr;
    dqm::reco::MonitorElement *h_nDaughters = nullptr;
    dqm::reco::MonitorElement *h_motherPdgId = nullptr;
    dqm::reco::MonitorElement *h_numAllSimSVs = nullptr;
    dqm::reco::MonitorElement *h_numSignalSimSVs = nullptr;
  } genericSimHistos_;
};

#endif  // Validation_RecoVertex_SecondaryVertexAnalyzerAlgo_h
