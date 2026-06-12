#ifndef Validation_RecoVertex_SVValidationStructs_h
#define Validation_RecoVertex_SVValidationStructs_h

// Package:    Validation/RecoVertex
//
/**\struct SVValidationStructs Validation/RecoVertex/interface/SVValidationStructs.h

 Description: Internal structs representing simulated and reconstructed secondary
              vertices as used by the SecondaryVertexAnalyzer. These are lightweight
              analysis-time objects built from EDM types and association maps; they
              are never stored in the event.

 Original Author: Jan Schulz
*/

#include <vector>
#include <optional>
#include <cmath>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

// =============================================================================
// SimSecondaryVertex
//
// Represents a TrackingVertex that is considered as a secondary vertex truth
// candidate. Built from a TrackingVertexRef in SecondaryVertexAnalyzer::getSimSVs().
// =============================================================================

struct SimSecondaryVertex {
  SimSecondaryVertex(double x1, double y1, double z1)
      : x(x1),
        y(y1),
        z(z1),
        r(std::sqrt(x1 * x1 + y1 * y1)),
        decayLength(-1.),
        nCharged(0),
        nReconstructable(0),
        num_matched_reco_vertices(0),
        average_match_quality(0.f),
        motherPdgId(0),
        isFromPileup(false) {}

  // Position
  double x, y, z;
  double r;  // transverse decay radius

  // Decay geometry — filled after PV association
  double decayLength;  // 3D decay length [cm]

  // Daughter track multiplicity
  int nCharged;          // number of charged daughter TrackingParticles
  int nReconstructable;  // number of daughters with sufficient hits to be reconstructable

  // Matching to reco
  int num_matched_reco_vertices;
  float average_match_quality;
  std::vector<float> matched_reco_shared_fractions;

  // Generator-level information
  int motherPdgId;    // PDG ID of the immediate decaying particle
  bool isFromPileup;  // true if this vertex comes from a pileup interaction

  // Event identification
  EncodedEventId eventId;

  // Reference to the underlying TrackingVertex
  TrackingVertexRef simVertex;
};

// =============================================================================
// RecoSecondaryVertex
//
// Represents a reconstructed secondary vertex. Templated on the underlying
// CMSSW type so the same struct can be used for both reco::Vertex (track-based)
// and reco::VertexCompositePtrCandidate (PF-based).
// =============================================================================

struct RecoSecondaryVertex {
  // Bitmask flags — consistent with PrimaryVertexAnalyzer4PUSlimmed conventions
  enum VertexProperties {
    NONE = 0,
    MATCHED = 1,
    DUPLICATE = 2,
    FAKE = 4,
    MERGED = 8,
  };

  RecoSecondaryVertex(double x1, double y1, double z1)
      : x(x1),
        y(y1),
        z(z1),
        r(std::sqrt(x1 * x1 + y1 * y1)),
        decayLength(-1.),
        decayLengthSignificance(-1.),
        chi2(-1.),
        ndof(-1.),
        nTracks(0),
        num_matched_sim_vertices(0),
        kind_of_vertex(NONE),
        mass(std::nullopt),
        motherPdgId(std::nullopt),
        isFromPileup(false) {}

  // Position
  double x, y, z;
  double r;  // transverse decay radius

  // Decay geometry — filled after PV association
  double decayLength;
  double decayLengthSignificance;

  // Fit quality
  double chi2;
  double ndof;
  double normalizedChi2() const { return (ndof > 0.) ? chi2 / ndof : -1.; }

  // Track multiplicity
  int nTracks;

  // Matching to sim
  int num_matched_sim_vertices;
  std::vector<const SimSecondaryVertex *> sim_vertices;
  std::vector<float> sim_vertices_shared_fraction;

  // Classification flags (bitmask of VertexProperties)
  int kind_of_vertex;

  // Optional fields — populated for VertexCompositePtrCandidate only
  std::optional<double> mass;

  // Optional fields — populated after MC truth matching
  std::optional<int> motherPdgId;
  bool isFromPileup;

  // Reference to the underlying reco vertex:
  // reco::VertexBaseRef or VertexCompositePtrCandidateRef.
  std::optional<edm::RefToBase<reco::Vertex>> recoVertexRef;
  std::optional<edm::RefToBase<reco::VertexCompositePtrCandidate>> recoVertexCPCRef;

  template <typename VertexType>
  edm::RefToBase<VertexType> recoVertex() const;
};

#endif  // Validation_RecoVertex_SVValidationStructs_h
