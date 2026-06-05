#ifndef SimTracker_VertexAssociation_calculateVertexSharedTracks_h
#define SimTracker_VertexAssociation_calculateVertexSharedTracks_h

#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

struct SharedTracksAndFractions {
  SharedTracksAndFractions(unsigned int nSharedTracks,
                           float sharedTracksFraction,
                           float sharedPt2Fraction,
                           float sharedDzErrFraction)
      : nSharedTracks_(nSharedTracks),
        sharedTracksFraction_(sharedTracksFraction),
        sharedPt2Fraction_(sharedPt2Fraction),
        sharedDzErrFraction_(sharedDzErrFraction) {}

  const unsigned int nSharedTracks_;
  const float sharedTracksFraction_;
  const float sharedPt2Fraction_;
  const float sharedDzErrFraction_;
};

// -----------------------------------------------------------------------------
// reco::Vertex overloads (original)
// -----------------------------------------------------------------------------

SharedTracksAndFractions calculateVertexSharedTracks(const reco::Vertex &recoV,
                                                     const TrackingVertex &simV,
                                                     const reco::RecoToSimCollection &trackRecoToSimAssociation);

SharedTracksAndFractions calculateVertexSharedTracks(const TrackingVertex &simV,
                                                     const reco::Vertex &recoV,
                                                     const reco::SimToRecoCollection &trackSimToRecoAssociation);

// -----------------------------------------------------------------------------
// reco::VertexCompositePtrCandidate overloads
//
// Track extraction from candidate daughters is done via dynamic_cast
// (reco::PFCandidate on RECO/AOD, pat::PackedCandidate on MiniAOD).
// Daughters from which no track can be recovered (neutral particles) are
// excluded from both numerator and denominator — they carry no tracking
// information and should not penalise the shared-track fraction.
//
// The fraction denominators follow the same conventions as the reco::Vertex
// overloads: sharedTracksFraction_ is relative to the number of reco
// daughters with a recoverable track; pt2 and dzError fractions are weighted
// sums over those same daughters.
// -----------------------------------------------------------------------------

SharedTracksAndFractions calculateVertexSharedTracks(const reco::VertexCompositePtrCandidate &recoV,
                                                     const TrackingVertex &simV,
                                                     const reco::RecoToSimCollection &trackRecoToSimAssociation);

SharedTracksAndFractions calculateVertexSharedTracks(const TrackingVertex &simV,
                                                     const reco::VertexCompositePtrCandidate &recoV,
                                                     const reco::SimToRecoCollection &trackSimToRecoAssociation);

#endif  // SimTracker_VertexAssociation_calculateVertexSharedTracks_h
