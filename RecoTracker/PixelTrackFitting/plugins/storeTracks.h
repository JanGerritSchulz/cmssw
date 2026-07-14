#ifndef RecoTracker_PixelTrackFitting_plugins_storeTracks_h
#define RecoTracker_PixelTrackFitting_plugins_storeTracks_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "RecoTracker/PixelTrackFitting/interface/TracksWithHits.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

template <typename Ev, typename TWH>
void storeTracks(Ev& ev,
                 TWH& tracksWithHits,
                 const TrackerTopology& ttopo,
                 const bool fillFullTrackExtra = false,
                 Propagator const* propagator = nullptr,
                 MagneticField const* magField = nullptr) {
  auto tracks = std::make_unique<reco::TrackCollection>();
  auto recHits = std::make_unique<TrackingRecHitCollection>();
  auto trackExtras = std::make_unique<reco::TrackExtraCollection>();

  int cc = 0, nTracks = tracksWithHits.size();

  trackExtras->resize(nTracks);
  tracks->reserve(nTracks);
  recHits->reserve(4 * nTracks);

  for (int i = 0; i < nTracks; i++) {
    reco::Track* track = tracksWithHits[i].first;
    auto& hits = tracksWithHits[i].second;

    if constexpr (std::is_same_v<TWH, pixeltrackfitting::TracksWithRecHits>) {
      // Sort hits by increasing global radius so the KF fitter sees them
      // in inside-out order (it assumes hits are ordered along the track).
      std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        return a->globalPosition().perp2() < b->globalPosition().perp2();
      });
    }

    for (unsigned int k = 0; k < hits.size(); k++) {
      auto* hit = hits[k]->clone();
      track->appendHitPattern(*hit, ttopo);
      recHits->push_back(hit);
    }
    tracks->push_back(*track);
    delete track;
  }

  LogDebug("TrackProducer") << "put the collection of TrackingRecHit in the event\n";
  edm::OrphanHandle<TrackingRecHitCollection> ohRH = ev.put(std::move(recHits));

  edm::RefProd<TrackingRecHitCollection> hitCollProd(ohRH);
  for (int k = 0; k < nTracks; k++) {
    auto& aTrackExtra = (*trackExtras)[k];
    const reco::Track& tk = (*tracks)[k];

    unsigned int nHits = tk.numberOfValidHits();
    aTrackExtra.setHits(hitCollProd, cc, nHits);

    AlgebraicVector5 v = AlgebraicVector5(0, 0, 0, 0, 0);
    reco::TrackExtra::TrajParams trajParams(nHits, LocalTrajectoryParameters(v, 1.));
    reco::TrackExtra::Chi2sFive chi2s(nHits, 0);
    aTrackExtra.setTrajParams(std::move(trajParams), std::move(chi2s));

    if (fillFullTrackExtra) {
      // DetIds of the innermost and outermost sorted hits
      const TrackingRecHit& innerHit = (*ohRH)[cc];
      const TrackingRecHit& outerHit = (*ohRH)[cc + nHits - 1];

      const auto& innerHitPos = innerHit.globalPosition();
      const auto& outerHitPos = outerHit.globalPosition();

      aTrackExtra.setInnerDetId(innerHit.geographicalId().rawId());
      aTrackExtra.setOuterDetId(outerHit.geographicalId().rawId());
      aTrackExtra.setSeedDirection(alongMomentum);

      // ---------------------------------------------------------------
      // Inner state: propagate the PCA FreeTrajectoryState to the
      // surface of the innermost hit.
      //
      // The track's reference point, momentum and covariance are defined
      // at the PCA and form a physically consistent state. Propagating
      // this to the innermost hit surface gives the correct starting
      // state for the KF fitter, which assumes the inner state lies on
      // the first hit's surface.
      //
      // Without propagation the position is at the PCA (~beamline) but
      // the fitter uses it directly on the detector surface, producing
      // an inconsistent TSOS that causes the smoother to return 0 hits.
      // ---------------------------------------------------------------
      if (propagator && magField && innerHit.det() && outerHit.det()) {
        // Build a FreeTrajectoryState at the PCA from the track parameters
        GlobalPoint pcaPos(tk.vx(), tk.vy(), tk.vz());
        GlobalVector pcaMom(tk.px(), tk.py(), tk.pz());
        GlobalTrajectoryParameters gtp(pcaPos, pcaMom, tk.charge(), magField);
        CurvilinearTrajectoryError cte(tk.covariance());
        FreeTrajectoryState pcaFTS(gtp, cte);

        // Propagate to the innermost hit surface
        TrajectoryStateOnSurface innerTSOS = propagator->propagate(pcaFTS, innerHit.det()->surface());

        if (innerTSOS.isValid()) {
          // const GlobalPoint& innerPos = innerTSOS.globalPosition();
          const GlobalVector& innerMom = innerTSOS.globalMomentum();
          aTrackExtra.setInnerPosition(math::XYZPoint(innerHitPos.x(), innerHitPos.y(), innerHitPos.z()));
          aTrackExtra.setInnerMomentum(math::XYZVector(innerMom.x(), innerMom.y(), innerMom.z()));
          reco::TrackExtra::CovarianceMatrix covInner = innerTSOS.curvilinearError().matrix();
          aTrackExtra.fillInner(covInner);
        } else {
          // Propagation failed — fall back to PCA state for the inner surface.
          // This is better than leaving the state default-constructed (zeros).
          std::cout << "storeTracks: inner propagation failed, falling back to PCA state";
          aTrackExtra.setInnerPosition(math::XYZPoint(innerHitPos.x(), innerHitPos.y(), innerHitPos.z()));
          aTrackExtra.setInnerMomentum(tk.momentum());
          reco::TrackExtra::CovarianceMatrix covInner = tk.covariance();
          aTrackExtra.fillInner(covInner);
        }

        // Propagate to the outermost hit surface
        TrajectoryStateOnSurface outerTSOS = propagator->propagate(pcaFTS, outerHit.det()->surface());

        if (outerTSOS.isValid()) {
          // const GlobalPoint& outerPos = outerTSOS.globalPosition();
          const GlobalVector& outerMom = outerTSOS.globalMomentum();
          aTrackExtra.setOuterPosition(math::XYZPoint(outerHitPos.x(), outerHitPos.y(), outerHitPos.z()));
          aTrackExtra.setOuterMomentum(math::XYZVector(-outerMom.x(), -outerMom.y(), -outerMom.z()));
          reco::TrackExtra::CovarianceMatrix covOuter = outerTSOS.curvilinearError().matrix();
          aTrackExtra.fillOuter(covOuter);
        } else {
          std::cout << "storeTracks: outer propagation failed, falling back to PCA state";
          aTrackExtra.setOuterPosition(math::XYZPoint(outerHitPos.x(), outerHitPos.y(), outerHitPos.z()));
          aTrackExtra.setOuterMomentum(tk.momentum());
          reco::TrackExtra::CovarianceMatrix covOuter = tk.covariance();
          aTrackExtra.fillOuter(covOuter);
        }

      } else {
        // No propagator available or hit has no det — use PCA state for both
        // ends. This avoids zeros in the covariance but is less accurate.
        std::cout << "storeTracks: no propagator/magField/det, using PCA state for inner/outer";

        aTrackExtra.setInnerPosition(tk.referencePoint());
        aTrackExtra.setInnerMomentum(tk.momentum());
        reco::TrackExtra::CovarianceMatrix covInner = tk.covariance();
        aTrackExtra.fillInner(covInner);

        aTrackExtra.setOuterPosition(tk.referencePoint());
        aTrackExtra.setOuterMomentum(tk.momentum());
        reco::TrackExtra::CovarianceMatrix covOuter = tk.covariance();
        aTrackExtra.fillOuter(covOuter);
      }
    }
    cc += nHits;
  }

  LogDebug("TrackProducer") << "put the collection of TrackExtra in the event"
                            << "\n";
  edm::OrphanHandle<reco::TrackExtraCollection> ohTE = ev.put(std::move(trackExtras));

  for (int k = 0; k < nTracks; k++) {
    const reco::TrackExtraRef theTrackExtraRef(ohTE, k);
    (*tracks)[k].setExtra(theTrackExtraRef);
  }

  ev.put(std::move(tracks));
}

#endif
