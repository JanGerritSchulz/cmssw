#ifndef RecoTracker_PixelTrackFitting_plugins_storeTracks_h
#define RecoTracker_PixelTrackFitting_plugins_storeTracks_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "RecoTracker/PixelTrackFitting/interface/TracksWithHits.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

template <typename Ev, typename TWH>
void storeTracks(Ev& ev, const TWH& tracksWithHits, const TrackerTopology& ttopo) {
  // std::cout << "Starting storeTracks " << std::endl;
  auto tracks = std::make_unique<reco::TrackCollection>();
  auto recHits = std::make_unique<TrackingRecHitCollection>();
  auto trackExtras = std::make_unique<reco::TrackExtraCollection>();

  int cc = 0, nTracks = tracksWithHits.size();
  // std::cout << "nTracks = " << nTracks << std::endl;

  trackExtras->resize(nTracks);
  tracks->reserve(nTracks);
  recHits->reserve(4 * nTracks);

  for (int i = 0; i < nTracks; i++) {
    // std::cout << "access track i = " << i << std::endl;
    reco::Track* track = tracksWithHits[i].first;
    // std::cout << "access hits of track i = " << i << std::endl;
    const auto& hits = tracksWithHits[i].second;

  // std::cout << "nHits in track = " << hits.size() << std::endl;
    for (unsigned int k = 0; k < hits.size(); k++) {
      // std::cout << "access hit k = " << k << ": hits[k] = " << hits[k] << std::endl;
      // auto hitRef = hits[k];
      // std::cout << "clone hit k = " << k << ": hits[k]->clone() = " << hits[k]->clone() << std::endl;
      auto* hit = hits[k]->clone();  // need to clone (at least if from SoA)
      // std::cout << "append hit k = " << k << std::endl;
      track->appendHitPattern(*hit, ttopo);
      recHits->push_back(hit);
    }
    // std::cout << "push track i = " << i << std::endl;
    tracks->push_back(*track);
    delete track;
  }

  // std::cout << "put the collection of TrackingRecHit in the event" << std::endl;
  LogDebug("TrackProducer") << "put the collection of TrackingRecHit in the event" << "\n";
  edm::OrphanHandle<TrackingRecHitCollection> ohRH = ev.put(std::move(recHits));

  edm::RefProd<TrackingRecHitCollection> hitCollProd(ohRH);
  for (int k = 0; k < nTracks; k++) {
    auto& aTrackExtra = (*trackExtras)[k];

    //fill the TrackExtra with TrackingRecHitRef
    unsigned int nHits = (*tracks)[k].numberOfValidHits();
    aTrackExtra.setHits(hitCollProd, cc, nHits);
    cc += nHits;
    AlgebraicVector5 v = AlgebraicVector5(0, 0, 0, 0, 0);
    reco::TrackExtra::TrajParams trajParams(nHits, LocalTrajectoryParameters(v, 1.));
    reco::TrackExtra::Chi2sFive chi2s(nHits, 0);
    aTrackExtra.setTrajParams(std::move(trajParams), std::move(chi2s));
  }

  // std::cout << "put the collection of TrackExtra in the event" << std::endl;
  LogDebug("TrackProducer") << "put the collection of TrackExtra in the event" << "\n";
  edm::OrphanHandle<reco::TrackExtraCollection> ohTE = ev.put(std::move(trackExtras));

  for (int k = 0; k < nTracks; k++) {
    const reco::TrackExtraRef theTrackExtraRef(ohTE, k);
    (*tracks)[k].setExtra(theTrackExtraRef);
  }

  ev.put(std::move(tracks));
  // std::cout << "Finished storeTracks " << std::endl;
}

#endif
