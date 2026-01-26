#ifndef PixelTrackFeaturesExtractorKernels_h
#define PixelTrackFeaturesExtractorKernels_h

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/TrackSoA/interface/TracksDevice.h"
#include "DataFormats/TrackSoA/interface/TracksHost.h"
#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"

#include "RecoTracker/FinalTrackSelectors/interface/PixelTrackFeaturesSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    void launchCAPreselectionKernel(
        Queue& queue,
        const int maxTracks,
        const int maxTracksPreselection,
        const int minNumberOfHits,
        const ::pixelTrack::Quality minQuality,
        const ::reco::TrackSoAConstView tracks,
        int* nKeptTracks,
        int* nKeptHits,
        int* oldIndex
    );

    void launchFeaturesExtractorKernel(
        Queue& queue,
        const int maxTracksPreselection,
        const int maxHitsPerTrack,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        const ::reco::TrackingRecHitConstView hits,
        PixelTrackFeaturesSoA::View trackFeatures,
        int* nKeptTracks,
        int* oldIndex
    );

    void launchPixelTrackFilterKernel(
        Queue& queue,
        const int maxTracksPreselection,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        int* oldIndex,
        const int* nKeptTracks,
        const int* nKeptHits,
        ::reco::TrackSoAView tracks_out,
        ::reco::TrackHitSoAView track_hits_out
    );

    reco::TracksSoACollection launchProduceOutputTracks(
        Queue& queue,
        const int maxTracksPreselection,
        const int avgHitsPerTrack,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        int* oldIndex,
        const int* nKeptTracks,
        const int* nKeptHits
    );

    void launchHitOffsetCompactKernel(
        Queue& queue,
        const int maxTracksPreselection,
        const int* nKeptTracks,
        int* nKeptHits
    );
}

#endif