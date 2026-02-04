#ifndef PixelTrackTorchHighPuritySelectorKernels_h
#define PixelTrackTorchHighPuritySelectorKernels_h

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
#include "RecoTracker/FinalTrackSelectors/interface/PixelRecHitFeaturesSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    void launchCAPreselection(
        Queue& queue,
        const int maxNumberOfTracks,
        const int minNumberOfHits,
        const ::pixelTrack::Quality minimumTrackQuality,
        const ::reco::TrackSoAConstView tracks,
        int* inputTrackIndices,
        int* preselectionOffsets,
        int* nSelectedTracks
    );

    void launchFeaturesExtractor(
        Queue& queue,
        const int maxPreselectedTracks,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        const ::reco::TrackingRecHitConstView hits,
        PixelTrackFeaturesSoA::View trackFeatures,
        RecHitFeatures::PixelRecHitFeaturesSoA::View hitFeatures,
        const int* nSelectedTracks,
        int * nKeptHits,
        int* inputTrackIndices
    );

    void launchScoreFilter(
        Queue& queue,
        const int maxPreselectedTracks,
        const double scoreThreshold,
        int* inputTrackIndices,
        int* nSelectedTracks,
        int* nKeptHits,
        const PixelTrackScoresSoA::View trackScores
    );

    void launchPixelTrackFilter(
        Queue& queue,
        const int maxTracksPreselection,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        int* inputTrackIndices,
        const int* nSelectedTracks,
        const int* nKeptHits,
        ::reco::TrackSoAView tracks_out,
        ::reco::TrackHitSoAView track_hits_out
    );

    reco::TracksSoACollection launchProduceOutputTracks(
        Queue& queue,
        const int maxPreselectedTracks,
        const int avgHitsPerTrack,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        int* inputTrackIndices,
        const int* nSelectedTracks,
        const int* nKeptHits
    );
}

#endif