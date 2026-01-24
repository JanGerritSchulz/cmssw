#ifndef PixelTrackFeaturesExtractorKernels_h
#define PixelTrackFeaturesExtractorKernels_h

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "RecoTracker/FinalTrackSelectors/interface/PixelTrackFeaturesSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    void launchTrackFeatureExtractorKernel(
        Queue& queue,
        const int maxTracks,
        const int minNumberOfHits,
        const ::pixelTrack::Quality minQuality,
        const ::reco::TrackSoAConstView& tracks,
        PixelTrackFeaturesSoA::View& trackFeatures,
        int* nKeptTracks,
        int* newIndex,
        int* nHitsPerKeptTrack
        );
    
    void launchCompactKernel(
        Queue& queue,
        const int maxTracks,
        const ::reco::TrackSoAConstView& tracks,
        ::reco::TrackSoAView& outTracks,
        const int* newIndex,
        const int* outHitOffsets
    );
}

#endif