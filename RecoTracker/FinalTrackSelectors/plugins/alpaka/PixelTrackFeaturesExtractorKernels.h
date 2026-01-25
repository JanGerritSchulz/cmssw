#ifndef PixelTrackFeaturesExtractorKernels_h
#define PixelTrackFeaturesExtractorKernels_h

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
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
        int* oldIndex
    );

    void launchFeaturesExtractorKernel(
        Queue& queue,
        const int maxTracksPreselection,
        const ::reco::TrackSoAConstView tracks,
        PixelTrackFeaturesSoA::View trackFeatures,
        int* oldIndex
    );
}

#endif