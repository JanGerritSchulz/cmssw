#ifndef PixelTrackFeaturesExtractorKernels_h
#define PixelTrackFeaturesExtractorKernels_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "RecoTracker/FinalTrackSelectors/interface/PixelTrackFeaturesSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    void launchTrackFeatureExtractorKernel(
        Queue& queue,
        const ::reco::TrackSoAConstView& tracks,
        PixelTrackFeaturesSoA::View& trackFeatures
        );
}

#endif