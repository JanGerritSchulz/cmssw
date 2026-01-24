#ifndef RecoTracker_FinalTrackSelectors_alpaka_PixelRecHitFeaturesDeviceCollection_h
#define RecoTracker_FinalTrackSelectors_alpaka_PixelRecHitFeaturesDeviceCollection_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "RecoTracker/FinalTrackSelectors/interface/PixelRecHitFeaturesSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using PixelRecHitFeaturesOnDevice = PortableCollection<RecHitFeatures::PixelRecHitFeaturesSoA>;
} 

#endif