#ifndef RecoTracker_FinalTrackSelectors_PixelRecHitFeaturesSoA_h
#define RecoTracker_FinalTrackSelectors_PixelRecHitFeaturesSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(PixelRecHitFeaturesSoALayout,
                    SOA_COLUMN(float, x),
                    SOA_COLUMN(float, y),
                    SOA_COLUMN(float, z),
                    SOA_COLUMN(float, xError),
                    SOA_COLUMN(float, yError),
                    SOA_COLUMN(float, zError),
                    SOA_COLUMN(float, r),
                    SOA_COLUMN(float, eta),
                    SOA_COLUMN(float, phi)
);

using PixelRecHitFeaturesSoA = PixelRecHitFeaturesSoALayout<>;

#endif