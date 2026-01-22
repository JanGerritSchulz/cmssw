#ifndef RecoTracker_FinalTrackSelectors_PixelTrackFeaturesSoA_h
#define RecoTracker_FinalTrackSelectors_PixelTrackFeaturesSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

GENERATE_SOA_LAYOUT(PixelTrackFeaturesSoALayout,
                    SOA_COLUMN(float, chi2), 
                    SOA_COLUMN(float, dzError),
                    SOA_COLUMN(float, dxyError),
                    SOA_COLUMN(float, eta),
                    SOA_COLUMN(int, ndof),
                    SOA_COLUMN(float, phi),
                    SOA_COLUMN(float, phiError),
                    SOA_COLUMN(float, pt),
                    SOA_COLUMN(float, ptError),
                    SOA_COLUMN(float, qoverp),
                    SOA_COLUMN(float, dzBS),
                    SOA_COLUMN(float, dxyBS)
);

using PixelTrackFeaturesSoA = PixelTrackFeaturesSoALayout<>;

// Define the SoA layout for track scores (output)
GENERATE_SOA_LAYOUT(PixelTrackScoresSoALayout, SOA_COLUMN(float, score))

using PixelTrackScoresSoA = PixelTrackScoresSoALayout<>;

#endif
//                    SOA_COLUMN(int, charge),
//                    SOA_COLUMN(float, dxy),
//                    SOA_COLUMN(float, dz),
//                    SOA_COLUMN(float, dsz),
//                    SOA_COLUMN(float, dszError),
//                    SOA_COLUMN(float, etaError),
//                    SOA_COLUMN(float, lambdaError),
//                    SOA_COLUMN(float, qoverpError),
//                    SOA_COLUMN(float, vx),
//                    SOA_COLUMN(float, vy),
//                    SOA_COLUMN(float, vz),