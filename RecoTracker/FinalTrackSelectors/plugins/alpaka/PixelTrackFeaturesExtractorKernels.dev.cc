#include <alpaka/alpaka.hpp>
#include <xtd/math/sqrt.h>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackFeaturesExtractorKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using PixelTrackFeaturesSoAView = PixelTrackFeaturesSoA::View;

    struct TrackFeatureKernel{ 
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const ::reco::TrackSoAConstView tracks,
            PixelTrackFeaturesSoAView trackFeatures) const {

        const int i = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
        if (i >= tracks.nTracks())
            return;
        
        constexpr int cPhi = 0;
        constexpr int cTip = 5;
        constexpr int cInvPt = 9;
        constexpr int cZip = 14;

        auto cov = tracks[i].covariance();
        auto state = tracks[i].state();
        float pt    = tracks[i].pt();
        int ndof  = nHits(tracks, i) * 2 - 5;
        float ptError = xtd::sqrt(cov(cInvPt)) * pt * pt;

        trackFeatures.chi2(i)     = tracks[i].chi2() * ndof;
        trackFeatures.dzError(i)  = xtd::sqrt(cov(cZip));
        trackFeatures.dxyError(i) = xtd::sqrt(cov(cTip));
        trackFeatures.eta(i)      = tracks[i].eta();
        trackFeatures.ndof(i)     = ndof;
        trackFeatures.phi(i)      = tracks[i].state()(0);
        trackFeatures.phiError(i) = xtd::sqrt(cov(cPhi));
        trackFeatures.pt(i)       = pt;
        trackFeatures.ptError(i)  = ptError;
        trackFeatures.qoverp(i)   = state(2);
        trackFeatures.dzBS(i)     = state(4);
        trackFeatures.dxyBS(i)    = state(1);
        }
    };

  void launchTrackFeatureExtractorKernel(
      Queue& queue,
      const ::reco::TrackSoAConstView& tracks,
      PixelTrackFeaturesSoAView& trackFeatures) {
    
    const int nTracks = tracks.nTracks();

    constexpr uint32_t threadsPerBlock = 256;
    const uint32_t blocks =
        cms::alpakatools::divide_up_by(nTracks, threadsPerBlock);

    const auto workDiv =
        cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

    alpaka::exec<Acc1D>(
        queue,
        workDiv,
        TrackFeatureKernel{},
        tracks,
        trackFeatures);
  }
}