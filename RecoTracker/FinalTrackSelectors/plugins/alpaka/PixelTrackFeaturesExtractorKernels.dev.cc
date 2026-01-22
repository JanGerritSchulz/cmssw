#include <alpaka/alpaka.hpp>

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
        
        enum CovIndex {
            cPhi = 0,
            cTip = 5,
            cInvPt = 9,
            cCotanTheta = 12,
            cZip = 14
        };

        auto cov = tracks[i].covariance();
        auto state = tracks[i].state();
        float pt    = tracks[i].pt();
        int ndof  = nHits(tracks, i) * 2 - 5;
        float ptError = alpaka::math::sqrt(acc, cov(cInvPt)) * pt * pt;

        trackFeatures.chi2(i)     = tracks[i].chi2() * ndof;
        trackFeatures.dzError(i)  = alpaka::math::sqrt(acc, cov(cZip));
        trackFeatures.dxyError(i) = alpaka::math::sqrt(acc, cov(cTip));
        trackFeatures.eta(i)      = tracks[i].eta();
        trackFeatures.ndof(i)     = ndof;
        trackFeatures.phi(i)      = tracks[i].state()(0);
        trackFeatures.phiError(i) = alpaka::math::sqrt(acc, cov(cPhi));
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