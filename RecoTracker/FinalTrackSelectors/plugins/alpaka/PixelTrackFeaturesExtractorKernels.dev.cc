#include <alpaka/alpaka.hpp>
#include <xtd/math/sqrt.h>
#include <limits>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackFeaturesExtractorKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using PixelTrackFeaturesSoAView = PixelTrackFeaturesSoA::View;

    struct CAPreselectionKernel{
        // Kernel used to preselect tracks which were flagged as good enough by the CA
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxTracks,
            const int maxTracksPreselection,
            const int minNumberOfHits,
            const ::pixelTrack::Quality minQuality,
            const ::reco::TrackSoAConstView tracks,
            int* nKeptTracks,
            int* oldIndex
        ) const {
            const int i = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
            if(i==0){
                printf("PixelTrackFeaturesExtractorKernels: nTracks=%d\n", tracks.nTracks());
                assert(tracks.nTracks()<maxTracks);
            }

            if(i<maxTracksPreselection)
                oldIndex[i] = -1;

            if(i < tracks.nTracks()){
                if(tracks[i].quality() >= minQuality && nHits(tracks, i)>=minNumberOfHits){
                    int idx = alpaka::atomicOp<alpaka::AtomicAdd>(acc, nKeptTracks, 1);
                    oldIndex[idx] = i;
                }
            }
        } 
    };


    struct FeaturesExtractorKernel{ 
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxTracksPreselection,
            const ::reco::TrackSoAConstView tracks,
            PixelTrackFeaturesSoAView trackFeatures,
            int* oldIndex
        ) const {

            const int i = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];

            if (i < maxTracksPreselection){
                int idx = oldIndex[i];
                if(idx >= 0){
                    constexpr int cPhi = 0;
                    constexpr int cTip = 5;
                    constexpr int cInvPt = 9;
                    constexpr int cZip = 14;
                    
                    const auto track = tracks[idx];
                    const auto cov   = track.covariance();
                    const auto state = track.state();
                    const float pt   = track.pt();
                    const int ndof   = nHits(tracks, idx) * 2 - 5;
                    const float ptError = xtd::sqrt(cov(cInvPt)) * pt * pt;

                    trackFeatures.chi2(i)     = track.chi2() * ndof;
                    trackFeatures.dzError(i)  = xtd::sqrt(cov(cZip));
                    trackFeatures.dxyError(i) = xtd::sqrt(cov(cTip));
                    trackFeatures.eta(i)      = track.eta();
                    trackFeatures.ndof(i)     = ndof;
                    trackFeatures.phi(i)      = track.state()(0);
                    trackFeatures.phiError(i) = xtd::sqrt(cov(cPhi));
                    trackFeatures.pt(i)       = pt;
                    trackFeatures.ptError(i)  = ptError;
                    trackFeatures.qoverp(i)   = state(2);
                    trackFeatures.dzBS(i)     = state(4);
                    trackFeatures.dxyBS(i)    = state(1);
                }
                else{
                    const auto NaN = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.chi2(i)     = NaN;
                    trackFeatures.dzError(i)  = NaN;
                    trackFeatures.dxyError(i) = NaN;
                    trackFeatures.eta(i)      = NaN;
                    trackFeatures.ndof(i)     = -1;
                    trackFeatures.phi(i)      = NaN;
                    trackFeatures.phiError(i) = NaN;
                    trackFeatures.pt(i)       = NaN;
                    trackFeatures.ptError(i)  = NaN;
                    trackFeatures.qoverp(i)   = NaN;
                    trackFeatures.dzBS(i)     = NaN;
                    trackFeatures.dxyBS(i)    = NaN;
                }
            }
        } 
    };

    void launchCAPreselectionKernel(
        Queue& queue,
        const int maxTracks,
        const int maxTracksPreselection,
        const int minNumberOfHits,
        const ::pixelTrack::Quality minQuality,
        const ::reco::TrackSoAConstView tracks,
        int* nKeptTracks,
        int* oldIndex)
    {
        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            CAPreselectionKernel{},
            maxTracks,
            maxTracksPreselection,
            minNumberOfHits,
            minQuality,
            tracks,
            nKeptTracks,
            oldIndex
        );
    }

    void launchFeaturesExtractorKernel(
        Queue& queue,
        const int maxTracksPreselection,
        const ::reco::TrackSoAConstView tracks,
        PixelTrackFeaturesSoA::View trackFeatures,
        int* oldIndex)
    {
        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxTracksPreselection, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            FeaturesExtractorKernel{},
            maxTracksPreselection,
            tracks,
            trackFeatures,
            oldIndex
        );
    }
}