#include <alpaka/alpaka.hpp>
#include <xtd/math/sqrt.h>
#include <limits>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackFeaturesExtractorKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using PixelTrackFeaturesSoAView = PixelTrackFeaturesSoA::View;

    struct TrackFeatureKernel{ 
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxTracks,
            const int minNumberOfHits,
            const ::pixelTrack::Quality minQuality,
            const ::reco::TrackSoAConstView tracks,
            PixelTrackFeaturesSoAView trackFeatures,
            int* nKeptTracks,
            int* newIndex,
            int* nHitsPerKeptTrack
        ) const {

            const int i = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
            if(i==0){
                printf("PixelTrackFeaturesExtractorKernels: nTracks=%d\n", tracks.nTracks());
                assert(tracks.nTracks()<maxTracks);
            }
                
            if (i < tracks.nTracks()){
                newIndex[i] = -1;
                nHitsPerKeptTrack[i] = 0;
                if(tracks[i].quality() >= minQuality && nHits(tracks, i)>=minNumberOfHits){
                    int idx = alpaka::atomicOp<alpaka::AtomicAdd>(acc, nKeptTracks, 1);
                    newIndex[i] = idx;
                    nHitsPerKeptTrack[i] = nHits(tracks, i);

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
                else{
                    trackFeatures.chi2(i)     = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.dzError(i)  = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.dxyError(i) = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.eta(i)      = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.ndof(i)     = -1;
                    trackFeatures.phi(i)      = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.phiError(i) = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.pt(i)       = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.ptError(i)  = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.qoverp(i)   = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.dzBS(i)     = std::numeric_limits<float>::quiet_NaN();
                    trackFeatures.dxyBS(i)    = std::numeric_limits<float>::quiet_NaN();
                }
            }

            else if(i < maxTracks){
                newIndex[i] = -1;
                nHitsPerKeptTrack[i] = 0;
                trackFeatures.chi2(i)     = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.dzError(i)  = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.dxyError(i) = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.eta(i)      = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.ndof(i)     = -1;
                trackFeatures.phi(i)      = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.phiError(i) = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.pt(i)       = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.ptError(i)  = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.qoverp(i)   = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.dzBS(i)     = std::numeric_limits<float>::quiet_NaN();
                trackFeatures.dxyBS(i)    = std::numeric_limits<float>::quiet_NaN();
            }
            else
                return;
        } 
    };

    struct CompactKernel{ 
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc, 
            const ::reco::TrackSoAConstView& tracks,
            ::reco::TrackSoAView& outTracks,
            const int* newIndex,
            const int* outHitOffsets
        ) const {
            const int i = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0]; 
            if (i >= tracks.nTracks()) return;
            int newTk = newIndex[i];
            if (newTk < 0) return; 
            
            outTracks[newTk] = tracks[i];
            outTracks[newTk].hitOffsets() = outHitOffsets[newTk];

            // copy hits
            uint32_t inStart  = (i == 0) ? 0 : tracks[i-1].hitOffsets();
            uint32_t inEnd    = tracks[i].hitOffsets();
            uint32_t outStart = (newTk == 0) ? 0 : outHitOffsets[newTk-1];

            for (uint32_t h = 0; h < (inEnd - inStart); ++h) {
                outHits[outStart + h].id()    = inHits[inStart + h].id();
                outHits[outStart + h].detId() = inHits[inStart + h].detId();
            }

        }
    };

    void launchTrackFeatureExtractorKernel(
        Queue& queue,
        const int maxTracks,
        const int minNumberOfHits,
        const ::pixelTrack::Quality minQuality,
        const ::reco::TrackSoAConstView& tracks,
        PixelTrackFeaturesSoAView& trackFeatures,
        int* nKeptTracks,
        int* newIndex,
        int* nHitsPerKeptTrack)
    {
        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            TrackFeatureKernel{},
            maxTracks,
            minNumberOfHits,
            minQuality,
            tracks,
            trackFeatures,
            nKeptTracks,
            newIndex,
            nHitsPerKeptTrack
        );
    }

    void launchCompactKernel(
        Queue& queue,
        const int maxTracks, 
        const ::reco::TrackSoAConstView& tracks,
        ::reco::TrackSoAView& outTracks,
        const int* newIndex,
        const int* outHitOffsets)
    {
        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            CompactKernel{},
            tracks,
            outTracks,
            newIndex,
            outHitOffsets
        );
    }
}