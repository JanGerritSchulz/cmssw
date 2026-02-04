#include <alpaka/alpaka.hpp>

#include <xtd/math/asinh.h>
#include <xtd/math/atan2.h>
#include <xtd/math/sqrt.h>
#include <limits>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/TrackSoA/interface/TracksDevice.h"
#include "DataFormats/TrackSoA/interface/TracksHost.h"
#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"

#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackTorchHighPuritySelectorKernels.h"

#ifndef KERNELS_DEBUG
#define KERNELS_DEBUG
#endif

namespace ALPAKA_ACCELERATOR_NAMESPACE {
    using PixelTrackFeaturesSoAView  = PixelTrackFeaturesSoA::View;
    using PixelRecHitFeaturesSoAView = RecHitFeatures::PixelRecHitFeaturesSoA::View;
    using TrackHitSoA  = ::reco::TrackHitSoA;
    using HitFeaturesIDX = RecHitFeatures::HitFeature;

    struct CAPreselectionKernel{
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxNumberOfTracks,
            const int minNumberOfHits,
            const ::pixelTrack::Quality minimumTrackQuality,
            const ::reco::TrackSoAConstView tracks,
            int* preselectionMask,
            int* trackIndexSoA) const 
        {
            const int trackLimit = std::min(maxNumberOfTracks, tracks.nTracks());
#ifdef KERNELS_DEBUG
            auto worker = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
            if(worker==0){
                printf("nTracks=%d\n", tracks.nTracks());
                if(tracks.nTracks() >= maxNumberOfTracks)
                    printf("PixelTrackTorchHighPuritySelectorKernels Warning: nTracks (%d) >= maxNumberOfTracks (%d)\n", tracks.nTracks(), maxNumberOfTracks);
            }
#endif
            for (auto i : alpaka::uniformElements(acc, trackLimit)){
                trackIndexSoA[i] = i;
                bool isGoodQuality = 
                    tracks[i].quality() >= minimumTrackQuality && 
                    nHits(tracks, i)>=minNumberOfHits;
                preselectionMask[i] = isGoodQuality ? 1 : 0;
            }
        } 
    };

    struct FeaturesExtractorKernel{ 
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxPreselectedTracks,
            const ::reco::TrackSoAConstView tracks,
            const ::reco::TrackHitSoAConstView track_hits,
            const ::reco::TrackingRecHitConstView hits,
            PixelTrackFeaturesSoAView trackFeatures,
            PixelRecHitFeaturesSoAView hitFeatures,
            const int* nKeptTracks,
            int* nKeptHits,
            int* originalTrackIndex
        ) const 
        {
            constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
            const uint32_t nPreselectedTracks = std::min(*nKeptTracks, maxPreselectedTracks);

            for (auto i : alpaka::uniformElements(acc, maxPreselectedTracks)) { 
                // Case 1: valid preselected track --> extract features
                if (i < nPreselectedTracks){
                    int idx = originalTrackIndex[i];
#ifdef KERNELS_DEBUG
                    if (idx<0) printf("PixelTrackTorchHighPuritySelectorKernels: Invalid originalTrackIndex for preselected idx %d\n", i);
#endif
                    // Indices into track covariance matrix
                    constexpr int cPhi = 0;
                    constexpr int cTip = 5;
                    constexpr int cInvPt = 9;
                    constexpr int cZip = 14;
                    
                    const auto track = tracks[idx];
                    const auto cov   = track.covariance();
                    const auto state = track.state();
                    const float pt   = track.pt();
                    const int n_hits = nHits(tracks, idx);
                    const int ndof   = n_hits * 2 - 5;
                    const float ptError = xtd::sqrt(cov(cInvPt)) * pt * pt;

                    nKeptHits[i] = n_hits;

#ifdef KERNELS_DEBUG
                    if(n_hits>RecHitFeatures::MaxHitsPerTrack) printf("PixelTrackTorchHighPuritySelectorKernels: Number of hits (%d) exceeds MaxHitsPerTrack (%d)\n", n_hits, RecHitFeatures::MaxHitsPerTrack);
#endif
                    trackFeatures.chi2(i)     = track.chi2() * ndof;
                    trackFeatures.dzError(i)  = xtd::sqrt(cov(cZip));
                    trackFeatures.dxyError(i) = xtd::sqrt(cov(cTip));
                    trackFeatures.eta(i)      = track.eta();
                    trackFeatures.ndof(i)     = ndof;
                    trackFeatures.phi(i)      = state(0);
                    trackFeatures.phiError(i) = xtd::sqrt(cov(cPhi));
                    trackFeatures.pt(i)       = pt;
                    trackFeatures.ptError(i)  = ptError;
                    trackFeatures.qoverp(i)   = state(2);
                    trackFeatures.dzBS(i)     = state(4);
                    trackFeatures.dxyBS(i)    = state(1);

                    //Prefill hit features:
                    auto hitMatrix = hitFeatures.hits(i);
                    for (int h = 0; h < RecHitFeatures::MaxHitsPerTrack; ++h) {
                        for (int f = 0; f < RecHitFeatures::HitFeatures; ++f) {
                            hitMatrix(h, f) = NaN;
                        }
                    }

                    uint32_t hitBegin = (idx == 0) ? 0 : tracks[idx - 1].hitOffsets();
                    uint32_t hitEnd   = track.hitOffsets();
                    uint32_t nHits   = hitEnd - hitBegin;

                    nHits = std::min(nHits, uint32_t(RecHitFeatures::MaxHitsPerTrack));

                    for (uint32_t h = 0; h < nHits; ++h) {
                        auto hit_id = track_hits[hitBegin + h].id();
                        const auto hit = hits[hit_id];
                        const float x = hit.xGlobal();
                        const float y = hit.yGlobal();
                        const float z = hit.zGlobal();
                        const float r = hit.rGlobal();

                        hitMatrix(h, HitFeaturesIDX::x) = x;
                        hitMatrix(h, HitFeaturesIDX::y) = y;
                        hitMatrix(h, HitFeaturesIDX::z) = z;

                        hitMatrix(h, HitFeaturesIDX::r) = r;                        
                        hitMatrix(h, HitFeaturesIDX::eta) =
                            (r > 0.f) ? xtd::asinh(z / r) : 0.f;
                        hitMatrix(h, HitFeaturesIDX::phi) =
                            xtd::atan2(y, x);

                        hitMatrix(h, HitFeaturesIDX::xErr) =
                            hit.xerrLocal();
                        hitMatrix(h, HitFeaturesIDX::yErr) =
                            hit.yerrLocal();
                    }
                }
                // Case 2: padding entries --> fill with NaNs for inference
                else if (i < (uint32_t) maxPreselectedTracks)
                {
                    originalTrackIndex[i] = -1; // mark as invalid
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

                    auto hitMatrix = hitFeatures.hits(i);
                    for (int h = 0; h < RecHitFeatures::MaxHitsPerTrack; ++h) {
                        for (int f = 0; f < RecHitFeatures::HitFeatures; ++f) {
                            hitMatrix(h, f) = NaN;
                        }
                    }
                }
            }
        }
    };

    struct PixelTrackFilterKernel{ 
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxPreselectedTracks,
            const ::reco::TrackSoAConstView tracks,
            const ::reco::TrackHitSoAConstView track_hits,
            int* originalTrackIndex,
            const int* nKeptTracks,
            const int* nKeptHits,
            ::reco::TrackSoAView tracks_out,
            ::reco::TrackHitSoAView track_hits_out
        ) const 
        {
            const int nTracks = std::min(*nKeptTracks, maxPreselectedTracks);
            if (alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0] == 0) 
                tracks_out.nTracks() = nTracks;

            for (auto i : alpaka::uniformElements(acc, nTracks)) { 
                int inputTrackIdx = originalTrackIndex[i];
                if(inputTrackIdx >= 0){
                    const auto track = tracks[inputTrackIdx];
                    tracks_out[i] = track;
                    tracks_out[i].hitOffsets() = nKeptHits[i];

                    //Access the hits associated to the track:
                    uint32_t hitBegin  = (inputTrackIdx == 0) ? 0 : tracks[inputTrackIdx-1].hitOffsets();
                    uint32_t hitEnd    = track.hitOffsets();
                    uint32_t outStart = (i == 0) ? 0 :  nKeptHits[i-1];
                    
                    for (uint32_t h = 0; h < (hitEnd - hitBegin); ++h) {
                        track_hits_out[outStart+h].id()    = track_hits[hitBegin + h].id();
                        track_hits_out[outStart+h].detId() = track_hits[hitBegin + h].detId();
                    }
                }
                else{
#ifdef KERNELS_DEBUG
                    printf("PixelTrackTorchHighPuritySelectorKernels: Error inputTrackIdx is negative");
#endif                    
                }
                
            } 
        }
    };

    struct HitOffsetCompactKernel{
        // TODO: implement the prefix scan provided my cmssw
        // NOTE:
        //  - originalTrackIndex[] is compacted in place
        //  - nKeptHits[] initially stores per-track hit counts
        //  - after this kernel, nKeptHits[] becomes an inclusive prefix sum

        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int  maxPreselectedTracks,
            int* originalTrackIndex,
            int* nKeptTracks,
            int* nKeptHits
        ) const 
        {
            int nTracks = 0;
            for (int j = 0; j < maxPreselectedTracks; j++){
                if (originalTrackIndex[j] != -1){
                    originalTrackIndex[nTracks]  = originalTrackIndex[j];
                    nKeptHits[nTracks] = nKeptHits[j];   
                    nTracks++;
                }
            }
            *nKeptTracks = nTracks;

            //loop over the tracks we kept and do the prefix sum
            for (int i = 1; i < nTracks; i++)
            {
                nKeptHits[i] += nKeptHits[i-1]; 
            }
            for (int i = nTracks; i < maxPreselectedTracks; i++){
                nKeptHits[i] = nKeptHits[nTracks-1];
                originalTrackIndex[i] = -1;
            }
        }
    };

    struct ScoreFilterKernel{
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            const int maxPreselectedTracks,
            const double scoreThreshold,
            int* originalTrackIndex,
            const PixelTrackScoresSoA::View trackScores,
            int* trackIndex,
            int* selectionMask,
            int* nKeptHits,
            int* nKeptHits_copy
        ) const 
        {
            for (auto i : alpaka::uniformElements(acc, maxPreselectedTracks)) {  
                nKeptHits_copy[i] = nKeptHits[i];
                trackIndex[i] = originalTrackIndex[i];  
                const float score = trackScores[i].score();  
                // Invalidate track slot if DNN score is below threshold
                bool isHP = score >= scoreThreshold && originalTrackIndex[i] != -1;
                selectionMask[i] = isHP ? 1 : 0; 
            }
        }
    };

    struct FilterArray{
        template <typename TAcc>
        ALPAKA_FN_ACC void operator()(
            TAcc const& acc,
            int* old_array,
            int* new_array,
            int* offsets,
            int old_size,
            int* new_size
        ) const 
        {
            auto worker = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
            const uint32_t new_size_local = offsets[old_size-1];
            if (worker == 0){
                *new_size = new_size_local;
            }
            for (auto i : alpaka::uniformElements(acc, old_size)){
                bool is_first = 
                    offsets[i] > 0 &&
                    ((i == 0) || (offsets[i] != offsets[i-1]));

                if (is_first){
                    new_array[offsets[i]-1] = old_array[i];
                }
            }
        }
    };

    void launchCAPreselectionKernel(
        Queue& queue,
        const int maxNumberOfTracks,
        const int minNumberOfHits,
        const ::pixelTrack::Quality minimumTrackQuality,
        const ::reco::TrackSoAConstView tracks,
        int* originalTrackIndex,
        int* preselectionOffsets,
        int* nKeptTracks)
    {
        auto trackIndexSoA    = cms::alpakatools::make_device_buffer<int[]>(queue, maxNumberOfTracks);
        auto preselectionMask = cms::alpakatools::make_device_buffer<int[]>(queue, maxNumberOfTracks);
        
        alpaka::memset(queue, trackIndexSoA, 0);
        alpaka::memset(queue, preselectionMask, 0);

        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxNumberOfTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);
        
        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            CAPreselectionKernel{},
            maxNumberOfTracks,
            minNumberOfHits,
            minimumTrackQuality,
            tracks,
            preselectionMask.data(),
            trackIndexSoA.data()
        );

        constexpr auto threadsPrefixScan = 1024;
        auto blocksPrefixScan = (maxNumberOfTracks + threadsPrefixScan - 1) / threadsPrefixScan;
        auto workDivPrefixScan = cms::alpakatools::make_workdiv<Acc1D>(blocksPrefixScan, threadsPrefixScan);
        auto bCounter = cms::alpakatools::make_device_buffer<int32_t>(queue);
        alpaka::memset(queue, bCounter, 0);

        alpaka::exec<Acc1D>(
            queue,
            workDivPrefixScan,
            cms::alpakatools::multiBlockPrefixScan<int>(),
            preselectionMask.data(),
            preselectionOffsets,
            maxNumberOfTracks,
            blocksPrefixScan,
            bCounter.data(),
            alpaka::getPreferredWarpSize(alpaka::getDev(queue))
        );

        alpaka::exec<Acc1D>(
            queue,
            workDivPrefixScan,
            FilterArray{},
            trackIndexSoA.data(),
            originalTrackIndex,
            preselectionOffsets,
            maxNumberOfTracks,
            nKeptTracks
        );

        /*
        validateSelectionsWithSorting(
            queue,
            maxPreselectedTracks,
            nKeptTracks,
            originalTrackIndex,
            originalTrackIndex_scan.data()
        );
        */
    }

    void launchFeaturesExtractorKernel(
        Queue& queue,
        const int maxPreselectedTracks,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        const ::reco::TrackingRecHitConstView hits,
        PixelTrackFeaturesSoAView trackFeatures,
        PixelRecHitFeaturesSoAView hitFeatures,
        const int* nKeptTracks,
        int* nKeptHits,
        int* originalTrackIndex)
    {
        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxPreselectedTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            FeaturesExtractorKernel{},
            maxPreselectedTracks,
            tracks,
            track_hits,
            hits,
            trackFeatures,
            hitFeatures,
            nKeptTracks,
            nKeptHits,
            originalTrackIndex
        );
    }

    void launchScoreFilterKernel(
        Queue& queue,
        const int maxPreselectedTracks,
        const double scoreThreshold,
        int* originalTrackIndex,
        int* nKeptTracks,
        int* nKeptHits,
        const PixelTrackScoresSoA::View trackScores)
    {
        // Produce a selection mask out of the DNN scores
        auto trackIndex    = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks);
        auto selectionMask = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks);
        auto selectionOffsets = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks);
        auto nKeptHits_copy = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks);

        alpaka::memset(queue, trackIndex, 0);
        alpaka::memset(queue, selectionMask, 0);
        alpaka::memset(queue, selectionOffsets, 0);
        alpaka::memset(queue, nKeptHits_copy, 0);

        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxPreselectedTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            ScoreFilterKernel{},
            maxPreselectedTracks,
            scoreThreshold,
            originalTrackIndex,
            trackScores,
            trackIndex.data(),
            selectionMask.data(),
            nKeptHits,
            nKeptHits_copy.data()
        );
        
        // Apply the selection mask to compact the originalTrackIndex array and produce the final list of selected tracks, while also counting the number of kept tracks
        constexpr auto threadsPrefixScan = 1024;
        auto blocksPrefixScan = (maxPreselectedTracks + threadsPrefixScan - 1) / threadsPrefixScan;
        auto workDivPrefixScan = cms::alpakatools::make_workdiv<Acc1D>(blocksPrefixScan, threadsPrefixScan);
        auto bCounter = cms::alpakatools::make_device_buffer<int32_t>(queue);
        alpaka::memset(queue, bCounter, 0);

        alpaka::exec<Acc1D>(
            queue,
            workDivPrefixScan,
            cms::alpakatools::multiBlockPrefixScan<int>(),
            selectionMask.data(),
            selectionOffsets.data(),
            maxPreselectedTracks,
            blocksPrefixScan,
            bCounter.data(),
            alpaka::getPreferredWarpSize(alpaka::getDev(queue))
        );

        alpaka::exec<Acc1D>(
            queue,
            workDivPrefixScan,
            FilterArray{},
            trackIndex.data(),
            originalTrackIndex,
            selectionOffsets.data(),
            maxPreselectedTracks,
            nKeptTracks
        );

        alpaka::exec<Acc1D>(
            queue,
            workDivPrefixScan,
            FilterArray{},
            nKeptHits,
            nKeptHits_copy.data(),
            selectionOffsets.data(),
            maxPreselectedTracks,
            nKeptTracks
        );

        alpaka::memset(queue, bCounter, 0);

        alpaka::exec<Acc1D>(
            queue,
            workDivPrefixScan,
            cms::alpakatools::multiBlockPrefixScan<int>(),
            nKeptHits_copy.data(),
            nKeptHits,
            maxPreselectedTracks,
            blocksPrefixScan,
            bCounter.data(),
            alpaka::getPreferredWarpSize(alpaka::getDev(queue))
        );
    }

    void launchHitOffsetCompactKernel(
        Queue& queue,
        const int  maxPreselectedTracks,
        int* originalTrackIndex,
        int* nKeptTracks,
        int* nKeptHits)
    {
        constexpr uint32_t threadsPerBlock = 1;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(1, 1);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);
        
        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            HitOffsetCompactKernel{},
            maxPreselectedTracks,
            originalTrackIndex,
            nKeptTracks,
            nKeptHits
        );
    }

    reco::TracksSoACollection launchProduceOutputTracks(
        Queue& queue,
        const int maxPreselectedTracks,
        const int avgHitsPerTrack,
        const ::reco::TrackSoAConstView tracks,
        const ::reco::TrackHitSoAConstView track_hits,
        int* originalTrackIndex,
        const int* nKeptTracks,
        const int* nKeptHits)
    {
        reco::TracksSoACollection tracks_out({{int(maxPreselectedTracks), int(maxPreselectedTracks * avgHitsPerTrack)}}, queue);
        
        constexpr uint32_t threadsPerBlock = 256;
        const uint32_t blocks =
            cms::alpakatools::divide_up_by(maxPreselectedTracks, threadsPerBlock);
        const auto workDiv =
            cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

        alpaka::exec<Acc1D>(
            queue,
            workDiv,
            PixelTrackFilterKernel{},
            maxPreselectedTracks,
            tracks,
            track_hits,
            originalTrackIndex,
            nKeptTracks,
            nKeptHits,
            tracks_out.view(),
            tracks_out.view<TrackHitSoA>()
        );

        return tracks_out;
    }
}