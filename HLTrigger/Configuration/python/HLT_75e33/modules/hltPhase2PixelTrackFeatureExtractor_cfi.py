import FWCore.ParameterSet.Config as cms

hltPhase2PixelTrackFeatureExtractor = cms.EDProducer('PixelTrackFeaturesExtractor@alpaka',
    pixelRecHitSrc = cms.InputTag('hltPhase2PixelRecHitsExtendedSoA'),
    pixelTrackSrc = cms.InputTag('hltPhase2PixelTracksSoA'),
    maxTracks = cms.int32(100000),
    minNumberOfHits = cms.int32(0),
    minQuality = cms.string('tight'),
)
