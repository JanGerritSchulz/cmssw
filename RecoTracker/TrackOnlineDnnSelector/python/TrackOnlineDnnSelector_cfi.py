import FWCore.ParameterSet.Config as cms

TrackOnlineDnnSelector = cms.EDProducer("trackOnlineDnnSelector",
    maxRecHits = cms.uint32(16),
    probThreshold = cms.double(0.095),
    tracksSrc = cms.InputTag("hltPhase2PixelTracksCAExtension"),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)
