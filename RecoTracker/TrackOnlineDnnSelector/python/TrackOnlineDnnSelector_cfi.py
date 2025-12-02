import FWCore.ParameterSet.Config as cms

TrackOnlineDnnSelector = cms.EDProducer("trackOnlineDnnSelector",
    skipNonExistingSrc = cms.bool(True),
    maxRecHits = cms.uint32(16),
    tracksSrc = cms.InputTag("hltPhase2PixelTracksCAExtension"),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)
