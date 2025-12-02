import FWCore.ParameterSet.Config as cms

hltTrackOnlineDnnSelector = cms.EDProducer("TrackOnlineDnnSelector",
    skipNonExistingSrc = cms.bool(True),
    maxRecHits = cms.uint32(16),
    tracksSrc = cms.InputTag("hltPhase2PixelTracksCAExtension"),
    beamSpot = cms.InputTag("hltOnlineBeamSpot")
)
