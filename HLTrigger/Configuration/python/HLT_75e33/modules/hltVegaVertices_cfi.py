import FWCore.ParameterSet.Config as cms

hltVegaVertices = cms.EDProducer('VegaVertexProducerFromSoAAlpaka',
    beamSpot = cms.InputTag('hltOnlineBeamSpot'),
    trackCollection = cms.InputTag('hltPhase2PixelTracks'),
    vertexSoACollection = cms.InputTag('hltVegaVerticesSoA'),
    minTrackWeight = cms.float(0.5),
    refitVerticesWithAVF = cms.bool(False),
    paramsAVF = cms.PSet(
      sigmacut = cms.double(3),
      tini = cms.double(256),
      ratio = cms.double(0.25)
    )
)
