import FWCore.ParameterSet.Config as cms

hltVegaVerticesSoA = cms.EDProducer('VegaVertexSoAProducer@alpaka',
    trackSoACollection = cms.InputTag('hltPhase2PixelTrackTorchHighPuritySelector'),
    maxVertices = cms.int32(4000),
    pairParams = cms.PSet(
      minTrackQuality = cms.string("tight"),
      minPt = cms.float(0.9),
      maxDPhi = cms.float(0.4),
      maxDEta = cms.float(0.4),
      maxDZ = cms.float(0.75),
      maxLinDistance = cms.float(0.1),
      maxChi2 = cms.float(5.0),
      minDXYSig = cms.float(2.0),
      maxDXYSig = cms.float(-1)
    ),
    tripletParams = cms.PSet(
      max3DDistance = cms.float(0.02)
    ),
    vertexParams = cms.PSet(
      maxChi2 = cms.float(5.0)
    ),
    alpaka = cms.untracked.PSet(
      backend = cms.untracked.string('')
    )
)