import FWCore.ParameterSet.Config as cms

hltVegaVerticesSoA = cms.EDProducer('VegaVertexSoAProducer@alpaka',
    trackSrc = cms.InputTag('hltPhase2PixelTrackTorchHighPuritySelector'),
    maxVertices = cms.int32(100),
    pairParams = cms.PSet(
      minPt = cms.float(0.9),
      maxDPhi = cms.float(0.4),
      maxDEta = cms.float(0.4),
      maxDZ = cms.float(0.75),
      maxLinDistance = cms.float(0.1),
      max3DDistance = cms.float(0.02)
    ),
    tripletParams = cms.PSet(
      max3DDistance = cms.float(0.02)
    ),
    alpaka = cms.untracked.PSet(
      backend = cms.untracked.string('')
    )
)