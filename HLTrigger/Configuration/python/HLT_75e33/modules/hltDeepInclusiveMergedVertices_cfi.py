import FWCore.ParameterSet.Config as cms

hltDeepInclusiveMergedVertices = cms.EDProducer("VertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("hltDeepTrackVertexArbitrator")
)
