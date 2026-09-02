import FWCore.ParameterSet.Config as cms

process = cms.Process("WRITE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("EmptySource")
process.maxEvents.input = 5

process.zvertexSoA = cms.EDProducer("TestWriteHostZVertexSoA",
    vertexSize = cms.uint32(2708)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(sys.argv[1])
)

process.path = cms.Path(process.zvertexSoA)
process.endPath = cms.EndPath(process.out)
