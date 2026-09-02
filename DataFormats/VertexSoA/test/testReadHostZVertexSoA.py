import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("READ")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:"+sys.argv[1]))

process.testReadHostZVertexSoA = cms.EDAnalyzer("TestReadHostZVertexSoA",
    input = cms.InputTag("zvertexSoA", "", "WRITE")
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('testZVertexSoAReader.root'),
    fastCloning = cms.untracked.bool(False)
)

process.path = cms.Path(process.testReadHostZVertexSoA)

process.endPath = cms.EndPath(process.out)

