import FWCore.ParameterSet.Config as cms

hltGenVertexProducer = cms.EDProducer("GenVertexCandidateProducer",
    genParticles = cms.InputTag("mergedGenParticles"),
    secondaryVertices = cms.InputTag("hltDeepInclusiveMergedVerticesPF"),
    pvSrc = cms.InputTag("hltOfflinePrimaryVertices"),
    nRequiredCommonTracks = cms.int32(2),        # number of tracks required to match the genDaughters
    dlenSigMin = cms.double(0.0),
    dR_max = cms.double(0.03),                   # dR between tracks and daughters to be considered matched
    relPt_max = cms.double(0.5)                  # dPt/pt between tracks and daughters to be considered matched
)

from PhysicsTools.SimpleVertexTable.mergedGenParticles import mergedGenParticles

def custom_GV_producer(process):
    print("Track collection is running")
    process.mergedGenParticles = mergedGenParticles
    process.hltGenVertexTable = hltGenVertexProducer
    process.NanoValTables.insert(-1, process.mergedGenParticles)
    process.NanoValTables.insert(-1, process.hltGenVertexTable)
    return process