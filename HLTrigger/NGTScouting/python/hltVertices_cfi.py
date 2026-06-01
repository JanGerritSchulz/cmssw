import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleSecondaryVertexFlatTableProducer_cfi import simpleSecondaryVertexFlatTableProducer

hltVertexTable = cms.EDProducer("HLTVertexTableProducer",
                                skipNonExistingSrc = cms.bool(True),
                                pvSrc = cms.InputTag("hltOfflinePrimaryVertices"),
                                svSrc = cms.InputTag("hltDeepInclusiveMergedVerticesPF"),
                                goodPvCut = cms.string("!isFake && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"), 
                                goodSvCut = cms.string(""), 
                                pfSrc = cms.InputTag("hltParticleFlowTmp"),
                                doSVs = cms.bool(True),
                                dlenMin = cms.double(0),
                                dlenSigMin = cms.double(3),
                                pvName = cms.string("hltPrimaryVertex"),
                                svName = cms.string("hltSecondaryVertex"),
                                svDoc  = cms.string("secondary vertices from IVF algorithm"))

hltPixelVertexTable = cms.EDProducer("HLTVertexTableProducer",
                                     skipNonExistingSrc = cms.bool(True),
                                     pvSrc = cms.InputTag("hltPhase2PixelVertices"),
                                     goodPvCut = cms.string(""),
                                     goodSvCut = cms.string(""), 
                                     usePF = cms.bool(False), # use directly the tracks from PV fit 
                                     pfSrc = cms.InputTag(""),
                                     doSVs = cms.bool(False),
                                     dlenMin = cms.double(0),
                                     dlenSigMin = cms.double(3),
                                     pvName = cms.string("hltPixelVertex"))
