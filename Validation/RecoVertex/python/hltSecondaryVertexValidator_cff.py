import FWCore.ParameterSet.Config as cms

from Validation.RecoVertex.secondaryVertexAnalyzerCPC_cfi import *
from Validation.RecoVertex.associators_cff import hltSVAssociationsTask

hltSecondaryVertexValidator = secondaryVertexAnalyzerCPC.clone(
    rootFolder = cms.untracked.string('HLT/SecondaryVertices/Validation'),
    verbose = cms.untracked.bool(False),
    doGenericSimPlots = cms.untracked.bool(True),
    doPerPdgPlots = cms.untracked.bool(True),
    recoVertexCollections = cms.VInputTag(["hltDeepInclusiveMergedVerticesPF"]),
    vertexAssociators = cms.VInputTag(["hltSVAssociatorByPositionAndTracks4GeneralTracks"]),
    primaryVertices = cms.InputTag('hltOfflinePrimaryVertices'),
    hepMCProduct = cms.InputTag('generatorSmeared'),
    simVertices = cms.InputTag('mix', 'MergedTrackTruth'),
    trackAssociation = cms.InputTag('tpToHLTGeneralTrackAssociation'),
    minDecayLength = cms.double(0.01),
    maxDecayLength = cms.double(20.),
    minReconstructableDaughters = cms.int32(2),
    minPtReconstructableDaughters = cms.double(0.5),
    signalPdgIds = cms.vint32(),
    bHadrons = cms.bool(True),
    cHadrons = cms.bool(True),
    sHadrons = cms.bool(True),
    taus = cms.bool(True),
    otherParticles = cms.bool(False),
)

hltSecondaryVertexValidation = cms.Sequence(
                                    hltSecondaryVertexValidator,
                                    hltSVAssociationsTask
                                    )
