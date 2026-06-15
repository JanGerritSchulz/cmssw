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
    simVertices = cms.InputTag('mix', 'MergedTrackTruth'),
    trackAssociation = cms.InputTag('trackingParticleRecoTrackAsssociation'),
    minDecayLength = cms.double(0.01),
    minReconstructableDaughters = cms.int32(2),
    absEtaMax = cms.double(2.5),
    signalPdgIds = cms.vint32(),
)

hltSecondaryVertexValidation = cms.Sequence(
                                    hltSecondaryVertexValidator,
                                    hltSVAssociationsTask
                                    )
