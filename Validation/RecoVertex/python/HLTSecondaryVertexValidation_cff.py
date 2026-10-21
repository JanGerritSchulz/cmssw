import FWCore.ParameterSet.Config as cms

from Validation.RecoVertex.secondaryVertexAnalyzer_cfi import *
from Validation.RecoVertex.secondaryVertexAnalyzerCPC_cfi import *
from Validation.RecoVertex.associators_cff import hltSVAssociationsTask

hltSecondaryVertexValidatorCPC = secondaryVertexAnalyzerCPC.clone(
    rootFolder = 'HLT/SecondaryVertices/Validation',
    verbose = False,
    ignoreMissingCollections = True,
    doGenericSimPlots = True,
    doPerPdgPlots = True,
    recoVertexCollections = ["hltDeepInclusiveMergedVerticesPF"],
    vertexAssociators = ["hltSVAssociatorByPositionAndTracks4GeneralTracks"],
    primaryVertices = 'hltOfflinePrimaryVertices',
    hepMCProduct = 'generatorSmeared',
    simVertices = cms.InputTag('mix', 'MergedTrackTruth'),
    trackAssociation = 'tpToHLTGeneralTrackAssociation',
    minDecayLength = 0.01,
    maxDecayLength = 20.,
    minPt = 10.,
    minReconstructableDaughters = 2,
    minPtReconstructableDaughters = 0.9,
    signalPdgIds = [],
    bHadrons = True,
    cHadrons = True,
    sHadrons = True,
    taus = True,
    otherParticles = False,
)

hltSecondaryVertexValidator = secondaryVertexAnalyzer.clone(
    rootFolder = 'HLT/SecondaryVertices/Validation',
    verbose = False,
    ignoreMissingCollections = True,
    doGenericSimPlots = True,
    doPerPdgPlots = True,
    recoVertexCollections = ["hltDeepInclusiveMergedVertices", "hltVegaVertices"],
    vertexAssociators = ["hltSVAssociatorByPositionAndTracks4PixelTracks", "hltSVAssociatorByPositionAndTracks4PixelTracks"],
    primaryVertices = 'hltOfflinePrimaryVertices',
    hepMCProduct = 'generatorSmeared',
    simVertices = cms.InputTag('mix', 'MergedTrackTruth'),
    trackAssociation = 'tpToHLTpixelTrackAssociation',
    minDecayLength = 0.01,
    maxDecayLength = 20.,
    minPt = 10.,
    minReconstructableDaughters = 2,
    minPtReconstructableDaughters = 0.9,
    signalPdgIds = [],
    bHadrons = True,
    cHadrons = True,
    sHadrons = True,
    taus = True,
    otherParticles = False,
)

HLTSecondaryVertexValidation = cms.Sequence(
                                    hltSecondaryVertexValidator
                                    +hltSecondaryVertexValidatorCPC,
                                    hltSVAssociationsTask
                                    )
