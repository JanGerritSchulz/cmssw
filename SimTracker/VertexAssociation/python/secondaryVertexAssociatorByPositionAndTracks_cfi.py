import FWCore.ParameterSet.Config as cms

# for SV of type reco::Vertex
secondaryVertexAssociatorByPositionAndTracks = cms.EDProducer("VertexAssociatorByPositionAndTracksProducer",
    absZ = cms.double(1.0),
    sigmaX = cms.double(5),
    sigmaY = cms.double(5),
    sigmaZ = cms.double(5),
    maxRecoZ = cms.double(1000),
    absT = cms.double(-1),
    sigmaT = cms.double(-1),
    maxRecoT = cms.double(-1),
    sharedTrackFraction = cms.double(2),
    filterSimVerticesForPVs = cms.bool(False),
    weightMethod = cms.string('nSharedTracks'),
    trackAssociation = cms.InputTag('trackingParticleRecoTrackAsssociation')
)

# for SV of type reco::VertexCompositeCandidate
secondaryVertexAssociatorByPositionAndTracksCPC = cms.EDProducer("VertexAssociatorByPositionAndTracksProducerCPC",
    absZ = cms.double(1.0),
    sigmaX = cms.double(5),
    sigmaY = cms.double(5),
    sigmaZ = cms.double(5),
    maxRecoZ = cms.double(1000),
    absT = cms.double(-1),
    sigmaT = cms.double(-1),
    maxRecoT = cms.double(-1),
    sharedTrackFraction = cms.double(2),
    filterSimVerticesForPVs = cms.bool(False),
    weightMethod = cms.string('nSharedTracks'),
    trackAssociation = cms.InputTag('trackingParticleRecoTrackAsssociation')
)
