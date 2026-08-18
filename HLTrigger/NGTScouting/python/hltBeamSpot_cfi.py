import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleBeamspotFlatTableProducer_cfi import simpleBeamspotFlatTableProducer

hltBeamSpotTable = simpleBeamspotFlatTableProducer.clone(
    src = cms.InputTag("hltOnlineBeamSpot"),
    name = cms.string("hltBeamSpot"),
    doc = cms.string("hltOnlineBeamSpot, the hlt reconstructed beamspot"),
    variables = cms.PSet(
       type = Var("type()","int16",doc="BeamSpot type (Unknown = -1, Fake = 0, LHC = 1, Tracker = 2)"),
       z = Var("position().z()",float,doc="BeamSpot center, z coordinate (cm)",precision=-1),
       zError = Var("z0Error()",float,doc="Error on BeamSpot center, z coordinate (cm)",precision=-1),
       sigmaZ = Var("sigmaZ()",float,doc="Width of BeamSpot in z (cm)",precision=-1),
       sigmaZError = Var("sigmaZ0Error()",float,doc="Error on width of BeamSpot in z (cm)",precision=-1),
    ),
)
