import FWCore.ParameterSet.Config as cms

hltESPChi2MeasurementEstimatorForRefit = cms.ESProducer(
    "Chi2MeasurementEstimatorESProducer",
    ComponentName    = cms.string('hltESPChi2MeasurementEstimatorForRefit'),
    MaxChi2          = cms.double(100000.),   # effectively open — accept all hits
    nSigma           = cms.double(100.),      # same
    MaxDisplacement  = cms.double(100.),      # loose geometric guard
    MaxSagitta       = cms.double(-1.),       # disabled
    MinimalTolerance = cms.double(10.),       # loose
    MinPtForHitRecoveryInGluedDet = cms.double(1e12),  # effectively disabled
)
