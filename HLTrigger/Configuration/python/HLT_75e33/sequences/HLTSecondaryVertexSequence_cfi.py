import FWCore.ParameterSet.Config as cms

from ..modules.hltDeepInclusiveMergedVertices_cfi import *
from ..modules.hltDeepInclusiveSecondaryVertices_cfi import *
from ..modules.hltDeepInclusiveVertexFinder_cfi import *
from ..modules.hltDeepTrackVertexArbitrator_cfi import *

HLTSecondaryVertexSequence = cms.Sequence(
    hltDeepInclusiveVertexFinder
    +hltDeepInclusiveSecondaryVertices
    +hltDeepTrackVertexArbitrator
    +hltDeepInclusiveMergedVertices
)
