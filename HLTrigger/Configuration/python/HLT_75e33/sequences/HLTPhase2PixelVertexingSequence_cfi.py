import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone

from ..modules.hltPhase2PixelVertices_cfi import *

HLTPhase2PixelVertexingSequence = cms.Sequence(
    hltPhase2PixelVertices
)

from ..modules.hltVegaVertices_cfi import *
from ..modules.hltVegaVerticesSoA_cfi import *

HLTPhase2PixelVertexingSequenceVega = cms.Sequence(
    hltPhase2PixelVertices
    +hltVegaVerticesSoA
    +hltVegaVertices
)

# Adding Vega vertices for performance measurements
from Configuration.ProcessModifiers.vegaVertices_cff import vegaVertices
vegaVertices.toReplaceWith(HLTPhase2PixelVertexingSequence,
    HLTPhase2PixelVertexingSequenceVega
)

# Empty sequence as a placeholder to be filled when alpakaValidationHLT is active
HLTPhase2PixelVertexingSequenceSerialSync = cms.Sequence()

hltPhase2PixelVerticesSerialSync = hltPhase2PixelVertices.clone(
    TrackCollection = "hltPhase2PixelTracksSerialSync"
)
# Serial sequence for CPU vs. GPU validation, to be kept in sync with default sequence
from Configuration.ProcessModifiers.alpakaValidationHLT_cff import alpakaValidationHLT
alpakaValidationHLT.toReplaceWith(HLTPhase2PixelVertexingSequenceSerialSync,
    cms.Sequence(hltPhase2PixelVerticesSerialSync)
)
