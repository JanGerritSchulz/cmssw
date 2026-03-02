import FWCore.ParameterSet.Config as cms


# list of layers to exclude from the CA (empty list doesn't exclude; [28, 29, 30] excludes the OT)
layersToExclude = []

# layers
layers = [
    #     0,        1,     2,       3,      4,     5,      6,       7
    # index, isBarrel, caDCA, caTheta, startR, DCurv, DCurv0, fishCut
    [     0,     True,  0.15,   0.002,   99.0,  99.0,   99.0, 0.99999],
    [     1,     True,  0.25,   0.002,   99.0,  99.0,   99.0, 0.99999],
    [     2,     True,  0.20,   0.002,   99.0,  99.0,   99.0, 0.99999],
    [     3,     True,  0.20,   0.002,   99.0,  99.0,   99.0, 0.99999],
    [     4,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [     5,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [     6,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [     7,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [     8,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [     9,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    10,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    11,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    12,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    13,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    14,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    15,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    16,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    17,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    18,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    19,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    20,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    21,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    22,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    23,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    24,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    25,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    26,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    27,    False,  0.25,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    28,     True,  0.10,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    29,     True,  0.10,   0.003,   99.0,  99.0,   99.0, 0.99999],
    [    30,     True,  0.10,   0.003,   99.0,  99.0,   99.0, 0.99999],
]

# layerPairs for doublet building including pair-specific cut values
layerPairs = [
    #  0,  1,     2,     3,      4,      5,      6,       7,       8,     9,     10,     11,     12
    #  i,  o, start,  skip, phiCut,  minIn,  maxIn,  minOut,  maxOut, maxDR,  minDZ,  maxDZ, ptCuts
    [  0,  1,  True, False,    350,  -17.0,   17.0,  -10000,   10000,   5.0,  -16.0,   16.0,  0.85],
    [  0,  2,  True, False,    600,  -14.0,   14.0,  -10000,   10000,  10.0,  -16.0,   16.0,  0.85],
    [  0,  4,  True, False,    450,    4.0,  10000,       0,    10.0,   8.0,    0.0,   25.0,  0.85],
    [  0,  5,  True, False,    522,    7.0,  10000,       0,   10000,   5.0,    0.0,   25.0,  0.85],
  # [  0,  6, False, False,    522,   11.0,  10000,       0,   10000,   5.0, -10000,  10000,  0.85],
    [  0, 16,  True, False,    450, -10000,   -4.0,       0,    10.0,   8.0,  -25.0,    0.0,  0.85],
    [  0, 17,  True, False,    522, -10000,   -7.0,       0,   10000,   5.0,  -25.0,    0.0,  0.85],
  # [  0, 18, False, False,    522, -10000,  -10.0,       0,   10000,   5.0, -10000,  10000,  0.85],
    [  1,  2,  True, False,    400,  -17.0,   17.0,  -10000,   10000,   7.0,  -13.0,   13.0,  0.85],
    [  1,  3, False, False,    650,  -15.0,   15.0,  -10000,   10000,  10.0,  -15.0,   15.0,  0.85],
    [  1,  4,  True, False,    500,    6.0,  10000,     6.5,   10000,   8.0,    0.0,   19.0,  0.85],
    [  1,  5, False, False,    730,    9.0,  10000,     6.5,   10000,  10.0,    0.0,   21.0,  0.85],
  # [  1,  6, False, False,    730,   13.0,  10000,     6.5,   10000,   8.0, -10000,  10000,  0.85],
    [  1, 16,  True, False,    500, -10000,   -6.0,     6.5,   10000,   8.0,  -19.0,    0.0,  0.85],
    [  1, 17, False, False,    730, -10000,   -9.0,     6.5,   10000,  10.0,  -21.0,    0.0,  0.85],
  # [  1, 18, False, False,    730, -10000,  -13.0,     6.5,   10000,   8.0, -10000,  10000,  0.85],
  # [  1, 28, False, False,   1300,    7.0,  10000,    30.0,    40.0, 10000,   19.0,   32.0,   1.0],
  # [  1, 28, False, False,   1300, -10000,   -7.0,   -40.0,   -30.0, 10000,  -32.0,  -19.0,   1.0],
    [  2,  3,  True, False,    350,  -18.0,   18.0,  -10000,   10000,   7.0,   -9.0,    9.0,  0.85],
    [  2,  4, False, False,    400,   11.0,  10000,    11.7,   10000,   7.0,    0.0,   13.0,  0.85],
    [  2, 16, False, False,    400, -10000,  -11.0,    11.7,   10000,   7.0,  -13.0,    0.0,  0.85],
    [  2, 28, False, False,   1200,    -10,     10,   -30.0,    30.0, 10000,  -15.0,   15.0,   2.0],
    [  2, 28, False, False,   1200,    -20,    -10,   -50.0,   -25.0, 10000,  -35.0,  -10.0,  0.85],
    [  2, 28, False, False,   1200,     10,     20,    25.0,    50.0, 10000,   10.0,   35.0,  0.85],
  # [  2, 28, False, False,   1200,    -20,     20,   -50.0,    50.0, 10000,  -35.0,   35.0,  0.85],
    [  3, 28, False, False,   1000,    -20,     20,   -45.0,    45.0, 10000,  -22.0,   22.0,  0.85],
  # [  3, 29, False, False,   1500,    -40,     40,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [  4,  5,  True, False,    300,      0,   14.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [  4,  6, False, False,    522,      0,   14.0,     3.5,   10000,   9.0, -10000,  10000,  0.85],
    [  4, 28, False, False,   1000,   11.6,  10000,    30.0,    57.5,  16.0,    5.0,   32.5,  0.85],
  # [  4, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [  5,  6,  True, False,    300,      0,   13.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [  5,  7, False, False,    522,      0,   13.0,     3.5,   10000,   9.0, -10000,  10000,  0.85],
    [  5, 28, False, False,   1000,   11.6,  10000,    40.0,    80.0,  16.0,  -10.0,   50.0,  0.85],
  # [  5, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [  6,  7,  True, False,    250,      0,   13.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [  6,  8, False, False,    522,      0,   13.0,     3.5,   10000,   9.0, -10000,  10000,  0.85],
    [  6, 28, False, False,   1000,   11.6,  10000,    55.0,    95.0,  16.0,    5.0,   50.0,  0.85],
  # [  6, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [  7,  8,  True, False,    250,      0,   13.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [  7,  9, False, False,    522,      0,   13.0,     3.5,   10000,   8.0, -10000,  10000,  0.85],
    [  7, 28, False, False,   1000,   11.8,  10000,    70.0,   110.0,  16.0,   15.0,   70.0,  0.85],
  # [  7, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [  8,  9,  True, False,    250,      0,   13.0,     3.5,   10000,   4.0, -10000,  10000,  0.85],
    [  8, 10, False, False,    522,      0,   13.0,     3.5,   10000,   8.0, -10000,  10000,  0.85],
    [  8, 28, False, False,    850,      0,  10000,    80.0,   10000,  14.0,   25.0,   70.0,  0.85],
  # [  8, 29, False, False,   1000,      0,  10000,       0,   10000, 10000, -10000,  10000,  0.85],
    [  9, 10,  True, False,    300,      0,   13.0,     4.0,   10000,   4.5, -10000,  10000,  0.85],
    [  9, 11, False, False,    522,      0,   13.0,     4.0,   10000,   8.0, -10000,  10000,  0.85],
  # [  9, 28, False, False,   1000,      0,  10000,       0,   10000, 10000, -10000,  10000,  0.85],
    [ 10, 11,  True, False,    240,      0,   13.0,     3.5,   10000,   4.0, -10000,  10000,  0.85],
    [ 10, 12, False, False,    650,   12.5,   16.5,    20.0,   10000,  10.0, -10000,  10000,  0.85],
    [ 11, 12, False, False,    300,      0,   16.5,     6.0,    21.0,   5.0, -10000,  10000,  0.85],
    [ 11, 13, False, False,    200,      0,    6.0,       0,     7.5,   3.0, -10000,  10000,  0.85],
    [ 11, 14, False, False,    220,      0,    4.6,       0,     7.5,   3.0, -10000,  10000,  0.85],
    [ 11, 15, False, False,    250,      0,    6.0,       0,   10000,   4.0, -10000,  10000,  0.85],
    [ 12, 13, False, False,    250,      0,   22.5,     7.0,   10000,   4.0, -10000,  10000,  0.85],
    [ 13, 14, False, False,    250,      0,   22.5,     7.0,   10000,   4.0, -10000,  10000,  0.85],
    [ 14, 15, False, False,    250,      0,   22.5,     7.0,   10000,   3.5, -10000,  10000,  0.85],
    [ 16, 17,  True, False,    300,      0,   14.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [ 16, 18, False, False,    522,      0,   14.0,     3.5,   10000,   9.0, -10000,  10000,  0.85],
    [ 16, 28, False, False,   1000,   11.6,  10000,   -57.5,   -30.0,  16.0,  -32.5,   -5.0,  0.85],
  # [ 16, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 17, 18,  True, False,    300,      0,   13.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [ 17, 19, False, False,    522,      0,   13.0,     3.5,   10000,   9.0, -10000,  10000,  0.85],
    [ 17, 28, False, False,   1000,   11.6,  10000,   -70.0,   -40.0,  16.0,  -50.0,  -10.0,  0.85],
  # [ 17, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 18, 19,  True, False,    250,      0,   13.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [ 18, 20, False, False,    522,      0,   13.0,     3.5,   10000,   9.0, -10000,  10000,  0.85],
    [ 18, 28, False, False,   1000,   11.6,  10000,   -95.0,   -55.0,  16.0,  -50.0,   -5.0,  0.85],
  # [ 18, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 19, 20,  True, False,    250,      0,   13.0,     3.5,   10000,   4.5, -10000,  10000,  0.85],
    [ 19, 21, False, False,    522,      0,   13.0,     3.5,   10000,   8.0, -10000,  10000,  0.85],
    [ 19, 28, False, False,   1000,   11.8,  10000,  -110.0,   -70.0,  16.0,  -70.0,  -15.0,  0.85],
  # [ 19, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 20, 21,  True, False,    250,      0,   13.0,     3.5,   10000,   4.0, -10000,  10000,  0.85],
    [ 20, 22, False, False,    522,      0,   13.0,     3.5,   10000,   8.0, -10000,  10000,  0.85],
    [ 20, 28, False, False,   1000,      0,  10000,  -10000,   -80.0,  14.0,  -70.0,  -25.0,  0.85],
  # [ 20, 29, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 21, 22,  True, False,    300,      0,   13.0,     4.0,   10000,   4.5, -10000,  10000,  0.85],
    [ 21, 23, False, False,    522,      0,   13.0,     4.0,   10000,   8.0, -10000,  10000,  0.85],
  # [ 21, 28, False, False,   1000,      0,  10000,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 22, 23,  True, False,    240,      0,   13.0,     3.5,   10000,   4.0, -10000,  10000,  0.85],
    [ 22, 24, False, False,    650,   12.5,   16.5,    20.0,   10000,  10.0, -10000,  10000,  0.85],
    [ 23, 24, False, False,    300,      0,   16.5,     6.0,    21.0,   5.0, -10000,  10000,  0.85],
    [ 23, 25, False, False,    200,      0,    6.0,       0,     7.5,   3.0, -10000,  10000,  0.85],
    [ 23, 26, False, False,    220,      0,    4.6,       0,     7.5,   3.0, -10000,  10000,  0.85],
    [ 23, 27, False, False,    250,      0,    6.0,       0,   10000,   4.0, -10000,  10000,  0.85],
    [ 24, 25, False, False,    250,      0,   22.5,     7.0,   10000,   4.0, -10000,  10000,  0.85],
    [ 25, 26, False, False,    250,      0,   22.5,     7.0,   10000,   4.0, -10000,  10000,  0.85],
    [ 26, 27, False, False,    250,      0,   22.5,     7.0,   10000,   3.5, -10000,  10000,  0.85],
    [ 28, 29, False, False,   1100,  -1200,   1200,  -10000,   10000, 10000,  -50.0,   50.0,  0.85],
  # [ 28, 30, False, False,   2000,    -40,     40,  -10000,   10000, 10000, -10000,  10000,  0.85],
    [ 29, 30, False, False,   1250,  -1200,   1200,  -10000,   10000, 10000,  -40.0,   40.0,  0.85],
]

# find the layerPairs that contain a layer that is excluded
excludeLayerPair = [any([(lp[0] == l) or (lp[1] == l) for l in layersToExclude]) for lp in layerPairs]
excludeCAExtension = [any([(lp[0] == l) or (lp[1] == l) for l in [28, 29, 30]]) for lp in layerPairs]

# exclude those layerPairs
layerPairsAlpaka = []
layerPairsCAExtension = []
for i, lp in enumerate(layerPairs):
    if (not excludeLayerPair[i]) and (not excludeCAExtension[i]):
        layerPairsAlpaka.append(lp)
    if not excludeLayerPair[i]:
        layerPairsCAExtension.append(lp)

# get startingPairs for Ntuplet building
startingPairsAlpaka = []
for i, lp in enumerate(layerPairsAlpaka):
    if lp[2]:
        startingPairsAlpaka.append(i)

startingPairsCAExtension = []
for i, lp in enumerate(layerPairsCAExtension):
    if lp[2]:
        startingPairsCAExtension.append(i)

hltPhase2PixelTracksSoA = cms.EDProducer('CAHitNtupletAlpakaPhase2@alpaka',
    pixelRecHitSrc = cms.InputTag('hltPhase2SiPixelRecHitsSoA'),
    ptmin = cms.double(0.9),
    hardCurvCut = cms.double(0.01425), # corresponds to 800 MeV in 3.8T.
    earlyFishbone = cms.bool(True),
    lateFishbone = cms.bool(False),
    onlySameLayersFishbone = cms.bool(False),
    fillStatistics = cms.bool(False),
    minHitsPerNtuplet = cms.uint32(4),
    maxNumberOfDoublets = cms.string(str(6*512*1024)),
    maxNumberOfTuples = cms.string(str(60*1024)),
    cellZ0Cut = cms.double(12.5), # it's half the BS width! It has nothing to do with the sample!!
    minYsizeB1 = cms.int32(20),
    minYsizeB2 = cms.int32(18),
    maxDYsize12 = cms.int32(12),
    maxDYsize = cms.int32(10),
    maxDYPred = cms.int32(24),
    avgHitsPerTrack = cms.double(7.0),
    avgCellsPerHit = cms.double(12),
    avgCellsPerCell = cms.double(0.151),
    avgTracksPerCell = cms.double(0.040),
    minHitsForSharingCut = cms.uint32(10),
    fitNas4 = cms.bool(False),
    useRiemannFit = cms.bool(False),
    doSharedHitCut = cms.bool(True),
    dupPassThrough = cms.bool(False),
    useSimpleTripletCleaner = cms.bool(True),
    disableTripletCleaner = cms.bool(False),
    disableFastDuplicateRemover = cms.bool(False),
    disableEarlyDuplicateRemover = cms.bool(False),
    trackQualityCuts = cms.PSet(
        maxChi2TripletsOrQuadruplets = cms.double(5.0),
        maxChi2Quintuplets = cms.double(5.0),
        maxChi2 = cms.double(5.0),
        minPt   = cms.double(0.9),
        maxTip  = cms.double(0.3),
        maxZip  = cms.double(12),
    ),
    geometry = cms.PSet(
        caDCACuts   = cms.vdouble([l[2] for l in layers[:28]]),
        caThetaCuts = cms.vdouble([l[3] for l in layers[:28]]),
        startingPairs = cms.vuint32(startingPairsAlpaka),
        startMaxInnerR = cms.vdouble([l[4] for l in layers[:28]]),
        caDCurvCuts    = cms.vdouble([l[5] for l in layers[:28]]),
        caDCurv0       = cms.vdouble([l[6] for l in layers[:28]]),
        fishboneCuts   = cms.vdouble([l[7] for l in layers[:28]]),
        pairGraph   = cms.vuint32(sum([[lp[0], lp[1]] for lp in layerPairsAlpaka], [])),
        skipsLayers = cms.vuint32( [int(lp[ 3]) for lp in layerPairsAlpaka]),
        phiCuts     = cms.vint32( [lp[ 4] for lp in layerPairsAlpaka]),
        minInner    = cms.vdouble([lp[ 5] for lp in layerPairsAlpaka]),
        maxInner    = cms.vdouble([lp[ 6] for lp in layerPairsAlpaka]),
        minOuter    = cms.vdouble([lp[ 7] for lp in layerPairsAlpaka]),
        maxOuter    = cms.vdouble([lp[ 8] for lp in layerPairsAlpaka]),
        maxDR       = cms.vdouble([lp[ 9] for lp in layerPairsAlpaka]),
        minDZ       = cms.vdouble([lp[10] for lp in layerPairsAlpaka]),
        maxDZ       = cms.vdouble([lp[11] for lp in layerPairsAlpaka]),
        ptCuts      = cms.vdouble([lp[12] for lp in layerPairsAlpaka]),
  ),
    # autoselect the alpaka backend
    alpaka = cms.untracked.PSet(backend = cms.untracked.string(''))
)

_hltPhase2PixelTracksSoA = cms.EDProducer("CAHitNtupletAlpakaPhase2OT@alpaka",
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')
    ),
    avgCellsPerCell = cms.double(0.5),
    avgCellsPerHit = cms.double(20),
    avgHitsPerTrack = cms.double(8.0),
    avgTracksPerCell = cms.double(0.09),
    cellZ0Cut = cms.double(13),
    disableEarlyDuplicateRemover = cms.bool(False),
    disableFastDuplicateRemover = cms.bool(False),
    disableTripletCleaner = cms.bool(False),
    doSharedHitCut = cms.bool(True),
    dupPassThrough = cms.bool(False),
    earlyFishbone = cms.bool(True),
    fillStatistics = cms.bool(False),
    fitNas4 = cms.bool(False),
    geometry = cms.PSet(
        caDCACuts = cms.vdouble(
            0.15, 0.25, 0.2, 0.2, 0.25,
            0.25, 0.25, 0.25, 0.25, 0.25,
            0.25, 0.25, 0.4, 0.4, 0.25,
            0.25, 0.25, 0.25, 0.25, 0.25,
            0.25, 0.25, 0.25, 0.25, 0.4,
            0.4, 0.25, 0.25, 0.5, 0.1,
            0.1
        ),
        caDCurv0 = cms.vdouble(
            99.0, 99.0, 99.0, 0.000854, 0.000913,
            0.00167, 0.00365, 0.00416, 0.00487, 0.00404,
            0.00401, 0.00415, 0.00415, 0.00433, 0.00433,
            0.00433, 0.000913, 0.00167, 0.00365, 0.00416,
            0.00487, 0.00404, 0.00401, 0.00415, 0.00415,
            0.00433, 0.00433, 0.00433, 0.000832, 0.000427,
            0.000276
        ),
        caDCurvCuts = cms.vdouble(
            99.0, 99.0, 99.0, 0.049, 0.0668,
            0.0873, 0.102, 0.106, 0.0524, 0.118,
            0.112, 0.0879, 0.0879, 0.105, 0.105,
            0.105, 0.0668, 0.0873, 0.102, 0.106,
            0.0524, 0.118, 0.112, 0.0879, 0.0879,
            0.105, 0.105, 0.105, 0.0816, 0.0482,
            0.0379
        ),
        caThetaCuts = cms.vdouble(
            0.002, 0.002, 0.002, 0.002, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.003, 0.003,
            0.003, 0.003, 0.003, 0.005, 0.005,
            0.003
        ),
        fishboneCuts = cms.vdouble(
            0.99999, 0.99999, 0.99999, 0.99999, 0.99999,
            0.99999, 0.99999, 0.99999, 0.999999, 0.999999,
            0.999999, 0.999999, 0.999999, 0.999999, 0.999999,
            0.999999, 0.99999, 0.99999, 0.99999, 0.99999,
            0.999999, 0.999999, 0.999999, 0.999999, 0.999999,
            0.999999, 0.999999, 0.999999, 0.99999, 0.99999,
            0.99999
        ),
        maxDR = cms.vdouble(
            5.0, 10.0, 8.5, 5.0, 8.5,
            5.0, 7.0, 10.0, 8.5, 10.0,
            8.5, 10.0, 7.0, 7.0, 7.0,
            10000, 10000, 10000, 10000, 4.5,
            9.0, 16.0, 4.5, 9.0, 16.0,
            4.5, 9.0, 16.0, 4.5, 8.0,
            16.0, 4.0, 8.0, 14.0, 4.5,
            8.0, 4.0, 10.0, 5.0, 3.0,
            3.0, 4.0, 4.0, 4.0, 3.5,
            4.5, 9.0, 16.0, 4.5, 9.0,
            16.0, 4.5, 9.0, 16.0, 4.5,
            8.0, 16.0, 4.0, 8.0, 14.0,
            4.5, 8.0, 4.0, 10.0, 5.0,
            3.0, 3.0, 4.0, 4.0, 4.0,
            3.5, 10000, 10000
        ),
        maxDZ = cms.vdouble(
            16.0, 16.0, 25.0, 25.0, 0.0,
            0.0, 13.0, 17.0, 19.0, 21.0,
            0.0, 0.0, 9.0, 13.0, 0.0,
            15.0, -10.0, 35.0, 22.0, 10000,
            10000, 32.5, 10000, 10000, 50.0,
            10000, 10000, 50.0, 10000, 10000,
            70.0, 10000, 10000, 70.0, 10000,
            10000, 10000, 10000, 10000, 10000,
            10000, 10000, 10000, 10000, 10000,
            10000, 10000, -5.0, 10000, 10000,
            -10.0, 10000, 10000, -5.0, 10000,
            10000, -15.0, 10000, 10000, -25.0,
            10000, 10000, 10000, 10000, 10000,
            10000, 10000, 10000, 10000, 10000,
            10000, 50.0, 40.0
        ),
        maxInner = cms.vdouble(
            17.0, 14.0, 10000, 10000, -3.0,
            -7.0, 17.0, 15.0, 10000, 10000,
            -6.0, -9.0, 18.0, 10000, -11.0,
            10, -10, 20, 20, 14.0,
            14.0, 10000, 13.0, 13.0, 10000,
            13.0, 13.0, 10000, 13.0, 13.0,
            10000, 13.0, 13.0, 10000, 13.0,
            13.0, 13.0, 16.5, 16.5, 6.0,
            4.6, 6.0, 22.5, 22.5, 22.5,
            14.0, 14.0, 10000, 13.0, 13.0,
            10000, 13.0, 13.0, 10000, 13.0,
            13.0, 10000, 13.0, 13.0, 10000,
            13.0, 13.0, 13.0, 16.5, 16.5,
            6.0, 4.6, 6.0, 22.5, 22.5,
            22.5, 1200, 1200
        ),
        maxOuter = cms.vdouble(
            10000, 10000, 12.0, 10000, 12.0,
            10000, 10000, 10000, 10000, 10000,
            10000, 10000, 10000, 10000, 10000,
            30.0, -25.0, 50.0, 45.0, 10000,
            10000, 57.5, 10000, 10000, 80.0,
            10000, 10000, 95.0, 10000, 10000,
            110.0, 10000, 10000, 10000, 10000,
            10000, 10000, 10000, 21.0, 7.5,
            7.5, 10000, 10000, 10000, 10000,
            10000, 10000, -30.0, 10000, 10000,
            -40.0, 10000, 10000, -55.0, 10000,
            10000, -70.0, 10000, 10000, -80.0,
            10000, 10000, 10000, 10000, 21.0,
            7.5, 7.5, 10000, 10000, 10000,
            10000, 10000, 10000
        ),
        minDZ = cms.vdouble(
            -16.0, -16.0, 0.0, 0.0, -25.0,
            -25.0, -13.0, -17.0, 0.0, 0.0,
            -19.0, -21.0, -9.0, 0.0, -13.0,
            -15.0, -35.0, 10.0, -22.0, -10000,
            -10000, 5.0, -10000, -10000, -10.0,
            -10000, -10000, 5.0, -10000, -10000,
            15.0, -10000, -10000, 25.0, -10000,
            -10000, -10000, -10000, -10000, -10000,
            -10000, -10000, -10000, -10000, -10000,
            -10000, -10000, -32.5, -10000, -10000,
            -50.0, -10000, -10000, -50.0, -10000,
            -10000, -70.0, -10000, -10000, -70.0,
            -10000, -10000, -10000, -10000, -10000,
            -10000, -10000, -10000, -10000, -10000,
            -10000, -50.0, -40.0
        ),
        minInner = cms.vdouble(
            -17.0, -14.0, 3.0, 7.0, -10000,
            -10000, -17.0, -15.0, 6.0, 9.0,
            -10000, -10000, -18.0, 11.0, -10000,
            -10, -20, 10, -20, 0,
            0, 11.6, 0, 0, 11.6,
            0, 0, 11.6, 0, 0,
            11.8, 0, 0, 0, 0,
            0, 0, 12.5, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 11.6, 0, 0,
            11.6, 0, 0, 11.6, 0,
            0, 11.8, 0, 0, 0,
            0, 0, 0, 12.5, 0,
            0, 0, 0, 0, 0,
            0, -1200, -1200
        ),
        minOuter = cms.vdouble(
            -10000, -10000, 0, 0, 0,
            0, -10000, -10000, 6.5, 6.5,
            6.5, 6.5, -10000, 11.7, 11.7,
            -30.0, -50.0, 25.0, -45.0, 3.5,
            3.5, 30.0, 3.5, 3.5, 40.0,
            3.5, 3.5, 55.0, 3.5, 3.5,
            70.0, 3.5, 3.5, 80.0, 4.0,
            4.0, 3.5, 20.0, 6.0, 0,
            0, 0, 7.0, 7.0, 7.0,
            3.5, 3.5, -57.5, 3.5, 3.5,
            -80.0, 3.5, 3.5, -95.0, 3.5,
            3.5, -110.0, 3.5, 3.5, -10000,
            4.0, 4.0, 3.5, 20.0, 6.0,
            0, 0, 0, 7.0, 7.0,
            7.0, -10000, -10000
        ),
        pairGraph = cms.vuint32(
            0, 1, 0, 2, 0,
            4, 0, 5, 0, 16,
            0, 17, 1, 2, 1,
            3, 1, 4, 1, 5,
            1, 16, 1, 17, 2,
            3, 2, 4, 2, 16,
            2, 28, 2, 28, 2,
            28, 3, 28, 4, 5,
            4, 6, 4, 28, 5,
            6, 5, 7, 5, 28,
            6, 7, 6, 8, 6,
            28, 7, 8, 7, 9,
            7, 28, 8, 9, 8,
            10, 8, 28, 9, 10,
            9, 11, 10, 11, 10,
            12, 11, 12, 11, 13,
            11, 14, 11, 15, 12,
            13, 13, 14, 14, 15,
            16, 17, 16, 18, 16,
            28, 17, 18, 17, 19,
            17, 28, 18, 19, 18,
            20, 18, 28, 19, 20,
            19, 21, 19, 28, 20,
            21, 20, 22, 20, 28,
            21, 22, 21, 23, 22,
            23, 22, 24, 23, 24,
            23, 25, 23, 26, 23,
            27, 24, 25, 25, 26,
            26, 27, 28, 29, 29,
            30
        ),
        phiCuts = cms.vint32(
            350, 600, 450, 522, 450,
            522, 400, 650, 500, 730,
            500, 730, 350, 400, 400,
            1200, 1200, 1200, 1000, 300,
            522, 1000, 300, 522, 1000,
            250, 522, 1000, 250, 522,
            1000, 250, 522, 850, 300,
            522, 240, 650, 300, 200,
            220, 250, 250, 250, 250,
            300, 522, 1000, 300, 522,
            1000, 250, 522, 1000, 250,
            522, 1000, 250, 522, 1000,
            300, 522, 240, 650, 300,
            200, 220, 250, 250, 250,
            250, 1100, 1250
        ),
        ptCuts = cms.vdouble(
            0.7, 0.8, 0.6, 0.85, 0.6,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            2.0, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85, 0.85, 0.85,
            0.85, 0.85, 0.85
        ),
        skipsLayers = cms.vuint32(
            0, 1, 0, 1, 0,
            1, 0, 1, 0, 1,
            0, 1, 0, 0, 0,
            0, 0, 0, 0, 0,
            1, 0, 0, 1, 0,
            0, 1, 0, 0, 1,
            0, 0, 1, 0, 0,
            1, 0, 1, 0, 0,
            0, 0, 0, 0, 0,
            0, 1, 0, 0, 1,
            0, 0, 1, 0, 0,
            1, 0, 0, 1, 0,
            0, 1, 0, 1, 0,
            0, 0, 0, 0, 0,
            0, 0, 0
        ),
        startMaxInnerR = cms.vdouble(
            99.0, 99.0, 99.0, 99.0, 99.0,
            99.0, 99.0, 99.0, 99.0, 99.0,
            99.0, 99.0, 99.0, 99.0, 99.0,
            99.0, 99.0, 99.0, 99.0, 99.0,
            99.0, 99.0, 99.0, 99.0, 99.0,
            99.0, 99.0, 99.0, 99.0, 99.0,
            99.0
        ),
        startingPairs = cms.vuint32(
            0, 1, 2, 3, 4,
            5, 6, 8, 10, 19,
            22, 25, 28, 31, 34,
            36, 45, 48, 51, 54,
            57, 60, 62
        )
    ),
    hardCurvCut = cms.double(0.02),
    lateFishbone = cms.bool(False),
    maxDYPred = cms.int32(24),
    maxDYsize = cms.int32(20),
    maxDYsize12 = cms.int32(15),
    maxNumberOfDoublets = cms.string('6000000.0'),
    maxNumberOfTuples = cms.string('200000.0'),
    minHitsForSharingCut = cms.uint32(1),
    minHitsPerNtuplet = cms.uint32(4),
    minYsizeB1 = cms.int32(15),
    minYsizeB2 = cms.int32(14),
    onlySameLayersFishbone = cms.bool(True),
    pixelRecHitSrc = cms.InputTag("hltPhase2PixelRecHitsExtendedSoA"),
    ptmin = cms.double(0.9),
    trackQualityCuts = cms.PSet(
        maxChi2 = cms.double(5.0),
        maxChi2Quintuplets = cms.double(3.0),
        maxChi2TripletsOrQuadruplets = cms.double(1.0),
        maxTip = cms.double(0.3),
        maxZip = cms.double(12),
        minPt = cms.double(0.9)
    ),
    useRiemannFit = cms.bool(False),
    useSimpleTripletCleaner = cms.bool(True)
)

def _exclude_OT_layers(hltPhase2PixelTracksSoA, layers_to_exclude = [28, 29, 30]):
    keep_indices = []
    num_pairs = len(hltPhase2PixelTracksSoA.geometry.pairGraph) // 2
    for i in range(num_pairs):
        a = hltPhase2PixelTracksSoA.geometry.pairGraph[2*i]
        b = hltPhase2PixelTracksSoA.geometry.pairGraph[2*i + 1]
        if a not in layers_to_exclude and b not in layers_to_exclude:
            keep_indices.append(i)
    # Now update in place
    # For pairGraph, build the new flat list from kept pairs
    new_pairGraph = []
    for i in keep_indices:
        new_pairGraph.extend([hltPhase2PixelTracksSoA.geometry.pairGraph[2*i], hltPhase2PixelTracksSoA.geometry.pairGraph[2*i+1]])

    hltPhase2PixelTracksSoA.geometry.pairGraph[:] = new_pairGraph
    # Update all other lists in place
    hltPhase2PixelTracksSoA.geometry.skipsLayers[:] = [hltPhase2PixelTracksSoA.geometry.skipsLayers[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.phiCuts[:] = [hltPhase2PixelTracksSoA.geometry.phiCuts[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.minInnerR[:] = [hltPhase2PixelTracksSoA.geometry.minInnerR[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.maxInnerR[:] = [hltPhase2PixelTracksSoA.geometry.maxInnerR[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.minOuterR[:] = [hltPhase2PixelTracksSoA.geometry.minOuterR[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.maxOuterR[:] = [hltPhase2PixelTracksSoA.geometry.maxOuterR[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.maxDR[:] = [hltPhase2PixelTracksSoA.geometry.maxDR[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.minInnerZ[:] = [hltPhase2PixelTracksSoA.geometry.minInnerZ[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.maxInnerZ[:] = [hltPhase2PixelTracksSoA.geometry.maxInnerZ[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.minOuterZ[:] = [hltPhase2PixelTracksSoA.geometry.minOuterZ[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.maxOuterZ[:] = [hltPhase2PixelTracksSoA.geometry.maxOuterZ[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.minDZ[:] = [hltPhase2PixelTracksSoA.geometry.minDZ[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.maxDZ[:] = [hltPhase2PixelTracksSoA.geometry.maxDZ[i] for i in keep_indices]
    hltPhase2PixelTracksSoA.geometry.ptCuts[:] = [hltPhase2PixelTracksSoA.geometry.ptCuts[i] for i in keep_indices]

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toReplaceWith(hltPhase2PixelTracksSoA, _hltPhase2PixelTracksSoA)

#print("Using {} pair connections: {}".format(len(hltPhase2PixelTracksSoA.geometry.pairGraph) // 2, hltPhase2PixelTracksSoA.geometry.pairGraph))
