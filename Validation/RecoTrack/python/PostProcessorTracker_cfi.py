import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from Configuration.Eras.Modifier_fastSim_cff import fastSim

def _addNoFlow(module):
    _noflowSeen = set()
    for eff in module.efficiency.value():
        tmp = eff.split(" ")
        if "cut" in tmp[0]:
            continue
        ind = -1
        if tmp[ind] == "fake" or tmp[ind] == "simpleratio":
            ind = -2
        if not tmp[ind] in _noflowSeen:
            module.noFlowDists.append(tmp[ind])
        if not tmp[ind-1] in _noflowSeen:
            module.noFlowDists.append(tmp[ind-1])

_defaultSubdirs = ["Tracking/Track/*", "Tracking/TrackTPPtLess09/*", "Tracking/TrackFromPV/*", "Tracking/TrackFromPVAllTP/*", "Tracking/TrackAllTPEffic/*", "Tracking/TrackBuilding/*","Tracking/TrackConversion/*", "Tracking/TrackGsf/*"]
_defaultSubdirsSummary = [e.replace("/*","") for e in _defaultSubdirs]

def makeEfficencyBundle(histoSuffix, label):
    """Produces list of the efficiencies (effic, duplicatesRate, fakerate, etc.) for a given variable.
    histoSuffix (str): suffix of the histograms used for calculating the efficencies, e.g. `eta` for pseudorapidity since histograms are named  like this `num_simul_eta`.
    label (str): label to put in the plot title for the quantity, e.g. `#eta` for the pseudorapidity.
    """
    suffix = "_vs_" + histoSuffix

    # eta and pT are historically named differently for no good reason...
    if histoSuffix == "eta":
        suffix = ""
    elif histoSuffix == "pT":
        suffix = "_Pt"

    strings = (suffix, label, histoSuffix, histoSuffix)
    # efficiency is again special for pT...
    stringsEff = (suffix if histoSuffix != "pT" else suffix[1:], label, histoSuffix, histoSuffix)

    effList = [
        "effic%s 'Efficiency vs %s' num_assoc(simToReco)_%s num_simul_%s" % stringsEff,
        "techEffic%s 'Technical efficiency vs %s' num_assoc(reconstructableSimToReco)_%s num_reconstructableSim_%s" % stringsEff,
        "duplicatesRate%s 'Duplicates Rate vs %s' num_duplicate_%s num_reco_%s" % strings,
        "chargeMisIdRate%s 'Charge MisID Rate vs %s' num_chargemisid_%s num_reco_%s" % strings,
        "pileuprate%s 'Pileup Rate vs %s' num_pileup_%s num_reco_%s" % strings,
        "fakerate%s 'Fake rate vs %s' num_assoc(recoToSim)_%s num_reco_%s fake" % strings,
    ]

    return effList

postProcessorTrack = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring(_defaultSubdirs),
    efficiency = cms.vstring(
        makeEfficencyBundle("eta", "#eta") +
        makeEfficencyBundle("pT", "p_{T}") +
        makeEfficencyBundle("hit", "# hits") +
        makeEfficencyBundle("layer", "# layers") +
        makeEfficencyBundle("pixellayer", "# pixel layers") +
        makeEfficencyBundle("3Dlayer", "# 3D layers") +
        makeEfficencyBundle("pu", "number of vertices/PU") +
        makeEfficencyBundle("phi", "#phi") +
        makeEfficencyBundle("dxy", "d_{xy}") +
        makeEfficencyBundle("dz", "d_{z}") +
        makeEfficencyBundle("dxypv", "d_{xy}(PV)") +
        makeEfficencyBundle("dzpv", "d_{z}(PV)") +
        makeEfficencyBundle("dxypv_zoomed", "d_{xy}(PV)") +
        makeEfficencyBundle("dzpv_zoomed", "d_{z}(PV)") +
        makeEfficencyBundle("vertpos", "r_{vertex} (vertical displacement)") +
        makeEfficencyBundle("zpos", "z_{vertex} (longitudinal displacement)") +
        makeEfficencyBundle("dr", "dr") +
        makeEfficencyBundle("drj", "dr (track,jet)") +
        [
        #    "efficPtvseta 'Efficiency in p_{T}-#eta plane' num_assoc(simToReco)_pTvseta num_simul_pTvseta",
        #   "technicalEfficPtvseta 'Technical efficiency in p_{T}-#eta plane' num_assoc(reconstructableSimToReco)_pTvseta num_reconstructableSim_pTvseta",
        #    "duplicatesRate_Ptvseta 'Duplicates Rate in (p_{T}-#eta) plane' num_duplicate_pTvseta num_reco_pTvseta",
        #    "chargeMisIdRate_Ptvseta 'Charge MisID Rate in (p_{T}-#eta) plane' num_chargemisid_pTvseta num_reco_pTvseta",
        #    "pileuprate_Ptvseta 'Pileup rate in (p_{T}-#eta) plane' num_pileup_pTvseta num_reco_pTvseta",
        #    "fakeratePtvseta 'Fake rate in (p_{T}-#eta) plane' num_assoc(recoToSim)_pTvseta num_reco_pTvseta fake",
            "duplicatesRate_chi2 'Duplicates Rate vs normalized #chi^{2}' num_duplicate_chi2 num_reco_chi2",
            "duplicatesRate_seedingLayerSet 'Duplicates rate vs. seedingLayerSet' num_duplicate_seedingLayerSet num_reco_seedingLayerSet",
            "chargeMisIdRate_chi2 'Charge MisID Rate vs normalized #chi^{2}' num_chargemisid_chi2 num_reco_chi2",
            "effic_vertcount_barrel 'efficiency in barrel vs N of pileup vertices' num_assoc(simToReco)_vertcount_barrel num_simul_vertcount_barrel",
            "effic_vertcount_fwdpos 'efficiency in endcap(+) vs N of pileup vertices' num_assoc(simToReco)_vertcount_fwdpos num_simul_vertcount_fwdpos",
            "effic_vertcount_fwdneg 'efficiency in endcap(-) vs N of pileup vertices' num_assoc(simToReco)_vertcount_fwdneg num_simul_vertcount_fwdneg",
            "effic_vertz_barrel 'efficiency in barrel vs z of primary interaction vertex' num_assoc(simToReco)_vertz_barrel num_simul_vertz_barrel",
            "effic_vertz_fwdpos 'efficiency in endcap(+) vs z of primary interaction vertex' num_assoc(simToReco)_vertz_fwdpos num_simul_vertz_fwdpos",
            "effic_vertz_fwdneg 'efficiency in endcap(-) vs z of primary interaction vertex' num_assoc(simToReco)_vertz_fwdneg num_simul_vertz_fwdneg",
            "technicalEffic_vertcount_barrel 'Technical efficiency in barrel vs N of pileup vertices' num_assoc(reconstructableSimToReco)_vertcount_barrel num_reconstructableSim_vertcount_barrel",
            "technicalEffic_vertcount_fwdpos 'Technical efficiency in endcap(+) vs N of pileup vertices' num_assoc(reconstructableSimToReco)_vertcount_fwdpos num_reconstructableSim_vertcount_fwdpos",
            "technicalEffic_vertcount_fwdneg 'Technical efficiency in endcap(-) vs N of pileup vertices' num_assoc(reconstructableSimToReco)_vertcount_fwdneg num_reconstructableSim_vertcount_fwdneg",
            "technicalEffic_vertz_barrel 'Technical efficiency in barrel vs z of primary interaction vertex' num_assoc(reconstructableSimToReco)_vertz_barrel num_reconstructableSim_vertz_barrel",
            "technicalEffic_vertz_fwdpos 'Technical efficiency in endcap(+) vs z of primary interaction vertex' num_assoc(reconstructableSimToReco)_vertz_fwdpos num_reconstructableSim_vertz_fwdpos",
            "technicalEffic_vertz_fwdneg 'Technical efficiency in endcap(-) vs z of primary interaction vertex' num_assoc(reconstructableSimToReco)_vertz_fwdneg num_reconstructableSim_vertz_fwdneg",
            "pileuprate_chi2 'Pileup rate vs normalized #chi^{2}' num_pileup_chi2 num_reco_chi2",
            "pileuprate_seedingLayerSet 'Pileup rate vs. seedingLayerSet' num_pileup_seedingLayerSet num_reco_seedingLayerSet",
            "fakerate_vs_chi2 'Fake rate vs normalized #chi^{2}' num_assoc(recoToSim)_chi2 num_reco_chi2 fake",
            "fakerate_vs_seedingLayerSet 'Fake rate vs. seedingLayerSet' num_assoc(recoToSim)_seedingLayerSet num_reco_seedingLayerSet fake",
            "fakerate_vertcount_barrel 'fake rate in barrel vs N of pileup vertices' num_assoc(recoToSim)_vertcount_barrel num_reco_vertcount_barrel fake",
            "fakerate_vertcount_fwdpos 'fake rate in endcap(+) vs N of pileup vertices' num_assoc(recoToSim)_vertcount_fwdpos num_reco_vertcount_fwdpos fake",
            "fakerate_vertcount_fwdneg 'fake rate in endcap(-) vs N of pileup vertices' num_assoc(recoToSim)_vertcount_fwdneg num_reco_vertcount_fwdneg fake",
            "fakerate_ootpu_entire 'fake rate from out of time pileup vs N of pileup vertices' num_assoc(recoToSim)_ootpu_entire num_reco_ootpu_entire",
            "fakerate_ootpu_barrel 'fake rate from out of time pileup in barrel vs N of pileup vertices' num_assoc(recoToSim)_ootpu_barrel num_reco_ootpu_barrel",
            "fakerate_ootpu_fwdpos 'fake rate from out of time pileup in endcap(+) vs N of pileup vertices' num_assoc(recoToSim)_ootpu_fwdpos num_reco_ootpu_fwdpos",
            "fakerate_ootpu_fwdneg 'fake rate from out of time pileup in endcap(-) vs N of pileup vertices' num_assoc(recoToSim)_ootpu_fwdneg num_reco_ootpu_fwdneg",

            "effic_vs_dzpvcut 'Efficiency vs. dz (PV)' num_assoc(simToReco)_dzpvcut num_simul_dzpvcut",
            "effic_vs_dzpvcut2 'Efficiency (tracking eff factorized out) vs. dz (PV)' num_assoc(simToReco)_dzpvcut num_simul2_dzpvcut",
            "technicalEffic_vs_dzpvcut 'Technical efficiency vs. dz (PV)' num_assoc(reconstructableSimToReco)_dzpvcut num_reconstructableSim_dzpvcut",
            "fakerate_vs_dzpvcut 'Fake rate vs. dz(PV)' num_assoc(recoToSim)_dzpvcut num_reco_dzpvcut fake",
            "pileuprate_dzpvcut 'Pileup rate vs. dz(PV)' num_pileup_dzpvcut num_reco_dzpvcut",

            "effic_vs_dzpvsigcut 'Efficiency vs. dz(PV)/dzError' num_assoc(simToReco)_dzpvsigcut num_simul_dzpvsigcut",
            "effic_vs_dzpvsigcut2 'Efficiency (tracking eff factorized out) vs. dz(PV)/dzError' num_assoc(simToReco)_dzpvsigcut num_simul2_dzpvsigcut",
            "fakerate_vs_dzpvsigcut 'Fake rate vs. dz(PV)/dzError' num_assoc(recoToSim)_dzpvsigcut num_reco_dzpvsigcut fake",
            "pileuprate_dzpvsigcut 'Pileup rate vs. dz(PV)/dzError' num_pileup_dzpvsigcut num_reco_dzpvsigcut",

            "effic_vs_simpvz 'Efficiency vs. sim PV z' num_assoc(simToReco)_simpvz num_simul_simpvz",
            "fakerate_vs_simpvz 'Fake rate vs. sim PV z' num_assoc(recoToSim)_simpvz num_reco_simpvz fake",
            "duplicatesRate_simpvz 'Duplicates Rate vs sim PV z' num_duplicate_simpvz num_reco_simpvz",
            "pileuprate_simpvz 'Pileup rate vs. sim PV z' num_pileup_simpvz num_reco_simpvz",

            "fakerate_vs_mva1 'Fake rate vs. MVA1' num_assoc(recoToSim)_mva1 num_reco_mva1 fake",
            "fakerate_vs_mva2 'Fake rate vs. MVA2' num_assoc(recoToSim)_mva2 num_reco_mva2 fake",
            "fakerate_vs_mva3 'Fake rate vs. MVA3' num_assoc(recoToSim)_mva3 num_reco_mva3 fake",

            "effic_vs_mva1cut 'Efficiency (tracking eff factorized out) vs. MVA1' num_assoc(simToReco)_mva1cut num_simul2_mva1cut",
            "fakerate_vs_mva1cut 'Fake rate vs. MVA1' num_assoc(recoToSim)_mva1cut num_reco_mva1cut fake",
            "effic_vs_mva2cut 'Efficiency (tracking eff factorized out) vs. MVA2' num_assoc(simToReco)_mva2cut num_simul2_mva2cut",
            "effic_vs_mva2cut_hp 'Efficiency (tracking eff factorized out) vs. MVA2' num_assoc(simToReco)_mva2cut_hp num_simul2_mva2cut_hp",
            "fakerate_vs_mva2cut 'Fake rate vs. MVA2' num_assoc(recoToSim)_mva2cut num_reco_mva2cut fake",
            "fakerate_vs_mva2cut_hp 'Fake rate vs. MVA2' num_assoc(recoToSim)_mva2cut_hp num_reco_mva2cut_hp fake",
            "effic_vs_mva3cut 'Efficiency (tracking eff factorized out) vs. MVA3' num_assoc(simToReco)_mva3cut num_simul2_mva3cut",
            "effic_vs_mva3cut_hp 'Efficiency (tracking eff factorized out) vs. MVA3' num_assoc(simToReco)_mva3cut_hp num_simul2_mva3cut_hp",
            "fakerate_vs_mva3cut 'Fake rate vs. MVA3' num_assoc(recoToSim)_mva3cut num_reco_mva3cut fake",
            "fakerate_vs_mva3cut_hp 'Fake rate vs. MVA3' num_assoc(recoToSim)_mva3cut_hp num_reco_mva3cut_hp fake",
        ]
    ),
    resolution = cms.vstring(
                             "cotThetares_vs_eta '#sigma(cot(#theta)) vs #eta' cotThetares_vs_eta",
                             "cotThetares_vs_pt '#sigma(cot(#theta)) vs p_{T}' cotThetares_vs_pt",
                             "h_dxypulleta 'd_{xy} Pull vs #eta' dxypull_vs_eta",
                             "h_dxypullpt 'd_{xy} Pull vs p_{T}' dxypull_vs_pt",
                             "dxyres_vs_eta '#sigma(d_{xy}) vs #eta' dxyres_vs_eta",
                             "dxyres_vs_phi '#sigma(d_{xy}) vs #phi' dxyres_vs_phi",
                             "dxyres_vs_pt '#sigma(d_{xy}) vs p_{T}' dxyres_vs_pt",
                             "h_dzpulleta 'd_{z} Pull vs #eta' dzpull_vs_eta",
                             "h_dzpullpt 'd_{z} Pull vs p_{T}' dzpull_vs_pt",
                             "dzres_vs_eta '#sigma(d_{z}) vs #eta' dzres_vs_eta",
                             "dzres_vs_phi '#sigma(d_{z}) vs #phi' dzres_vs_phi",
                             "dzres_vs_pt '#sigma(d_{z}) vs p_{T}' dzres_vs_pt",
                             "etares_vs_eta '#sigma(#eta) vs #eta' etares_vs_eta",
                             "h_phipulleta '#phi Pull vs #eta' phipull_vs_eta",
                             "h_phipullpt '#phi Pull vs p_{T}' phipull_vs_pt",
                             "h_phipullphi '#phi Pull vs #phi' phipull_vs_phi",
                             "phires_vs_eta '#sigma(#phi) vs #eta' phires_vs_eta",
                             "phires_vs_phi '#sigma(#phi) vs #phi' phires_vs_phi",
                             "phires_vs_pt '#sigma(#phi) vs p_{T}' phires_vs_pt",
                             "h_ptpulleta 'p_{T} Pull vs #eta' ptpull_vs_eta",
                             "h_ptpullpt 'p_{T} Pull vs p_{T}' ptpull_vs_pt",
                             "h_ptpullphi 'p_{T} Pull vs #phi' ptpull_vs_phi",
                             "ptres_vs_eta '#sigma(p_{T}) vs #eta' ptres_vs_eta",
                             "ptres_vs_phi '#sigma(p_{T}) vs #phi' ptres_vs_phi",
                             "ptres_vs_pt '#sigma(p_{T}) vs p_{T}' ptres_vs_pt",
                             "h_thetapulleta '#theta Pull vs #eta' thetapull_vs_eta",
                             "h_thetapullpt '#theta Pull vs p_{T}' thetapull_vs_pt",
                             "h_thetapullphi '#theta Pull vs #phi' thetapull_vs_phi"
                             ),
    cumulativeDists = cms.untracked.vstring(
        "num_reco_dzpvcut",
        "num_assoc(recoToSim)_dzpvcut",
        "num_assoc(simToReco)_dzpvcut",
        "num_simul_dzpvcut",
        "num_simul2_dzpvcut",
        "num_pileup_dzpvcut",
        "num_reco_dzpvsigcut",
        "num_assoc(recoToSim)_dzpvsigcut",
        "num_assoc(simToReco)_dzpvsigcut",
        "num_simul_dzpvsigcut",
        "num_simul2_dzpvsigcut",
        "num_pileup_dzpvsigcut",
        "num_reco_mva1cut descending",
        "num_reco_mva2cut descending",
        "num_reco_mva2cut_hp descending",
        "num_reco_mva3cut descending",
        "num_reco_mva3cut_hp descending",
        "num_assoc(recoToSim)_mva1cut descending",
        "num_assoc(recoToSim)_mva2cut descending",
        "num_assoc(recoToSim)_mva2cut_hp descending",
        "num_assoc(recoToSim)_mva3cut descending",
        "num_assoc(recoToSim)_mva3cut_hp descending",
        "num_assoc(simToReco)_mva1cut descending",
        "num_assoc(simToReco)_mva2cut descending",
        "num_assoc(simToReco)_mva2cut_hp descending",
        "num_assoc(simToReco)_mva3cut descending",
        "num_assoc(simToReco)_mva3cut_hp descending",
        "num_simul2_mva1cut descending",
        "num_simul2_mva2cut descending",
        "num_simul2_mva2cut_hp descending",
        "num_simul2_mva3cut descending",
        "num_simul2_mva3cut_hp descending",
    ),
    noFlowDists = cms.untracked.vstring(),
    outputFileName = cms.untracked.string("")
)

_addNoFlow(postProcessorTrack)

postProcessorTrack2D = DQMEDHarvester("DQMGenericClient",
    makeGlobalEffienciesPlot = cms.untracked.bool(False),
    subDirs = cms.untracked.vstring(_defaultSubdirs),
    efficiency = cms.vstring(
    "efficPtvseta 'Efficiency in p_{T}-#eta plane' num_assoc(simToReco)_pTvseta num_simul_pTvseta",
    "duplicatesRate_Ptvseta 'Duplicates Rate in (p_{T}-#eta) plane' num_duplicate_pTvseta num_reco_pTvseta",
    "chargeMisIdRate_Ptvseta 'Charge MisID Rate in (p_{T}-#eta) plane' num_chargemisid_pTvseta num_reco_pTvseta",
    "pileuprate_Ptvseta 'Pileup rate in (p_{T}-#eta) plane' num_pileup_pTvseta num_reco_pTvseta",
    "fakeratePtvseta 'Fake rate in (p_{T}-#eta) plane' num_assoc(recoToSim)_pTvseta num_reco_pTvseta fake",
    ),
    resolution = cms.vstring(),
    noFlowDists = cms.untracked.vstring(),
    outputFileName = cms.untracked.string("")
)

_addNoFlow(postProcessorTrack2D)

# nrec/nsim makes sense only for
# - all tracks vs. all in-time TrackingParticles
# - PV tracks vs. signal TrackingParticles
postProcessorTrackNrecVsNsim = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring("Tracking/TrackFromPV/*", "Tracking/TrackAllTPEffic/*"),
    efficiency = cms.vstring(
        "nrecPerNsim 'Tracks/TrackingParticles vs #eta' num_reco2_eta num_simul_eta simpleratio",
        "nrecPerNsimPt 'Tracks/TrackingParticles vs p_{T}' num_reco2_pT num_simul_pT simpleratio",
#        "nrecPerNsimPtvseta 'Tracks/TrackingParticles in (p_{T}-#eta) plane' num_reco2_pTvseta num_simul_pTvseta simpleratio",
        "nrecPerNsim_vs_pu 'Tracks/TrackingParticles vs pu' num_reco2_pu num_simul_pu simpleratio",
    ),
    resolution = cms.vstring(),
    noFlowDists = cms.untracked.vstring(),
)
_addNoFlow(postProcessorTrackNrecVsNsim)
postProcessorTrackNrecVsNsim2D = DQMEDHarvester("DQMGenericClient",
    makeGlobalEffienciesPlot = cms.untracked.bool(False),
    subDirs = cms.untracked.vstring("Tracking/TrackFromPV/*", "Tracking/TrackAllTPEffic/*"),
    efficiency = cms.vstring(
        "nrecPerNsimPtvseta 'Tracks/TrackingParticles in (p_{T}-#eta) plane' num_reco2_pTvseta num_simul_pTvseta simpleratio",
    ),
    resolution = cms.vstring(),
    noFlowDists = cms.untracked.vstring(),
)
_addNoFlow(postProcessorTrackNrecVsNsim2D)


postProcessorTrackSummary = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring(_defaultSubdirsSummary),
    efficiency = cms.vstring(
    "effic_vs_coll 'Efficiency vs track collection' num_assoc(simToReco)_coll num_simul_coll",
    "techEffic_vs_coll 'Technical efficiency vs track collection' num_assoc(reconstructableSimToReco)_coll num_reconstructableSim_coll",
    "effic_vs_coll_allPt 'Efficiency vs track collection' num_assoc(simToReco)_coll_allPt num_simul_coll_allPt",
    "duplicatesRate_coll 'Duplicates Rate vs track collection' num_duplicate_coll num_reco_coll",
    "pileuprate_coll 'Pileup rate vs track collection' num_pileup_coll num_reco_coll",
    "fakerate_vs_coll 'Fake rate vs track collection' num_assoc(recoToSim)_coll num_reco_coll fake",
    ),
    resolution = cms.vstring(),
    noFlowDists = cms.untracked.vstring(),
)
_addNoFlow(postProcessorTrackSummary)

postProcessorTrackSequence = cms.Sequence(
    postProcessorTrack+
    postProcessorTrackNrecVsNsim+
    postProcessorTrackSummary
)

from Configuration.ProcessModifiers.seedingDeepCore_cff import seedingDeepCore
postProcessorTrackDeepCore = postProcessorTrack.clone()
postProcessorTrackDeepCore.subDirs.extend(["Tracking/JetCore/*"])
seedingDeepCore.toReplaceWith(postProcessorTrack,postProcessorTrackDeepCore)
postProcessorTrackSummaryDeepCore = postProcessorTrackSummary.clone()
postProcessorTrackSummaryDeepCore.subDirs.extend(["Tracking/JetCore/*"])
seedingDeepCore.toReplaceWith(postProcessorTrackSummary,postProcessorTrackSummaryDeepCore)
postProcessorTrack2DDeepCore = postProcessorTrack2D.clone()
postProcessorTrack2DDeepCore.subDirs.extend(["Tracking/JetCore/*"])
seedingDeepCore.toReplaceWith(postProcessorTrack2D,postProcessorTrack2DDeepCore)


fastSim.toModify(postProcessorTrack, subDirs = [e for e in _defaultSubdirs if e not in ["Tracking/TrackGsf/*","Tracking/TrackConversion/*"]])
fastSim.toModify(postProcessorTrackSummary, subDirs = [e for e in _defaultSubdirsSummary if e not in ["Tracking/TrackGsf","Tracking/TrackConversion"]])

#######
# Define a standalone seuquence to support the Standalone harvesting mode
# see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMultiTrackValidator#cmsDriver_MTV_alone_i_e_standalone for more information
########

postProcessorTrackStandalone = postProcessorTrack.clone(
    subDirs = _defaultSubdirs+["Tracking/TrackBHadron/*"]
)
postProcessorTrackSummaryStandalone = postProcessorTrackSummary.clone(
    subDirs = _defaultSubdirs+["Tracking/TrackBHadron"]
)

postProcessorTrackSequenceStandalone = cms.Sequence(
    postProcessorTrackStandalone+
    postProcessorTrackNrecVsNsim+
    postProcessorTrackSummaryStandalone
)



postProcessorTrackPhase2 = postProcessorTrack.clone()
postProcessorTrackPhase2.subDirs.extend(["Tracking/TrackTPEtaGreater2p7/*"])
postProcessorTrackSummaryPhase2 = postProcessorTrackSummary.clone()
postProcessorTrackSummaryPhase2.subDirs.extend(["Tracking/TrackTPEtaGreater2p7/*"])

from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
phase2_tracker.toReplaceWith(postProcessorTrack,postProcessorTrackPhase2)
phase2_tracker.toReplaceWith(postProcessorTrackSummary,postProcessorTrackSummaryPhase2)

from Configuration.ProcessModifiers.displacedTrackValidation_cff import displacedTrackValidation
postProcessorTrackDisplaced = postProcessorTrack.clone()
postProcessorTrackDisplaced.subDirs.extend(["Tracking/TrackDisplaced/*"])
postProcessorTrackSummaryDisplaced = postProcessorTrackSummary.clone()
postProcessorTrackSummaryDisplaced.subDirs.extend(["Tracking/TrackDisplaced/*"])
displacedTrackValidation.toReplaceWith(postProcessorTrack,postProcessorTrackDisplaced)
displacedTrackValidation.toReplaceWith(postProcessorTrackSummary,postProcessorTrackSummaryDisplaced)

from Configuration.ProcessModifiers.pp_on_AA_cff import pp_on_AA

_defaultSubdirsHIon =  _defaultSubdirs + ["Tracking/HIPixelTrack/*"]

(pp_on_AA & ~phase2_tracker).toModify(postProcessorTrack,subDirs = _defaultSubdirsHIon)
(pp_on_AA & ~phase2_tracker).toModify(postProcessorTrack2D,subDirs = _defaultSubdirsHIon)
(pp_on_AA & ~phase2_tracker).toModify(postProcessorTrackSummary,subDirs = _defaultSubdirsHIon)


postProcessorTrackTrackingOnly = postProcessorTrack.clone()
postProcessorTrackTrackingOnly.subDirs.extend(["Tracking/TrackBHadron/*", "Tracking/TrackSeeding/*", "Tracking/PixelTrack/*", "Tracking/PixelTrackFromPV/*", "Tracking/PixelTrackFromPVAllTP/*", "Tracking/PixelTrackBHadron/*"])
postProcessorTrackSummaryTrackingOnly = postProcessorTrackSummary.clone()
postProcessorTrackSummaryTrackingOnly.subDirs.extend(["Tracking/TrackBHadron", "Tracking/TrackSeeding", "Tracking/PixelTrack", "Tracking/PixelTrackFromPV", "Tracking/PixelTrackFromPVAllTP", "Tracking/PixelTrackBHadron"])

postProcessorTrackTrackingOnlyHIon = postProcessorTrackTrackingOnly.clone()
postProcessorTrackTrackingOnlyHIon.subDirs.extend(["Tracking/HIPixelTrack/*"])

postProcessorTrackSummaryTrackingOnlyHIon = postProcessorTrackSummaryTrackingOnly.clone()
postProcessorTrackSummaryTrackingOnlyHIon.subDirs.extend(["Tracking/HIPixelTrack"])

(pp_on_AA & ~phase2_tracker).toReplaceWith(postProcessorTrackTrackingOnly,postProcessorTrackTrackingOnlyHIon)
(pp_on_AA & ~phase2_tracker).toReplaceWith(postProcessorTrackSummaryTrackingOnly,postProcessorTrackSummaryTrackingOnlyHIon)

postProcessorTrackSequenceTrackingOnly = cms.Sequence(
    postProcessorTrackTrackingOnly+
    postProcessorTrackNrecVsNsim+
    postProcessorTrackSummaryTrackingOnly
)

fastSim.toModify(postProcessorTrackTrackingOnly,subDirs = [e for e in _defaultSubdirs if e not in ["Tracking/TrackGsf/*","Tracking/TrackConversion/*","Tracking/TrackBHadron/*"]])
fastSim.toModify(postProcessorTrackSummaryTrackingOnly,subDirs = [e for e in _defaultSubdirsSummary if e not in ["Tracking/TrackGsf","Tracking/TrackConversion","Tracking/TrackBHadron"]])
