#include "Validation/RecoTrack/interface/MTVHistoProducerAlgoForTracker.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/trackFromSeedFitFailed.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"

#include "TMath.h"
#include <TF1.h>
#include <cassert>

using namespace std;

namespace {
  typedef dqm::reco::DQMStore DQMStore;

  void setBinLabels(dqm::reco::MonitorElement*& h, const std::vector<std::string>& labels) {
    for (size_t i = 0; i < labels.size(); ++i) {
      h->setBinLabel(i + 1, labels[i]);
    }
    h->disableAlphanumeric();
  }

  void setBinLabelsAlgo(dqm::reco::MonitorElement*& h, int axis = 1) {
    for (size_t i = 0; i < reco::TrackBase::algoSize; ++i) {
      h->setBinLabel(i + 1, reco::TrackBase::algoName(static_cast<reco::TrackBase::TrackAlgorithm>(i)), axis);
    }
    h->disableAlphanumeric();
  }

  void fillMVAHistos(const std::vector<dqm::reco::MonitorElement*>& h_mva,
                     const std::vector<dqm::reco::MonitorElement*>& h_mvacut,
                     const std::vector<dqm::reco::MonitorElement*>& h_mva_hp,
                     const std::vector<dqm::reco::MonitorElement*>& h_mvacut_hp,
                     const std::vector<float>& mvas,
                     unsigned int selectsLoose,
                     unsigned int selectsHP) {
    // Fill MVA1 histos with all tracks, MVA2 histos only with tracks
    // not selected by MVA1, etc.
    for (size_t i = 0; i < mvas.size(); ++i) {
      if (i <= selectsLoose) {
        h_mva[i]->Fill(mvas[i]);
        h_mvacut[i]->Fill(mvas[i]);
      }
      if (i >= 1 && i <= selectsHP) {
        h_mva_hp[i]->Fill(mvas[i]);
        h_mvacut_hp[i]->Fill(mvas[i]);
      }
    }
  }

  void fillMVAHistos(double xval,
                     const std::vector<dqm::reco::MonitorElement*>& h_mva,
                     const std::vector<dqm::reco::MonitorElement*>& h_mva_hp,
                     const std::vector<float>& mvas,
                     unsigned int selectsLoose,
                     unsigned int selectsHP) {
    // Fill MVA1 histos with all tracks, MVA2 histos only with tracks
    // not selected by MVA1, etc.
    for (size_t i = 0; i < mvas.size(); ++i) {
      if (i <= selectsLoose) {
        h_mva[i]->Fill(xval, mvas[i]);
      }
      if (i >= 1 && i <= selectsHP) {
        h_mva_hp[i]->Fill(xval, mvas[i]);
      }
    }
  }
}  // namespace

MTVHistoProducerAlgoForTracker::MTVHistoProducerAlgoForTracker(const edm::ParameterSet& pset, const bool doSeedPlots)
    : doSeedPlots_(doSeedPlots),
      doMTDPlots_(pset.getUntrackedParameter<bool>("doMTDPlots")),
      doDzPVcutPlots_(pset.getUntrackedParameter<bool>("doDzPVcutPlots")) {
  //parameters for _vs_eta plots
  minEta = pset.getParameter<double>("minEta");
  maxEta = pset.getParameter<double>("maxEta");
  nintEta = pset.getParameter<int>("nintEta");
  useFabsEta = pset.getParameter<bool>("useFabsEta");

  //parameters for _vs_pt plots
  minPt = pset.getParameter<double>("minPt");
  maxPt = pset.getParameter<double>("maxPt");
  nintPt = pset.getParameter<int>("nintPt");
  useInvPt = pset.getParameter<bool>("useInvPt");
  useLogPt = pset.getUntrackedParameter<bool>("useLogPt", false);

  //parameters for _vs_Hit plots
  minHit = pset.getParameter<double>("minHit");
  maxHit = pset.getParameter<double>("maxHit");
  nintHit = pset.getParameter<int>("nintHit");

  //parameters for _vs_Pu plots
  minPu = pset.getParameter<double>("minPu");
  maxPu = pset.getParameter<double>("maxPu");
  nintPu = pset.getParameter<int>("nintPu");

  //parameters for _vs_Layer plots
  minLayers = pset.getParameter<double>("minLayers");
  maxLayers = pset.getParameter<double>("maxLayers");
  nintLayers = pset.getParameter<int>("nintLayers");

  //parameters for _vs_phi plots
  minPhi = pset.getParameter<double>("minPhi");
  maxPhi = pset.getParameter<double>("maxPhi");
  nintPhi = pset.getParameter<int>("nintPhi");

  //parameters for _vs_Dxy plots
  minDxy = pset.getParameter<double>("minDxy");
  maxDxy = pset.getParameter<double>("maxDxy");
  nintDxy = pset.getParameter<int>("nintDxy");

  //parameters for _vs_Dz plots
  minDz = pset.getParameter<double>("minDz");
  maxDz = pset.getParameter<double>("maxDz");
  nintDz = pset.getParameter<int>("nintDz");

  dxyDzZoom = pset.getParameter<double>("dxyDzZoom");

  //parameters for _vs_ProductionVertexTransvPosition plots
  minVertpos = pset.getParameter<double>("minVertpos");
  maxVertpos = pset.getParameter<double>("maxVertpos");
  nintVertpos = pset.getParameter<int>("nintVertpos");
  useLogVertpos = pset.getUntrackedParameter<bool>("useLogVertpos");

  //parameters for _vs_ProductionVertexZPosition plots
  minZpos = pset.getParameter<double>("minZpos");
  maxZpos = pset.getParameter<double>("maxZpos");
  nintZpos = pset.getParameter<int>("nintZpos");

  //parameters for _vs_dR plots
  mindr = pset.getParameter<double>("mindr");
  maxdr = pset.getParameter<double>("maxdr");
  nintdr = pset.getParameter<int>("nintdr");

  //parameters for _vs_dR_jet plots
  mindrj = pset.getParameter<double>("mindrj");
  maxdrj = pset.getParameter<double>("maxdrj");
  nintdrj = pset.getParameter<int>("nintdrj");

  // paramers for _vs_chi2 plots
  minChi2 = pset.getParameter<double>("minChi2");
  maxChi2 = pset.getParameter<double>("maxChi2");
  nintChi2 = pset.getParameter<int>("nintChi2");

  //parameters for dE/dx plots
  minDeDx = pset.getParameter<double>("minDeDx");
  maxDeDx = pset.getParameter<double>("maxDeDx");
  nintDeDx = pset.getParameter<int>("nintDeDx");

  //parameters for Pileup plots
  minVertcount = pset.getParameter<double>("minVertcount");
  maxVertcount = pset.getParameter<double>("maxVertcount");
  nintVertcount = pset.getParameter<int>("nintVertcount");

  //parameters for number of tracks plots
  minTracks = pset.getParameter<double>("minTracks");
  maxTracks = pset.getParameter<double>("maxTracks");
  nintTracks = pset.getParameter<int>("nintTracks");

  //parameters for vs. PV z plots
  minPVz = pset.getParameter<double>("minPVz");
  maxPVz = pset.getParameter<double>("maxPVz");
  nintPVz = pset.getParameter<int>("nintPVz");

  //parameters for vs. MVA plots
  minMVA = pset.getParameter<double>("minMVA");
  maxMVA = pset.getParameter<double>("maxMVA");
  nintMVA = pset.getParameter<int>("nintMVA");

  //parameters for resolution plots
  ptRes_rangeMin = pset.getParameter<double>("ptRes_rangeMin");
  ptRes_rangeMax = pset.getParameter<double>("ptRes_rangeMax");
  ptRes_nbin = pset.getParameter<int>("ptRes_nbin");

  phiRes_rangeMin = pset.getParameter<double>("phiRes_rangeMin");
  phiRes_rangeMax = pset.getParameter<double>("phiRes_rangeMax");
  phiRes_nbin = pset.getParameter<int>("phiRes_nbin");

  cotThetaRes_rangeMin = pset.getParameter<double>("cotThetaRes_rangeMin");
  cotThetaRes_rangeMax = pset.getParameter<double>("cotThetaRes_rangeMax");
  cotThetaRes_nbin = pset.getParameter<int>("cotThetaRes_nbin");

  dxyRes_rangeMin = pset.getParameter<double>("dxyRes_rangeMin");
  dxyRes_rangeMax = pset.getParameter<double>("dxyRes_rangeMax");
  dxyRes_nbin = pset.getParameter<int>("dxyRes_nbin");

  dzRes_rangeMin = pset.getParameter<double>("dzRes_rangeMin");
  dzRes_rangeMax = pset.getParameter<double>("dzRes_rangeMax");
  dzRes_nbin = pset.getParameter<int>("dzRes_nbin");

  maxDzpvCum = pset.getParameter<double>("maxDzpvCumulative");
  nintDzpvCum = pset.getParameter<int>("nintDzpvCumulative");

  //--- tracking particle selectors for efficiency measurements
  using namespace edm;
  using namespace reco::modules;
  auto initTPselector = [&](auto& sel, auto& name) {
    sel = std::make_unique<TrackingParticleSelector>(
        ParameterAdapter<TrackingParticleSelector>::make(pset.getParameter<ParameterSet>(name)));
  };
  auto initTrackSelector = [&](auto& sel, auto& name) {
    sel = makeRecoTrackSelectorFromTPSelectorParameters(pset.getParameter<ParameterSet>(name));
  };
  auto initGPselector = [&](auto& sel, auto& name) {
    sel = std::make_unique<GenParticleCustomSelector>(
        ParameterAdapter<GenParticleCustomSelector>::make(pset.getParameter<ParameterSet>(name)));
  };

  initTPselector(generalTpSelector, "generalTpSelector");
  initTPselector(TpSelectorForTechnicalEfficiency, "TpSelectorForTechnicalEfficiency");
  initTPselector(TpSelectorForEfficiencyVsEta, "TpSelectorForEfficiencyVsEta");
  initTPselector(TpSelectorForEfficiencyVsPhi, "TpSelectorForEfficiencyVsPhi");
  initTPselector(TpSelectorForEfficiencyVsPt, "TpSelectorForEfficiencyVsPt");
  initTPselector(TpSelectorForEfficiencyVsVTXR, "TpSelectorForEfficiencyVsVTXR");
  initTPselector(TpSelectorForEfficiencyVsVTXZ, "TpSelectorForEfficiencyVsVTXZ");

  initTrackSelector(trackSelectorVsEta, "TpSelectorForEfficiencyVsEta");
  initTrackSelector(trackSelectorVsPhi, "TpSelectorForEfficiencyVsPhi");
  initTrackSelector(trackSelectorVsPt, "TpSelectorForEfficiencyVsPt");

  initGPselector(generalGpSelector, "generalGpSelector");
  initGPselector(GpSelectorForTechnicalEfficiency, "GpSelectorForTechnicalEfficiency");
  initGPselector(GpSelectorForEfficiencyVsEta, "GpSelectorForEfficiencyVsEta");
  initGPselector(GpSelectorForEfficiencyVsPhi, "GpSelectorForEfficiencyVsPhi");
  initGPselector(GpSelectorForEfficiencyVsPt, "GpSelectorForEfficiencyVsPt");
  initGPselector(GpSelectorForEfficiencyVsVTXR, "GpSelectorForEfficiencyVsVTXR");
  initGPselector(GpSelectorForEfficiencyVsVTXZ, "GpSelectorForEfficiencyVsVTXZ");

  // SeedingLayerSets
  // If enabled, use last bin to denote other or unknown cases
  seedingLayerSetNames = pset.getParameter<std::vector<std::string>>("seedingLayerSets");
  std::vector<std::pair<SeedingLayerSetId, std::string>> stripPairSets;
  if (!seedingLayerSetNames.empty()) {
    std::vector<std::vector<std::string>> layerSets = SeedingLayerSetsBuilder::layerNamesInSets(seedingLayerSetNames);
    for (size_t i = 0; i < layerSets.size(); ++i) {
      const auto& layerSet = layerSets[i];
      if (layerSet.size() > std::tuple_size<SeedingLayerSetId>::value) {
        throw cms::Exception("Configuration")
            << "Got seedingLayerSet " << seedingLayerSetNames[i] << " with " << layerSet.size()
            << " elements, but I have a hard-coded maximum of " << std::tuple_size<SeedingLayerSetId>::value
            << ". Please increase the maximum in MTVHistoProducerAlgoForTracker.h";
      }
      SeedingLayerSetId setId;
      for (size_t j = 0; j < layerSet.size(); ++j) {
        // SeedingLayerSetsBuilder::fillDescriptions() kind-of
        // suggests that the 'M' prefix stands for strip mono hits
        // (maybe it should force), so making the assumption here is
        // (still) a bit ugly. But, this is the easiest way.
        bool isStripMono = !layerSet[j].empty() && layerSet[j][0] == 'M';
        setId[j] = std::make_tuple(SeedingLayerSetsBuilder::nameToEnumId(layerSet[j]), isStripMono);
      }
      // Account for the fact that strip triplet seeding may give pairs
      if (layerSet.size() == 3 && isTrackerStrip(std::get<GeomDetEnumerators::SubDetector>(std::get<0>(setId[0])))) {
        SeedingLayerSetId pairId;
        pairId[0] = setId[0];
        pairId[1] = setId[1];
        stripPairSets.emplace_back(pairId, layerSet[0] + "+" + layerSet[1]);
      }

      auto inserted = seedingLayerSetToBin.insert(std::make_pair(setId, i));
      if (!inserted.second)
        throw cms::Exception("Configuration") << "SeedingLayerSet " << seedingLayerSetNames[i]
                                              << " is specified twice, while the set list should be unique.";
    }

    // Add the "strip pairs from strip triplets" if they don't otherwise exist
    for (const auto& setIdName : stripPairSets) {
      auto inserted = seedingLayerSetToBin.insert(std::make_pair(setIdName.first, seedingLayerSetNames.size()));
      if (inserted.second)
        seedingLayerSetNames.push_back(setIdName.second);
    }

    seedingLayerSetNames.emplace_back("Other/Unknown");
  }

  // fix for the LogScale by Ryan
  if (useLogPt) {
    maxPt = log10(maxPt);
    if (minPt > 0) {
      minPt = log10(minPt);
    } else {
      edm::LogWarning("MultiTrackValidator")
          << "minPt = " << minPt << " <= 0 out of range while requesting log scale.  Using minPt = 0.1.";
      minPt = log10(0.1);
    }
  }
  if (useLogVertpos) {
    maxVertpos = std::log10(maxVertpos);
    if (minVertpos > 0) {
      minVertpos = std::log10(minVertpos);
    } else {
      edm::LogWarning("MultiTrackValidator")
          << "minVertpos = " << minVertpos << " <= 0 out of range while requesting log scale.  Using minVertpos = 0.1.";
      minVertpos = -1;
    }
  }
}

MTVHistoProducerAlgoForTracker::~MTVHistoProducerAlgoForTracker() {}

std::unique_ptr<RecoTrackSelectorBase> MTVHistoProducerAlgoForTracker::makeRecoTrackSelectorFromTPSelectorParameters(
    const edm::ParameterSet& pset) {
  edm::ParameterSet psetTrack;
  psetTrack.copyForModify(pset);
  psetTrack.eraseSimpleParameter("minHit");
  psetTrack.eraseSimpleParameter("minLayer");
  psetTrack.eraseSimpleParameter("signalOnly");
  psetTrack.eraseSimpleParameter("intimeOnly");
  psetTrack.eraseSimpleParameter("chargedOnly");
  psetTrack.eraseSimpleParameter("stableOnly");
  psetTrack.addParameter("maxChi2", 1e10);
  psetTrack.addParameter("minHit", 0);
  psetTrack.addParameter("minPixelHit", 0);
  psetTrack.addParameter("maxPixelHit", 99);
  psetTrack.addParameter("minLayer", 0);
  psetTrack.addParameter("min3DLayer", 0);
  psetTrack.addParameter("quality", std::vector<std::string>{});
  psetTrack.addParameter("algorithm", std::vector<std::string>{});
  psetTrack.addParameter("originalAlgorithm", std::vector<std::string>{});
  psetTrack.addParameter("algorithmMaskContains", std::vector<std::string>{});
  psetTrack.addParameter("invertRapidityCut", false);
  psetTrack.addParameter("minPhi", -3.2);
  psetTrack.addParameter("maxPhi", 3.2);
  return std::make_unique<RecoTrackSelectorBase>(psetTrack);
}

void MTVHistoProducerAlgoForTracker::pushbackNewMTVMonitoringBundles(Histograms& histograms) {
  // monitoring histograms for reco and sim tracks
  histograms.hs_eta.push_back(MTVMonitoringBundle());
  histograms.hs_pT.push_back(MTVMonitoringBundle());
  histograms.hs_pTvseta.push_back(MTVMonitoringBundle());
  histograms.hs_hit.push_back(MTVMonitoringBundle());
  histograms.hs_layer.push_back(MTVMonitoringBundle());
  histograms.hs_pixellayer.push_back(MTVMonitoringBundle());
  histograms.hs_3Dlayer.push_back(MTVMonitoringBundle());
  histograms.hs_pu.push_back(MTVMonitoringBundle());
  histograms.hs_phi.push_back(MTVMonitoringBundle());
  histograms.hs_dxy.push_back(MTVMonitoringBundle());
  histograms.hs_dz.push_back(MTVMonitoringBundle());
  histograms.hs_dxypv.push_back(MTVMonitoringBundle());
  histograms.hs_dzpv.push_back(MTVMonitoringBundle());
  histograms.hs_dxypvzoomed.push_back(MTVMonitoringBundle());
  histograms.hs_dzpvzoomed.push_back(MTVMonitoringBundle());
  histograms.hs_vertpos.push_back(MTVMonitoringBundle());
  histograms.hs_zpos.push_back(MTVMonitoringBundle());
  histograms.hs_dr.push_back(MTVMonitoringBundle());
  histograms.hs_drj.push_back(MTVMonitoringBundle());
  histograms.hs_dzpvcut.push_back(MTVMonitoringBundle());
  histograms.hs_dzpvsigcut.push_back(MTVMonitoringBundle());
  histograms.hs_simpvz.push_back(MTVMonitoringBundle());
  histograms.hs_chi2.push_back(MTVMonitoringBundle());
  histograms.hs_chi2prob.push_back(MTVMonitoringBundle());
  histograms.hs_seedingLayerSet.push_back(MTVMonitoringBundle());
  // monitoring histograms for resolutions
  histograms.hr_dxy.push_back(MTVResolutionBundle());
  histograms.hr_dz.push_back(MTVResolutionBundle());
  histograms.hr_phi.push_back(MTVResolutionBundle());
  histograms.hr_pt.push_back(MTVResolutionBundle());
  histograms.hr_cotTheta.push_back(MTVResolutionBundle());
  histograms.hr_theta.push_back(MTVResolutionBundle());
}

void MTVHistoProducerAlgoForTracker::bookSimHistos(DQMStore::IBooker& ibook, Histograms& histograms) {
  histograms.h_ptSIM = make1DIfLogX(ibook, useLogPt, "ptSIM", "generated p_{T}", nintPt, minPt, maxPt);
  histograms.h_etaSIM = ibook.book1D("etaSIM", "generated pseudorapidity #eta", nintEta, minEta, maxEta);
  histograms.h_phiSIM = ibook.book1D("phiSIM", "generated #phi", nintPhi, minPhi, maxPhi);
  histograms.h_tracksSIM =
      ibook.book1D("tracksSIM", "number of simulated tracks", nintTracks, minTracks, maxTracks * 10);
  histograms.h_vertposSIM =
      ibook.book1D("vertposSIM", "Transverse position of sim vertices", nintVertpos, minVertpos, maxVertpos);
  histograms.h_bunchxSIM = ibook.book1D("bunchxSIM", "bunch crossing", 21, -15.5, 5.5);
}

void MTVHistoProducerAlgoForTracker::bookSimAndRecoTrackHistos(DQMStore::IBooker& ibook,
                                                               Histograms& histograms,
                                                               const bool doSimTrackPlots,
                                                               const bool doRecoTrackPlots,
                                                               const bool doResolutionPlots) {
  if (doRecoTrackPlots) {
    histograms.h_tracks.push_back(
        ibook.book1D("tracks", "number of reconstructed tracks", nintTracks, minTracks, maxTracks));
    histograms.h_fakes.push_back(ibook.book1D("fakes", "number of fake reco tracks", nintTracks, minTracks, maxTracks));
    histograms.h_charge.push_back(ibook.book1D("charge", "charge", 3, -1.5, 1.5));

    histograms.h_hits.push_back(ibook.book1D("hits", "number of hits per track", nintHit, minHit, maxHit));
    histograms.h_losthits.push_back(ibook.book1D("losthits", "number of lost hits per track", nintHit, minHit, maxHit));
    histograms.h_nchi2.push_back(ibook.book1D("chi2", "normalized #chi^{2}", 200, 0, 20));
    histograms.h_nchi2_prob.push_back(ibook.book1D("chi2_prob", "normalized #chi^{2} probability", 100, 0, 1));

    histograms.h_nmisslayers_inner.push_back(
        ibook.book1D("missing_inner_layers", "number of missing inner layers", nintLayers, minLayers, maxLayers));
    histograms.h_nmisslayers_outer.push_back(
        ibook.book1D("missing_outer_layers", "number of missing outer layers", nintLayers, minLayers, maxLayers));

    histograms.h_algo.push_back(
        ibook.book1D("h_algo", "Tracks by algo", reco::TrackBase::algoSize, 0., double(reco::TrackBase::algoSize)));
    for (size_t ibin = 0; ibin < reco::TrackBase::algoSize - 1; ibin++)
      histograms.h_algo.back()->setBinLabel(ibin + 1, reco::TrackBase::algoNames[ibin]);
    histograms.h_algo.back()->disableAlphanumeric();
  }

  histograms.hs_eta.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "eta", "Pseudorapidity #eta", "", nintEta, minEta, maxEta);
  histograms.hs_pT.back().book1DIfLogX(
      ibook, useLogPt, doSimTrackPlots, doRecoTrackPlots, "pT", "p_{T}", "", nintPt, minPt, maxPt);
  histograms.hs_pTvseta.back().book2DIfLogY(ibook,
                                            useLogPt,
                                            doSimTrackPlots,
                                            doRecoTrackPlots,
                                            "pTvseta",
                                            "Pseudorapidity #eta",
                                            "p_{T}",
                                            nintEta,
                                            minEta,
                                            maxEta,
                                            nintPt,
                                            minPt,
                                            maxPt);

  histograms.hs_hit.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "hit", "Number of hits", "", nintHit, minHit, maxHit);
  histograms.hs_layer.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "layer", "Number of layers", "", nintLayers, minLayers, maxLayers);
  histograms.hs_pixellayer.back().book1D(ibook,
                                         doSimTrackPlots,
                                         doRecoTrackPlots,
                                         "pixellayer",
                                         "Number of pixel layers",
                                         "",
                                         nintLayers,
                                         minLayers,
                                         maxLayers);
  histograms.hs_3Dlayer.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "3Dlayer", "Number of 3D layers", "", nintLayers, minLayers, maxLayers);
  histograms.hs_pu.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "pu", "Number of Primary Vertices / Pileup", "", nintPu, minPu, maxPu);
  histograms.hs_phi.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "phi", "#phi angle", "", nintPhi, minPhi, maxPhi);
  histograms.hs_dxy.back().book1D(ibook,
                                  doSimTrackPlots,
                                  doRecoTrackPlots,
                                  "dxy",
                                  "Transverse impact parameter d_{xy} [cm]",
                                  "",
                                  nintDxy,
                                  minDxy,
                                  maxDxy);
  histograms.hs_dz.back().book1D(ibook,
                                 doSimTrackPlots,
                                 doRecoTrackPlots,
                                 "dz",
                                 "Longitudinal impact parameter d_{z} [cm]",
                                 "",
                                 nintDz,
                                 minDz,
                                 maxDz);
  histograms.hs_vertpos.back().book1DIfLogX(ibook,
                                            useLogVertpos,
                                            doSimTrackPlots,
                                            doRecoTrackPlots,
                                            "vertpos",
                                            "Radial displacement of production vertex [cm]",
                                            "",
                                            nintVertpos,
                                            minVertpos,
                                            maxVertpos);
  histograms.hs_zpos.back().book1D(ibook,
                                   doSimTrackPlots,
                                   doRecoTrackPlots,
                                   "zpos",
                                   "z coordinate of production vertex [cm]",
                                   "",
                                   nintZpos,
                                   minZpos,
                                   maxZpos);
  histograms.hs_dr.back().book1DLogX(
      ibook, doSimTrackPlots, doRecoTrackPlots, "dr", "dR", "", nintdr, log10(mindr), log10(maxdr));
  histograms.hs_drj.back().book1DLogX(
      ibook, doSimTrackPlots, doRecoTrackPlots, "drj", "dR(TP,jet)", "", nintdrj, log10(mindrj), log10(maxdrj));
  histograms.hs_simpvz.back().book1D(
      ibook, doSimTrackPlots, doRecoTrackPlots, "simpvz", "z of the simulated PV", "", nintPVz, minPVz, maxPVz);
  histograms.hs_chi2.back().book1D(
      ibook, false, doRecoTrackPlots, "chi2", "Normalized #chi^{2} / ndof", "", nintChi2, minChi2, maxChi2);
  histograms.hs_chi2prob.back().book1D(
      ibook, false, doRecoTrackPlots, "chi2prob", "Probability for given #chi^{2}", "", 100, 0., 1.);
  if (!seedingLayerSetNames.empty()) {
    const auto size = seedingLayerSetNames.size();
    histograms.hs_seedingLayerSet.back().book1D(
        ibook, false, doRecoTrackPlots, "seedingLayerSet", "Seeding layer set", "", size, 0, size);
    histograms.hs_seedingLayerSet.back().modifyHistograms(setBinLabels, seedingLayerSetNames);
  }

  if (doRecoTrackPlots) {
    auto bookResolutionPlots1D = [&](std::vector<dqm::reco::MonitorElement*>& vec, auto&&... params) {
      vec.push_back(doResolutionPlots ? ibook.book1D(std::forward<decltype(params)>(params)...) : nullptr);
    };
    auto bookResolutionPlots2D = [&](std::vector<dqm::reco::MonitorElement*>& vec, bool logx, auto&&... params) {
      vec.push_back(doResolutionPlots ? make2DIfLogX(ibook, logx, std::forward<decltype(params)>(params)...) : nullptr);
    };
    auto bookResolutionPlotsProfile2D = [&](std::vector<dqm::reco::MonitorElement*>& vec, auto&&... params) {
      vec.push_back(doResolutionPlots ? ibook.bookProfile2D(std::forward<decltype(params)>(params)...) : nullptr);
    };

    auto bookResolutionBundle = [&](std::vector<MTVResolutionBundle>& vec, auto&&... params) {
      vec.back().bookResolutions(ibook,
                                 nintEta,
                                 minEta,
                                 maxEta,
                                 nintPhi,
                                 minPhi,
                                 maxPhi,
                                 nintPt,
                                 minPt,
                                 maxPt,
                                 useLogPt,
                                 std::forward<decltype(params)>(params)...);
    };

    bookResolutionPlots1D(histograms.h_eta, "eta", "pseudorapidity residue", 1000, -0.1, 0.1);
    bookResolutionPlots1D(histograms.h_pt, "pullPt", "pull of p_{T}", 100, -10, 10);
    bookResolutionPlots1D(histograms.h_pullTheta, "pullTheta", "pull of #theta parameter", 250, -25, 25);
    bookResolutionPlots1D(histograms.h_pullPhi, "pullPhi", "pull of #phi parameter", 250, -25, 25);
    bookResolutionPlots1D(histograms.h_pullDxy, "pullDxy", "pull of d_{xy} parameter", 250, -25, 25);
    bookResolutionPlots1D(histograms.h_pullDz, "pullDz", "pull of d_{z} parameter", 250, -25, 25);
    bookResolutionPlots1D(histograms.h_pullQoverp, "pullQoverp", "pull of qoverp parameter", 250, -25, 25);

    /* TO BE FIXED -----------
    if (associators[ww]=="TrackAssociatorByChi2"){
      histograms.h_assocSimToReco_chi2.push_back( ibook.book1D("assocChi2","track association #chi^{2}",1000000,0,100000) );
      histograms.h_assocSimToReco_chi2_prob.push_back(ibook.book1D("assocChi2_prob","probability of association #chi^{2}",100,0,1));
    } else if (associators[ww]=="quickTrackAssociatorByHits"){
      histograms.h_assocSimToReco_Fraction.push_back( ibook.book1D("assocFraction","fraction of shared hits",200,0,2) );
      histograms.h_assocSimToReco_SharedHit.push_back(ibook.book1D("assocSharedHit","number of shared hits",20,0,20));
    }
    */
    histograms.h_assocSimToReco_Fraction.push_back(ibook.book1D("assocFraction", "fraction of shared hits", 200, 0, 2));
    histograms.h_assocSimToReco_SharedHit.push_back(
        ibook.book1D("assocSharedHit", "number of shared hits", 41, -0.5, 40.5));
    // ----------------------

    // use the standard error of the mean as the errors in the profile
    histograms.chi2_vs_nhits.push_back(
        ibook.bookProfile("chi2mean_vs_nhits", "mean #chi^{2} vs nhits", nintHit, minHit, maxHit, 100, 0, 10, " "));

    bookResolutionPlots2D(
        histograms.etares_vs_eta, false, "etares_vs_eta", "etaresidue vs eta", nintEta, minEta, maxEta, 200, -0.1, 0.1);
    bookResolutionPlots2D(
        histograms.nrec_vs_nsim,
        false,
        "nrec_vs_nsim",
        "Number of selected reco tracks vs. number of selected sim tracks;TrackingParticles;Reco tracks",
        nintTracks,
        minTracks,
        maxTracks,
        nintTracks,
        minTracks,
        maxTracks);

    histograms.chi2_vs_eta.push_back(
        ibook.bookProfile("chi2mean", "mean #chi^{2} vs #eta", nintEta, minEta, maxEta, 200, 0, 20, " "));
    histograms.chi2_vs_phi.push_back(
        ibook.bookProfile("chi2mean_vs_phi", "mean #chi^{2} vs #phi", nintPhi, minPhi, maxPhi, 200, 0, 20, " "));
    histograms.chi2_vs_pt.push_back(
        makeProfileIfLogX(ibook, useLogPt, "chi2mean_vs_pt", "mean #chi^{2} vs p_{T}", nintPt, minPt, maxPt, 0, 20));
    histograms.chi2_vs_drj.push_back(makeProfileIfLogX(
        ibook, true, "chi2mean_vs_drj", "mean #chi^{2} vs dR(track,jet)", nintdrj, log10(mindrj), log10(maxdrj), 0, 20));

    histograms.assoc_chi2_vs_eta.push_back(
        ibook.bookProfile("assoc_chi2mean", "mean #chi^{2} vs #eta", nintEta, minEta, maxEta, 200, 0., 20., " "));
    histograms.assoc_chi2prob_vs_eta.push_back(ibook.bookProfile(
        "assoc_chi2prob_vs_eta", "mean #chi^{2} probability vs #eta", nintEta, minEta, maxEta, 100, 0., 1., " "));
    histograms.assoc_chi2_vs_pt.push_back(makeProfileIfLogX(
        ibook, useLogPt, "assoc_chi2mean_vs_pt", "mean #chi^{2} vs p_{T}", nintPt, minPt, maxPt, 0., 20.));
    histograms.assoc_chi2prob_vs_pt.push_back(makeProfileIfLogX(
        ibook, useLogPt, "assoc_chi2prob_vs_pt", "mean #chi^{2} probability vs p_{T}", nintPt, minPt, maxPt, 0., 1.));
    histograms.assoc_chi2_vs_drj.push_back(makeProfileIfLogX(ibook,
                                                             true,
                                                             "assoc_chi2mean_vs_drj",
                                                             "mean #chi^{2} vs dR(track,jet)",
                                                             nintdrj,
                                                             log10(mindrj),
                                                             log10(maxdrj),
                                                             0.,
                                                             20));
    histograms.assoc_chi2prob_vs_drj.push_back(makeProfileIfLogX(ibook,
                                                                 true,
                                                                 "assoc_chi2prob_vs_drj",
                                                                 "mean #chi^{2} probability vs dR(track,jet)",
                                                                 nintdrj,
                                                                 log10(mindrj),
                                                                 log10(maxdrj),
                                                                 0.,
                                                                 1.));

    histograms.nhits_vs_eta.push_back(
        ibook.bookProfile("hits_eta", "mean hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nPXBhits_vs_eta.push_back(ibook.bookProfile(
        "PXBhits_vs_eta", "mean # PXB its vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nPXFhits_vs_eta.push_back(ibook.bookProfile(
        "PXFhits_vs_eta", "mean # PXF hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nPXLhits_vs_eta.push_back(ibook.bookProfile(
        "PXLhits_vs_eta", "mean # PXL hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nTIBhits_vs_eta.push_back(ibook.bookProfile(
        "TIBhits_vs_eta", "mean # TIB hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nTIDhits_vs_eta.push_back(ibook.bookProfile(
        "TIDhits_vs_eta", "mean # TID hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nTOBhits_vs_eta.push_back(ibook.bookProfile(
        "TOBhits_vs_eta", "mean # TOB hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nTEChits_vs_eta.push_back(ibook.bookProfile(
        "TEChits_vs_eta", "mean # TEC hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    histograms.nSTRIPhits_vs_eta.push_back(ibook.bookProfile(
        "STRIPhits_vs_eta", "mean # STRIP hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));

    histograms.nLayersWithMeas_vs_eta.push_back(ibook.bookProfile("LayersWithMeas_eta",
                                                                  "mean # Layers with measurement vs eta",
                                                                  nintEta,
                                                                  minEta,
                                                                  maxEta,
                                                                  nintLayers,
                                                                  minLayers,
                                                                  maxLayers,
                                                                  " "));
    histograms.nPXLlayersWithMeas_vs_eta.push_back(ibook.bookProfile("PXLlayersWithMeas_vs_eta",
                                                                     "mean # PXL Layers with measurement vs eta",
                                                                     nintEta,
                                                                     minEta,
                                                                     maxEta,
                                                                     nintLayers,
                                                                     minLayers,
                                                                     maxLayers,
                                                                     " "));
    histograms.nSTRIPlayersWithMeas_vs_eta.push_back(ibook.bookProfile("STRIPlayersWithMeas_vs_eta",
                                                                       "mean # STRIP Layers with measurement vs eta",
                                                                       nintEta,
                                                                       minEta,
                                                                       maxEta,
                                                                       nintLayers,
                                                                       minLayers,
                                                                       maxLayers,
                                                                       " "));
    histograms.nSTRIPlayersWith1dMeas_vs_eta.push_back(
        ibook.bookProfile("STRIPlayersWith1dMeas_vs_eta",
                          "mean # STRIP Layers with 1D measurement vs eta",
                          nintEta,
                          minEta,
                          maxEta,
                          nintLayers,
                          minLayers,
                          maxLayers,
                          " "));
    histograms.nSTRIPlayersWith2dMeas_vs_eta.push_back(
        ibook.bookProfile("STRIPlayersWith2dMeas_vs_eta",
                          "mean # STRIP Layers with 2D measurement vs eta",
                          nintEta,
                          minEta,
                          maxEta,
                          nintLayers,
                          minLayers,
                          maxLayers,
                          " "));

    if (doMTDPlots_) {
      histograms.nMTDhits_vs_eta.push_back(ibook.bookProfile(
          "MTDhits_vs_eta", "mean # MTD hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));

      histograms.nBTLhits_vs_eta.push_back(ibook.bookProfile(
          "BTLhits_vs_eta", "mean # BTL hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));

      histograms.nETLhits_vs_eta.push_back(ibook.bookProfile(
          "ETLhits_vs_eta", "mean # ETL hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));
    }

    histograms.nhits_vs_phi.push_back(
        ibook.bookProfile("hits_phi", "mean # hits vs #phi", nintPhi, minPhi, maxPhi, nintHit, minHit, maxHit, " "));

    histograms.nlosthits_vs_eta.push_back(ibook.bookProfile(
        "losthits_vs_eta", "mean # lost hits vs eta", nintEta, minEta, maxEta, nintHit, minHit, maxHit, " "));

    //resolution of track parameters
    //                       dPt/Pt    cotTheta        Phi            TIP            LIP
    // log10(pt)<0.5        100,0.1    240,0.08     100,0.015      100,0.1000    150,0.3000
    // 0.5<log10(pt)<1.5    100,0.1    120,0.01     100,0.003      100,0.0100    150,0.0500
    // >1.5                 100,0.3    100,0.005    100,0.0008     100,0.0060    120,0.0300

    if (doResolutionPlots) {
      bookResolutionBundle(histograms.hr_pt, "pt", ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax);
      bookResolutionBundle(histograms.hr_phi, "phi", phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax);
      bookResolutionBundle(histograms.hr_dxy, "dxy", dxyRes_nbin, dxyRes_rangeMin, dxyRes_rangeMax);
      bookResolutionBundle(histograms.hr_dz, "dz", dzRes_nbin, dzRes_rangeMin, dzRes_rangeMax);
      // FIXME: those are random values for theta here. Check the histograms and reset them!!!
      bookResolutionBundle(histograms.hr_theta, "theta", 100, -1, 1);
      bookResolutionBundle(
          histograms.hr_cotTheta, "cotTheta", cotThetaRes_nbin, cotThetaRes_rangeMin, cotThetaRes_rangeMax);
    }

    bookResolutionPlotsProfile2D(histograms.ptmean_vs_eta_phi,
                                 "ptmean_vs_eta_phi",
                                 "mean p_{t} vs #eta and #phi",
                                 nintPhi,
                                 minPhi,
                                 maxPhi,
                                 nintEta,
                                 minEta,
                                 maxEta,
                                 1000,
                                 0,
                                 1000);
    bookResolutionPlotsProfile2D(histograms.phimean_vs_eta_phi,
                                 "phimean_vs_eta_phi",
                                 "mean #phi vs #eta and #phi",
                                 nintPhi,
                                 minPhi,
                                 maxPhi,
                                 nintEta,
                                 minEta,
                                 maxEta,
                                 nintPhi,
                                 minPhi,
                                 maxPhi);

    //      histograms.h_ptshiftetamean.push_back( ibook.book1D("h_ptshifteta_Mean","<#deltapT/pT>[%] vs #eta",nintEta,minEta,maxEta) );

    bookResolutionPlots2D(histograms.nrecHit_vs_nsimHit_rec2sim,
                          false,
                          "nrecHit_vs_nsimHit_rec2sim",
                          "nrecHit vs nsimHit (Rec2simAssoc)",
                          nintHit,
                          minHit,
                          maxHit,
                          nintHit,
                          minHit,
                          maxHit);
  }

  if (doSimTrackPlots) {
    histograms.nrecHit_vs_nsimHit_sim2rec.push_back(doResolutionPlots
                                                        ? ibook.book2D("nrecHit_vs_nsimHit_sim2rec",
                                                                       "nrecHit vs nsimHit (Sim2RecAssoc)",
                                                                       nintHit,
                                                                       minHit,
                                                                       maxHit,
                                                                       nintHit,
                                                                       minHit,
                                                                       maxHit)
                                                        : nullptr);

    // TODO: use the dynamic track algo priority order also here
    constexpr auto nalgos = reco::TrackBase::algoSize;
    histograms.h_duplicates_oriAlgo_vs_oriAlgo.push_back(ibook.book2D("duplicates_oriAlgo_vs_oriAlgo",
                                                                      "Duplicate tracks: originalAlgo vs originalAlgo",
                                                                      nalgos,
                                                                      0,
                                                                      nalgos,
                                                                      nalgos,
                                                                      0,
                                                                      nalgos));
    setBinLabelsAlgo(histograms.h_duplicates_oriAlgo_vs_oriAlgo.back(), 1);
    setBinLabelsAlgo(histograms.h_duplicates_oriAlgo_vs_oriAlgo.back(), 2);
  }
}

void MTVHistoProducerAlgoForTracker::bookSimTrackPVAssociationHistos(DQMStore::IBooker& ibook, Histograms& histograms) {
  histograms.hs_dxypv.back().book1D(
      ibook, true, false, "dxypv", "Transverse impact parameter d_{xy} wrt PV", "", nintDxy, minDxy, maxDxy);
  histograms.hs_dzpv.back().book1D(
      ibook, true, false, "dzpv", "Longitudinal impact parameter d_{z} wrt PV", "", nintDz, minDz, maxDz);
  histograms.hs_dxypvzoomed.back().book1D(ibook,
                                          true,
                                          false,
                                          "dxypvzoomed",
                                          "Transverse impact parameter d_{xy} wrt P",
                                          "",
                                          nintDxy,
                                          minDxy / dxyDzZoom,
                                          maxDxy / dxyDzZoom);
  histograms.hs_dzpvzoomed.back().book1D(ibook,
                                         true,
                                         false,
                                         "dzpvzoomed",
                                         "Longitudinal impact parameter d_{z} wrt PV",
                                         "",
                                         nintDz,
                                         minDz / dxyDzZoom,
                                         maxDz / dxyDzZoom);

  if (doDzPVcutPlots_) {
    histograms.hs_dzpvcut.back().book1D(
        ibook, true, false, "dzpvcut", "Longitudinal impact parameter d_{z} wrt PV", "", nintDzpvCum, 0, maxDzpvCum);
    histograms.h_simul2_dzpvcut.push_back(ibook.book1D("num_simul2_dzpvcut",
                                                       "N of simulated tracks (associated to any track) from sim PV",
                                                       nintDzpvCum,
                                                       0,
                                                       maxDzpvCum));
  }
}

void MTVHistoProducerAlgoForTracker::bookRecoPVAssociationHistos(DQMStore::IBooker& ibook, Histograms& histograms) {
  histograms.hs_dxypv.back().book1D(
      ibook, false, true, "dxypv", "Transverse impact parameter dxy wrt PV", "", nintDxy, minDxy, maxDxy);
  histograms.hs_dzpv.back().book1D(
      ibook, false, true, "dzpv", "Longitudinal impact parameter dz wrt PV", "", nintDz, minDz, maxDz);
  histograms.hs_dxypvzoomed.back().book1D(ibook,
                                          true,
                                          false,
                                          "dxypvzoomed",
                                          "Transverse impact parameter dxy wrt P",
                                          "",
                                          nintDxy,
                                          minDxy / dxyDzZoom,
                                          maxDxy / dxyDzZoom);
  histograms.hs_dzpvzoomed.back().book1D(ibook,
                                         true,
                                         false,
                                         "dzpvzoomed",
                                         "Longitudinal impact parameter dz wrt PV",
                                         "",
                                         nintDz,
                                         minDz / dxyDzZoom,
                                         maxDz / dxyDzZoom);

  if (doDzPVcutPlots_) {
    histograms.hs_dzpvcut.back().book1D(
        ibook, false, true, "dzpvcut", "Longitudinal impact parameter dz wrt PV", "", nintDzpvCum, 0, maxDzpvCum);
    histograms.h_simul2_dzpvcut.push_back(ibook.book1D("num_simul2_dzpvcut",
                                                       "N of simulated tracks (associated to any track) from sim PV",
                                                       nintDzpvCum,
                                                       0,
                                                       maxDzpvCum));
  }
}

void MTVHistoProducerAlgoForTracker::bookRecodEdxHistos(DQMStore::IBooker& ibook, Histograms& histograms) {
  // dE/dx stuff
  histograms.h_dedx_estim.emplace_back();
  histograms.h_dedx_estim.back().push_back(
      ibook.book1D("h_dedx_estim1", "dE/dx estimator 1", nintDeDx, minDeDx, maxDeDx));
  histograms.h_dedx_estim.back().push_back(
      ibook.book1D("h_dedx_estim2", "dE/dx estimator 2", nintDeDx, minDeDx, maxDeDx));

  histograms.h_dedx_nom.emplace_back();
  histograms.h_dedx_nom.back().push_back(
      ibook.book1D("h_dedx_nom1", "dE/dx number of measurements", nintHit, minHit, maxHit));
  histograms.h_dedx_nom.back().push_back(
      ibook.book1D("h_dedx_nom2", "dE/dx number of measurements", nintHit, minHit, maxHit));

  histograms.h_dedx_sat.emplace_back();
  histograms.h_dedx_sat.back().push_back(
      ibook.book1D("h_dedx_sat1", "dE/dx number of measurements with saturation", nintHit, minHit, maxHit));
  histograms.h_dedx_sat.back().push_back(
      ibook.book1D("h_dedx_sat2", "dE/dx number of measurements with saturation", nintHit, minHit, maxHit));
}

void MTVHistoProducerAlgoForTracker::bookSeedHistos(DQMStore::IBooker& ibook, Histograms& histograms) {
  histograms.h_seedsFitFailed.push_back(
      ibook.book1D("seeds_fitFailed", "Number of seeds for which the fit failed", nintTracks, minTracks, maxTracks));
  histograms.h_seedsFitFailedFraction.push_back(
      ibook.book1D("seeds_fitFailedFraction", "Fraction of seeds for which the fit failed", 100, 0, 1));
}

void MTVHistoProducerAlgoForTracker::bookMVAHistos(DQMStore::IBooker& ibook, Histograms& histograms, size_t nMVAs) {
  histograms.h_reco_mva.emplace_back();
  histograms.h_assocRecoToSim_mva.emplace_back();

  histograms.h_reco_mvacut.emplace_back();
  histograms.h_assocSimToReco_mvacut.emplace_back();
  histograms.h_assocRecoToSim_mvacut.emplace_back();
  histograms.h_simul2_mvacut.emplace_back();

  histograms.h_reco_mva_hp.emplace_back();
  histograms.h_assocRecoToSim_mva_hp.emplace_back();

  histograms.h_reco_mvacut_hp.emplace_back();
  histograms.h_assocSimToReco_mvacut_hp.emplace_back();
  histograms.h_assocRecoToSim_mvacut_hp.emplace_back();
  histograms.h_simul2_mvacut_hp.emplace_back();

  histograms.h_assocRecoToSim_mva_vs_pt.emplace_back();
  histograms.h_fake_mva_vs_pt.emplace_back();
  histograms.h_assocRecoToSim_mva_vs_pt_hp.emplace_back();
  histograms.h_fake_mva_vs_pt_hp.emplace_back();
  histograms.h_assocRecoToSim_mva_vs_eta.emplace_back();
  histograms.h_fake_mva_vs_eta.emplace_back();
  histograms.h_assocRecoToSim_mva_vs_eta_hp.emplace_back();
  histograms.h_fake_mva_vs_eta_hp.emplace_back();

  for (size_t i = 1; i <= nMVAs; ++i) {
    auto istr = std::to_string(i);
    std::string pfix;

    if (i == 1) {
      histograms.h_reco_mva_hp.back().emplace_back();
      histograms.h_assocRecoToSim_mva_hp.back().emplace_back();

      histograms.h_reco_mvacut_hp.back().emplace_back();
      histograms.h_assocSimToReco_mvacut_hp.back().emplace_back();
      histograms.h_assocRecoToSim_mvacut_hp.back().emplace_back();
      histograms.h_simul2_mvacut_hp.back().emplace_back();

      histograms.h_assocRecoToSim_mva_vs_pt_hp.back().emplace_back();
      histograms.h_fake_mva_vs_pt_hp.back().emplace_back();
      histograms.h_assocRecoToSim_mva_vs_eta_hp.back().emplace_back();
      histograms.h_fake_mva_vs_eta_hp.back().emplace_back();
    } else {
      pfix = " (not loose-selected)";
      std::string pfix2 = " (not HP-selected)";

      histograms.h_reco_mva_hp.back().push_back(ibook.book1D(
          "num_reco_mva" + istr + "_hp", "N of reco track after vs MVA" + istr + pfix2, nintMVA, minMVA, maxMVA));
      histograms.h_assocRecoToSim_mva_hp.back().push_back(
          ibook.book1D("num_assoc(recoToSim)_mva" + istr + "_hp",
                       "N of associated tracks (recoToSim) vs MVA" + istr + pfix2,
                       nintMVA,
                       minMVA,
                       maxMVA));

      histograms.h_reco_mvacut_hp.back().push_back(ibook.book1D("num_reco_mva" + istr + "cut" + "_hp",
                                                                "N of reco track vs cut on MVA" + istr + pfix2,
                                                                nintMVA,
                                                                minMVA,
                                                                maxMVA));
      histograms.h_assocSimToReco_mvacut_hp.back().push_back(
          ibook.book1D("num_assoc(simToReco)_mva" + istr + "cut_hp",
                       "N of associated tracks (simToReco) vs cut on MVA" + istr + pfix2,
                       nintMVA,
                       minMVA,
                       maxMVA));
      histograms.h_assocRecoToSim_mvacut_hp.back().push_back(
          ibook.book1D("num_assoc(recoToSim)_mva" + istr + "cut_hp",
                       "N of associated tracks (recoToSim) vs cut on MVA" + istr + pfix2,
                       nintMVA,
                       minMVA,
                       maxMVA));
      histograms.h_simul2_mvacut_hp.back().push_back(
          ibook.book1D("num_simul2_mva" + istr + "cut_hp",
                       "N of simulated tracks (associated to any track) vs cut on MVA" + istr + pfix2,
                       nintMVA,
                       minMVA,
                       maxMVA));

      histograms.h_assocRecoToSim_mva_vs_pt_hp.back().push_back(
          makeProfileIfLogX(ibook,
                            useLogPt,
                            ("mva_assoc(recoToSim)_mva" + istr + "_pT_hp").c_str(),
                            ("MVA" + istr + " of associated tracks (recoToSim) vs. track p_{T}" + pfix2).c_str(),
                            nintPt,
                            minPt,
                            maxPt,
                            minMVA,
                            maxMVA));
      histograms.h_fake_mva_vs_pt_hp.back().push_back(
          makeProfileIfLogX(ibook,
                            useLogPt,
                            ("mva_fake_mva" + istr + "pT_hp").c_str(),
                            ("MVA" + istr + " of non-associated tracks (recoToSim) vs. track p_{T}" + pfix2).c_str(),
                            nintPt,
                            minPt,
                            maxPt,
                            minMVA,
                            maxMVA));
      histograms.h_assocRecoToSim_mva_vs_eta_hp.back().push_back(
          ibook.bookProfile("mva_assoc(recoToSim)_mva" + istr + "_eta_hp",
                            "MVA" + istr + " of associated tracks (recoToSim) vs. track #eta" + pfix2,
                            nintEta,
                            minEta,
                            maxEta,
                            nintMVA,
                            minMVA,
                            maxMVA));
      histograms.h_fake_mva_vs_eta_hp.back().push_back(
          ibook.bookProfile("mva_fake_mva" + istr + "eta_hp",
                            "MVA" + istr + " of non-associated tracks (recoToSim) vs. track #eta" + pfix2,
                            nintEta,
                            minEta,
                            maxEta,
                            nintMVA,
                            minMVA,
                            maxMVA));
    }

    histograms.h_reco_mva.back().push_back(
        ibook.book1D("num_reco_mva" + istr, "N of reco track vs MVA" + istr + pfix, nintMVA, minMVA, maxMVA));
    histograms.h_assocRecoToSim_mva.back().push_back(
        ibook.book1D("num_assoc(recoToSim)_mva" + istr,
                     "N of associated tracks (recoToSim) vs MVA" + istr + pfix,
                     nintMVA,
                     minMVA,
                     maxMVA));

    histograms.h_reco_mvacut.back().push_back(ibook.book1D(
        "num_reco_mva" + istr + "cut", "N of reco track vs cut on MVA" + istr + pfix, nintMVA, minMVA, maxMVA));
    histograms.h_assocSimToReco_mvacut.back().push_back(
        ibook.book1D("num_assoc(simToReco)_mva" + istr + "cut",
                     "N of associated tracks (simToReco) vs cut on MVA" + istr + pfix,
                     nintMVA,
                     minMVA,
                     maxMVA));
    histograms.h_assocRecoToSim_mvacut.back().push_back(
        ibook.book1D("num_assoc(recoToSim)_mva" + istr + "cut",
                     "N of associated tracks (recoToSim) vs cut on MVA" + istr + pfix,
                     nintMVA,
                     minMVA,
                     maxMVA));
    histograms.h_simul2_mvacut.back().push_back(
        ibook.book1D("num_simul2_mva" + istr + "cut",
                     "N of simulated tracks (associated to any track) vs cut on MVA" + istr + pfix,
                     nintMVA,
                     minMVA,
                     maxMVA));

    histograms.h_assocRecoToSim_mva_vs_pt.back().push_back(
        makeProfileIfLogX(ibook,
                          useLogPt,
                          ("mva_assoc(recoToSim)_mva" + istr + "_pT").c_str(),
                          ("MVA" + istr + " of associated tracks (recoToSim) vs. track p_{T}" + pfix).c_str(),
                          nintPt,
                          minPt,
                          maxPt,
                          minMVA,
                          maxMVA));
    histograms.h_fake_mva_vs_pt.back().push_back(
        makeProfileIfLogX(ibook,
                          useLogPt,
                          ("mva_fake_mva" + istr + "_pT").c_str(),
                          ("MVA" + istr + " of non-associated tracks (recoToSim) vs. track p_{T}" + pfix).c_str(),
                          nintPt,
                          minPt,
                          maxPt,
                          minMVA,
                          maxMVA));
    histograms.h_assocRecoToSim_mva_vs_eta.back().push_back(
        ibook.bookProfile("mva_assoc(recoToSim)_mva" + istr + "_eta",
                          "MVA" + istr + " of associated tracks (recoToSim) vs. track #eta" + pfix,
                          nintEta,
                          minEta,
                          maxEta,
                          nintMVA,
                          minMVA,
                          maxMVA));
    histograms.h_fake_mva_vs_eta.back().push_back(
        ibook.bookProfile("mva_fake_mva" + istr + "_eta",
                          "MVA" + istr + " of non-associated tracks (recoToSim) vs. track #eta" + pfix,
                          nintEta,
                          minEta,
                          maxEta,
                          nintMVA,
                          minMVA,
                          maxMVA));
  }
}

void MTVHistoProducerAlgoForTracker::fill_generic_simTrack_histos(const Histograms& histograms,
                                                                  const TrackingParticle::Vector& momentumTP,
                                                                  const TrackingParticle::Point& vertexTP,
                                                                  int bx) const {
  if (bx == 0) {
    histograms.h_ptSIM->Fill(sqrt(momentumTP.perp2()));
    histograms.h_etaSIM->Fill(momentumTP.eta());
    histograms.h_phiSIM->Fill(momentumTP.phi());
    histograms.h_vertposSIM->Fill(sqrt(vertexTP.perp2()));
  }
  histograms.h_bunchxSIM->Fill(bx);
}

void MTVHistoProducerAlgoForTracker::fill_recoAssociated_simTrack_histos(
    const Histograms& histograms,
    int count,
    const TrackingParticle& tp,
    const TrackingParticle::Vector& momentumTP,
    const TrackingParticle::Point& vertexTP,
    double dxySim,
    double dzSim,
    double dxyPVSim,
    double dzPVSim,
    int nSimHits,
    int nSimLayers,
    int nSimPixelLayers,
    int nSimStripMonoAndStereoLayers,
    const reco::Track* track,
    int numVertices,
    double dR,
    double dRJet,
    const math::XYZPoint* pvPosition,
    const TrackingVertex::LorentzVector* simPVPosition,
    const math::XYZPoint& bsPosition,
    const std::vector<float>& mvas,
    unsigned int selectsLoose,
    unsigned int selectsHP) const {
  const bool isMatched = track;
  const bool isReconstructable = (*TpSelectorForTechnicalEfficiency)(tp);
  const auto eta = getEta(momentumTP.eta());
  const auto phi = momentumTP.phi();
  const auto pt = getPt(sqrt(momentumTP.perp2()));
  const auto nSim3DLayers = nSimPixelLayers + nSimStripMonoAndStereoLayers;

  const auto vertexTPwrtBS = vertexTP - bsPosition;
  const auto vertxy = std::sqrt(vertexTPwrtBS.perp2());
  const auto vertz = vertexTPwrtBS.z();

  //efficiency vs. cut on MVA
  //
  // Note that this includes also pileup TPs, as "signalOnly"
  // selection is applied only in the TpSelector*. Have to think if
  // this is really what we want.
  if (isMatched) {
    for (size_t i = 0; i < mvas.size(); ++i) {
      if (i <= selectsLoose) {
        histograms.h_simul2_mvacut[count][i]->Fill(maxMVA);
        histograms.h_assocSimToReco_mvacut[count][i]->Fill(mvas[i]);
      }
      if (i >= 1 && i <= selectsHP) {
        histograms.h_simul2_mvacut_hp[count][i]->Fill(maxMVA);
        histograms.h_assocSimToReco_mvacut_hp[count][i]->Fill(mvas[i]);
      }
    }
  }

  if ((*TpSelectorForEfficiencyVsEta)(tp))
    //effic vs eta
    histograms.hs_eta[count].fillSimTrackHistos(isMatched, isReconstructable, eta);

  if ((*TpSelectorForEfficiencyVsPhi)(tp)) {
    histograms.hs_phi[count].fillSimTrackHistos(isMatched, isReconstructable, phi);
    histograms.hs_hit[count].fillSimTrackHistos(isMatched, isReconstructable, nSimHits);
    histograms.hs_layer[count].fillSimTrackHistos(isMatched, isReconstructable, nSimLayers);
    histograms.hs_pixellayer[count].fillSimTrackHistos(isMatched, isReconstructable, nSimPixelLayers);
    histograms.hs_3Dlayer[count].fillSimTrackHistos(isMatched, isReconstructable, nSim3DLayers);
    histograms.hs_pu[count].fillSimTrackHistos(isMatched, isReconstructable, numVertices);
    histograms.hs_dr[count].fillSimTrackHistos(isMatched, isReconstructable, dR);
    histograms.hs_drj[count].fillSimTrackHistos(isMatched, isReconstructable, dRJet);

    if (isMatched && histograms.nrecHit_vs_nsimHit_sim2rec[count])
      histograms.nrecHit_vs_nsimHit_sim2rec[count]->Fill(track->numberOfValidHits(), nSimHits);
  }

  if ((*TpSelectorForEfficiencyVsPt)(tp))
    histograms.hs_pT[count].fillSimTrackHistos(isMatched, isReconstructable, pt);

  if ((*TpSelectorForEfficiencyVsVTXR)(tp)) {
    histograms.hs_vertpos[count].fillSimTrackHistos(isMatched, isReconstructable, vertxy);
    histograms.hs_dxy[count].fillSimTrackHistos(isMatched, isReconstructable, dxySim);
    if (pvPosition) {
      histograms.hs_dxypv[count].fillSimTrackHistos(isMatched, isReconstructable, dxyPVSim);
      histograms.hs_dxypvzoomed[count].fillSimTrackHistos(isMatched, isReconstructable, dxyPVSim);
    }
  }

  if ((*TpSelectorForEfficiencyVsVTXZ)(tp)) {
    histograms.hs_zpos[count].fillSimTrackHistos(isMatched, isReconstructable, vertz);
    histograms.hs_dz[count].fillSimTrackHistos(isMatched, isReconstructable, dzSim);
    if (pvPosition) {
      histograms.hs_dzpv[count].fillSimTrackHistos(isMatched, isReconstructable, dzPVSim);
      histograms.hs_dzpvzoomed[count].fillSimTrackHistos(isMatched, isReconstructable, dzPVSim);
      if (doDzPVcutPlots_) {
        const double dzpvcut = std::abs(track->dz(*pvPosition));
        histograms.hs_dzpvcut[count].fillSimTrackHistos(isMatched, isReconstructable, dzpvcut);
        if (isMatched)
          histograms.h_simul2_dzpvcut[count]->Fill(dzpvcut);
      }
    }
    if (simPVPosition) {
      const auto simpvz = simPVPosition->z();
      histograms.hs_simpvz[count].fillSimTrackHistos(isMatched, isReconstructable, simpvz);
    }
  }
}

void MTVHistoProducerAlgoForTracker::fill_duplicate_histos(const Histograms& histograms,
                                                           int count,
                                                           const reco::Track& track1,
                                                           const reco::Track& track2) const {
  histograms.h_duplicates_oriAlgo_vs_oriAlgo[count]->Fill(track1.originalAlgo(), track2.originalAlgo());
}

void MTVHistoProducerAlgoForTracker::fill_simTrackBased_histos(const Histograms& histograms, int numSimTracks) const {
  histograms.h_tracksSIM->Fill(numSimTracks);
}

// dE/dx
void MTVHistoProducerAlgoForTracker::fill_dedx_recoTrack_histos(
    const Histograms& histograms,
    int count,
    const edm::RefToBase<reco::Track>& trackref,
    const std::vector<const edm::ValueMap<reco::DeDxData>*>& v_dEdx) const {
  for (unsigned int i = 0; i < v_dEdx.size(); i++) {
    const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[i]);
    const reco::DeDxData& dedx = dEdxTrack[trackref];
    histograms.h_dedx_estim[count][i]->Fill(dedx.dEdx());
    histograms.h_dedx_nom[count][i]->Fill(dedx.numberOfMeasurements());
    histograms.h_dedx_sat[count][i]->Fill(dedx.numberOfSaturatedMeasurements());
  }
}

void MTVHistoProducerAlgoForTracker::fill_generic_recoTrack_histos(const Histograms& histograms,
                                                                   int count,
                                                                   const reco::Track& track,
                                                                   const TrackerTopology& ttopo,
                                                                   const math::XYZPoint& bsPosition,
                                                                   const math::XYZPoint* pvPosition,
                                                                   const TrackingVertex::LorentzVector* simPVPosition,
                                                                   bool isMatched,
                                                                   bool isSigMatched,
                                                                   bool isChargeMatched,
                                                                   int numAssocRecoTracks,
                                                                   int numVertices,
                                                                   int nSimHits,
                                                                   double sharedFraction,
                                                                   double dR,
                                                                   double dRJet,
                                                                   const std::vector<float>& mvas,
                                                                   unsigned int selectsLoose,
                                                                   unsigned int selectsHP) const {
  //Fill track algo histogram
  histograms.h_algo[count]->Fill(track.algo());
  int sharedHits = sharedFraction * track.numberOfValidHits();

  //Compute fake rate vs eta
  const auto eta = getEta(track.momentum().eta());
  const auto phi = track.momentum().phi();
  const auto pt = getPt(sqrt(track.momentum().perp2()));
  const auto dxy = track.dxy(bsPosition);
  const auto dz = track.dz(bsPosition);
  const auto dxypv = pvPosition ? track.dxy(*pvPosition) : 0.0;
  const auto dzpv = pvPosition ? track.dz(*pvPosition) : 0.0;
  const auto nhits = track.found();
  const auto nlayers = track.hitPattern().trackerLayersWithMeasurement();
  const auto nPixelLayers = track.hitPattern().pixelLayersWithMeasurement();
  const auto n3DLayers = nPixelLayers + track.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
  const auto refPointWrtBS = track.referencePoint() - bsPosition;
  const auto vertxy = std::sqrt(refPointWrtBS.perp2());
  const auto vertz = refPointWrtBS.z();
  const auto chi2 = track.normalizedChi2();
  const auto chi2prob = TMath::Prob(track.chi2(), (int)track.ndof());
  const bool fillSeedingLayerSets = !seedingLayerSetNames.empty();
  const unsigned int seedingLayerSetBin = fillSeedingLayerSets ? getSeedingLayerSetBin(track, ttopo) : 0;
  const auto simpvz = simPVPosition ? simPVPosition->z() : 0.0;

  const bool paramsValid = !trackFromSeedFitFailed(track);
  const bool isSelected = (*trackSelectorVsPhi)(track, bsPosition);
  const bool isDuplicate = (numAssocRecoTracks > 1);
  const bool isPileup = isMatched && (!isSigMatched);
  isChargeMatched = isChargeMatched || (!doSeedPlots_);

  if (paramsValid) {
    histograms.hs_eta[count].fillRecoHistos(
        isMatched, (*trackSelectorVsEta)(track, bsPosition), isDuplicate, isPileup, isChargeMatched, eta);
    histograms.hs_phi[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, phi);
    histograms.hs_pT[count].fillRecoHistos(
        isMatched, (*trackSelectorVsPt)(track, bsPosition), isDuplicate, isPileup, isChargeMatched, pt);
    histograms.hs_pTvseta[count].fillRecoHistos(
        isMatched, (*trackSelectorVsPt)(track, bsPosition), isDuplicate, isPileup, isChargeMatched, eta, pt);
    histograms.hs_chi2[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, chi2);
    histograms.hs_chi2prob[count].fillRecoHistos(
        isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, chi2prob);
    histograms.hs_hit[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, nhits);
    histograms.hs_layer[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, nlayers);
    histograms.hs_pixellayer[count].fillRecoHistos(
        isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, nPixelLayers);
    histograms.hs_3Dlayer[count].fillRecoHistos(
        isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, n3DLayers);
    histograms.hs_pu[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, numVertices);
    histograms.hs_dr[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dR);
    histograms.hs_drj[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dRJet);
    histograms.hs_vertpos[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, vertxy);
    histograms.hs_dxy[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dxy);
    histograms.hs_zpos[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, vertz);
    histograms.hs_dz[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dz);
    if (pvPosition) {
      histograms.hs_dxypv[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dxypv);
      histograms.hs_dxypvzoomed[count].fillRecoHistos(
          isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dxypv);
      histograms.hs_dzpv[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dzpv);
      histograms.hs_dzpvzoomed[count].fillRecoHistos(
          isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, dzpv);
      if (doDzPVcutPlots_)
        histograms.hs_dzpvcut[count].fillRecoHistos(
            isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, std::abs(dzpv));
    }
    if (simPVPosition)
      histograms.hs_simpvz[count].fillRecoHistos(isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, simpvz);

    if (dRJet <= 99999)  //dRJet can be set to numeric_limits max^2, this is a protection
      histograms.chi2_vs_drj[count]->Fill(dRJet, chi2);
    if (fillSeedingLayerSets)
      histograms.hs_seedingLayerSet[count].fillRecoHistos(
          isMatched, isSelected, isDuplicate, isPileup, isChargeMatched, seedingLayerSetBin);

    if (isMatched) {
      histograms.assoc_chi2_vs_eta[count]->Fill(eta, chi2);
      histograms.assoc_chi2prob_vs_eta[count]->Fill(eta, chi2prob);
      histograms.assoc_chi2_vs_pt[count]->Fill(pt, chi2);
      histograms.assoc_chi2prob_vs_pt[count]->Fill(pt, chi2prob);
      if (dRJet <= 99999) {  //dRJet can be set to numeric_limits max^2, this is a protection
        histograms.assoc_chi2_vs_drj[count]->Fill(dRJet, chi2);
        histograms.assoc_chi2prob_vs_drj[count]->Fill(dRJet, chi2prob);
      }
    }
  }

  if (!mvas.empty()) {
    assert(histograms.h_reco_mva.size() > static_cast<size_t>(count));
    assert(histograms.h_reco_mvacut.size() > static_cast<size_t>(count));
    assert(histograms.h_reco_mva_hp.size() > static_cast<size_t>(count));
    assert(histograms.h_reco_mvacut_hp.size() > static_cast<size_t>(count));

    fillMVAHistos(histograms.h_reco_mva[count],
                  histograms.h_reco_mvacut[count],
                  histograms.h_reco_mva_hp[count],
                  histograms.h_reco_mvacut_hp[count],
                  mvas,
                  selectsLoose,
                  selectsHP);

    if (isMatched) {
      fillMVAHistos(histograms.h_assocRecoToSim_mva[count],
                    histograms.h_assocRecoToSim_mvacut[count],
                    histograms.h_assocRecoToSim_mva_hp[count],
                    histograms.h_assocRecoToSim_mvacut_hp[count],
                    mvas,
                    selectsLoose,
                    selectsHP);
      assert(histograms.h_assocRecoToSim_mva_vs_pt.size() > static_cast<size_t>(count));
      assert(histograms.h_assocRecoToSim_mva_vs_pt_hp.size() > static_cast<size_t>(count));
      fillMVAHistos(pt,
                    histograms.h_assocRecoToSim_mva_vs_pt[count],
                    histograms.h_assocRecoToSim_mva_vs_pt_hp[count],
                    mvas,
                    selectsLoose,
                    selectsHP);
      assert(histograms.h_assocRecoToSim_mva_vs_eta.size() > static_cast<size_t>(count));
      assert(histograms.h_assocRecoToSim_mva_vs_eta_hp.size() > static_cast<size_t>(count));
      fillMVAHistos(eta,
                    histograms.h_assocRecoToSim_mva_vs_eta[count],
                    histograms.h_assocRecoToSim_mva_vs_eta_hp[count],
                    mvas,
                    selectsLoose,
                    selectsHP);
    } else {  // !isMatched
      assert(histograms.h_fake_mva_vs_pt.size() > static_cast<size_t>(count));
      assert(histograms.h_fake_mva_vs_pt_hp.size() > static_cast<size_t>(count));
      assert(histograms.h_fake_mva_vs_eta.size() > static_cast<size_t>(count));
      assert(histograms.h_fake_mva_vs_eta_hp.size() > static_cast<size_t>(count));
      fillMVAHistos(
          pt, histograms.h_fake_mva_vs_pt[count], histograms.h_fake_mva_vs_pt_hp[count], mvas, selectsLoose, selectsHP);
      fillMVAHistos(eta,
                    histograms.h_fake_mva_vs_eta[count],
                    histograms.h_fake_mva_vs_eta_hp[count],
                    mvas,
                    selectsLoose,
                    selectsHP);
    }
  }

  if (isMatched) {
    if (histograms.nrecHit_vs_nsimHit_rec2sim[count])
      histograms.nrecHit_vs_nsimHit_rec2sim[count]->Fill(track.numberOfValidHits(), nSimHits);
    histograms.h_assocSimToReco_Fraction[count]->Fill(sharedFraction);
    histograms.h_assocSimToReco_SharedHit[count]->Fill(sharedHits);
  }
}

void MTVHistoProducerAlgoForTracker::fill_simAssociated_recoTrack_histos(const Histograms& histograms,
                                                                         int count,
                                                                         const reco::Track& track) const {
  //nchi2 and hits global distributions
  histograms.h_hits[count]->Fill(track.numberOfValidHits());
  histograms.h_losthits[count]->Fill(track.numberOfLostHits());
  histograms.h_nmisslayers_inner[count]->Fill(
      track.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
  histograms.h_nmisslayers_outer[count]->Fill(
      track.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
  if (trackFromSeedFitFailed(track))
    return;

  histograms.h_nchi2[count]->Fill(track.normalizedChi2());
  histograms.h_nchi2_prob[count]->Fill(TMath::Prob(track.chi2(), (int)track.ndof()));
  histograms.chi2_vs_nhits[count]->Fill(track.numberOfValidHits(), track.normalizedChi2());
  histograms.h_charge[count]->Fill(track.charge());

  //chi2 and #hit vs eta: fill 2D histos
  const auto eta = getEta(track.eta());
  histograms.chi2_vs_eta[count]->Fill(eta, track.normalizedChi2());
  histograms.nhits_vs_eta[count]->Fill(eta, track.numberOfValidHits());
  const auto pt = getPt(sqrt(track.momentum().perp2()));
  histograms.chi2_vs_pt[count]->Fill(pt, track.normalizedChi2());
  const auto pxbHits = track.hitPattern().numberOfValidPixelBarrelHits();
  const auto pxfHits = track.hitPattern().numberOfValidPixelEndcapHits();
  const auto tibHits = track.hitPattern().numberOfValidStripTIBHits();
  const auto tidHits = track.hitPattern().numberOfValidStripTIDHits();
  const auto tobHits = track.hitPattern().numberOfValidStripTOBHits();
  const auto tecHits = track.hitPattern().numberOfValidStripTECHits();
  histograms.nPXBhits_vs_eta[count]->Fill(eta, pxbHits);
  histograms.nPXFhits_vs_eta[count]->Fill(eta, pxfHits);
  histograms.nPXLhits_vs_eta[count]->Fill(eta, pxbHits + pxfHits);
  histograms.nTIBhits_vs_eta[count]->Fill(eta, tibHits);
  histograms.nTIDhits_vs_eta[count]->Fill(eta, tidHits);
  histograms.nTOBhits_vs_eta[count]->Fill(eta, tobHits);
  histograms.nTEChits_vs_eta[count]->Fill(eta, tecHits);
  histograms.nSTRIPhits_vs_eta[count]->Fill(eta, tibHits + tidHits + tobHits + tecHits);
  histograms.nLayersWithMeas_vs_eta[count]->Fill(eta, track.hitPattern().trackerLayersWithMeasurement());
  histograms.nPXLlayersWithMeas_vs_eta[count]->Fill(eta, track.hitPattern().pixelLayersWithMeasurement());
  if (doMTDPlots_) {
    //  const auto mtdHits = track.hitPattern().numberOfValidTimingHits();
    const auto btlHits = track.hitPattern().numberOfValidTimingBTLHits();
    const auto etlHits = track.hitPattern().numberOfValidTimingETLHits();
    histograms.nMTDhits_vs_eta[count]->Fill(eta, btlHits + etlHits);
    histograms.nBTLhits_vs_eta[count]->Fill(eta, btlHits);
    histograms.nETLhits_vs_eta[count]->Fill(eta, etlHits);
  }
  int LayersAll = track.hitPattern().stripLayersWithMeasurement();
  int Layers2D = track.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
  int Layers1D = LayersAll - Layers2D;
  histograms.nSTRIPlayersWithMeas_vs_eta[count]->Fill(eta, LayersAll);
  histograms.nSTRIPlayersWith1dMeas_vs_eta[count]->Fill(eta, Layers1D);
  histograms.nSTRIPlayersWith2dMeas_vs_eta[count]->Fill(eta, Layers2D);

  histograms.nlosthits_vs_eta[count]->Fill(eta, track.numberOfLostHits());
}

void MTVHistoProducerAlgoForTracker::fill_trackBased_histos(const Histograms& histograms,
                                                            int count,
                                                            int assTracks,
                                                            int numRecoTracks,
                                                            int numRecoTracksSelected,
                                                            int numSimTracksSelected) const {
  histograms.h_tracks[count]->Fill(assTracks);
  histograms.h_fakes[count]->Fill(numRecoTracks - assTracks);
  if (histograms.nrec_vs_nsim[count])
    histograms.nrec_vs_nsim[count]->Fill(numSimTracksSelected, numRecoTracksSelected);
}

void MTVHistoProducerAlgoForTracker::fill_ResoAndPull_recoTrack_histos(const Histograms& histograms,
                                                                       int count,
                                                                       const TrackingParticle::Vector& momentumTP,
                                                                       const TrackingParticle::Point& vertexTP,
                                                                       int chargeTP,
                                                                       const reco::Track& track,
                                                                       const math::XYZPoint& bsPosition) const {
  if (trackFromSeedFitFailed(track))
    return;

  // evaluation of TP parameters
  double qoverpSim = chargeTP / sqrt(momentumTP.x() * momentumTP.x() + momentumTP.y() * momentumTP.y() +
                                     momentumTP.z() * momentumTP.z());
  double lambdaSim = M_PI / 2 - momentumTP.theta();
  double phiSim = momentumTP.phi();
  double dxySim = TrackingParticleIP::dxy(vertexTP, momentumTP, bsPosition);
  double dzSim = TrackingParticleIP::dz(vertexTP, momentumTP, bsPosition);

  //  reco::Track::ParameterVector rParameters = track.parameters(); // UNUSED

  double qoverpRec(0);
  double qoverpErrorRec(0);
  double ptRec(0);
  double ptErrorRec(0);
  double lambdaRec(0);
  double lambdaErrorRec(0);
  double phiRec(0);
  double phiErrorRec(0);

  /* TO BE FIXED LATER  -----------
  //loop to decide whether to take gsfTrack (utilisation of mode-function) or common track
  const GsfTrack* gsfTrack(0);
  if(useGsf){
    gsfTrack = dynamic_cast<const GsfTrack*>(&(*track));
    if (gsfTrack==0) edm::LogInfo("TrackValidator") << "Trying to access mode for a non-GsfTrack";
  }

  if (gsfTrack) {
    // get values from mode
    getRecoMomentum(*gsfTrack, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec,
		    lambdaRec,lambdaErrorRec, phiRec, phiErrorRec);
  }

  else {
    // get values from track (without mode)
    getRecoMomentum(*track, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec,
		    lambdaRec,lambdaErrorRec, phiRec, phiErrorRec);
  }
  */
  getRecoMomentum(track, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec, lambdaRec, lambdaErrorRec, phiRec, phiErrorRec);
  // -------------

  double ptError = ptErrorRec;
  double ptres = ptRec - sqrt(momentumTP.perp2());
  double etares = track.eta() - momentumTP.Eta();

  double dxyRec = track.dxy(bsPosition);
  double dzRec = track.dz(bsPosition);

  const auto phiRes = phiRec - phiSim;
  const auto dxyRes = dxyRec - dxySim;
  const auto dzRes = dzRec - dzSim;
  const auto cotThetaRes = 1 / tan(M_PI * 0.5 - lambdaRec) - 1 / tan(M_PI * 0.5 - lambdaSim);
  const double thetaRes = (lambdaRec - lambdaSim);

  // eta residue; pt, k, theta, phi, dxy, dz pulls
  double qoverpPull = (qoverpRec - qoverpSim) / qoverpErrorRec;
  double thetaPull = (lambdaRec - lambdaSim) / lambdaErrorRec;
  double phiPull = phiRes / phiErrorRec;
  double dxyPull = dxyRes / track.dxyError();
  double dzPull = dzRes / track.dzError();

#ifdef EDM_ML_DEBUG
  double contrib_Qoverp = ((qoverpRec - qoverpSim) / qoverpErrorRec) * ((qoverpRec - qoverpSim) / qoverpErrorRec) / 5;
  double contrib_dxy = ((dxyRec - dxySim) / track.dxyError()) * ((dxyRec - dxySim) / track.dxyError()) / 5;
  double contrib_dz = ((dzRec - dzSim) / track.dzError()) * ((dzRec - dzSim) / track.dzError()) / 5;
  double contrib_theta = ((lambdaRec - lambdaSim) / lambdaErrorRec) * ((lambdaRec - lambdaSim) / lambdaErrorRec) / 5;
  double contrib_phi = ((phiRec - phiSim) / phiErrorRec) * ((phiRec - phiSim) / phiErrorRec) / 5;

  LogTrace("TrackValidatorTEST")
      //<< "assocChi2=" << tp.begin()->second << "\n"
      << ""
      << "\n"
      << "ptREC=" << ptRec << "\n"
      << "etaREC=" << track.eta() << "\n"
      << "qoverpREC=" << qoverpRec << "\n"
      << "dxyREC=" << dxyRec << "\n"
      << "dzREC=" << dzRec << "\n"
      << "thetaREC=" << track.theta() << "\n"
      << "phiREC=" << phiRec << "\n"
      << ""
      << "\n"
      << "qoverpError()=" << qoverpErrorRec << "\n"
      << "dxyError()=" << track.dxyError() << "\n"
      << "dzError()=" << track.dzError() << "\n"
      << "thetaError()=" << lambdaErrorRec << "\n"
      << "phiError()=" << phiErrorRec << "\n"
      << ""
      << "\n"
      << "ptSIM=" << sqrt(momentumTP.perp2()) << "\n"
      << "etaSIM=" << momentumTP.Eta() << "\n"
      << "qoverpSIM=" << qoverpSim << "\n"
      << "dxySIM=" << dxySim << "\n"
      << "dzSIM=" << dzSim << "\n"
      << "thetaSIM=" << M_PI / 2 - lambdaSim << "\n"
      << "phiSIM=" << phiSim << "\n"
      << ""
      << "\n"
      << "contrib_Qoverp=" << contrib_Qoverp << "\n"
      << "contrib_dxy=" << contrib_dxy << "\n"
      << "contrib_dz=" << contrib_dz << "\n"
      << "contrib_theta=" << contrib_theta << "\n"
      << "contrib_phi=" << contrib_phi << "\n"
      << ""
      << "\n"
      << "chi2PULL=" << contrib_Qoverp + contrib_dxy + contrib_dz + contrib_theta + contrib_phi << "\n";
#endif

  histograms.h_pullQoverp[count]->Fill(qoverpPull);
  histograms.h_pullTheta[count]->Fill(thetaPull);
  histograms.h_pullPhi[count]->Fill(phiPull);
  histograms.h_pullDxy[count]->Fill(dxyPull);
  histograms.h_pullDz[count]->Fill(dzPull);

  const auto etaSim = getEta(momentumTP.eta());
  const auto ptSim = getPt(sqrt(momentumTP.perp2()));

  histograms.h_pt[count]->Fill(ptres / ptError);
  histograms.h_eta[count]->Fill(etares);
  //histograms.etares_vs_eta[count]->Fill(getEta(track.eta()),etares);
  histograms.etares_vs_eta[count]->Fill(etaSim, etares);

  //resolution of track params: fill 2D histos of residuals + pulls
  histograms.hr_dxy[count].fill(etaSim, ptSim, phiSim, dxyRes, dxyPull);
  histograms.hr_pt[count].fill(etaSim, ptSim, phiSim, ptres / ptRec, ptres / ptError);
  histograms.hr_dz[count].fill(etaSim, ptSim, phiSim, dzRes, dzPull);
  histograms.hr_phi[count].fill(etaSim, ptSim, phiSim, phiRes, phiPull);
  histograms.hr_cotTheta[count].fill(etaSim, ptSim, phiSim, cotThetaRes, -1);
  histograms.hr_cotTheta[count].fill(etaSim, ptSim, phiSim, thetaRes, thetaPull);

  //plots vs phi
  histograms.nhits_vs_phi[count]->Fill(phiRec, track.numberOfValidHits());
  histograms.chi2_vs_phi[count]->Fill(phiRec, track.normalizedChi2());
  histograms.ptmean_vs_eta_phi[count]->Fill(phiRec, getEta(track.eta()), ptRec);
  histograms.phimean_vs_eta_phi[count]->Fill(phiRec, getEta(track.eta()), phiRec);
}

void MTVHistoProducerAlgoForTracker::getRecoMomentum(const reco::Track& track,
                                                     double& pt,
                                                     double& ptError,
                                                     double& qoverp,
                                                     double& qoverpError,
                                                     double& lambda,
                                                     double& lambdaError,
                                                     double& phi,
                                                     double& phiError) const {
  pt = track.pt();
  ptError = track.ptError();
  qoverp = track.qoverp();
  qoverpError = track.qoverpError();
  lambda = track.lambda();
  lambdaError = track.lambdaError();
  phi = track.phi();
  phiError = track.phiError();
  //   cout <<"test1" << endl;
}

void MTVHistoProducerAlgoForTracker::getRecoMomentum(const reco::GsfTrack& gsfTrack,
                                                     double& pt,
                                                     double& ptError,
                                                     double& qoverp,
                                                     double& qoverpError,
                                                     double& lambda,
                                                     double& lambdaError,
                                                     double& phi,
                                                     double& phiError) const {
  pt = gsfTrack.ptMode();
  ptError = gsfTrack.ptModeError();
  qoverp = gsfTrack.qoverpMode();
  qoverpError = gsfTrack.qoverpModeError();
  lambda = gsfTrack.lambdaMode();
  lambdaError = gsfTrack.lambdaModeError();
  phi = gsfTrack.phiMode();
  phiError = gsfTrack.phiModeError();
  //   cout <<"test2" << endl;
}

double MTVHistoProducerAlgoForTracker::getEta(double eta) const {
  if (useFabsEta)
    return fabs(eta);
  else
    return eta;
}

double MTVHistoProducerAlgoForTracker::getPt(double pt) const {
  if (useInvPt && pt != 0)
    return 1 / pt;
  else
    return pt;
}

unsigned int MTVHistoProducerAlgoForTracker::getSeedingLayerSetBin(const reco::Track& track,
                                                                   const TrackerTopology& ttopo) const {
  if (track.seedRef().isNull() || !track.seedRef().isAvailable())
    return seedingLayerSetNames.size() - 1;

  const TrajectorySeed& seed = *(track.seedRef());
  SeedingLayerSetId searchId;
  const int nhits = seed.nHits();
  if (nhits > static_cast<int>(std::tuple_size<SeedingLayerSetId>::value)) {
    LogDebug("TrackValidator") << "Got seed with " << nhits << " hits, but I have a hard-coded maximum of "
                               << std::tuple_size<SeedingLayerSetId>::value
                               << ", classifying the seed as 'unknown'. Please increase the maximum in "
                                  "MTVHistoProducerAlgoForTracker.h if needed.";
    return seedingLayerSetNames.size() - 1;
  }
  int i = 0;
  for (auto const& recHit : seed.recHits()) {
    DetId detId = recHit.geographicalId();

    if (detId.det() != DetId::Tracker) {
      throw cms::Exception("LogicError") << "Encountered seed hit detId " << detId.rawId() << " not from Tracker, but "
                                         << detId.det();
    }

    GeomDetEnumerators::SubDetector subdet;
    bool subdetStrip = false;
    switch (detId.subdetId()) {
      case PixelSubdetector::PixelBarrel:
        subdet = GeomDetEnumerators::PixelBarrel;
        break;
      case PixelSubdetector::PixelEndcap:
        subdet = GeomDetEnumerators::PixelEndcap;
        break;
      case StripSubdetector::TIB:
        subdet = GeomDetEnumerators::TIB;
        subdetStrip = true;
        break;
      case StripSubdetector::TID:
        subdet = GeomDetEnumerators::TID;
        subdetStrip = true;
        break;
      case StripSubdetector::TOB:
        subdet = GeomDetEnumerators::TOB;
        subdetStrip = true;
        break;
      case StripSubdetector::TEC:
        subdet = GeomDetEnumerators::TEC;
        subdetStrip = true;
        break;
      default:
        throw cms::Exception("LogicError") << "Unknown subdetId " << detId.subdetId();
    };

    TrackerDetSide side = static_cast<TrackerDetSide>(ttopo.side(detId));

    // Even with the recent addition of
    // SeedingLayerSetsBuilder::fillDescription() this assumption is a
    // bit ugly.
    const bool isStripMono = subdetStrip && trackerHitRTTI::isSingle(recHit);
    searchId[i] =
        SeedingLayerId(SeedingLayerSetsBuilder::SeedingLayerId(subdet, side, ttopo.layer(detId)), isStripMono);
    ++i;
  }
  auto found = seedingLayerSetToBin.find(searchId);
  if (found == seedingLayerSetToBin.end()) {
    return seedingLayerSetNames.size() - 1;
  }
  return found->second;
}

void MTVHistoProducerAlgoForTracker::fill_recoAssociated_simTrack_histos(const Histograms& histograms,
                                                                         int count,
                                                                         const reco::GenParticle& tp,
                                                                         const TrackingParticle::Vector& momentumTP,
                                                                         const TrackingParticle::Point& vertexTP,
                                                                         double dxySim,
                                                                         double dzSim,
                                                                         int nSimHits,
                                                                         const reco::Track* track,
                                                                         int numVertices) const {
  bool isMatched = track;
  bool isReconstructable = (*GpSelectorForTechnicalEfficiency)(tp);

  if ((*GpSelectorForEfficiencyVsEta)(tp))
    //effic vs eta
    histograms.hs_eta[count].fillSimTrackHistos(isMatched, isReconstructable, getEta(momentumTP.eta()));

  if ((*GpSelectorForEfficiencyVsPhi)(tp)) {
    histograms.hs_phi[count].fillSimTrackHistos(isMatched, isReconstructable, momentumTP.phi());
    histograms.hs_hit[count].fillSimTrackHistos(isMatched, isReconstructable, (int)nSimHits);
    histograms.hs_pu[count].fillSimTrackHistos(isMatched, isReconstructable, numVertices);
    if (isMatched && histograms.nrecHit_vs_nsimHit_sim2rec[count])
      histograms.nrecHit_vs_nsimHit_sim2rec[count]->Fill(track->numberOfValidHits(), nSimHits);
  }

  if ((*GpSelectorForEfficiencyVsPt)(tp)) {
    histograms.hs_pT[count].fillSimTrackHistos(isMatched, isReconstructable, getPt(sqrt(momentumTP.perp2())));
    histograms.hs_pTvseta[count].fillSimTrackHistos(
        isMatched, isReconstructable, getEta(momentumTP.eta()), getPt(sqrt(momentumTP.perp2())));
  }

  if ((*GpSelectorForEfficiencyVsVTXR)(tp)) {
    histograms.hs_vertpos[count].fillSimTrackHistos(isMatched, isReconstructable, sqrt(vertexTP.perp2()));
    histograms.hs_dxy[count].fillSimTrackHistos(isMatched, isReconstructable, dxySim);
  }

  if ((*GpSelectorForEfficiencyVsVTXZ)(tp)) {
    histograms.hs_zpos[count].fillSimTrackHistos(isMatched, isReconstructable, vertexTP.z());
    histograms.hs_dz[count].fillSimTrackHistos(isMatched, isReconstructable, dzSim);
  }
}

void MTVHistoProducerAlgoForTracker::fill_seed_histos(const Histograms& histograms,
                                                      int count,
                                                      int seedsFitFailed,
                                                      int seedsTotal) const {
  histograms.h_seedsFitFailed[count]->Fill(seedsFitFailed);
  histograms.h_seedsFitFailedFraction[count]->Fill(static_cast<double>(seedsFitFailed) / seedsTotal);
}
