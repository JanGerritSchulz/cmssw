#ifndef Validation_RecoTrack_MTVResolutionBundle_h
#define Validation_RecoTrack_MTVResolutionBundle_h

// Package:    Validation/RecoTrack
// Class:      MTVResolutionBundle
//
/**\class MTVResolutionBundle Validation/RecoTrack/MTVResolutionBundle.h
 Description: Bundle of monitoring elements for the MultiTrackValidator to have all histograms 
              used for resolutions of the same variable at one unique place 
              (e.g. the various resolution histograms of pT: ptres_vs_phi, ptres_vs_eta, ptres_vs_pt)
*/
//
// Original Author:  Jan Schulz (2026)

#include "Validation/RecoTrack/interface/MTVBookers.h"

// -----------------------------------------------------------------------------------------------------------------

class MTVResolutionBundle {
public:
  MTVResolutionBundle() {}
  virtual ~MTVResolutionBundle() = default;

  void fill(const double eta, const double pt, const double phi, const double residual, const double pull) {
    h_res_vs_eta->Fill(eta, residual);
    h_res_vs_pt->Fill(pt, residual);
    h_res_vs_phi->Fill(phi, residual);
    h_pull_vs_eta->Fill(eta, pull);
    h_pull_vs_pt->Fill(pt, pull);
    h_pull_vs_phi->Fill(phi, pull);
  }

  // booking all histograms for resolution plots
  void bookResolutions(DQMStore::IBooker& ibooker,
                       const int etaNBins,
                       const double etaValMin,
                       const double etaValMax,
                       const int phiNBins,
                       const double phiValMin,
                       const double phiValMax,
                       const int ptNBins,
                       const double ptValMin,
                       const double ptValMax,
                       const bool ptLog,
                       const std::string& name,
                       const int resNBins,
                       const double resValMin,
                       const double resValMax) {
    h_res_vs_eta = make2DIfLogX(ibooker,
                                false,
                                (name + "res_vs_eta").c_str(),
                                (name + " residuals vs #eta;Pseudorapidity #eta;Residuals").c_str(),
                                etaNBins,
                                etaValMin,
                                etaValMax,
                                resNBins,
                                resValMin,
                                resValMax);
    h_res_vs_pt = make2DIfLogX(ibooker,
                               ptLog,
                               (name + "res_vs_pt").c_str(),
                               (name + " residuals vs p_{T};p_{T};Residuals").c_str(),
                               ptNBins,
                               ptValMin,
                               ptValMax,
                               resNBins,
                               resValMin,
                               resValMax);
    h_res_vs_phi = make2DIfLogX(ibooker,
                                false,
                                (name + "res_vs_phi").c_str(),
                                (name + " residuals vs #phi;#phi;Residuals").c_str(),
                                phiNBins,
                                phiValMin,
                                phiValMax,
                                resNBins,
                                resValMin,
                                resValMax);
    h_pull_vs_eta = make2DIfLogX(ibooker,
                                 false,
                                 (name + "pull_vs_eta").c_str(),
                                 (name + " pulls vs #eta;Pseudorapidity #eta;Pulls").c_str(),
                                 etaNBins,
                                 etaValMin,
                                 etaValMax,
                                 100,
                                 -10,
                                 10);
    h_pull_vs_pt = make2DIfLogX(ibooker,
                                ptLog,
                                (name + "pull_vs_pt").c_str(),
                                (name + " pulls vs p_{T};p_{T};Pulls").c_str(),
                                ptNBins,
                                ptValMin,
                                ptValMax,
                                100,
                                -10,
                                10);
    h_pull_vs_phi = make2DIfLogX(ibooker,
                                 false,
                                 (name + "pull_vs_phi").c_str(),
                                 (name + " pulls vs #phi;#phi;Pulls").c_str(),
                                 phiNBins,
                                 phiValMin,
                                 phiValMax,
                                 100,
                                 -10,
                                 10);
  }

  // allow modification of the histograms,
  // e.g. for setting the labels on the x-ticks (for collection summary plots)
  template <typename ModificationFunc, typename... Args>
  void modifyHistograms(ModificationFunc modify, Args&&... args) {
    for (auto h : {h_res_vs_eta, h_res_vs_pt, h_res_vs_phi, h_pull_vs_eta, h_pull_vs_pt, h_pull_vs_phi}) {
      if (h)
        modify(h, std::forward<Args>(args)...);
    }
  }

private:
  // histograms
  dqm::reco::MonitorElement* h_res_vs_eta = nullptr;
  dqm::reco::MonitorElement* h_res_vs_pt = nullptr;
  dqm::reco::MonitorElement* h_res_vs_phi = nullptr;
  dqm::reco::MonitorElement* h_pull_vs_eta = nullptr;
  dqm::reco::MonitorElement* h_pull_vs_pt = nullptr;
  dqm::reco::MonitorElement* h_pull_vs_phi = nullptr;
};

#endif
