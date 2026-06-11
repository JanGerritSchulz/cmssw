#ifndef Validation_RecoVertex_SVResolutionBundle_h
#define Validation_RecoVertex_SVResolutionBundle_h

// Package:    Validation/RecoVertex
// Class:      SVResolutionBundle
//
/**\class SVResolutionBundle Validation/RecoVertex/interface/SVResolutionBundle.h

 Description: Bundle of monitoring elements for the SecondaryVertexAnalyzer,
              grouping residual and pull histograms for a single vertex
              quantity at one place (e.g. x-position residuals and pulls vs
              decay length, eta, and nTracks).

              The SV case differs from the track case (MTVResolutionBundle) in
              the choice of kinematic axes: rather than eta/pT/phi as for tracks,
              the natural axes for SV resolution studies are decay length (the
              most physically relevant SV variable), transverse radius, and
              track multiplicity. These are configurable at booking time.

              A dedicated decay-length resolution bundle (SVDecayLengthBundle)
              is provided separately since decay length and its significance are
              the most important SV-specific quantities and deserve dedicated
              significance and pull treatment.

 Original Author: Jan Schulz
*/

#include "DQMServices/Core/interface/DQMBookingHelpers.h"

// =============================================================================
// SVResolutionBundle
//
// Generic resolution bundle: residuals and pulls of any vertex quantity
// vs decay length, transverse radius, and nTracks.
// =============================================================================

class SVResolutionBundle {
public:
  using IBooker = dqm::reco::DQMStore::IBooker;

  SVResolutionBundle() = default;
  virtual ~SVResolutionBundle() = default;

  // Fill all six residual/pull histograms for one matched reco-sim pair.
  void fill(const double decayLength, const double r, const double nTracks, const double residual, const double pull) {
    if (h_res_vs_decayLength)
      h_res_vs_decayLength->Fill(decayLength, residual);
    if (h_res_vs_r)
      h_res_vs_r->Fill(r, residual);
    if (h_res_vs_nTracks)
      h_res_vs_nTracks->Fill(nTracks, residual);
    if (h_pull_vs_decayLength)
      h_pull_vs_decayLength->Fill(decayLength, pull);
    if (h_pull_vs_r)
      h_pull_vs_r->Fill(r, pull);
    if (h_pull_vs_nTracks)
      h_pull_vs_nTracks->Fill(nTracks, pull);
  }

  // Book all histograms for this resolution bundle.
  void bookResolutions(IBooker &ibooker,
                       // decay length axis
                       const int lxyNBins,
                       const double lxyMin,
                       const double lxyMax,
                       // transverse radius axis
                       const int rNBins,
                       const double rMin,
                       const double rMax,
                       // track multiplicity axis
                       const int nTracksNBins,
                       const double nTracksMin,
                       const double nTracksMax,
                       // residual axis
                       const std::string &name,
                       const int resNBins,
                       const double resMin,
                       const double resMax) {
    h_res_vs_decayLength =
        dqm::booking::book2DIfLogX(ibooker,
                                   false,
                                   (name + "res_vs_decayLength").c_str(),
                                   (name + " residuals vs decay length;Decay length L_{3D} [cm];Residuals").c_str(),
                                   lxyNBins,
                                   lxyMin,
                                   lxyMax,
                                   resNBins,
                                   resMin,
                                   resMax);
    h_res_vs_r = dqm::booking::book2DIfLogX(ibooker,
                                            false,
                                            (name + "res_vs_r").c_str(),
                                            (name + " residuals vs transverse radius;r_{T} [cm];Residuals").c_str(),
                                            rNBins,
                                            rMin,
                                            rMax,
                                            resNBins,
                                            resMin,
                                            resMax);
    h_res_vs_nTracks =
        dqm::booking::book2DIfLogX(ibooker,
                                   false,
                                   (name + "res_vs_nTracks").c_str(),
                                   (name + " residuals vs track multiplicity;N tracks;Residuals").c_str(),
                                   nTracksNBins,
                                   nTracksMin,
                                   nTracksMax,
                                   resNBins,
                                   resMin,
                                   resMax);
    h_pull_vs_decayLength =
        dqm::booking::book2DIfLogX(ibooker,
                                   false,
                                   (name + "pull_vs_decayLength").c_str(),
                                   (name + " pulls vs decay length;Decay length L_{3D} [cm];Pulls").c_str(),
                                   lxyNBins,
                                   lxyMin,
                                   lxyMax,
                                   100,
                                   -10.,
                                   10.);
    h_pull_vs_r = dqm::booking::book2DIfLogX(ibooker,
                                             false,
                                             (name + "pull_vs_r").c_str(),
                                             (name + " pulls vs transverse radius;r_{T} [cm];Pulls").c_str(),
                                             rNBins,
                                             rMin,
                                             rMax,
                                             100,
                                             -10.,
                                             10.);
    h_pull_vs_nTracks = dqm::booking::book2DIfLogX(ibooker,
                                                   false,
                                                   (name + "pull_vs_nTracks").c_str(),
                                                   (name + " pulls vs track multiplicity;N tracks;Pulls").c_str(),
                                                   nTracksNBins,
                                                   nTracksMin,
                                                   nTracksMax,
                                                   100,
                                                   -10.,
                                                   10.);
  }

  template <typename ModificationFunc, typename... Args>
  void modifyHistograms(ModificationFunc modify, Args &&...args) {
    for (auto h :
         {h_res_vs_decayLength, h_res_vs_r, h_res_vs_nTracks, h_pull_vs_decayLength, h_pull_vs_r, h_pull_vs_nTracks}) {
      if (h)
        modify(h, std::forward<Args>(args)...);
    }
  }

private:
  dqm::reco::MonitorElement *h_res_vs_decayLength = nullptr;
  dqm::reco::MonitorElement *h_res_vs_r = nullptr;
  dqm::reco::MonitorElement *h_res_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_pull_vs_decayLength = nullptr;
  dqm::reco::MonitorElement *h_pull_vs_r = nullptr;
  dqm::reco::MonitorElement *h_pull_vs_nTracks = nullptr;
};

// =============================================================================
// SVDecayLengthBundle
//
// Dedicated bundle for decay length resolution and significance.
// This is separated from SVResolutionBundle because:
//   - Decay length significance (L/sigma_L) is an independent observable
//     that needs its own residual/pull treatment
//   - The relevant axes are different: nTracks and eta of the SV, not
//     decay length itself (which is the quantity being resolved here)
//   - This bundle is the most b-tagging-relevant one and deserves prominence
// =============================================================================

class SVDecayLengthBundle {
public:
  using IBooker = dqm::reco::DQMStore::IBooker;

  SVDecayLengthBundle() = default;
  virtual ~SVDecayLengthBundle() = default;

  // Fill decay length and significance residuals/pulls together.
  void fill(const double nTracks,
            const double eta,
            const double lResidual,
            const double lPull,
            const double lSigResidual,
            const double lSigPull) {
    if (h_lRes_vs_nTracks)
      h_lRes_vs_nTracks->Fill(nTracks, lResidual);
    if (h_lRes_vs_eta)
      h_lRes_vs_eta->Fill(eta, lResidual);
    if (h_lPull_vs_nTracks)
      h_lPull_vs_nTracks->Fill(nTracks, lPull);
    if (h_lPull_vs_eta)
      h_lPull_vs_eta->Fill(eta, lPull);
    if (h_lSigRes_vs_nTracks)
      h_lSigRes_vs_nTracks->Fill(nTracks, lSigResidual);
    if (h_lSigRes_vs_eta)
      h_lSigRes_vs_eta->Fill(eta, lSigResidual);
    if (h_lSigPull_vs_nTracks)
      h_lSigPull_vs_nTracks->Fill(nTracks, lSigPull);
    if (h_lSigPull_vs_eta)
      h_lSigPull_vs_eta->Fill(eta, lSigPull);
  }

  void bookResolutions(IBooker &ibooker,
                       const int nTracksNBins,
                       const double nTracksMin,
                       const double nTracksMax,
                       const int etaNBins,
                       const double etaMin,
                       const double etaMax,
                       const int lResNBins,
                       const double lResMin,
                       const double lResMax,
                       const int lSigResNBins,
                       const double lSigResMin,
                       const double lSigResMax) {
    // Decay length residuals and pulls
    h_lRes_vs_nTracks = dqm::booking::book2DIfLogX(ibooker,
                                                   false,
                                                   "decayLengthRes_vs_nTracks",
                                                   "L_{3D} residuals vs N tracks;N tracks;L_{3D} residuals [cm]",
                                                   nTracksNBins,
                                                   nTracksMin,
                                                   nTracksMax,
                                                   lResNBins,
                                                   lResMin,
                                                   lResMax);
    h_lRes_vs_eta = dqm::booking::book2DIfLogX(ibooker,
                                               false,
                                               "decayLengthRes_vs_eta",
                                               "L_{3D} residuals vs #eta;#eta;L_{3D} residuals [cm]",
                                               etaNBins,
                                               etaMin,
                                               etaMax,
                                               lResNBins,
                                               lResMin,
                                               lResMax);
    h_lPull_vs_nTracks = dqm::booking::book2DIfLogX(ibooker,
                                                    false,
                                                    "decayLengthPull_vs_nTracks",
                                                    "L_{3D} pulls vs N tracks;N tracks;L_{3D} pulls",
                                                    nTracksNBins,
                                                    nTracksMin,
                                                    nTracksMax,
                                                    100,
                                                    -10.,
                                                    10.);
    h_lPull_vs_eta = dqm::booking::book2DIfLogX(ibooker,
                                                false,
                                                "decayLengthPull_vs_eta",
                                                "L_{3D} pulls vs #eta;#eta;L_{3D} pulls",
                                                etaNBins,
                                                etaMin,
                                                etaMax,
                                                100,
                                                -10.,
                                                10.);
    // Decay length significance residuals and pulls
    h_lSigRes_vs_nTracks =
        dqm::booking::book2DIfLogX(ibooker,
                                   false,
                                   "decayLengthSigRes_vs_nTracks",
                                   "L_{3D}/#sigma_{L} residuals vs N tracks;N tracks;L_{3D}/#sigma_{L} residuals",
                                   nTracksNBins,
                                   nTracksMin,
                                   nTracksMax,
                                   lSigResNBins,
                                   lSigResMin,
                                   lSigResMax);
    h_lSigRes_vs_eta =
        dqm::booking::book2DIfLogX(ibooker,
                                   false,
                                   "decayLengthSigRes_vs_eta",
                                   "L_{3D}/#sigma_{L} residuals vs #eta;#eta;L_{3D}/#sigma_{L} residuals",
                                   etaNBins,
                                   etaMin,
                                   etaMax,
                                   lSigResNBins,
                                   lSigResMin,
                                   lSigResMax);
    h_lSigPull_vs_nTracks =
        dqm::booking::book2DIfLogX(ibooker,
                                   false,
                                   "decayLengthSigPull_vs_nTracks",
                                   "L_{3D}/#sigma_{L} pulls vs N tracks;N tracks;L_{3D}/#sigma_{L} pulls",
                                   nTracksNBins,
                                   nTracksMin,
                                   nTracksMax,
                                   100,
                                   -10.,
                                   10.);
    h_lSigPull_vs_eta = dqm::booking::book2DIfLogX(ibooker,
                                                   false,
                                                   "decayLengthSigPull_vs_eta",
                                                   "L_{3D}/#sigma_{L} pulls vs #eta;#eta;L_{3D}/#sigma_{L} pulls",
                                                   etaNBins,
                                                   etaMin,
                                                   etaMax,
                                                   100,
                                                   -10.,
                                                   10.);
  }

  template <typename ModificationFunc, typename... Args>
  void modifyHistograms(ModificationFunc modify, Args &&...args) {
    for (auto h : {h_lRes_vs_nTracks,
                   h_lRes_vs_eta,
                   h_lPull_vs_nTracks,
                   h_lPull_vs_eta,
                   h_lSigRes_vs_nTracks,
                   h_lSigRes_vs_eta,
                   h_lSigPull_vs_nTracks,
                   h_lSigPull_vs_eta}) {
      if (h)
        modify(h, std::forward<Args>(args)...);
    }
  }

private:
  // Decay length residuals and pulls
  dqm::reco::MonitorElement *h_lRes_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_lRes_vs_eta = nullptr;
  dqm::reco::MonitorElement *h_pull_vs_nTracks = nullptr;  // kept for completeness
  dqm::reco::MonitorElement *h_lPull_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_lPull_vs_eta = nullptr;

  // Decay length significance residuals and pulls
  dqm::reco::MonitorElement *h_lSigRes_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_lSigRes_vs_eta = nullptr;
  dqm::reco::MonitorElement *h_lSigPull_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_lSigPull_vs_eta = nullptr;
};

#endif  // Validation_RecoVertex_SVResolutionBundle_h
