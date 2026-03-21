#ifndef Validation_RecoTrack_MTVHistoProducerAlgoForTracker_h
#define Validation_RecoTrack_MTVHistoProducerAlgoForTracker_h

/* \author B.Mangano, UCSD
 *
 * Concrete class implementing the MTVHistoProducerAlgo interface.
 * To be used within the MTV to fill histograms for Tracker tracks.
 */

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/CandAlgos/interface/GenParticleCustomSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "Validation/RecoTrack/interface/MTVMonitoringElement.h"
#include "Validation/RecoTrack/interface/MTVResolutionBundle.h"

struct MTVHistoProducerAlgoForTrackerHistograms {
  //sim
  using METype = dqm::reco::MonitorElement*;
  METype h_ptSIM, h_etaSIM, h_phiSIM, h_tracksSIM, h_vertposSIM, h_bunchxSIM;

  //1D
  std::vector<METype> h_tracks, h_fakes, h_hits, h_charge, h_algo, h_seedsFitFailed, h_seedsFitFailedFraction;
  mutable std::vector<MTVMonitoringElement> hs_eta, hs_pT, hs_pTvseta, hs_hit, hs_layer, hs_pixellayer, hs_3Dlayer,
      hs_pu, hs_phi, hs_dxy, hs_dz, hs_dxypv, hs_dzpv, hs_dxypvzoomed, hs_dzpvzoomed, hs_vertpos, hs_zpos, hs_dr,
      hs_drj, hs_dzpvcut, hs_dzpvsigcut, hs_simpvz, hs_chi2, hs_chi2prob, hs_seedingLayerSet;
  std::vector<METype> h_pt, h_eta, h_pullTheta, h_pullPhi, h_pullDxy, h_pullDz, h_pullQoverp;
  std::vector<METype> h_assocRecoToSim_itpu_eta, h_assocRecoToSim_itpu_sig_eta, h_assocRecoToSim_eta_sig;
  std::vector<METype> h_assocRecoToSim_itpu_vertcount, h_assocRecoToSim_itpu_sig_vertcount;
  std::vector<METype> h_assocRecoToSim_ootpu_eta, h_assocRecoToSim_ootpu_vertcount;
  std::vector<METype> h_reco_ootpu_eta, h_reco_ootpu_vertcount;
  std::vector<METype> h_con_eta, h_con_vertcount, h_con_zpos;
  std::vector<METype> h_simul2_dzpvcut;

  std::vector<std::vector<METype>> h_reco_mva, h_assocRecoToSim_mva;
  std::vector<std::vector<METype>> h_reco_mvacut, h_assocSimToReco_mvacut, h_assocRecoToSim_mvacut, h_simul2_mvacut;
  std::vector<std::vector<METype>> h_reco_mva_hp, h_assocRecoToSim_mva_hp;
  std::vector<std::vector<METype>> h_reco_mvacut_hp, h_assocSimToReco_mvacut_hp, h_assocRecoToSim_mvacut_hp,
      h_simul2_mvacut_hp;

  std::vector<std::vector<METype>> h_assocRecoToSim_mva_vs_pt, h_fake_mva_vs_pt, h_assocRecoToSim_mva_vs_pt_hp,
      h_fake_mva_vs_pt_hp;
  std::vector<std::vector<METype>> h_assocRecoToSim_mva_vs_eta, h_fake_mva_vs_eta, h_assocRecoToSim_mva_vs_eta_hp,
      h_fake_mva_vs_eta_hp;

  // dE/dx
  // in the future these might become an array
  std::vector<std::vector<METype>> h_dedx_estim;
  std::vector<std::vector<METype>> h_dedx_nom;
  std::vector<std::vector<METype>> h_dedx_sat;

  //2D
  std::vector<METype> nrec_vs_nsim;
  std::vector<METype> nrecHit_vs_nsimHit_sim2rec;
  std::vector<METype> nrecHit_vs_nsimHit_rec2sim;
  std::vector<METype> h_duplicates_oriAlgo_vs_oriAlgo;

  //assoc hits
  std::vector<METype> h_assocSimToReco_Fraction, h_assocSimToReco_SharedHit;

  //#hit vs eta: to be used with doProfileX
  std::vector<METype> nhits_vs_eta, nPXBhits_vs_eta, nPXFhits_vs_eta, nPXLhits_vs_eta, nTIBhits_vs_eta, nTIDhits_vs_eta,
      nTOBhits_vs_eta, nTEChits_vs_eta, nSTRIPhits_vs_eta, nLayersWithMeas_vs_eta, nPXLlayersWithMeas_vs_eta,
      nSTRIPlayersWithMeas_vs_eta, nSTRIPlayersWith1dMeas_vs_eta, nSTRIPlayersWith2dMeas_vs_eta, nMTDhits_vs_eta,
      nBTLhits_vs_eta, nETLhits_vs_eta;

  //---- second set of histograms (originally not used by the SeedGenerator)
  //1D
  std::vector<METype> h_nchi2, h_nchi2_prob, h_losthits, h_nmisslayers_inner, h_nmisslayers_outer;

  //2D
  std::vector<METype> chi2_vs_nhits, etares_vs_eta;
  std::vector<METype> h_ptshifteta;
  std::vector<METype> chi2_vs_phi, nhits_vs_phi;

  //Profile2D
  std::vector<METype> ptmean_vs_eta_phi, phimean_vs_eta_phi;

  //assoc chi2
  std::vector<METype> h_assocSimToReco_chi2, h_assocSimToReco_chi2_prob;

  //chi2 and # lost hits vs eta: to be used with doProfileX
  std::vector<METype> chi2_vs_eta, chi2_vs_pt, chi2_vs_drj, nlosthits_vs_eta;
  std::vector<METype> assoc_chi2_vs_eta, assoc_chi2_vs_pt, assoc_chi2_vs_drj, assoc_chi2prob_vs_eta,
      assoc_chi2prob_vs_pt, assoc_chi2prob_vs_drj;

  //   resolution of track params: to be used with fitslicesytool
  // + pulls of track params vs eta: to be used with fitslicesytool
  mutable std::vector<MTVResolutionBundle> hr_dxy, hr_pt, hr_dz, hr_phi, hr_cotTheta, hr_theta;

  //std::vector<METype> dxyres_vs_eta, ptres_vs_eta, dzres_vs_eta, phires_vs_eta, cotThetares_vs_eta;
  //std::vector<METype> dxyres_vs_pt, ptres_vs_pt, dzres_vs_pt, phires_vs_pt, cotThetares_vs_pt;
  //std::vector<METype> dxypull_vs_eta, ptpull_vs_eta, dzpull_vs_eta, phipull_vs_eta, thetapull_vs_eta;
  //std::vector<METype> dxypull_vs_pt, ptpull_vs_pt, dzpull_vs_pt, phipull_vs_pt, thetapull_vs_pt;
  //std::vector<METype> ptpull_vs_phi, phipull_vs_phi, thetapull_vs_phi;
};

class MTVHistoProducerAlgoForTracker {
public:
  typedef dqm::reco::DQMStore DQMStore;

  MTVHistoProducerAlgoForTracker(const edm::ParameterSet& pset, const bool doSeedPlots);
  ~MTVHistoProducerAlgoForTracker();

  static std::unique_ptr<RecoTrackSelectorBase> makeRecoTrackSelectorFromTPSelectorParameters(
      const edm::ParameterSet& pset);

  bool tpIsReconstructable(const TrackingParticle& tp) { return (*TpSelectorForTechnicalEfficiency)(tp); }

  using Histograms = MTVHistoProducerAlgoForTrackerHistograms;
  void pushbackNewMTVMonitoringElements(Histograms& histograms);
  void bookSimHistos(DQMStore::IBooker& ibook, Histograms& histograms);
  void bookSimAndRecoTrackHistos(DQMStore::IBooker& ibook,
                                 Histograms& histograms,
                                 const bool doSimTrackPlots,
                                 const bool doRecoTrackPlots,
                                 const bool doResolutionPlots);
  void bookSimTrackPVAssociationHistos(DQMStore::IBooker& ibook, Histograms& histograms);
  void bookRecoPVAssociationHistos(DQMStore::IBooker& ibook, Histograms& histograms);
  void bookRecodEdxHistos(DQMStore::IBooker& ibook, Histograms& histograms);
  void bookSeedHistos(DQMStore::IBooker& ibook, Histograms& histograms);
  void bookMVAHistos(DQMStore::IBooker& ibook, Histograms& histograms, size_t nMVAs);

  void fill_generic_simTrack_histos(const Histograms& histograms,
                                    const TrackingParticle::Vector&,
                                    const TrackingParticle::Point& vertex,
                                    int bx) const;
  void fill_simTrackBased_histos(const Histograms& histograms, int numSimTracks) const;

  void fill_recoAssociated_simTrack_histos(const Histograms& histograms,
                                           int count,
                                           const TrackingParticle& tp,
                                           const TrackingParticle::Vector& momentumTP,
                                           const TrackingParticle::Point& vertexTP,
                                           double dxy,
                                           double dz,
                                           double dxyPV,
                                           double dzPV,
                                           int nSimHits,
                                           int nSimLayers,
                                           int nSimPixelLayers,
                                           int nSimStripMonoAndStereoLayers,
                                           const reco::Track* track,
                                           int numVertices,
                                           double dR,
                                           double dR_jet,
                                           const math::XYZPoint* pvPosition,
                                           const TrackingVertex::LorentzVector* simPVPosition,
                                           const math::XYZPoint& bsPosition,
                                           const std::vector<float>& mvas,
                                           unsigned int selectsLoose,
                                           unsigned int selectsHP) const;

  void fill_recoAssociated_simTrack_histos(const Histograms& histograms,
                                           int count,
                                           const reco::GenParticle& tp,
                                           const TrackingParticle::Vector& momentumTP,
                                           const TrackingParticle::Point& vertexTP,
                                           double dxy,
                                           double dz,
                                           int nSimHits,
                                           const reco::Track* track,
                                           int numVertices) const;

  void fill_duplicate_histos(const Histograms& histograms,
                             int count,
                             const reco::Track& track1,
                             const reco::Track& track2) const;

  void fill_generic_recoTrack_histos(const Histograms& histograms,
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
                                     double dR_jet,
                                     const std::vector<float>& mvas,
                                     unsigned int selectsLoose,
                                     unsigned int selectsHP) const;

  void fill_dedx_recoTrack_histos(const Histograms& histograms,
                                  int count,
                                  const edm::RefToBase<reco::Track>& trackref,
                                  const std::vector<const edm::ValueMap<reco::DeDxData>*>& v_dEdx) const;

  void fill_simAssociated_recoTrack_histos(const Histograms& histograms, int count, const reco::Track& track) const;

  void fill_trackBased_histos(const Histograms& histograms,
                              int count,
                              int assTracks,
                              int numRecoTracks,
                              int numRecoTracksSelected,
                              int numSimTracksSelected) const;

  void fill_ResoAndPull_recoTrack_histos(const Histograms& histograms,
                                         int count,
                                         const TrackingParticle::Vector& momentumTP,
                                         const TrackingParticle::Point& vertexTP,
                                         int chargeTP,
                                         const reco::Track& track,
                                         const math::XYZPoint& bsPosition) const;

  void fill_seed_histos(const Histograms& histograms, int count, int seedsFitFailed, int seedsTotal) const;

private:
  /// retrieval of reconstructed momentum components from reco::Track (== mean values for GSF)
  void getRecoMomentum(const reco::Track& track,
                       double& pt,
                       double& ptError,
                       double& qoverp,
                       double& qoverpError,
                       double& lambda,
                       double& lambdaError,
                       double& phi,
                       double& phiError) const;
  /// retrieval of reconstructed momentum components based on the mode of a reco::GsfTrack
  void getRecoMomentum(const reco::GsfTrack& gsfTrack,
                       double& pt,
                       double& ptError,
                       double& qoverp,
                       double& qoverpError,
                       double& lambda,
                       double& lambdaError,
                       double& phi,
                       double& phiError) const;

  double getEta(double eta) const;

  double getPt(double pt) const;

  unsigned int getSeedingLayerSetBin(const reco::Track& track, const TrackerTopology& ttopo) const;

  //private data members
  std::unique_ptr<TrackingParticleSelector> generalTpSelector;
  std::unique_ptr<TrackingParticleSelector> TpSelectorForTechnicalEfficiency;
  std::unique_ptr<TrackingParticleSelector> TpSelectorForEfficiencyVsEta;
  std::unique_ptr<TrackingParticleSelector> TpSelectorForEfficiencyVsPhi;
  std::unique_ptr<TrackingParticleSelector> TpSelectorForEfficiencyVsPt;
  std::unique_ptr<TrackingParticleSelector> TpSelectorForEfficiencyVsVTXR;
  std::unique_ptr<TrackingParticleSelector> TpSelectorForEfficiencyVsVTXZ;

  std::unique_ptr<RecoTrackSelectorBase> trackSelectorVsEta;
  std::unique_ptr<RecoTrackSelectorBase> trackSelectorVsPhi;
  std::unique_ptr<RecoTrackSelectorBase> trackSelectorVsPt;

  std::unique_ptr<GenParticleCustomSelector> generalGpSelector;
  std::unique_ptr<GenParticleCustomSelector> GpSelectorForTechnicalEfficiency;
  std::unique_ptr<GenParticleCustomSelector> GpSelectorForEfficiencyVsEta;
  std::unique_ptr<GenParticleCustomSelector> GpSelectorForEfficiencyVsPhi;
  std::unique_ptr<GenParticleCustomSelector> GpSelectorForEfficiencyVsPt;
  std::unique_ptr<GenParticleCustomSelector> GpSelectorForEfficiencyVsVTXR;
  std::unique_ptr<GenParticleCustomSelector> GpSelectorForEfficiencyVsVTXZ;

  double minEta, maxEta;
  int nintEta;
  bool useFabsEta;
  double minPt, maxPt;
  int nintPt;
  bool useInvPt;
  bool useLogPt;
  double minHit, maxHit;
  int nintHit;
  double minPu, maxPu;
  int nintPu;
  double minLayers, maxLayers;
  int nintLayers;
  double minPhi, maxPhi;
  int nintPhi;
  double minDxy, maxDxy;
  int nintDxy;
  double minDz, maxDz;
  int nintDz;
  double dxyDzZoom;
  double minVertpos, maxVertpos;
  int nintVertpos;
  bool useLogVertpos;
  double minZpos, maxZpos;
  int nintZpos;
  double mindr, maxdr;
  int nintdr;
  double mindrj, maxdrj;
  int nintdrj;
  double minChi2, maxChi2;
  int nintChi2;
  double minDeDx, maxDeDx;
  int nintDeDx;
  double minVertcount, maxVertcount;
  int nintVertcount;
  double minTracks, maxTracks;
  int nintTracks;
  double minPVz, maxPVz;
  int nintPVz;
  double minMVA, maxMVA;
  int nintMVA;

  const bool doSeedPlots_;
  const bool doMTDPlots_;
  const bool doDzPVcutPlots_;

  //
  double ptRes_rangeMin, ptRes_rangeMax;
  int ptRes_nbin;
  double phiRes_rangeMin, phiRes_rangeMax;
  int phiRes_nbin;
  double cotThetaRes_rangeMin, cotThetaRes_rangeMax;
  int cotThetaRes_nbin;
  double dxyRes_rangeMin, dxyRes_rangeMax;
  int dxyRes_nbin;
  double dzRes_rangeMin, dzRes_rangeMax;
  int dzRes_nbin;

  double maxDzpvCum;
  int nintDzpvCum;
  double maxDzpvsigCum;
  int nintDzpvsigCum;

  std::vector<std::string> seedingLayerSetNames;
  using SeedingLayerId =
      std::tuple<SeedingLayerSetsBuilder::SeedingLayerId, bool>;  // last bool for strip mono (true) or not (false)
  using SeedingLayerSetId = std::array<SeedingLayerId, 4>;
  std::map<SeedingLayerSetId, unsigned int> seedingLayerSetToBin;
};

#endif
