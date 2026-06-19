#ifndef Validation_RecoVertex_SVTrackQualityBundle_h
#define Validation_RecoVertex_SVTrackQualityBundle_h

// Package:    Validation/RecoVertex
// Class:      SVTrackQualityBundle
//
/**\class SVTrackQualityBundle Validation/RecoVertex/interface/SVTrackQualityBundle.h

 Description: Bundle of monitoring elements for the SecondaryVertexAnalyzer,
              grouping all histograms for the track content quality at one
              place (e.g. track purity or the number of shared tracks).

 Original Author: Jan Schulz
*/

#include "DQMServices/Core/interface/DQMBookingHelpers.h"
// SVTrackQualityBundle.h
class SVTrackQualityBundle {
  using IBooker = dqm::reco::DQMStore::IBooker;

public:
  void fill(double nTracksInRecoSV, double nMatchedTracksInSimSV, double purity, double efficiency, double nShared) {
    h_purity_vs_nTracks->Fill(nTracksInRecoSV, purity);
    h_efficiency_vs_nTracks->Fill(nMatchedTracksInSimSV, efficiency);
    h_nShared_vs_nTracks->Fill(nMatchedTracksInSimSV, nShared);
    h_purity->Fill(purity);
    h_efficiency->Fill(efficiency);
    h_nShared->Fill(nShared);
  }

  void bookHistograms(IBooker &ibooker, int nTracksNBins, double nTracksMin, double nTracksMax) {
    h_purity_vs_nTracks = dqm::booking::book2DIfLogX(
        ibooker,
        false,
        "trackPurity_vs_nTracks",
        "Track purity vs N(tracks in RecoSV);N(tracks in RecoSV);Track purity = nSharedTracks(RecoSV, SimSV)",
        nTracksNBins,
        nTracksMin,
        nTracksMax,
        50,
        0.,
        1.0001);
    h_efficiency_vs_nTracks = dqm::booking::book2DIfLogX(
        ibooker,
        false,
        "trackEfficiency_vs_nTracks",
        "Track efficiency vs N(matched RecoTracks in SimSV);N(matched RecoTracks in SimSV);Track efficiency = "
        "nSharedTracks(RecoSV, SimSV) / nMatchedRecoTracks(SimSV)",
        nTracksNBins,
        nTracksMin,
        nTracksMax,
        50,
        0.,
        1.0001);
    h_nShared_vs_nTracks = dqm::booking::book2DIfLogX(
        ibooker,
        false,
        "nSharedTracks_vs_nTracks",
        "N(shared tracks) vs N(matched RecoTracks in SimSV);N(matched RecoTracks in SimSV);N shared tracks",
        nTracksNBins,
        nTracksMin,
        nTracksMax,
        20,
        -0.5,
        19.5);

    h_purity = ibooker.book1D(
        "trackPurity",
        "Track purity per RecoSV;Purity = nSharedTracks(RecoSV, SimSV) / nTracks(RecoSV);Sim-matched RecoSVs",
        50,
        0.,
        1.);
    h_efficiency = ibooker.book1D("trackEfficiency",
                                  "Track efficiency per SimSV;Efficiency = nSharedTracks(RecoSV, SimSV) / "
                                  "nMatchedRecoTracks(SimSV);Reco-matched SimSVs",
                                  50,
                                  0.,
                                  1.);
    h_nShared = ibooker.book1D("nSharedTracks", "N shared tracks;N shared tracks;Entries", 20, -0.5, 19.5);
  }

private:
  dqm::reco::MonitorElement *h_purity_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_efficiency_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_nShared_vs_nTracks = nullptr;
  dqm::reco::MonitorElement *h_purity = nullptr;
  dqm::reco::MonitorElement *h_efficiency = nullptr;
  dqm::reco::MonitorElement *h_nShared = nullptr;
};

#endif  // Validation_RecoVertex_SVTrackQualityBundle_h
