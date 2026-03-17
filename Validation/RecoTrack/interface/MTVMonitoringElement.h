#ifndef Validation_RecoTracker_MTVMonitoringElement_h
#define Validation_RecoTracker_MTVMonitoringElement_h

// Package:    Validation/RecoTracker
// Class:      MTVMonitoringElement
//
/**\class MTVMonitoringElement Validation/RecoTracker/MTVMonitoringElement.h
 Description: Bundle of monitoring elements for the MultiTrackValidator to have all histograms 
              for the same variable at one unique place 
              (e.g. the various plots over eta like reco, sim, assocSimToReco, assocRecoToSim, ...)
*/
//
// Original Author:  Jan Schulz (2026)

#include "DQMServices/Core/interface/MonitorElement.h"
#include <TF1.h>

// Wrapper functions taken from original MTV code
namespace {
  typedef dqm::reco::DQMStore DQMStore;

  void BinLogX(TH1* h) {
    TAxis* axis = h->GetXaxis();
    int bins = axis->GetNbins();

    float from = axis->GetXmin();
    float to = axis->GetXmax();
    float width = (to - from) / bins;
    std::vector<float> new_bins(bins + 1, 0);

    for (int i = 0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
    }
    axis->Set(bins, new_bins.data());
  }

  void BinLogY(TH1* h) {
    TAxis* axis = h->GetYaxis();
    int bins = axis->GetNbins();

    float from = axis->GetXmin();
    float to = axis->GetXmax();
    float width = (to - from) / bins;
    std::vector<float> new_bins(bins + 1, 0);

    for (int i = 0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
    }
    axis->Set(bins, new_bins.data());
  }

  template <typename... Args>
  dqm::reco::MonitorElement* make1DIfLogX(DQMStore::IBooker& ibook, bool logx, Args&&... args) {
    auto h = std::make_unique<TH1F>(std::forward<Args>(args)...);
    if (logx)
      BinLogX(h.get());
    const auto& name = h->GetName();
    return ibook.book1D(name, h.release());
  }

  template <typename... Args>
  dqm::reco::MonitorElement* makeProfileIfLogX(DQMStore::IBooker& ibook, bool logx, Args&&... args) {
    auto h = std::make_unique<TProfile>(std::forward<Args>(args)...);
    if (logx)
      BinLogX(h.get());
    const auto& name = h->GetName();
    return ibook.bookProfile(name, h.release());
  }

  template <typename... Args>
  dqm::reco::MonitorElement* make2DIfLogX(DQMStore::IBooker& ibook, bool logx, Args&&... args) {
    auto h = std::make_unique<TH2F>(std::forward<Args>(args)...);
    if (logx)
      BinLogX(h.get());
    const auto& name = h->GetName();
    return ibook.book2D(name, h.release());
  }

  template <typename... Args>
  dqm::reco::MonitorElement* make2DIfLogY(DQMStore::IBooker& ibook, bool logy, Args&&... args) {
    auto h = std::make_unique<TH2F>(std::forward<Args>(args)...);
    if (logy)
      BinLogY(h.get());
    const auto& name = h->GetName();
    return ibook.book2D(name, h.release());
  }
}  // namespace

// -----------------------------------------------------------------------------------------------------------------

class MTVMonitoringElement {
public:
  MTVMonitoringElement() {}
  virtual ~MTVMonitoringElement() = default;

  template <typename... Args>
  void fillSimTrackHistos(const bool isMatched, const bool isReconstructable, Args... args) {
    h_sim->Fill(args...);
    if (isMatched)
      h_assocSimToReco->Fill(args...);
    if (isReconstructable) {
      h_reconstructableSim->Fill(args...);
      if (isMatched)
        h_assocReconstructableSimToReco->Fill(args...);
    }
  }

  template <typename... Args>
  void fillRecoHistos(const bool isMatched,
                      const bool isSelected,
                      const bool isDuplicate,
                      const bool isPileup,
                      const bool isChargeMatched,
                      Args... args) {
    h_reco->Fill(args...);
    if (isMatched)
      h_assocRecoToSim->Fill(args...);
    if (isSelected)
      h_selectedReco->Fill(args...);
    if (isDuplicate)
      h_duplicate->Fill(args...);
    if (isPileup)
      h_pileup->Fill(args...);
    if (!isChargeMatched)
      h_chargemisid->Fill(args...);
  }

  // booking a generic type of histograms
  template <typename BookFunc, typename... Args>
  void bookGeneric(BookFunc bookFunc,
                   DQMStore::IBooker& ibooker,
                   const bool logScale,
                   const bool bookSimTrackHistos,
                   const bool bookRecoHistos,
                   const std::string& name,
                   const std::string& xlabel,
                   const std::string& ylabel,
                   const int nBins,
                   const double valMin,
                   const double valMax,
                   Args... args) {
    const std::string& xylabels = "; " + xlabel + "; " + ylabel;

    if (bookSimTrackHistos) {
      h_sim = bookFunc(ibooker,
                       logScale,
                       ("num_simul_" + name).c_str(),
                       ("N of simulated tracks (for efficiency); " + xylabels).c_str(),
                       nBins,
                       valMin,
                       valMax,
                       args...);
      h_assocSimToReco =
          bookFunc(ibooker,
                   logScale,
                   ("num_assoc(simToReco)_" + name).c_str(),
                   ("N of simulated tracks associated to a reco track (for efficiency); " + xylabels).c_str(),
                   nBins,
                   valMin,
                   valMax,
                   args...);
      h_reconstructableSim =
          bookFunc(ibooker,
                   logScale,
                   ("num_reconstructableSim_" + name).c_str(),
                   ("N of reconstructable simulated tracks (for technical efficiency); " + xylabels).c_str(),
                   nBins,
                   valMin,
                   valMax,
                   args...);
      h_assocReconstructableSimToReco = bookFunc(
          ibooker,
          logScale,
          ("num_assoc(reconstructableSimToReco)_" + name).c_str(),
          ("N of reconstructable simulated tracks associated to a reco track (for technical efficiency); " + xylabels)
              .c_str(),
          nBins,
          valMin,
          valMax,
          args...);
    }

    if (bookRecoHistos) {
      h_reco = bookFunc(ibooker,
                        logScale,
                        ("num_reco_" + name).c_str(),
                        ("N of reconstructed tracks; " + xylabels).c_str(),
                        nBins,
                        valMin,
                        valMax,
                        args...);
      h_selectedReco = bookFunc(ibooker,
                                logScale,
                                ("num_reco2_" + name).c_str(),
                                ("N of selected reconstructed tracks; " + xylabels).c_str(),
                                nBins,
                                valMin,
                                valMax,
                                args...);
      h_assocRecoToSim =
          bookFunc(ibooker,
                   logScale,
                   ("num_assoc(recoToSim)_" + name).c_str(),
                   ("N of reconstructed tracks matched to a simulated track (for fake rate); " + xylabels).c_str(),
                   nBins,
                   valMin,
                   valMax,
                   args...);
      h_duplicate = bookFunc(
          ibooker,
          logScale,
          ("num_duplicate_" + name).c_str(),
          ("N of reconstructed duplicated tracks (matched to multi-assoc simulated track); " + xylabels).c_str(),
          nBins,
          valMin,
          valMax,
          args...);
      h_chargemisid = bookFunc(
          ibooker,
          logScale,
          ("num_chargemisid_" + name).c_str(),
          ("N of reconstructed tracks matched to a simulated track but with mis-identified charge; " + xylabels).c_str(),
          nBins,
          valMin,
          valMax,
          args...);
      h_pileup = bookFunc(ibooker,
                          logScale,
                          ("num_pileup_" + name).c_str(),
                          ("N of reconstructed tracks matched to a simulated pileup track; " + xylabels).c_str(),
                          nBins,
                          valMin,
                          valMax,
                          args...);
    }
  }

  // booking a 1D histogram with option to set log scale on x
  template <typename... Args>
  void book1DIfLogX(DQMStore::IBooker& ibooker,
                    const bool logScale,
                    const bool bookSimTrackHistos,
                    const bool bookRecoHistos,
                    Args&&... args) {
    bookGeneric(
        [](auto& ib, bool log, auto&&... innerArgs) {
          return make1DIfLogX(ib, log, std::forward<decltype(innerArgs)>(innerArgs)...);
        },
        ibooker,
        logScale,
        bookSimTrackHistos,
        bookRecoHistos,
        std::forward<Args>(args)...);
  }

  // booking a standard 1D histogram
  template <typename... Args>
  void book1D(DQMStore::IBooker& ibooker, Args... args) {
    book1DIfLogX(ibooker, false, std::forward<Args>(args)...);
  }

  // booking a 1D histogram with log scale on x
  template <typename... Args>
  void book1DLogX(DQMStore::IBooker& ibooker, Args... args) {
    book1DIfLogX(ibooker, true, std::forward<Args>(args)...);
  }

  // booking a 2D histogram with option to set log scale on x
  template <typename... Args>
  void book2DIfLogX(DQMStore::IBooker& ibooker,
                    const bool logScale,
                    const bool bookSimTrackHistos,
                    const bool bookRecoHistos,
                    Args&&... args) {
    bookGeneric(
        [](auto& ib, bool log, auto&&... innerArgs) {
          return make2DIfLogX(ib, log, std::forward<decltype(innerArgs)>(innerArgs)...);
        },
        ibooker,
        logScale,
        bookSimTrackHistos,
        bookRecoHistos,
        std::forward<Args>(args)...);
  }

  // booking a standard 2D histogram
  template <typename... Args>
  void book2D(DQMStore::IBooker& ibooker, Args... args) {
    book2DIfLogX(ibooker, false, std::forward<Args>(args)...);
  }

  // booking a 2D histogram with log scale on x
  template <typename... Args>
  void book2DLogX(DQMStore::IBooker& ibooker, Args... args) {
    book2DIfLogX(ibooker, true, std::forward<Args>(args)...);
  }

  // booking a 2D histogram with option to set log scale on y
  template <typename... Args>
  void book2DIfLogY(DQMStore::IBooker& ibooker,
                    const bool logScale,
                    const bool bookSimTrackHistos,
                    const bool bookRecoHistos,
                    Args&&... args) {
    bookGeneric(
        [](auto& ib, bool log, auto&&... innerArgs) {
          return make2DIfLogY(ib, log, std::forward<decltype(innerArgs)>(innerArgs)...);
        },
        ibooker,
        logScale,
        bookSimTrackHistos,
        bookRecoHistos,
        std::forward<Args>(args)...);
  }

  // booking a 2D histogram with log scale on y
  template <typename... Args>
  void book2DLogY(DQMStore::IBooker& ibooker, Args... args) {
    book2DIfLogY(ibooker, true, std::forward<Args>(args)...);
  }

private:
  // histograms
  dqm::reco::MonitorElement* h_reco = nullptr;
  dqm::reco::MonitorElement* h_selectedReco = nullptr;
  dqm::reco::MonitorElement* h_sim = nullptr;
  dqm::reco::MonitorElement* h_reconstructableSim = nullptr;
  dqm::reco::MonitorElement* h_assocSimToReco = nullptr;
  dqm::reco::MonitorElement* h_assocRecoToSim = nullptr;
  dqm::reco::MonitorElement* h_assocReconstructableSimToReco = nullptr;
  dqm::reco::MonitorElement* h_duplicate = nullptr;
  dqm::reco::MonitorElement* h_chargemisid = nullptr;
  dqm::reco::MonitorElement* h_pileup = nullptr;
};

#endif