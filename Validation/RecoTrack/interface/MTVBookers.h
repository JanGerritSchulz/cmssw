#ifndef Validation_RecoTrack_MTVBookers_h
#define Validation_RecoTrack_MTVBookers_h

// Package:    Validation/RecoTrack
//
/**
 Description: Bunch of wrapper functions to easily book histograms in the MTV
*/

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

#endif
