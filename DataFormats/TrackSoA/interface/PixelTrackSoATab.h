#ifndef DataFormats_TrackSoA_PixelTrackSoATab_h
#define DataFormats_TrackSoA_PixelTrackSoATab_h

#include <Eigen/Core>
#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include <cmath>

class PixelTrackSoATab {
    // Indices in the state vector
    static constexpr int iPhi = 0;
    static constexpr int iTip = 1;
    static constexpr int iQoverP = 2;
    static constexpr int iCotanTheta = 3;
    static constexpr int iZip = 4;

    // Indices in the covariance matrix
    static constexpr int cPhi = 0;
    static constexpr int cTip = 5;
    static constexpr int cInvPt = 9;
    static constexpr int cZip = 14;

public:
    PixelTrackSoATab() = default;
    PixelTrackSoATab(const reco::TrackSoAConstView& tracks, int idx)
        : tracks_(&tracks), idx_(idx) {}

    // ---- kinematics ----
    float pt() const { return tracks_->pt(idx_); }
    float eta() const { return tracks_->eta(idx_); }
    float phi() const { return tracks_->state(idx_)(iPhi); }

    int charge() const {
        float qop = tracks_->state(idx_)(iQoverP);
        return (qop > 0.f) - (qop < 0.f);
    }

    // ---- impact parameters (beam spot frame) ----
    float dxy() const { return tracks_->state(idx_)(iTip); }
    float dz()  const { return tracks_->state(idx_)(iZip); }

    // ---- errors (from SoA covariance layout) ----
    float dxyError() const { return std::sqrt(tracks_->covariance(idx_)(cTip)); }
    float dzError()  const { return std::sqrt(tracks_->covariance(idx_)(cZip)); }
    float phiError() const { return std::sqrt(tracks_->covariance(idx_)(cPhi)); }

    float ptError() const {
        float pt = tracks_->pt(idx_);
        return std::sqrt(tracks_->covariance(idx_)(cInvPt)) * pt * pt;
    }

    // ---- fit quality ----
    float chi2() const {
        return tracks_->chi2(idx_) * ndof();
    }

    int ndof() const {
        int nh = nHits();
        return 2 * nh - 5;
    }

    int nHits() const {
        return reco::nHits(*tracks_, idx_);
    }

private:
    const reco::TrackSoAConstView* tracks_ = nullptr;
    int idx_ = -1;
};

#endif
