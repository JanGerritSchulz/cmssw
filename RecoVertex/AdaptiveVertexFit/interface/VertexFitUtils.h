#ifndef RecoVertex_AdaptiveVertexFitter_interface_VertexFitUtils_h
#define RecoVertex_AdaptiveVertexFitter_interface_VertexFitUtils_h

namespace vertexfit {

  enum CovIdx : int {
    kVarPhi = 0,
    kCovPhiDxy = 1,
    kCovPhiQOverPt = 2,
    kVarDxy = 5,
    kCovDxyQOverPt = 6,
    kVarQOverPt = 9,
    kVarCotTheta = 12,
    kCovCotThetaDz = 13,
    kVarDz = 14
    // slots 3,4,7,8,10,11 are the exact cross-block zeros for pixel tracks from broken line
  };

  template <size_t NTRKS>
  struct TrackParams {
    std::array<float, NTRKS> phi, dxy, qOverPt, cotTheta, dz;
  };

  template <size_t NTRKS>
  struct TrackCov {
    std::array<std::array<float, 15>, NTRKS> packed;
  };

  template <size_t NTRKS>
  struct TrackGeomParams {
    std::array<float, NTRKS> cx, cy, r, q, refAngle;
  };
}  // namespace vertexfit

#endif
