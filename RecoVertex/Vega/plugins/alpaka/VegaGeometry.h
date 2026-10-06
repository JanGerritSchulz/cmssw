#ifndef RecoVertex_Vega_plugins_alpaka_VegaGeometry_h
#define RecoVertex_Vega_plugins_alpaka_VegaGeometry_h

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"

// Helix geometry for track pairs: circle-circle closest approach in the
// transverse plane plus z from propagating each track along its own helix.
//
// Conventions (must match KernelBuildTrackExtra / TrackExtraSoA):
//   - circle center (cx, cy), unsigned radius r, charge sign q, and
//     refAngle = angular position of the track's reference point (PCA to
//     beamspot) on its own circle, reduced to (-pi, pi]
//   - refZ = z of the reference point (TrackSoA's zip/dz)
//   - cotTheta = dz/dl along the helix; identical to sinh(eta) for the SoA
//     parametrization, so eta itself is never needed here
//
// Empirical facts from the Python reference that motivate the design:
//   - ~99.6% of true same-vertex pairs have *intersecting* circles in x-y,
//     so that branch is the hot path and the disjoint regimes mostly matter
//     for background
//   - for intersecting circles the xy separation is zero by construction;
//     the discriminating power is carried by z-consistency, which is also
//     how the two candidate crossings are disambiguated (nearestPoints)

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::geometry {

  // Per-track helix geometry, bundled so the functions below never need to
  // know about SoA layouts. Fill from TrackSoA + TrackExtraSoA per track:
  //   cx, cy, r, q, refAngle  <- trksExtra[it]   (KernelBuildTrackExtra)
  //   refZ                    <- reco::zip(trks, it)
  //   cotTheta                <- trks[it].state()(3)
  struct HelixGeom {
    float cx;
    float cy;
    float r;
    float q;         // only the sign is ever used
    float refAngle;  // angular position of the reference point on the circle
    float refZ;      // z of the reference point
    float cotTheta;  // dz per unit signed transverse arc length
  };

  // Signed transverse arc length from trk's reference point to (targetX,
  // targetY), both assumed to lie on trk's own circle:
  //   l = -sign(q) * r * dTheta
  // with dTheta the shortest signed angle from refAngle to the target's
  // angular position on the circle. Inverse of the position map
  // theta(l) = refAngle - sign(q) * l / r.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float arcLengthTo(TAcc const& acc,
                                                   HelixGeom const& trk,
                                                   float targetX,
                                                   float targetY) {
    const float targetAngle = alpaka::math::atan2(acc, targetY - trk.cy, targetX - trk.cx);
    const float dTheta = cms::alpakatools::reducePhiRange(acc, targetAngle - trk.refAngle);
    return -trk.q * trk.r * dTheta;
  }

  // z reached by propagating trk along its own helix to (targetX, targetY).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float propagateZ(TAcc const& acc, HelixGeom const& trk, float targetX, float targetY) {
    return trk.refZ + arcLengthTo(acc, trk, targetX, targetY) * trk.cotTheta;
  }

  // Result of the circle-circle closest-approach search in x-y.
  struct CircleCirclePoints {
    bool intersecting;  // true:  (ax, ay)/(bx, by) are the two crossings
                        // false: (x1, y1)/(x2, y2) are the per-track nearest points
    float ax, ay;
    float bx, by;
    float x1, y1;
    float x2, y2;
  };

  // The two circles' closest-approach points in x-y, handling all three
  // geometric regimes (intersecting / externally separate / one contains the
  // other). In the intersecting case BOTH crossings are returned unresolved; 
  // disambiguation needs z-propagation and happens in nearestPoints.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE CircleCirclePoints closestPointsOnCircles(TAcc const& acc,
                                                                           HelixGeom const& t1,
                                                                           HelixGeom const& t2) {
    const float dx = t2.cx - t1.cx;
    const float dy = t2.cy - t1.cy;
    const float d = alpaka::math::sqrt(acc, dx * dx + dy * dy) + 1e-12f;  // guard against concentric circles
    const float ux = dx / d;
    const float uy = dy / d;

    CircleCirclePoints out;
    out.intersecting = (d >= alpaka::math::abs(acc, t1.r - t2.r)) && (d <= t1.r + t2.r);

    if (out.intersecting) {
      // standard two-circle intersection: chord at distance a from center 1, half-length h
      const float a = (d * d - t2.r * t2.r + t1.r * t1.r) / (2.f * d);
      const float h = alpaka::math::sqrt(acc, alpaka::math::max(acc, t1.r * t1.r - a * a, 0.f));
      const float mx = t1.cx + a * ux;
      const float my = t1.cy + a * uy;
      out.ax = mx - h * uy;
      out.ay = my + h * ux;
      out.bx = mx + h * uy;
      out.by = my - h * ux;
    } else if (d > t1.r + t2.r) {
      // externally separate: nearest points face each other along the center line
      out.x1 = t1.cx + t1.r * ux;
      out.y1 = t1.cy + t1.r * uy;
      out.x2 = t2.cx - t2.r * ux;
      out.y2 = t2.cy - t2.r * uy;
    } else {
      // one circle contains the other: nearest points on the same side
      const float sign = (t1.r >= t2.r) ? 1.f : -1.f;
      out.x1 = t1.cx + sign * t1.r * ux;
      out.y1 = t1.cy + sign * t1.r * uy;
      out.x2 = t2.cx + sign * t2.r * ux;
      out.y2 = t2.cy + sign * t2.r * uy;
    }
    return out;
  }

  // Each track's own (x, y, z) near-point for the circle-circle
  // closest-approach estimate - port of geometry._nearestPoints. In the
  // intersecting regime the two candidate crossings are disambiguated by
  // z-consistency: propagate both tracks to both candidates and keep the
  // one with smaller |z1 - z2|.
  struct TrackNearPoints {
    float x1, y1, z1;
    float x2, y2, z2;
  };

  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE TrackNearPoints nearestPoints(TAcc const& acc, HelixGeom const& t1, HelixGeom const& t2) {
    const auto cc = closestPointsOnCircles(acc, t1, t2);

    TrackNearPoints np;
    if (cc.intersecting) {
      const float zA1 = propagateZ(acc, t1, cc.ax, cc.ay);
      const float zA2 = propagateZ(acc, t2, cc.ax, cc.ay);
      const float zB1 = propagateZ(acc, t1, cc.bx, cc.by);
      const float zB2 = propagateZ(acc, t2, cc.bx, cc.by);
      const bool preferA = alpaka::math::abs(acc, zA1 - zA2) <= alpaka::math::abs(acc, zB1 - zB2);
      np.x1 = preferA ? cc.ax : cc.bx;
      np.y1 = preferA ? cc.ay : cc.by;
      np.z1 = preferA ? zA1 : zB1;
      np.x2 = np.x1;  // crossings lie on both circles: xy coincides
      np.y2 = np.y1;
      np.z2 = preferA ? zA2 : zB2;
    } else {
      np.x1 = cc.x1;
      np.y1 = cc.y1;
      np.z1 = propagateZ(acc, t1, cc.x1, cc.y1);
      np.x2 = cc.x2;
      np.y2 = cc.y2;
      np.z2 = propagateZ(acc, t2, cc.x2, cc.y2);
    }
    return np;
  }

  // Geometric (non-iterative) common-vertex estimate for a track pair -
  // port of geometry.estimatedVertex: xy-midpoint and z-average of the two
  // tracks' near-points, plus their full 3D separation (the quantity a
  // max3DDistance-style cut applies to).
  struct VertexEstimate {
    float x, y, z;
    float distance3D;
  };

  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexEstimate estimatedVertex(TAcc const& acc, HelixGeom const& t1, HelixGeom const& t2) {
    const auto np = nearestPoints(acc, t1, t2);

    VertexEstimate vtx;
    vtx.x = 0.5f * (np.x1 + np.x2);
    vtx.y = 0.5f * (np.y1 + np.y2);
    vtx.z = 0.5f * (np.z1 + np.z2);
    const float dx = np.x1 - np.x2;
    const float dy = np.y1 - np.y2;
    const float dz = np.z1 - np.z2;
    vtx.distance3D = alpaka::math::sqrt(acc, dx * dx + dy * dy + dz * dz);
    return vtx;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::geometry

#endif  // RecoVertex_Vega_plugins_alpaka_VegaGeometry_h
