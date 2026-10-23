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
    float sinPhi;
    float cosPhi;
    float pt;
  };

  // Result of the circle-circle closest-approach search in x-y.
  struct CircleCirclePoints {
    bool intersecting;  // true:  (ax, ay)/(bx, by) are the two crossings
                        // false: (x1, y1)/(x2, y2) are the per-track nearest points
    float ax, ay;
    float bx, by;
    float x1, y1;
    float x2, y2;
  };

  // Each track's own (x, y, z) near-point for the circle-circle
  // closest-approach estimate - port of geometry._nearestPoints. In the
  // intersecting regime the two candidate crossings are disambiguated by
  // z-consistency: propagate both tracks to both candidates and keep the
  // one with smaller |z1 - z2|.
  struct TrackNearPoints {
    float x1, y1, z1;
    float x2, y2, z2;
  };

  // Geometric (non-iterative) common-vertex estimate for a track pair -
  // port of geometry.estimatedVertex: xy-midpoint and z-average of the two
  // tracks' near-points, plus their full 3D separation (the quantity a
  // max3DDistance-style cut applies to).
  struct VertexEstimate {
    float x, y, z;
    float distance3D;
  };

  // Minimal 3-vector for the refinement math. Host+device (unlike the
  // acc-templated functions), so tests can use it anywhere.
  struct Point3D {
    float x, y, z;
  };

  // Signed transverse arc length from trk's reference point to (targetX,
  // targetY), both assumed to lie on trk's own circle:
  //   l = -sign(q) * r * dTheta
  // with dTheta the shortest signed angle from refAngle to the target's
  // angular position on the circle. Inverse of the position map
  // alpha(l) = refAngle - sign(q) * l / r.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float arcLengthTo(TAcc const& acc, HelixGeom const& trk, float targetX, float targetY) {
    const float targetAngle = alpaka::math::atan2(acc, targetY - trk.cy, targetX - trk.cx);
    const float dTheta = cms::alpakatools::reducePhiRange(acc, targetAngle - trk.refAngle);
    return -trk.q * trk.r * dTheta;
  }

  // z reached by propagating trk along its own helix to (targetX, targetY).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float propagateZ(TAcc const& acc, HelixGeom const& trk, float targetX, float targetY) {
    return trk.refZ + arcLengthTo(acc, trk, targetX, targetY) * trk.cotTheta;
  }

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

  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE TrackNearPoints nearestPoints(TAcc const& acc,
                                                               HelixGeom const& t1,
                                                               HelixGeom const& t2) {
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

  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexEstimate estimatedVertex(TAcc const& acc,
                                                                HelixGeom const& t1,
                                                                HelixGeom const& t2) {
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

  // -----------------------------------------------------------------------
  // Helix-helix closest-approach refinement (Gauss-Newton + Levenberg-
  // Marquardt)
  //
  // The LM loop itself (iteration budget, convergence tolerance, step
  // accept/reject, damping adaptation, dual starting point) is deliberately
  // NOT here: that is control-flow/policy and belongs to the calling kernel.
  // What lives here are the pure math primitives the loop is built from,
  // each individually testable against the Python reference:
  //
  //   helixPositionAndDerivatives  position map x(l) and derivatives dx/dl + d2x/dl2
  //   pairNormalEquationsSolve     one damped 2x2 normal-equations solve
  //   clampArcLength               the maxArcLength bound on |l|
  //   smartStartArcLengths         seed (l1, l2) from the geometric estimate
  //   refinedVertex                epilogue: converged (l1, l2) -> VertexEstimate
  //
  // Per-iteration math (mirrors refineVertex): with residual r = p1 - p2 and
  // j_i = dr/dl_i, solve the LM-damped normal equations
  //
  //   [ a(1+lam)      b     ] [dl1]   [ j1 . r ]
  //   [    b       c(1+lam) ] [dl2] = [ j2 . r ]
  //
  // with a = |j1|^2, c = |j2|^2, b = j1.j2, and step l_i -> l_i - dl_i (note
  // the sign convention: the solve returns the gradient-direction step,
  // which the caller subtracts from l).
  // -----------------------------------------------------------------------

  // LM policy defaults, taken from the Python reference (maxIter=4 validated
  // against 30 on real data with indistinguishable results). The kernel owns
  // the loop and may override these; they live here only to keep the magic
  // numbers in one place.
  inline constexpr int defaultMaxRefineIter = 4;
  inline constexpr float defaultLmDamping = 1e-3f;
  inline constexpr float defaultLmGrow = 5.f;
  inline constexpr float defaultLmShrink = 3.f;

  // Bound on |l| (default: the CMS tracker's outer radius scale). No real
  // track's vertex displacement can exceed it, but the unconstrained
  // least-squares problem has genuine far-away minima for near-tangent
  // high-R circle pairs - see refineVertex's docstring in the Python
  // reference for the confirmed failure mode this guards against.
  inline constexpr float defaultMaxArcLength = 100.f;

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE Point3D operator-(Point3D a, Point3D b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  }

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE float dot(Point3D a, Point3D b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE float norm2(Point3D v) { return dot(v, v); }

  // Position on trk's helix at signed transverse arc length l, and the
  // derivative d(position)/dl there - one trig evaluation shared between the
  // two. Position map: alpha(l) = refAngle - q * l / r (inverse of
  // arcLengthTo), x = cx + r cos(alpha), y = cy + r sin(alpha),
  // z = refZ + l * cotTheta.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void helixPositionAndDerivatives(
      TAcc const& acc, HelixGeom const& trk, float l, Point3D& pos, Point3D& dpos, Point3D& ddpos) {
    const float alpha = trk.refAngle - trk.q * l / trk.r;
    const float sinAlpha = alpaka::math::sin(acc, alpha);
    const float cosAlpha = alpaka::math::cos(acc, alpha);
    pos = {trk.cx + trk.r * cosAlpha, trk.cy + trk.r * sinAlpha, trk.refZ + l * trk.cotTheta};
    dpos = {trk.q * sinAlpha, -trk.q * cosAlpha, trk.cotTheta};
    ddpos = {-cosAlpha / trk.r, -sinAlpha / trk.r, 0.f};
  }

  // Position only (epilogue use, when the derivative is not needed).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE Point3D helixPosition(TAcc const& acc, HelixGeom const& trk, float l) {
    Point3D pos, dpos, ddpos;
    helixPositionAndDerivatives(acc, trk, l, pos, dpos, ddpos);
    return pos;
  }

  // First derivative d(x)/dl of the helix position map.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE Point3D helixFirstDerivative(TAcc const& acc, HelixGeom const& trk, float l) {
    Point3D pos, dpos, ddpos;
    helixPositionAndDerivatives(acc, trk, l, pos, dpos, ddpos);
    return dpos;
  }

  // Second derivative d^2(x)/dl^2 of the helix position map. The transverse
  // part is the position rotated by 90 degrees scaled by 1/r (no new trig);
  // z is linear in l so its second derivative vanishes.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE Point3D helixSecondDerivative(TAcc const& acc, HelixGeom const& trk, float l) {
    Point3D pos, dpos, ddpos;
    helixPositionAndDerivatives(acc, trk, l, pos, dpos, ddpos);
    return ddpos;
  }

  // One LM-damped 2x2 normal-equations solve for the pair refinement. Takes
  // the two tracks' position derivatives and the residual r = pos1 - pos2,
  // returns the step (dl1, dl2) to SUBTRACT from (l1, l2). j2 = -dpos2 is
  // handled internally (residual convention r = pos1 - pos2).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void pairNormalEquationsSolve(
      TAcc const& acc, Point3D dpos1, Point3D dpos2, Point3D residual, float lam, float& dl1, float& dl2) {
    const float a = norm2(dpos1) * (1.f + lam);
    const float c = norm2(dpos2) * (1.f + lam);
    const float b = -dot(dpos1, dpos2);
    const float rhs0 = dot(dpos1, residual);
    const float rhs1 = -dot(dpos2, residual);
    // Cauchy-Schwarz guarantees a*c - b*b >= 0; it vanishes only for exactly
    // parallel local directions (genuinely degenerate pairs), where the
    // guard plus lam > 0 keep the step finite instead of overshooting.
    const float det = alpaka::math::max(acc, a * c - b * b, 1e-12f);
    dl1 = (c * rhs0 - b * rhs1) / det;
    dl2 = (a * rhs1 - b * rhs0) / det;
  }

  // Clip an arc length to [-maxArcLength, maxArcLength] - applied to every
  // starting and trial point before it is evaluated or accepted.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float clampArcLength(TAcc const& acc,
                                                      float l,
                                                      float maxArcLength = defaultMaxArcLength) {
    return alpaka::math::min(acc, alpaka::math::max(acc, l, -maxArcLength), maxArcLength);
  }

  // Seed (l1, l2) from the geometric circle-circle estimate - refineVertex's
  // "smart" start. (The "dumb" start is simply l1 = l2 = 0 and needs no
  // helper.)
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void smartStartArcLengths(TAcc const& acc,
                                                           HelixGeom const& t1,
                                                           HelixGeom const& t2,
                                                           float& l1,
                                                           float& l2,
                                                           float maxArcLength = defaultMaxArcLength) {
    const auto np = nearestPoints(acc, t1, t2);
    l1 = clampArcLength(acc, arcLengthTo(acc, t1, np.x1, np.y1), maxArcLength);
    l2 = clampArcLength(acc, arcLengthTo(acc, t2, np.x2, np.y2), maxArcLength);
  }

  // Epilogue of the refinement: converged arc lengths -> VertexEstimate
  // (midpoint + true 3D separation). Unlike estimatedVertex's circle
  // crossings, the two tracks' refined points need not coincide in xy.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexEstimate
  refinedVertex(TAcc const& acc, HelixGeom const& t1, HelixGeom const& t2, float l1, float l2) {
    const Point3D p1 = helixPosition(acc, t1, l1);
    const Point3D p2 = helixPosition(acc, t2, l2);
    VertexEstimate vtx;
    vtx.x = 0.5f * (p1.x + p2.x);
    vtx.y = 0.5f * (p1.y + p2.y);
    vtx.z = 0.5f * (p1.z + p2.z);
    vtx.distance3D = alpaka::math::sqrt(acc, norm2(p1 - p2));
    return vtx;
  }

  // -----------------------------------------------------------------------
  // Unweighted LM-damped Gauss-Newton helix-helix refinement
  //
  // Minimizes f(l1, l2) = |x1(l1) - x2(l2)|^2 (squared 3D distance between
  // the two helices) over the two signed arc lengths, starting from (l1, l2).
  //
  // Gauss-Newton: the Hessian is approximated by J^T J (the residual-times-
  // curvature term is dropped), so no second derivatives are needed and the
  // 2x2 normal matrix is always positive semi-definite. Levenberg-Marquardt:
  // the diagonal is scaled by (1 + lambda), and lambda is adapted by
  // accept/reject of every trial step:
  //   - trial point has smaller |r|^2  -> accept, lambda /= lmShrink
  //   - otherwise                      -> reject, lambda *= lmGrow, retry
  // Every trial arc length is clamped to |l| <= maxArcLength (see
  // defaultMaxArcLength for the failure mode this prevents).
  //
  // One iteration costs one pairNormalEquationsSolve and one evaluation of
  // both helices at the trial point (position + derivative share one trig
  // evaluation per track); an accepted trial's values are reused for the
  // next iteration, so there is no redundant evaluation.
  // -----------------------------------------------------------------------

  // Tunable policy; defaults mirror the Python reference (see the constants
  // above). tolArcLength is the |dl| below which an accepted step counts as
  // converged (cm; 1e-4 cm = 1 micron, below tracker resolution).
  struct PairRefineConfig {
    int maxIter = defaultMaxRefineIter;
    float lmDamping = defaultLmDamping;
    float lmGrow = defaultLmGrow;
    float lmShrink = defaultLmShrink;
    float maxArcLength = defaultMaxArcLength;
    float tolArcLength = 1e-4f;
    float lmDampingMax = 1e4f;  // give up (stalled) when rejects push lambda beyond this
  };

  struct PairRefineResult {
    float l1, l2;     // refined arc lengths
    float distance2;  // |x1(l1) - x2(l2)|^2 at the returned point
    int iterations;   // trial steps taken (accepted + rejected)
    bool converged;   // an accepted step fell below tolArcLength
  };

  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE PairRefineResult refinePairGaussNewton(TAcc const& acc,
                                                                        HelixGeom const& t1,
                                                                        HelixGeom const& t2,
                                                                        float l1Start,
                                                                        float l2Start,
                                                                        PairRefineConfig const& cfg = {}) {
    float l1 = clampArcLength(acc, l1Start, cfg.maxArcLength);
    float l2 = clampArcLength(acc, l2Start, cfg.maxArcLength);

    Point3D p1, d1, dd1, p2, d2, dd2;  // dd* (second derivative) is unused by Gauss-Newton
    helixPositionAndDerivatives(acc, t1, l1, p1, d1, dd1);
    helixPositionAndDerivatives(acc, t2, l2, p2, d2, dd2);
    Point3D r = p1 - p2;
    float cost = norm2(r);

    float lam = cfg.lmDamping;
    int iter = 0;
    bool converged = false;

    while (iter < cfg.maxIter) {
      ++iter;

      // damped step, to be subtracted from (l1, l2)
      float dl1, dl2;
      pairNormalEquationsSolve(acc, d1, d2, r, lam, dl1, dl2);

      const float n1 = clampArcLength(acc, l1 - dl1, cfg.maxArcLength);
      const float n2 = clampArcLength(acc, l2 - dl2, cfg.maxArcLength);

      Point3D np1, nd1, ndd1, np2, nd2, ndd2;
      helixPositionAndDerivatives(acc, t1, n1, np1, nd1, ndd1);
      helixPositionAndDerivatives(acc, t2, n2, np2, nd2, ndd2);
      const Point3D nr = np1 - np2;
      const float nCost = norm2(nr);

      if (nCost < cost) {
        // accept; judge convergence on the step actually applied (after clamping)
        const float applied1 = alpaka::math::abs(acc, n1 - l1);
        const float applied2 = alpaka::math::abs(acc, n2 - l2);
        l1 = n1;
        l2 = n2;
        p1 = np1;
        d1 = nd1;
        p2 = np2;
        d2 = nd2;
        r = nr;
        cost = nCost;
        lam = lam / cfg.lmShrink;
        if (applied1 < cfg.tolArcLength && applied2 < cfg.tolArcLength) {
          converged = true;
          break;
        }
      } else {
        // reject: more damping, retry from the same point
        lam = lam * cfg.lmGrow;
        if (lam > cfg.lmDampingMax)
          break;
      }
    }
    return {l1, l2, cost, iter, converged};
  }

  // Two-start variant: run from the geometric ("smart") seed and from
  // l1 = l2 = 0 (the reference points, "dumb" start), keep the lower final
  // distance. ~2x cost; use where wrong-basin pairs matter.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE PairRefineResult
  refinePairDualStart(TAcc const& acc, HelixGeom const& t1, HelixGeom const& t2, PairRefineConfig const& cfg = {}) {
    float s1, s2;
    smartStartArcLengths(acc, t1, t2, s1, s2, cfg.maxArcLength);
    const auto smart = refinePairGaussNewton(acc, t1, t2, s1, s2, cfg);
    const auto dumb = refinePairGaussNewton(acc, t1, t2, 0.f, 0.f, cfg);
    return (dumb.distance2 < smart.distance2) ? dumb : smart;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::geometry

#endif  // RecoVertex_Vega_plugins_alpaka_VegaGeometry_h
