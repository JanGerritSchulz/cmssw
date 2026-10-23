#ifndef RecoVertex_Vega_plugins_alpaka_VegaVertexFit_h
#define RecoVertex_Vega_plugins_alpaka_VegaVertexFit_h

#include <alpaka/alpaka.hpp>

#include "VegaGeometry.h"

#define VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE 1

// Weighted, covariance-aware two-track vertex fit - structured
// for one-thread-per-pair GPU execution (no shared memory, no sync, all
// state in registers).
//
// Formulation:
// the vertex position V is an explicit free parameter; each track carries
// its own nuisance parameter l_i (signed transverse arc length along its
// helix). Given V the tracks decouple completely, so each iteration is:
//   1. per track, independently: Newton+LM solve for the l_i minimizing the
//      weighted residual against the current V (profileArcLength), with the
//      position-space weight from covariance transport (positionWeight).
//   2. closed-form weighted least squares for V (weightedPoint).
// After convergence the per-track weights get the profiled-arc-length
// (Schur-complement) correction (correctWeightForProfiledArcLength) before
// the final vertex covariance and chi2 are computed - skipping this was the
// dominant source of a ~40x transverse-covariance underestimate in the
// Python reference.
//
// Everything here is per-track or per-pair scalar math on fixed-size
// objects (3-vectors, symmetric 3x3 matrices) - the N-track generalization
// reuses every per-track helper unchanged and only widens the sums in
// weightedPoint / the final covariance.
//
// Covariance conventions (must match the SoA track parametrization):
// the track's 5-parameter covariance is block-diagonal: a 3x3 block for
// {phi, dxy, qOverPt} (circle fit) and a 2x2 block for {cotTheta, dz}
// (line fit), with NO cross-block terms (confirmed on real data in soa.py).
// TrackCovariance stores the 8 independent entries; positionWeight
// exploits the block structure instead of assembling a dense 5x5.

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::vertexfit {

  using geometry::dot;
  using geometry::HelixGeom;
  using geometry::norm2;
  using geometry::Point3D;

  // -----------------------------------------------------------------------
  // Small fixed-size linear algebra (symmetric 3x3, in registers)
  // -----------------------------------------------------------------------

  // Symmetric 3x3 matrix, stored as the 6 independent entries.
  struct SymMat3 {
    float xx, xy, xz, yy, yz, zz;
  };
  struct SymMat3d {
    double xx, xy, xz, yy, yz, zz;
  };

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE SymMat3 operator+(SymMat3 a, SymMat3 b) {
    return {a.xx + b.xx, a.xy + b.xy, a.xz + b.xz, a.yy + b.yy, a.yz + b.yz, a.zz + b.zz};
  }

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE Point3D mul(SymMat3 const& m, Point3D v) {
    return {m.xx * v.x + m.xy * v.y + m.xz * v.z,
            m.xy * v.x + m.yy * v.y + m.yz * v.z,
            m.xz * v.x + m.yz * v.y + m.zz * v.z};
  }

  // Outer product v v^T as a symmetric matrix.
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE SymMat3 outerProduct(Point3D v) {
    return {v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.y, v.y * v.z, v.z * v.z};
  }

  // Inverse of a symmetric 3x3 via cofactors. The weights handled here are
  // sums of positive-semidefinite terms and well-conditioned except for
  // genuinely degenerate (near-parallel) track pairs, which the LM damping
  // in profileArcLength already flags; a small pivot guard keeps the
  // determinant away from exact zero.
  // returns false for singular / non-finite input (NaN fails the comparisons)
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE bool inverseSym(TAcc const& acc, SymMat3 const& m, SymMat3& out) {
    const float s = (m.xx + m.yy + m.zz) * (1.0 / 3.0);
    if (!(s > 0.0))
      return false;
    const float is = 1.0 / s;
    const float xx = m.xx * is, xy = m.xy * is, xz = m.xz * is, yy = m.yy * is, yz = m.yz * is, zz = m.zz * is;
    const float cxx = yy * zz - yz * yz, cxy = xz * yz - xy * zz, cxz = xy * yz - xz * yy;
    const float cyy = xx * zz - xz * xz, cyz = xy * xz - xx * yz, czz = xx * yy - xy * xy;
    const float det = xx * cxx + xy * cxy + xz * cxz;
    if (!(det > 1e-20))
      return false;
    const float k = is / det;
    out = {float(cxx * k), float(cxy * k), float(cxz * k), float(cyy * k), float(cyz * k), float(czz * k)};
    return true;
  }

  // -----------------------------------------------------------------------
  // Track covariance (block-diagonal 5-parameter helix covariance)
  // -----------------------------------------------------------------------

  // Independent entries of the block-diagonal track covariance:
  //   {phi, dxy, qOverPt}  3x3 block (circle fit)
  //   {cotTheta, dz}       2x2 block (line fit)
  // No cross-block terms exist (see module docstring).
  struct TrackCovariance {
    float varPhi, covPhiDxy, covPhiQOverPt, varDxy, covDxyQOverPt, varQOverPt;
    float varCotTheta, covCotThetaDz, varDz;
  };

  // -----------------------------------------------------------------------
  // Covariance transport: helix-parameter covariance -> position-space
  // weight at the current arc length
  // -----------------------------------------------------------------------

  // Position-space covariance Sigma = J C J^T of the track's helix position
  // at arc length l, and its inverse (the weight W used in the fit).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE SymMat3
  positionWeight(TAcc const& acc, HelixGeom const& trk, TrackCovariance const& cov, float l) {
    const float alpha = trk.refAngle - trk.q * l / trk.r;
    const float sinAlpha = alpaka::math::sin(acc, alpha);
    const float cosAlpha = alpaka::math::cos(acc, alpha);

    // Jacobian rows (x, y, z) x cols (phi, dxy, qOverPt, cotTheta, dz).
    const float j[3][5] = {{/* x */ -trk.cy - trk.r * sinAlpha,
                            trk.sinPhi,
                            trk.pt * (-trk.r * trk.sinPhi - trk.q * trk.r * cosAlpha + l * sinAlpha),
                            0.f,
                            0.f},
                           {/* y */ trk.cx + trk.r * cosAlpha,
                            -trk.cosPhi,
                            trk.pt * (trk.r * trk.cosPhi - trk.q * trk.r * sinAlpha - l * cosAlpha),
                            0.f,
                            0.f},
                           {/* z */ 0.f, 0.f, 0.f, l, 1.f}};

    // Sigma = J C J^T, exploiting the block-diagonal C: the {phi,dxy,qOverPt}
    // block contributes to x,y only, the {cotTheta,dz} block to z only.
    SymMat3 sigma{};
    // circle-fit block (indices 0..2)
    const float c3[3][3] = {{cov.varPhi, cov.covPhiDxy, cov.covPhiQOverPt},
                            {cov.covPhiDxy, cov.varDxy, cov.covDxyQOverPt},
                            {cov.covPhiQOverPt, cov.covDxyQOverPt, cov.varQOverPt}};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        sigma.xx += j[0][a] * c3[a][b] * j[0][b];
        sigma.xy += j[0][a] * c3[a][b] * j[1][b];
        // sigma.xz += j[0][a] * c3[a][b] * j[2][b]; // j[2][0...2] = 0
        sigma.yy += j[1][a] * c3[a][b] * j[1][b];
        // sigma.yz += j[1][a] * c3[a][b] * j[2][b]; // j[2][0...2] = 0
        // sigma.zz += j[2][a] * c3[a][b] * j[2][b]; // j[2][0...2] = 0
      }
    // line-fit block (indices 3..4)
    // Full math:
    // const float c2[2][2] = {{cov.varCotTheta, cov.covCotThetaDz}, {cov.covCotThetaDz, cov.varDz}};
    // for (int a = 0; a < 2; ++a)
    //   for (int b = 0; b < 2; ++b) {
    //     sigma.xx += j[0][a + 3] * c2[a][b] * j[0][b + 3];
    //     sigma.xy += j[0][a + 3] * c2[a][b] * j[1][b + 3];
    //     sigma.xz += j[0][a + 3] * c2[a][b] * j[2][b + 3];
    //     sigma.yy += j[1][a + 3] * c2[a][b] * j[1][b + 3];
    //     sigma.yz += j[1][a + 3] * c2[a][b] * j[2][b + 3];
    //     sigma.zz += j[2][a + 3] * c2[a][b] * j[2][b + 3];
    //   }
    // However, doing the math explicitly, you will just find
    sigma.zz = j[2][3] * j[2][3] * cov.varCotTheta + 2 * j[2][3] * cov.covCotThetaDz + cov.varDz;

    SymMat3 weights;
    inverseSym(acc, sigma, weights);
    // if (!inverseValid)
    //   printf("VegaVertexFit: Warning invalid weights matrix due to failing matrix inversion!");
    return weights;
  }

  // -----------------------------------------------------------------------
  // Per-track arc-length profiling (the inner 1-D Newton+LM solve)
  // -----------------------------------------------------------------------

  // Find the arc length minimizing track i's weighted residual against the
  // current vertex estimate V: a few Newton steps with light LM damping,
  // warm-started from lStart (the previous iteration's value, or the
  // geometric-estimate seed on the first call).
  //
  //   g(l) = d^T W (x(l) - V)
  //   H(l) = d^T W d + (x(l) - V)^T W d'(l)
  //
  // Do NOT replace this with the closed-form transverse-only projection
  // (alpha = atan2(Vy-Cy, Vx-Cx)): the Python reference measured up to ~1 cm
  // arc-length error and ~94% chi2 inflation from that shortcut on displaced
  // tracks. It is a fine seed, not a solution.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float profileArcLength(TAcc const& acc,
                                                        HelixGeom const& trk,
                                                        SymMat3 const& w,
                                                        Point3D v,
                                                        float lStart,
                                                        int maxIter = 3,
                                                        float lmDamping = 1e-3f) {
    float l = lStart;
    for (int iter = 0; iter < maxIter; ++iter) {
      Point3D pos, dpos, ddpos;
      geometry::helixPositionAndDerivatives(acc, trk, l, pos, dpos, ddpos);
      const Point3D residual = pos - v;

      const float grad = dot(dpos, mul(w, residual));
      const float hess1 = dot(dpos, mul(w, dpos));
      const float hess2 = dot(residual, mul(w, ddpos));
      float hess = hess1 + hess2;
      // guard against negative/zero hess (use hess1 only instead)
      if (hess <= 0)
        hess = hess1;
      const float step = grad / alpaka::math::max(acc, hess * (1.f + lmDamping), 1e-12f);
      l = geometry::clampArcLength(acc, l - step);
      if (alpaka::math::abs(acc, step) < 1e-6f)
        break;
    }
    return l;
  }

  // -----------------------------------------------------------------------
  // Profiled-arc-length (Schur-complement) weight correction
  // -----------------------------------------------------------------------

  // Correct a track's converged weight for the fact that its arc length l is
  // itself fitted (a profiled nuisance parameter), not independently known.
  // The direction along the track's own tangent is unmeasured by that track,
  // so its weight in that direction is removed (partially restored by the
  // nuisance uncertainty in l):
  //
  //   W_eff = W - (W d)(W d)^T / (d^T W d + 1/var_l)
  //
  // with d = dx/dl at the converged l and var_l the nuisance variance from
  // the profiling solve (1/H at convergence). Applied ONCE to the final
  // weight, not inside the iteration - see vertexfit.py /
  // _correctWeightForProfiledArcLength for the derivation and the ~40x
  // covariance underestimate this fixes.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE SymMat3
  correctWeightForProfiledArcLength(TAcc const& acc, HelixGeom const& trk, SymMat3 const& w, float l, float varL) {
    Point3D pos, dpos, ddpos;
    geometry::helixPositionAndDerivatives(acc, trk, l, pos, dpos, ddpos);
    const Point3D wd = mul(w, dpos);
    const float denom = alpaka::math::max(acc, dot(dpos, wd) + 1.f / alpaka::math::max(acc, varL, 1e-12f), 1e-12f);
    const SymMat3 corr = outerProduct(wd);
    const float s = 1.f / denom;
    return {w.xx - s * corr.xx,
            w.xy - s * corr.xy,
            w.xz - s * corr.xz,
            w.yy - s * corr.yy,
            w.yz - s * corr.yz,
            w.zz - s * corr.zz};
  }

  // -----------------------------------------------------------------------
  // Closed-form weighted-least-squares vertex update
  // -----------------------------------------------------------------------

  // Exact WLS vertex position given both tracks' optimal points and weights:
  //   V = (W1 + W2)^-1 (W1 x1 + W2 x2)
  // One 3x3 solve, no iteration. The N-track generalization just widens the
  // two sums.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE Point3D
  weightedPoint(TAcc const& acc, Point3D x1, SymMat3 const& w1, Point3D x2, SymMat3 const& w2) {
    const SymMat3 a = w1 + w2;
    const Point3D b{w1.xx * x1.x + w1.xy * x1.y + w1.xz * x1.z + w2.xx * x2.x + w2.xy * x2.y + w2.xz * x2.z,
                    w1.xy * x1.x + w1.yy * x1.y + w1.yz * x1.z + w2.xy * x2.x + w2.yy * x2.y + w2.yz * x2.z,
                    w1.xz * x1.x + w1.yz * x1.y + w1.zz * x1.z + w2.xz * x2.x + w2.yz * x2.y + w2.zz * x2.z};
    SymMat3 ia;
    inverseSym(acc, a, ia);
    return mul(ia, b);
  }

  // -----------------------------------------------------------------------
  // Fit result
  // -----------------------------------------------------------------------

  struct VertexFitResult {
    Point3D position;
    SymMat3 covariance;
    float chi2;
    float l1, l2;  // converged arc lengths (kept for downstream use / debugging)
    bool valid;
  };

  // -----------------------------------------------------------------------
  // Two-track fit driver: alternates per-track profiling with the
  // closed-form vertex update, then applies the profiled-arc-length
  // correction and computes the final covariance and chi2.
  //
  // Seeded from the geometric estimate (estimatedVertex) for V and
  // smartStartArcLengths for l1/l2 - both already close for genuine pairs,
  // so a small fixed iteration budget suffices (mirrors maxIter=4 in the
  // Python reference).
  // -----------------------------------------------------------------------
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexFitResult fitTwoTrackVertex(TAcc const& acc,
                                                                   HelixGeom const& t1,
                                                                   HelixGeom const& t2,
                                                                   TrackCovariance const& cov1,
                                                                   TrackCovariance const& cov2,
                                                                   int maxIter = 4,
                                                                   int maxIterProfiling = 3,
                                                                   float lmDampingProfiling = 1e-3f) {
    // seed
    const auto seed = geometry::estimatedVertex(acc, t1, t2);  // FIXME: avoid conversion maybe?
    Point3D v{seed.x, seed.y, seed.z};
    float l1, l2;
    geometry::smartStartArcLengths(acc, t1, t2, l1, l2);  // FIXME: this completely repeats the lines above!!!

#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
    printf(
        "fitTwoTrackVertex: track1 (cx, cy, r, q, refAngle, refZ, cotTheta, sinPhi, cosPhi, pt) = "
        "(%f, %f, %f, %f, %f, %f, %f, %f, %f, %f)\n",
        t1.cx, t1.cy, t1.r, t1.q, t1.refAngle, t1.refZ, t1.cotTheta, t1.sinPhi, t1.cosPhi, t1.pt);
    printf(
        "fitTwoTrackVertex: track2 (cx, cy, r, q, refAngle, refZ, cotTheta, sinPhi, cosPhi, pt) = "
        "(%f, %f, %f, %f, %f, %f, %f, %f, %f, %f)\n",
        t2.cx, t2.cy, t2.r, t2.q, t2.refAngle, t2.refZ, t2.cotTheta, t2.sinPhi, t2.cosPhi, t2.pt);
    printf(
        "fitTwoTrackVertex: cov1 (varPhi, covPhiDxy, covPhiQOverPt, varDxy, covDxyQOverPt, varQOverPt, "
        "varCotTheta, covCotThetaDz, varDz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)\n",
        cov1.varPhi, cov1.covPhiDxy, cov1.covPhiQOverPt, cov1.varDxy, cov1.covDxyQOverPt, cov1.varQOverPt,
        cov1.varCotTheta, cov1.covCotThetaDz, cov1.varDz);
    printf(
        "fitTwoTrackVertex: cov2 (varPhi, covPhiDxy, covPhiQOverPt, varDxy, covDxyQOverPt, varQOverPt, "
        "varCotTheta, covCotThetaDz, varDz) = (%g, %g, %g, %g, %g, %g, %g, %g, %g)\n",
        cov2.varPhi, cov2.covPhiDxy, cov2.covPhiQOverPt, cov2.varDxy, cov2.covDxyQOverPt, cov2.varQOverPt,
        cov2.varCotTheta, cov2.covCotThetaDz, cov2.varDz);
    printf("fitTwoTrackVertex: seed V = (%f, %f, %f), (l1, l2) = (%f, %f)\n", v.x, v.y, v.z, l1, l2);
#endif

    // alternate per-track profiling and the closed-form vertex update
    for (int iter = 0; iter < maxIter; ++iter) {
      const SymMat3 w1 = positionWeight(acc, t1, cov1, l1);
      const SymMat3 w2 = positionWeight(acc, t2, cov2, l2);
      l1 = profileArcLength(acc, t1, w1, v, l1, maxIterProfiling, lmDampingProfiling);
      l2 = profileArcLength(acc, t2, w2, v, l2, maxIterProfiling, lmDampingProfiling);
      const Point3D x1 = geometry::helixPosition(acc, t1, l1);
      const Point3D x2 = geometry::helixPosition(acc, t2, l2);
      v = weightedPoint(acc, x1, w1, x2, w2);

#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
      printf(
          "fitTwoTrackVertex: iter %d  (l1, l2) = (%f, %f)  x1 = (%f, %f, %f)  x2 = (%f, %f, %f)  "
          "V = (%f, %f, %f)\n",
          iter, l1, l2, x1.x, x1.y, x1.z, x2.x, x2.y, x2.z, v.x, v.y, v.z);
#endif
    }

    // final weights with the profiled-arc-length correction
    const SymMat3 w1 = positionWeight(acc, t1, cov1, l1);
    const SymMat3 w2 = positionWeight(acc, t2, cov2, l2);
    // nuisance variance from the profiling Hessian at convergence
    Point3D pos1, dpos1, ddpos1, pos2, dpos2, ddpos2;
    geometry::helixPositionAndDerivatives(acc, t1, l1, pos1, dpos1, ddpos1);
    geometry::helixPositionAndDerivatives(acc, t2, l2, pos2, dpos2, ddpos2);
    const float varL1 = 1.f / alpaka::math::max(acc, dot(dpos1, mul(w1, dpos1)), 1e-12f);
    const float varL2 = 1.f / alpaka::math::max(acc, dot(dpos2, mul(w2, dpos2)), 1e-12f);
    const SymMat3 w1eff = correctWeightForProfiledArcLength(acc, t1, w1, l1, varL1);
    const SymMat3 w2eff = correctWeightForProfiledArcLength(acc, t2, w2, l2, varL2);

    // final vertex position, covariance, chi2
    const Point3D x1 = geometry::helixPosition(acc, t1, l1);
    const Point3D x2 = geometry::helixPosition(acc, t2, l2);
    v = weightedPoint(acc, x1, w1eff, x2, w2eff);
    SymMat3 covariance;
    inverseSym(acc, w1eff + w2eff, covariance);
    const Point3D r1 = x1 - v;
    const Point3D r2 = x2 - v;
    const float chi2 = dot(r1, mul(w1eff, r1)) + dot(r2, mul(w2eff, r2));

#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
    printf(
        "fitTwoTrackVertex: done  V = (%f, %f, %f)  chi2 = %f  (l1, l2) = (%f, %f)  "
        "cov (xx, yy, zz) = (%g, %g, %g)\n",
        v.x, v.y, v.z, chi2, l1, l2, covariance.xx, covariance.yy, covariance.zz);
#endif

    return {v, covariance, chi2, l1, l2, true};
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::vertexfit

#endif  // RecoVertex_Vega_plugins_alpaka_VegaVertexFit_h
