#ifndef RecoVertex_Vega_plugins_alpaka_VegaVertexFit_h
#define RecoVertex_Vega_plugins_alpaka_VegaVertexFit_h

#include <alpaka/alpaka.hpp>

#include "VegaGeometry.h"
#include "VegaMath.h"

#ifndef VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
#define VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE 0
#endif

// Weighted, covariance-aware two-track vertex fit - structured for
// one-thread-per-pair GPU execution (no shared memory, no sync, all state in
// registers).
//
// Formulation:
// the vertex position V is an explicit free parameter; each track carries its
// own nuisance parameter l_i (signed transverse arc length along its helix).
// Given V the tracks decouple, so each iteration is:
//   1. per track, independently: 1-D Newton solve for the l_i minimizing the
//      weighted residual against the current V (profileArcLength), with the
//      position-space weight W_i from covariance transport (positionWeight).
//   2. closed-form weighted least squares for V (solveVertex), using the
//      Schur-corrected weights W_i^eff (correctWeightForProfiledArcLength).
//
// LESSONS BUILT IN (found by reproducing the first GPU convergence logs):
//
//  (a) The position weight W_i = (J C J^T)^-1 is evaluated at FIXED l, so it
//      contains no uncertainty along the track: 0.006 - 0.25 micron versus
//      tens of microns transverse. Using that W in the V update pins V to
//      each track's foot point along the tangent, so the alternation cannot
//      travel along the tracks at all (it stalled at l ~ 1e-4 cm where the
//      answer was l = 2.85 cm, in double precision). The V update therefore
//      uses W^eff = W - (W d)(W d)^T / h, which has no stiffness along d.
//      W itself is only used for the l-profiling (scale invariant there).
//  (b) h is the full Hessian d^T W d + r^T W d'. The earlier
//      1/varL = d^T W d form halves the denominator and removes only half of
//      the tangential weight (chi2 inflated by 5 - 450x).
//  (c) W^eff has rank 2; the sum of two of them is nearly singular for
//      near-parallel tracks (0.9 mrad crossing is real). All 3x3 algebra is
//      therefore done in double; in float the final solve returned garbage
//      (V ~ 1e6 cm, NaN) once (b) was fixed. Every inversion reports failure
//      and the fit then returns valid = false instead of garbage.
//  (d) The alternation is a polish, not a search: seed it with the Gauss-
//      Newton helix-helix refinement (pairfit::refinePairGaussNewton) and
//      call the seeded overload (VegaPairFit.h has the pair recipe).
//
// Covariance conventions (must match the SoA track parametrization, which
// stores {phi, tip, q/pT, cotTheta, zip}; tip = -dxy of reco::Track, see the
// Jacobian note in positionWeight): the track's 5-parameter covariance is
// block-diagonal: a 3x3 block for {phi, tip, qOverPt} (circle fit) and a 2x2
// block for {cotTheta, zip} (line fit), with NO cross-block terms. All 9
// independent entries are stored in HelixCovariance.
//
// Per-track helpers are per-track on purpose: the N-track generalization
// reuses them unchanged and only widens the sums in solveVertex.

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::vertexfit {

  using geometry::HelixCovariance;
  using geometry::HelixGeom;
  using math::SymMat3d;
  using math::Vec3d;

  // The helix evaluation at the fit precision (the geometry stays float).
  using HelixEvald = geometry::HelixEvalT<double>;

  // Floor added to the xy position variance (cm^2): (10 nm)^2. Physically
  // irrelevant (far below any resolution; and a variance added along the
  // direction that l absorbs does not change W^eff anyway) but it keeps the
  // 2x2 position block invertible when the track passes exactly through the
  // beamspot at l = 0, where J C J^T is exactly rank 1.
  inline constexpr double positionVarianceFloor = 1e-12;

  // -----------------------------------------------------------------------
  // Fit configuration and result
  // -----------------------------------------------------------------------

  struct VertexFitConfig {
    int maxIter = 6;           // outer alternations; a GN-seeded fit converges in ~4
    int maxIterProfiling = 3;  // Newton steps per track per alternation
    double lmDamping = 1e-3;   // LM damping of the profiling Newton step
    double maxArcLength = geometry::defaultMaxArcLength;
    double tolPosition = 1e-6;  // cm: stop when the vertex moves less than this
  };

  // Effective constraints: each W_eff has rank 2, so 2 tracks x 2 - 3
  // vertex parameters = 1. (N tracks: 2N - 3.)
  inline constexpr int twoTrackNdf = 1;

  // Fit result, templated on precision: the fit itself produces the double
  // (double) version; the public VertexFitResult (float, consistent with the
  // float HelixGeom / HelixCovariance inputs) is produced from it by a single
  // narrowing at the very end (toFloat). A float 3x3 covariance keeps ~1e-7
  // relative precision, enough for an output; a consumer that must re-invert
  // it or use it as a weight in a further fit should widen it first
  // (math::cast<double>) or, better, refit from the tracks.
  template <typename T>
  struct VertexFitResultT {
    math::Vec3<T> position;
    math::SymMat3<T> covariance;
    T chi2;
    int ndf;
    T l1, l2;        // converged arc lengths (seeds for downstream refits / debugging)
    int iterations;  // outer alternations actually performed
    bool valid;      // false: a matrix inversion failed; the other fields are not meaningful
  };

  using VertexFitResult = VertexFitResultT<float>;

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE VertexFitResult toFloat(VertexFitResultT<double> const& r) {
    return {math::cast<float>(r.position),
            math::cast<float>(r.covariance),
            float(r.chi2),
            r.ndf,
            float(r.l1),
            float(r.l2),
            r.iterations,
            r.valid};
  }

  // -----------------------------------------------------------------------
  // Covariance transport: helix-parameter covariance -> position-space
  // weight at the current arc length
  // -----------------------------------------------------------------------

  // Position-space covariance Sigma = J C J^T of the helix position at arc
  // length l (at fixed l!) and its inverse W. Returns false if Sigma is not
  // invertible (never for sane tracks, thanks to positionVarianceFloor).
  //
  // Jacobian columns (phi, tip, qOverPt, cotTheta, zip), derived from
  //   alpha = alpha0 - q l/R, alpha0 = phi + q pi/2, R = pT/(cB),
  //   C = (tip + q R) (sin phi, -cos phi),  x = C + R (cos alpha, sin alpha),
  //   z = zip + cotTheta l
  // (checked against finite differences in the unit test):
  //   dx/dphi   = -cy - R sin(alpha)        dy/dphi   = cx + R cos(alpha)
  //   dx/dtip   = sin(phi)                  dy/dtip   = -cos(phi)
  //   dx/dk     = pT (-R sin(phi) - q R cos(alpha) + l sin(alpha))
  //   dy/dk     = pT ( R cos(phi) - q R sin(alpha) - l cos(alpha))
  //   dz/dcot   = l,  dz/dzip = 1
  // Sigma is block-diagonal: the circle block feeds (x,y) only, the line
  // block feeds z only, so W is block-diagonal too and is inverted as a 2x2
  // block plus a scalar (no 3x3 cofactor cancellation).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE bool positionWeight(
      TAcc const&, HelixGeom const& trk, HelixCovariance const& cov, HelixEvald const& he, double l, SymMat3d& w) {
    const double q = trk.q, r = trk.r, pt = trk.pt;
    const double sinPhi = trk.sinPhi, cosPhi = trk.cosPhi;
    const double jx[3] = {-double(trk.cy) - r * he.sinA, sinPhi, pt * (-r * sinPhi - q * r * he.cosA + l * he.sinA)};
    const double jy[3] = {double(trk.cx) + r * he.cosA, -cosPhi, pt * (r * cosPhi - q * r * he.sinA - l * he.cosA)};

    const double c3[3][3] = {{cov.varPhi, cov.covPhiDxy, cov.covPhiQOverPt},
                             {cov.covPhiDxy, cov.varDxy, cov.covDxyQOverPt},
                             {cov.covPhiQOverPt, cov.covDxyQOverPt, cov.varQOverPt}};
    double cx[3], cy[3];  // C3 * jx, C3 * jy
    for (int a = 0; a < 3; ++a) {
      cx[a] = c3[a][0] * jx[0] + c3[a][1] * jx[1] + c3[a][2] * jx[2];
      cy[a] = c3[a][0] * jy[0] + c3[a][1] * jy[1] + c3[a][2] * jy[2];
    }
    const double sxx = jx[0] * cx[0] + jx[1] * cx[1] + jx[2] * cx[2] + positionVarianceFloor;
    const double sxy = jx[0] * cy[0] + jx[1] * cy[1] + jx[2] * cy[2];
    const double syy = jy[0] * cy[0] + jy[1] * cy[1] + jy[2] * cy[2] + positionVarianceFloor;
    const double szz = l * l * cov.varCotTheta + 2.0 * l * cov.covCotThetaDz + cov.varDz;

    const double det2 = sxx * syy - sxy * sxy;
    const double tr = sxx + syy;
    w = {0, 0, 0, 0, 0, 0};
    if (!(det2 > 1e-12 * tr * tr * 0.25) || !(szz > 0.0))
      return false;
    w = {syy / det2, -sxy / det2, 0.0, sxx / det2, 0.0, 1.0 / szz};
    return true;
  }

  // -----------------------------------------------------------------------
  // Per-track arc-length profiling (the inner 1-D Newton solve)
  // -----------------------------------------------------------------------

  // Find the arc length minimizing track i's weighted residual against the
  // current vertex estimate V: a few Newton steps with light LM damping,
  // warm-started from lStart.
  //
  //   g(l) = d^T W (x(l) - V)
  //   H(l) = d^T W d + (x(l) - V)^T W d'(l)
  //
  // If the exact Hessian is not positive, the Gauss-Newton part d^T W d is
  // used (always > 0). The step is invariant under a common scale of W, so
  // the huge tangential stiffness of W does not matter here.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE double profileArcLength(TAcc const& acc,
                                                         HelixGeom const& trk,
                                                         SymMat3d const& w,
                                                         Vec3d v,
                                                         double lStart,
                                                         int maxIter,
                                                         double lmDamping,
                                                         double maxArcLength) {
    double l = lStart;
    for (int iter = 0; iter < maxIter; ++iter) {
      const HelixEvald e = geometry::helixEval<double>(acc, trk, l);
      const Vec3d residual = e.pos - v;
      const double grad = dot(e.d, mul(w, residual));
      const double hess1 = dot(e.d, mul(w, e.d));
      double hess = hess1 + dot(residual, mul(w, e.dd));
      if (!(hess > 0.0))
        hess = hess1;
      if (!(hess > 0.0))
        break;  // degenerate weight: keep l
      const double step = grad / (hess * (1.0 + lmDamping));
      l = math::clamp(acc, l - step, maxArcLength);
      if (alpaka::math::abs(acc, step) < 1e-9)
        break;
    }
    return l;
  }

  // -----------------------------------------------------------------------
  // Profiled-arc-length (Schur-complement) weight correction
  // -----------------------------------------------------------------------

  // Eliminate the nuisance parameter l from the track's quadratic form
  // (x(l) - V)^T W (x(l) - V) exactly (Schur complement of the l-l entry):
  //
  //   W_eff = W - (W d)(W d)^T / h,   h = d^T W d + (x - V)^T W d'
  //
  // with d, d' at the converged l. h is the full Hessian of the profiling
  // problem (NOT d^T W d + 1/var_l with var_l = 1/(d^T W d), which gives
  // 2 d^T W d and removes only half the tangential weight). W_eff d = 0 up
  // to the small residual term: the track constrains V only across its own
  // direction. h is clamped to [d^T W d / 2, inf) so that an unlucky
  // residual term can never over-subtract.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE SymMat3d
  correctWeightForProfiledArcLength(TAcc const&, SymMat3d const& w, HelixEvald const& he, Vec3d v) {
    const Vec3d wd = mul(w, he.d);
    const double h1 = dot(he.d, wd);
    double h = h1 + dot(he.pos - v, mul(w, he.dd));
    if (!(h > 0.5 * h1))
      h = h1;
    if (!(h > 0.0))
      return w;
    const SymMat3d o = outerProduct(wd);
    const double s = 1.0 / h;
    return {w.xx - s * o.xx, w.xy - s * o.xy, w.xz - s * o.xz, w.yy - s * o.yy, w.yz - s * o.yz, w.zz - s * o.zz};
  }

  // -----------------------------------------------------------------------
  // Closed-form weighted-least-squares vertex update
  // -----------------------------------------------------------------------

  // Exact WLS vertex given both tracks' points and weights:
  //   V = (W1 + W2)^-1 (W1 x1 + W2 x2),   cov = (W1 + W2)^-1
  // One 3x3 inversion, no iteration. The N-track generalization just widens
  // the two sums. Returns false if (W1 + W2) is not invertible (exactly
  // parallel tracks).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE bool solveVertex(
      TAcc const& acc, Vec3d x1, SymMat3d const& w1, Vec3d x2, SymMat3d const& w2, Vec3d& v, SymMat3d& cov) {
    const Vec3d a1 = mul(w1, x1), a2 = mul(w2, x2);
    const Vec3d b{a1.x + a2.x, a1.y + a2.y, a1.z + a2.z};
    if (!inverseSym(acc, w1 + w2, cov))
      return false;
    v = mul(cov, b);
    return true;
  }

  // -----------------------------------------------------------------------
  // Two-track fit driver (seeded)
  //
  // l1Seed/l2Seed should come from pairfit::refinePairGaussNewton (or the
  // dual-start variant). The fit alternates profiling with the Schur-
  // corrected closed-form vertex update, then returns V, its covariance and
  // the chi2 of the converged state.
  //
  // fitTwoTrackVertexInternal returns the double (double) result, for callers
  // that still compare or combine fits (e.g. the pair stage picking the better
  // of two seeds); fitTwoTrackVertex below is the float-in / float-out
  // interface.
  // -----------------------------------------------------------------------
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexFitResultT<double> fitTwoTrackVertexInternal(TAcc const& acc,
                                                                                    HelixGeom const& t1,
                                                                                    HelixGeom const& t2,
                                                                                    HelixCovariance const& cov1,
                                                                                    HelixCovariance const& cov2,
                                                                                    double l1Seed,
                                                                                    double l2Seed,
                                                                                    VertexFitConfig const& cfg = {}) {
    VertexFitResultT<double> bad{{0, 0, 0}, {0, 0, 0, 0, 0, 0}, 0.0, twoTrackNdf, l1Seed, l2Seed, 0, false};

    double l1 = math::clamp(acc, l1Seed, cfg.maxArcLength);
    double l2 = math::clamp(acc, l2Seed, cfg.maxArcLength);

    HelixEvald e1 = geometry::helixEval<double>(acc, t1, l1);
    HelixEvald e2 = geometry::helixEval<double>(acc, t2, l2);
    Vec3d v{0.5 * (e1.pos.x + e2.pos.x), 0.5 * (e1.pos.y + e2.pos.y), 0.5 * (e1.pos.z + e2.pos.z)};

    SymMat3d w1, w2, eff1, eff2, vcov;
    if (!positionWeight(acc, t1, cov1, e1, l1, w1) || !positionWeight(acc, t2, cov2, e2, l2, w2)) {
#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
      printf("fitTwoTrackVertex: FAILED propagating covariance to seed position!!!\n");
#endif
      return bad;
    }

#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
    printf("fitTwoTrackVertex: track1 (cx, cy, r, q, refAngle, refZ, cotTheta) = (%f, %f, %f, %f, %f, %f, %f)\n",
           t1.cx,
           t1.cy,
           t1.r,
           t1.q,
           t1.refAngle,
           t1.refZ,
           t1.cotTheta);
    printf("fitTwoTrackVertex: track2 (cx, cy, r, q, refAngle, refZ, cotTheta) = (%f, %f, %f, %f, %f, %f, %f)\n",
           t2.cx,
           t2.cy,
           t2.r,
           t2.q,
           t2.refAngle,
           t2.refZ,
           t2.cotTheta);
    printf(
        "fitTwoTrackVertex: cov1 (varPhi, covPhiDxy, covPhiQOverPt, varDxy, covDxyQOverPt, varQOverPt, varCotTheta, "
        "covCotThetaDz, varDz) = (%f, %f, %f, %f, %f, %f, %f, %f, %f)\n",
        cov1.varPhi,
        cov1.covPhiDxy,
        cov1.covPhiQOverPt,
        cov1.varDxy,
        cov1.covDxyQOverPt,
        cov1.varQOverPt,
        cov1.varCotTheta,
        cov1.covCotThetaDz,
        cov1.varDz);
    printf(
        "fitTwoTrackVertex: cov2 (varPhi, covPhiDxy, covPhiQOverPt, varDxy, covDxyQOverPt, varQOverPt, varCotTheta, "
        "covCotThetaDz, varDz) = (%f, %f, %f, %f, %f, %f, %f, %f, %f)\n",
        cov2.varPhi,
        cov2.covPhiDxy,
        cov2.covPhiQOverPt,
        cov2.varDxy,
        cov2.covDxyQOverPt,
        cov2.varQOverPt,
        cov2.varCotTheta,
        cov2.covCotThetaDz,
        cov2.varDz);
    printf("fitTwoTrackVertex: seed (l1, l2) = (%f, %f)  V = (%f, %f, %f)\n", l1, l2, v.x, v.y, v.z);
#endif

    int iter = 0;
    for (; iter < cfg.maxIter; ++iter) {
      // 1. profile each track against the current V (decoupled)
      l1 = profileArcLength(acc, t1, w1, v, l1, cfg.maxIterProfiling, cfg.lmDamping, cfg.maxArcLength);
      l2 = profileArcLength(acc, t2, w2, v, l2, cfg.maxIterProfiling, cfg.lmDamping, cfg.maxArcLength);
      e1 = geometry::helixEval<double>(acc, t1, l1);
      e2 = geometry::helixEval<double>(acc, t2, l2);
      // weights at the refined arc lengths (also the profiling weights of the next iteration)
      if (!positionWeight(acc, t1, cov1, e1, l1, w1) || !positionWeight(acc, t2, cov2, e2, l2, w2)) {
#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
        printf("fitTwoTrackVertex: FAILED propagating covariance to new position!!!\n");
#endif
        return bad;
      }

      // 2. closed-form vertex update with the Schur-corrected weights
      eff1 = correctWeightForProfiledArcLength(acc, w1, e1, v);
      eff2 = correctWeightForProfiledArcLength(acc, w2, e2, v);
      const Vec3d vOld = v;
      if (!solveVertex(acc, e1.pos, eff1, e2.pos, eff2, v, vcov)) {
#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
        printf("fitTwoTrackVertex: FAILED solving for the vertex position update!!!\n");
#endif
        return bad;
      }

#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
      printf("fitTwoTrackVertex: iter %d  (l1, l2) = (%f, %f)  V = (%f, %f, %f)\n", iter, l1, l2, v.x, v.y, v.z);
#endif
      const Vec3d dv = v - vOld;
      if (dot(dv, dv) < cfg.tolPosition * cfg.tolPosition) {
        ++iter;
        break;
      }
    }

    // final, self-consistent state: W_eff at the final V, one more solve,
    // then covariance and chi2 from the same weights
    eff1 = correctWeightForProfiledArcLength(acc, w1, e1, v);
    eff2 = correctWeightForProfiledArcLength(acc, w2, e2, v);
    if (!solveVertex(acc, e1.pos, eff1, e2.pos, eff2, v, vcov)) {
#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
      printf("fitTwoTrackVertex: FAILED solving for the final vertex position!!!\n");
#endif
      return bad;
    }
    const Vec3d r1 = e1.pos - v;
    const Vec3d r2 = e2.pos - v;
    const double chi2 = dot(r1, mul(eff1, r1)) + dot(r2, mul(eff2, r2));
    if (!(chi2 >= 0.0)) {
#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
      printf("fitTwoTrackVertex: FAILED due to negative chi2!!!\n");
#endif
      return bad;
    }

#if VEGA_VERTEXFIT_MONITOR_FIT_CONVERGENCE
    printf("fitTwoTrackVertex: done after %d iters  V = (%f, %f, %f)  chi2 = %g  (l1, l2) = (%f, %f)\n",
           iter,
           v.x,
           v.y,
           v.z,
           chi2,
           l1,
           l2);
#endif
    return {v, vcov, chi2, twoTrackNdf, l1, l2, iter, true};
  }

  // Public two-track fit: float in, float out.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexFitResult fitTwoTrackVertex(TAcc const& acc,
                                                                   HelixGeom const& t1,
                                                                   HelixGeom const& t2,
                                                                   HelixCovariance const& cov1,
                                                                   HelixCovariance const& cov2,
                                                                   float l1Seed,
                                                                   float l2Seed,
                                                                   VertexFitConfig const& cfg = {}) {
    return toFloat(fitTwoTrackVertexInternal(acc, t1, t2, cov1, cov2, double(l1Seed), double(l2Seed), cfg));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::vertexfit

#endif  // RecoVertex_Vega_plugins_alpaka_VegaVertexFit_h
