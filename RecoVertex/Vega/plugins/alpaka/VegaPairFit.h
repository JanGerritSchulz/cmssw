#ifndef RecoVertex_Vega_plugins_alpaka_VegaPairFit_h
#define RecoVertex_Vega_plugins_alpaka_VegaPairFit_h

#include <alpaka/alpaka.hpp>

#include "VegaGeometry.h"
#include "VegaMath.h"
#include "VegaVertexFit.h"

#ifndef VEGA_PAIRFIT_MONITOR_FIT_CONVERGENCE
#define VEGA_PAIRFIT_MONITOR_FIT_CONVERGENCE 0
#endif

// Pair stage: seeding and the one-thread-per-pair recipe.
//
//   VegaMath.h       small linear algebra (Vec3, SymMat3, inverseSym)
//   VegaGeometry.h   helix geometry (HelixGeom, helixEval, circle-circle)
//   VegaVertexFit.h  weighted fit with explicit vertex (reusable for N tracks)
//   VegaPairFit.h    this file: everything that is specific to PAIRS
//                    - geometric and Gauss-Newton seeding of (l1, l2)
//                    - the recipe  GN seed -> fitTwoTrackVertex
//
// Gauss-Newton (float) finds the vertex "quickly"; the weighted fit (double
// internally) then provides the covariance and the chi2.

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::pairfit {

  using geometry::HelixCovariance;
  using geometry::HelixGeom;
  using geometry::HelixEval;
  using math::Vec3f;

  // Geometric / refined common-vertex estimate for a track pair: xy-midpoint
  // and z-average of the two tracks' near-points, plus their full 3D
  // separation (the quantity a max3DDistance-style cut applies to).
  struct VertexSeed {
    float x, y, z;
    float distance3D;
  };

  // Geometric (non-iterative) estimate from the circle-circle closest approach.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexSeed estimatedVertex(TAcc const& acc, HelixGeom const& t1, HelixGeom const& t2) {
    const auto np = geometry::nearestPoints(acc, t1, t2);

    VertexSeed vtx;
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
  // The building blocks (pairNormalEquationsSolve, smartStartArcLengths,
  // refinedVertex) are kept separate and individually testable against the
  // Python reference; refinePairGaussNewton is the loop built from them.
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
  // against 30 on real data with indistinguishable results).
  inline constexpr int defaultMaxRefineIter = 4;
  inline constexpr float defaultLmDamping = 1e-3f;
  inline constexpr float defaultLmGrow = 5.f;
  inline constexpr float defaultLmShrink = 3.f;

  // One LM-damped 2x2 normal-equations solve for the pair refinement. Takes
  // the two tracks' position derivatives and the residual r = pos1 - pos2,
  // returns the step (dl1, dl2) to SUBTRACT from (l1, l2). j2 = -dpos2 is
  // handled internally (residual convention r = pos1 - pos2).
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void pairNormalEquationsSolve(
      TAcc const& acc, Vec3f dpos1, Vec3f dpos2, Vec3f residual, float lam, float& dl1, float& dl2) {
    const float a = math::norm2(dpos1) * (1.f + lam);
    const float c = math::norm2(dpos2) * (1.f + lam);
    const float b = -math::dot(dpos1, dpos2);
    const float rhs0 = math::dot(dpos1, residual);
    const float rhs1 = -math::dot(dpos2, residual);
    // Cauchy-Schwarz guarantees a*c - b*b >= 0; it vanishes only for exactly
    // parallel local directions (genuinely degenerate pairs), where the
    // guard plus lam > 0 keep the step finite instead of overshooting.
    const float det = alpaka::math::max(acc, a * c - b * b, 1e-12f);
    dl1 = (c * rhs0 - b * rhs1) / det;
    dl2 = (a * rhs1 - b * rhs0) / det;
  }

  // Seed (l1, l2) from the geometric circle-crossing estimate -
  // refineVertex's "smart" start. (The "dumb" start is simply l1 = l2 = 0 and
  // needs no helper.)
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void smartStartArcLengths(TAcc const& acc,
                                                           HelixGeom const& t1,
                                                           HelixGeom const& t2,
                                                           float& l1,
                                                           float& l2,
                                                           float maxArcLength = geometry::defaultMaxArcLength) {
    const auto np = geometry::nearestPoints(acc, t1, t2);
    l1 = math::clamp(acc, geometry::arcLengthTo(acc, t1, np.x1, np.y1), maxArcLength);
    l2 = math::clamp(acc, geometry::arcLengthTo(acc, t2, np.x2, np.y2), maxArcLength);
  }

  // Epilogue of the refinement: converged arc lengths -> VertexSeed
  // (midpoint + true 3D separation). Unlike estimatedVertex's circle
  // crossings, the two tracks' refined points need not coincide in xy.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE VertexSeed
  refinedVertex(TAcc const& acc, HelixGeom const& t1, HelixGeom const& t2, float l1, float l2) {
    const Vec3f p1 = geometry::helixPosition(acc, t1, l1);
    const Vec3f p2 = geometry::helixPosition(acc, t2, l2);
    VertexSeed vtx;
    vtx.x = 0.5f * (p1.x + p2.x);
    vtx.y = 0.5f * (p1.y + p2.y);
    vtx.z = 0.5f * (p1.z + p2.z);
    vtx.distance3D = alpaka::math::sqrt(acc, math::norm2(p1 - p2));
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
  // geometry::defaultMaxArcLength for the failure mode this prevents).
  //
  // One iteration costs one pairNormalEquationsSolve and one evaluation of
  // both helices at the trial point; an accepted trial's values are reused
  // for the next iteration, so there is no redundant evaluation.
  // -----------------------------------------------------------------------

  // Tunable policy; defaults mirror the Python reference (see the constants
  // above). tolArcLength is the |dl| below which an accepted step counts as
  // converged (cm; 1e-4 cm = 1 micron, below tracker resolution).
  struct PairRefineConfig {
    int maxIter = defaultMaxRefineIter;
    float lmDamping = defaultLmDamping;
    float lmGrow = defaultLmGrow;
    float lmShrink = defaultLmShrink;
    float maxArcLength = geometry::defaultMaxArcLength;
    float tolArcLength = 1e-4f;
    float lmDampingMax = 1e4f;  // give up (stalled) when rejects push lambda beyond this
  };

  struct PairRefineResult {
    float l1, l2;     // refined arc lengths
    float distance2;  // |x1(l1) - x2(l2)|^2 at the returned point
    int iterations;   // trial steps taken (accepted + rejected)
    bool converged;   // an accepted step fell below tolArcLength
  };

  // GN policy for the pair stage: displaced pairs need many more iterations
  // than prompt ones (11-13 LM iterations from l = 0 for a 2.8 cm
  // displacement with a 0.9 mrad crossing angle; the library default of 4 was
  // validated on prompt pairs only).
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE PairRefineConfig defaultPairGaussNewtonConfig() {
    PairRefineConfig c;
    c.maxIter = 15;
    return c;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE PairRefineResult refinePairGaussNewton(TAcc const& acc,
                                                                        HelixGeom const& t1,
                                                                        HelixGeom const& t2,
                                                                        float l1Start,
                                                                        float l2Start,
                                                                        PairRefineConfig const& cfg = {}) {
    float l1 = math::clamp(acc, l1Start, cfg.maxArcLength);
    float l2 = math::clamp(acc, l2Start, cfg.maxArcLength);

#if VEGA_PAIRFIT_MONITOR_FIT_CONVERGENCE
    printf("refinePairGaussNewton: track1 (cx, cy, r, q, refAngle, refZ, cotTheta) = (%f, %f, %f, %f, %f, %f, %f)\n",
           t1.cx,
           t1.cy,
           t1.r,
           t1.q,
           t1.refAngle,
           t1.refZ,
           t1.cotTheta);
    printf("refinePairGaussNewton: track2 (cx, cy, r, q, refAngle, refZ, cotTheta) = (%f, %f, %f, %f, %f, %f, %f)\n",
           t2.cx,
           t2.cy,
           t2.r,
           t2.q,
           t2.refAngle,
           t2.refZ,
           t2.cotTheta);
    printf("refinePairGaussNewton: seed (l1, l2) = (%f, %f)\n", l1, l2);
#endif

    HelixEval e1 = geometry::helixEval(acc, t1, l1);  // e*.dd (second derivative) is unused by Gauss-Newton
    HelixEval e2 = geometry::helixEval(acc, t2, l2);
    Vec3f r = e1.pos - e2.pos;
    float cost = math::norm2(r);

    float lam = cfg.lmDamping;
    int iter = 0;
    bool converged = false;

    while (iter < cfg.maxIter) {
      ++iter;

      // damped step, to be subtracted from (l1, l2)
      float dl1, dl2;
      pairNormalEquationsSolve(acc, e1.d, e2.d, r, lam, dl1, dl2);

      const float n1 = math::clamp(acc, l1 - dl1, cfg.maxArcLength);
      const float n2 = math::clamp(acc, l2 - dl2, cfg.maxArcLength);

      const HelixEval trial1 = geometry::helixEval(acc, t1, n1);
      const HelixEval trial2 = geometry::helixEval(acc, t2, n2);
      const Vec3f nr = trial1.pos - trial2.pos;
      const float nCost = math::norm2(nr);

      if (nCost < cost) {
        // accept; judge convergence on the step actually applied (after clamping)
        const float applied1 = alpaka::math::abs(acc, n1 - l1);
        const float applied2 = alpaka::math::abs(acc, n2 - l2);
        l1 = n1;
        l2 = n2;
        e1 = trial1;
        e2 = trial2;
        r = nr;
        cost = nCost;
        lam = lam / cfg.lmShrink;
#if VEGA_PAIRFIT_MONITOR_FIT_CONVERGENCE
        printf(
            "refinePairGaussNewton: iter %d ACCEPT  (l1, l2) = (%f, %f)  p1 = (%f, %f, %f)  p2 = (%f, %f, %f)  "
            "dist = %f  lam = %g\n",
            iter,
            l1,
            l2,
            e1.pos.x,
            e1.pos.y,
            e1.pos.z,
            e2.pos.x,
            e2.pos.y,
            e2.pos.z,
            alpaka::math::sqrt(acc, cost),
            lam);
#endif
        if (applied1 < cfg.tolArcLength && applied2 < cfg.tolArcLength) {
          converged = true;
          break;
        }
      } else {
        // reject: more damping, retry from the same point
        lam = lam * cfg.lmGrow;
#if VEGA_PAIRFIT_MONITOR_FIT_CONVERGENCE
        printf("refinePairGaussNewton: iter %d REJECT  (l1, l2) = (%f, %f)  trial dist = %f  lam -> %g\n",
               iter,
               l1,
               l2,
               alpaka::math::sqrt(acc, nCost),
               lam);
#endif
        if (lam > cfg.lmDampingMax)
          break;
      }
    }
#if VEGA_PAIRFIT_MONITOR_FIT_CONVERGENCE
    printf("refinePairGaussNewton: done after %d iters, converged = %d, final (l1, l2) = (%f, %f), dist = %f\n",
           iter,
           converged,
           l1,
           l2,
           alpaka::math::sqrt(acc, cost));
#endif
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

  // -----------------------------------------------------------------------
  // The pair-stage recipe
  //   1. Gauss-Newton helix-helix refinement from the zero start (l1 = l2 = 0,
  //      i.e. the reference points = PCA to the beamspot).
  //   2. The weighted fit seeded with the GN result.
  //
  // tryGeometricStart = true additionally runs GN from the circle-crossing
  // ("smart") start, fits that too and keeps the valid fit with the lower
  // weighted chi2. It is OFF by default because in the toy test it did not
  // help: for near-collinear pairs (dphi < 60 mrad) the circle crossings of
  // near-tangent circles are ill-conditioned and the smart start left 7-8 %
  // of such pairs in a wrong basin (maha2 > 30), versus 1.5 % for the zero
  // start alone; for prompt and displaced pairs (0.5 - 3 cm) all choices
  // were equally good (< 0.02 % outliers). Selecting by the weighted chi2 or
  // by the GN distance made no real difference. A fit seeded at the true
  // vertex (oracle) reaches 0.1 %, so a better second seed would still pay.
  // -----------------------------------------------------------------------
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE vertexfit::VertexFitResult fitTwoTrackVertexGaussNewtonSeeded(
      TAcc const& acc,
      HelixGeom const& t1,
      HelixGeom const& t2,
      HelixCovariance const& cov1,
      HelixCovariance const& cov2,
      vertexfit::VertexFitConfig const& cfg = {},
      PairRefineConfig const& gnCfg = defaultPairGaussNewtonConfig(),
      bool tryGeometricStart = false) {
    using vertexfit::VertexFitResultT;

    const auto gnDumb = refinePairGaussNewton(acc, t1, t2, 0.f, 0.f, gnCfg);
    const VertexFitResultT<double> fitDumb =
        vertexfit::fitTwoTrackVertexInternal(acc, t1, t2, cov1, cov2, double(gnDumb.l1), double(gnDumb.l2), cfg);
    if (!tryGeometricStart)
      return vertexfit::toFloat(fitDumb);

    float s1, s2;
    smartStartArcLengths(acc, t1, t2, s1, s2, gnCfg.maxArcLength);
    const auto gnSmart = refinePairGaussNewton(acc, t1, t2, s1, s2, gnCfg);
    const double dl1 = double(gnSmart.l1) - double(gnDumb.l1);
    const double dl2 = double(gnSmart.l2) - double(gnDumb.l2);
    if (alpaka::math::abs(acc, dl1) < 1e-3 && alpaka::math::abs(acc, dl2) < 1e-3)
      return vertexfit::toFloat(fitDumb);  // same basin

    const VertexFitResultT<double> fitSmart =
        vertexfit::fitTwoTrackVertexInternal(acc, t1, t2, cov1, cov2, double(gnSmart.l1), double(gnSmart.l2), cfg);
    if (!fitDumb.valid)
      return vertexfit::toFloat(fitSmart);
    if (!fitSmart.valid)
      return vertexfit::toFloat(fitDumb);
    return vertexfit::toFloat((fitSmart.chi2 < fitDumb.chi2) ? fitSmart : fitDumb);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::pairfit

#endif  // RecoVertex_Vega_plugins_alpaka_VegaPairFit_h
