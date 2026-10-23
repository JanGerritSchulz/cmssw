#ifndef RecoVertex_Vega_plugins_alpaka_VegaMath_h
#define RecoVertex_Vega_plugins_alpaka_VegaMath_h

#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::math {

  // Minimal 3-vector, used for points and directions alike. Host+device
  // (unlike the acc-templated functions), so tests can use it anywhere.
  template <typename T>
  struct Vec3 {
    T x, y, z;
  };

  // Symmetric 3x3 matrix, stored as the 6 independent entries.
  template <typename T>
  struct SymMat3 {
    T xx, xy, xz, yy, yz, zz;
  };

  using Vec3f = Vec3<float>;
  using Vec3d = Vec3<double>;
  using SymMat3f = SymMat3<float>;
  using SymMat3d = SymMat3<double>;

  // Precision conversion, e.g. cast<float>(Vec3d) at the end of a double fit or
  // cast<double>(SymMat3f) before inverting a stored float covariance.
  template <typename To, typename From>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE Vec3<To> cast(Vec3<From> v) {
    return {To(v.x), To(v.y), To(v.z)};
  }

  template <typename To, typename From>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE SymMat3<To> cast(SymMat3<From> m) {
    return {To(m.xx), To(m.xy), To(m.xz), To(m.yy), To(m.yz), To(m.zz)};
  }

  template <typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE Vec3<T> operator-(Vec3<T> a, Vec3<T> b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  }

  template <typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE T dot(Vec3<T> a, Vec3<T> b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  template <typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE T norm2(Vec3<T> v) {
    return dot(v, v);
  }

  template <typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE SymMat3<T> operator+(SymMat3<T> a, SymMat3<T> b) {
    return {a.xx + b.xx, a.xy + b.xy, a.xz + b.xz, a.yy + b.yy, a.yz + b.yz, a.zz + b.zz};
  }

  template <typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE Vec3<T> mul(SymMat3<T> const& m, Vec3<T> v) {
    return {m.xx * v.x + m.xy * v.y + m.xz * v.z,
            m.xy * v.x + m.yy * v.y + m.yz * v.z,
            m.xz * v.x + m.yz * v.y + m.zz * v.z};
  }

  template <typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE SymMat3<T> outerProduct(Vec3<T> v) {
    return {v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.y, v.y * v.z, v.z * v.z};
  }

  // Clip a value to [-maxAbsValue, maxAbsValue] (used for the arc-length bound
  // |l| <= maxArcLength in both precisions).
  template <typename TAcc, typename T>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE T clamp(TAcc const& acc, T value, T maxAbsValue) {
    return alpaka::math::min(acc, alpaka::math::max(acc, value, -maxAbsValue), maxAbsValue);
  }

  // Smallest accepted normalized determinant of a matrix to be inverted
  // (det of m / trace^3-scaled, see inverseSym). Exactly parallel track pairs
  // end up below it and are reported as invalid instead of garbage.
  inline constexpr double minNormalizedDet = 1e-12;

  // Inverse of a symmetric 3x3 via cofactors, in double, with a RELATIVE
  // pivot guard: m is normalized by its mean diagonal first, so the guard
  // does not depend on units (an absolute guard on cm^2 quantities was
  // clipping every real covariance). Returns false for singular or
  // non-finite input (NaN fails the comparisons); out is then all zeros.
  // Explicitly uses double since float proved to be insufficient.
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE bool inverseSym(TAcc const&, SymMat3d const& m, SymMat3d& out) {
    out = {0, 0, 0, 0, 0, 0};
    const double s = (m.xx + m.yy + m.zz) * (1.0 / 3.0);
    if (!(s > 0.0))
      return false;
    const double is = 1.0 / s;
    const double xx = m.xx * is, xy = m.xy * is, xz = m.xz * is, yy = m.yy * is, yz = m.yz * is, zz = m.zz * is;
    const double cxx = yy * zz - yz * yz, cxy = xz * yz - xy * zz, cxz = xy * yz - xz * yy;
    const double cyy = xx * zz - xz * xz, cyz = xy * xz - xx * yz, czz = xx * yy - xy * xy;
    const double det = xx * cxx + xy * cxy + xz * cxz;
    if (!(det > minNormalizedDet))
      return false;
    const double k = is / det;
    out = {cxx * k, cxy * k, cxz * k, cyy * k, cyz * k, czz * k};
    return true;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega::math

#endif
