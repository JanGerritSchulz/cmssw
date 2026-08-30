#ifndef RecoVertex_AdaptiveVertexFitter_interface_VertexFitResult_h
#define RecoVertex_AdaptiveVertexFitter_interface_VertexFitResult_h

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace vertexfit {

  using Vector3d = Eigen::Vector3d;
  using Matrix3d = Eigen::Matrix3d;

  struct VertexFitResult {
    // position vector: (vx, vy, vz)
    Vector3d position;
    /* covariance matrix:
       | cov(x,x) | cov(y,x) | cov(z,x) |      
       | cov(x,y) | cov(y,y) | cov(z,y) |
       | cov(x,z) | cov(y,z) | cov(z,z) |
    */
    Matrix3d covariances;
    float chi2;
    int ndof;
  };
}  // namespace vertexFit

#endif
