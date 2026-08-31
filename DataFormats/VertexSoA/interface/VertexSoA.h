#ifndef DataFormats_VertexSoA_interface_VertexSoA_h
#define DataFormats_VertexSoA_interface_VertexSoA_h

#include <Eigen/Core>
#include <Eigen/Dense>
#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace reco {
  // covariance elements
  using CovarianceMatrix = Eigen::Matrix<float, 6, 1>;
  enum VertexSoACovIdx : size_t { kVarX = 0, kCovXY, kCovXZ, kVarY, kCovYZ, kVarZ };

  GENERATE_SOA_LAYOUT(VertexSoALayout,
                      SOA_COLUMN(float, x),
                      SOA_COLUMN(float, y),
                      SOA_COLUMN(float, z),
                      SOA_COLUMN(float, t),

                      SOA_COLUMN(float, chi2),
                      SOA_COLUMN(float, ndof),
                      SOA_EIGEN_COLUMN(CovarianceMatrix, covariance),

                      SOA_COLUMN(uint8_t, nTracks),
                      SOA_COLUMN(uint32_t, trackOffsets),

                      SOA_SCALAR(uint32_t, nVertices))

  GENERATE_SOA_LAYOUT(VertexTracksSoALayout, SOA_COLUMN(uint32_t, id))

  GENERATE_SOA_BLOCKS(VertexSoABlocksLayout,
                      SOA_BLOCK(vertex, VertexSoALayout),
                      SOA_BLOCK(tracks, VertexTracksSoALayout))

  using VertexSoA = VertexSoALayout<>;
  using VertexSoAView = VertexSoA::View;
  using VertexSoAConstView = VertexSoA::ConstView;

  using VertexTracksSoA = VertexTracksSoALayout<>;
  using VertexTracksSoAView = VertexTracksSoA::View;
  using VertexTracksSoAConstView = VertexTracksSoA::ConstView;

  // SoABlocks Layout that combines vertices and associated tracks
  using VertexSoABlocks = VertexSoABlocksLayout<>;
  using VertexSoABlocksView = VertexSoABlocks::View;
  using VertexSoABlocksConstView = VertexSoABlocks::ConstView;
}  // namespace reco

#endif  // DataFormats_VertexSoA_interface_VertexSoA_h
