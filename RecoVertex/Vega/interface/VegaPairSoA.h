#ifndef RecoVertex_Vega_interface_VegaPairsSoA_h
#define RecoVertex_Vega_interface_VegaPairsSoA_h

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"

namespace vega {

  // pair of tracks (two vertex-compatible tracks)
  GENERATE_SOA_LAYOUT(VegaPairSoALayout,
                      SOA_COLUMN(uint32_t, idx1),  // index of track 1
                      SOA_COLUMN(uint32_t, idx2),  // index of track 2
                      SOA_COLUMN(float, arcLen1),  // arc length of track 1 at closest approach
                      SOA_COLUMN(float, arcLen2),  // arc length of track 2 at closest approach

                      SOA_COLUMN(float, x),                                  // vertex fit x
                      SOA_COLUMN(float, y),                                  // vertex fit y
                      SOA_COLUMN(float, z),                                  // vertex fit z
                      SOA_COLUMN(float, chi2),                               // chi2 of the vertex fit
                      SOA_EIGEN_COLUMN(reco::CovarianceMatrix, covariance),  // vertex fit covariances

                      // pair can be used as starting point for candidate building
                      // because it is not used as a second pair in a triplet somewhere
                      SOA_COLUMN(bool, start),
                      SOA_COLUMN(bool, used),  // used in larger object (flag to mark isolated pairs)

                      SOA_SCALAR(uint32_t, nPairs))  // number of pairs

  using VegaPairSoA = VegaPairSoALayout<>;
  using VegaPairSoAView = VegaPairSoA::View;
  using VegaPairSoAConstView = VegaPairSoA::ConstView;

}  // namespace vega

#endif  // RecoVertex_Vega_interface_VegaPairsSoA_h
