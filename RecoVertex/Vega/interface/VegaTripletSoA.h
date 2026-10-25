#ifndef RecoVertex_Vega_interface_VegaTripletSoA_h
#define RecoVertex_Vega_interface_VegaTripletSoA_h

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"

namespace vega {

  // pair of tracks (two vertex-compatible tracks)
  GENERATE_SOA_LAYOUT(VegaTripletSoALayout,
                      SOA_COLUMN(uint32_t, idx1),  // index of track 1
                      SOA_COLUMN(uint32_t, idx2),  // index of track 2
                      SOA_COLUMN(uint32_t, idx3),  // index of track 3
                      SOA_COLUMN(float, arcLen1),  // arc length of track 1 at closest approach
                      SOA_COLUMN(float, arcLen2),  // arc length of track 2 at closest approach
                      SOA_COLUMN(float, arcLen3),  // arc length of track 3 at closest approach

                      SOA_COLUMN(float, x),                                  // vertex fit x
                      SOA_COLUMN(float, y),                                  // vertex fit y
                      SOA_COLUMN(float, z),                                  // vertex fit z
                      SOA_COLUMN(float, chi2),                               // chi2 of the vertex fit
                      SOA_EIGEN_COLUMN(reco::CovarianceMatrix, covariance),  // vertex fit covariances

                      SOA_SCALAR(uint32_t, nTriplets))  // number of pairs

  using VegaTripletSoA = VegaTripletSoALayout<>;
  using VegaTripletSoAView = VegaTripletSoA::View;
  using VegaTripletSoAConstView = VegaTripletSoA::ConstView;

}  // namespace vega

#endif  // RecoVertex_Vega_interface_VegaTripletSoA_h
