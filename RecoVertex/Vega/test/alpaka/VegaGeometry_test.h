#ifndef RecoVertex_Vega_test_alpaka_VegaGeometry_test_h
#define RecoVertex_Vega_test_alpaka_VegaGeometry_test_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace vega::geometry::test {
    void runKernels(Queue& queue);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoVertex_Vega_test_alpaka_VegaGeometry_test_h
