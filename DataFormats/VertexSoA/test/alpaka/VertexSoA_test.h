#ifndef DataFormats_VertexSoA_test_alpaka_VertexSoA_test_h
#define DataFormats_VertexSoA_test_alpaka_VertexSoA_test_h

#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT {

  void runKernels(::reco::VertexSoABlocksView view, Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::testVertexSoAT

#endif  // DataFormats_VertexSoA_test_alpaka_VertexSoA_test_h
