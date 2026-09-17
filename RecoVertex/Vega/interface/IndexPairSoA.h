#ifndef RecoVertex_Vega_interface_IndexPairSoA_h
#define RecoVertex_Vega_interface_IndexPairSoA_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace reco {

  // pair of indices
  GENERATE_SOA_LAYOUT(IndexPairSoALayout, SOA_COLUMN(uint32_t, idx1), SOA_COLUMN(uint32_t, idx2))

  using IndexPairSoA = IndexPairSoALayout<>;
  using IndexPairSoAView = IndexPairSoA::View;
  using IndexPairSoAConstView = IndexPairSoA::ConstView;

}  // namespace reco

#endif  // RecoVertex_Vega_interface_IndexPairSoA_h
