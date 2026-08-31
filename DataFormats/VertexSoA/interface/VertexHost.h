#ifndef DataFormats_VertexSoA_VertexHost_H
#define DataFormats_VertexSoA_VertexHost_H

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace reco {
  using VertexHost = PortableHostCollection<reco::VertexSoABlocks>;
}  // namespace reco

#endif  // DataFormats_VertexSoA_VertexHost_H
