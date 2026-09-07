#ifndef RecoVertex_Vega_interface_VegaStructures_h
#define RecoVertex_Vega_interface_VegaStructures_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/OneToManyAssoc.h"
#include "HeterogeneousCore/AlpakaInterface/interface/SimpleVector.h"
#include "RecoVertex/Vega/interface/TrackExtraSoA.h"
#include "RecoVertex/Vega/plugins/VegaParams.h"

namespace vega::structs {

  // Indices for tracks and pairs
  using TrkIndex = uint16_t;
  using PairIndex = uint16_t;

  using VtxSoAView = ::reco::VertexSoAView;
  using VtxTrkSoAView = ::reco::VertexTracksSoAView;
  using TrkSoAConstView = ::reco::TrackSoAConstView;
  using TrkExtraSoAView = ::vega::TrackExtraSoAView;
  using Params = ::vega::VegaParams;
  using PairParams = ::vega::PairParams;
  using TripletParams = ::vega::TripletParams;

  // generic one-to-many association map with random access/filling
  using OneToMany = cms::alpakatools::
      OneToManyAssocRandomAccess<TrkIndex, cms::alpakatools::kDynamicSize, cms::alpakatools::kDynamicSize>;
  using OneToManyStorage = typename OneToMany::value_type;
  using OneToManyOffsets = typename OneToMany::Counter;
  using OneToManyView = typename OneToMany::View;

  // generic one-to-many association map with sequential access/filling
  using OneToManySeq = cms::alpakatools::
      OneToManyAssocSequential<TrkIndex, cms::alpakatools::kDynamicSize, cms::alpakatools::kDynamicSize>;
  using OneToManySeqStorage = typename OneToManySeq::value_type;
  using OneToManySeqOffsets = typename OneToManySeq::Counter;
  using OneToManySeqView = typename OneToManySeq::View;

  // specific association maps for vega
  using PairToPair = vega::structs::OneToMany;
  using TrackToPair = vega::structs::OneToMany;

}  // namespace vega::structs

#endif  // RecoVertex_Vega_interface_VegaStructures_h
