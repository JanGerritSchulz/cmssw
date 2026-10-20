#ifndef RecoVertex_Vega_interface_VegaStructures_h
#define RecoVertex_Vega_interface_VegaStructures_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackSoA/interface/TracksSoA.h"
#include "DataFormats/VertexSoA/interface/VertexSoA.h"
#include "DataFormats/VertexSoA/interface/alpaka/VertexSoACollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/OneToManyAssoc.h"
#include "HeterogeneousCore/AlpakaInterface/interface/SimpleVector.h"
#include "RecoVertex/Vega/interface/IndexPairSoA.h"
#include "RecoVertex/Vega/interface/TrackExtraSoA.h"
#include "RecoVertex/Vega/interface/VegaPairSoA.h"
#include "RecoVertex/Vega/plugins/VegaParams.h"

namespace vega::structs {

  // enum for setting the kernel type/mode:
  // Count: only count the number of pairs/triplets/etc, do not fill them
  //        (used for allocating the exact memory needed)
  // Form: fill the pairs/triplets/etc into the allocated memory
  enum class KernelMode { Count, Form };

  // empty tag type for counting mode kernel
  struct NoInput {};

  // Indices for tracks and pairs
  using TrkIdx = uint16_t;
  using PairIdx = uint16_t;

  using VtxSoAView = ::reco::VertexSoAView;
  using VtxTrkSoAView = ::reco::VertexTracksSoAView;
  using TrkSoAConstView = ::reco::TrackSoAConstView;
  using TrkExtraSoAView = ::vega::TrackExtraSoAView;
  using TrkExtraSoAConstView = ::vega::TrackExtraSoAConstView;

  // track pairs
  using PairCandSoAView = ::vega::IndexPairSoAView;
  using PairCandSoAConstView = ::vega::IndexPairSoAConstView;
  using VegaPairSoAView = ::vega::VegaPairSoAView;
  using VegaPairSoAConstView = ::vega::VegaPairSoAConstView;

  // Vega parameters
  using Params = ::vega::VegaParams;
  using PairParams = ::vega::PairParams;
  using TripletParams = ::vega::TripletParams;
  using VertexParams = ::vega::VertexParams;

  // generic one-to-many association map with random access/filling
  using OneToMany = cms::alpakatools::
      OneToManyAssocRandomAccess<TrkIdx, cms::alpakatools::kDynamicSize, cms::alpakatools::kDynamicSize>;

  // generic one-to-many association map with sequential access/filling
  using OneToManySeq =
      cms::alpakatools::OneToManyAssocSequential<TrkIdx, cms::alpakatools::kDynamicSize, cms::alpakatools::kDynamicSize>;

  // specific association maps for vega
  using PairToPair = OneToMany;
  using TrackToPair = OneToMany;
  using VtxToTrack = OneToManySeq;

}  // namespace vega::structs

#endif  // RecoVertex_Vega_interface_VegaStructures_h
