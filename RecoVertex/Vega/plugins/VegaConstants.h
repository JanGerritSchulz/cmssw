#ifndef RecoVertex_Vega_plugins_VegaConstants_h
#define RecoVertex_Vega_plugins_VegaConstants_h

#include <cstdint>

namespace vega {

  // Upper bound on tracks considered per event for vertexing (pairing/tiling/sorting).
  // Fixed at compile time: sizes static shared-memory buffers in the sorting and
  // pair-finding kernels. Unrelated to TrackerTraits::maxNumberOfTuples, which bounds
  // an earlier, unrelated reconstruction stage.
  constexpr int32_t maxTracksForVertexing = 4096;

  // If true, PairFinder copies the device-computed track count (TrackSoA's
  // nTracks) back to host and uses that exact value to size buffers and work
  // division for the vertexing stage. This forces a device-to-host
  // synchronization at that point in the pipeline, but avoids allocating and
  // launching kernels against the full maxTracksForVertexing capacity when
  // the real per-event track count is much smaller.
  //
  // If false, all vertexing buffers/kernels are sized against the fixed
  // maxTracksForVertexing upper bound instead, and the real nTracks is only
  // read where needed inside device kernels. This keeps the pipeline fully
  // asynchronous (no host readback), at the cost of some wasted allocation
  // and work-division headroom on typical, lower-multiplicity events.
  //
  // FIXME: check which option is to be preferred (timing, memory footprint, ...)
  constexpr bool useExactNTracksOnHost = true;

}  // namespace vega

#endif
