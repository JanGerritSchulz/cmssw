#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "RecoVertex/Vega/interface/alpaka/IndexPairSoACollection.h"
#include "RecoVertex/Vega/plugins/alpaka/VegaGeometry.h"
#include "RecoVertex/Vega/plugins/alpaka/VegaPairFit.h"
#include "RecoVertex/Vega/plugins/VegaConstants.h"

#include "VegaPairFinder.h"

#define VEGA_PAIRS_DEBUG 0
#define VEGA_PAIRS_WARNINGS 0
#define VEGA_PAIRS_PRINTALL 0

namespace ALPAKA_ACCELERATOR_NAMESPACE::vega {

  // -------------------------------------------------------------
  // TRACK SORTING KERNELS
  // -------------------------------------------------------------

  class KernelSortTracksByDz {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, TrkIdx* sortIdx, TrkIdx const nTracks) const {
      for (auto it : cms::alpakatools::uniform_elements(acc, nTracks)) {
        sortIdx[it] = it;
      }

      auto trkStateRow0 = trks[0].state();
      // FIXME: index 4 is hardcoded for dz, should be replaced with a proper accessor
      float const* dz = trkStateRow0.data() + 4 * trkStateRow0.stride();

      auto& sortWorkSpace = alpaka::declareSharedVar<TrkIdx[::vega::maxTracksForVertexing], __COUNTER__>(acc);
      // sort using only 24 bits (meaning it might happen that very close tracks in dz are ordered incorrectly,
      // but this is not a problem for the pair finding, also rare-ish: ~10 in 1k tracks)
      cms::alpakatools::radixSort<Acc1D, float, 3>(acc, dz, sortIdx, sortWorkSpace, nTracks);
    }
  };

  class KernelPrintSortedTracks {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, TrkSoAConstView trks, TrkIdx* sortIdx, TrkIdx const nTracks) const {
      float dzPrev = -9999.f;
      TrkIdx counterWrongOrder = 0;
      for (TrkIdx i{0}; i < nTracks; ++i) {
        TrkIdx it = sortIdx[i];
        float dz = ::reco::zip(trks, it);
        bool wrongOrder = (dz < dzPrev);
        if (wrongOrder) {
          ++counterWrongOrder;
        }
        printf("vega::KernelPrintSortedTracks: (< than prev)=%d, i=%d, it=%d, dz=%f\n", wrongOrder, i, it, dz);
        dzPrev = dz;
      }
      printf("vega::KernelPrintSortedTracks: Total wrong orders found: %d out of %d\n", counterWrongOrder, nTracks);
    }
  };

  // -------------------------------------------------------------
  // PAIR FINDING KERNELS
  // -------------------------------------------------------------

  template <KernelMode Mode>
  class KernelFindPairCandidates {
  public:
    using PairsArg = std::conditional_t<Mode == KernelMode::Form, PairCandSoAView, NoInput>;
    static constexpr TrkIdx T = 32u;  // tile size

    ALPAKA_FN_ACC void operator()(Acc2D const& acc,
                                  PairParams const& params,
                                  TrkSoAConstView trks,
                                  TrkIdx const* sortIdx,
                                  TrkIdx const nTracks,
                                  uint32_t* nPairs,  // Count writes this, Form writes this as well
                                  PairsArg pairs     // Form writes this
    ) const {
      for (TrkIdx i : cms::alpakatools::uniform_elements_y(acc, nTracks)) {
        const TrkIdx si = sortIdx[i];
        auto itrk = trks[si];

        if (itrk.quality() < params.minTrkQuality)
          continue;

        const float idxy = ::reco::tip(trks, si);
        const float iz = ::reco::zip(trks, si);
        const float iphi = ::reco::phi(trks, si);
        const float ieta = itrk.eta();
        const float icottheta = itrk.state()(3);  // FIXME hardcoded index!
        const float icosphi = alpaka::math::cos(acc, iphi);
        const float isinphi = alpaka::math::sin(acc, iphi);
        const float ix = -idxy * isinphi;
        const float iy = idxy * icosphi;

        for (TrkIdx j : cms::alpakatools::uniform_elements_x(acc, nTracks)) {
          // only upper triangle
          if (i >= j)
            continue;

          // find pairs by applying compatibility cuts
          const TrkIdx sj = sortIdx[j];
          auto jtrk = trks[sj];

          if (jtrk.quality() < params.minTrkQuality)
            continue;

          // dz cut
          const float jz = ::reco::zip(trks, sj);
          if (alpaka::math::abs(acc, iz - jz) > params.maxDZ)
            continue;

          // phi cut
          const float jphi = ::reco::phi(trks, sj);
          if (alpaka::math::abs(acc, iphi - jphi) > params.maxDPhi)
            continue;

          // eta cut
          const float jeta = jtrk.eta();
          if (alpaka::math::abs(acc, ieta - jeta) > params.maxDEta)
            continue;

          // linear distance cut
          const float jdxy = ::reco::tip(trks, sj);
          const float jcottheta = jtrk.state()(3);
          const float jcosphi = alpaka::math::cos(acc, jphi);
          const float jsinphi = alpaka::math::sin(acc, jphi);
          const float jx = -jdxy * jsinphi;
          const float jy = jdxy * jcosphi;

          const float nx = isinphi * jcottheta - jsinphi * icottheta;
          const float ny = icottheta * jcosphi - jcottheta * icosphi;
          const float nz = icosphi * jsinphi - jcosphi * isinphi;
          const float dx = ix - jx;
          const float dy = iy - jy;
          const float dz = iz - jz;

          const float maxLinDistance2 = (params.maxLinDistance * params.maxLinDistance) * (nx * nx + ny * ny + nz * nz);
          const float linearApproxDistance = dx * nx + dy * ny + dz * nz;
          if (linearApproxDistance * linearApproxDistance > maxLinDistance2)
            continue;

          // count or fill the pair
          if constexpr (Mode == KernelMode::Count) {
            alpaka::atomicAdd(acc, nPairs, 1u, alpaka::hierarchy::Blocks{});
          } else {
            // fill the pair
            auto idx = alpaka::atomicAdd(acc, nPairs, 1u, alpaka::hierarchy::Blocks{});
            if (idx >= static_cast<uint32_t>(pairs.metadata().size())) {
#if VEGA_PAIRS_WARNINGS
              printf("vega::KernelFindPairCandidates: Warning!!!! Too many pairs (maxNumOfPairs = %d)!\n",
                     pairs.metadata().size());
#endif
              alpaka::atomicSub(acc, nPairs, 1u, alpaka::hierarchy::Blocks{});
              continue;
            }

            pairs[idx].idx1() = alpaka::math::min(acc, si, sj);
            pairs[idx].idx2() = alpaka::math::max(acc, si, sj);

#if VEGA_PAIRS_PRINTALL
            printf("vega::KernelFindPairCandidates: new track pair %d: (%d, %d)\n",
                   idx,
                   alpaka::math::min(acc, si, sj),
                   alpaka::math::max(acc, si, sj));
#endif
          }
        }
      }
    }
  };

  template <KernelMode Mode>
  class KernelFindPairs {
  public:
    using PairsArg = std::conditional_t<Mode == KernelMode::Form, VegaPairSoAView, NoInput>;
    using MapArg = std::conditional_t<Mode == KernelMode::Form, TrkToPair*, NoInput>;

    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  PairParams const& params,
                                  TrkSoAConstView trks,
                                  TrkExtraSoAConstView trksExtra,
                                  uint32_t const* nPairCands,
                                  PairCandSoAView const pairCands,
                                  uint32_t* nPairs,  // Count writes this, Form writes this as well
                                  PairsArg pairs,    // Form writes this
                                  MapArg trk1ToPair  // Form writes this
    ) const {
      for (auto c : cms::alpakatools::uniform_elements(acc, *nPairCands)) {
        auto cand = pairCands[c];
        const TrkIdx i1 = cand.idx1();
        const TrkIdx i2 = cand.idx2();
        auto tk1 = trks[i1];
        auto tk2 = trks[i2];
        auto te1 = trksExtra[i1];
        auto te2 = trksExtra[i2];
        auto covTrk1 = tk1.covariance();
        auto covTrk2 = tk2.covariance();

        const geometry::HelixGeom t1 = {.cx = te1.cx(),
                                        .cy = te1.cy(),
                                        .r = te1.r(),
                                        .q = te1.q(),
                                        .refAngle = te1.refAngle(),
                                        .refZ = ::reco::zip(trks, i1),
                                        .cotTheta = tk1.state()(3),
                                        .sinPhi = te1.sinPhi(),
                                        .cosPhi = te1.cosPhi(),
                                        .pt = tk1.pt()};
        const geometry::HelixGeom t2 = {.cx = te2.cx(),
                                        .cy = te2.cy(),
                                        .r = te2.r(),
                                        .q = te2.q(),
                                        .refAngle = te2.refAngle(),
                                        .refZ = ::reco::zip(trks, i2),
                                        .cotTheta = tk2.state()(3),
                                        .sinPhi = te2.sinPhi(),
                                        .cosPhi = te2.cosPhi(),
                                        .pt = tk2.pt()};

        const geometry::HelixCovariance cov1 = {.varPhi = covTrk1(0),
                                                .covPhiDxy = covTrk1(1),
                                                .covPhiQOverPt = covTrk1(2),
                                                .varDxy = covTrk1(5),
                                                .covDxyQOverPt = covTrk1(6),
                                                .varQOverPt = covTrk1(9),
                                                .varCotTheta = covTrk1(12),
                                                .covCotThetaDz = covTrk1(13),
                                                .varDz = covTrk1(14)};

        const geometry::HelixCovariance cov2 = {.varPhi = covTrk2(0),
                                                .covPhiDxy = covTrk2(1),
                                                .covPhiQOverPt = covTrk2(2),
                                                .varDxy = covTrk2(5),
                                                .covDxyQOverPt = covTrk2(6),
                                                .varQOverPt = covTrk2(9),
                                                .varCotTheta = covTrk2(12),
                                                .covCotThetaDz = covTrk2(13),
                                                .varDz = covTrk2(14)};

        const vertexfit::VertexFitConfig vertexFitConfig = {
            .maxIter = 1,
            .maxIterProfiling = 6,
            .lmDamping = 1e-3,
            .maxArcLength = 100.,
            .tolPosition = 1e-6,
        };
        const pairfit::PairRefineConfig seedingConfig{
            .maxIter = 10,
            .lmDamping = 1e-3f,
            .lmGrow = 5.f,
            .lmShrink = 3.f,
            .maxArcLength = 100.f,
            .tolArcLength = 1e-4f,
            .lmDampingMax = 1e4f,
        };
        auto fit = pairfit::fitTwoTrackVertexGaussNewtonSeeded(acc, t1, t2, cov1, cov2, vertexFitConfig, seedingConfig);

        // reject failing fits
        if (!fit.valid) {
#if VEGA_PAIRS_WARNINGS
          printf("vega::KernelFindPairs: Warning!!!! Fit of track pair failed (track pair %d)!\n", c);
#endif
          continue;
        }

        // cut on the fit quality
        if (fit.chi2 > params.maxChi2)
          continue;

        // cuts on radial displacement significance
        const auto& pos = fit.position;
        const auto& cov = fit.covariance;
        const float x2 = pos.x * pos.x;
        const float y2 = pos.y * pos.y;
        const float dxySigNum2 = (x2 + y2) * (x2 + y2);
        const float dxySigDenom2 = x2 * cov.xx + 2 * pos.x * pos.y * cov.xy + y2 * cov.yy;
        if ((params.minDXYSig >= 0.f) && (dxySigNum2 < params.minDXYSig * params.minDXYSig * dxySigDenom2))
          continue;
        if ((params.maxDXYSig >= 0.f) && (dxySigNum2 > params.maxDXYSig * params.maxDXYSig * dxySigDenom2))
          continue;

        // count or fill the pair
        if constexpr (Mode == KernelMode::Count) {
          alpaka::atomicAdd(acc, nPairs, 1u, alpaka::hierarchy::Blocks{});
        } else {
          // fill the pair
          auto idx = alpaka::atomicAdd(acc, nPairs, 1u, alpaka::hierarchy::Blocks{});
          if (idx >= static_cast<uint32_t>(pairs.metadata().size())) {
#if VEGA_PAIRS_WARNINGS
            printf("vega::KernelFindPairs: Warning!!!! Too many pairs (maxNumOfPairs = %d)!\n",
                   pairs.metadata().size());
#endif
            alpaka::atomicSub(acc, nPairs, 1u, alpaka::hierarchy::Blocks{});
            continue;
          }
          auto pair = pairs[idx];
          pair.idx1() = i1;
          pair.idx2() = i2;
          pair.arcLen1() = fit.l1;
          pair.arcLen2() = fit.l2;
          pair.x() = pos.x;
          pair.y() = pos.y;
          pair.z() = pos.z;
          pair.chi2() = fit.chi2;
          pair.covariance() = math::toEigenCovariance(fit.covariance);
          pair.used() = false;
          pair.start() = true;

          trk1ToPair->count(acc, i1);
        }
      }
    }
  };

  // -------------------------------------------------------------
  // MAP FILLING KERNELS
  // -------------------------------------------------------------

  class KernelFillTrk1ToPair {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  uint32_t* nPairs,
                                  VegaPairSoAConstView pairs,
                                  TrkToPair* trk1ToPair) const {
      for (auto pairIdx : cms::alpakatools::uniform_elements(acc, *nPairs)) {
#if VEGA_PAIRS_PRINTALL
        printf("vega::KernelFillTrk1ToPair: fill %d (trkIdx) -> %d (pairIdx)\n", pairs[pairIdx].idx1(), pairIdx);
#endif
        trk1ToPair->fill(acc, pairs[pairIdx].idx1(), pairIdx);
      }
    }
  };

  // -------------------------------------------------------------
  // MAIN FUNCTION
  // -------------------------------------------------------------

  PairProducts PairFinder::find(TrkSoAConstView trks, TrkExtraSoAConstView trksExtra) const {
    // Get sorted track indices by dz (z of closest approach)
    auto sortedTrackIndices = cms::alpakatools::make_device_buffer<TrkIdx[]>(queue_, nTracks_);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: created sorted track indices\n");
    alpaka::wait(queue_);
#endif

    const auto workDivSortTracksByDz = cms::alpakatools::make_workdiv<Acc1D>(1u, 256u);
    alpaka::exec<Acc1D>(
        queue_, workDivSortTracksByDz, KernelSortTracksByDz{}, trks, sortedTrackIndices.data(), nTracks_);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: sorted tracks by dz\n");
    alpaka::wait(queue_);
#endif

#if VEGA_PAIRS_PRINTALL
    const auto workDivPrintSortedTracks = cms::alpakatools::make_workdiv<Acc1D>(1u, 1u);
    alpaka::exec<Acc1D>(
        queue_, workDivPrintSortedTracks, KernelPrintSortedTracks{}, trks, sortedTrackIndices.data(), nTracks_);

    printf("PairFinder::find: done sorting and printing\n");
    alpaka::wait(queue_);
#endif

    // Find pairs of tracks compatible with vertexing
    auto nPairCands = cms::alpakatools::make_device_buffer<uint32_t>(queue_);
    alpaka::memset(queue_, nPairCands, 0u);
    auto nPairs = cms::alpakatools::make_device_buffer<uint32_t>(queue_);
    alpaka::memset(queue_, nPairs, 0u);

    constexpr TrkIdx T = KernelFindPairCandidates<KernelMode::Count>::T;  // tile size
    const uint32_t nTiles = cms::alpakatools::divide_up_by(nTracks_ + T - 1u, T);

    using Vec2D = alpaka::Vec<alpaka::DimInt<2u>, uint32_t>;
    Vec2D blocksPerGrid{nTiles, nTiles};  // (row tiles, col tiles) -- Alpaka convention: index 0 = outer/y
    Vec2D threadsPerBlock{T, T};
    const auto workDivFindPairCands = cms::alpakatools::make_workdiv<Acc2D>(blocksPerGrid, threadsPerBlock);

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: made work division\n");
    alpaka::wait(queue_);
#endif

    // Size the pair collection: exact count from the Count kernel (requires a
    // device-to-host copy + sync), or the heuristic upper bound.
    uint32_t pairCandCapacity = nTracks_ * 10u;
    if constexpr (::vega::useExactNPairsOnHost) {
      // const auto workDivFindPairCands = cms::alpakatools::make_workdiv<Acc2D>(1u, 256u);
      alpaka::exec<Acc2D>(queue_,
                          workDivFindPairCands,
                          KernelFindPairCandidates<KernelMode::Count>{},
                          params_,
                          trks,
                          sortedTrackIndices.data(),
                          nTracks_,
                          nPairCands.data(),
                          NoInput{});  // Count mode, no pairs are filled

      auto nPairCandsHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
      alpaka::memcpy(queue_, nPairCandsHost, nPairCands);
      alpaka::wait(queue_);  // host sync: pair count is now valid
      pairCandCapacity = *nPairCandsHost.data();

#if VEGA_PAIRS_DEBUG
      printf("PairFinder::find: counted %d pair candidates\n", *nPairCandsHost.data());
      alpaka::wait(queue_);
#endif
    }

    auto pairCands = IndexPairSoACollection(queue_, pairCandCapacity);

    // reset nPairCands counter for the filling
    alpaka::memset(queue_, nPairCands, 0u);

    alpaka::exec<Acc2D>(queue_,
                        workDivFindPairCands,
                        KernelFindPairCandidates<KernelMode::Form>{},
                        params_,
                        trks,
                        sortedTrackIndices.data(),
                        nTracks_,
                        nPairCands.data(),
                        pairCands.view());

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: formed pair candidates\n");
    alpaka::wait(queue_);
#endif

    const uint32_t blockSize = 128u;
    const uint32_t numberOfBlocks = cms::alpakatools::divide_up_by(pairCandCapacity + blockSize - 1u, blockSize);
    const auto workDivFindPairs = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocks, blockSize);

    // Size the pair collection: exact count from the Count kernel (requires a
    // device-to-host copy + sync), or the heuristic upper bound.
    uint32_t pairCapacity = pairCandCapacity;
    if constexpr (false) {  // currently disabled since expensive, could be ::vega::useExactNPairsOnHost again
      // const auto workDivFindPairCands = cms::alpakatools::make_workdiv<Acc2D>(1u, 256u);
      alpaka::exec<Acc1D>(queue_,
                          workDivFindPairs,
                          KernelFindPairs<KernelMode::Count>{},
                          params_,
                          trks,
                          trksExtra,
                          nPairCands.data(),
                          pairCands.view(),
                          nPairs.data(),
                          NoInput{},
                          NoInput{});  // Count mode, no pairs are filled

      auto nPairHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
      alpaka::memcpy(queue_, nPairHost, nPairs);
      alpaka::wait(queue_);  // host sync: pair count is now valid
      pairCapacity = *nPairHost.data();

#if VEGA_PAIRS_DEBUG
      printf("PairFinder::find: counted %d pairs\n", *nPairHost.data());
      alpaka::wait(queue_);
#endif
    }

    // pairs: collection of track pairs
    // trk1ToPair: map from track to pairs in which this track fills idx1
    auto pairs = VegaPairSoACollection(queue_, pairCapacity);
    auto trk1ToPair = PortableOneToManyAssoc<TrkToPair>(queue_, nTracks_, pairCapacity);

    // reset nPairs counter for the filling
    alpaka::memset(queue_, nPairs, 0u);

    alpaka::exec<Acc1D>(queue_,
                        workDivFindPairs,
                        KernelFindPairs<KernelMode::Form>{},
                        params_,
                        trks,
                        trksExtra,
                        nPairCands.data(),
                        pairCands.view(),
                        nPairs.data(),
                        pairs.view(),
                        trk1ToPair.data());

#if VEGA_PAIRS_DEBUG
    auto nPairHost = cms::alpakatools::make_host_buffer<uint32_t, Platform>();
    alpaka::memcpy(queue_, nPairHost, nPairs);
    alpaka::wait(queue_);  // host sync: pair count is now valid
    pairCapacity = *nPairHost.data();

    printf("PairFinder::find: formed %d pairs\n", *nPairHost.data());
    alpaka::wait(queue_);
#endif

    // publish the actual pair count into the SoA (device-to-device copy)
    auto nPairsView = cms::alpakatools::make_device_view(alpaka::getDev(queue_), pairs.view().nPairs());
    alpaka::memcpy(queue_, nPairsView, nPairs);

    const uint32_t blockSizeFill = 128u;
    const uint32_t numberOfBlocksFill =
        cms::alpakatools::divide_up_by(pairCapacity + blockSizeFill - 1u, blockSizeFill);
    const auto workDivFill = cms::alpakatools::make_workdiv<Acc1D>(numberOfBlocksFill, blockSizeFill);

    // fill the trk1ToPair map
    TrkToPair::template launchFinalize<Acc1D>(trk1ToPair.view(), queue_);
    alpaka::exec<Acc1D>(queue_, workDivFill, KernelFillTrk1ToPair{}, nPairs.data(), pairs.view(), trk1ToPair.data());

#if VEGA_PAIRS_DEBUG
    printf("PairFinder::find: done\n");
    alpaka::wait(queue_);
#endif

    return PairProducts{pairCapacity, std::move(pairs), std::move(trk1ToPair)};
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vega
