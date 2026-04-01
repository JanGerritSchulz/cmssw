#ifndef SimDataFormats_TrackingAnalysis_SimPixelTrack_h
#define SimDataFormats_TrackingAnalysis_SimPixelTrack_h

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitFwd.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

/** @brief Semi-Monte Carlo truth information used for pixel-tracking opimization.
 *
 * SimPixelTracks hold references to all pixel RecHits of a simulated TrackingParticle.
 * Ones those RecHits are sorted according to their position relative to the particle vertex
 * by the method sortRecHits(), you can create the true doublets of RecHits that the 
 * TrackingParticle left in the detector. These SimPixelTrack::Doublet objects can be used  
 * to optimize the doublet creation in the reconstruction.
 *
 * The Doublets are generated as the RecHit pairs between two consecutively hit layers.
 * I.e., if a TrackingParticle produces
 *  - 1 hit (A) in 1st layer
 *  - 2 hits (B, C) in 3rd layer
 *  - 1 hit (D) in 4th layer
 * then, the true Doublets are:
 *  (A-B), (A-C), (B-D) and (C-D).
 * So, neither does it matter that the 2nd layer got "skipped" as there are no hits,
 * nor is the Doublet of (A-D) formed since there is a layer with hits in between.
 * Doublets are not created between hits within the same layer.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date January 2025
 */
class SimPixelTrack {
public:
  // types for SimPixelTrack properties (also used for CA cuts in the SimPixelTracksAnalyzer)
  using float_type = double;
  using int_type = int16_t;
  using layer_type = uint16_t;
  using status_type = uint8_t;

  inline static layer_type kInvalidLayerId = std::numeric_limits<layer_type>::max();

  struct Hits {
    // vectors of usual hit properties
    std::vector<unsigned int> detIds;          // detector Ids of the RecHits
    std::vector<int> moduleIds;                // module Ids of the RecHits
    std::vector<GlobalPoint> globalPositions;  // global positions of the RecHits (corrected by beamspot)
    std::vector<layer_type> layerIds;          // layer IDs corresponding to the RecHits
    std::vector<int_type> clusterYSizes;       // cluster sizes (local y) corresponding to the RecHits

    // vectors of vectorHit properties
    std::vector<float_type> dPhiDrs;     // dPhi / dR of the stub, -1 for non-stub hit
    std::vector<float_type> dPhiDrErrs;  // error of dPhi / dR of the stub, -1 for non-stub hit
  };

  /**
    * Sub-class for true doublets of RecHits
    *  - first hit = inner RecHit
    *  - second hit = outer RecHit
    */
  class Doublet {
  public:
    // possible states of the doublet (could be set by an analyzer according to doublet cuts)
    enum class Status : status_type { undef, alive, killedByCuts, killedByMissingLayerPair };

    struct Triplet {
      Triplet(size_t innerDoubletIndex, size_t nInnerTriplets)
          : index_(innerDoubletIndex), status_(Status::undef), quadrupletIsKilled_(nInnerTriplets, false) {}

      size_t innerDoubletIndex() const { return index_; }

      // methods to set status to undef, alive or killed
      void setUndef() { status_ = Status::undef; }
      void setAlive() { status_ = Status::alive; }
      void setKilled() { status_ = Status::killedByCuts; }

      // methods to check if status is undef, alive or killed
      bool isUndef() const { return status_ == Status::undef; }
      bool isAlive() const { return status_ == Status::alive; }
      bool isKilled() const { return status_ == Status::killedByCuts; }

      void setCurvature(float_type const curvature) { curvature_ = curvature; }
      float_type curvature() const { return curvature_; }

      void setKilledQuadruplet(size_t i) { quadrupletIsKilled_.at(i) = true; }
      bool isKilledQuadruplet(size_t i) { return quadrupletIsKilled_.at(i); }
      std::vector<bool> const& quadruplets() const { return quadrupletIsKilled_; }

      size_t index_;                            // index of the inner doublet of the triplet
      Status status_;                           // status of the triplet
      float_type curvature_{-99999};            // curvature of the triplet
      std::vector<bool> quadrupletIsKilled_{};  // status of the quadruplets with this triplet as the the outer triplet
    };

    // default constructor
    Doublet() = default;

    // constructor
    Doublet(SimPixelTrack const&, size_t const, size_t const, const TrackerTopology*, std::vector<size_t> const&);

    // method to get the layer pair
    std::pair<layer_type, layer_type> layerIds() const { return layerIds_; }

    // method to get the number of skipped layers
    int_type numSkippedLayers() const { return numSkippedLayers_; }

    // method to get the layer pair ID
    layer_type layerPairId() const { return layerPairId_; }

    // methods to get the inner/outer layerId
    layer_type innerLayerId() const { return layerIds_.first; }
    layer_type outerLayerId() const { return layerIds_.second; }

    // methods to get the cluster size of the inner/outer RecHit
    int_type innerClusterYSize() const { return clusterYSizes_.first; }
    int_type outerClusterYSize() const { return clusterYSizes_.second; }

    // methods to get the module ids of the inner/outer RecHit
    unsigned int innerModuleId() const { return moduleIds_.first; }
    unsigned int outerModuleId() const { return moduleIds_.second; }

    // methods to get the global position of the inner/outer RecHit
    GlobalPoint innerGlobalPos() const { return globalPositions_.first; };
    GlobalPoint outerGlobalPos() const { return globalPositions_.second; };

    // methods to set status to undef, alive or killed
    void setUndef() { status_ = Status::undef; }
    void setAlive() { status_ = Status::alive; }
    void setKilledByCuts() { status_ = Status::killedByCuts; }
    void setKilledByMissingLayerPair() { status_ = Status::killedByMissingLayerPair; }
    void setValidStart() { validStart_ = true; }

    // methods to check if status is undef, alive or killed
    bool isUndef() const { return status_ == Status::undef; }
    bool isAlive() const { return status_ == Status::alive; }
    bool isKilledByCuts() const { return status_ == Status::killedByCuts; }
    bool isKilledByMissingLayerPair() const { return status_ == Status::killedByMissingLayerPair; }
    bool isKilled() const { return isKilledByCuts() || isKilledByMissingLayerPair(); }
    bool isValidStart() const { return validStart_; }

    // methods to get the vector of inner triplets
    std::vector<Triplet>& innerTriplets() { return innerTriplets_; }
    std::vector<Triplet> const& innerTripletsView() const { return innerTriplets_; }
    size_t innerNeighborIndex(size_t i) const { return innerTriplets_.at(i).innerDoubletIndex(); }
    // method to get the number of Triplets
    size_t numInnerTriplets() const { return innerTriplets_.size(); }
    // method to get the inner layer ID of the Triplets
    layer_type innerTripletsInnerLayerId() const { return innerTripletsInnerLayerId_; }

  private:
    std::pair<int, int> moduleIds_;                        // module Ids of the RecHits of the Doublet
    std::pair<GlobalPoint, GlobalPoint> globalPositions_;  // global position of the RecHits of the Doublet
                                                           // (corrected by beamspot)
    std::pair<layer_type, layer_type> layerIds_;           // pair of layer IDs corresponding to the RecHits
    std::pair<int_type, int_type> clusterYSizes_;          // pair of cluster sizes corresponding to the RecHits
    Status status_;                                        // status of the doublet
    bool validStart_{false};                               // doublet passes cuts for starting Ntuplets
    int_type numSkippedLayers_;                            // number of layers skipped by the Doublet
    layer_type layerPairId_;                // ID of the layer pair as defined in the reconstruction for the doublets
    std::vector<Triplet> innerTriplets_{};  // indices of inner triplets and its status
    layer_type innerTripletsInnerLayerId_{99};  // layer ID of the inner RecHit of the triplets
  };

  /**
    * Sub-class for true Ntuplets of the Tracking Particle
    * - keep track of length == number of doublets
    * - first and last layer
    * - whether the Ntuplet is actually created (survives all cuts)
    */
  class Ntuplet {
  public:
    // flags indicating qualities of Ntuplet (depending on its constituents)
    // The order is chosen in such a way that a smaller status value means that the Ntuplet get farther
    // in the reconstruction chain. Hence, a value of 0 corresponds to the Ntuplet surviving reconstruction.
    enum class StatusBit : status_type {
      isTooShort = 1,
      hasMissingLayerPair = 1 << 1,
      hasUndefDoubletCuts = 1 << 2,
      hasKilledDoublets = 1 << 3,
      hasUndefTripletCuts = 1 << 4,
      hasKilledTriplets = 1 << 5,
      hasKilledQuadruplets = 1 << 6,
      invalidStart = 1 << 7
    };

    // default constructor
    Ntuplet() = default;

    // constructor
    Ntuplet(size_t const numDoublets,
            status_type const status,
            layer_type const firstLayerId,
            layer_type const secondLayerId,
            layer_type const lastLayerId,
            int_type const numSkippedLayers)
        : numDoublets_(numDoublets),
          status_(status),
          firstLayerId_(firstLayerId),
          secondLayerId_(secondLayerId),
          lastLayerId_(lastLayerId),
          numSkippedLayers_(numSkippedLayers) {};

    // accessing the different members
    size_t numDoublets() const { return numDoublets_; }
    size_t numRecHits() const { return (numDoublets_ + 1); }
    layer_type firstLayerId() const { return firstLayerId_; }
    layer_type secondLayerId() const { return secondLayerId_; }
    layer_type lastLayerId() const { return lastLayerId_; }
    int_type numSkippedLayers() const { return numSkippedLayers_; }

    // method to update an external status
    static status_type updateStatus(status_type status,
                                    bool const hasUndefDoubletCuts,
                                    bool const hasMissingLayerPair,
                                    bool const hasKilledDoublets,
                                    bool const hasUndefTripletCuts,
                                    bool const hasKilledTriplets,
                                    bool const hasKilledQuadruplets,
                                    bool const isTooShort = false,
                                    bool const invalidStart = false) {
      return status | (status_type(hasUndefDoubletCuts) * status_type(StatusBit::hasUndefDoubletCuts) +
                       status_type(hasMissingLayerPair) * status_type(StatusBit::hasMissingLayerPair) +
                       status_type(hasKilledDoublets) * status_type(StatusBit::hasKilledDoublets) +
                       status_type(hasUndefTripletCuts) * status_type(StatusBit::hasUndefTripletCuts) +
                       status_type(hasKilledTriplets) * status_type(StatusBit::hasKilledTriplets) +
                       status_type(hasKilledQuadruplets) * status_type(StatusBit::hasKilledQuadruplets) +
                       status_type(isTooShort) * status_type(StatusBit::isTooShort) +
                       status_type(invalidStart) * status_type(StatusBit::invalidStart));
    }

    // methods to set status to alive, undef or killed
    void setUndefDoubletCuts() { status_ |= status_type(StatusBit::hasUndefDoubletCuts); }
    void setUndefTripletCuts() { status_ |= status_type(StatusBit::hasUndefTripletCuts); }
    void setMissingLayerPair() { status_ |= status_type(StatusBit::hasMissingLayerPair); }
    void setKilledDoublets() { status_ |= status_type(StatusBit::hasKilledDoublets); }
    void setKilledTriplet() { status_ |= status_type(StatusBit::hasKilledTriplets); }
    void setKilledQuadruplet() { status_ |= status_type(StatusBit::hasKilledQuadruplets); }
    void setTooShort() { status_ |= status_type(StatusBit::isTooShort); }
    void setInvalidStart() { status_ |= status_type(StatusBit::invalidStart); }

    // methods to check if status is undef, alive or killed
    bool hasUndefDoubletCuts() const { return status_ & status_type(StatusBit::hasUndefDoubletCuts); }
    bool hasUndefTripletCuts() const { return status_ & status_type(StatusBit::hasUndefTripletCuts); }
    bool hasUndef() const { return hasUndefDoubletCuts() || hasUndefTripletCuts(); }
    bool hasMissingLayerPair() const { return status_ & status_type(StatusBit::hasMissingLayerPair); }
    bool hasKilledDoublets() const { return status_ & status_type(StatusBit::hasKilledDoublets); }
    bool hasKilledTriplets() const { return status_ & status_type(StatusBit::hasKilledTriplets); }
    bool hasKilledQuadruplets() const { return status_ & status_type(StatusBit::hasKilledQuadruplets); }
    bool isKilled() const {
      return hasMissingLayerPair() || hasKilledDoublets() || hasKilledTriplets() || hasKilledQuadruplets();
    }
    bool isTooShort() const { return status_ & status_type(StatusBit::isTooShort); }
    bool invalidStart() const { return status_ & status_type(StatusBit::invalidStart); }
    bool isAlive() const { return !(status_); }  // if nothing is set (no undef and no kills) the tuplet is alive

    // method to get the leading digit of the status (first non-zero one),
    // e.g. status=00110100 -> failingRecoStep()=00000100
    // This represents the first step of the reco chain the given Ntuplet fails.
    // For an alive Ntuplet, return the max value 11111111.
    status_type failingRecoStep() const {
      if (isAlive())
        return 0b11111111;
      else
        return (status_ & ((~status_) + 1));
    }

    // method to compare the own status to a given reference and check which one gets farther in the reconstruction chain
    bool getsFartherInRecoChainThanReference(Ntuplet const& referenceNtuplet) const {
      return failingRecoStep() > referenceNtuplet.failingRecoStep();
    }
    // method to check if the own Ntuplet gets exactly as far in the reco chain as a given reference
    bool getsAsFarInRecoChainAsReference(Ntuplet const& referenceNtuplet) const {
      return failingRecoStep() == referenceNtuplet.failingRecoStep();
    }

  private:
    size_t numDoublets_;  // number of doublets in the Ntuplet
    status_type status_;  // status flags of the Ntuplet (missing layer pairs, undefined cuts, killed doublets, etc.)
    layer_type firstLayerId_;    // index of the first layer of the Ntuplet
    layer_type secondLayerId_;   // index of the second layer of the Ntuplet
    layer_type lastLayerId_;     // index of the last layer of the Ntuplet
    int_type numSkippedLayers_;  // number of skipped layers over the full Ntuplet (sum of skips by doublets)
  };

  // default contructor
  SimPixelTrack() = default;

  // constructor
  SimPixelTrack(TrackingParticleRef const trackingParticleRef, reco::BeamSpot const& beamSpot)
      : trackingParticleRef_(trackingParticleRef), beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}
  SimPixelTrack(reco::TrackBaseRef const trackRef, reco::BeamSpot const& beamSpot)
      : trackRef_(trackRef), beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}
  SimPixelTrack(reco::BeamSpot const& beamSpot) : beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}

  // method to add a RecHit to the SimPixelTrack
  void addRecHit(TrackingRecHit const& recHit,
                 layer_type const layerId,
                 int_type const clusterYSize,
                 unsigned int const detId,
                 int const moduleId);

  // method to get the reference to the TrackingParticle
  TrackingParticleRef trackingParticle() const { return trackingParticleRef_; }
  // method to get the reference to the track
  reco::TrackBaseRef track() const { return trackRef_; }

  // method to get the detector id vector
  std::vector<unsigned int> detIds() const { return hits_.detIds; }
  // method to get the detector id at index i
  unsigned int detIds(size_t const i) const { return hits_.detIds[i]; }

  // method to get the module id vector
  std::vector<int> moduleIds() const { return hits_.moduleIds; }
  // method to get the module id at index i
  int moduleIds(size_t const i) const { return hits_.moduleIds[i]; }

  // method to get the global position vector of the RecHits
  std::vector<GlobalPoint> globalPositions() const { return hits_.globalPositions; }
  // method to get the global position of the RecHit at index i
  GlobalPoint globalPositions(size_t const i) const { return hits_.globalPositions[i]; }

  // method to get the layer id vector
  std::vector<layer_type> layerIds() const { return hits_.layerIds; }
  // method to get the layer id at index i
  layer_type layerIds(size_t const i) const { return hits_.layerIds[i]; }

  // method to get the cluster size vector
  std::vector<int_type> clusterYSizes() const { return hits_.clusterYSizes; }
  // method to get the cluster size at index i
  int_type clusterYSizes(size_t const i) const { return hits_.clusterYSizes[i]; }

  // method to get the beam spot position
  GlobalVector beamSpotPosition() const { return beamSpotPosition_; }

  // method to get the number of layers
  size_t numLayers() const { return numLayers_; }
  // method to get number of RecHits in the SimPixelTrack
  size_t numRecHits() const { return hits_.layerIds.size(); }
  // method to get the number of SimDoublets
  size_t numDoublets() const { return doublets_.size(); }

  // method to sort the RecHits according to the position (either a given reference point or the TP vertex)
  void sortRecHits();
  void sortRecHits(float const, float const, float const);

  // method to produce the SimDoublets from the RecHits
  void buildSimDoublets(const TrackerTopology* trackerTopology) const;
  // method to access the SimDoublets
  std::vector<Doublet>& getSimDoublets() const { return doublets_; }
  // method to build and access the SimDoublets
  std::vector<Doublet>& buildAndGetSimDoublets(const TrackerTopology* trackerTopology) const {
    buildSimDoublets(trackerTopology);
    return doublets_;
  }
  // method to access a single SimDoublet
  Doublet const& getSimDoublet(size_t const index) const { return doublets_.at(index); }

  // method to build the SimNtuplets
  // minNumDoubletsToPass = the number of doublets required for the Ntuplet to not be considered too short
  void buildSimNtuplets(size_t const minNumDoubletsToPass = 0) const;
  // method to access the SimNtuplets
  std::vector<Ntuplet>& getSimNtuplets() const { return ntuplets_; };
  // method to build and access the SimNtuplets in one go
  // minNumDoubletsToPass = the number of doublets required for the Ntuplet to not be considered too short
  std::vector<Ntuplet>& buildAndGetSimNtuplets(size_t const minNumDoubletsToPass = 0) const {
    buildSimNtuplets(minNumDoubletsToPass);
    return ntuplets_;
  };

  // method to check if there are SimNtuplets
  bool hasSimNtuplet() const { return longestNtupletIndex_.has_value(); }
  // method to check if there are alive SimNtuplet
  bool hasAliveSimNtuplet() const { return longestAliveNtupletIndex_.has_value(); }

  // method to access the longest SimNtuplet
  Ntuplet const& longestSimNtuplet() const { return ntuplets_.at(*longestNtupletIndex_); }
  // method to access the longest alive SimNtuplet
  Ntuplet const& longestAliveSimNtuplet() const { return ntuplets_.at(*longestAliveNtupletIndex_); }
  // method to access the best SimNtuplet
  Ntuplet const& bestSimNtuplet() const { return ntuplets_.at(*bestNtupletIndex_); }

  // method to get fishbone alignments
  std::vector<std::pair<layer_type, float_type>> fishboneScores() const;

  // method to clear the mutable vectors once you finished using them
  void clearMutables() const {
    doublets_.clear();
    ntuplets_.clear();
    longestNtupletIndex_.reset();
    longestAliveNtupletIndex_.reset();
    bestNtupletIndex_.reset();
  }

private:
  // function for recursive building of Ntuplets
  void buildSimNtuplets(Doublet const& doublet,
                        std::vector<bool> const& quadruplets,
                        size_t numSimDoublets,
                        layer_type const lastLayerId,
                        status_type const status,
                        int_type const numSkippedLayers,
                        size_t const minNumDoubletsToPass) const;

  // class members
  TrackingParticleRef trackingParticleRef_;  // reference to the TrackingParticle (if SimPixelTrack is based on a TP)
  reco::TrackBaseRef trackRef_;              // referency to the track (if SimPixelTrack is based on a track)
  Hits hits_;                                // RecHits associated to the TP
  GlobalVector beamSpotPosition_;  // global position of the beam spot (needed to correct the global RecHit position)
  bool recHitsAreSorted_{false};   // true if RecHits were sorted
  size_t numLayers_{0};            // number of layers hit by the TrackingParticle

  // non-persistent, mutable members:
  // vector of true doublets
  mutable std::vector<Doublet> doublets_{};
  // vector of true Ntuplets
  mutable std::vector<Ntuplet> ntuplets_{};
  // index of the longest SimNtuplet
  mutable std::optional<size_t> longestNtupletIndex_{-1};
  // index of the longest SimNtuplet that survives
  mutable std::optional<size_t> longestAliveNtupletIndex_{-1};
  // index of the SimNtuplet that gets the farthest in the reco chain
  mutable std::optional<size_t> bestNtupletIndex_{-1};
  // vector of length NrecHits that holds for each RecHit references to all the
  // doublets that have this hit as an outer hit
  mutable std::vector<std::vector<size_t>> innerDoubletsOfRecHit_{};
};

// collection of SimPixelTrack
typedef std::vector<SimPixelTrack> SimPixelTrackCollection;

#endif
