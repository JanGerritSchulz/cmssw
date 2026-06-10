#ifndef SimTracker_VertexAssociation_VertexAssociatorByPositionAndTracks_h
#define SimTracker_VertexAssociation_VertexAssociatorByPositionAndTracks_h

#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociatorBaseImpl.h"

/**
 * This class associates reco vertex collections and TrackingVertices by their
 * position (maximum distance in Z should be smaller than absZ and
 * sigmaZ*zError of VertexType), and (optionally) by the fraction of
 * tracks shared by VertexType and TrackingVertex divided by the
 * number of tracks in VertexType. This fraction is always used as
 * the quality in the association, i.e. multiple associations are
 * sorted by it in descending order. The fraction supports different
 * weighting methods, selectable via the string argument `weightMethod`:
 *   "none" -> Use the raw fraction:
 *            fraction = nSharedTracks / nTracksInRecoVertex
 *   "pt2" -> Use the pt^2-weighted fraction:
 *            fraction = SumPt2(sharedTracks) / SumPt2(tracksInRecoVertex)
 *   "dzError" -> Use the 1/dzError^2-weighted fraction:
 *            fraction = Sum1OverDzErr2(sharedTracks) / Sum1OverDzErr2(tracksInRecoVertex)
 *   "nSharedTracks" -> Use the raw number of shared tracks:
 *            fraction = nSharedTracks
 *
 * Supported vertex collection types:
 *   std::vector<reco::Vertex>                      (track-based PVs and SVs)
 *   std::vector<reco::VertexCompositePtrCandidate>  (PF-candidate-based SVs)
 *
 * The SimVertex filter applied in the PV case (keeping only the first
 * TrackingVertex per event at BX=0) is controlled via the filterSimVerticesForPVs
 * flag. It should be enabled for PV association and disabled for SV
 * association, where all TrackingVertices are candidates for matching.
 */
template <typename VertexCollection>
class VertexAssociatorByPositionAndTracks : public reco::VertexToTrackingVertexAssociatorBaseImpl<VertexCollection> {
public:
  using VertexType = typename VertexCollection::value_type;
  using SimToRecoCollection =
      typename reco::VertexToTrackingVertexAssociatorBaseImpl<VertexCollection>::SimToRecoCollection;
  using RecoToSimCollection =
      typename reco::VertexToTrackingVertexAssociatorBaseImpl<VertexCollection>::RecoToSimCollection;

  static constexpr double kCheckDisabled = std::numeric_limits<double>::max();

  /// Full constructor including timing parameters.
  VertexAssociatorByPositionAndTracks(const edm::EDProductGetter *productGetter,
                                      double sigmaX,
                                      double sigmaY,
                                      double sigmaZ,
                                      double absZ,
                                      double maxRecoZ,
                                      double sigmaT,
                                      double absT,
                                      double maxRecoT,
                                      double sharedTrackFraction,
                                      const reco::RecoToSimCollection *trackRecoToSimAssociation,
                                      const reco::SimToRecoCollection *trackSimToRecoAssociation,
                                      const std::string &weightMethod,
                                      bool filterSimVerticesForPVs = true);

  /// Constructor without timing parameters (timing disabled).
  VertexAssociatorByPositionAndTracks(const edm::EDProductGetter *productGetter,
                                      double sigmaX,
                                      double sigmaY,
                                      double sigmaZ,
                                      double absZ,
                                      double maxRecoZ,
                                      double sharedTrackFraction,
                                      const reco::RecoToSimCollection *trackRecoToSimAssociation,
                                      const reco::SimToRecoCollection *trackSimToRecoAssociation,
                                      const std::string &weightMethod,
                                      bool filterSimVerticesForPVs = true);

  ~VertexAssociatorByPositionAndTracks() override = default;

  RecoToSimCollection associateRecoToSim(const edm::Handle<edm::View<VertexType>> &vCH,
                                         const edm::Handle<TrackingVertexCollection> &tVCH) const override;

  SimToRecoCollection associateSimToReco(const edm::Handle<edm::View<VertexType>> &vCH,
                                         const edm::Handle<TrackingVertexCollection> &tVCH) const override;

private:
  // Returns kCheckDisabled for negative values and the original value otherwise.
  double getValueIfEnable(const double value) const { return value < 0 ? kCheckDisabled : value; }

  // Returns true if the reco vertex should be skipped entirely.
  // Specialised per VertexType in the .cc file.
  bool isRecoVertexInvalid(const VertexType &vertex) const;

  // Returns the X,Y,Z position of the reco vertex.
  // Specialised per VertexType in the .cc file.
  double recoVertexX(const VertexType &vertex) const;
  double recoVertexY(const VertexType &vertex) const;
  double recoVertexZ(const VertexType &vertex) const;

  // Returns the X,Y,Z error of the reco vertex, used for maxSigmaZ cut.
  // Specialised per VertexType in the .cc file.
  double recoVertexXError(const VertexType &vertex) const;
  double recoVertexYError(const VertexType &vertex) const;
  double recoVertexZError(const VertexType &vertex) const;

  // Returns the T (time) of the reco vertex.
  // Specialised per VertexType in the .cc file.
  double recoVertexT(const VertexType &vertex) const;

  // Returns the T error of the reco vertex.
  // Specialised per VertexType in the .cc file.
  double recoVertexTError(const VertexType &vertex) const;

  // Computes the shared-track fraction between a reco vertex and a sim vertex.
  // Specialised per VertexType in the .cc file to handle the different
  // track-access patterns of reco::Vertex and reco::VertexCompositePtrCandidate.
  float sharedTrackFractionForVertex(const VertexType &recoVertex, const TrackingVertex &simVertex) const;

  // Returns different ref types depending on the reco vertex collection which is required by the association map.
  auto makeVertexRef(const edm::Handle<edm::View<VertexType>> &handle, size_t index) const;

  // ----- member data -----
  const edm::EDProductGetter *productGetter_;

  const double sigmaX_;
  const double sigmaY_;
  const double sigmaZ_;
  const double absZ_;
  const double maxRecoZ_;
  const double sigmaT_;
  const double absT_;
  const double maxRecoT_;
  const double sharedTrackFraction_;

  const reco::RecoToSimCollection *trackRecoToSimAssociation_;
  const reco::SimToRecoCollection *trackSimToRecoAssociation_;

  bool useWeightPtSum2_;
  bool useWeightDzErr_;
  bool useNSharedTracks_;

  // When true, only the first TrackingVertex per event at BX=0 is considered
  // as a sim vertex candidate. Should be true for PV association, false for SV.
  bool filterSimVerticesForPVs_;
};

#endif  // SimTracker_VertexAssociation_VertexAssociatorByPositionAndTracks_h
