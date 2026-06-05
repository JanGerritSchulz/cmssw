#ifndef SimTracker_VertexAssociation_VertexAssociatorByPositionAndTracks_h
#define SimTracker_VertexAssociation_VertexAssociatorByPositionAndTracks_h

#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociatorBaseImpl.h"

/**
 * This class associates reco::Vertices and TrackingVertices by their
 * position (maximum distance in Z should be smaller than absZ and
 * sigmaZ*zError of VertexType), and (optionally) by the fraction of
 * tracks shared by VertexType and TrackingVertex divided by the
 * number of tracks in VertexType. This fraction is always used as
 * the quality in the association, i.e. multiple associations are
 * sorted by it in descending order.
 */
template <typename VertexCollection>
class VertexAssociatorByPositionAndTracks : public reco::VertexToTrackingVertexAssociatorBaseImpl<VertexCollection> {
public:
  using VertexType = typename VertexCollection::value_type;
  using SimToRecoCollection = reco::VertexToTrackingVertexAssociatorBaseImpl<VertexCollection>::SimToRecoCollection;
  using RecoToSimCollection = reco::VertexToTrackingVertexAssociatorBaseImpl<VertexCollection>::RecoToSimCollection;

  VertexAssociatorByPositionAndTracks(const edm::EDProductGetter *productGetter,
                                      double absZ,
                                      double sigmaZ,
                                      double maxRecoZ,
                                      double absT,
                                      double sigmaT,
                                      double maxRecoT,
                                      double sharedTrackFraction,
                                      const reco::RecoToSimCollection *trackRecoToSimAssociation,
                                      const reco::SimToRecoCollection *trackSimToRecoAssociation,
                                      const std::string &weightMethod);

  VertexAssociatorByPositionAndTracks(const edm::EDProductGetter *productGetter,
                                      double absZ,
                                      double sigmaZ,
                                      double maxRecoZ,
                                      double sharedTrackFraction,
                                      const reco::RecoToSimCollection *trackRecoToSimAssociation,
                                      const reco::SimToRecoCollection *trackSimToRecoAssociation,
                                      const std::string &weightMethod);

  ~VertexAssociatorByPositionAndTracks() override = default;

  /* Associate TrackingVertex to RecoVertex By Hits */
  RecoToSimCollection associateRecoToSim(const edm::Handle<edm::View<VertexType>> &vCH,
                                         const edm::Handle<TrackingVertexCollection> &tVCH) const override;

  SimToRecoCollection associateSimToReco(const edm::Handle<edm::View<VertexType>> &vCH,
                                         const edm::Handle<TrackingVertexCollection> &tVCH) const override;

private:
  // ----- member data
  const edm::EDProductGetter *productGetter_;

  const double absZ_;
  const double sigmaZ_;
  const double maxRecoZ_;
  const double absT_;
  const double sigmaT_;
  const double maxRecoT_;
  const double sharedTrackFraction_;

  const reco::RecoToSimCollection *trackRecoToSimAssociation_;
  const reco::SimToRecoCollection *trackSimToRecoAssociation_;

  bool useWeightPtSum2_;
  bool useWeightDzErr_;
};

#endif
