#include <limits>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/Common/interface/Handle.h"

#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "SimTracker/VertexAssociation/interface/VertexAssociatorByPositionAndTracks.h"

#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"

template <typename VertexCollection>
class VertexAssociatorByPositionAndTracksProducerBase : public edm::global::EDProducer<> {
public:
  explicit VertexAssociatorByPositionAndTracksProducerBase(const edm::ParameterSet &);
  ~VertexAssociatorByPositionAndTracksProducerBase() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  // ----------member data ---------------------------
  const double sigmaX_;
  const double sigmaY_;
  const double sigmaZ_;
  const double absZ_;
  const double maxRecoZ_;
  const double sigmaT_;
  const double absT_;
  const double maxRecoT_;
  const double sharedTrackFraction_;
  const bool filterSimVerticesForPVs_;

  edm::EDGetTokenT<reco::RecoToSimCollection> trackRecoToSimAssociationToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> trackSimToRecoAssociationToken_;
  const std::string weightMethod_;
};

template <typename VertexCollection>
VertexAssociatorByPositionAndTracksProducerBase<VertexCollection>::VertexAssociatorByPositionAndTracksProducerBase(
    const edm::ParameterSet &config)
    : sigmaX_(config.getParameter<double>("sigmaX")),
      sigmaY_(config.getParameter<double>("sigmaY")),
      sigmaZ_(config.getParameter<double>("sigmaZ")),
      absZ_(config.getParameter<double>("absZ")),
      maxRecoZ_(config.getParameter<double>("maxRecoZ")),
      sigmaT_(config.getParameter<double>("sigmaT")),
      absT_(config.getParameter<double>("absT")),
      maxRecoT_(config.getParameter<double>("maxRecoT")),
      sharedTrackFraction_(config.getParameter<double>("sharedTrackFraction")),
      filterSimVerticesForPVs_(config.getParameter<bool>("filterSimVerticesForPVs")),
      trackRecoToSimAssociationToken_(
          consumes<reco::RecoToSimCollection>(config.getParameter<edm::InputTag>("trackAssociation"))),
      trackSimToRecoAssociationToken_(
          consumes<reco::SimToRecoCollection>(config.getParameter<edm::InputTag>("trackAssociation"))),
      weightMethod_(config.getParameter<std::string>("weightMethod")) {
  produces<reco::VertexToTrackingVertexAssociator<VertexCollection>>();
}

template <typename VertexCollection>
VertexAssociatorByPositionAndTracksProducerBase<VertexCollection>::~VertexAssociatorByPositionAndTracksProducerBase() {}

template <typename VertexCollection>
void VertexAssociatorByPositionAndTracksProducerBase<VertexCollection>::fillDescriptions(
    edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;

  // Matching conditions
  desc.add<double>("sigmaX", -1.0);
  desc.add<double>("sigmaY", -1.0);
  desc.add<double>("sigmaZ", 3.0);
  desc.add<double>("absZ", 0.1);
  desc.add<double>("maxRecoZ", 1000.0);
  desc.add<double>("sigmaT", -1.0);
  desc.add<double>("absT", -1.0);
  desc.add<double>("maxRecoT", -1.0);
  desc.add<double>("sharedTrackFraction", -1.0);
  desc.add<std::string>("weightMethod", "none");
  desc.add<bool>("filterSimVerticesForPVs", true);

  // Track-TrackingParticle association
  desc.add<edm::InputTag>("trackAssociation", edm::InputTag("trackingParticleRecoTrackAsssociation"));

  descriptions.addWithDefaultLabel(desc);
}

template <typename VertexCollection>
void VertexAssociatorByPositionAndTracksProducerBase<VertexCollection>::produce(edm::StreamID,
                                                                            edm::Event &iEvent,
                                                                            const edm::EventSetup &) const {
  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  iEvent.getByToken(trackRecoToSimAssociationToken_, recotosimCollectionH);

  edm::Handle<reco::SimToRecoCollection> simtorecoCollectionH;
  iEvent.getByToken(trackSimToRecoAssociationToken_, simtorecoCollectionH);

  std::unique_ptr<VertexAssociatorByPositionAndTracks<VertexCollection>> impl;

  edm::LogWarning("VertexAssociatorByPositionAndTracksProducer") << "XXXXXXX VertexAssociatorByPositionAndTracksProducerBase we are producing! XXXXXXX";

  if (!recotosimCollectionH.isValid() || !simtorecoCollectionH.isValid()) {
    if (!recotosimCollectionH.isValid())
      edm::LogWarning("VertexAssociatorByPositionAndTracksProducer") << "trackRecoToSimAssociation is not available in the event";
    if (!simtorecoCollectionH.isValid())
      edm::LogWarning("VertexAssociatorByPositionAndTracksProducer") << "trackSimToRecoAssociation is not available in the event";
    return;
  }
  if (sigmaT_ < 0.0) {
    impl = std::make_unique<VertexAssociatorByPositionAndTracks<VertexCollection>>(&(iEvent.productGetter()),
                                                                                   sigmaX_,
                                                                                   sigmaY_,
                                                                                   sigmaZ_,
                                                                                   absZ_,
                                                                                   maxRecoZ_,
                                                                                   sharedTrackFraction_,
                                                                                   recotosimCollectionH.product(),
                                                                                   simtorecoCollectionH.product(),
                                                                                   weightMethod_,
                                                                                   filterSimVerticesForPVs_);
  } else {
    impl = std::make_unique<VertexAssociatorByPositionAndTracks<VertexCollection>>(&(iEvent.productGetter()),
                                                                                   sigmaX_,
                                                                                   sigmaY_,
                                                                                   sigmaZ_,
                                                                                   absZ_,
                                                                                   maxRecoZ_,
                                                                                   sigmaT_,
                                                                                   absT_,
                                                                                   maxRecoT_,
                                                                                   sharedTrackFraction_,
                                                                                   recotosimCollectionH.product(),
                                                                                   simtorecoCollectionH.product(),
                                                                                   weightMethod_,
                                                                                   filterSimVerticesForPVs_);
  }

  auto toPut = std::make_unique<reco::VertexToTrackingVertexAssociator<VertexCollection>>(std::move(impl));
  iEvent.put(std::move(toPut));
}

using VertexAssociatorByPositionAndTracksProducer = VertexAssociatorByPositionAndTracksProducerBase<std::vector<reco::Vertex>>;
using VertexAssociatorByPositionAndTracksProducerCPC = VertexAssociatorByPositionAndTracksProducerBase<std::vector<reco::VertexCompositePtrCandidate>>;
DEFINE_FWK_MODULE(VertexAssociatorByPositionAndTracksProducer);
DEFINE_FWK_MODULE(VertexAssociatorByPositionAndTracksProducerCPC);
