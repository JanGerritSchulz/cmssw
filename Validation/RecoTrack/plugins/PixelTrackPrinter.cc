#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/DynArray.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <cstddef>

namespace simdoublets {
  // function that determines the layerId from the detId for Phase 2
  unsigned int getLayerId(DetId const& detId, const TrackerTopology* trackerTopology) {
    // number of barrel layers
    constexpr unsigned int numBarrelLayers{4};
    // number of disks per endcap
    constexpr unsigned int numEndcapDisks = 12;

    // set default to 999 (invalid)
    unsigned int layerId{99};

    if (detId.subdetId() == PixelSubdetector::PixelBarrel) {
      // subtract 1 in the barrel to get, e.g. for Phase 2, from (1,4) to (0,3)
      layerId = trackerTopology->pxbLayer(detId) - 1;
    } else if (detId.subdetId() == PixelSubdetector::PixelEndcap) {
      if (trackerTopology->pxfSide(detId) == 1) {
        // add offset in the backward endcap to get, e.g. for Phase 2, from (1,12) to (16,27)
        layerId = trackerTopology->pxfDisk(detId) + numBarrelLayers + numEndcapDisks - 1;
      } else {
        // add offest in the forward endcap to get, e.g. for Phase 2, from (1,12) to (4,15)
        layerId = trackerTopology->pxfDisk(detId) + numBarrelLayers - 1;
      }
    }
    // return the determined Id
    return layerId;
  }
}  // namespace simdoublets

// -------------------------------------------------------------------------------------------------------------
// class declaration
// -------------------------------------------------------------------------------------------------------------

class PixelTrackPrinter : public DQMEDAnalyzer {
public:
  explicit PixelTrackPrinter(const edm::ParameterSet&);
  ~PixelTrackPrinter() override;

  void dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) override {};
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  // ------------ member data ------------

  const TrackerTopology* trackerTopology_ = nullptr;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topology_getToken_;
  const edm::EDGetTokenT<edm::View<reco::Track>> tracks_getToken_;
};

void PixelTrackPrinter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("hltPhase2PixelTracks"));

  descriptions.addWithDefaultLabel(desc);
}

PixelTrackPrinter::PixelTrackPrinter(const edm::ParameterSet& pSet)
    : topology_getToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      tracks_getToken_(consumes<edm::View<reco::Track>>(pSet.getParameter<edm::InputTag>("trackSrc"))) {}

PixelTrackPrinter::~PixelTrackPrinter() {}

void PixelTrackPrinter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::cout << "\n================================================" << std::endl;
  std::cout << " Event " << iEvent.id().event() << std::endl;

  // get tracker topology
  trackerTopology_ = &iSetup.getData(topology_getToken_);
  // get tracks
  edm::Handle<View<reco::Track>> trackCollectionHandle;
  if (!iEvent.getByToken(tracks_getToken_, trackCollectionHandle)) {
    return;
  }
  const edm::View<reco::Track>& trackCollection = *trackCollectionHandle;

  // loop over tracks and print information
  for (View<reco::Track>::size_type i = 0; i < trackCollection.size(); ++i) {
    auto track = trackCollection.refAt(i);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << " Track " << i << std::endl;

    std::cout << "   - RecHit layers: ";
    for (size_t j{0}; j < track->recHitsSize(); j++) {
      std::cout << simdoublets::getLayerId(track->recHit(j)->geographicalId(), trackerTopology_) << ", ";
    }
    std::cout << "\n";

    std::cout << "   - RecHit z: ";
    for (size_t j{0}; j < track->recHitsSize(); j++) {
      std::cout << track->recHit(j)->globalPosition().z() << ", ";
    }
    std::cout << "\n";

    std::cout << "   - RecHit r: ";
    for (size_t j{0}; j < track->recHitsSize(); j++) {
      std::cout << track->recHit(j)->globalPosition().perp() << ", ";
    }
    std::cout << "\n";
  }
}

DEFINE_FWK_MODULE(PixelTrackPrinter);
