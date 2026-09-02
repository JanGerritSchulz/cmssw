#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexSoA/interface/VertexHost.h"

namespace edmtest {

  class TestReadHostVertexSoA : public edm::global::EDAnalyzer<> {
  public:
    TestReadHostVertexSoA(edm::ParameterSet const&);
    void analyze(edm::StreamID, edm::Event const&, edm::EventSetup const&) const override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    edm::EDGetTokenT<reco::VertexHost> getToken_;
  };

  TestReadHostVertexSoA::TestReadHostVertexSoA(edm::ParameterSet const& iPSet)
      : getToken_(consumes(iPSet.getParameter<edm::InputTag>("input"))) {}

  void TestReadHostVertexSoA::analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const&) const {
    auto const& Vertexs = iEvent.get(getToken_);
    auto VertexsView = Vertexs.view();

    std::cout << "VertexsView.vertex().metadata().size() " << VertexsView.vertex().metadata().size() << std::endl;
    std::cout << "VertexsView.vertex()[10].chi2() " << VertexsView.vertex()[10].chi2() << std::endl;
    for (int i = 0; i < VertexsView.vertex().metadata().size(); ++i) {
      if (VertexsView.vertex()[i].chi2() != float(i)) {
        throw cms::Exception("TestReadHostVertexSoA Failure") << "TestReadHostVertexSoA::analyze, entry. i = " << i;
      }
    }
  }

  void TestReadHostVertexSoA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("input");
    descriptions.addDefault(desc);
  }
}  // namespace edmtest

using edmtest::TestReadHostVertexSoA;
DEFINE_FWK_MODULE(TestReadHostVertexSoA);
