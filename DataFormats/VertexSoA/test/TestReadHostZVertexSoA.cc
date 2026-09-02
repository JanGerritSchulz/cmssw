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

#include "DataFormats/VertexSoA/interface/ZVertexHost.h"

namespace edmtest {

  class TestReadHostZVertexSoA : public edm::global::EDAnalyzer<> {
  public:
    TestReadHostZVertexSoA(edm::ParameterSet const&);
    void analyze(edm::StreamID, edm::Event const&, edm::EventSetup const&) const override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    edm::EDGetTokenT<reco::ZVertexHost> getToken_;
  };

  TestReadHostZVertexSoA::TestReadHostZVertexSoA(edm::ParameterSet const& iPSet)
      : getToken_(consumes(iPSet.getParameter<edm::InputTag>("input"))) {}

  void TestReadHostZVertexSoA::analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const&) const {
    auto const& ZVertexs = iEvent.get(getToken_);
    auto ZVertexsView = ZVertexs.view();

    std::cout << "ZVertexsView.zvertex().metadata().size() " << ZVertexsView.zvertex().metadata().size() << std::endl;
    std::cout << "ZVertexsView.zvertex()[10].chi2() " << ZVertexsView.zvertex()[10].chi2() << std::endl;
    for (int i = 0; i < ZVertexsView.zvertex().metadata().size(); ++i) {
      if (ZVertexsView.zvertex()[i].chi2() != float(i)) {
        throw cms::Exception("TestReadHostZVertexSoA Failure") << "TestReadHostZVertexSoA::analyze, entry. i = " << i;
      }
    }
  }

  void TestReadHostZVertexSoA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("input");
    descriptions.addDefault(desc);
  }
}  // namespace edmtest

using edmtest::TestReadHostZVertexSoA;
DEFINE_FWK_MODULE(TestReadHostZVertexSoA);
