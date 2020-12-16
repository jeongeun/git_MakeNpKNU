#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "TTree.h"
#include "TClonesArray.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TH1F.h>
#include <TDirectoryFile.h>

class MakeNpKNUHist : public edm::EDAnalyzer {
   public:
      explicit MakeNpKNUHist(const edm::ParameterSet&);
      ~MakeNpKNUHist();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // Pileup
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupToken;   bool DoPileup ;
      void ProdPileup(edm::Event& iEvent, const edm::EventSetup& iSetup);

		TH1F*	h1_NumberOfEvents         ;
		TH1F*	h1_TrueNumInteractions    ;
		TH1F*	h1_TrueNumInteractions60  ;
		TH1F*	h1_TrueNumInteractions70  ;
		TH1F*	h1_TrueNumInteractions80  ;
		TH1F*	h1_TrueNumInteractions90  ;
		TH1F*	h1_TrueNumInteractions100 ;
		TH1F*	h1_PU_NumInteractions     ;

		//int debugMax;
		int LoopEvent;
		//bool debug ; 

};
MakeNpKNUHist::MakeNpKNUHist(const edm::ParameterSet& iConfig) {
	//debugMax = iConfig.getParameter<int> ("debugMax"); 

   pileupToken = (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup")))   ;
//   vertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex")) ;

	edm::Service<TFileService> fs;
	h1_NumberOfEvents = fs->make<TH1F>("h1_NumberOfEvents","h1_NumberOfEvents",2,0,2);
	TFileDirectory pileupDir = fs->mkdir("pileupHist");
	h1_TrueNumInteractions = pileupDir.make<TH1F>("h1_TrueNumInteractions","h1_TrueNumInteractions", 50,0,50);
	h1_TrueNumInteractions60 = pileupDir.make<TH1F>("h1_TrueNumInteractions60","h1_TrueNumInteractions60", 60,0,60);
	h1_TrueNumInteractions70 = pileupDir.make<TH1F>("h1_TrueNumInteractions70","h1_TrueNumInteractions70", 70,0,70);
	h1_TrueNumInteractions80 = pileupDir.make<TH1F>("h1_TrueNumInteractions80","h1_TrueNumInteractions80", 80,0,80);
	h1_TrueNumInteractions90 = pileupDir.make<TH1F>("h1_TrueNumInteractions90","h1_TrueNumInteractions90", 90,0,90);
	h1_TrueNumInteractions100 = pileupDir.make<TH1F>("h1_TrueNumInteractions100","h1_TrueNumInteractions100", 100,0,100);
	h1_PU_NumInteractions  = pileupDir.make<TH1F>("h1_PU_NumInteractions" ,"h1_PU_NumInteractions", 100,0,100);

	LoopEvent = 0;
}
MakeNpKNUHist::~MakeNpKNUHist() {
}

void MakeNpKNUHist::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
	using std::cout;
	using std::endl;
	
	h1_NumberOfEvents->Fill(1);

	if(iEvent.isRealData()) return;

	LoopEvent++;
	//debug = (debugMax<0 or LoopEvent <= debugMax);

	edm::Handle<std::vector<PileupSummaryInfo> > PileupInfos;
	iEvent.getByToken(pileupToken, PileupInfos);
	for (const PileupSummaryInfo& it : *PileupInfos) {
		if(it.getBunchCrossing() == 0) {
			h1_TrueNumInteractions->Fill(it.getTrueNumInteractions());
			h1_TrueNumInteractions60 ->Fill(it.getTrueNumInteractions()); 
			h1_TrueNumInteractions70 ->Fill(it.getTrueNumInteractions()); 
			h1_TrueNumInteractions80 ->Fill(it.getTrueNumInteractions()); 
			h1_TrueNumInteractions90 ->Fill(it.getTrueNumInteractions()); 
			h1_TrueNumInteractions100->Fill(it.getTrueNumInteractions()); 
			h1_PU_NumInteractions->Fill(it.getPU_NumInteractions());
			break;
		}
	}
}

void MakeNpKNUHist::beginJob() {}
void MakeNpKNUHist::endJob() {}
void MakeNpKNUHist::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(MakeNpKNUHist);
