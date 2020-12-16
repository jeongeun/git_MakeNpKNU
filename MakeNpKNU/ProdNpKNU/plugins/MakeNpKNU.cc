#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "TTree.h"
#include "TClonesArray.h"


class MakeNpKNU : public edm::EDAnalyzer {
   public:
      explicit MakeNpKNU(const edm::ParameterSet&);
      ~MakeNpKNU();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

		TTree* outTree;
		TClonesArray* evtTCA    ;

		bool doElectron;
		edm::EDGetTokenT<npknu::ElectronCollection>    electronToken     ;
		TClonesArray* electronTCA    ;

		bool doMuon;
		edm::EDGetTokenT<npknu::MuonCollection>    muonToken     ;
		TClonesArray* muonTCA    ;

		bool doJet;
		edm::EDGetTokenT<npknu::JetCollection>    jetToken     ;
		TClonesArray* jetTCA    ;

		bool doMET;
		edm::EDGetTokenT<npknu::METCollection>    metToken     ;
		TClonesArray* metTCA    ;

		bool doPhoton;
		edm::EDGetTokenT<npknu::PhotonCollection>   photonToken     ;
		TClonesArray* photonTCA    ;

		bool doGenInfo     ;
		edm::EDGetTokenT<npknu::GenInfoCollection>   genInfoToken     ;
		TClonesArray* genInfoTCA    ;

		bool doGenParticle ;
		edm::EDGetTokenT<npknu::GenParticleCollection>   genParticleToken     ;
		TClonesArray* genParticleTCA    ;

		bool doPileup      ;
		edm::EDGetTokenT<npknu::PileupCollection>   pileupToken     ;
		TClonesArray* pileupTCA    ;

		bool doTrigger     ;
		edm::EDGetTokenT<npknu::TriggerCollection>   triggerToken     ;
		TClonesArray* triggerTCA    ;
		edm::EDGetTokenT<npknu::TriggerObjectCollection>   triggerObjectToken     ;
		TClonesArray* triggerObjectTCA    ;

		bool doVertex      ;
		edm::EDGetTokenT<npknu::VertexCollection>   vertexToken     ;
		TClonesArray* vertexTCA    ;

		bool doMakeHist;
		std::string histFileName;
		TFile* histFile;

		npknu::ElectronHist* eleHistEB ;	
		npknu::ElectronHist* eleHistEE ;	

		npknu::MuonHist* muHistMB      ;
		npknu::MuonHist* muHistOL      ;
		npknu::MuonHist* muHistME      ;

		npknu::METHist* metHist        ;

		npknu::PhotonHist* phoHistEB   ;
		npknu::PhotonHist* phoHistEE   ;

		npknu::EtcHist* etcHist        ;
		int LoopNumber;

};
MakeNpKNU::MakeNpKNU(const edm::ParameterSet& iConfig) {
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree> ("NpKNU","NpKNUTree");

	evtTCA = new TClonesArray("npknu::Evt");
	outTree->Branch("evt",&evtTCA);

	doElectron = false;
	if(iConfig.existsAs<edm::InputTag>("electronToken")) {
		doElectron = true;
		electronToken = (consumes<npknu::ElectronCollection>    (iConfig.getParameter<edm::InputTag>("electronToken"))    ) ;
		electronTCA = new TClonesArray("npknu::Electron");
		outTree->Branch("electron",&electronTCA);
	}

	doMuon = false;
	if(iConfig.existsAs<edm::InputTag>("muonToken")) {
		doMuon = true;
		muonToken = (consumes<npknu::MuonCollection>    (iConfig.getParameter<edm::InputTag>("muonToken"))    ) ;
		muonTCA = new TClonesArray("npknu::Muon");
		outTree->Branch("muon",&muonTCA);
	}

	doJet = false;
	if(iConfig.existsAs<edm::InputTag>("jetToken")) {
		doJet = true;
		jetToken = (consumes<npknu::JetCollection>    (iConfig.getParameter<edm::InputTag>("jetToken"))    ) ;
		jetTCA = new TClonesArray("npknu::Jet");
		outTree->Branch("jet",&jetTCA);
	}

	doMET = false;
	if(iConfig.existsAs<edm::InputTag>("metToken")) {
		doMET = true;
		metToken = (consumes<npknu::METCollection>    (iConfig.getParameter<edm::InputTag>("metToken"))    ) ;
		metTCA = new TClonesArray("npknu::MET");
		outTree->Branch("met",&metTCA);
	}

	doPhoton = false;
	if(iConfig.existsAs<edm::InputTag>("photonToken")) {
		doPhoton = true;
		photonToken = (consumes<npknu::PhotonCollection>    (iConfig.getParameter<edm::InputTag>("photonToken"))    ) ;
		photonTCA = new TClonesArray("npknu::Photon");
		outTree->Branch("photon",&photonTCA);
	}

	doGenInfo = false;
	if(iConfig.existsAs<edm::InputTag>("genInfoToken")) {
		doGenInfo = true;
		genInfoToken = (consumes<npknu::GenInfoCollection>    (iConfig.getParameter<edm::InputTag>("genInfoToken"))    ) ;
		genInfoTCA = new TClonesArray("npknu::GenInfo");
		outTree->Branch("genInfo",&genInfoTCA);
	}

	doGenParticle = false;
	if(iConfig.existsAs<edm::InputTag>("genParticleToken")) {
		doGenParticle = true;
		genParticleToken = (consumes<npknu::GenParticleCollection>    (iConfig.getParameter<edm::InputTag>("genParticleToken"))    ) ;
		genParticleTCA = new TClonesArray("npknu::GenParticle");
		outTree->Branch("genParticle",&genParticleTCA);
	}

	doTrigger = false;
	if(iConfig.existsAs<edm::InputTag>("triggerToken")) {
		doTrigger = true;
		triggerToken = (consumes<npknu::TriggerCollection>    (iConfig.getParameter<edm::InputTag>("triggerToken"))    ) ;
		triggerTCA = new TClonesArray("npknu::Trigger");
		outTree->Branch("trigger",&triggerTCA);
		triggerObjectToken = (consumes<npknu::TriggerObjectCollection>    (iConfig.getParameter<edm::InputTag>("triggerObjectToken"))    ) ;
		triggerObjectTCA = new TClonesArray("npknu::TriggerObject");
		outTree->Branch("triggerObject",&triggerObjectTCA);
	}

	doPileup = false;
	if(iConfig.existsAs<edm::InputTag>("pileupToken")) {
		doPileup = true;
		pileupToken = (consumes<npknu::PileupCollection>    (iConfig.getParameter<edm::InputTag>("pileupToken"))    ) ;
		pileupTCA = new TClonesArray("npknu::Pileup");
		outTree->Branch("pileup",&pileupTCA);
	}

	doVertex = false;
	if(iConfig.existsAs<edm::InputTag>("vertexToken")) {
		doVertex = true;
		vertexToken = (consumes<npknu::VertexCollection>    (iConfig.getParameter<edm::InputTag>("vertexToken"))    ) ;
		vertexTCA = new TClonesArray("npknu::Vertex");
		outTree->Branch("vertex",&vertexTCA);
	}

	doMakeHist = false;
	if(iConfig.existsAs<std::string>("HistFileName")) {
		doMakeHist = true;
		histFileName = iConfig.getParameter<std::string>("HistFileName");
		histFile  = new TFile(histFileName.c_str(), "recreate");
		eleHistEB = new npknu::ElectronHist(histFile, "eleHistEB", true );	
		eleHistEE = new npknu::ElectronHist(histFile, "eleHistEE", false);
		muHistMB  = new npknu::MuonHist(histFile, "muHistMB", true );
		muHistOL  = new npknu::MuonHist(histFile, "muHistOL", false );
		muHistME  = new npknu::MuonHist(histFile, "muHistME", false );
		metHist   = new npknu::METHist(histFile, "metHist", true );
		phoHistEB = new npknu::PhotonHist(histFile, "phoHistEB", true );
		phoHistEE = new npknu::PhotonHist(histFile, "phoHistEE", false);
		etcHist   = new npknu::EtcHist(histFile, "etcHist" );
	}

	LoopNumber = 0;
}
MakeNpKNU::~MakeNpKNU() {
}

void MakeNpKNU::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
   //std::cout << "MakeNpKNU::analyze Started " << std::endl;
	LoopNumber++;

   evtTCA->Clear() ;
   npknu::Evt* evtPtr = new ((*evtTCA)[0]) npknu::Evt();
   evtPtr->run   = iEvent.id().run();
   evtPtr->lumi  = iEvent.id().luminosityBlock();
   evtPtr->event = iEvent.id().event();
   evtPtr->isRealData = iEvent.isRealData();

	if(doElectron) {
		electronTCA->Clear();
   	edm::Handle<npknu::ElectronCollection> electrons; iEvent.getByToken(electronToken, electrons);
	   for(const auto& it : *electrons) { 
			npknu::Electron* elePtr = new ((*electronTCA)[(int)electronTCA->GetEntries()]) npknu::Electron(it); 
			if(doMakeHist) {
				if(elePtr->isCutPassHEEP and elePtr->isEB()) eleHistEB->Fill(elePtr);
				if(elePtr->isCutPassHEEP and elePtr->isEE()) eleHistEE->Fill(elePtr);
			}
		}
	}

	if(doMuon) {
		muonTCA->Clear();
   	edm::Handle<npknu::MuonCollection> muons; iEvent.getByToken(muonToken, muons);
	   for(const auto& it : *muons) { 
		        npknu::Muon* muPtr = new ((*muonTCA)[(int)muonTCA->GetEntries()]) npknu::Muon(it); 
			if(doMakeHist) {
				if( muPtr->isMB() ){
					if((muPtr->statusGV > 110) or (muPtr->statusGV>10 and muPtr->statusGV< 20) ) muHistMB->Fill(muPtr);//highPTmuon
				}
				if( muPtr->isOverlap() ){
				        if((muPtr->statusGV > 110) or (muPtr->statusGV>10 and muPtr->statusGV< 20) ) muHistOL->Fill(muPtr);//highPTmuon
				}

				if( muPtr->isME() ){
				        if((muPtr->statusGV > 110) or (muPtr->statusGV>10 and muPtr->statusGV< 20) ) muHistME->Fill(muPtr);//highPTmuon
				}
			}
		}
	}

	if(doJet) {
		jetTCA->Clear();
   	edm::Handle<npknu::JetCollection> jets; iEvent.getByToken(jetToken, jets);
	   for(const auto& it : *jets) { new ((*jetTCA)[(int)jetTCA->GetEntries()]) npknu::Jet(it); }
	}

	if(doMET) {
		metTCA->Clear();
   	edm::Handle<npknu::METCollection> mets; iEvent.getByToken(metToken, mets);
	   for(const auto& it : *mets) { 
				npknu::MET* metPtr = new ((*metTCA)[(int)metTCA->GetEntries()]) npknu::MET(it);
			if(doMakeHist) {
				metHist->Fill(metPtr);
			}
		 }
	}

	if(doPhoton) {
		photonTCA->Clear();
   	edm::Handle<npknu::PhotonCollection> photons; iEvent.getByToken(photonToken, photons);
	   for(const auto& it : *photons) { 
			npknu::Photon* ptr = new ((*photonTCA)[(int)photonTCA->GetEntries()]) npknu::Photon(it); 
			if(doMakeHist) {
				if(ptr->isCutPassTight and ptr->isEB()) phoHistEB->Fill(ptr);
				if(ptr->isCutPassTight and ptr->isEE()) phoHistEE->Fill(ptr);
			}
		}
	}

	if(doGenInfo and !iEvent.isRealData()) {
		genInfoTCA->Clear();
   	edm::Handle<npknu::GenInfoCollection> genInfos; iEvent.getByToken(genInfoToken, genInfos);
	   for(const auto& it : *genInfos) { new ((*genInfoTCA)[(int)genInfoTCA->GetEntries()]) npknu::GenInfo(it); }
	}

	if(doGenParticle and !iEvent.isRealData()) {
		genParticleTCA->Clear();
   	edm::Handle<npknu::GenParticleCollection> genParticles; iEvent.getByToken(genParticleToken, genParticles);
	   for(const auto& it : *genParticles) { new ((*genParticleTCA)[(int)genParticleTCA->GetEntries()]) npknu::GenParticle(it); }
	}

	if(doPileup and !iEvent.isRealData()) {
		pileupTCA->Clear();
   	edm::Handle<npknu::PileupCollection> pileups; iEvent.getByToken(pileupToken, pileups);
	   for(const auto& it : *pileups) {  
                             npknu::Pileup* PUptr = new ((*pileupTCA)[(int)pileupTCA->GetEntries()]) npknu::Pileup(it); 
                             
                             if(doMakeHist){
                                    if(PUptr->BunchCrossing != 0){
                                        double pileupTrue =(PUptr->BunchCrossing == 0)? -999 : PUptr->TrueNumInteractions ; 
 				        int    pileupPU   =(PUptr->BunchCrossing == 0)? -999 : PUptr->PU_NumInteractions  ; 
				        etcHist->Fill(pileupTrue, pileupPU);
				     }
                             }
		}
	}

	if(doTrigger) {
		triggerTCA->Clear();
   	edm::Handle<npknu::TriggerCollection> triggers; iEvent.getByToken(triggerToken, triggers);
	   for(const auto& it : *triggers) { new ((*triggerTCA)[(int)triggerTCA->GetEntries()]) npknu::Trigger(it); }
		triggerObjectTCA->Clear();
   	edm::Handle<npknu::TriggerObjectCollection> triggerObjects; iEvent.getByToken(triggerObjectToken, triggerObjects);
	   for(const auto& it : *triggerObjects) { new ((*triggerObjectTCA)[(int)triggerObjectTCA->GetEntries()]) npknu::TriggerObject(it); }
	}

	if(doVertex) {
		vertexTCA->Clear();
   	edm::Handle<npknu::VertexCollection> vertexs; iEvent.getByToken(vertexToken, vertexs);
	   for(const auto& it : *vertexs) { new ((*vertexTCA)[(int)vertexTCA->GetEntries()]) npknu::Vertex(it); }
	}

	outTree->Fill();

	if(LoopNumber < 10) {
		std::cout << "------- number " << LoopNumber << std::endl;
   	if(doElectron){ std::cout << "@MakeNpKNU ElectronSize  " << electronTCA->GetEntries() << std::endl; }	
   	if(doMuon)    { std::cout << "@MakeNpKNU MuonSize      " << muonTCA->GetEntries()     << std::endl; }
   	if(doMET)     { std::cout << "@MakeNpKNU METSize       " << metTCA->GetEntries()      << std::endl; }
   	if(doJet)     { std::cout << "@MakeNpKNU JetSize       " << jetTCA->GetEntries()      << std::endl; }
   	if(doPhoton)  { std::cout << "@MakeNpKNU PhotonSize    " << photonTCA->GetEntries()   << std::endl; }
   	if(doTrigger) { std::cout << "@MakeNpKNU TriggerSize   " << triggerTCA->GetEntries()  << std::endl; }
   	if(doTrigger) { std::cout << "@MakeNpKNU TriggerObjSize" << triggerObjectTCA->GetEntries() << std::endl; }
   	if(doVertex)  { std::cout << "@MakeNpKNU VertexSize    " << vertexTCA->GetEntries()   << std::endl; }
   	if(doGenInfo and !iEvent.isRealData())    { std::cout << "@MakeNpKNU GenInfoSize     "<< genInfoTCA->GetEntries()     << std::endl; }
   	if(doGenParticle and !iEvent.isRealData()){ std::cout << "@MakeNpKNU GenParticleSize "<< genParticleTCA->GetEntries() << std::endl; }
   	if(doPileup and !iEvent.isRealData())     { std::cout << "@MakeNpKNU PileUpSize      "<< pileupTCA->GetEntries()      << std::endl; } 
	}
 
}

void MakeNpKNU::beginJob() { }
void MakeNpKNU::endJob() { 
	if(doMakeHist) {
		std::cout << "endJob Write HistFile " << std::endl;
		eleHistEB->Write(); 
		eleHistEE->Write(); 
		muHistMB->Write(); 
		muHistOL->Write(); 
		muHistME->Write(); 
		metHist  ->Write(); 
		//phoHistEB->Write(); 
		//phoHistEE->Write(); 
		etcHist  ->Write(); 
		histFile ->Write();
	}
}
void MakeNpKNU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(MakeNpKNU);
