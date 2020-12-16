#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

class ProdNpGenerator : public edm::EDProducer {
   public:
      explicit ProdNpGenerator(const edm::ParameterSet&);
      ~ProdNpGenerator();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

		// GenParticles
      void FillDaughters(const reco::GenParticle& gen, int mIdx = -1, int printDepth = 0);
      bool DoDaughters(const int pdgId);
      void PrintGenParticle(const reco::GenParticle& ref);
      std::vector<std::pair<reco::GenParticle, int> > genVec;
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  genParticleToken ; bool DoGenParticle;
      void ProdGenParticle(edm::Event& iEvent, const edm::EventSetup& iSetup);
		bool PrintParticle;
		bool PrintGenParticleNum;
		unsigned prdGenParticleNum;

		// GenInfo
      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken ; bool DoGenInfo ;
      void ProdGenInfo(edm::Event& iEvent, const edm::EventSetup& iSetup);

		void Warning(std::string war);
		std::vector<std::string> WarningVec;
};

ProdNpGenerator::ProdNpGenerator(const edm::ParameterSet& iConfig) {

   DoGenParticle = (iConfig.getUntrackedParameter<bool>("DoGenParticle",false));
   if(DoGenParticle and iConfig.existsAs<edm::InputTag>("genParticles")) {
      genParticleToken = (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))) ;
      produces <std::vector<npknu::GenParticle> > ("GenParticle");
		PrintParticle = (iConfig.getUntrackedParameter<bool>("PrintParticle",false));
		PrintGenParticleNum = (iConfig.getUntrackedParameter<bool>("PrintGenParticleNum",false));
		prdGenParticleNum = 0;
   }

   DoGenInfo = (iConfig.getUntrackedParameter<bool>("DoGenInfo",false));
   if(DoGenInfo and iConfig.existsAs<edm::InputTag>("genInfoToken")) {
      genInfoToken = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoToken")) ;
      produces <std::vector<npknu::GenInfo> > ("GenInfo");
   }
}
ProdNpGenerator::~ProdNpGenerator() {}

void ProdNpGenerator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

	if(iEvent.isRealData()) return;

	if(DoGenParticle) ProdGenParticle(iEvent, iSetup);
	if(DoGenInfo) ProdGenInfo(iEvent, iSetup);

}

void ProdNpGenerator::ProdGenParticle(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   genVec.clear();
   edm::Handle<edm::View<reco::GenParticle> > genParticles; iEvent.getByToken(genParticleToken, genParticles);
   std::auto_ptr<std::vector<npknu::GenParticle> > genParticleProd ( new std::vector<npknu::GenParticle>() );
   //for(edm::View<reco::GenParticle>::const_iterator gen=genParticles->begin(); gen != genParticles->end(); ++gen) {
   // npknu::GenParticle genParti;
   // genParti.pt     = gen->pt();
   // genParti.eta    = gen->eta();
   // genParti.phi    = gen->phi();
   // genParti.energy = gen->energy();
   // genParti.pdgId  = gen->pdgId();
   // genParti.status = gen->status();
   // genParti.mass   = gen->mass();
   // genParti.nDau   = gen->numberOfDaughters();
   // genParticleProd->push_back(genParti);
   // prdGenParticleNum++;
   //}

   for(edm::View<reco::GenParticle>::const_iterator gen=genParticles->begin(); gen != genParticles->end(); ++gen) {
      if((std::abs(gen->pdgId()) == 22) and (gen->status() == 23)) { FillDaughters(*gen); }
      if((std::abs(gen->pdgId()) == 23) and (gen->status() == 22)) { FillDaughters(*gen); }
      if((std::abs(gen->pdgId()) == 24) and (gen->status() == 22)) { FillDaughters(*gen); }
      if((std::abs(gen->pdgId()) == 34) and (gen->status() == 22)) { FillDaughters(*gen); }
   }
	prdGenParticleNum = 0;
   for(std::vector<std::pair<reco::GenParticle, int> >::const_iterator gen=genVec.begin(); gen != genVec.end(); ++gen) {
      npknu::GenParticle genParti;
      if(PrintParticle) PrintGenParticle(gen->first);
      genParti.pt     = (gen->first).pt();
      genParti.eta    = (gen->first).eta();
      genParti.phi    = (gen->first).phi();
      genParti.energy = (gen->first).energy();
      genParti.pdgId  = (gen->first).pdgId();
      genParti.status = (gen->first).status();
      genParti.mass   = (gen->first).mass();
      genParti.nDau   = (gen->first).numberOfDaughters();
      genParti.mIdx   = gen->second;
      genParticleProd->push_back(genParti);
      prdGenParticleNum++;
   }
   iEvent.put(genParticleProd,"GenParticle");
   if(PrintGenParticleNum) { std::cout << "ProdNpGenerator NumberOfGenParticles " << genParticles->size() << " -> " << prdGenParticleNum << std::endl; }
}

bool ProdNpGenerator::DoDaughters(const int pdgId) {
   int ApdgId = std::abs(pdgId);
   if(ApdgId == 22 or ApdgId == 23 or ApdgId == 24) return true;
   if(ApdgId == 6) return true;
   if(ApdgId == 34) return true;
   if(ApdgId > 10 and ApdgId < 17) return true;
   return false;
}

void ProdNpGenerator::FillDaughters(const reco::GenParticle& gen, int mIdx, int printDepth) {
   int thisIndex = genVec.size();
   if(mIdx < 0) mIdx = thisIndex;
   if(gen.status() != 44) genVec.push_back(std::make_pair(gen,mIdx));

   if(PrintParticle) {
      if(gen.status() != 44) {
         for(int i=0; i<printDepth; i++) std::cout << "   " ;
         std::cout << "GenParticle " << (genVec.size() -1) << " " << gen.pdgId() << " " << gen.status() << " " << gen.numberOfDaughters() << " " << mIdx << std::endl;
      }
   }
   if(gen.status() != 44) printDepth++;
   for(size_t dIdx=0; dIdx<gen.numberOfDaughters(); dIdx++) {
      const reco::Candidate* dau = gen.daughter(dIdx);
      if(DoDaughters(gen.pdgId())) FillDaughters((reco::GenParticle&)(*dau), thisIndex, printDepth);
   }
}

void ProdNpGenerator::PrintGenParticle(const reco::GenParticle& ref) {
   std::cout << "GenParticle " << "status "<< ref.status() << " pdgId " << ref.pdgId() << " noM " << ref.numberOfMothers() << " noD " << ref.numberOfDaughters() << " " ;
   std::cout << "pt " << ref.pt() << " eta " << ref.eta() << " phi " << ref.phi() << " mass " << ref.mass() << std::endl;
}

void ProdNpGenerator::Warning(std::string war) {
	if(std::find(WarningVec.begin(), WarningVec.end(), war) == WarningVec.end()) {
		WarningVec.push_back(war);
	}
}

void ProdNpGenerator::ProdGenInfo(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   std::auto_ptr<std::vector<npknu::GenInfo> > genInfoProd ( new std::vector<npknu::GenInfo>() );
   edm::Handle<GenEventInfoProduct> GenEventInfoProducts;
   iEvent.getByToken(genInfoToken, GenEventInfoProducts);
	if(GenEventInfoProducts->binningValues().size()>0){
	   npknu::GenInfo genInfo;
      genInfo.scalePDF = GenEventInfoProducts->pdf()->scalePDF     ;
   	genInfo.pdg1     = GenEventInfoProducts->pdf()->id.first     ;
	   genInfo.pdg2     = GenEventInfoProducts->pdf()->id.second    ;
      genInfo.x1       = GenEventInfoProducts->pdf()->x.first      ;
   	genInfo.x2       = GenEventInfoProducts->pdf()->x.second     ;
	   genInfo.xpdf1    = GenEventInfoProducts->pdf()->xPDF.first   ;
      genInfo.xpdf2    = GenEventInfoProducts->pdf()->xPDF.second  ;
   	genInfoProd->push_back(genInfo);
	} else {
		std::stringstream buffer;
		buffer << "GenEventInfoProduct->binningValues().size() = " << GenEventInfoProducts->binningValues().size() << std::endl;
		Warning(buffer.str());
		//std::cout << "GenEventInfoProduct->binningValues().size() " << GenEventInfoProducts->binningValues().size() << std::endl;
	}
  	iEvent.put(genInfoProd,"GenInfo");
}

void ProdNpGenerator::beginJob() { }
void ProdNpGenerator::endJob() { 
	for(std::vector<std::string>::iterator it=WarningVec.begin(); it!=WarningVec.end(); it++) std::cout << "### Final Warning " << *it <<std::endl;
}
void ProdNpGenerator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpGenerator);
