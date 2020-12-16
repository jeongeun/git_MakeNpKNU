#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "DataFormats/PatCandidates/interface/Jet.h"

class ProdNpJet : public edm::EDProducer {
   public:
      explicit ProdNpJet(const edm::ParameterSet&);
      ~ProdNpJet();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
		double MinPtCut;
		bool PFJetId(const pat::Jet& pfjet, int idType);
};

ProdNpJet::ProdNpJet(const edm::ParameterSet& iConfig) {
   jetToken = (consumes<edm::View<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))) ;
	produces <std::vector<npknu::Jet> > ("Jet");
	MinPtCut = (iConfig.getUntrackedParameter<double>("MinPtCut",0.0));
	std::cout << "ProdNpJet MinPtCut " << MinPtCut << std::endl;
}

ProdNpJet::~ProdNpJet() { }
void ProdNpJet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

	std::auto_ptr<std::vector<npknu::Jet> > jetProd ( new std::vector<npknu::Jet>() );
   edm::Handle<edm::View<pat::Jet> > jets;
   iEvent.getByToken(jetToken, jets);
   for(edm::View<pat::Jet>::const_iterator jetPtr = jets->begin(); jetPtr != jets->end(); ++jetPtr){
      if(jetPtr->pt() < MinPtCut) continue;

		npknu::Jet jet;
		jet.pt     = jetPtr->pt();
		jet.eta    = jetPtr->eta();
		jet.phi    = jetPtr->phi();
		jet.energy = jetPtr->energy();
		jet.neutralHadronEnergyFraction   = jetPtr->neutralHadronEnergyFraction() ;
		jet.neutralEmEnergyFraction       = jetPtr->neutralEmEnergyFraction()     ;
		jet.chargedMultiplicity           = jetPtr->chargedMultiplicity()         ;
		jet.neutralMultiplicity           = jetPtr->neutralMultiplicity();
		jet.muonEnergyFraction            = jetPtr->muonEnergyFraction();
		jet.chargedHadronEnergyFraction   = jetPtr->chargedHadronEnergyFraction();
		jet.chargedEmEnergyFraction       = jetPtr->chargedEmEnergyFraction();

		jet.passPFLooseId  = PFJetId(*jetPtr, 0);
		jet.passPFMediumId = PFJetId(*jetPtr, 1);
		jet.passPFTightId  = PFJetId(*jetPtr, 2);

	
		jetProd->push_back(jet);
   }
	iEvent.put(jetProd,"Jet");

}

bool ProdNpJet::PFJetId(const pat::Jet& pfjet, int idType){ // idType==0 for Loose,  1 for Medium, 2 for Tight
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID Loose(Recommended)

	float neuHadFrac = pfjet.neutralHadronEnergyFraction();
	float neuEMFrac  = pfjet.neutralEmEnergyFraction();
	float numConst   = pfjet.chargedMultiplicity() + pfjet.neutralMultiplicity();
	float muFrac     = pfjet.muonEnergyFraction();

	float chHadFrac = pfjet.chargedHadronEnergyFraction();
	float chMul     = pfjet.chargedMultiplicity(); 
	float chEMFrac  = pfjet.chargedEmEnergyFraction();

	float CutNeuHadFrac[3] = {0.99, 0.95, 0.90}; 
	float CutNeuEMFrac[3]  = {0.99, 0.95, 0.90}; 
	float CutNumConst[3]   = {0, 0, 0}; 
	float CutMuFrac[3]     = {0.8, 0.8, 0.8}; 
	float CutChHadFrac[3]  = {1, 1, 1}; 
	float CutChMul[3]      = {0, 0, 0}; 
	float CutChEMFrac[3]   = {0.99, 0.99, 0.99}; 

   bool passPFJetId = false;

   if (      (neuHadFrac < CutNeuHadFrac[idType] )
         and (neuEMFrac  < CutNeuEMFrac[idType] )
         and (numConst   < CutNumConst[idType] )
         and (muFrac     < CutMuFrac[idType] )
         and (fabs(pfjet.eta()) > 2.4 || chHadFrac > CutChHadFrac[idType] )
         and (fabs(pfjet.eta()) > 2.4 || chMul     > CutChMul[idType]     )
         and (fabs(pfjet.eta()) > 2.4 || chEMFrac  < CutChEMFrac[idType]  )
   ) passPFJetId = true;
   return passPFJetId;


// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopJME
}




void ProdNpJet::beginJob() { }
void ProdNpJet::endJob() { } 
void ProdNpJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpJet);
