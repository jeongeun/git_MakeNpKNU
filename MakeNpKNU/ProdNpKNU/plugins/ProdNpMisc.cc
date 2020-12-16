#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


class ProdNpMisc : public edm::EDProducer {
   public:
      explicit ProdNpMisc(const edm::ParameterSet&);
      ~ProdNpMisc();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

		// Pileup
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupToken;   bool DoPileup ;
      void ProdPileup(edm::Event& iEvent, const edm::EventSetup& iSetup);

		//Vertex
      edm::EDGetTokenT<reco::VertexCollection> vertexToken ; bool DoVertex ;
      void ProdVertex(edm::Event& iEvent, const edm::EventSetup& iSetup) ;

		// Trigger
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
      std::vector<std::string> triggerIdentifiers; bool DoTrigger;
		void ProdTrigger(edm::Event& iEvent, const edm::EventSetup& iSetup);
		bool PrintTriggerResult;

      unsigned LoopNumber;
};

ProdNpMisc::ProdNpMisc(const edm::ParameterSet& iConfig) {

   DoPileup = (iConfig.getUntrackedParameter<bool>("DoPileup",false));
   if(DoPileup and iConfig.existsAs<edm::InputTag>("pileupToken")) {
      pileupToken = (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupToken")))   ;
      produces <std::vector<npknu::Pileup> > ("Pileup");
		std::cout << "Doing ProdPileup " << std::endl;
   }

   DoVertex = (iConfig.getUntrackedParameter<bool>("DoVertex",false));
   if(DoVertex and iConfig.existsAs<edm::InputTag>("vertexToken")) {
      vertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexToken")) ;
      produces <std::vector<npknu::Vertex> > ("Vertex");
		std::cout << "Doing ProdVertex " << std::endl;
   }

	DoTrigger = (iConfig.getUntrackedParameter<bool>("DoTrigger",false));;
	if(DoTrigger and iConfig.existsAs<std::vector<std::string > >("triggerIdentifiers")) {
	   triggerBitsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBitsToken"));
	   triggerObjectsToken = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjectsToken"));
	   triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescalesToken"));
	   triggerIdentifiers = iConfig.getParameter<std::vector<std::string> >("triggerIdentifiers");
		PrintTriggerResult = (iConfig.getUntrackedParameter<bool>("PrintTriggerResult",false));
   	produces <std::vector<npknu::Trigger> > ("Trigger");
   	produces <std::vector<npknu::TriggerObject> > ("TriggerObject");
		std::cout << "Doing ProdTrigger " << std::endl;
	}

   LoopNumber = 0;
 
}
ProdNpMisc::~ProdNpMisc() {}

void ProdNpMisc::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
   LoopNumber++;

	if(DoPileup && (!iEvent.isRealData()))  ProdPileup(iEvent, iSetup);
	if(DoVertex)  ProdVertex(iEvent, iSetup);
	if(DoTrigger) ProdTrigger(iEvent, iSetup);
}

void ProdNpMisc::ProdTrigger(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerBitsToken, triggerBits);
   iEvent.getByToken(triggerObjectsToken, triggerObjects);
   iEvent.getByToken(triggerPrescalesToken, triggerPrescales);

   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   if(LoopNumber < 5) { //for printing trigger path at only first event
      std::cout << "\n === TRIGGER PATHS === " << std::endl;
      for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
      }

      std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
           obj.unpackPathNames(names);
           std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
           // Print trigger object collection and type
           std::cout << "\t   Collection: " << obj.collection() << std::endl;
           std::cout << "\t   Type IDs:   ";
           for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] << std::endl;
           std::cout << "\t   Filters:    ";  // Print associated trigger filters  
           for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h]<< std::endl;
           std::vector<std::string> pathNamesAll  = obj.pathNames(false);
           std::vector<std::string> pathNamesLast = obj.pathNames(true);
           // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
           // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
           // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
           std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
           for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
               bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
               bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
               bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
               bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
               std::cout << "   " << pathNamesAll[h];
               if (isBoth) std::cout << "(L,3)";
               if (isL3 && !isBoth) std::cout << "(*,3)";
               if (isLF && !isBoth) std::cout << "(L,*)";
               if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
           }
           std::cout << std::endl;
       }//triggerObject loop over


   }//1 loop over
	using std::cout;
	using std::endl;

   std::auto_ptr<std::vector<npknu::Trigger> > triggerProd ( new std::vector<npknu::Trigger>() );
   for (unsigned int j = 0; j < triggerIdentifiers.size(); ++j) {
      std::string idName = triggerIdentifiers[j];
      std::string idNameUnstarred = idName;
      bool isStarred = (idName.find("*")!=std::string::npos);
      if( isStarred ) idNameUnstarred.erase( idName.find("*"), 1 );
      for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
         if( (isStarred && names.triggerName(i).find(idNameUnstarred)!=std::string::npos ) || (!isStarred && names.triggerName(i)==idName)) {
            if(PrintTriggerResult) std::cout << "TriggerName " << names.triggerName(i) << " Accept " << triggerBits->accept(i) << " PrescaleIndex " << triggerPrescales->getPrescaleForIndex(i) << std::endl;
            npknu::Trigger trigger;
            trigger.name = names.triggerName(i);
            trigger.accept = triggerBits->accept(i);
            trigger.prescaleForIndex = triggerPrescales->getPrescaleForIndex(i);
            //triggerProd->push_back(trigger);
            if(triggerBits->accept(i) ) triggerProd->push_back(trigger);
         }
      }
   }
   iEvent.put(triggerProd,"Trigger");

   std::auto_ptr<std::vector<npknu::TriggerObject> > triggerObjectProd ( new std::vector<npknu::TriggerObject>() );
	for(pat::TriggerObjectStandAlone obj : *triggerObjects) { 
	   obj.unpackPathNames(names);
	   npknu::TriggerObject triObj;
	   triObj.pt = obj.pt();
	   triObj.eta = obj.eta();
	   triObj.phi = obj.phi();
	
	   std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	   std::vector<std::string> pathNamesLast = obj.pathNames(true);
		for(unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
         bool isBoth = obj.hasPathName( pathNamesLast[h], true, true );
         bool isL3   = obj.hasPathName( pathNamesLast[h], false, true );
         bool isLF   = obj.hasPathName( pathNamesLast[h], true, false );
         bool isNone = obj.hasPathName( pathNamesLast[h], false, false );
         triObj.pathNameVec.push_back( pathNamesLast[h]);
         triObj.isBoth.push_back(isBoth);
         triObj.isL3.push_back(isL3);
         triObj.isLF.push_back(isLF);
         triObj.isNone.push_back(isNone);
		}
      triggerObjectProd->push_back(triObj);
   }
	
   iEvent.put(triggerObjectProd,"TriggerObject");
}


void ProdNpMisc::ProdPileup(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   std::auto_ptr<std::vector<npknu::Pileup> > pileupProd ( new std::vector<npknu::Pileup>() );
   edm::Handle<std::vector<PileupSummaryInfo> > PileupInfos;
   iEvent.getByToken(pileupToken, PileupInfos);
   for (const PileupSummaryInfo& it : *PileupInfos) {
      npknu::Pileup pileup;
      pileup.BunchCrossing       = it.getBunchCrossing();
      pileup.TrueNumInteractions = it.getTrueNumInteractions();
      pileup.PU_NumInteractions  = it.getPU_NumInteractions();
      pileupProd->push_back(pileup);
   }
   iEvent.put(pileupProd,"Pileup");
}

void ProdNpMisc::ProdVertex(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   std::auto_ptr<std::vector<npknu::Vertex> > vertexProd ( new std::vector<npknu::Vertex>() );
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken, vertices);
   for(const reco::Vertex& it : *vertices) {
      npknu::Vertex vtx;
      vtx.x              = it.x()              ;
      vtx.y              = it.y()              ;
      vtx.z              = it.z()              ;
      vtx.xError         = it.xError()         ;
      vtx.yError         = it.yError()         ;
      vtx.zError         = it.zError()         ;
      vtx.tracksSize     = it.tracksSize()     ;
      vtx.nTracks        = it.nTracks()        ;
      vtx.isFake         = it.isFake()         ;
      vtx.ndof           = it.ndof()           ;
      vtx.position_rho   = it.position().rho() ;
      vtx.chi2           = it.chi2()           ;
      vtx.normalizedChi2 = it.normalizedChi2() ;
      vertexProd->push_back(vtx);
   }
   iEvent.put(vertexProd,"Vertex");
}



void ProdNpMisc::beginJob() { }
void ProdNpMisc::endJob() { }
void ProdNpMisc::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpMisc);
