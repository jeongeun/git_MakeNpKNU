#include <memory>
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

#include "Math/VectorUtil.h"

#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//these header files give us easy to use shortcuts for which cut
////corresponds to which which cutnr
////this is fixed for a given ID (and can be different for each ID)
////hence its hard coded
////also these headerfiles are intentionally completely standalone 
////so you can easily include them in your analysis if you find them
////useful
#include "HEEP/VID/interface/CutNrs.h"
#include "HEEP/VID/interface/VIDCutCodes.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"

#include "TH1F.h"
//2015-July27
//https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_747/ElectronNtupler/plugins/SimpleElectronNtupler.cc
//2016-Feb28
//https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_electrons.cc
//2016-Nov10
//https://github.com/Sam-Harper/usercode/blob/HEEPV70Example/HEEPAnalyzer/plugins/HEEPV70Example.cc
//https://github.com/Sam-Harper/usercode/blob/HEEPV70Example/HEEPAnalyzer/test/HEEPV70Example_cfg.py


//**********************************************************
////
//// class: HEEPV70Example
////
//// author: Sam Harper (RAL)
////
//// this class is meant to serve as an example of how to 
//// access the HEEP V7.0 ID + new tracker isolation 
//// using the value maps calculated by VID

// the following value maps are of interest
// //   edm::ValueMap<vid::CutFlowResult>  egmGsfElectronIDs:heepElectronID-HEEPV70  :
// //        the full VID result with lots of info
// //   edm::ValueMap<unsigned int>  egmGsfElectronIDs:heepElectronID-HEEPV70Bitmap  :
// //        bitmap of which cuts passed, bitNr X = cut X and 0=fail, 1 =pass
// //   edm::ValueMap<bool>  egmGsfElectronIDs:heepElectronID-HEEPV70 
// //        global pass/fail bool,  0=fail HEEPV70, 1 =pass HEEPV70
// //   edm::ValueMap<int> heepIDVarValueMaps:eleNrSaturateIn5x5  :
// //        specific to HEEEP ID (not technically part of VID but run by VID automatically)
// //        nr of staturated crystals in 5x5 centred on the seed crystal
// //   edm::ValueMap<int> heepIDVarValueMaps:eleNrSaturateIn5x5  :
// //        specific to HEEEP ID (not technically part of VID but run by VID automatically)
// //        the new tracker isolation used in the ID
// //

//namespace{
//  struct NrPassFail {
//    NrPassFail():nrPass(0),nrFail(0){}
//    mutable std::atomic<int> nrPass;
//    mutable std::atomic<int> nrFail;
//  };
//}


class ProdNpElectron : public edm::EDProducer {
   public:
      explicit ProdNpElectron(const edm::ParameterSet&);
      ~ProdNpElectron();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

		enum ElectronMatchType {UNMATCHED = 0,
           TRUE_PROMPT_ELECTRON,
           TRUE_ELECTRON_FROM_TAU,
           TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

		int matchToTruth(const edm::Ptr<reco::GsfElectron> el, const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
		void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
		void printCutFlowResult(vid::CutFlowResult &cutflow);

                edm::EDGetTokenT<edm::View<pat::Electron> > electronToken ; 
                //edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronToken;
                edm::EDGetTokenT<reco::ConversionCollection> conversionToken ;
                edm::EDGetTokenT<reco::BeamSpot> beamSpotToken ;
                edm::EDGetTokenT<edm::SortedCollection<EcalRecHit> > ebReducedRecHitToken;
                edm::EDGetTokenT<edm::SortedCollection<EcalRecHit> > eeReducedRecHitToken;
                edm::EDGetTokenT<edm::SortedCollection<EcalRecHit> > esReducedRecHitToken;
                edm::EDGetTokenT<reco::VertexCollection> vertexToken ;
		edm::EDGetTokenT<double> rhoToken;

		double MinPtCut;
		int LoopEvent;    

		bool doGenMatch;
		edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleToken;

                //edm::EDGetTokenT<edm::ValueMap<bool> > 	mvaMediumIdMapToken   ; 
                //edm::EDGetTokenT<edm::ValueMap<bool> > 	mvaTightIdMapToken    ; 
                //edm::EDGetTokenT<edm::ValueMap<float> >	mvaValuesMapToken     ; 
                //edm::EDGetTokenT<edm::ValueMap<int> >  	mvaCategoriesMapToken ; 

		//edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenTight ;
		//edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenMedium;
		//edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenLoose ;
		//edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenVeto  ;
		//edm::EDGetTokenT<edm::ValueMap<bool> > HEEP70IdToken  ; //heep70
                //edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > HEEP60IdFullInfoMapToken ; //heepcutflow
                edm::EDGetTokenT<edm::ValueMap<bool> > HEEP70IdToken  ;
                edm::EDGetTokenT<edm::ValueMap<unsigned int> > HEEP70IdBitmapToken;
                edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> >  HEEP70IdResultToken ;
                
                edm::EDGetTokenT<edm::ValueMap<int> > nrSatCrysMapToken; 
                edm::EDGetTokenT<edm::ValueMap<float> > trkIsolMapToken;   

};

ProdNpElectron::ProdNpElectron(const edm::ParameterSet& iConfig) {
	electronToken        = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronToken"));
	//electronToken        = consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronToken"));
	conversionToken      = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversionToken")) ;
	beamSpotToken        = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotToken")) ;
	ebReducedRecHitToken = consumes<edm::SortedCollection<EcalRecHit> > (iConfig.getParameter<edm::InputTag>("ebReducedRecHitToken"));
	eeReducedRecHitToken = consumes<edm::SortedCollection<EcalRecHit> > (iConfig.getParameter<edm::InputTag>("eeReducedRecHitToken"));
	esReducedRecHitToken = consumes<edm::SortedCollection<EcalRecHit> > (iConfig.getParameter<edm::InputTag>("esReducedRecHitToken"));
	vertexToken          = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexToken")) ;
	rhoToken             = (consumes<double> (iConfig.getParameter<edm::InputTag>("rhoToken")))                 ;

	MinPtCut = (iConfig.getUntrackedParameter<double>("MinPtCut",10.0));	

	doGenMatch = false;
   if(iConfig.existsAs<edm::InputTag>("genParticles")) {
		doGenMatch = true;
        genParticleToken = (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))) ;
	}
	
	//mvaMediumIdMapToken   = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaMediumIdMap"))) ;
	//mvaTightIdMapToken    = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaTightIdMap")))  ;
	//mvaValuesMapToken     = (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap")))  ;
	//mvaCategoriesMapToken = (consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))) ;

	//cutBasedIdTokenTight  = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenTight"))) ;
	//cutBasedIdTokenMedium = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenMedium"))) ;
	//cutBasedIdTokenLoose  = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenLoose"))) ;
	//cutBasedIdTokenVeto   = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenVeto"))) ;

//        //HEEP60IdFullInfoMapToken = consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("HEEP60IdFullInfoMap"));
//        trkIsolMapToken       = (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trkIsolMap")));
	HEEP70IdToken         = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("HEEP70Id"))) ;
        HEEP70IdBitmapToken= (consumes<edm::ValueMap<unsigned int> >(iConfig.getParameter<edm::InputTag>("HEEP70IdBitmap"))) ;
        HEEP70IdResultToken= (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("HEEP70Id"))) ;
        nrSatCrysMapToken  = (consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("nrSatCrysMap"))) ;
        trkIsolMapToken    = (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trkIsolMap"))) ;

	produces <npknu::ElectronCollection> ("Electron");
	std::cout << "ProdNpElectron MinPt " << MinPtCut << std::endl;

	LoopEvent = 0;

}

ProdNpElectron::~ProdNpElectron(){}

void ProdNpElectron::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
   using namespace std;
   using namespace reco;

	using std::cout;
	using std::endl;

	LoopEvent++;

   edm::Handle<double> rhoHandle;
	iEvent.getByToken(rhoToken, rhoHandle); 
   double rhoValue = (*(rhoHandle.product()));

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken, vertices);
   if (vertices->empty()) return; //skip the event if no PV found
   reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
   int firstGoodVertexIdx = 0;
   for(reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
       // Replace isFake() for miniAOD because 
       // it requires tracks and miniAOD vertices don't have tracks:
       // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
       // bool isFake = vtx->isFake();
       bool isFake =  (vtx->chi2()==0 && vtx->ndof()==0);
       if ( !isFake &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
	      firstGoodVertex = vtx;
   	   break;
           }
    }
   if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs
   //std::cout << "good PVs , firstGoodVertexIdx " << firstGoodVertexIdx << std::endl;

   edm::Handle<reco::ConversionCollection> conversions;
   iEvent.getByToken(conversionToken, conversions);

   edm::Handle<reco::BeamSpot> beamSpot_h;
   iEvent.getByToken(beamSpotToken, beamSpot_h);
   const reco::BeamSpot& beamSpot = *(beamSpot_h.product());

   noZS::EcalClusterLazyTools* lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebReducedRecHitToken, eeReducedRecHitToken, esReducedRecHitToken); 

   edm::Handle<edm::View<reco::GenParticle> > genParticles; 
	if((!iEvent.isRealData()) and doGenMatch) {
	   iEvent.getByToken(genParticleToken, genParticles);
	} else {
		doGenMatch = false;
	}

	//edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsTight ;
	//edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsMedium;
	//edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsLoose ;
	//edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsVeto  ;

        edm::Handle<edm::ValueMap<bool> >  heep70Id ;
        edm::Handle<edm::ValueMap<unsigned int> > heep70IdBitmap;
        edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep70IdResult;

        edm::Handle<edm::ValueMap<int> > nrSatCrysMap;
        edm::Handle<edm::ValueMap<float> > trkIsolMap;


        iEvent.getByToken(HEEP70IdToken, heep70Id);
        iEvent.getByToken(HEEP70IdBitmapToken,heep70IdBitmap);
        iEvent.getByToken(HEEP70IdResultToken,heep70IdResult);
        iEvent.getByToken(trkIsolMapToken,trkIsolMap);
        iEvent.getByToken(nrSatCrysMapToken,nrSatCrysMap);


	//iEvent.getByToken(cutBasedIdTokenTight  ,cutIdDecisionsTight );
	//iEvent.getByToken(cutBasedIdTokenMedium ,cutIdDecisionsMedium);
	//iEvent.getByToken(cutBasedIdTokenLoose  ,cutIdDecisionsLoose );
	//iEvent.getByToken(cutBasedIdTokenVeto   ,cutIdDecisionsVeto  );

	//edm::Handle<edm::ValueMap<bool> >  mvaMediumIdDecisions ;
	//edm::Handle<edm::ValueMap<bool> >  mvaTightIdDecisions ;
	//edm::Handle<edm::ValueMap<float> > mvaValues           ;
	//edm::Handle<edm::ValueMap<int> >   mvaCategories       ;
	//iEvent.getByToken(mvaMediumIdMapToken   ,mvaMediumIdDecisions);
	//iEvent.getByToken(mvaTightIdMapToken    ,mvaTightIdDecisions) ;
	//iEvent.getByToken(mvaValuesMapToken     ,mvaValues)           ;
	//iEvent.getByToken(mvaCategoriesMapToken ,mvaCategories)       ;
        // Full cut flow info for one of the working points:
        //edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep60_Id_cutflow_data;
        //iEvent.getByToken(HEEP60IdFullInfoMapToken,heep60_Id_cutflow_data);


	std::auto_ptr<npknu::ElectronCollection> electronProd ( new npknu::ElectronCollection() );
   edm::Handle<edm::View<pat::Electron> > electrons;
   //edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(electronToken, electrons);
//   for(edm::View<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); el++) {
   for(size_t eleNr=0; eleNr<electrons->size(); eleNr++){
      edm::Ptr<pat::Electron> el(electrons,eleNr);

      if(el->pt() < MinPtCut) continue;
      npknu::Electron ele;
      /////https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
      /////https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
      ele.pt      = el->pt()     ;      
      ele.eta     = el->eta()    ;      
      ele.phi     = el->phi()    ;      
      ele.energy  = el->energy() ;      
      ele.gsfTrack_Px      = el->gsfTrack()->px()  ;  
      ele.gsfTrack_Py      = el->gsfTrack()->py()  ;  
      ele.gsfTrack_Pz      = el->gsfTrack()->pz()  ;  
      ele.gsfTrack_Pt      = el->gsfTrack()->pt()  ;  
      ele.superCluster_x   = el->superCluster()->x()   ; 
      ele.superCluster_y   = el->superCluster()->y()   ; 
      ele.superCluster_z   = el->superCluster()->z()   ; 
      ele.superCluster_eta = el->superCluster()->eta() ; 
      ele.superCluster_phi = el->superCluster()->phi() ; 
      //add preshower
      ele.superCluster_ESenergy = el->superCluster()->preshowerEnergy();
      ele.superCluster_ESplane1 = el->superCluster()->preshowerEnergyPlane1();
      ele.superCluster_ESplane2 = el->superCluster()->preshowerEnergyPlane2();

      ele.charge              = el->charge()                 ;
      ele.caloEnergy          = el->caloEnergy()             ;
      ele.ecalEnergy          = el->ecalEnergy()             ;
      ele.superCluster_energy = el->superCluster()->energy() ;

      ele.scE1x5          = el->scE1x5()         ;
      ele.scE2x5Max       = el->scE2x5Max()      ;
      ele.scE5x5          = el->scE5x5()         ;
      ele.scPixCharge     = el->scPixCharge()    ;
      ele.scSigmaEtaEta   = el->scSigmaEtaEta()  ;
      ele.scSigmaIEtaIEta = el->scSigmaIEtaIEta();

      ele.e2x5Max = el->e2x5Max() ;
      ele.e5x5    = el->e5x5()    ;
      ele.e1x5    = el->e1x5()    ;
      ele.full5x5_sigmaIetaIeta   = el->full5x5_sigmaIetaIeta() ;
      ele.sigmaIetaIphi           = el->sigmaIetaIphi()         ;
      ele.sigmaIphiIphi           = el->sigmaIphiIphi()         ;

      ele.gsfTrack_dxy      = el->gsfTrack()->dxy()                ;
      ele.gsfTrack_dz       = el->gsfTrack()->dz()                 ;

      ele.gsfTrack_dxyPVtx  = (!vertices->empty()) ? el->gsfTrack()->dxy(vertices->front().position()) : el->gsfTrack()->dxy()  ;
      ele.gsfTrack_dzPVtx   = (!vertices->empty()) ? el->gsfTrack()->dz(vertices->front().position())  : el->gsfTrack()->dz()  ;

      ele.dr03EcalRecHitSumEt              = el->dr03EcalRecHitSumEt()            ;
      ele.dr03HcalDepth1TowerSumEt         = el->dr03HcalDepth1TowerSumEt()       ;
      //ele.dr03HcalDepth2TowerSumEt         = el->dr03HcalDepth2TowerSumEt()       ;
      ele.dr03HcalTowerSumEt               = el->dr03HcalTowerSumEt()             ;
      ele.dr03TkSumPt                      = el->dr03TkSumPt()                    ;

      ele.ecalDrivenSeed                   = (bool)el->ecalDrivenSeed()           ;
      ele.deltaEtaSuperClusterTrackAtVtx   = el->deltaEtaSuperClusterTrackAtVtx() ;
      ele.deltaEtaSeedClusterTrackAtCalo   = el->deltaEtaSeedClusterTrackAtCalo() ;
      ele.deltaEtaSeedClusterTrackAtVtx    = el->deltaEtaSeedClusterTrackAtVtx() ;
      ele.deltaPhiSuperClusterTrackAtVtx   = el->deltaPhiSuperClusterTrackAtVtx() ;

      //HEEP Id variable delta Eta seed at vertex
      double eledEtaSeedAtVtx = el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() ?
      el->deltaEtaSuperClusterTrackAtVtx() - el->superCluster()->eta() + el->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

      ele.dEtaSeedAtVtx = eledEtaSeedAtVtx;

      ele.hcalOverEcal                     = el->hcalOverEcal()                   ;
      ele.hadronicOverEm                   = el->hadronicOverEm()                 ;
      ele.nMissingHits            = el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
      ele.passConversionVeto      = el->passConversionVeto()   ;
      ele.fbrem                   = el->fbrem()                ;

      ele.caloIso                 = el->caloIso()              ;
      ele.ecalIso                 = el->ecalIso()              ;
      ele.hcalIso                 = el->hcalIso()              ;
      ele.trackIso                = el->trackIso()             ;

      ele.isPF                    = el->isPF()                 ;
      ele.r9                      = el->r9()                   ;

      ele.eSuperClusterOverP =  el->eSuperClusterOverP() ;
      //ele.ooEmooP            = (1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy());

    if( el->ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ele.ooEmooP = ( 1e30 );
    }else if( !std::isfinite(el->ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ele.ooEmooP = ( 1e30 );
    }else{
      ele.ooEmooP = ( fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
    }


      ele.vtxFitConversion = ConversionTools::hasMatchedConversion(*el, conversions, beamSpot.position());

      // effective area for isolation
      ele.effArea = npknu::EffAreaElectron_Spring15_25ns(el->superCluster()->eta());

      // see TestElectronID/ElectronIDAnalyzer/plugins/ElectronIDAnalyzer.cc
      reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
      ele.pfIso_sumChargedHadronPt = pfIso.sumChargedHadronPt;
      ele.pfIso_sumNeutralHadronEt = pfIso.sumNeutralHadronEt;
      ele.pfIso_sumPhotonEt        = pfIso.sumPhotonEt       ;
      ele.pfIso_sumPUPt            = pfIso.sumPUPt           ;
      ele.absIsoWithDBeta          = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );

      ele.effAreaPFIso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rhoValue * npknu::EffAreaElectron_Spring15_25ns(el->superCluster()->eta()));

      // lazyTools  see 
      // lazyTools https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer.cc
      // EgammaAnalysis/ElectronTools/plugins/ElectronIDValueMapProducer.cc
      std::vector<float> vCov = lazyToolnoZS->localCovariances(*(el->superCluster()->seed()));
      float lazyToolnoZS_eleFull5x5SigmaIEtaIEta = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
      float lazyToolnoZS_eleFull5x5SigmaIEtaIPhi = vCov[1];
      float lazyToolnoZS_eleFull5x5R9            = lazyToolnoZS->e3x3(*(el->superCluster()->seed()) ) / el->superCluster()->rawEnergy() ;
      float lazyToolnoZS_e1x5                    = lazyToolnoZS->e1x5(*(el->superCluster()->seed()) );
      float lazyToolnoZS_e2x5Max                 = lazyToolnoZS->e2x5Max(*(el->superCluster()->seed()) );
      float lazyToolnoZS_e5x5                    = lazyToolnoZS->e5x5(*(el->superCluster()->seed()) );
      float lazyToolnoZS_circularity             = (lazyToolnoZS_e5x5 != 0.) ? 1.-lazyToolnoZS_e1x5/lazyToolnoZS_e5x5 : -1;
      float lazyToolnoZS_R_e1x5_e5x5             = (lazyToolnoZS_e5x5 != 0.) ? lazyToolnoZS_e1x5/lazyToolnoZS_e5x5 : 0;
      float lazyToolnoZS_R_e2x5_e5x5             = (lazyToolnoZS_e5x5 != 0.) ? lazyToolnoZS_e2x5Max/lazyToolnoZS_e5x5 : 0;

      ele.lazyToolnoZS_eleFull5x5SigmaIEtaIEta = lazyToolnoZS_eleFull5x5SigmaIEtaIEta ;
      ele.lazyToolnoZS_eleFull5x5SigmaIEtaIPhi = lazyToolnoZS_eleFull5x5SigmaIEtaIPhi ;
      ele.lazyToolnoZS_eleFull5x5R9            = lazyToolnoZS_eleFull5x5R9            ;
      ele.lazyToolnoZS_e1x5                    = lazyToolnoZS_e1x5                    ;
      ele.lazyToolnoZS_e2x5Max                 = lazyToolnoZS_e2x5Max                 ;
      ele.lazyToolnoZS_e5x5                    = lazyToolnoZS_e5x5                    ;
      ele.lazyToolnoZS_circularity             = lazyToolnoZS_circularity             ;
      ele.lazyToolnoZS_R_e1x5_e5x5             = lazyToolnoZS_R_e1x5_e5x5             ;
      ele.lazyToolnoZS_R_e2x5_e5x5             = lazyToolnoZS_R_e2x5_e5x5             ;

      // rho
      ele.rho = rhoValue;

      // Electrin ID
      //const edm::Ptr<pat::Electron> elPtr(electrons, el - electrons->begin());
      //ele.isCutPassTight  = (*cutIdDecisionsTight )[el];
      //ele.isCutPassMedium = (*cutIdDecisionsMedium)[el];
      //ele.isCutPassLoose  = (*cutIdDecisionsLoose )[el];
      //ele.isCutPassVeto   = (*cutIdDecisionsVeto  )[el];
      const bool passHEEPV70 = (*heep70Id )[el];

      ele.isCutPassHEEP   = passHEEPV70 ;

      float trkIsol   = (*trkIsolMap)[el];
      unsigned int heepV70Bitmap = (*heep70IdBitmap)[el];
      ele.heepV70Bitmap  = heepV70Bitmap ;
      ele.trkIsol         = trkIsol ;
      ele.matchToTruth = (doGenMatch) ? matchToTruth(el, genParticles) : -1 ;

      using HEEPV70 = VIDCutCodes<cutnrs::HEEPV70>; 
      const bool passEtShowerShapeHE = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ET, HEEPV70::SIGMAIETAIETA, HEEPV70::E2X5OVER5X5, HEEPV70::HADEM});
      const bool passN1TrkIso        = HEEPV70::pass(heepV70Bitmap,HEEPV70::TRKISO,HEEPV70::IGNORE);
      ele.passEtShowerShapeHE  = passEtShowerShapeHE ;  
      ele.passN1TrkIso         = passN1TrkIso        ;  
       //access # saturated crystals in the 5x5
      int nrSatCrys=(*nrSatCrysMap)[el];
      ele.nrSatCrys       = nrSatCrys ;
      if(nrSatCrys!=0) std::cout <<"###  nrSatCrys exist !!!  "<<nrSatCrys<<std::endl;

      const vid::CutFlowResult& heepCutFlowResult = (*heep70IdResult)[el];

      const bool passHEEPV70VID = heepCutFlowResult.cutFlowPassed();
      ele.passHEEPV70VID = passHEEPV70VID ;


      if(passHEEPV70!=passHEEPV70VID) std::cout <<"### error in HEEP VID " << passHEEPV70VID << ", HEEP ID result "<< passHEEPV70 << std::endl;
    

      //Fill MVA Info
      //ele.mvaIsPassMedium = (*mvaMediumIdDecisions)[el] ;
      //ele.mvaIsPassTight  = (*mvaTightIdDecisions)[el]  ;
      //ele.mvaValue        = (*mvaValues)[el] ;
      //ele.mvaCategory     = (*mvaCategories)[el] ;

    //how to get the track isolation from VID
    //note, this works for all cuts except E2x5/E5x5 although this is a feature which is
    //not often used so safer to use the standard accessors
    //trk isolation value is confirmed to be okay though
    const float trkIsoVID = heepCutFlowResult.getValueCutUpon(HEEPV70::TRKISO);
    ele.trkIsoVID = trkIsoVID ;
    if(trkIsol!=trkIsoVID) std::cout <<"### error in VID " << trkIsoVID << "  trk isol " << trkIsol << std::endl;

    const bool passEtShowerShapeHEVID = heepCutFlowResult.getCutResultByIndex(HEEPV70::ET)
      && heepCutFlowResult.getCutResultByIndex(HEEPV70::SIGMAIETAIETA) 
      && heepCutFlowResult.getCutResultByIndex(HEEPV70::E2X5OVER5X5)
      && heepCutFlowResult.getCutResultByIndex(HEEPV70::HADEM);
    //now for track isolation
    //now this may not be the fastest function as it makes a new cut flow...
   const bool passN1TrkIsoVID = heepCutFlowResult.getCutFlowResultMasking(HEEPV70::TRKISO).cutFlowPassed();
   ele.passEtShowerShapeHEVID = passEtShowerShapeHEVID ;
   ele.passN1TrkIsoVID        = passN1TrkIsoVID ;


 //   if(passN1TrkIso && !passHEEPV70) std::cout <<" trk isol "<<trkIso<<std::endl;

    if(passEtShowerShapeHE != passEtShowerShapeHEVID) std::cout <<"### error in VID showershape cuts"<<std::endl;
    if(passN1TrkIso != passN1TrkIsoVID) std::cout <<"### error in VID trk iso cuts"<<std::endl;

    bool VIDvsID = ((trkIsol!=trkIsoVID) || (passEtShowerShapeHE != passEtShowerShapeHEVID) || (passN1TrkIso != passN1TrkIsoVID)) ;
    ele.CheckVIDvsID = VIDvsID ;
	
      electronProd->push_back(ele);
   }
   iEvent.put(electronProd,"Electron");

   delete lazyToolnoZS;
}

//START From https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc
int ProdNpElectron::matchToTruth(const edm::Ptr<reco::GsfElectron> el, const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ProdNpElectron::findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus){
  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }
  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  return;
}

void ProdNpElectron::printCutFlowResult(vid::CutFlowResult &cutflow){
  printf("    CutFlow name= %s    decision is %d\n", cutflow.cutFlowName().c_str(), (int)cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut<ncuts; icut++){
		printf("  %2d      %50s    %d        %f          %d\n", icut,
	   	cutflow.getNameAtIndex(icut).c_str(),
		   (int)cutflow.isCutMasked(icut),
		   cutflow.getValueCutUpon(icut),
	   	(int)cutflow.getCutResultByIndex(icut));
  }
}

//void ProdNpElectron::CutFlowResult(vid::CutFlowResult &cutflow, std::vector<std::string>& nameVec, std::vector<float>& valueVec){
//	nameVec.clear();
//	valueVec.clear();
//  int ncuts = cutflow.cutFlowSize();
//  for(int icut = 0; icut<ncuts; icut++){
//			nameVec.push_back(cutflow.getNameAtIndex(icut).c_str());
//			valueVec.push_back(cutflow.getValueCutUpon(icut));
//  }
//}
// END From https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc


void ProdNpElectron::beginJob() { }
void ProdNpElectron::endJob() { }
void ProdNpElectron::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpElectron);
