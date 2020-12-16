#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"

#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"

class ProdNpPhoton : public edm::EDProducer {
   public:
      explicit ProdNpPhoton(const edm::ParameterSet&);
      ~ProdNpPhoton();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		enum PhotonMatchType {UNMATCHED = 0, 
						MATCHED_FROM_GUDSCB,
						MATCHED_FROM_PI0,
						MATCHED_FROM_OTHER_SOURCES};

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      bool hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol,
                           const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,
                           float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);

		int matchToTruth(const reco::Photon &pho, const edm::Handle<edm::View<reco::GenParticle> > &genParticles);
		void findFirstNonPhotonMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
		void printCutFlowResult(vid::CutFlowResult &cutflow); 

      edm::EDGetTokenT<edm::View<pat::Photon> > photonToken;
      edm::EDGetTokenT<double> rhoToken;
      edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken;
      edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken;
      edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken;
      edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken;

		edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken;
		bool DoGenMatch;

      edm::EDGetTokenT<edm::View<pat::Electron> > electronToken ;

      // ID decision objects
		edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenTight  ;
		edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenMedium ;
		edm::EDGetTokenT<edm::ValueMap<bool> > cutBasedIdTokenLoose  ;

      edm::EDGetTokenT<edm::ValueMap<bool> >               mvaMediumIdMapToken        ;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > mvaMediumIdFullInfoMapToken;
		edm::EDGetTokenT<edm::ValueMap<float> >              mvaValuesMapToken          ;
		edm::EDGetTokenT<edm::ValueMap<int> >                mvaCategoriesMapToken      ;

		double MinPtCut ;
		double MaxAEtaCut ;
		int LoopEvent;
};

ProdNpPhoton::ProdNpPhoton(const edm::ParameterSet& iConfig) {
   photonToken      = (consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))) ;
   rhoToken         = (consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))                     ;
   full5x5SigmaIEtaIEtaMapToken   = (consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap")))   ;
   phoChargedIsolationToken       = (consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoChargedIsolation")))       ;
   phoNeutralHadronIsolationToken = (consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))) ;
   phoPhotonIsolationToken        = (consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoPhotonIsolation")))        ;
   electronToken    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons")))       ;

	DoGenMatch = false;
	if(iConfig.existsAs<edm::InputTag>("genParticles")) {
		genParticlesToken = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
		DoGenMatch = true;
	}

   cutBasedIdTokenTight  = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenTight"))) ;
	cutBasedIdTokenMedium = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenMedium"))) ;
	cutBasedIdTokenLoose  = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("cutBasedIdTokenLoose"))) ;


   mvaMediumIdMapToken         = (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaMediumIdMap")))                       ;
   mvaMediumIdFullInfoMapToken = (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("mvaMediumIdFullInfoMap"))) ;
	mvaValuesMapToken           = (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap")))                        ;
	mvaCategoriesMapToken       = (consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))                      ;

	produces <std::vector<npknu::Photon> > ("Photon");

	MinPtCut = (iConfig.getUntrackedParameter<double>("MinPtCut",0.0));
	MaxAEtaCut = (iConfig.getUntrackedParameter<double>("MaxAEtaCut",100.0));
	std::cout << "ProdNtPhoton MinPt " << MinPtCut << " MaxAEta " << MaxAEtaCut << std::endl;

	LoopEvent=0;
}

ProdNpPhoton::~ProdNpPhoton() {}

void ProdNpPhoton::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
	using std::cout;
	using std::endl;

	LoopEvent++;
	bool IsMC = !iEvent.isRealData();

	std::auto_ptr<std::vector<npknu::Photon> > photonProd ( new std::vector<npknu::Photon>() );

   edm::Handle<edm::View<pat::Photon> > photons;   iEvent.getByToken(photonToken, photons);
   edm::Handle< double > rhoH;                     iEvent.getByToken(rhoToken,rhoH); float rho = *rhoH;
   edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap      ;     iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken, full5x5SigmaIEtaIEtaMap);
   edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap       ;     iEvent.getByToken(phoChargedIsolationToken, phoChargedIsolationMap);
   edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap ;     iEvent.getByToken(phoNeutralHadronIsolationToken, phoNeutralHadronIsolationMap);
   edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap        ;     iEvent.getByToken(phoPhotonIsolationToken, phoPhotonIsolationMap);

	edm::Handle<edm::View<reco::GenParticle> > genParticles         ;     if(IsMC and DoGenMatch) iEvent.getByToken(genParticlesToken,genParticles);

   edm::Handle<edm::View<pat::Electron> > electrons; iEvent.getByToken(electronToken, electrons);

   edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsTight ;
	edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsMedium;
	edm::Handle<edm::ValueMap<bool> >  cutIdDecisionsLoose ;
	iEvent.getByToken(cutBasedIdTokenTight  ,cutIdDecisionsTight );
	iEvent.getByToken(cutBasedIdTokenMedium ,cutIdDecisionsMedium);
	iEvent.getByToken(cutBasedIdTokenLoose  ,cutIdDecisionsLoose );

	edm::Handle<edm::ValueMap<bool> >               mvaMediumIdDecisions  ;
	edm::Handle<edm::ValueMap<vid::CutFlowResult> > mvaMediumIdCutflowData;
	edm::Handle<edm::ValueMap<float> >              mvaValues             ;
	edm::Handle<edm::ValueMap<int> >                mvaCategories         ;

	iEvent.getByToken(mvaMediumIdMapToken        ,mvaMediumIdDecisions  );
	iEvent.getByToken(mvaMediumIdFullInfoMapToken,mvaMediumIdCutflowData);
	iEvent.getByToken(mvaValuesMapToken          ,mvaValues             );
	iEvent.getByToken(mvaCategoriesMapToken      ,mvaCategories         );

   for(edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
      if(pho->pt() < MinPtCut) continue;
		npknu::Photon photon;

      const edm::Ptr<pat::Photon> phoPtr(photons, pho - photons->begin());
      float full5x5_sigmaIetaIetaMap = (*full5x5SigmaIEtaIEtaMap)[phoPtr]  ;
      float isoChargedHadronsMap = (*phoChargedIsolationMap)[phoPtr]       ;
      float isoNeutralHadronsMap = (*phoNeutralHadronIsolationMap)[phoPtr] ;
      float isoPhotonsMap        = (*phoPhotonIsolationMap)[phoPtr]        ;

      float isoChargedHadronsWithEA = std::max((float)0.0, isoChargedHadronsMap  - rho * npknu::EffAreaPhotonCha(pho->superCluster()->eta()) );
      float isoNeutralHadronsWithEA = std::max((float)0.0, isoNeutralHadronsMap  - rho * npknu::EffAreaPhotonNeu(pho->superCluster()->eta()) ) ;
      float isoPhotonsWithEA        = std::max((float)0.0, isoPhotonsMap         - rho * npknu::EffAreaPhotonPho(pho->superCluster()->eta()) ) ;
		bool passEleVeto = pho->passElectronVeto();
		bool hasPixelSeed = pho->hasPixelSeed();

      photon.pt                       = pho->pt()     ;
      photon.eta                      = pho->eta()    ;
      photon.phi                      = pho->phi()    ;
      photon.energy                   = pho->energy() ;
		photon.superCluster_eta         = pho->superCluster()->eta()  ;
		photon.superCluster_phi         = pho->superCluster()->phi()  ;
      photon.hadTowOverEm             = pho->hadTowOverEm()         ;
      photon.full5x5_sigmaIetaIeta    = pho->full5x5_sigmaIetaIeta();
      photon.full5x5_sigmaIetaIetaMap = full5x5_sigmaIetaIetaMap    ;
      photon.isoPhotonsMap               = isoPhotonsMap            ;
      photon.isoChargedHadronsMap        = isoChargedHadronsMap     ;
      photon.isoNeutralHadronsMap        = isoNeutralHadronsMap     ;
      photon.isoPhotonsWithEA         = isoPhotonsWithEA            ;
      photon.isoNeutralHadronsWithEA  = isoNeutralHadronsWithEA     ;
      photon.isoChargedHadronsWithEA  = isoChargedHadronsWithEA     ;
      photon.passEleVeto              = passEleVeto                 ;
		photon.hasPixelSeed             = hasPixelSeed                ;
		photon.matchToTruth             = (IsMC and DoGenMatch) ? matchToTruth(*pho, genParticles) : -1 ;  


		photon.isCutPassTight  = (*cutIdDecisionsTight)[phoPtr];
		photon.isCutPassMedium = (*cutIdDecisionsTight)[phoPtr];
		photon.isCutPassLoose  = (*cutIdDecisionsTight)[phoPtr];

		vid::CutFlowResult mvaMediumCutflowResult = (*mvaMediumIdCutflowData)[phoPtr];
		//printCutFlowResult(mvaMediumCutflowResult);
		photon.mvaIsMedium   = (*mvaMediumIdDecisions  )[phoPtr];
		photon.mvaValues     = (*mvaValues             )[phoPtr];
		photon.mvaCategories = (*mvaCategories         )[phoPtr];


		photonProd->push_back(photon);
   }
	iEvent.put(photonProd,"Photon");
}

bool ProdNpPhoton::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol, const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,  float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {
   if (sc.isNull()) return false;
   for (edm::View<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
      if (it->superCluster()!=sc) continue;
      if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
      if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
      return true;
   }
   return false;
}

// START From https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/PhotonNtupler/plugins/PhotonNtuplerVIDDemo.cc
int ProdNpPhoton::matchToTruth(const reco::Photon &pho, const edm::Handle<edm::View<reco::GenParticle> > &genParticles) {
  // 
  // Explicit loop and geometric matching method 
  //

  // Find the closest status 1 gen photon to the reco photon
  double dR = 999;
  const reco::Candidate *closestPhoton = 0;
  for(size_t i=0; i<genParticles->size();i++){
    const reco::Candidate *particle = &(*genParticles)[i];
    // Drop everything that is not photon or not status 1
    if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestPhoton = particle;
    }
  }
  // See if the closest photon (if it exists) is close enough.
  // If not, no match found.
  if( !(closestPhoton != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // Find ID of the parent of the found generator level photon match
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

  // Allowed parens: quarks pdgId 1-5, or a gluon 21
  std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
  if( !(std::find(allowedParents.begin(), 
		 allowedParents.end(), ancestorPID)
	!= allowedParents.end()) ){
    // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not. 
    if( abs(ancestorPID) == 111 )
      return MATCHED_FROM_PI0;
    else
      return MATCHED_FROM_OTHER_SOURCES;
  }
  return MATCHED_FROM_GUDSCB;
   
}

void ProdNpPhoton::findFirstNonPhotonMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus) {
  if( particle == 0 ){
    printf("PhotonNtuplerVIDDemo: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}

void ProdNpPhoton::printCutFlowResult(vid::CutFlowResult &cutflow) {
  printf("    CutFlow name= %s    decision is %d\n", cutflow.cutFlowName().c_str(),	 (int) cutflow.cutFlowPassed());
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut< (int)cutflow.cutFlowSize(); icut++){
    printf("  %d       %50s    %d        %f          %d\n", icut, cutflow.getNameAtIndex(icut).c_str(), (int)cutflow.isCutMasked(icut), cutflow.getValueCutUpon(icut), (int)cutflow.getCutResultByIndex(icut) );
  }
  
}

// END From https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/PhotonNtupler/plugins/PhotonNtuplerVIDDemo.cc


void ProdNpPhoton::beginJob() {}
void ProdNpPhoton::endJob() { }
void ProdNpPhoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpPhoton);
