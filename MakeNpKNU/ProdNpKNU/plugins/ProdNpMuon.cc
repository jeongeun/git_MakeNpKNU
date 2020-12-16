#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"


class ProdNpMuon : public edm::EDProducer {
   public:
      explicit ProdNpMuon(const edm::ParameterSet&);
      ~ProdNpMuon();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<edm::View<pat::Muon> > muonToken;
      double MinPtCut;
      
      edm::EDGetTokenT<reco::VertexCollection> vertexToken ;
    //  edm::EDGetTokenT<reco::Candidate> tokenPFCandidates ;

      edm::EDGetTokenT<edm::View<pat::PackedCandidate> >   tokenPFCandidates; 

};

ProdNpMuon::ProdNpMuon(const edm::ParameterSet& iConfig) {
   muonToken  = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("muons")))               ;
   MinPtCut = (iConfig.getUntrackedParameter<double>("MinPtCut",0.0));
   
   vertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexToken")) ;
   tokenPFCandidates= consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag> ("PFCandidates")  );

	produces <npknu::MuonCollection> ("Muon");
   std::cout << "ProdNpMuon MinPt " << MinPtCut << std::endl;
}

ProdNpMuon::~ProdNpMuon() { }
void ProdNpMuon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
   using namespace reco;
   using namespace pat;
   using namespace std;

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken, vertices);

   int idxGV = -1;
   int idx=0;
   for(const reco::Vertex& it : *vertices) {
      if(idxGV < 0) {
         if((it.ndof() > 4) and (std::fabs(it.z() <= 24)) and (std::fabs(it.position().rho()) <= 2)) idxGV = idx;
      }
      idx++;
   }
   idxGV = (idxGV < 0) ? 0 : idxGV;

   math::XYZPoint PV(vertices->at(0).x(), vertices->at(0).y(), vertices->at(0).z());
   math::XYZPoint GV(vertices->at(idxGV).x(), vertices->at(idxGV).y(), vertices->at(idxGV).z());

   std::auto_ptr<npknu::MuonCollection> muonProd ( new npknu::MuonCollection() );
   edm::Handle<edm::View<pat::Muon> > muons;
   iEvent.getByToken(muonToken, muons);

   edm::Handle<edm::View<pat::PackedCandidate> > pfCandidates;
   iEvent.getByToken(tokenPFCandidates, pfCandidates);

	for(edm::View<pat::Muon>::const_iterator mu = muons->begin(); mu!=muons->end(); ++mu) {
	    if(mu->pt() < MinPtCut) continue;
            //if(!(mu->isPFMuon()) || !(mu->isGlobalMuon()) ) continue; 

            using pat::Muon; // avoid hiding the base implementation
    	    npknu::Muon muon;
	    muon.isPFMuon = mu->isPFMuon() ;
	    muon.isGlobalMuon = mu->isGlobalMuon();
    	    muon.pt_P4    = mu->pt();
    	    muon.eta_P4   = mu->eta();
    	    muon.phi_P4   = mu->phi();
    	    muon.type     = mu->type() ;
    	    muon.energy   = mu->energy();

            //const reco::TrackRef& tunePMuonBestTk = mu->tunePMuonBestTrack() ;
            //const reco::TrackRef& muonBestTk      = mu->muonBestTrack()      ;
	    
            //const reco::TrackRef tpfmsTk         = mu->tpfmsTrack()         ;
            //const reco::TrackRef pickyTk         = mu->pickyTrack()         ;
            //const reco::TrackRef dytTk           = mu->dytTrack()           ;

	    if( mu->tunePMuonBestTrack().isNull() ){
	        muon.hastunePMuonBestTrack = false ;
	        muon.pt       =  -999.0 ;
	        muon.pt_TPq  =  -999.0 ;
	        muon.ptErr    =  -999.0 ;
	        muon.eta      =  -999.0 ;
	        muon.phi      =  -999.0 ;
                muon.charge   =  -999.0 ;
                muon.dxy      =  -999.0 ;
                muon.dz       =  -999.0 ;
                muon.dxy_GV   =  -999.0 ;
                muon.dz_GV    =  -999.0 ;
                muon.edxy     =  -999.0 ;
                muon.edz      =  -999.0 ;
            } else {
	        muon.hastunePMuonBestTrack = true ;
	        muon.pt       = mu->tunePMuonBestTrack()->pt()           ;
	        muon.pt_TPq  = mu->tunePMuonBestTrack()->pt() * mu->tunePMuonBestTrack()->charge()  ;
	        muon.ptErr    = mu->tunePMuonBestTrack()->ptError()      ;
	        muon.eta      = mu->tunePMuonBestTrack()->eta()          ;
	        muon.phi      = mu->tunePMuonBestTrack()->phi()          ;
                muon.charge   = mu->tunePMuonBestTrack()->charge()       ;
                muon.dxy      = mu->tunePMuonBestTrack()->dxy(vertices->front().position()  )  ;
                muon.dz       = mu->tunePMuonBestTrack()->dz( vertices->front().position()  )  ;
                muon.dxy_GV   = mu->tunePMuonBestTrack()->dxy(vertices->at(idxGV).position())  ;
                muon.dz_GV    = mu->tunePMuonBestTrack()->dz( vertices->at(idxGV).position())  ;
                muon.edxy     = mu->tunePMuonBestTrack()->dxyError()     ;
                muon.edz      = mu->tunePMuonBestTrack()->dzError()      ;
	    }

	   if( mu->muonBestTrack().isNull() ){
	       muon.hasmuonBestTrack     = false ;
               muon.pt_BTq     = -999.0 ;   
               muon.pt_BT       = -999.0 ;   
               muon.ptErr_BT    = -999.0 ;   
               muon.eta_BT      = -999.0 ;   
               muon.phi_BT      = -999.0 ;   
               muon.dxy_BT      =-999.0 ; 
               muon.dz_BT       =-999.0 ; 
               muon.dxy_GV_BT   =-999.0 ; 
               muon.dz_GV_BT    =-999.0 ; 
               muon.edxy_BT     =-999.0 ; 
               muon.edz_BT      =-999.0 ; 
               muon.momQuality  =-999.0 ; 
           } else {
               muon.hasmuonBestTrack = true ;
               muon.pt_BTq   = mu->muonBestTrack()->pt() * mu->muonBestTrack()->charge()  ;
               muon.pt_BT     = mu->muonBestTrack()->pt()  ;
               muon.ptErr_BT  = mu->muonBestTrack()->ptError()  ;
               muon.eta_BT    = mu->muonBestTrack()->eta()    ;
               muon.phi_BT    = mu->muonBestTrack()->phi()    ;
               muon.dxy_BT      = mu->muonBestTrack()->dxy(vertices->front().position()  )  ;
               muon.dz_BT       = mu->muonBestTrack()->dz( vertices->front().position()  )  ;
               muon.dxy_GV_BT   = mu->muonBestTrack()->dxy(vertices->at(idxGV).position())  ;
               muon.dz_GV_BT    = mu->muonBestTrack()->dz( vertices->at(idxGV).position())  ;
               muon.edxy_BT     = mu->muonBestTrack()->dxyError()     ;
               muon.edz_BT      = mu->muonBestTrack()->dzError()      ;
               muon.momQuality  = mu->muonBestTrack()->ptError()/mu->muonBestTrack()->pt() ;
           }
  
  			   //if( mu->tpfmsTrack().isNull() ) {
  				// muon.pt_TPFMS = -999.0 ;
  			   //} else{
  				// muon.pt_TPFMS = mu->tpfmsTrack()->pt() * mu->tpfmsTrack()->charge()  ;
            //}
			   //if( mu->pickyTrack().isNull() ) {
				// muon.pt_PICKY = -999.0 ;
			   //} else{
				// muon.pt_PICKY = mu->pickyTrack()->pt() * mu->pickyTrack()->charge()  ;
            //}
			   //if( mu->dytTrack().isNull() ) {
				// muon.pt_DYT = -999.0 ;
			   //} else{
				// muon.pt_DYT = mu->dytTrack()->pt() * mu->dytTrack()->charge()  ;
            //}

         //   muon.pt_TPFMS = (mu->tpfmsTrack().isNull()) ? -999.0 : mu->tpfmsTrack()->pt() * mu->tpfmsTrack()->charge()  ;    
         //   muon.pt_PICKY = (mu->pickyTrack().isNull()) ? -999.0 : mu->pickyTrack()->pt() * mu->pickyTrack()->charge()  ;   
         //   muon.pt_DYT   = (mu->dytTrack().isNull()  ) ? -999.0 : mu->dytTrack()->pt()   * mu->dytTrack()->charge()    ;     

        	//const reco::TrackRef globalTk  = mu->globalTrack()         ;
            if(mu->globalTrack().isNull())    {
	       muon.hasglobalTrack = false ;
               muon.pt_GB               =  -999.0 ;
               muon.pt_GBq              =  -999.0 ;
               muon.dxy_GB              =  -999.0 ;
               muon.dz_GB               =  -999.0 ;
               muon.dxy_GV_GB           =  -999.0 ;
               muon.dz_GV_GB            =  -999.0 ;
               muon.vx_GB               =  -999.0 ;
               muon.vy_GB               =  -999.0 ;
               muon.vz_GB               =  -999.0 ;
               muon.normalizedChi2_GB   =  -999.0 ;
               muon.nValidTrkHits       =  -999.0 ;
               muon.nValidMuonHits      =  -999.0 ;
            } else {
	       muon.hasglobalTrack = true ;
               muon.pt_GB                 = mu->globalTrack()->pt()   ;
               muon.pt_GBq                = mu->globalTrack()->pt() * mu->globalTrack()->charge()   ;
               muon.dxy_GB                = mu->globalTrack()->dxy(PV)                                    ;
               muon.dz_GB                 = mu->globalTrack()->dz(PV)                                     ;
               muon.dxy_GV_GB             = mu->globalTrack()->dxy(GV)                                    ;
               muon.dz_GV_GB              = mu->globalTrack()->dz(GV)                                     ;
               muon.vx_GB                 = mu->globalTrack()->vx()                                       ;
               muon.vy_GB                 = mu->globalTrack()->vy()                                       ;
               muon.vz_GB                 = mu->globalTrack()->vz()                                       ;
               muon.normalizedChi2_GB     = mu->globalTrack()->normalizedChi2()                           ;
               muon.nValidTrkHits         = mu->globalTrack()->hitPattern().numberOfValidTrackerHits()    ;
               muon.nValidMuonHits        = mu->globalTrack()->hitPattern().numberOfValidMuonHits()       ;
            }

            if(mu->innerTrack().isNull()) {
	       muon.hasinnerTrack = false;
               muon.pt_IN               =  -999.0  ;
               muon.pt_INq              =  -999.0  ;
               muon.TrkQuality_IN       =  -999.0  ;
               muon.dxy_IN              =  -999.0  ;
               muon.dz_IN               =  -999.0  ;
               muon.dxy_GV_IN           =  -999.0  ;
               muon.dz_GV_IN            =  -999.0  ;
               muon.ptErr_IN            =  -999.0  ;
               muon.normalizedChi2_IN   =  -999.0  ;
               muon.nValidPixelHits     =  -999.0  ;
               muon.nTrackerLayers      =  -999.0  ;
               muon.npixelLayers        =  -999.0  ;
               muon.momQuality_IN       =  -999.0  ;
               muon.Algo_IN             =  -999.0  ;
               muon.originalAlgo_IN     =  -999.0  ;
            } else {
	       muon.hasinnerTrack = true;
               muon.pt_IN                =  mu->innerTrack()->pt()                                        ;
               muon.pt_INq               =  mu->innerTrack()->pt()  * mu->innerTrack()->charge()          ;
               muon.TrkQuality_IN        =  mu->innerTrack()->quality(reco::TrackBase::highPurity)        ;
               muon.dxy_IN               =  mu->innerTrack()->dxy(vertices->front().position()  )         ;
               muon.dz_IN                =  mu->innerTrack()->dz( vertices->front().position()  )         ;
               muon.dxy_GV_IN            =  mu->innerTrack()->dxy(vertices->at(idxGV).position())         ;
               muon.dz_GV_IN             =  mu->innerTrack()->dz( vertices->at(idxGV).position())         ;
               muon.ptErr_IN             =  mu->innerTrack()->ptError()                                   ;
               muon.normalizedChi2_IN    =  mu->innerTrack()->normalizedChi2()                            ;
               muon.nValidPixelHits      =  mu->innerTrack()->hitPattern().numberOfValidPixelHits()       ;
               muon.nTrackerLayers       =  mu->innerTrack()->hitPattern().trackerLayersWithMeasurement() ;
               muon.npixelLayers         =  mu->innerTrack()->hitPattern().pixelLayersWithMeasurement()   ;
               muon.momQuality_IN        =  mu->innerTrack()->ptError()/mu->innerTrack()->pt() ;
               muon.Algo_IN              =  mu->innerTrack()->algo()  ;
               muon.originalAlgo_IN      =  mu->innerTrack()->originalAlgo()  ;  //Consider only muons from muonSeededStepOutIn algo
            }
            muon.nMatchedStations  = mu->numberOfMatchedStations() ;
            muon.nMatches          = mu->numberOfMatches()         ;
            muon.TrkIsoR03_sumPt   = mu->isolationR03().sumPt    ;
            muon.relTrkIsoR03_BT          = (mu->muonBestTrack()->pt() == 0.0 || !(muon.hasmuonBestTrack)) ? -999 : ((mu->isolationR03().sumPt)/(mu->muonBestTrack()->pt() ))    ;
            muon.relTrkIsoR03_TP          = (mu->tunePMuonBestTrack()->pt() == 0.0) ? -999 : ((mu->isolationR03().sumPt)/(mu->tunePMuonBestTrack()->pt() ))    ;
            muon.trackIso                 = mu->trackIso();
            muon.caloIso                  = mu->caloIso() ;
            muon.ecalIso                  = mu->ecalIso() ;
            muon.hcalIso                  = mu->hcalIso() ;
            muon.pfIsoR03_sumChargedHadronPt              = mu->pfIsolationR03().sumChargedHadronPt              ;
            muon.pfIsoR03_sumNeutralHadronEt              = mu->pfIsolationR03().sumNeutralHadronEt              ;
            muon.pfIsoR03_sumPhotonEt                     = mu->pfIsolationR03().sumPhotonEt                     ;
            muon.pfIsoR03_sumPUPt                         = mu->pfIsolationR03().sumPUPt                         ;
            muon.pfIsoR03_sumChargedParticlePt            = mu->pfIsolationR03().sumChargedParticlePt            ;
            muon.pfIsoR03_sumNeutralHadronEtHighThreshold = mu->pfIsolationR03().sumNeutralHadronEtHighThreshold ;
            muon.pfIsoR03_sumPhotonEtHighThreshold        = mu->pfIsolationR03().sumPhotonEtHighThreshold        ;
            muon.pfIsoR04_sumChargedHadronPt              = mu->pfIsolationR04().sumChargedHadronPt              ;
            muon.pfIsoR04_sumNeutralHadronEt              = mu->pfIsolationR04().sumNeutralHadronEt              ;
            muon.pfIsoR04_sumPhotonEt                     = mu->pfIsolationR04().sumPhotonEt                     ;
            muon.pfIsoR04_sumPUPt                         = mu->pfIsolationR04().sumPUPt                         ;
            muon.pfIsoR04_sumChargedParticlePt            = mu->pfIsolationR04().sumChargedParticlePt            ;
            muon.pfIsoR04_sumNeutralHadronEtHighThreshold = mu->pfIsolationR04().sumNeutralHadronEtHighThreshold ;
            muon.pfIsoR04_sumPhotonEtHighThreshold        = mu->pfIsolationR04().sumPhotonEtHighThreshold        ;

            const reco::MuonPFIsolation& pfiso04 = mu->pfIsolationR04();
            muon.pfIsoVar       = (mu->tunePMuonBestTrack()->pt() == 0.0) ? -999.0 : 
                               (pfiso04.sumChargedHadronPt + std::max(0.,pfiso04.sumNeutralHadronEt+pfiso04.sumPhotonEt-0.5*pfiso04.sumPUPt))/mu->tunePMuonBestTrack()->pt();

            muon.isLooseMuon    = mu->isLooseMuon()    ;
            muon.isTightMuonPV  = mu->isTightMuon(vertices->front())  ;
            muon.isSoftMuonPV   = mu->isSoftMuon(vertices->front())   ;
            muon.isHighPtMuonPV = mu->isHighPtMuon(vertices->front()) ;
            muon.isTightMuonGV  = mu->isTightMuon(vertices->at(idxGV))  ;
            muon.isSoftMuonGV   = mu->isSoftMuon(vertices->at(idxGV))   ;
            muon.isHighPtMuonGV = mu->isHighPtMuon(vertices->at(idxGV)) ;

            muon.dB_PV             = mu->dB();
            muon.edB_PV            = mu->edB();
            muon.dB_BS             = mu->dB(pat::Muon::BS2D);
            muon.edB_BS            = mu->edB(pat::Muon::BS2D);
            

       	   int xPV = (mu->isTightMuon(vertices->front())) ? 2 : 1;
       	   int yPV = (mu->isHighPtMuon(vertices->front())) ? 1 : 0;
       	   int zPV = (mu->isPFMuon()) ? 1 : 0;
       
       	   int xGV = (mu->isTightMuon( vertices->at(idxGV))) ? 2 : 1;
       	   int yGV = (mu->isHighPtMuon(vertices->at(idxGV))) ? 1 : 0;
       	   int zGV = (mu->isPFMuon()) ? 1 : 0;

            muon.statusPV  =  100*zPV+10*yPV+xPV   ; //isPF isHighPt isTight
            muon.statusGV  =  100*zGV+10*yGV+xGV   ; 

          for ( unsigned j=0; j < pfCandidates->size(); ++j ) {
               const pat::PackedCandidate & pfCandidate = (*pfCandidates)[j];
               // look for pf muon
               bool isPFMuon = ( (abs(pfCandidate.pdgId()) == 13) and (pfCandidate.pt() > 100.0) );
               if( isPFMuon ){ 
                    float dr = deltaR( muon.eta, muon.phi, pfCandidate.eta(), pfCandidate.phi() );
                    muon.dr_muPf = dr ;
                    if(dr < 0.001){
                       bool isBadPFMuon = true ;
                       muon.isBadPFMuon = isBadPFMuon ;
                    }
                 }
               
               bool isChHad  = (abs(pfCandidate.pdgId()) == 211); 
               bool isInnerTrack = (not (mu->innerTrack().isNull()) ) ;
               // require small DR
               if( isChHad and isInnerTrack){ 
                   float dr_ch = deltaR( mu->innerTrack()->eta(), mu->innerTrack()->phi(), pfCandidate.eta(), pfCandidate.phi());
                   float dpt = ( pfCandidate.pt() -(mu->innerTrack()->pt()))/(0.5*((mu->innerTrack()->pt()) + pfCandidate.pt()));
                    muon.dr_chPf = dr_ch ;
                    muon.dpt_chPf = dpt ;
                    if( (dr_ch < 0.00001) and (dpt < 0.00001) ){
                       bool isBadChHadron = true ;
                       muon.isBadChHadron = isBadChHadron ;
                    }
                 }
              }



		muonProd->push_back(muon);
	}
	iEvent.put(muonProd,"Muon");
}

void ProdNpMuon::beginJob() { }
void ProdNpMuon::endJob() { } 
void ProdNpMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpMuon);
