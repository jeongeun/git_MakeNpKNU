#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "FWCore/Utilities/interface/InputTag.h"


class ProdNpMET : public edm::EDProducer {
   public:
      explicit ProdNpMET(const edm::ParameterSet&);
      ~ProdNpMET();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

		edm::EDGetTokenT<pat::METCollection> metToken              ;
                edm::EDGetTokenT<pat::METCollection> metEGCleanToken       ;
                edm::EDGetTokenT<pat::METCollection> metMuEGCleanToken     ;
       //         edm::EDGetTokenT<pat::METCollection> metMuEGCleanCorrToken ;
                edm::EDGetTokenT<pat::METCollection> metUncorrectedToken   ;
      		edm::EDGetTokenT<pat::METCollection> metPUPPIToken         ;
		edm::EDGetTokenT<bool> BadChCandFilterToken;
		edm::EDGetTokenT<bool> BadPFMuonFilterToken;
       //         edm::EDGetTokenT<edm::PtrVector<reco::Muon> > badGlobalMuonFilterToken;
       //         edm::EDGetTokenT<edm::PtrVector<reco::Muon> > duplicateMuonFilterToken;      
		//edm::EDGetTokenT<bool> BadGlobalMuonTaggerFilterToken;      

		edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken;
		edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_data;
                edm::EDGetTokenT<bool> dupECALClustersToken ; //  particleFlowEGammaGSFixed;                
                edm::EDGetTokenT<edm::EDCollection<DetId>> hitsNotReplacedToken ; //ecalMultiAndGSGlobalRecHitEB;

	//	edm::EDGetTokenT<bool> hbheNoiseFilterToken_;
	//	edm::EDGetTokenT<bool> hbheTightNoiseFilterToken_;
	//	edm::EDGetTokenT<bool> hbheIsoNoiseFilterToken_;

		double MinPtCut ;
		bool PrintNum   ;
                bool IsMC = false;
		int edmNum      ;
		int prdNum     ;
		//int prdNum1     ;
		int prdNum2     ;
		int prdNum3     ;
		int prdNum4     ;
		int prdNum5     ;
		int prdNum6     ;
		int prdNum7     ;
		int prdNum8     ;

};

ProdNpMET::ProdNpMET(const edm::ParameterSet& iConfig) {
	metToken             = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metToken")) ;
        //Re-miniaod data (2017/Feb/03)
        metEGCleanToken      = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metEGCleanToken"))  ; 
        metMuEGCleanToken    = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metMuEGCleanToken"))  ; 
   //     metMuEGCleanCorrToken= consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metMuEGCleanCorrToken"))  ; 
        metUncorrectedToken  = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metUncorrectedToken"))  ; 
	metPUPPIToken        = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metPUPPIToken")) ;
	BadPFMuonFilterToken = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter")); 
	BadChCandFilterToken = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter")); 
	//BadGlobalMuonTaggerFilterToken = consumes<bool>(iConfig.getParameter<edm::InputTag>("badGlobalMuonTagger")); 

        //badGlobalMuonFilterToken    = consumes<edm::PtrVector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("badGlobalMuonFilterToken")) ;
        //duplicateMuonFilterToken    = consumes<edm::PtrVector<reco::Muon> >(iConfig.getParameter<edm::InputTag>("duplicateMuonFilterToken")) ;      
        metFilterBitsToken = mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"));
        metFilterBitsToken_data = mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits_data"));
        //Re-miniaod data (2017/Feb/03)
        dupECALClustersToken = consumes<bool>(iConfig.getParameter<edm::InputTag>("particleFlowEGammaGSFixed")) ;
        hitsNotReplacedToken = consumes<edm::EDCollection<DetId> >(iConfig.getParameter<edm::InputTag>("ecalMultiAndGSGlobalRecHitEB")) ;

	produces <std::vector<npknu::MET> > ("MET");

	MinPtCut = (iConfig.getUntrackedParameter<double>("MinPtCut",0.0));
	PrintNum = (iConfig.getUntrackedParameter<bool>("PrintNum",false));
	edmNum = 0;
	prdNum = 0;
	prdNum2 = 0;
	prdNum3 = 0;
	prdNum4 = 0;
	prdNum5 = 0;
	prdNum6 = 0;
	prdNum7 = 0;
	prdNum8 = 0;
	std::cout << "ProdNtMET MinPt " << MinPtCut << " Print " << PrintNum << std::endl;
	

}
ProdNpMET::~ProdNpMET() {}

void ProdNpMET::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

	bool IsMC = !iEvent.isRealData();

	std::auto_ptr<std::vector<npknu::MET> > metProd ( new std::vector<npknu::MET>() );
	//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/MET.h
	// double shiftedPt(METUncertainty shift, METCorrectionLevel level=Type1) const { return shiftedP2(shift,level).pt(); };

	edm::Handle<pat::METCollection> mets ;
        edm::Handle<pat::METCollection> metsEGClean;
        edm::Handle<pat::METCollection> metsMuEGClean;
        edm::Handle<pat::METCollection> metsUncorr;

	iEvent.getByToken(metToken, mets);
	const pat::MET &it = mets->front();

	iEvent.getByToken(metEGCleanToken, metsEGClean);
	const pat::MET &itEG = metsEGClean->front();

	iEvent.getByToken(metMuEGCleanToken, metsMuEGClean);
	const pat::MET &itMuEG = metsMuEGClean->front();

	iEvent.getByToken(metUncorrectedToken, metsUncorr);
	const pat::MET &itUncorr = metsUncorr->front();

        edm::Handle<pat::METCollection> metPUPPI;
        iEvent.getByToken(metPUPPIToken, metPUPPI);
	const pat::MET &itpuppi = metPUPPI->front();        
      
     //   edm::Handle<edm::PtrVector<reco::Muon> > badGlobalMuonFilter;
     //   edm::Handle<edm::PtrVector<reco::Muon> > duplicateMuonFilter;


	edmNum = mets->size(); //prdNum=0;prdNum2=0;prdNum3=0;prdNum4=0;prdNum5=0;
           npknu::MET met;
           met.isCaloMET           = it.isCaloMET() ;
           met.isPFMET             = it.isPFMET()   ;
           met.isRecoMET           = it.isRecoMET() ;

           met.pt                  = it.pt()                           ;
           met.phi                 = it.phi()                          ;
           met.sumEt               = it.sumEt()                        ;
           met.metSignificance     = it.metSignificance() ;

           met.EGClean_pt         = (IsMC) ? -999.0 : itEG.pt()                           ;
           met.EGClean_phi        = (IsMC) ? -999.0 : itEG.phi()                          ;
           met.EGClean_sumEt      = (IsMC) ? -999.0 : itEG.sumEt()                        ;
           met.EGClean_metSig     = (IsMC) ? -999.0 : itEG.metSignificance() ;

           met.MuEGClean_pt       = (IsMC) ? -999.0 : itMuEG.pt()                           ;
           met.MuEGClean_phi      = (IsMC) ? -999.0 : itMuEG.phi()                          ;
           met.MuEGClean_sumEt    = (IsMC) ? -999.0 : itMuEG.sumEt()                        ;
           met.MuEGClean_metSig   = (IsMC) ? -999.0 : itMuEG.metSignificance() ;

           met.Uncorr_pt          = (IsMC) ? -999.0 : itUncorr.pt()                           ;
           met.Uncorr_phi         = (IsMC) ? -999.0 : itUncorr.phi()                          ;
           met.Uncorr_sumEt       = (IsMC) ? -999.0 : itUncorr.sumEt()                        ;
           met.Uncorr_metSig      = (IsMC) ? -999.0 : itUncorr.metSignificance() ;

           met.genMET_pt           = (IsMC) ? it.genMET()->pt()    : -999.0 ;
           met.genMET_phi          = (IsMC) ? it.genMET()->phi()   : -999.0 ;
           met.genMET_sumEt        = (IsMC) ? it.genMET()->sumEt() : -999.0 ;

           met.corPt_Type01          = it.corPt(pat::MET::METCorrectionLevel::Type01 )  ;
           met.corPhi_Type01         = it.corPhi(pat::MET::METCorrectionLevel::Type01 )  ;
           met.corSumEt_Type01       = it.corSumEt(pat::MET::METCorrectionLevel::Type01 )  ;
           met.corPt_TypeXY          = it.corPt(pat::MET::METCorrectionLevel::TypeXY )  ;
           met.corPhi_TypeXY         = it.corPhi(pat::MET::METCorrectionLevel::TypeXY )  ;
           met.corSumEt_TypeXY       = it.corSumEt(pat::MET::METCorrectionLevel::TypeXY )  ;
           met.corPt_Type1XY         = it.corPt(pat::MET::METCorrectionLevel::Type1XY )  ;
           met.corPhi_Type1XY        = it.corPhi(pat::MET::METCorrectionLevel::Type1XY )  ;
           met.corSumEt_Type1XY      = it.corSumEt(pat::MET::METCorrectionLevel::Type1XY )  ;
           met.corPt_Type01XY        = it.corPt(pat::MET::METCorrectionLevel::Type01XY )  ;
           met.corPhi_Type01XY       = it.corPhi(pat::MET::METCorrectionLevel::Type01XY )  ;
           met.corSumEt_Type01XY     = it.corSumEt(pat::MET::METCorrectionLevel::Type01XY )  ;
           met.corPt_Type1Smear      = it.corPt(pat::MET::METCorrectionLevel::Type1Smear )  ;
           met.corPhi_Type1Smear     = it.corPhi(pat::MET::METCorrectionLevel::Type1Smear )  ;
           met.corSumEt_Type1Smear   = it.corSumEt(pat::MET::METCorrectionLevel::Type1Smear )  ;
           met.corPt_Type01Smear     = it.corPt(pat::MET::METCorrectionLevel::Type01Smear )  ;
           met.corPhi_Type01Smear    = it.corPhi(pat::MET::METCorrectionLevel::Type01Smear )  ;         
           met.corSumEt_Type01Smear  = it.corSumEt(pat::MET::METCorrectionLevel::Type01Smear )  ;        
           met.corPt_Type1SmearXY    = it.corPt(pat::MET::METCorrectionLevel::Type1SmearXY)  ;          
           met.corPhi_Type1SmearXY   = it.corPhi(pat::MET::METCorrectionLevel::Type1SmearXY)  ;          
           met.corSumEt_Type1SmearXY = it.corSumEt(pat::MET::METCorrectionLevel::Type1SmearXY)  ;           
           met.corPt_Type01SmearXY   = it.corPt(pat::MET::METCorrectionLevel::Type01SmearXY)  ;
           met.corPhi_Type01SmearXY  = it.corPhi(pat::MET::METCorrectionLevel::Type01SmearXY)  ;
           met.corSumEt_Type01SmearXY= it.corSumEt(pat::MET::METCorrectionLevel::Type01SmearXY)  ;      
           
           met.shiftedPt_Type1XY           = it.shiftedPt(pat::MET::NoShift,pat::MET::Type1XY) ;
           met.shiftedPhi_Type1XY           = it.shiftedPhi(pat::MET::NoShift,pat::MET::Type1XY) ;

           met.shiftedPt_JetResUp           =  it.shiftedPt(pat::MET::JetResUp          );
           met.shiftedPt_JetResDown         =  it.shiftedPt(pat::MET::JetResDown        );
           met.shiftedPt_JetEnUp            =  it.shiftedPt(pat::MET::JetEnUp           );
           met.shiftedPt_JetEnDown          =  it.shiftedPt(pat::MET::JetEnDown         );
           //met.shiftedPt_MuonEnUp           = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::MuonEnUp          );
           //met.shiftedPt_MuonEnDown         = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::MuonEnDown        );
           //met.shiftedPt_ElectronEnUp       = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::ElectronEnUp      );
           //met.shiftedPt_ElectronEnDown     = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::ElectronEnDown    );
           //met.shiftedPt_TauEnUp	    = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::TauEnUp           );
           //met.shiftedPt_TauEnDown          = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::TauEnDown         );
           met.shiftedPt_UnclusteredEnUp    =  it.shiftedPt(pat::MET::UnclusteredEnUp   );
           met.shiftedPt_UnclusteredEnDown  =  it.shiftedPt(pat::MET::UnclusteredEnDown );
           //met.shiftedPt_PhotonEnUp         = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::PhotonEnUp        );
           //met.shiftedPt_PhotonEnDown       = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::PhotonEnDown      );
           //met.shiftedPt_JetResUpSmear      = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::JetResUpSmear     );
           //met.shiftedPt_JetResDownSmear    = (!IsMC) ? -999.0 : it.shiftedPt(pat::MET::JetResDownSmear   );

           met.METPuppi_pt         = itpuppi.pt();
           met.METPuppi_phi        = itpuppi.phi();

           met.shiftedPt_Type1XY_puppi        = itpuppi.shiftedPt(pat::MET::NoShift,pat::MET::Type1XY) ;
           met.shiftedPhi_Type1XY_puppi       = itpuppi.shiftedPhi(pat::MET::NoShift,pat::MET::Type1XY) ;
           met.METPuppi_corPt_Type01          = itpuppi.corPt(pat::MET::METCorrectionLevel::Type01 )  ;
           met.METPuppi_corPhi_Type01         = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type01 )  ;
           met.METPuppi_corSumEt_Type01       = itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type01 )  ;
           met.METPuppi_corPt_TypeXY          = itpuppi.corPt(pat::MET::METCorrectionLevel::TypeXY )  ;
           met.METPuppi_corPhi_TypeXY         = itpuppi.corPhi(pat::MET::METCorrectionLevel::TypeXY )  ;
           met.METPuppi_corSumEt_TypeXY       = itpuppi.corSumEt(pat::MET::METCorrectionLevel::TypeXY )  ;
           met.METPuppi_corPt_Type1XY         = itpuppi.corPt(pat::MET::METCorrectionLevel::Type1XY )  ;
           met.METPuppi_corPhi_Type1XY        = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type1XY )  ;
           met.METPuppi_corSumEt_Type1XY      = itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type1XY )  ;
           met.METPuppi_corPt_Type01XY        = itpuppi.corPt(pat::MET::METCorrectionLevel::Type01XY )  ;
           met.METPuppi_corPhi_Type01XY       = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type01XY )  ;
           met.METPuppi_corSumEt_Type01XY     = itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type01XY )  ;
           met.METPuppi_corPt_Type1Smear      = itpuppi.corPt(pat::MET::METCorrectionLevel::Type1Smear )  ;
           met.METPuppi_corPhi_Type1Smear     = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type1Smear )  ;
           met.METPuppi_corSumEt_Type1Smear   = itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type1Smear )  ;
           met.METPuppi_corPt_Type01Smear     = itpuppi.corPt(pat::MET::METCorrectionLevel::Type01Smear )  ;
           met.METPuppi_corPhi_Type01Smear    = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type01Smear )  ;         
           met.METPuppi_corSumEt_Type01Smear  = itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type01Smear )  ;        
           met.METPuppi_corPt_Type1SmearXY    = itpuppi.corPt(pat::MET::METCorrectionLevel::Type1SmearXY)  ;          
           met.METPuppi_corPhi_Type1SmearXY   = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type1SmearXY)  ;          
           met.METPuppi_corSumEt_Type1SmearXY = itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type1SmearXY)  ;           
           met.METPuppi_corPt_Type01SmearXY   = itpuppi.corPt(pat::MET::METCorrectionLevel::Type01SmearXY)  ;
           met.METPuppi_corPhi_Type01SmearXY  = itpuppi.corPhi(pat::MET::METCorrectionLevel::Type01SmearXY)  ;
           met.METPuppi_corSumEt_Type01SmearXY= itpuppi.corSumEt(pat::MET::METCorrectionLevel::Type01SmearXY)  ;      
     edm::Handle<edm::TriggerResults> metFilterBits;
     iEvent.getByToken(metFilterBitsToken, metFilterBits);
     if( !metFilterBits.isValid() ){
         iEvent.getByToken(metFilterBitsToken_data, metFilterBits);
     }
     const edm::TriggerNames &metFilterNames = iEvent.triggerNames(*metFilterBits);

    for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i){

      if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_trackingFailureFilter") == 0)
	met.Flag_trackingFailureFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0)
        met.Flag_HBHENoiseFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0)
         met.Flag_HBHENoiseIsoFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_CSCTightHaloFilter") == 0)
	met.Flag_CSCTightHaloFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_CSCTightHaloTrkMuUnvetoFilter") == 0)
	met.Flag_CSCTightHaloTrkMuUnvetoFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_CSCTightHalo2015Filter") == 0)
	met.Flag_CSCTightHalo2015Filter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)
	met.Flag_goodVertices = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_globalTightHalo2016Filter") == 0)
	met.Flag_globalTightHalo2016Filter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_globalSuperTightHalo2016Filter") == 0)
	met.Flag_globalSuperTightHalo2016Filter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_HcalStripHaloFilter") == 0)
	met.Flag_HcalStripHaloFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_hcalLaserEventFilter") == 0)
	met.Flag_hcalLaserEventFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0)
	met.Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_EcalDeadCellBoundaryEnergyFilter") == 0)
	met.Flag_EcalDeadCellBoundaryEnergyFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0)
	met.Flag_eeBadScFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_ecalLaserCorrFilter") == 0)
	met.Flag_ecalLaserCorrFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_trkPOGFilters") == 0)
	met.Flag_trkPOGFilters = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_chargedHadronTrackResolutionFilter") == 0)
	met.Flag_chargedHadronTrackResolutionFilter = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_muonBadTrackFilter") == 0)
	met.Flag_muonBadTrackFilter = metFilterBits->accept(i);
      //# and the sub-filters
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_trkPOG_manystripclus53X") == 0)
	met.Flag_trkPOG_manystripclus53X = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_trkPOG_toomanystripclus53X") == 0)
	met.Flag_trkPOG_toomanystripclus53X = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_trkPOG_logErrorTooManyClusters") == 0)
	met.Flag_trkPOG_logErrorTooManyClusters = metFilterBits->accept(i);
      //# and the summary = HBHENoise primaryVertex CSCTightHalo eeBadSc chargedHadronTrackResol muonBadTrack
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_METFilters") == 0)
	met.Flag_METFilters = metFilterBits->accept(i);
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_badMuons") == 0) {
	met.Flag_badGlobalMuonFilter = metFilterBits->accept(i);  
        met.IsCorrectedBadGlobalMuon = metFilterBits->accept(i);
	//std::cout << "found bad muon flag : " << std::endl;
      }
      else if(strcmp(metFilterNames.triggerName(i).c_str(), "Flag_duplicateMuons") == 0)
	met.Flag_duplicateMuonFilter = metFilterBits->accept(i);  
        met.IsCorrectedDuplicateMuon = metFilterBits->accept(i);
    } //loop over met filters

	//https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
        edm::Handle<bool> BadChargedCandidateFilter;
        iEvent.getByToken(BadChCandFilterToken, BadChargedCandidateFilter);

        edm::Handle<bool> BadPFMuonFilter;
        iEvent.getByToken(BadPFMuonFilterToken, BadPFMuonFilter);

        //edm::Handle<bool> ifilterbadGlobalMuonTagger;
        //iEvent.getByToken(BadGlobalMuonTaggerFilterToken, ifilterbadGlobalMuonTagger);
	//met.filterbadGlobalMuonTagger = *ifilterbadGlobalMuonTagger;

        met.IsNoBadPFMuon           = *BadPFMuonFilter         ;
        met.IsNoBadChargedCandidate = *BadChargedCandidateFilter;


        met.Flag_BadPFMuonFilter           = *BadPFMuonFilter;
        met.Flag_BadChargedCandidateFilter = *BadChargedCandidateFilter;

	if (*BadPFMuonFilter           == true) { prdNum2++;  }
        if (*BadChargedCandidateFilter == true) { prdNum3++;  }
	if ( met.IsCorrectedBadGlobalMuon == true) { prdNum4++;  }
        if ( met.IsCorrectedDuplicateMuon == true) { prdNum5++;  }


        //ECAL Slew correction check
        edm::Handle<bool> dupECALClusters ;           
        iEvent.getByToken(dupECALClustersToken, dupECALClusters);
        const bool IsdupECALCluster = (dupECALClusters.isValid() and *(dupECALClusters) == true);
        met.IsDuplicateECALClusters = IsdupECALCluster;
        

        edm::Handle<edm::EDCollection<DetId> > hitsNotReplaced ;
        iEvent.getByToken(hitsNotReplacedToken,hitsNotReplaced);
        const bool IshitNotReplaced = ( /*hitsNotReplaced.isValid() and */ !(hitsNotReplaced->empty())); 
        met.IsNotReplatedHitsECAL = IshitNotReplaced ;
        
        const bool IsBadECALSlew = (met.IsDuplicateECALClusters ==true or met.IsNotReplatedHitsECAL == true);
        met.IsBadECALSlew = IsBadECALSlew ;
        //std::cout << "isDuplicateECAL " << met.IsDuplicateECALClusters << " notRelatedHit " << met.IsNotReplatedHitsECAL << " BadECALSlew " << met.IsBadECALSlew << std::endl;
        if(met.IsDuplicateECALClusters == true){ prdNum6++;} 
        if(met.IsNotReplatedHitsECAL == true){ prdNum7++;} 
        if(met.IsBadECALSlew    == true){ prdNum8++;} 


	metProd->push_back(met);
		prdNum++;
	iEvent.put(metProd, "MET");
        //if(PrintNum) { std::cout << "ProdNpMET NumberOf METs =" << edmNum   <<  "=" << prdNum << std::endl;
        //               std::cout << "-> GoodMET isNoBadPFMu  =" << prdNum2  << " BadPFMu  " << prdNum - prdNum2 <<std::endl;
        //               std::cout << "-> GoodMET isNoBadChCan =" << prdNum3  << " BadChCan " << prdNum - prdNum3 <<std::endl; 
        //               std::cout << "-> CorrectedMET BadGBMu =" << prdNum4  << std::endl;
        //               std::cout << "-> CorrectedMET DupliMu =" << prdNum5  << std::endl; 
        //               std::cout << "Re-miniaod NumberOf Bad ECAL Slew Events " << std::endl; 
        //               std::cout << "-> Bad hasDuplicateECALClusters =" << prdNum6  << std::endl; 
        //               std::cout << "-> Bad NotReplacedHitsECAL      =" << prdNum7  << std::endl; 
        //               std::cout << "-> BadECALSlew =" << prdNum8  << std::endl; 
        //  }


// Example (metFilters_cff.py)
//    metFilters = cms.Sequence(
//  >>     HBHENoiseFilterResultProducer *
//  >>     HBHENoiseFilter *
//  >>     primaryVertexFilter*
// 	   #   HBHENoiseIsoFilter*
//   	 #   HcalStripHaloFilter *
//  >>     CSCTightHaloFilter *
//	    #   hcalLaserEventFilter *
//      	 #Various proposals for updated halo filters.
//	       ##2015 proposals: 
//    		   #CSCTightHaloTrkMuUnvetoFilter *
//	    	   #CSCTightHalo2015Filter *
//     		  ##2016 proposals
//	     	  #globalTightHalo2016Filter*
//     		  #globalSuperTightHalo2016Filter*
//    >>   EcalDeadCellTriggerPrimitiveFilter* 
//	    #   *goodVertices * trackingFailureFilter *
//    >>   eeBadScFilter*
//	    #   ecalLaserCorrFilter *
//	    #   trkPOGFilters
//    >>   chargedHadronTrackResolutionFilter *
//    >>   muonBadTrackFilter
//    )


}

void ProdNpMET::beginJob() {}
void ProdNpMET::endJob() { 
   if(PrintNum) { std::cout << "ProdNpMET NumberOf METs =" << edmNum   <<  "=" << prdNum << std::endl;
                  std::cout << "-> GoodMET isNoBadPFMu  =" << prdNum2  << " BadPFMu  " << prdNum - prdNum2 <<std::endl;
                  std::cout << "-> GoodMET isNoBadChCan =" << prdNum3  << " BadChCan " << prdNum - prdNum3 <<std::endl; 
                  std::cout << "-> CorrectedMET BadGBMu =" << prdNum4  << std::endl;
                  std::cout << "-> CorrectedMET DupliMu =" << prdNum5  << std::endl; 
                  std::cout << "Re-miniaod NumberOf Bad ECAL Slew Events " << std::endl; 
                  std::cout << "-> Bad hasDuplicateECALClusters =" << prdNum6  << std::endl; 
                  std::cout << "-> Bad NotReplacedHitsECAL      =" << prdNum7  << std::endl; 
                  std::cout << "-> BadECALSlew =" << prdNum8  << std::endl; 
     }
 }
void ProdNpMET::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ProdNpMET);
