#ifndef NpKNU_HH
#define NpKNU_HH

#include <iostream>
#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>
#include "TClonesArray.h"
#include <string>
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

namespace npknu {

int CompareWild(const char *string, const char *wild);
TString AddStringToFileName(TString _orig, TString _addString) ;

float EffAreaElectron_Spring15_25ns(double scEta);
float EffAreaElectron_Spring15_50ns(double scEta);
float EffAreaPhotonCha(double scEta, bool is50ns=false);
float EffAreaPhotonNeu(double scEta, bool is50ns=false);
float EffAreaPhotonPho(double scEta, bool is50ns=false);


class Evt : public TObject {
public:
   int run    ; 
   int lumi   ; 
   double event  ; 
	bool isRealData ; 
   Evt() {};
   Evt(int _run, int _lumi, double _event, bool _isRealData) : run(_run), lumi(_lumi), event(_event), isRealData(_isRealData) {};
   ~Evt() {}; 
	friend std::ostream& operator<<(std::ostream& os, const Evt& evt);
   ClassDef(Evt,1);
};

class P4 : public TObject {
public:
	double pt;
	double eta;
	double phi;
	double energy;
	P4() {};
	 ~P4() {};
	double deltaPhi(const P4* ptr);
	double deltaPhi(const P4 ptr);
	double deltaPhi(double PHI);
	double etRatio(const P4* ptr);
	double etRatio(const P4 obj);
	double etRatio(double PT);
	double deltaR(const P4* ptr);
	double deltaR(const P4 obj);
   TLorentzVector P4Vector();
   TLorentzVector P4TVector();
	double GetMt(const P4* ptr);
	double GetMt(const P4 ptr);
	double GetMt(double PT, double PHI);
	double GetM(const P4 ptr);
	double GetM(const P4* ptr);
	void PrintP4();
	friend std::ostream& operator<<(std::ostream& os, const P4& obj);
	ClassDef(P4,1);

};


class GenParticle : public P4 {
public:
	int    status ;
	int    pdgId  ;
	double mass   ;
   unsigned int nDau  ; 
   int mIdx ;

   GenParticle() {};
    ~GenParticle() {}; 
   ClassDef(GenParticle,1);
};

class Photon : public P4 {
public:
	double superCluster_eta         ;
	double superCluster_phi         ;
	double hadTowOverEm             ;
	double full5x5_sigmaIetaIeta    ;
	double full5x5_sigmaIetaIetaMap ;
	double isoChargedHadronsMap     ;
	double isoNeutralHadronsMap     ;
	double isoPhotonsMap            ;
	double isoPhotonsWithEA         ;
	double isoNeutralHadronsWithEA  ;
	double isoChargedHadronsWithEA  ;
	bool   passEleVeto              ;
	bool   hasPixelSeed             ;
	int    matchToTruth             ;
//	unsigned int idBits             ;
//	bool passTightCutId             ;
	std::vector<bool> isID          ; 
	std::vector<std::vector<float> > idValuesVec;

   bool isCutPassTight  ;
   bool isCutPassMedium ;
   bool isCutPassLoose  ;

   bool   mvaIsMedium              ; 
   float  mvaValues                ;
   int    mvaCategories            ;

	std::vector<float> GetIdValues(int idIndex) { return idValuesVec.at(idIndex);	};
	void PrintIdValues(int idIndex) {
		std::vector<float> thisVec = GetIdValues(idIndex);
		int thisI = 0;
		for(std::vector<float>::iterator it=thisVec.begin(); it!=thisVec.end(); it++, thisI++) {
			std::cout << "PhoIdValue " << idIndex << " " << thisI << " " << *it << std::endl;
		}
	};	
	double scEta() { return superCluster_eta; };
	double scAEta() { return std::abs(superCluster_eta); };
	bool   isGap()  { return ((scAEta() > 1.4442) and (scAEta() < 1.566)) ; };
	bool   isEB()   { return ( scAEta() <= 1.479) ; };
	bool   isEE()   { return ((scAEta() > 1.479 ) and (scAEta() < 2.5)); };

	float CutThresholer(std::string cutName, int idIndex=2, bool is50ns=false) {
		if(is50ns) {
			float cut_Spring50nsEB_hOverE[3]    = {0.05   , 0.05   , 0.05   };                 		float cut_Spring50nsEE_hOverE[3]    = {0.05   , 0.05   , 0.05   };               
			float cut_Spring50nsEB_SIetaIeta[3] = {0.0103 , 0.0100 , 0.0100 };                 		float cut_Spring50nsEE_SIetaIeta[3] = {0.0277 , 0.0267 , 0.0267 };
			float cut_Spring50nsEB_PFchaIso[3]  = {2.44   , 1.31   , 0.91   };                 		float cut_Spring50nsEE_PFchaIso[3]  = {1.84   , 1.25   , 0.65   };
			float cut_Spring50nsEB_PFneuIso[3];                                                		float cut_Spring50nsEE_PFneuIso[3]; 
					cut_Spring50nsEB_PFneuIso[0] = 2.57 + exp( 0.0044 * pt + 0.5809) ;     						cut_Spring50nsEE_PFneuIso[0] = 4.00 + exp(0.0040 * pt + 0.9402)   ; 
					cut_Spring50nsEB_PFneuIso[1] = 0.60 + exp( 0.0044 * pt + 0.5809) ;     						cut_Spring50nsEE_PFneuIso[1] = 1.65 + exp(0.0040 * pt + 0.9402)   ; 
					cut_Spring50nsEB_PFneuIso[2] = 0.33 + exp( 0.0044 * pt + 0.5809) ;     						cut_Spring50nsEE_PFneuIso[2] = 0.93 + exp(0.0040 * pt + 0.9402)    ;
			float cut_Spring50nsEB_PFphoIso[3];                                                		float cut_Spring50nsEE_PFphoIso[3];
					cut_Spring50nsEB_PFphoIso[0] = 1.92 + 0.0043 * pt ;                    				      cut_Spring50nsEE_PFphoIso[0] = 2.15 + 0.0041 * pt  ;
					cut_Spring50nsEB_PFphoIso[1] = 1.33 + 0.0043 * pt ;                    						cut_Spring50nsEE_PFphoIso[1] = 1.02 + 0.0041 * pt  ;
					cut_Spring50nsEB_PFphoIso[2] = 0.61 + 0.0043 * pt ;                    						cut_Spring50nsEE_PFphoIso[2] = 0.54 + 0.0041 * pt  ;
			if(cutName.compare("hOverE"       ) == 0) return ( (isEB()) ?  cut_Spring50nsEB_hOverE[idIndex]    : cut_Spring50nsEE_hOverE[idIndex]      );
			if(cutName.compare("sigmaIetaIeta") == 0) return ( (isEB()) ?  cut_Spring50nsEB_SIetaIeta[idIndex] : cut_Spring50nsEE_SIetaIeta[idIndex]   );
			if(cutName.compare("PFChaIsoEA"   ) == 0) return ( (isEB()) ?  cut_Spring50nsEB_PFchaIso[idIndex]  : cut_Spring50nsEE_PFchaIso[idIndex]    );
			if(cutName.compare("PFNeuIsoEA"   ) == 0) return ( (isEB()) ?  cut_Spring50nsEB_PFneuIso[idIndex]  : cut_Spring50nsEE_PFneuIso[idIndex]    );
			if(cutName.compare("PFPhoIsoEA"   ) == 0) return ( (isEB()) ?  cut_Spring50nsEB_PFphoIso[idIndex]  : cut_Spring50nsEE_PFphoIso[idIndex]    );
		} else {
			float cut_Spring25nsEB_hOverE[3]    = {0.05   , 0.05   , 0.05   };               		float cut_Spring25nsEE_hOverE[3]    = {0.05   , 0.05   , 0.05   };                                
			float cut_Spring25nsEB_SIetaIeta[3] = {0.0102 , 0.0102 , 0.0100 };               		float cut_Spring25nsEE_SIetaIeta[3] = {0.0274 , 0.0268 , 0.0268 };
			float cut_Spring25nsEB_PFchaIso[3]  = {3.32   , 1.37   , 0.76   };               		float cut_Spring25nsEE_PFchaIso[3]  = {1.97   , 1.10   , 0.56   };
			float cut_Spring25nsEB_PFneuIso[3];                                              		float cut_Spring25nsEE_PFneuIso[3]; 
					cut_Spring25nsEB_PFneuIso[0] = 1.92 + 0.014 * pt + 0.000019 * pt * pt;     				cut_Spring25nsEE_PFneuIso[0] = 11.86 + 0.0139 * pt + 0.000025 * pt * pt; 
					cut_Spring25nsEB_PFneuIso[1] = 1.06 + 0.014 * pt + 0.000019 * pt * pt;     				cut_Spring25nsEE_PFneuIso[1] = 2.69  + 0.0139 * pt + 0.000025 * pt * pt; 
					cut_Spring25nsEB_PFneuIso[2] = 0.97 + 0.014 * pt + 0.000019 * pt * pt;     				cut_Spring25nsEE_PFneuIso[2] = 2.09  + 0.0139 * pt + 0.000025 * pt * pt;
			float cut_Spring25nsEB_PFphoIso[3];                                              		float cut_Spring25nsEE_PFphoIso[3];
					cut_Spring25nsEB_PFphoIso[0] = 0.81 + 0.0053 * pt ;                        				cut_Spring25nsEE_PFphoIso[0] = 0.83 + 0.0034 * pt ;
					cut_Spring25nsEB_PFphoIso[1] = 0.28 + 0.0053 * pt ;                        				cut_Spring25nsEE_PFphoIso[1] = 0.39 + 0.0034 * pt ;
					cut_Spring25nsEB_PFphoIso[2] = 0.08 + 0.0053 * pt ;                        				cut_Spring25nsEE_PFphoIso[2] = 0.16 + 0.0034 * pt ;
			if(cutName.compare("hOverE"       ) == 0) return ( (isEB()) ?  cut_Spring25nsEB_hOverE[idIndex]    : cut_Spring25nsEE_hOverE[idIndex]      );
			if(cutName.compare("sigmaIetaIeta") == 0) return ( (isEB()) ?  cut_Spring25nsEB_SIetaIeta[idIndex] : cut_Spring25nsEE_SIetaIeta[idIndex]   );
			if(cutName.compare("PFChaIsoEA"   ) == 0) return ( (isEB()) ?  cut_Spring25nsEB_PFchaIso[idIndex]  : cut_Spring25nsEE_PFchaIso[idIndex]    );
			if(cutName.compare("PFNeuIsoEA"   ) == 0) return ( (isEB()) ?  cut_Spring25nsEB_PFneuIso[idIndex]  : cut_Spring25nsEE_PFneuIso[idIndex]    );
			if(cutName.compare("PFPhoIsoEA"   ) == 0) return ( (isEB()) ?  cut_Spring25nsEB_PFphoIso[idIndex]  : cut_Spring25nsEE_PFphoIso[idIndex]    );
		}
		return -999;
	};
	bool CutBasedId(int idIndex = 2, bool is50ns = false) {
		float cut_hOverE     = CutThresholer("hOverE"       , idIndex, is50ns) ;  
		float cut_SIetaIeta  = CutThresholer("sigmaIetaIeta", idIndex, is50ns) ;  
		float cut_PFchaIso   = CutThresholer("PFChaIsoEA"   , idIndex, is50ns) ;  
		float cut_PFneuIso   = CutThresholer("PFNeuIsoEA"   , idIndex, is50ns) ;  
		float cut_PFphoIso   = CutThresholer("PFPhoIsoEA"   , idIndex, is50ns) ;  

		if(isEB() or isEE()                         ) 
		if(hadTowOverEm             < cut_hOverE    )
		if(full5x5_sigmaIetaIeta    < cut_SIetaIeta )
		if(isoChargedHadronsWithEA  < cut_PFchaIso  )
		if(isoNeutralHadronsWithEA  < cut_PFneuIso  )
		if(isoPhotonsWithEA         < cut_PFphoIso  ) return true;
		return false;
	};
	bool isTight()  { return isCutPassTight;                   };
	bool isMedium() { return (!isTight()  and isCutPassMedium); };
	bool isLoose()  { return (!isMedium() and isCutPassLoose); };

	bool isDeepTight(bool is50ns=false)  { return CutBasedId(2, is50ns);                   };
	bool isDeepMedium(bool is50ns=false) { return (!isDeepTight(is50ns)  and CutBasedId(1, is50ns)); };
	bool isDeepLoose(bool is50ns=false)  { return (!isDeepMedium(is50ns) and CutBasedId(0, is50ns)); };


	void PrintCutValue() {
		std::cout <<"@@ CutValue hadTowOverEm             " << hadTowOverEm              << std::endl;
		std::cout <<"@@ CutValue full5x5_sigmaIetaIeta    " << full5x5_sigmaIetaIeta     << std::endl;
		std::cout <<"@@ CutValue isoChargedHadronsWithEA  " << isoChargedHadronsWithEA   << std::endl;
		std::cout <<"@@ CutValue isoNeutralHadronsWithEA  " << isoNeutralHadronsWithEA   << std::endl;
		std::cout <<"@@ CutValue isoPhotonsWithEA         " << isoPhotonsWithEA          << std::endl;
	};


	Photon() {};
	 ~Photon() {};
	ClassDef(Photon,1);
};


class Muon : public P4 {
public:

   //TunePBestTrack(=tpbt) : reference to the Track chosen to assign the momentum value to the muon by PF 
   //int jetRef;  ///< referxence to the jet containing it (in the JET branch)
   bool isPFMuon  ;
   bool isGlobalMuon ;

   bool hastunePMuonBestTrack = false ;
   bool hasmuonBestTrack = false ;
   bool hasglobalTrack = false ;
   bool hasinnerTrack = false;

   int    type   ; 
   double pt_P4  ; 
   double eta_P4 ; 
   double phi_P4 ; 
   double pt_TPq ; 
   double ptErr  ; 
   double charge ; 
   double dxy    ; 
   double dz     ; 
   double dxy_GV ; 
   double dz_GV  ; 
   double edxy   ; 
   double edz    ; 

   double pt_BTq     ;
   double pt_BT      ;
   double ptErr_BT   ;
   double eta_BT     ;
   double phi_BT     ;
   double dxy_BT     ;
   double dz_BT      ;
   double dxy_GV_BT  ;
   double dz_GV_BT   ;
   double edxy_BT    ;
   double edz_BT     ;
   double momQuality ;

   double pt_GB              ; 
   double pt_GBq             ; 
   double dxy_GB             ; 
   double dz_GB              ; 
   double dxy_GV_GB          ; 
   double dz_GV_GB           ; 
   double vx_GB              ; 
   double vy_GB              ; 
   double vz_GB              ; 
   double normalizedChi2_GB  ; 
   int nValidTrkHits         ; 
   int nValidMuonHits        ; 

   double pt_IN                ;
   double pt_INq               ;
   double dxy_IN               ;
   double dz_IN                ;
   double dxy_GV_IN            ;
   double dz_GV_IN             ;
   double ptErr_IN             ;
   double normalizedChi2_IN    ;
   double dB_PV ; 
   double edB_PV ; 
   double dB_BS ; 
   double edB_BS ; 
   int    nValidPixelHits      ;
   int    nTrackerLayers       ;
   int    npixelLayers         ;

   int    nMatchedStations; 
   int    nMatches        ; 
   double TrkIsoR03_sumPt ; 
   double relTrkIsoR03_BT ; 
   double relTrkIsoR03_TP ; 
   double trackIso        ; 
   double caloIso         ; 
   double ecalIso         ; 
   double hcalIso         ; 
   double pfIsoR03_sumChargedHadronPt             ; 
   double pfIsoR03_sumNeutralHadronEt             ; 
   double pfIsoR03_sumPhotonEt                    ; 
   double pfIsoR03_sumPUPt                        ; 
   double pfIsoR03_sumChargedParticlePt           ; 
   double pfIsoR03_sumNeutralHadronEtHighThreshold; 
   double pfIsoR03_sumPhotonEtHighThreshold       ; 
   double pfIsoR04_sumChargedHadronPt             ; 
   double pfIsoR04_sumNeutralHadronEt             ; 
   double pfIsoR04_sumPhotonEt                    ; 
   double pfIsoR04_sumPUPt                        ; 
   double pfIsoR04_sumChargedParticlePt           ; 
   double pfIsoR04_sumNeutralHadronEtHighThreshold; 
   double pfIsoR04_sumPhotonEtHighThreshold       ; 
   double pfIsoVar ; 
   bool   isLooseMuon     ;
   bool   isTightMuonPV   ;
   bool   isSoftMuonPV    ;
   bool   isHighPtMuonPV  ;
   bool   isTightMuonGV   ;
   bool   isSoftMuonGV    ;
   bool   isHighPtMuonGV  ;

   float  Algo_IN             =  -999.0  ;
   float  originalAlgo_IN     =  -999.0  ;
   float  momQuality_IN       =  -999.0  ;
   float  TrkQuality_IN       =  -999.0  ;
   bool   isBadPFMuon  = false ;
   float dr_muPf  = -999.0 ;
   bool isBadChHadron  = false ;
   float dr_chPf  = -999.0 ;
   float dpt_chPf = -999.0 ;

   int statusPV           ; /// status = 100*z+10*y+x
   int statusGV           ; ///        x is general ID --> noTight: x>=1; Tight: x==2;
                             ///        HighPt: y==1;
	double Eta() { return eta; };
	double AEta() { return std::abs(eta); };
	bool   isMB() { return (AEta() <= 0.9) ; };
	bool   isOverlap() { return ((AEta() <= 1.2) and (AEta()>0.9)) ; };
	bool   isME() { return ((AEta() > 1.2) and (AEta() < 2.5)); };

	Muon() {};
	 ~Muon() {};
   void PrintMuon() {
        std::cout << "Muon pt " << pt << " eta " << eta << " phi " << phi << " monQuality " <<  ptErr_BT/pt_BT  << " relTrkIsoR03 " << relTrkIsoR03_BT << " nValidPixelHits "  << nValidPixelHits << " nTrkLayers " << nTrackerLayers << " nMatchedStations " << nMatchedStations << " nValidMuonHits " <<  nValidMuonHits << " dz " << dz_GV_BT << "dxy " << dxy_GV_BT <<  " @@ momQuIN " << momQuality_IN << " TrkQual " << TrkQuality_IN << " dr_muPf " << dr_muPf << " dr_chPf " << dr_chPf << " dpt_chPf " << dpt_chPf << std::endl;
   };

	ClassDef(Muon,1);
};

class Cut : public TObject {
public:
	std::string name;
	bool isPass;
	std::vector<std::string> nameVec;
	std::vector<float> valueVec;
	std::vector<bool> isVec;
	Cut() { name = "NoName" ; };
	Cut(std::string _name) { name = _name;};
	~Cut() {};
	void Fill(std::string cutname, float value, bool isOk) {
		nameVec.push_back(cutname);
		valueVec.push_back(value);
		isVec.push_back(isOk);
	};
	int Size() { return (int)nameVec.size(); };
	float GetValue(std::string str) {
		for(int i=0; i < (int)nameVec.size(); i++) {
			if(str.compare(nameVec.at(i)) == 0) return valueVec.at(i);
		}
		std::cout << "Error in Cut::GetValue " << str << std::endl;
		return 999.0;
	};
	void PrintCut() {
		std::cout << "CutResult " << name << " results " << isPass << std::endl;
		for(int i=0; i < (int)nameVec.size(); i++) {
			std::cout << "Cut " << name << " " << i << " " << nameVec.at(i) << " " << valueVec.at(i) << " " << isVec.at(i) << std::endl;
		}
	};
	ClassDef(Cut,1);
};


class Electron : public P4 {
public:

	double gsfTrack_Px;
	double gsfTrack_Py;
	double gsfTrack_Pz;
	double gsfTrack_Pt;
	double superCluster_x;
	double superCluster_y;
	double superCluster_z;
	double superCluster_eta;
	double superCluster_phi;
	double superCluster_energy   ; 
	double ecalEnergy            ;
	double caloEnergy            ; 

	double superCluster_ESenergy   ;//
	double superCluster_ESplane1   ;//
	double superCluster_ESplane2   ;//

	int    charge          ;
	double scE1x5          ;
	double scE2x5Max       ;
	double scE5x5          ;
	double scPixCharge     ;
	double scSigmaEtaEta   ;
	double scSigmaIEtaIEta ;

	bool   ecalDrivenSeed                   ; 
	double deltaEtaSuperClusterTrackAtVtx   ; 
	double deltaPhiSuperClusterTrackAtVtx   ; 
	double deltaEtaSeedClusterTrackAtCalo   ;//
	double deltaEtaSeedClusterTrackAtVtx    ;//
	double dEtaSeedAtVtx                    ;//
	double hcalOverEcal                     ;
	double hadronicOverEm                   ; 
	double e2x5Max                          ; 
	double e5x5                             ; 
	double e1x5                             ; 
	double dr03EcalRecHitSumEt              ; 
	double dr03HcalDepth1TowerSumEt         ;
	double dr03HcalTowerSumEt               ;  
	double dr03TkSumPt                      ; 
	double gsfTrack_dxy                     ; 
	double gsfTrack_dz                      ; 
	double gsfTrack_dxyPVtx                 ; 
	double gsfTrack_dzPVtx                  ; 
	int    nMissingHits                     ; 

	double caloIso               ;  
	double ecalIso               ;  
	double hcalIso               ;  
	double trackIso              ; 
	double full5x5_sigmaIetaIeta ; 
	double sigmaIetaIphi         ;  
	double sigmaIphiIphi         ;  
	bool   isPF                  ;  
	bool   passConversionVeto    ;  
	double r9                    ;  
	double fbrem                 ;

	double eSuperClusterOverP    ;
	double ooEmooP               ;

	double vtxFitConversion      ;
	double effArea               ; 

        double pfIso_sumChargedHadronPt ; 
        double pfIso_sumNeutralHadronEt ; 
        double pfIso_sumPhotonEt        ; 
        double pfIso_sumPUPt            ; 
        double absIsoWithDBeta          ; 
     	double effAreaPFIso             ;
     
        double lazyToolnoZS_eleFull5x5SigmaIEtaIEta ; 
        double lazyToolnoZS_eleFull5x5SigmaIEtaIPhi ; 
        double lazyToolnoZS_eleFull5x5R9            ; 
        double lazyToolnoZS_e1x5                    ; 
        double lazyToolnoZS_e2x5Max                 ; 
        double lazyToolnoZS_e5x5                    ; 
        double lazyToolnoZS_circularity             ; 
        double lazyToolnoZS_R_e1x5_e5x5             ; 
        double lazyToolnoZS_R_e2x5_e5x5             ; 

	double rho               ; 
	std::vector<bool> isID   ; 
	int matchToTruth         ;

        bool isCutPassTight  ;
        bool isCutPassMedium ;
        bool isCutPassLoose  ;
        bool isCutPassVeto   ;
        bool isCutPassHEEP   ;//
	bool passHEEPV70VID  ;
	int nrSatCrys  ;
        float trkIsol    ;
	float trkIsoVID  ;

	unsigned int heepV70Bitmap  ;
	
	bool passEtShowerShapeHE  ;
	bool passEtShowerShapeHEVID  ;
	bool passN1TrkIso  ;
	bool passN1TrkIsoVID  ;
	
	bool CheckVIDvsID  ;

	bool  mvaIsPassMedium; 
	bool  mvaIsPassTight ; 
	float mvaValue       ; 
	int   mvaCategory    ; 

	std::vector<std::vector<float> > idValuesVec;
   Electron() {} ;
    ~Electron() {};
	double scEta() { return superCluster_eta; };
	double scAEta() { return std::abs(superCluster_eta); };
	bool   isGap() { return ((scAEta() > 1.4442) and (scAEta() < 1.566)) ; };
	bool   isEB() { return (scAEta() <= 1.479) ; };
	bool   isEE() { return ((scAEta() > 1.479) and (scAEta() < 2.5)); };

	void PrintObj() {
		std::cout << "Electron pt " << pt << " eta " << eta << " phi " << phi << " energy " << energy
                << " charge " << charge << " scEta " << scEta()
                << " passHEEP " << isCutPassHEEP 
                << " H/E " << hcalOverEcal 
                << " sigIetaIeta " << scSigmaIEtaIEta 
                << " deltaEtain " << dEtaSeedAtVtx 
                << " deltaPhi " << deltaPhiSuperClusterTrackAtVtx   
                << " nrSatCrys " << nrSatCrys << " trkIsol " << trkIsol << " passEtShowerShapeHE " <<  passEtShowerShapeHE 
                << " passConversionVeto " << passConversionVeto 
                << " matchMCTruth " << matchToTruth 
                << " nMissingHit " << nMissingHits << " ";
      std::cout <<  std::endl;            
	}

	std::vector<float> GetIdValues(int idIndex) { return idValuesVec.at(idIndex);	};
	void PrintIdValues(int idIndex) {
		std::vector<float> thisVec = GetIdValues(idIndex);
		int thisI = 0;
		for(std::vector<float>::iterator it=thisVec.begin(); it!=thisVec.end(); it++, thisI++) {
			std::cout << "EleIdValue " << idIndex << " " << thisI << " " << *it << std::endl;
		}
	};	

	double EffArea(bool is50ns = false) { return ( (is50ns) ? EffAreaElectron_Spring15_50ns(std::abs(superCluster_eta)) : EffAreaElectron_Spring15_25ns(std::abs(superCluster_eta)) ) ; };
	double PFIsoWithEA(bool is50ns = false) { return pfIso_sumChargedHadronPt + std::max(0.0 , pfIso_sumNeutralHadronEt + pfIso_sumPhotonEt - rho * EffArea(is50ns) ); };

	float CutThresholer(std::string cutName, int idIndex =3, bool is50ns=false) {
		if(is50ns) {
			float cut_Spring50nsEB_full5x5_sigmaIetaIeta[4]    = { 0.012  , 0.0105  , 0.0101 , 0.0101  }; float cut_Spring50nsEE_full5x5_sigmaIetaIeta[4]    = { 0.0339,  0.0318  , 0.0287  , 0.0287   };
			float cut_Spring50nsEB_abs_dEtaIn[4]               = { 0.0126 , 0.00976 , 0.0094 , 0.00864 }; float cut_Spring50nsEE_abs_dEtaIn[4]               = { 0.0109,  0.00952 , 0.00773 , 0.00762  };
			float cut_Spring50nsEB_abs_dPhiIn[4]               = { 0.107  , 0.0929  , 0.0296 , 0.0286  }; float cut_Spring50nsEE_abs_dPhiIn[4]               = { 0.219 ,  0.181   , 0.148   , 0.0439   };
			float cut_Spring50nsEB_hOverE[4]                   = { 0.186  , 0.0765  , 0.0372 , 0.0342  }; float cut_Spring50nsEE_hOverE[4]                   = { 0.0962,  0.0824  , 0.0546  , 0.0544   };
			float cut_Spring50nsEB_relIsoWithEA[4]             = { 0.161  , 0.118   , 0.0987 , 0.0591  }; float cut_Spring50nsEE_relIsoWithEA[4]             = { 0.193 ,  0.118   , 0.0902  , 0.0759   };
			float cut_Spring50nsEB_ooEmooP[4]                  = { 0.239  , 0.184   , 0.118  , 0.0116  }; float cut_Spring50nsEE_ooEmooP[4]                  = { 0.141 ,  0.125   , 0.104   , 0.01     };
			float cut_Spring50nsEB_abs_d0[4]                   = { 0.0621 , 0.0227  , 0.0151 , 0.0103  }; float cut_Spring50nsEE_abs_d0[4]                   = { 0.279 ,  0.242   , 0.0535  , 0.0377   };
			float cut_Spring50nsEB_abs_dz[4]                   = { 0.613  , 0.379   , 0.238  , 0.170   }; float cut_Spring50nsEE_abs_dz[4]                   = { 0.947 ,  0.921   , 0.572   , 0.571    };
			float cut_Spring50nsEB_expectedMissingInnerHits[4] = { 2      , 2       , 2      , 2       }; float cut_Spring50nsEE_expectedMissingInnerHits[4] = { 3 	 ,  1 	   , 1 	    , 1         };  
			if(cutName.compare("fullSigma") == 0) return ( (isEB()) ? cut_Spring50nsEB_full5x5_sigmaIetaIeta[idIndex]    : cut_Spring50nsEE_full5x5_sigmaIetaIeta[idIndex]   );
			if(cutName.compare("dEtaIn"   ) == 0) return ( (isEB()) ? cut_Spring50nsEB_abs_dEtaIn[idIndex]               : cut_Spring50nsEE_abs_dEtaIn[idIndex]              );
			if(cutName.compare("dPhiIn"   ) == 0) return ( (isEB()) ? cut_Spring50nsEB_abs_dPhiIn[idIndex]               : cut_Spring50nsEE_abs_dPhiIn[idIndex]              );
			if(cutName.compare("hOverE"   ) == 0) return ( (isEB()) ? cut_Spring50nsEB_hOverE[idIndex]                   : cut_Spring50nsEE_hOverE[idIndex]                  );
			if(cutName.compare("relIsoEA" ) == 0) return ( (isEB()) ? cut_Spring50nsEB_relIsoWithEA[idIndex]             : cut_Spring50nsEE_relIsoWithEA[idIndex]            );
			if(cutName.compare("ooEmooP"  ) == 0) return ( (isEB()) ? cut_Spring50nsEB_ooEmooP[idIndex]                  : cut_Spring50nsEE_ooEmooP[idIndex]                 );
			if(cutName.compare("d0"       ) == 0) return ( (isEB()) ? cut_Spring50nsEB_abs_d0[idIndex]                   : cut_Spring50nsEE_abs_d0[idIndex]                  );
			if(cutName.compare("dz"       ) == 0) return ( (isEB()) ? cut_Spring50nsEB_abs_dz[idIndex]                   : cut_Spring50nsEE_abs_dz[idIndex]                  );
			if(cutName.compare("missHit"  ) == 0) return ( (isEB()) ? cut_Spring50nsEB_expectedMissingInnerHits[idIndex] : cut_Spring50nsEE_expectedMissingInnerHits[idIndex]);
		} else {                                          
			float cut_Spring25nsEB_full5x5_sigmaIetaIeta[4]    = { 0.0114 , 0.0103  , 0.0101 , 0.0101  }; float cut_Spring25nsEE_full5x5_sigmaIetaIeta[4]    = { 0.0352   , 0.0301   , 0.0283   , 0.0279   };
			float cut_Spring25nsEB_abs_dEtaIn[4]               = { 0.0152 , 0.0105  , 0.0103 , 0.00926 }; float cut_Spring25nsEE_abs_dEtaIn[4]               = { 0.0113   , 0.00814  , 0.00733  , 0.00724  };
			float cut_Spring25nsEB_abs_dPhiIn[4]               = { 0.216  , 0.115 	, 0.0336 , 0.0336  }; float cut_Spring25nsEE_abs_dPhiIn[4]               = { 0.237    , 0.182    , 0.114    , 0.0918   };
			float cut_Spring25nsEB_hOverE[4]                   = { 0.181  , 0.104 	, 0.0876 , 0.0597  }; float cut_Spring25nsEE_hOverE[4]                   = { 0.116    , 0.0897   , 0.0678   , 0.0615   };
			float cut_Spring25nsEB_relIsoWithEA[4]             = { 0.126  , 0.0893  , 0.0766 , 0.0354  }; float cut_Spring25nsEE_relIsoWithEA[4]             = { 0.144    , 0.121    , 0.0678   , 0.0646   };
			float cut_Spring25nsEB_ooEmooP[4]                  = { 0.207  , 0.102 	, 0.0174 , 0.012   }; float cut_Spring25nsEE_ooEmooP[4]                  = { 0.174    , 0.126    , 0.0898   , 0.00999  };
			float cut_Spring25nsEB_abs_d0[4]                   = { 0.0564 , 0.0261  , 0.0118 , 0.0111  }; float cut_Spring25nsEE_abs_d0[4]                   = { 0.222    , 0.118    , 0.0739   , 0.0351   };
			float cut_Spring25nsEB_abs_dz[4]                   = { 0.472  , 0.41 	, 0.373  , 0.0466  }; float cut_Spring25nsEE_abs_dz[4]                   = { 0.921    , 0.822    , 0.602    , 0.417    };
			float cut_Spring25nsEB_expectedMissingInnerHits[4] = { 2 	  , 2 	   , 2 	   , 2       }; float cut_Spring25nsEE_expectedMissingInnerHits[4] = { 3        , 1        , 1        , 1        };  
			if(cutName.compare("fullSigma") == 0) return ( (isEB()) ? cut_Spring25nsEB_full5x5_sigmaIetaIeta[idIndex]    : cut_Spring25nsEE_full5x5_sigmaIetaIeta[idIndex]   );
			if(cutName.compare("dEtaIn"   ) == 0) return ( (isEB()) ? cut_Spring25nsEB_abs_dEtaIn[idIndex]               : cut_Spring25nsEE_abs_dEtaIn[idIndex]              );
			if(cutName.compare("dPhiIn"   ) == 0) return ( (isEB()) ? cut_Spring25nsEB_abs_dPhiIn[idIndex]               : cut_Spring25nsEE_abs_dPhiIn[idIndex]              );
			if(cutName.compare("hOverE"   ) == 0) return ( (isEB()) ? cut_Spring25nsEB_hOverE[idIndex]                   : cut_Spring25nsEE_hOverE[idIndex]                  );
			if(cutName.compare("relIsoEA" ) == 0) return ( (isEB()) ? cut_Spring25nsEB_relIsoWithEA[idIndex]             : cut_Spring25nsEE_relIsoWithEA[idIndex]            );
			if(cutName.compare("ooEmooP"  ) == 0) return ( (isEB()) ? cut_Spring25nsEB_ooEmooP[idIndex]                  : cut_Spring25nsEE_ooEmooP[idIndex]                 );
			if(cutName.compare("d0"       ) == 0) return ( (isEB()) ? cut_Spring25nsEB_abs_d0[idIndex]                   : cut_Spring25nsEE_abs_d0[idIndex]                  );
			if(cutName.compare("dz"       ) == 0) return ( (isEB()) ? cut_Spring25nsEB_abs_dz[idIndex]                   : cut_Spring25nsEE_abs_dz[idIndex]                  );
			if(cutName.compare("missHit"  ) == 0) return ( (isEB()) ? cut_Spring25nsEB_expectedMissingInnerHits[idIndex] : cut_Spring25nsEE_expectedMissingInnerHits[idIndex]);
		}
		return -999;
	};

	bool CutBasedId(int idIndex = 3, bool is50ns = false) {
		//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
		float cut_full5x5_sigmaIetaIeta    = CutThresholer("fullSigma", idIndex, is50ns) ;  
		float cut_abs_dEtaIn               = CutThresholer("dEtaIn"   , idIndex, is50ns) ;  
		float cut_abs_dPhiIn               = CutThresholer("dPhiIn"   , idIndex, is50ns) ;  
		float cut_hOverE                   = CutThresholer("hOverE"   , idIndex, is50ns) ;  
		float cut_relIsoWithEA             = CutThresholer("relIsoEA" , idIndex, is50ns) ;  
		float cut_ooEmooP                  = CutThresholer("ooEmooP"  , idIndex, is50ns) ;  
		float cut_abs_d0                   = CutThresholer("d0"       , idIndex, is50ns) ;  
		float cut_abs_dz                   = CutThresholer("dz"       , idIndex, is50ns) ;  
		float cut_expectedMissingInnerHits = CutThresholer("missHit"  , idIndex, is50ns) ;  

		if(isEB() or isEE()                                                        )
		if(full5x5_sigmaIetaIeta                    < cut_full5x5_sigmaIetaIeta    )
		if(std::abs(deltaEtaSuperClusterTrackAtVtx) < cut_abs_dEtaIn               )
		if(std::abs(deltaPhiSuperClusterTrackAtVtx) < cut_abs_dPhiIn               )
		if(hadronicOverEm                           < cut_hOverE                   )
		if(relIsoWithEA(is50ns)                     < cut_relIsoWithEA             )
		if(std::abs(ooEmooP)                        < cut_ooEmooP                  )
		if(std::abs(gsfTrack_dxyPVtx)               < cut_abs_d0                   )
		if(std::abs(gsfTrack_dzPVtx)                < cut_abs_dz                   )
		if(nMissingHits                            <= cut_expectedMissingInnerHits ) 
		if(passConversionVeto                         ) return true;

		return false;
	};
	bool isTight()  { return isCutPassTight;                    };
	bool isMedium() { return (!isTight()  and isCutPassMedium); };
	bool isLoose()  { return (!isMedium() and isCutPassLoose);  };
	bool isVeto()   { return (!isLoose()  and isCutPassVeto);   };
	bool isHEEP()   { return isCutPassHEEP ;                   };

	bool isDeepTight(bool is50ns=false)  { return CutBasedId(3, is50ns);                   };
	bool isDeepMedium(bool is50ns=false) { return (!isDeepTight(is50ns)  and CutBasedId(2, is50ns)); };
	bool isDeepLoose(bool is50ns=false)  { return (!isDeepMedium(is50ns) and CutBasedId(1, is50ns)); };
	bool isDeepVeto(bool is50ns=false)   { return (!isDeepLoose(is50ns)  and CutBasedId(0, is50ns)); };

	double relIsoWithEA(bool is50ns=false) { return (PFIsoWithEA(is50ns) / pt); };

	void PrintCutValue(bool is50ns = false) {
		std::cout <<"@@ CutValue full5x5_sigmaIetaIeta                    " << full5x5_sigmaIetaIeta                     << std::endl;
		std::cout <<"@@ CutValue std::abs(deltaEtaSuperClusterTrackAtVtx) " << std::abs(deltaEtaSuperClusterTrackAtVtx)  << std::endl;
		std::cout <<"@@ CutValue std::abs(deltaPhiSuperClusterTrackAtVtx) " << std::abs(deltaPhiSuperClusterTrackAtVtx)  << std::endl;
		std::cout <<"@@ CutValue hadronicOverEm                           " << hadronicOverEm                            << std::endl;
		std::cout <<"@@ CutValue (PFIsoWithEA(is50ns) / pt)               " << relIsoWithEA(is50ns)                      << std::endl;
		std::cout <<"@@ CutValue std::abs(ooEmooP)                        " << std::abs(ooEmooP)                         << std::endl;
		std::cout <<"@@ CutValue std::abs(gsfTrack_dxyPVtx)               " << std::abs(gsfTrack_dxyPVtx)                << std::endl;
		std::cout <<"@@ CutValue std::abs(gsfTrack_dzPVtx)                " << std::abs(gsfTrack_dzPVtx)                 << std::endl;
		std::cout <<"@@ CutValue nMissingHits                             " << nMissingHits                              << std::endl;
		std::cout <<"@@ CutValue passConversionVeto                       " << passConversionVeto                        << std::endl;
	};

   ClassDef(Electron,1);
};

class Jet : public P4 {
public:
   float neutralHadronEnergyFraction ; 
   float neutralEmEnergyFraction     ; 
   float chargedMultiplicity         ; 
   float neutralMultiplicity         ; 
   float muonEnergyFraction          ; 
   float chargedHadronEnergyFraction ; 
   float chargedEmEnergyFraction     ; 

	bool passPFLooseId  ;
	bool passPFMediumId ;
	bool passPFTightId  ;
	
	Jet() {};
	 ~Jet() {};
	ClassDef(Jet,1)	
};

class MET : public P4 {
public:
   bool  isCaloMET ;
   bool  isPFMET   ;
   bool  isRecoMET ;
   double sumEt                    ; 
   double metSignificance          ;

   double MuClean_pt        = -999.0        ; 
   double MuClean_phi       = -999.0        ; 
   double MuClean_sumEt     = -999.0        ; 
   double MuClean_metSig    = -999.0        ;

   double EGClean_pt        = -999.0        ; 
   double EGClean_phi       = -999.0        ; 
   double EGClean_sumEt     = -999.0        ; 
   double EGClean_metSig    = -999.0        ;

   double MuEGClean_pt      = -999.0        ; 
   double MuEGClean_phi     = -999.0        ; 
   double MuEGClean_sumEt   = -999.0        ; 
   double MuEGClean_metSig  = -999.0        ;

   double Uncorr_pt         = -999.0        ; 
   double Uncorr_phi        = -999.0        ; 
   double Uncorr_sumEt      = -999.0        ; 
   double Uncorr_metSig     = -999.0        ;


   double genMET_pt                ; 
   double genMET_phi               ; 
   double genMET_sumEt             ; 
   double METPuppi_pt              ;
   double METPuppi_phi             ;

   double corPt_Type01           ; 
   double corPhi_Type01          ; 
   double corSumEt_Type01        ; 
   double corPt_TypeXY           ; 
   double corPhi_TypeXY          ; 
   double corSumEt_TypeXY        ; 
   double corPt_Type1XY          ; 
   double corPhi_Type1XY         ; 
   double corSumEt_Type1XY       ; 
   double corPt_Type01XY         ; 
   double corPhi_Type01XY        ; 
   double corSumEt_Type01XY      ; 
   double corPt_Type1Smear       ; 
   double corPhi_Type1Smear      ; 
   double corSumEt_Type1Smear    ; 
   double corPt_Type01Smear      ; 
   double corPhi_Type01Smear     ; 
   double corSumEt_Type01Smear   ; 
   double corPt_Type1SmearXY     ; 
   double corPhi_Type1SmearXY    ; 
   double corSumEt_Type1SmearXY  ; 
   double corPt_Type01SmearXY    ; 
   double corPhi_Type01SmearXY   ; 
   double corSumEt_Type01SmearXY ;

   double shiftedPt_JetResUp           ;
   double shiftedPt_JetResDown         ;
   double shiftedPt_JetEnUp            ;
   double shiftedPt_JetEnDown          ;
   //double shiftedPt_MuonEnUp           ;
   //double shiftedPt_MuonEnDown         ;
   //double shiftedPt_ElectronEnUp       ;
   //double shiftedPt_ElectronEnDown     ;
   //double shiftedPt_TauEnUp            ;
   //double shiftedPt_TauEnDown          ;
   double shiftedPt_UnclusteredEnUp    ;
   double shiftedPt_UnclusteredEnDown  ;
   //double shiftedPt_PhotonEnUp         ;
   //double shiftedPt_PhotonEnDown       ;
   double shiftedPt_JetResUpSmear      ;
   double shiftedPt_JetResDownSmear    ;

   double shiftedPt_Type1XY    ;
   double shiftedPhi_Type1XY    ;
   double shiftedPt_Type1XY_puppi    ;
   double shiftedPhi_Type1XY_puppi    ;


   double METPuppi_corPt_Type01            ; 
   double METPuppi_corPhi_Type01           ; 
   double METPuppi_corSumEt_Type01         ; 
   double METPuppi_corPt_TypeXY            ; 
   double METPuppi_corPhi_TypeXY           ; 
   double METPuppi_corSumEt_TypeXY         ; 
   double METPuppi_corPt_Type1XY           ; 
   double METPuppi_corPhi_Type1XY          ; 
   double METPuppi_corSumEt_Type1XY        ; 
   double METPuppi_corPt_Type01XY          ; 
   double METPuppi_corPhi_Type01XY         ; 
   double METPuppi_corSumEt_Type01XY       ; 
   double METPuppi_corPt_Type1Smear        ; 
   double METPuppi_corPhi_Type1Smear       ; 
   double METPuppi_corSumEt_Type1Smear     ; 
   double METPuppi_corPt_Type01Smear       ; 
   double METPuppi_corPhi_Type01Smear      ; 
   double METPuppi_corSumEt_Type01Smear    ; 
   double METPuppi_corPt_Type1SmearXY      ; 
   double METPuppi_corPhi_Type1SmearXY     ; 
   double METPuppi_corSumEt_Type1SmearXY   ; 
   double METPuppi_corPt_Type01SmearXY     ; 
   double METPuppi_corPhi_Type01SmearXY    ; 
   double METPuppi_corSumEt_Type01SmearXY  ;

   //double METPuppi_shiftedPt_JetResUp           ;
   //double METPuppi_shiftedPt_JetResDown         ;
   //double METPuppi_shiftedPt_JetEnUp            ;
   //double METPuppi_shiftedPt_JetEnDown          ;
   //double METPuppi_shiftedPt_MuonEnUp           ;
   //double METPuppi_shiftedPt_MuonEnDown         ;
   //double METPuppi_shiftedPt_ElectronEnUp       ;
   //double METPuppi_shiftedPt_ElectronEnDown     ;
   //double METPuppi_shiftedPt_TauEnUp            ;
   //double METPuppi_shiftedPt_TauEnDown          ;
   //double METPuppi_shiftedPt_UnclusteredEnUp    ;
   //double METPuppi_shiftedPt_UnclusteredEnDown  ;
   //double METPuppi_shiftedPt_PhotonEnUp         ;
   //double METPuppi_shiftedPt_PhotonEnDown       ;
   //double METPuppi_shiftedPt_JetResUpSmear      ;
   //double METPuppi_shiftedPt_JetResDownSmear    ;

   bool Flag_HBHENoiseFilter = false   ;
   bool Flag_HBHENoiseIsoFilter = false   ;
   bool Flag_CSCTightHaloFilter = false  ;
   bool Flag_hcalLaserEventFilter = false   ;
   bool Flag_EcalDeadCellTriggerPrimitiveFilter = false   ;
   bool Flag_EcalDeadCellBoundaryEnergyFilter = false   ;
   bool Flag_goodVertices = false   ;
   bool Flag_trackingFailureFilter = false   ;
   bool Flag_eeBadScFilter = false   ;
   bool Flag_ecalLaserCorrFilter = false   ;
   bool Flag_trkPOGFilters = false   ;  
   bool Flag_trkPOG_manystripclus53X = false   ;
   bool Flag_trkPOG_toomanystripclus53X = false  ; 
   bool Flag_trkPOG_logErrorTooManyClusters = false  ;
   bool Flag_HcalStripHaloFilter                = false ;
   bool Flag_CSCTightHaloTrkMuUnvetoFilter      = false ;
   bool Flag_CSCTightHalo2015Filter             = false ;
   bool Flag_globalSuperTightHalo2016Filter     = false ;
   bool Flag_globalTightHalo2016Filter          = false ;
   bool Flag_chargedHadronTrackResolutionFilter = false ;
   bool Flag_muonBadTrackFilter                 = false ;
   bool Flag_METFilters = false  ;
   bool IsNoBadPFMuon            = false; // true is goodEvent 
   bool IsNoBadChargedCandidate  = false; // true is goodEvent
   bool Flag_BadPFMuonFilter     = false; // true is goodEvent
   bool Flag_BadChargedCandidateFilter = false; // true is goodEvent

   bool Flag_badGlobalMuonFilter = false; // true is CorrectedMuonEvent 
   bool IsCorrectedBadGlobalMuon = false; // true is CorrectedMuonEvent 
   bool Flag_duplicateMuonFilter = false; // true is CorrectedMuonEvent 
   bool IsCorrectedDuplicateMuon = false; // true is CorrectedMuonEvent 

   //ECAL slew (Re-miniAOD)
   bool IsDuplicateECALClusters = false; // true is BadEvent
   bool IsNotReplatedHitsECAL   = false; // true is BadEvent
   bool IsBadECALSlew           = false; // true is BadEvent

   bool filterbadPFMuon = false;
   bool filterbadChCandidate  =false;

   MET() {};
    ~MET() {};

   void PrintCompareMET() {
       std::cout << "Comparing MET:MuEGClean:MuClean:EGClean:Uncorr= pt " << pt << ":" << MuEGClean_pt << ":" << MuClean_pt << ":" << EGClean_pt << ":" << Uncorr_pt
       << " phi " << phi << ":" << MuEGClean_phi << ":" << MuClean_phi << ":" << EGClean_phi << ":" << Uncorr_phi <<  std::endl;
   };
   void PrintMET() {
       std::cout << "MET pt " << pt << " phi " << phi << " sumEt " << sumEt << " shiftedPt_JetEnup " << shiftedPt_JetEnUp << " down " << shiftedPt_JetEnDown 
       <<" Significance " << metSignificance<< " isNotBadPFMu " << IsNoBadPFMuon << " isNotBadChCand " << IsNoBadChargedCandidate << " IsCorrectedBadGlobalMuon " << IsCorrectedBadGlobalMuon << " IsCorrectedDuplicateMuon " << IsCorrectedDuplicateMuon << " ECALSlew :IsDuplicateECAL:IsNotRepHits= IsBadECALSlew " << IsDuplicateECALClusters << IsNotReplatedHitsECAL << IsBadECALSlew << std::endl;
   };
   ClassDef(MET,1)
};


class Trigger : public TObject {
public:
   std::string name;
   bool accept;
   double prescaleForIndex;
   Trigger() {};
    ~Trigger() {};
	bool Pass() { return accept; };
	bool Pass(std::vector<std::string> strVec) {
		for(std::vector<std::string>::iterator it=strVec.begin(); it!=strVec.end(); it++) { if(Pass(*it)) return true; }
		return false;
	};
	bool Pass(std::string str) {
		//std::cout << "Debug : comapre " << name.c_str() << " | " << str.c_str() << std::endl;
		return (CompareWild(name.c_str(), str.c_str()) == 1 and accept) ;
	};
	void PrintResult(const char* space="") {
		std::cout << space ;
		std::cout << "Trigger " << name << " accept " << accept << " prescaleForIndex " << prescaleForIndex << std::endl;
	};
   ClassDef(Trigger,1)
};

class TriggerObject : public TObject {
public:
	double pt;
	double eta;
	double phi;
	std::vector<std::string> pathNameVec;
	std::vector<bool> isBoth;
	std::vector<bool> isL3;
	std::vector<bool> isLF;
	std::vector<bool> isNone;
	bool hasPATH(std::string name){
		for (std::vector<std::string>::iterator it=pathNameVec.begin(); it!=pathNameVec.end(); it++) {
			if(CompareWild(it->c_str(), name.c_str()) == 1)  return true;
			//if(it->compare(name) == 0) return true;
		}
   	return false;
	};
	void PrintObj(const char* space = "") {
		std::cout << space ;
		std::cout << "TriggerObject " << pt << " " << eta << " " << phi << " " ;
		std::vector<bool>::iterator itBoth = isBoth.begin();
		for(std::vector<std::string>::iterator it=pathNameVec.begin(); it!=pathNameVec.end(); it++, itBoth++) {
			std::cout << " [ " << *it << " " << *itBoth << " ] " ;
		}
		std::cout << std::endl;
	};

   TriggerObject() {};
    ~TriggerObject() {};
   ClassDef(TriggerObject,1)
};


class GenInfo : public TObject {
public:
   double scalePDF;
   int pdg1;
   int pdg2;
   double x1;
   double x2;
   double xpdf1;
   double xpdf2;
	GenInfo() {};
	 ~GenInfo() {};
	ClassDef(GenInfo,1);
};

class PDFxfx : public GenInfo {
public:
	std::vector<std::pair<int, std::pair<double, double> > > xfxVec; //from LHAPDF
	PDFxfx(const npknu::GenInfo& a) : GenInfo(a) {} ;
	PDFxfx() {} ;
	 ~PDFxfx() {};
	int GetIndex(int i) { return xfxVec.at(i).first ;};
	double Getxfx1(int i) { return xfxVec.at(i).second.first ;};
	double Getxfx2(int i) { return xfxVec.at(i).second.second ;};
	double Getxpdf1(int i) { return Getxfx1(i) / x1 ;};
	double Getxpdf2(int i) { return Getxfx2(i) / x2 ;};
	double GetVal(int i) { return  Getxpdf1(i) * Getxpdf2(i) ; };
	double Weight(int i) { return GetVal(i) / GetVal(0); };
	ClassDef(PDFxfx,1);
};

class Pileup : public TObject {
public:
	int    BunchCrossing;
	double TrueNumInteractions;
	int    PU_NumInteractions;
   Pileup() {};
    ~Pileup() {}; 
   ClassDef(Pileup,1);
};

class Vertex : public TObject {
public:
	double x              ;  
	double y              ;  
	double z              ;  
	double xError         ;  
	double yError         ;  
	double zError         ;  
	int    tracksSize     ;  
	int    nTracks        ;  
	bool   isFake         ;  
	double ndof           ;  
	double position_rho   ;  
	double chi2           ;  
	double normalizedChi2 ;  

   Vertex() {};
    ~Vertex() {}; 
	bool isGoodVertex() { //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_7.4.12/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc
		return ((!isFake) and (ndof >= 4.0) and (position_rho <= 2.0) and (std::abs(z) <= 24.0));
   };
   ClassDef(Vertex,1);
};

class PUReweight {
private:
   long int events;
	double sumWeight;
   TH1D* PUweightHist;
public:
	TH1D* h1_data;
	TH1F* h1_mc;
	TH1D* h1_ratio;
   PUReweight(TString fileName, TString fileNameMC, TString mcHistName = "pileup") {
		events = 0;
		sumWeight = 0;
		TFile* pileupFile = new TFile(fileName,"READ");
		PUweightHist = (TH1D*)pileupFile->Get("pileup");
		PUweightHist->SetDirectory(0);
		pileupFile->Close();
		double PUweightInt = PUweightHist->Integral();
	
		TH1F* mcPU=NULL;
		TFile* mcFile = new TFile(fileNameMC,"READ");
		if( mcPU==NULL) mcPU = (TH1F*)mcFile->Get(mcHistName);
//		else mcPU->Add((TH1F*)mcFile->Get(mcHistName));
		mcPU->SetDirectory(0);
		mcFile->Close();
		
		PUweightHist->Scale(1.0/PUweightInt);
		mcPU->Scale(1.0/mcPU->Integral());
		PUweightHist->Divide(mcPU);

		delete mcPU;
	};
   ~PUReweight() {
	};
   double getWeight(double nPU) {
		events++;
		double PUweight = 0.0;
		PUweight = PUweightHist->GetBinContent(PUweightHist->GetXaxis()->FindBin(nPU));
		sumWeight+= PUweight;
		return PUweight;
	};
	double getAvgWeight() { return (sumWeight / events); };
};

class P4Hist {
public:
	TDirectory* dir;
	TH1F* h1_pt    ;
	TH1F* h1_eta   ;
	TH1F* h1_phi   ;
	TH1F* h1_energy;

	P4Hist(TDirectory* histFile, TString dirName) {
		dir = histFile->mkdir(dirName,dirName);
		dir->cd();
		h1_pt     = new TH1F("h1_pt","h1_pt",1000,0,2000)        ;
		h1_eta    = new TH1F("h1_eta","h1_eta",100,-2.5,2.5)     ;
		h1_phi    = new TH1F("h1_phi","h1_phi",64,-3.2,3.2)      ;
		h1_energy = new TH1F("h1_energy","h1_energy",1000,0,2000);
	};
	~P4Hist() {};
	void Fill(npknu::P4 ptr, double weight = 1.0) {
		h1_pt     ->Fill(ptr.pt     , weight) ; 
		h1_eta    ->Fill(ptr.eta    , weight) ; 
		h1_phi    ->Fill(ptr.phi    , weight) ; 
		h1_energy ->Fill(ptr.energy , weight) ; 
	};
	void Fill(npknu::P4* ptr, double weight = 1.0) { Fill(*ptr, weight); };
	void Write() {
		dir->cd();
		h1_pt     ->Write() ; 
		h1_eta    ->Write() ; 
		h1_phi    ->Write() ; 
		h1_energy ->Write() ; 
	};
	void Sumw2(bool isSumw2 = true) {
		h1_pt     ->Sumw2(isSumw2); 
		h1_eta    ->Sumw2(isSumw2); 
		h1_phi    ->Sumw2(isSumw2); 
		h1_energy ->Sumw2(isSumw2); 
	};
	void Print(TString outName) {
		Sumw2(false);
		TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
		c1->Divide(2,2);
		c1->cd(1); h1_pt    ->Draw();
		c1->cd(2); h1_eta   ->Draw();
		c1->cd(3); h1_phi   ->Draw();
		c1->cd(4); h1_energy->Draw();
		c1->Print(outName);
		delete c1;
		Sumw2(true);
	};
	TDirectory* GetDir() { return dir;};
};

class P4DHist {
public:
	TDirectory* dir   ;
	TH1F* h1_mass     ;
	TH1F* h1_mt       ;
	TH1F* h1_dPhi     ;
	TH1F* h1_etRatio  ;
	P4Hist* p4Hist1   ;
	P4Hist* p4Hist2   ;
	P4DHist(TDirectory* histFile, TString dirName) {
		dir = histFile->mkdir(dirName,dirName);
		dir->cd();
		h1_mass = new TH1F("h1_mass","h1_mass",1000,0,5000);
		h1_mt = new TH1F("h1_mt","h1_mt",1000,0,5000);
		h1_dPhi = new TH1F("h1_dPhi","h1_dPhi",100,0,3.2);
		h1_etRatio = new TH1F("h1_etRatio","h1_etRatio",100,0,5);
		p4Hist1 = new P4Hist(dir, "p4Hist1");
		p4Hist2 = new P4Hist(dir, "p4Hist2");
		Sumw2(true);
	};
	~P4DHist() {};
	void Fill(npknu::P4 ptr1, npknu::P4 ptr2, double weight = 1.0) {
		h1_mass->Fill(ptr1.GetM(ptr2), weight);
		h1_mt->Fill(ptr1.GetMt(ptr2), weight);
		h1_dPhi->Fill(ptr1.deltaPhi(ptr2), weight);
		h1_etRatio->Fill(ptr1.etRatio(ptr2), weight);
		p4Hist1->Fill(ptr1, weight);
		p4Hist2->Fill(ptr2, weight);
	};
	void Fill(npknu::P4* ptr1, npknu::P4* ptr2, double weight = 1.0) { Fill(*ptr1, *ptr2, weight); };
	void Write() {
		dir->cd();
		h1_mass->Write();
		h1_mt     ->Write();
		h1_dPhi   ->Write();
		h1_etRatio->Write();
		p4Hist1->Write();
		p4Hist2->Write();
	};
	void Sumw2(bool isSumw2 = true) {
		p4Hist1->Sumw2(isSumw2);
		p4Hist2->Sumw2(isSumw2);
		h1_mass->Sumw2(isSumw2);
		h1_mt->Sumw2(isSumw2); 
		h1_dPhi->Sumw2(isSumw2);
		h1_etRatio->Sumw2(isSumw2);

	};
	void Print(TString outName) {
		Sumw2(false);
		TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
		c1->Divide(2,2);
		c1->cd(1); h1_mass  ->Draw();  
		c1->cd(2); h1_mt     ->Draw(); 
		c1->cd(3); h1_dPhi   ->Draw(); 
		c1->cd(4); h1_etRatio->Draw(); 
		c1->Print(outName);
		delete c1;
		Sumw2(true);
	};
	void PrintAll(TString outName) {
		Print(outName);
		p4Hist1->Print(AddStringToFileName(outName, (TString)"_" + (TString)p4Hist1->GetDir()->GetName()));
		p4Hist2->Print(AddStringToFileName(outName, (TString)"_" + (TString)p4Hist2->GetDir()->GetName()));
	};
	TDirectory* GetDir() { return dir;};
};

class ElectronHist {
public:
	TDirectory* dir;
	TH1F* h1_scEta        ;
	TH1F* h1_fullSigma    ; 
	TH1F* h1_dEtaIn       ; 
	TH1F* h1_dPhiIn       ; 
	TH1F* h1_hOverE       ; 
	TH1F* h1_relIsoEA     ; 
	TH1F* h1_ooEmooP      ; 
	TH1F* h1_d0           ; 
	TH1F* h1_dz           ; 
	TH1F* h1_missHit      ; 
	TH1F* h1_passConvVeto ; 
	P4Hist* p4Hist;
	ElectronHist(TDirectory* mDir, TString dirName, bool isBarrel = true) {
		dir = mDir->mkdir(dirName);
		p4Hist = new P4Hist(dir, "p4Hist");
		dir->cd();
		if(isBarrel) {
			h1_scEta        = new TH1F("h1_scEta"       ,"h1_scEta"       ,  50, -2.5     , 2.5     );
			h1_fullSigma    = new TH1F("h1_fullSigma"   ,"h1_fullSigma"   , 100, 0        , 0.0101  );
			h1_dEtaIn       = new TH1F("h1_dEtaIn"      ,"h1_dEtaIn"      , 100, -0.00926 , 0.00926 );
			h1_dPhiIn       = new TH1F("h1_dPhiIn"      ,"h1_dPhiIn"      , 100, -0.0336  , 0.0336  );
			h1_hOverE       = new TH1F("h1_hOverE"      ,"h1_hOverE"      , 100, 0        , 0.0597  );
			h1_relIsoEA     = new TH1F("h1_relIsoEA"    ,"h1_relIsoEA"    , 100, 0        , 0.0354  );
			h1_ooEmooP      = new TH1F("h1_ooEmooP"     ,"h1_ooEmooP"     , 100, -0.012   , 0.012   );
			h1_d0           = new TH1F("h1_d0"          ,"h1_d0"          , 100, -0.0111  , 0.0111  );
			h1_dz           = new TH1F("h1_dz"          ,"h1_dz"          , 100, -0.0466  , 0.0466  );
			h1_missHit      = new TH1F("h1_missHit"     ,"h1_missHit"     , 4  , 0        , 3       );
			h1_passConvVeto = new TH1F("h1_passConvVeto","h1_passConvVeto", 4  , 0        , 3       );
		} else {
			h1_scEta        = new TH1F("h1_scEta"       ,"h1_scEta"       ,  50, -2.5     , 2.5     );
			h1_fullSigma    = new TH1F("h1_fullSigma"   ,"h1_fullSigma"   , 100, 0        , 0.0279  );
			h1_dEtaIn       = new TH1F("h1_dEtaIn"      ,"h1_dEtaIn"      , 100, -0.00724 , 0.00724 );
			h1_dPhiIn       = new TH1F("h1_dPhiIn"      ,"h1_dPhiIn"      , 100, -0.0918  , 0.0918  );
			h1_hOverE       = new TH1F("h1_hOverE"      ,"h1_hOverE"      , 100, 0        , 0.0615  );
			h1_relIsoEA     = new TH1F("h1_relIsoEA"    ,"h1_relIsoEA"    , 100, 0        , 0.0646  );
			h1_ooEmooP      = new TH1F("h1_ooEmooP"     ,"h1_ooEmooP"     , 100, -0.00999 , 0.00999 );
			h1_d0           = new TH1F("h1_d0"          ,"h1_d0"          , 100, -0.0351  , 0.0351  );
			h1_dz           = new TH1F("h1_dz"          ,"h1_dz"          , 100, -0.417   , 0.417   );
			h1_missHit      = new TH1F("h1_missHit"     ,"h1_missHit"     , 4  , 0        , 3       );
			h1_passConvVeto = new TH1F("h1_passConvVeto","h1_passConvVeto", 4  , 0        , 3       );
		}
	};

	void Fill(Electron& ele, double weight = 1.0, bool is50ns = false) {
		h1_scEta        ->Fill(ele.superCluster_eta               , weight) ;
		h1_fullSigma    ->Fill(ele.full5x5_sigmaIetaIeta          , weight) ; 
		h1_dEtaIn       ->Fill(ele.deltaEtaSuperClusterTrackAtVtx , weight) ; 
		h1_dPhiIn       ->Fill(ele.deltaPhiSuperClusterTrackAtVtx , weight) ; 
		h1_hOverE       ->Fill(ele.hadronicOverEm                 , weight) ; 
		h1_relIsoEA     ->Fill(ele.relIsoWithEA(is50ns)           , weight) ; 
		h1_ooEmooP      ->Fill(ele.ooEmooP                        , weight) ; 
		h1_d0           ->Fill(ele.gsfTrack_dxyPVtx               , weight) ; 
		h1_dz           ->Fill(ele.gsfTrack_dzPVtx                , weight) ; 
		h1_missHit      ->Fill(ele.nMissingHits                   , weight) ; 
		h1_passConvVeto ->Fill(ele.passConversionVeto             , weight) ; 
		p4Hist->Fill(ele, weight);
	};
	void Fill(Electron* elePtr, double weight = 1.0, bool is50ns = false) { Fill(*elePtr, weight, is50ns); };
	void Write() {
		dir->cd();
      h1_scEta        ->Write();  
      h1_fullSigma    ->Write();  
      h1_dEtaIn       ->Write();  
      h1_dPhiIn       ->Write();  
      h1_hOverE       ->Write();  
      h1_relIsoEA     ->Write();  
      h1_ooEmooP      ->Write();  
      h1_d0           ->Write();  
      h1_dz           ->Write();  
      h1_missHit      ->Write();  
      h1_passConvVeto ->Write();  
      p4Hist          ->Write();  
	};
	TDirectory* GetDir() { return dir; };
};

class MuonHist   {
public:
        TDirectory* dir;
        TH1F* h1_pt                   ;
        TH1F* h1_pt_BT                   ;
        TH1F* h1_pt_GB                   ;
        TH1F* h1_pt_IN                   ;
        TH1F* h1_eta                   ;
        TH1F* h1_phi                   ;
    	TH1F* h1_nValidPixelHits              ;
    	TH1F* h1_nTrackerLayers          ;
    	TH1F* h1_nValidMuonHits          ;
    	TH1F* h1_nMatchedStations ;
    	TH1F* h1_dxy                     ;    
    	TH1F* h1_dz                      ;
    	TH1F* h1_dxy_GV                  ;    
    	TH1F* h1_dz_GV                   ;       
    	TH1F* h1_dxy_GV_BT                  ;    
    	TH1F* h1_dz_GV_BT                   ;       
    	TH1F* h1_relTrkIsoR03_TP          ;
    	TH1F* h1_relTrkIsoR03_BT          ;
    	TH1F* h1_TrkIsoR03_sumPt                ; 
    	TH1F* h1_isHighPtMuonGV          ;     
    	TH1F* h1_statusGV                ;    

    	  P4Hist* p4Hist;
    	  MuonHist(TDirectory* mDir, TString dirName, bool isMB = true) {
    	 	dir = mDir->mkdir(dirName);
    	 	p4Hist = new P4Hist(dir, "p4Hist");
    	 	dir->cd();
    	 	if(isMB) {
                   h1_pt                     = new TH1F("h1_pt",                      "h1_pt", 500,0,5000 )   ; 
                   h1_pt_BT                     = new TH1F("h1_pt_BT",                      "h1_pt_BT", 500,0,5000 )   ; 
                   h1_pt_GB                     = new TH1F("h1_pt_GB",                      "h1_pt_GB", 500,0,5000 )   ; 
                   h1_pt_IN                     = new TH1F("h1_pt_IN",                      "h1_pt_IN", 500,0,5000 )   ; 
                   h1_eta                    = new TH1F("h1_eta",                     "h1_eta", 50,-2.5,2.5 )   ; 
                   h1_phi                    = new TH1F("h1_phi",                     "h1_phi", 50,-3.2,3.2 )   ; 
                   h1_nValidPixelHits             = new TH1F("h1_nValidPixelHits",              "h1_nValidPixelHits",20,0,20)    ; 
                   h1_nTrackerLayers         = new TH1F("h1_nTrackerLayers",          "h1_nTrackerLayers",50,0,50); 
                   h1_nValidMuonHits         = new TH1F("h1_nValidMuonHits",          "h1_nValidMuonHits",60,0,60); 
                   h1_nMatchedStations= new TH1F("h1_nMatchedStations", "h1_nMatchedStations",20,0,20) ; 
                   h1_dxy                    = new TH1F("h1_dxy",                     "h1_dxy",   100, -5,5)      ;     
                   h1_dz                     = new TH1F("h1_dz",                      "h1_dz",    100, -50,50)      ; 
                   h1_dxy_GV                 = new TH1F("h1_dxy_GV",                  "h1_dxy_GV",100, -5,5)      ;     
                   h1_dz_GV                  = new TH1F("h1_dz_GV",                   "h1_dz_GV", 100, -50,50)      ;      
                   h1_dxy_GV_BT                 = new TH1F("h1_dxy_GV_BT",                  "h1_dxy_GV_BT",100, -5,5)      ;     
                   h1_dz_GV_BT                  = new TH1F("h1_dz_GV_BT",                   "h1_dz_GV_BT", 100, -50,50)      ;      
                   h1_relTrkIsoR03_TP         = new TH1F("h1_relTrkIsoR03_TP",          "h1_relTrkIsoR03_TP",1000,0,1)  ; 
                   h1_relTrkIsoR03_BT         = new TH1F("h1_relTrkIsoR03_BT",          "h1_relTrkIsoR03_BT",1000,0,1)  ; 
                   h1_TrkIsoR03_sumPt               = new TH1F("h1_TrkIsoR03_sumPt",                "h1_TrkIsoR03_sumPt",      500,0,1000)  ;  
    	           h1_isHighPtMuonGV         = new TH1F("h1_isHighPtMuonGV",          "h1_isHighPtMuonGV",4,0,3)     ;     
                   h1_statusGV               = new TH1F("h1_statusGV",                "h1_statusGV",240,0,120)       ;     
		  }else {
                   h1_pt                     = new TH1F("h1_pt",                      "h1_pt", 500,0,5000 )   ; 
                   h1_pt_BT                     = new TH1F("h1_pt_BT",                      "h1_pt_BT", 500,0,5000 )   ; 
                   h1_pt_GB                     = new TH1F("h1_pt_GB",                      "h1_pt_GB", 500,0,5000 )   ; 
                   h1_pt_IN                     = new TH1F("h1_pt_IN",                      "h1_pt_IN", 500,0,5000 )   ; 
                   h1_eta                    = new TH1F("h1_eta",                     "h1_eta", 50,-2.5,2.5 )   ; 
                   h1_phi                    = new TH1F("h1_phi",                     "h1_phi", 50,-3.2,3.2 )   ; 
                   h1_nValidPixelHits             = new TH1F("h1_nValidPixelHits",              "h1_nValidPixelHits",20,0,20)    ; 
                   h1_nTrackerLayers         = new TH1F("h1_nTrackerLayers",          "h1_nTrackerLayers",50,0,50); 
                   h1_nValidMuonHits         = new TH1F("h1_nValidMuonHits",          "h1_nValidMuonHits",60,0,60); 
                   h1_nMatchedStations= new TH1F("h1_nMatchedStations", "h1_nMatchedStations",20,0,20) ; 
                   h1_dxy                    = new TH1F("h1_dxy",                     "h1_dxy",   400, -2,2)      ;     
                   h1_dz                     = new TH1F("h1_dz",                      "h1_dz",    400, -2,2)      ; 
                   h1_dxy_GV                 = new TH1F("h1_dxy_GV",                  "h1_dxy_GV",400, -2,2)      ;     
                   h1_dz_GV                  = new TH1F("h1_dz_GV",                   "h1_dz_GV", 400, -2,2)      ;      
                   h1_dxy_GV_BT                 = new TH1F("h1_dxy_GV_BT",                  "h1_dxy_GV_BT",100, -5,5)      ;     
                   h1_dz_GV_BT                  = new TH1F("h1_dz_GV_BT",                   "h1_dz_GV_BT", 100, -50,50)      ;      
                   h1_relTrkIsoR03_TP         = new TH1F("h1_relTrkIsoR03_TP",          "h1_relTrkIsoR03_TP",1000,0,1)  ; 
                   h1_relTrkIsoR03_BT         = new TH1F("h1_relTrkIsoR03_BT",          "h1_relTrkIsoR03_BT",1000,0,1)  ; 
                   h1_TrkIsoR03_sumPt               = new TH1F("h1_TrkIsoR03_sumPt",                "h1_TrkIsoR03_sumPt",      500,0,1000)  ;  
    	           h1_isHighPtMuonGV         = new TH1F("h1_isHighPtMuonGV",          "h1_isHighPtMuonGV",4,0,3)     ;     
                   h1_statusGV               = new TH1F("h1_statusGV",                "h1_statusGV",240,0,120)       ;     
                  }
	};
	void Fill(Muon mu, double weight = 1.0) {
		 h1_pt                     ->Fill(mu.pt                   , weight) ;
		 h1_pt_BT                     ->Fill(mu.pt_BT                   , weight) ;
		 h1_pt_GB                     ->Fill(mu.pt_GB                   , weight) ;
		 h1_pt_IN                     ->Fill(mu.pt_IN                   , weight) ;
		 h1_eta                    ->Fill(mu.eta                   , weight) ;
		 h1_phi                    ->Fill(mu.phi                   , weight) ; 
		 h1_nValidPixelHits             ->Fill(mu.nValidPixelHits              , weight) ; 
		 h1_nTrackerLayers         ->Fill(mu.nTrackerLayers          , weight) ; 
		 h1_nValidMuonHits         ->Fill(mu.nValidMuonHits          , weight) ; 
		 h1_nMatchedStations->Fill(mu.nMatchedStations , weight) ; 
		 h1_dxy                    ->Fill(mu.dxy                     , weight) ; 
		 h1_dz                     ->Fill(mu.dz                      , weight) ; 
		 h1_dxy_GV                 ->Fill(mu.dxy_GV                  , weight) ;
		 h1_dz_GV                  ->Fill(mu.dz_GV                   , weight) ; 
		 h1_dxy_GV_BT                 ->Fill(mu.dxy_GV_BT                  , weight) ;
		 h1_dz_GV_BT                  ->Fill(mu.dz_GV_BT                   , weight) ; 
		 h1_relTrkIsoR03_TP         ->Fill(mu.relTrkIsoR03_TP          , weight) ; 
		 h1_relTrkIsoR03_BT         ->Fill(mu.relTrkIsoR03_BT          , weight) ; 
		 h1_TrkIsoR03_sumPt               ->Fill(mu.TrkIsoR03_sumPt                , weight) ; 
	         h1_isHighPtMuonGV         ->Fill(mu.isHighPtMuonGV          , weight) ;
		 h1_statusGV               ->Fill(mu.statusGV                , weight) ; 
		 p4Hist->Fill(mu, weight);
	};
	void Fill(Muon* muPtr, double weight = 1.0) { Fill(*muPtr, weight); };
	void Write() {
		dir->cd();
		 h1_pt                  ->Write()   ;
		 h1_pt_BT                  ->Write()   ;
		 h1_pt_GB                  ->Write()   ;
		 h1_pt_IN                  ->Write()   ;
		 h1_eta                  ->Write()   ;
		 h1_phi                  ->Write()   ; 
		 h1_nValidPixelHits             ->Write()   ; 
		 h1_nTrackerLayers         ->Write()   ; 
		 h1_nValidMuonHits         ->Write()   ; 
		 h1_nMatchedStations->Write()   ; 
		 h1_dxy                    ->Write()   ; 
		 h1_dz                     ->Write()   ; 
		 h1_dxy_GV                 ->Write()   ;
		 h1_dz_GV                  ->Write()   ; 
		 h1_dxy_GV_BT                 ->Write()   ;
		 h1_dz_GV_BT                  ->Write()   ; 
		 h1_relTrkIsoR03_TP         ->Write()   ; 
		 h1_relTrkIsoR03_BT         ->Write()   ; 
		 h1_TrkIsoR03_sumPt               ->Write()   ; 
	         h1_isHighPtMuonGV         ->Write()   ;
		 h1_statusGV               ->Write()   ; 
	  	 p4Hist                    ->Write()  ;
	};
	TDirectory* GetDir() { return dir; };
};


class METHist   {
public:
        TDirectory* dir;
        TH1F* h1_isCaloMET           ;
        TH1F* h1_isPFMET             ;
        TH1F* h1_isRecoMET           ;
    	TH1F* h1_sumEt               ;
        TH1F* h1_corPt_Type01                 ; 
        TH1F* h1_corPhi_Type01                ; 
        TH1F* h1_corSumEt_Type01              ; 
        TH1F* h1_corPt_TypeXY                 ; 
        TH1F* h1_corPhi_TypeXY                ; 
        TH1F* h1_corSumEt_TypeXY              ; 
        TH1F* h1_corPt_Type1XY                ; 
        TH1F* h1_corPhi_Type1XY               ; 
        TH1F* h1_corSumEt_Type1XY             ; 
        TH1F* h1_corPt_Type01XY               ; 
        TH1F* h1_corPhi_Type01XY              ; 
        TH1F* h1_corSumEt_Type01XY            ; 
        TH1F* h1_corPt_Type1Smear             ; 
        TH1F* h1_corPhi_Type1Smear            ; 
        TH1F* h1_corSumEt_Type1Smear          ; 
        TH1F* h1_corPt_Type01Smear            ; 
        TH1F* h1_corPhi_Type01Smear           ; 
        TH1F* h1_corSumEt_Type01Smear         ; 
        TH1F* h1_corPt_Type1SmearXY           ; 
        TH1F* h1_corPhi_Type1SmearXY          ; 
        TH1F* h1_corSumEt_Type1SmearXY        ; 
        TH1F* h1_corPt_Type01SmearXY          ; 
        TH1F* h1_corPhi_Type01SmearXY         ; 
        TH1F* h1_corSumEt_Type01SmearXY       ;
        TH1F* h1_shiftedPt_JetResUp           ;
        TH1F* h1_shiftedPt_JetResDown         ;
        TH1F* h1_shiftedPt_JetEnUp            ;
        TH1F* h1_shiftedPt_JetEnDown          ;
        TH1F* h1_shiftedPt_UnclusteredEnUp    ;
        TH1F* h1_shiftedPt_UnclusteredEnDown  ;
        TH1F* h1_shiftedPt_JetResUpSmear      ;
        TH1F* h1_shiftedPt_JetResDownSmear    ;
    	TH1F* h1_METPuppi_pt           ;
    	TH1F* h1_METPuppi_phi          ;
    	TH1F* h1_genMET_pt           ;
    	TH1F* h1_genMET_phi          ;
    	TH1F* h1_genMET_sumEt        ;
    	TH1F* h1_metSignificance     ;    
    	TH1F* h1_metFilter_badPFMuon    ;    
    	TH1F* h1_metFilter_badChCan     ;    
    	P4Hist* p4Hist;

    	  METHist(TDirectory* mDir, TString dirName, bool isMET = true) {
    	 	dir = mDir->mkdir(dirName);
    	 	p4Hist = new P4Hist(dir, "p4Hist");
    	 	dir->cd();
    	 	if(isMET) {
                   h1_isCaloMET          = new TH1F("isCaloMET          ",  "h1_isCaloMET          ",4 ,0 ,3) ; 
                   h1_isPFMET            = new TH1F("isPFMET            ",  "h1_isPFMET            ",4 ,0 ,3) ;
                   h1_isRecoMET          = new TH1F("isRecoMET          ",  "h1_isRecoMET          ",4 ,0 ,3) ;
                   h1_sumEt              = new TH1F("sumEt              ",  "h1_sumEt              ",500 ,0 ,5000) ;
                   h1_corPt_Type01                 = new TH1F("h1_corPt_Type01"                ,"corPt_Type01"       ,500 ,0 ,5000 ) ; 
                   h1_corPhi_Type01                = new TH1F("h1_corPhi_Type01"               ,"corPhi_Type01"      ,50 ,-3.2,3.2 ) ; 
                   h1_corSumEt_Type01              = new TH1F("h1_corSumEt_Type01"             ,"corSumEt_Type01"    ,500 ,0 ,5000 )              ; 
                   h1_corPt_TypeXY                 = new TH1F("h1_corPt_TypeXY"                ,"corPt_TypeXY"       ,500 ,0 ,5000  )             ; 
                   h1_corPhi_TypeXY                = new TH1F("h1_corPhi_TypeXY"               ,"corPhi_TypeXY"      ,50 ,-3.2,3.2 )            ; 
                   h1_corSumEt_TypeXY              = new TH1F("h1_corSumEt_TypeXY"             ,"corSumEt_TypeXY"    ,500 ,0 ,5000  )            ; 
                   h1_corPt_Type1XY                = new TH1F("h1_corPt_Type1XY"               ,"corPt_Type1XY"      ,500 ,0 ,5000  )            ; 
                   h1_corPhi_Type1XY               = new TH1F("h1_corPhi_Type1XY"              ,"corPhi_Type1XY"     ,50 ,-3.2,3.2 )            ; 
                   h1_corSumEt_Type1XY             = new TH1F("h1_corSumEt_Type1XY"            ,"corSumEt_Type1XY"   ,500 ,0 ,5000  )            ; 
                   h1_corPt_Type01XY               = new TH1F("h1_corPt_Type01XY"              ,"corPt_Type01XY"     ,500 ,0 ,5000  )            ; 
                   h1_corPhi_Type01XY              = new TH1F("h1_corPhi_Type01XY"             ,"corPhi_Type01XY"    ,50 ,-3.2,3.2 )           ; 
                   h1_corSumEt_Type01XY            = new TH1F("h1_corSumEt_Type01XY"           ,"corSumEt_Type01XY"  ,500 ,0 ,5000  )           ; 
                   h1_corPt_Type1Smear             = new TH1F("h1_corPt_Type1Smear"            ,"corPt_Type1Smear"   ,500 ,0 ,5000  )           ; 
                   h1_corPhi_Type1Smear            = new TH1F("h1_corPhi_Type1Smear"           ,"corPhi_Type1Smear"  ,50 ,-3.2,3.2 )           ; 
                   h1_corSumEt_Type1Smear          = new TH1F("h1_corSumEt_Type1Smear"         ,"corSumEt_Type1Smear",500 ,0 ,5000     )          ; 
                   h1_corPt_Type01Smear            = new TH1F("h1_corPt_Type01Smear"           ,"corPt_Type01Smear"    ,500 ,0 ,5000   )          ; 
                   h1_corPhi_Type01Smear           = new TH1F("h1_corPhi_Type01Smear"          ,"corPhi_Type01Smear"   ,50 ,-3.2,3.2  )        ; 
                   h1_corSumEt_Type01Smear         = new TH1F("h1_corSumEt_Type01Smear"        ,"corSumEt_Type01Smear" ,500 ,0 ,5000   )        ; 
                   h1_corPt_Type1SmearXY           = new TH1F("h1_corPt_Type1SmearXY"          ,"corPt_Type1SmearXY"   ,500 ,0 ,5000   )         ; 
                   h1_corPhi_Type1SmearXY          = new TH1F("h1_corPhi_Type1SmearXY"         ,"corPhi_Type1SmearXY"  ,50 ,-3.2,3.2  )         ; 
                   h1_corSumEt_Type1SmearXY        = new TH1F("h1_corSumEt_Type1SmearXY"       ,"corSumEt_Type1SmearXY",500 ,0 ,5000 )       ; 
                   h1_corPt_Type01SmearXY          = new TH1F("h1_corPt_Type01SmearXY"         ,"corPt_Type01SmearXY"   ,500 ,0 ,5000   )       ; 
                   h1_corPhi_Type01SmearXY         = new TH1F("h1_corPhi_Type01SmearXY"        ,"corPhi_Type01SmearXY"  ,50 ,-3.2,3.2  )       ; 
                   h1_corSumEt_Type01SmearXY       = new TH1F("h1_corSumEt_Type01SmearXY"      ,"corSumEt_Type01SmearXY",500 ,0 ,5000   )       ;
                   h1_shiftedPt_JetResUp           = new TH1F("h1_shiftedPt_JetResUp"          ,"shiftedPt_JetResUp"     ,500 ,0 ,5000  )        ;
                   h1_shiftedPt_JetResDown         = new TH1F("h1_shiftedPt_JetResDown"        ,"shiftedPt_JetResDown"   ,500 ,0 ,5000  )        ;
                   h1_shiftedPt_JetEnUp            = new TH1F("h1_shiftedPt_JetEnUp"           ,"shiftedPt_JetEnUp"     ,500 ,0 ,5000  )         ;
                   h1_shiftedPt_JetEnDown          = new TH1F("h1_shiftedPt_JetEnDown"         ,"shiftedPt_JetEnDown"   ,500 ,0 ,5000  )         ;
                   h1_shiftedPt_UnclusteredEnUp    = new TH1F("h1_shiftedPt_UnclusteredEnUp"   ,"shiftedPt_UnclusteredEnUp"   ,500 ,0 ,5000 )  ;
                   h1_shiftedPt_UnclusteredEnDown  = new TH1F("h1_shiftedPt_UnclusteredEnDown" ,"shiftedPt_UnclusteredEnDown" ,500 ,0 ,5000 )  ;
                   h1_shiftedPt_JetResUpSmear      = new TH1F("h1_shiftedPt_JetResUpSmear"     ,"shiftedPt_JetResUpSmear"    ,500 ,0 ,5000 )   ;
                   h1_shiftedPt_JetResDownSmear    = new TH1F("h1_shiftedPt_JetResDownSmear"   ,"shiftedPt_JetResDownSmear"  ,500 ,0 ,5000 )   ;
                   h1_METPuppi_pt          = new TH1F("METPuppi_pt          ",  "h1_METPuppi_pt          ",500 ,0 ,5000) ;
                   h1_METPuppi_phi         = new TH1F("METPuppi_phi         ",  "h1_METPuppi_phi         ",50,-3.2 ,3.2) ;
                   h1_genMET_pt          = new TH1F("genMET_pt          ",  "h1_genMET_pt          ",500 ,0 ,5000) ;
                   h1_genMET_phi         = new TH1F("genMET_phi         ",  "h1_genMET_phi         ",50,-3.2 ,3.2) ;
                   h1_genMET_sumEt       = new TH1F("genMET_sumEt       ",  "h1_genMET_sumEt       ",500 ,0 ,5000) ;
                   h1_metSignificance    = new TH1F("metSignificance    ",  "h1_metSignificance    ",100 ,-10 ,10  )  ;
                   h1_metFilter_badPFMuon    = new TH1F("metFilter_badPFMuon",  "h1_metFilter_badPFMuon",4 ,0 ,3 ) ;
                   h1_metFilter_badChCan    = new TH1F("metFilter_badChCan  ",  "h1_metFilter_badChCan ",4 ,0 ,3 ) ;
		  }
	};
	void Fill(MET met, double weight = 1.0) {
		 h1_isCaloMET          ->Fill(met.isCaloMET           , weight) ;
		 h1_isPFMET            ->Fill(met.isPFMET             , weight) ;
		 h1_isRecoMET          ->Fill(met.isRecoMET           , weight) ; 
		 h1_sumEt              ->Fill(met.sumEt               , weight) ; 
                 h1_corPt_Type01               ->Fill(met.corPt_Type01                 , weight)  ; 
                 h1_corPhi_Type01              ->Fill(met.corPhi_Type01                , weight)  ; 
                 h1_corSumEt_Type01            ->Fill(met.corSumEt_Type01              , weight)  ; 
                 h1_corPt_TypeXY               ->Fill(met.corPt_TypeXY                 , weight)  ; 
                 h1_corPhi_TypeXY              ->Fill(met.corPhi_TypeXY                , weight)  ; 
                 h1_corSumEt_TypeXY            ->Fill(met.corSumEt_TypeXY              , weight)  ; 
                 h1_corPt_Type1XY              ->Fill(met.corPt_Type1XY                , weight)  ; 
                 h1_corPhi_Type1XY             ->Fill(met.corPhi_Type1XY               , weight)  ; 
                 h1_corSumEt_Type1XY           ->Fill(met.corSumEt_Type1XY             , weight)  ; 
                 h1_corPt_Type01XY             ->Fill(met.corPt_Type01XY               , weight)  ; 
                 h1_corPhi_Type01XY            ->Fill(met.corPhi_Type01XY              , weight)  ; 
                 h1_corSumEt_Type01XY          ->Fill(met.corSumEt_Type01XY            , weight)  ; 
                 h1_corPt_Type1Smear           ->Fill(met.corPt_Type1Smear             , weight)  ; 
                 h1_corPhi_Type1Smear          ->Fill(met.corPhi_Type1Smear            , weight)  ; 
                 h1_corSumEt_Type1Smear        ->Fill(met.corSumEt_Type1Smear          , weight)  ; 
                 h1_corPt_Type01Smear          ->Fill(met.corPt_Type01Smear            , weight)  ; 
                 h1_corPhi_Type01Smear         ->Fill(met.corPhi_Type01Smear           , weight)  ; 
                 h1_corSumEt_Type01Smear       ->Fill(met.corSumEt_Type01Smear         , weight)  ; 
                 h1_corPt_Type1SmearXY         ->Fill(met.corPt_Type1SmearXY           , weight)  ; 
                 h1_corPhi_Type1SmearXY        ->Fill(met.corPhi_Type1SmearXY          , weight)  ; 
                 h1_corSumEt_Type1SmearXY      ->Fill(met.corSumEt_Type1SmearXY        , weight)  ; 
                 h1_corPt_Type01SmearXY        ->Fill(met.corPt_Type01SmearXY          , weight)  ; 
                 h1_corPhi_Type01SmearXY       ->Fill(met.corPhi_Type01SmearXY         , weight)  ; 
                 h1_corSumEt_Type01SmearXY     ->Fill(met.corSumEt_Type01SmearXY       , weight)  ;
                 h1_shiftedPt_JetResUp         ->Fill(met.shiftedPt_JetResUp           , weight)  ;
                 h1_shiftedPt_JetResDown       ->Fill(met.shiftedPt_JetResDown         , weight)  ;
                 h1_shiftedPt_JetEnUp          ->Fill(met.shiftedPt_JetEnUp            , weight)  ;
                 h1_shiftedPt_JetEnDown        ->Fill(met.shiftedPt_JetEnDown          , weight)  ;
                 h1_shiftedPt_UnclusteredEnUp  ->Fill(met.shiftedPt_UnclusteredEnUp    , weight)  ;
                 h1_shiftedPt_UnclusteredEnDown->Fill(met.shiftedPt_UnclusteredEnDown  , weight)  ;
                 h1_shiftedPt_JetResUpSmear    ->Fill(met.shiftedPt_JetResUpSmear      , weight)  ;
                 h1_shiftedPt_JetResDownSmear  ->Fill(met.shiftedPt_JetResDownSmear    , weight)  ;
		 h1_METPuppi_pt          ->Fill(met.METPuppi_pt           , weight) ; 
		 h1_METPuppi_phi         ->Fill(met.METPuppi_phi          , weight) ; 
		 h1_genMET_pt          ->Fill(met.genMET_pt           , weight) ; 
		 h1_genMET_phi         ->Fill(met.genMET_phi          , weight) ; 
		 h1_genMET_sumEt       ->Fill(met.genMET_sumEt        , weight) ; 
		 h1_metSignificance    ->Fill(met.metSignificance     , weight) ;
		 //h1_metFilter_badPFMuon->Fill(met.filterbadPFMuon     , weight) ;
		 //h1_metFilter_badChCan ->Fill(met.filterbadChCandidate, weight) ;
		 p4Hist->Fill(met, weight);
	};
	void Fill(MET* metPtr, double weight = 1.0) { Fill(*metPtr, weight); };
	void Write() {
		dir->cd();
		 h1_isCaloMET          ->Write()   ;
		 h1_isPFMET            ->Write()   ;
		 h1_isRecoMET          ->Write()   ; 
		 h1_sumEt              ->Write()   ; 
                 h1_corPt_Type01               ->Write()  ; 
                 h1_corPhi_Type01              ->Write()  ; 
                 h1_corSumEt_Type01            ->Write()  ; 
                 h1_corPt_TypeXY               ->Write()  ; 
                 h1_corPhi_TypeXY              ->Write()  ; 
                 h1_corSumEt_TypeXY            ->Write()  ; 
                 h1_corPt_Type1XY              ->Write()  ; 
                 h1_corPhi_Type1XY             ->Write()  ; 
                 h1_corSumEt_Type1XY           ->Write()  ; 
                 h1_corPt_Type01XY             ->Write()  ; 
                 h1_corPhi_Type01XY            ->Write()  ; 
                 h1_corSumEt_Type01XY          ->Write()  ; 
                 h1_corPt_Type1Smear           ->Write()  ; 
                 h1_corPhi_Type1Smear          ->Write()  ; 
                 h1_corSumEt_Type1Smear        ->Write()  ; 
                 h1_corPt_Type01Smear          ->Write()  ; 
                 h1_corPhi_Type01Smear         ->Write()  ; 
                 h1_corSumEt_Type01Smear       ->Write()  ; 
                 h1_corPt_Type1SmearXY         ->Write()  ; 
                 h1_corPhi_Type1SmearXY        ->Write()  ; 
                 h1_corSumEt_Type1SmearXY      ->Write()  ; 
                 h1_corPt_Type01SmearXY        ->Write()  ; 
                 h1_corPhi_Type01SmearXY       ->Write()  ; 
                 h1_corSumEt_Type01SmearXY     ->Write()  ;
                 h1_shiftedPt_JetResUp         ->Write()  ;
                 h1_shiftedPt_JetResDown       ->Write()  ;
                 h1_shiftedPt_JetEnUp          ->Write()  ;
                 h1_shiftedPt_JetEnDown        ->Write()  ;
                 h1_shiftedPt_UnclusteredEnUp  ->Write()  ;
                 h1_shiftedPt_UnclusteredEnDown->Write()  ;
                 h1_shiftedPt_JetResUpSmear    ->Write()  ;
                 h1_shiftedPt_JetResDownSmear  ->Write()  ;
		 h1_METPuppi_pt          ->Write()   ; 
		 h1_METPuppi_phi         ->Write()   ; 
		 h1_genMET_pt          ->Write()   ; 
		 h1_genMET_phi         ->Write()   ; 
		 h1_genMET_sumEt       ->Write()   ; 
		 h1_metSignificance    ->Write()   ;
		 //h1_metFilter_badPFMuon->Write()   ;
		 //h1_metFilter_badChCan ->Write()   ;
	  	 p4Hist                ->Write()  ;
	};
	TDirectory* GetDir() { return dir; };
};


class PhotonHist {
public:
	TDirectory* dir;
	TH1F* h1_scEta                    ;
	TH1F* h1_hadTowOverEm             ; 
	TH1F* h1_full5x5_sigmaIetaIeta    ; 
	TH1F* h1_isoChargedHadronsWithEA  ; 
	TH1F* h1_isoNeutralHadronsWithEA  ; 
	TH1F* h1_isoPhotonsWithEA         ;
	TH1F* h1_passEleVeto              ; 
	TH1F* h1_hasPixelSeed             ; 
	P4Hist* p4Hist;
	PhotonHist(TDirectory* mDir, TString dirName, bool isBarrel = true) {
		dir = mDir->mkdir(dirName);
		dir->cd();
		if(isBarrel) {
			h1_scEta                    = new TH1F("h1_scEta"                   ,"h1_scEta"                   ,  50, -2.5, 2.5    );
			h1_hadTowOverEm             = new TH1F("h1_hadTowOverEm"            ,"h1_hadTowOverEm"            , 100, 0   , 0.05   );
			h1_full5x5_sigmaIetaIeta    = new TH1F("h1_full5x5_sigmaIetaIeta"   ,"h1_full5x5_sigmaIetaIeta"   , 100, 0   , 0.01   );
			h1_isoChargedHadronsWithEA  = new TH1F("h1_isoChargedHadronsWithEA" ,"h1_isoChargedHadronsWithEA" , 100, 0   , 0.76   );
			h1_isoNeutralHadronsWithEA  = new TH1F("h1_isoNeutralHadronsWithEA" ,"h1_isoNeutralHadronsWithEA" , 100, 0   , 1.0    );
			h1_isoPhotonsWithEA         = new TH1F("h1_isoPhotonsWithEA"        ,"h1_isoPhotonsWithEA"        , 100, 0   , 1.0    );
			h1_passEleVeto              = new TH1F("h1_passEleVeto"             ,"h1_passEleVeto"             ,   3, 0   , 2      );
			h1_hasPixelSeed             = new TH1F("h1_hasPixelSeed"            ,"h1_hasPixelSeed"            ,   3, 0   , 2      );
		} else {
			h1_scEta                    = new TH1F("h1_scEta"                   ,"h1_scEta"                   ,  50, -2.5, 2.5    );
			h1_hadTowOverEm             = new TH1F("h1_hadTowOverEm"            ,"h1_hadTowOverEm"            , 100, 0   , 0.05   );
			h1_full5x5_sigmaIetaIeta    = new TH1F("h1_full5x5_sigmaIetaIeta"   ,"h1_full5x5_sigmaIetaIeta"   , 100, 0   , 0.0268 );
			h1_isoChargedHadronsWithEA  = new TH1F("h1_isoChargedHadronsWithEA" ,"h1_isoChargedHadronsWithEA" , 100, 0   , 0.56   );
			h1_isoNeutralHadronsWithEA  = new TH1F("h1_isoNeutralHadronsWithEA" ,"h1_isoNeutralHadronsWithEA" , 100, 0   , 1.0    );
			h1_isoPhotonsWithEA         = new TH1F("h1_isoPhotonsWithEA"        ,"h1_isoPhotonsWithEA"        , 100, 0   , 1.0    );
			h1_passEleVeto              = new TH1F("h1_passEleVeto"             ,"h1_passEleVeto"             ,   3, 0   , 2      );
			h1_hasPixelSeed             = new TH1F("h1_hasPixelSeed"            ,"h1_hasPixelSeed"            ,   3, 0   , 2      );
		}
		p4Hist = new P4Hist(dir, "p4Hist");
	};
	void Fill(Photon pho, double weight = 1.0) {
		float rel_PFneuIso   = pho.isoNeutralHadronsWithEA / pho.CutThresholer("PFNeuIsoEA"   , 2, false) ;  
		float rel_PFphoIso   = pho.isoPhotonsWithEA        / pho.CutThresholer("PFPhoIsoEA"   , 2, false) ;  
		h1_scEta                    ->Fill(pho.superCluster_eta        , weight) ;
		h1_hadTowOverEm             ->Fill(pho.hadTowOverEm            , weight) ; 
		h1_full5x5_sigmaIetaIeta    ->Fill(pho.full5x5_sigmaIetaIeta   , weight) ; 
		h1_isoChargedHadronsWithEA  ->Fill(pho.isoChargedHadronsWithEA , weight) ; 
		h1_isoNeutralHadronsWithEA  ->Fill(rel_PFneuIso                , weight) ; 
		h1_isoPhotonsWithEA         ->Fill(rel_PFphoIso                , weight) ; 
		h1_passEleVeto              ->Fill(pho.passEleVeto             , weight) ; 
		h1_hasPixelSeed             ->Fill(pho.hasPixelSeed            , weight) ; 
		p4Hist->Fill(pho, weight);
	};
	void Fill(Photon* phoPtr, double weight = 1.0) { Fill(*phoPtr, weight); };
	void Write() {
		dir->cd();
		h1_scEta                    ->Write(); 
		h1_hadTowOverEm             ->Write(); 
		h1_full5x5_sigmaIetaIeta    ->Write(); 
		h1_isoChargedHadronsWithEA  ->Write(); 
		h1_isoNeutralHadronsWithEA  ->Write(); 
		h1_isoPhotonsWithEA         ->Write(); 
		h1_passEleVeto              ->Write(); 
		h1_hasPixelSeed             ->Write(); 
		p4Hist                      ->Write();
	};
	TDirectory* GetDir() { return dir; };
};

class EtcHist {
public:
	TDirectory* dir;
	//TH1F* h1_nVtx  ;
	//TH1F* h1_nGVtx ;
	TH1F* h1_pileupTrue ; 
	TH1F* h1_pileupPU ;
	TH1F* h1_puWeight ; 

	EtcHist(TDirectory* mDir, TString dirName) {
		dir = mDir->mkdir(dirName);
		dir->cd();
		//h1_nVtx  = new TH1F("h1_nVtx" ,"h1_nVtx" ,50,0,50);
		//h1_nGVtx = new TH1F("h1_nGVtx","h1_nGVtx",50,0,50);
		h1_pileupTrue = new TH1F("h1_pileupTrue","h1_pileupTrue",100,0,100);
		h1_pileupPU   = new TH1F("h1_pileupPU"  ,"h1_pileupPU"  ,100,0,100);
		h1_puWeight = new TH1F("h1_puWeight","h1_puWeight",10000,0,10);
	}
	~EtcHist() {};
	//void Fill(int nVtx, int nGVtx, double pileupTrue, int pileupPU, double weight = 1.0) {
	//	h1_nVtx ->Fill(nVtx, weight);
	//	h1_nGVtx->Fill(nGVtx, weight);
	//	h1_pileupTrue->Fill(pileupTrue, weight);
	//	h1_pileupPU  ->Fill(pileupPU, weight);
	//	h1_puWeight  ->Fill(weight) ;
	void Fill(double pileupTrue, int pileupPU, double weight = 1.0) {
		//h1_nVtx ->Fill(nVtx, weight);
		//h1_nGVtx->Fill(nGVtx, weight);
		h1_pileupTrue->Fill(pileupTrue, weight);
		h1_pileupPU  ->Fill(pileupPU, weight);
		h1_puWeight  ->Fill(weight) ;
	};
	void Write() {
		dir->cd();
		//h1_nVtx      ->Write() ; 
		//h1_nGVtx     ->Write() ; 
		h1_pileupTrue->Write() ;
		h1_pileupPU  ->Write() ;
		h1_puWeight  ->Write() ;
	};
	TDirectory* GetDir() {return dir;};
};


typedef std::vector<npknu::GenInfo>      GenInfoCollection       ;
typedef std::vector<npknu::PDFxfx>       PDFxfxCollection        ;
typedef std::vector<npknu::GenParticle>  GenParticleCollection   ;
typedef std::vector<npknu::Pileup>       PileupCollection        ;
typedef std::vector<npknu::Vertex>       VertexCollection        ;
typedef std::vector<npknu::Muon>         MuonCollection          ;
typedef std::vector<npknu::Electron>     ElectronCollection      ;
typedef std::vector<npknu::Photon>       PhotonCollection        ;
typedef std::vector<npknu::MET>          METCollection           ;
typedef std::vector<npknu::Jet>          JetCollection           ;
typedef std::vector<npknu::Trigger>      TriggerCollection       ;
typedef std::vector<npknu::TriggerObject> TriggerObjectCollection       ;


int GetElectronIdTCA(TClonesArray* inTCA, TClonesArray* outTCA, int idIndex=3, double minPt=0.0, double maxAEta = 100.0, bool doMCTruth = false);
int NumElectronIdTCA(TClonesArray* inTCA, int idIndex=3, double minPt=0.0, double maxAEta = 100.0, bool doMCTruth = false);
int GetPhotonIdTCA(TClonesArray* inTCA, TClonesArray* outTCA, int idIndex=2, double minPt=0.0, double maxAEta = 100.0, bool doMCTruth = false);
int NumPhotonIdTCA(TClonesArray* inTCA, int idIndex=2, double minPt=0.0, double maxAEta = 100.0, bool doMCTruth = false);
int GetMuonCutTCA(TClonesArray* inTCA, TClonesArray* outTCA, double minPt=0.0, double minAEta = 0.0, double maxAEta = 100.0);
int GetMuonIdTCA(TClonesArray* inTCA, TClonesArray* outTCA, int idIndex=3, double minPt=0.0, double minAEta = 0.0, double maxAEta = 100.0);
int NumMuonIdTCA(TClonesArray* inTCA, int idIndex=3, double minPt=0.0, double maxAEta = 100.0);
int GetGoodVertex(TClonesArray* inTCA, TClonesArray* outTCA);
int NumGoodVertex(TClonesArray* inTCA);
bool GetElectronPreselection(TClonesArray* inTCA, double minPt=0.0, double minAEta = 0.0, double maxAEta = 100.0);
bool GetMuonPreselection(TClonesArray* inTCA, double minPt=0.0, double minAEta = 0.0, double maxAEta = 100.0);
bool GetMuonVeto2nd(TClonesArray* inTCA,double minPt = 25.0);


int GetPassTriggerTCA(TClonesArray* triggerTCA, TClonesArray* outTCA);
int NumPassTriggerTCA(TClonesArray* triggerTCA);
int GetPassTriggerTCA(TClonesArray* triggerTCA, TClonesArray* outTCA, std::string name);
int NumPassTriggerTCA(TClonesArray* triggerTCA, std::string name);
int GetPassTriggerTCA(TClonesArray* triggerTCA, TClonesArray* outTCA, std::vector<std::string> nameVec);
bool isPassTriggerTCA(TClonesArray* triggerTCA, TClonesArray* outTCA, std::vector<std::string> nameVec);
int NumPassTriggerTCA(TClonesArray* triggerTCA, std::vector<std::string> nameVec);

int GetPassTriggerObjTCA(TClonesArray* triggerObjectTCA, TClonesArray* outTCA, std::string name = "");
bool isPassTriggerObjTCA(TClonesArray* triggerObjectTCA, TClonesArray* outTCA, std::string name = "");
int NumPassTriggerObjTCA(TClonesArray* triggerObjectTCA, std::string name ="");

int NumElectronMVATCA(TClonesArray* inTCA, unsigned int idBit = 3, double minPt=0.0, double maxAEta=100.0, bool doMCTruth =false); 
int GetElectronMVATCA(TClonesArray* inTCA, TClonesArray* outTCA, unsigned int idBit = 3, double minPt=0.0, double maxAEta=100.0, bool doMCTruth = false);

int GetDirEvents(TString fileNames="./*.root", TString treename = "MakeNpKNU/NpKNU" );


double DiLepMass(TClonesArray* lepTCA);
int GetPileupPU(TClonesArray* inTCA);
double GetPileupTrue(TClonesArray* inTCA);
void PrintCMD(int argc, char** argv);
void PrintElectronTCA(TClonesArray* electronTCA);
void PrintBit(unsigned int bit);

}; // End namespace npknu 

void WelcomeNpKNU();

#endif
