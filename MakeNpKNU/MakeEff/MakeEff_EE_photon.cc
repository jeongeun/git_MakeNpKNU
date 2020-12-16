#include <Riostream.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "./ProdNpKNU/src/NpKNU.hh"

   using namespace std;
   using namespace TMath;
   npknu::PUReweight* puWeight ;

   const char * ee_total[8] = {"ee_pt_tot", "ee_et_tot", "ee_eta_tot",  "ee_nPV_tot",  "ee_sig_tot" ,"ee_hoe_tot" ,"ee_deta_tot","ee_dphi_tot"} ; 
   const char * ee_pass1[8] = {"ee_pt_hlt0","ee_et_hlt0","ee_eta_hlt0", "ee_nPV_hlt0", "ee_sig_hlt0","ee_hoe_hlt0","ee_deta_hlt0","ee_dphi_hlt0"} ; 
   const char * ee_pass2[8] = {"ee_pt_hlt2","ee_et_hlt2","ee_eta_hlt2", "ee_nPV_hlt2", "ee_sig_hlt2","ee_hoe_hlt2","ee_deta_hlt2","ee_dphi_hlt2"} ; 

  const Int_t ptbins = 12;
  Double_t edges[ptbins+1] = {80, 90, 100, 120, 130, 140, 150, 200, 300, 400, 500, 800, 2000};

   Double_t  elebin[8]  =  {  20,  400,   16, 100, 50,   50 ,   200, 200 };
   Double_t  elexmin[8]  = { 0.0,  0.0, -2.5, 0.0, 0.0, 0.0 , -0.01,-0.1};
   Double_t  elexmax[8]  = {2000, 2000,  2.5, 100, 0.03,0.15,  0.01, 0.1};
   TString  elexTitle[8] = {"Electron p_{T} [GeV]", "Electron p_{T} [GeV]","Electron #eta_{SC}", "nPV", "#sigma_{i#etai#eta}", "H/E", "#delta#eta_{in}^{seed}", "#delta#phi_{in}"};

   TH1F* hee_total[8] ;  
   TH1F* hee_pass1[8] ;  
   TH1F* hee_pass2[8] ;  
        
int main(int argc, char* argv[]) {
     WelcomeNpKNU();
     int type;
     double lumi;
     double xsec;
     double gen ;
     type = atoi(argv[1]);
     lumi = atof(argv[2]);
     xsec = atof(argv[3]);
     gen  = atof(argv[4]);
     TString fileName = argv[5];
     double lumiWeight ;

     TChain* inChain = new TChain("MakeNpKNU/NpKNU");
     inChain->Add(argv[5]);
     TClonesArray* evtTCA         = new TClonesArray("npknu::Evt")            ;         inChain->SetBranchAddress("evt"          , &evtTCA         );
     TClonesArray* electronTCA    = new TClonesArray("npknu::Electron")       ;         inChain->SetBranchAddress("electron"     , &electronTCA    );
     TClonesArray* muonTCA        = new TClonesArray("npknu::Muon")           ;         inChain->SetBranchAddress("muon"         , &muonTCA        );
/*     TClonesArray* photonTCA      = new TClonesArray("npknu::Photon")         ;       inChain->SetBranchAddress("photon"       , &photonTCA      ); */
     TClonesArray* jetTCA         = new TClonesArray("npknu::Jet")            ;         inChain->SetBranchAddress("jet"          , &jetTCA         );
     TClonesArray* metTCA         = new TClonesArray("npknu::MET")            ;         inChain->SetBranchAddress("met"          , &metTCA         );
     TClonesArray* triggerTCA     = new TClonesArray("npknu::Trigger")        ;         inChain->SetBranchAddress("trigger"      , &triggerTCA     );
     TClonesArray* trigObjTCA     = new TClonesArray("npknu::TriggerObject")  ;         inChain->SetBranchAddress("triggerObject", &trigObjTCA     );
     TClonesArray* vertexTCA      = new TClonesArray("npknu::Vertex")         ;         inChain->SetBranchAddress("vertex"       , &vertexTCA      );
     TClonesArray* pileupTCA      = new TClonesArray("npknu::Pileup")         ;         inChain->SetBranchAddress("pileup"       , &pileupTCA      );
     TClonesArray* genInfoTCA     = new TClonesArray("npknu::GenInfo")        ;         inChain->SetBranchAddress("genInfo"      , &genInfoTCA     );
//    TClonesArray* genParticleTCA = new TClonesArray("npknu::GenParticle")    ;         inChain->SetBranchAddress("genParticle"  , &genParticleTCA );

     if(type == 1) {
        lumiWeight = 1.0;
     } else {
        lumiWeight = lumi *1000 * xsec / gen;
        puWeight = new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_FullRun_69p2mb.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
     }
     cout << "Type " << type << ": lumiWeight " << lumiWeight << "=" << lumi << "*" << xsec << "/" << gen << "  : " << fileName << endl;

     for(int j=0; j<8; j++) { //Draw Variables
        if(j==0){
           hee_total[j] = new TH1F(ee_total[j], ee_total[j], ptbins, edges )   ; hee_total[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   hee_total[j]->Sumw2();    
           hee_pass1[j] = new TH1F(ee_pass1[j], ee_pass1[j], ptbins, edges )   ; hee_pass1[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   hee_pass1[j]->Sumw2();    
           hee_pass2[j] = new TH1F(ee_pass2[j], ee_pass2[j], ptbins, edges )   ; hee_pass2[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   hee_pass2[j]->Sumw2();    
         }if(j>0){
         hee_total[j] = new TH1F(ee_total[j], ee_total[j], elebin[j], elexmin[j], elexmax[j]) ; hee_total[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   hee_total[j]->Sumw2();   
         hee_pass1[j] = new TH1F(ee_pass1[j], ee_pass1[j], elebin[j], elexmin[j], elexmax[j]) ; hee_pass1[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   hee_pass1[j]->Sumw2();   
         hee_pass2[j] = new TH1F(ee_pass2[j], ee_pass2[j], elebin[j], elexmin[j], elexmax[j]) ; hee_pass2[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   hee_pass2[j]->Sumw2();   
        }
     }

     TFile* outFile = new TFile(argv[6],"recreate");
     int TotalN = (int)inChain->GetEntries();
     int per99 = (TotalN > 10) ? TotalN / 9 : 1;
     int per100 = 0;

     std::vector<std::string> ReftrigNameVec;
     ReftrigNameVec.push_back("HLT_Mu*");
     ReftrigNameVec.push_back("HLT_IsoMu*");
     ReftrigNameVec.push_back("HLT_TkMu*");

     std::vector<std::string> SignalpathNameVec;
     SignalpathNameVec.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");  //reference HLT muon

     std::vector<std::string> PhotonNameVec;
     PhotonNameVec.push_back("HLT_Photon175_v*");

     TClonesArray* eleIdTCA         = new TClonesArray("npknu::Electron")       ;
     TClonesArray* referencePassTCA = new TClonesArray("npknu::Trigger");
     TClonesArray* signalPassTCA    = new TClonesArray("npknu::Trigger");
     TClonesArray* ECALPassTCA      = new TClonesArray("npknu::Trigger");
     TClonesArray* PhotonPassTCA    = new TClonesArray("npknu::Trigger");
	int denomenator =0;
	int nomenator =0;
	int nomenatorch =0;


     int NBadECALSlewden =0;
     int NDupECALden =0;
     int NunRepHitden =0;
     int NBadECALSlewnom =0;
     int NDupECALnom =0;
     int NunRepHitnom =0;
     int NBadECALSlewnomch =0;
     int NDupECALnomch =0;
     int NunRepHitnomch =0;

     for(int eventLoop=0; eventLoop<TotalN ; eventLoop++) {
        inChain->GetEntry(eventLoop);
        if((eventLoop%per99) == 0) cout << "Running " << (per100++ * 10)<< " % " << eventLoop << " / " << TotalN << endl;
        npknu::Evt* evtPtr = (npknu::Evt*)evtTCA->At(0);
        npknu::MET* metPtr = (npknu::MET*)metTCA->At(0);

        int NsignalPass =0;
        int NPhotonPass =0;
        int nVtx = 0;

        eleIdTCA->Clear();
        signalPassTCA->Clear();   referencePassTCA->Clear();
        PhotonPassTCA->Clear(); 

        bool isElePreselect = npknu::GetElectronPreselection(electronTCA, 35.0, 1.566, 2.5 ); /* ptCut, etaCut*/
        if(!isElePreselect) continue;

        npknu::GetPassTriggerTCA(triggerTCA, referencePassTCA, ReftrigNameVec);
        if(referencePassTCA->GetEntries() < 1) continue;

        npknu::GetElectronIdTCA(electronTCA, eleIdTCA, 4, 35.0, 2.5); /*HEEP, ptCut, etaCut*/
        nVtx = vertexTCA->GetEntries();

        if( eleIdTCA->GetEntries() != 1  ) continue; //HEEP 130Pt one Ele selection_2

        double pileup       = (evtPtr->isRealData) ? -999 : npknu::GetPileupTrue(pileupTCA);
        double pileupWeight = (evtPtr->isRealData) ?  1.0 : puWeight->getWeight(pileup);
        double tWeight      = pileupWeight * lumiWeight; // Notice:We need additional weight for gen failed event (Ngen/NGencut)

        npknu::Electron* heepPtr = (npknu::Electron*)eleIdTCA->At(0) ;
        if(heepPtr->isEE() != 1) continue;

        if(metPtr->IsDuplicateECALClusters == true){NDupECALden++;}
        if(metPtr->IsNotReplatedHitsECAL == true){NunRepHitden++;}
        if(metPtr->IsBadECALSlew == true){NBadECALSlewden++;}


        hee_total[0]->Fill( heepPtr->pt ,tWeight  );  
        hee_total[1]->Fill( heepPtr->pt ,tWeight  );  

        npknu::GetPassTriggerTCA(triggerTCA, signalPassTCA, SignalpathNameVec);
        npknu::GetPassTriggerTCA(triggerTCA, PhotonPassTCA, PhotonNameVec);

        bool isHLT115     = ( signalPassTCA->GetEntries() > 0) ;
        bool isPhotonpass = ( PhotonPassTCA->GetEntries() > 0) ;

        bool isBothpass   = (isHLT115 && isPhotonpass) ;
        bool isNopass     = (!isHLT115 && !isPhotonpass) ;

        if(isHLT115 ){ hee_pass1[0]->Fill( heepPtr->pt ,tWeight  );  hee_pass1[1]->Fill( heepPtr->pt ,tWeight  );  }
        if(isHLT115 || isPhotonpass){ hee_pass2[0]->Fill( heepPtr->pt ,tWeight  ); hee_pass2[1]->Fill( heepPtr->pt ,tWeight  ); }                            
        if(heepPtr->pt < 130.0) continue;
	denomenator++;
        hee_total[2]->Fill( heepPtr->superCluster_eta  ,tWeight  );  
        hee_total[3]->Fill( nVtx   ,tWeight );  
        hee_total[4]->Fill( heepPtr->full5x5_sigmaIetaIeta   ,tWeight );  
        hee_total[5]->Fill( heepPtr->hadronicOverEm    ,tWeight );  
        hee_total[6]->Fill( heepPtr->dEtaSeedAtVtx     ,tWeight );  
        hee_total[7]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx    ,tWeight );  
        if(isHLT115){
	    nomenator++;
            hee_pass1[2]->Fill( heepPtr->superCluster_eta  ,tWeight  );                 
            hee_pass1[3]->Fill( nVtx   ,tWeight );                                      
            hee_pass1[4]->Fill( heepPtr->full5x5_sigmaIetaIeta   ,tWeight );  
            hee_pass1[5]->Fill( heepPtr->hadronicOverEm    ,tWeight );  
            hee_pass1[6]->Fill( heepPtr->dEtaSeedAtVtx     ,tWeight );  
            hee_pass1[7]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );  

            if(metPtr->IsDuplicateECALClusters == true) {NDupECALnom++;}
            if(metPtr->IsNotReplatedHitsECAL == true) {NunRepHitnom++;}
            if(metPtr->IsBadECALSlew == true) {NBadECALSlewnom++;}

        }if(isHLT115 || isPhotonpass) {
	     nomenatorch++;
             hee_pass2[2]->Fill( heepPtr->superCluster_eta  ,tWeight  );                 
             hee_pass2[3]->Fill( nVtx   ,tWeight );                                      
             hee_pass2[4]->Fill( heepPtr->full5x5_sigmaIetaIeta   ,tWeight );  
             hee_pass2[5]->Fill( heepPtr->hadronicOverEm    ,tWeight );  
             hee_pass2[6]->Fill( heepPtr->dEtaSeedAtVtx     ,tWeight );  
             hee_pass2[7]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );  

            if(metPtr->IsDuplicateECALClusters == true) {NDupECALnomch++;}
            if(metPtr->IsNotReplatedHitsECAL == true) {NunRepHitnomch++;}
            if(metPtr->IsBadECALSlew == true) {NBadECALSlewnomch++;}
           }
        if(heepPtr->pt < 400.0) continue;
        if(!isHLT115 &&  isPhotonpass){
             cout << "@@ EE Only PhoPass run:lumi:evt " << evtPtr->run << ":" << evtPtr->lumi << ":" << evtPtr->event << " Ele115 " <<  isHLT115 << " Pho175 " << isPhotonpass << "===>";
	     heepPtr->PrintObj();
         }
        if(!isHLT115 && !isPhotonpass ){
             cout << "## EE All 2HLTFail run:lumi:evt " << evtPtr->run << ":" << evtPtr->lumi << ":" << evtPtr->event << " Ele115 " <<  isHLT115 << " Pho175 " << isPhotonpass << "===>";
	     heepPtr->PrintObj();
         }


  }
        delete eleIdTCA ;
        delete referencePassTCA ;
        delete signalPassTCA ;
//        delete ECALPassTCA ;
        delete PhotonPassTCA ;
        delete evtTCA          ;
        delete electronTCA     ;
        delete muonTCA         ;
        /*delete photonTCA       ;*/
        delete jetTCA          ;
        delete metTCA          ;
        delete triggerTCA      ;
        delete trigObjTCA      ;
        delete vertexTCA       ;
        delete pileupTCA       ;
        delete genInfoTCA      ;
//        delete genParticleTCA  ;

   hee_total[0]->Write() ;
   hee_pass1[0]->Write() ;
   hee_pass2[0]->Write() ;

   hee_total[1]->Write() ; 
   hee_pass1[1]->Write() ; 
   hee_pass2[1]->Write() ; 

   hee_total[2]->Write() ; 
   hee_pass1[2]->Write() ; 
   hee_pass2[2]->Write() ; 

   hee_total[3]->Write() ;  
   hee_pass1[3]->Write() ;  
   hee_pass2[3]->Write() ;  

   hee_total[4]->Write() ; 
   hee_pass1[4]->Write() ; 
   hee_pass2[4]->Write() ; 

   hee_total[5]->Write() ; 
   hee_pass1[5]->Write() ; 
   hee_pass2[5]->Write() ; 

   hee_total[6]->Write() ;  
   hee_pass1[6]->Write() ;  
   hee_pass2[6]->Write() ;  

   hee_total[7]->Write() ; 
   hee_pass1[7]->Write() ; 
   hee_pass2[7]->Write() ; 

   outFile->Write();

     cout << "1.(HEEP130) NBadECALSlewden  "<<   NBadECALSlewden  ;
     cout << " NDupECALden      "<<   NDupECALden      ;
     cout << " NunRepHitden     "<<   NunRepHitden     << endl;

     cout << "2.(Ele115)  NBadECALSlewnom  "<<   NBadECALSlewnom  ;
     cout << " NDupECALnom      "<<   NDupECALnom      ;
     cout << " NunRepHitnom     "<<   NunRepHitnom     <<  endl;

     cout << "3.(ECALpho) NBadECALSlewnomch"<<   NBadECALSlewnomch;
     cout << " NDupECALnomch    "<<   NDupECALnomch   ; 
     cout << " NunRepHitnomch   "<<   NunRepHitnomch   <<   endl;

        cout << "## " << fileName << " totalN " << TotalN <<endl;
        cout << " EE Total ele    = " << denomenator << "   " << int(hee_total[1]->Integral()) << endl; //eta
        cout << " EE HLT pass     = " << nomenator << "   " << int(hee_pass1[1]->Integral()) << endl;
        cout << " EE HLT+pho = " << nomenatorch << "   " <<int(hee_pass2[1]->Integral()) << endl;
        cout << "  # hlt/Tot = " <<  (double(nomenator)/double(denomenator)) << " %    " << ((double(hee_pass1[1]->Integral())/double(hee_total[1]->Integral()))*100) <<  " %" << endl;
        cout << "  # hlt+pho/Tot = " << (double(nomenatorch)/double(denomenator)) << " %    " <<  ((double(hee_pass2[1]->Integral())/double(hee_total[1]->Integral()))*100) <<  " %" << endl;

        outFile->Close();
        cout << "*** File Closed " << endl;

        return 0;
}