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

   const char * eb_total[8] = {"eb_pt_tot", "eb_et_tot", "eb_eta_tot" ,  "eb_nPV_tot",  "eb_sig_tot" ,"eb_hoe_tot" ,"eb_deta_tot","eb_dphi_tot"} ; 
   const char * eb_passE[8] = {"eb_pt_hlt0","eb_et_hlt0","eb_eta_hlt0", "eb_nPV_hlt0", "eb_sig_hlt0","eb_hoe_hlt0","eb_deta_hlt0","eb_dphi_hlt0"} ; 
   const char * eb_passP[8] = {"eb_pt_hltp","eb_et_hltp","eb_eta_hltp", "eb_nPV_hltp", "eb_sig_hltp","eb_hoe_hltp","eb_deta_hltp","eb_dphi_hltp"} ; 

  const Int_t ptbins = 12;
//  Double_t edges[ptbins+1] = { 80, 90, 100, 120, 130, 140, 150, 200, 300, 400, 500, 800, 2000};
  Double_t edges[ptbins+1] = {150 ,155,160, 170, 180, 190, 200, 260, 300, 400, 500, 800, 2000};

   Double_t  elebin[8]  =  {  20,  400,   16, 100, 50,   50 ,   200, 200 };
   Double_t  elexmin[8]  = { 0.0,  0.0, -2.5, 0.0, 0.0, 0.0 , -0.01,-0.1};
   Double_t  elexmax[8]  = {2000, 2000,  2.5, 100, 0.03,0.15,  0.01, 0.1};
   TString  elexTitle[8] = {"Electron p_{T} [GeV]", "Electron p_{T} [GeV]","Electron #eta_{SC}", "nPV", "#sigma_{i#etai#eta}", "H/E", "#delta#eta_{in}^{seed}", "#delta#phi_{in}"};

   TH1F* heb_total[8] ;  
   TH1F* heb_passE[8] ;  
   TH1F* heb_passP[8] ;  
        
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
//     TClonesArray* muonTCA        = new TClonesArray("npknu::Muon")           ;         inChain->SetBranchAddress("muon"         , &muonTCA        );
///*     TClonesArray* photonTCA      = new TClonesArray("npknu::Photon")         ;       inChain->SetBranchAddress("photon"       , &photonTCA      ); */
//     TClonesArray* jetTCA         = new TClonesArray("npknu::Jet")            ;         inChain->SetBranchAddress("jet"          , &jetTCA         );
//     TClonesArray* metTCA         = new TClonesArray("npknu::MET")            ;         inChain->SetBranchAddress("met"          , &metTCA         );
     TClonesArray* triggerTCA     = new TClonesArray("npknu::Trigger")        ;         inChain->SetBranchAddress("trigger"      , &triggerTCA     );
     TClonesArray* trigObjTCA     = new TClonesArray("npknu::TriggerObject")  ;         inChain->SetBranchAddress("triggerObject", &trigObjTCA     );
     TClonesArray* vertexTCA      = new TClonesArray("npknu::Vertex")         ;         inChain->SetBranchAddress("vertex"       , &vertexTCA      );
     TClonesArray* pileupTCA      = new TClonesArray("npknu::Pileup")         ;         inChain->SetBranchAddress("pileup"       , &pileupTCA      );
     TClonesArray* genInfoTCA     = new TClonesArray("npknu::GenInfo")        ;         inChain->SetBranchAddress("genInfo"      , &genInfoTCA     );
//     TClonesArray* genParticleTCA = new TClonesArray("npknu::GenParticle")    ;         inChain->SetBranchAddress("genParticle"  , &genParticleTCA );

     if(type == 1) {
        lumiWeight = 1.0;
     } else {
        lumiWeight = lumi *1000 * xsec / gen;
        puWeight = new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_FullRun_69p2mb.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
     }
     cout << "Type " << type << ": lumiWeight " << lumiWeight << "=" << lumi << "*" << xsec << "/" << gen << "  : " << fileName << endl;

     for(int j=0; j<8; j++) { //Draw Variables
        if(j==0){
           heb_total[j] = new TH1F(eb_total[j], eb_total[j], ptbins, edges )   ; heb_total[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_total[j]->Sumw2();    
           heb_passE[j] = new TH1F(eb_passE[j], eb_passE[j], ptbins, edges )   ; heb_passE[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_passE[j]->Sumw2();    
           heb_passP[j] = new TH1F(eb_passP[j], eb_passP[j], ptbins, edges )   ; heb_passP[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_passP[j]->Sumw2();    
         }if(j>0){
         heb_total[j] = new TH1F(eb_total[j], eb_total[j], elebin[j], elexmin[j], elexmax[j]) ; heb_total[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_total[j]->Sumw2();   
         heb_passE[j] = new TH1F(eb_passE[j], eb_passE[j], elebin[j], elexmin[j], elexmax[j]) ; heb_passE[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_passE[j]->Sumw2();   
         heb_passP[j] = new TH1F(eb_passP[j], eb_passP[j], elebin[j], elexmin[j], elexmax[j]) ; heb_passP[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_passP[j]->Sumw2();   
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
     TClonesArray* PhotonPassTCA    = new TClonesArray("npknu::Trigger");
	int denomenator =0;
	int nomenator =0;
	int nomenatorch =0;

     for(int eventLoop=0; eventLoop<TotalN ; eventLoop++) {
        inChain->GetEntry(eventLoop);
        if((eventLoop%per99) == 0) cout << "Running " << (per100++ * 10)<< " % " << eventLoop << " / " << TotalN << endl;
        npknu::Evt* evtPtr = (npknu::Evt*)evtTCA->At(0);

        int NsignalPass =0;
        int NPhotonPass =0;
        int nVtx = 0;

        eleIdTCA->Clear();
        signalPassTCA->Clear();   referencePassTCA->Clear();
        PhotonPassTCA->Clear(); 

        bool isElePreselect = npknu::GetElectronPreselection(electronTCA, 80.0, 0.0, 1.4442 ); /* ptCut, etaCut*/
        if(!isElePreselect) continue;

        npknu::GetPassTriggerTCA(triggerTCA, referencePassTCA, ReftrigNameVec);
        if(referencePassTCA->GetEntries() < 1) continue;

        npknu::GetElectronIdTCA(electronTCA, eleIdTCA, 4, 35.0, 1.4442); /*HEEP, ptCut, etaCut*/
        nVtx = vertexTCA->GetEntries();

        if( eleIdTCA->GetEntries() != 1  ) continue; //HEEP 130Pt one Ele selection_2

        double pileup       = (evtPtr->isRealData) ? -999 : npknu::GetPileupTrue(pileupTCA);
        double pileupWeight = (evtPtr->isRealData) ?  1.0 : puWeight->getWeight(pileup);
        double tWeight      = pileupWeight * lumiWeight; // Notice:We need additional weight for gen failed event (Ngen/NGencut)

        npknu::Electron* heepPtr = (npknu::Electron*)eleIdTCA->At(0) ;
        if(heepPtr->isEB() != 1) continue;
        heb_total[0]->Fill( heepPtr->pt ,tWeight  );  
        heb_total[1]->Fill( heepPtr->pt ,tWeight  );  

        npknu::GetPassTriggerTCA(triggerTCA, signalPassTCA, SignalpathNameVec);
        npknu::GetPassTriggerTCA(triggerTCA, PhotonPassTCA, PhotonNameVec);

        bool isHLT115     = ( signalPassTCA->GetEntries() > 0) ;
        bool isPhotonpass = ( PhotonPassTCA->GetEntries() > 0) ;

        bool isBothpass   = (isHLT115 && isPhotonpass) ;
        bool isNopass     = (!isHLT115 && !isPhotonpass) ;

        if(isHLT115 ){ heb_passE[0]->Fill( heepPtr->pt ,tWeight  );  heb_passE[1]->Fill( heepPtr->pt ,tWeight  );  }
        if(isPhotonpass){ heb_passP[0]->Fill( heepPtr->pt ,tWeight  ); heb_passP[1]->Fill( heepPtr->pt ,tWeight  ); }                            
        if(heepPtr->pt < 130.0) continue;
	denomenator++;
        heb_total[2]->Fill( heepPtr->superCluster_eta  ,tWeight  );  
        heb_total[3]->Fill( nVtx   ,tWeight );  
        heb_total[4]->Fill( heepPtr->full5x5_sigmaIetaIeta   ,tWeight );  
        heb_total[5]->Fill( heepPtr->hadronicOverEm    ,tWeight );  
        heb_total[6]->Fill( heepPtr->dEtaSeedAtVtx     ,tWeight );  
        heb_total[7]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx    ,tWeight );  
        if(isHLT115){
	    nomenator++;
            heb_passE[2]->Fill( heepPtr->superCluster_eta  ,tWeight  );                 
            heb_passE[3]->Fill( nVtx   ,tWeight );                                      
            heb_passE[4]->Fill( heepPtr->full5x5_sigmaIetaIeta   ,tWeight );  
            heb_passE[5]->Fill( heepPtr->hadronicOverEm    ,tWeight );  
            heb_passE[6]->Fill( heepPtr->dEtaSeedAtVtx     ,tWeight );  
            heb_passE[7]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );  
        }
        if(isPhotonpass) {
	     nomenatorch++;
             heb_passP[2]->Fill( heepPtr->superCluster_eta  ,tWeight  );                 
             heb_passP[3]->Fill( nVtx   ,tWeight );                                      
             heb_passP[4]->Fill( heepPtr->full5x5_sigmaIetaIeta   ,tWeight );  
             heb_passP[5]->Fill( heepPtr->hadronicOverEm    ,tWeight );  
             heb_passP[6]->Fill( heepPtr->dEtaSeedAtVtx     ,tWeight );  
             heb_passP[7]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );  
           }
  }
        delete eleIdTCA ;
        delete referencePassTCA ;
        delete signalPassTCA ;
        delete PhotonPassTCA ;
        delete evtTCA          ;
        delete electronTCA     ;
        //delete muonTCA         ;
        //delete jetTCA          ;
        //delete metTCA          ;
        delete triggerTCA      ;
        delete trigObjTCA      ;
        delete vertexTCA       ;
        delete pileupTCA       ;
        delete genInfoTCA      ;
        //delete genParticleTCA  ;

   heb_total[0]->Write() ;
   heb_passE[0]->Write() ;
   heb_passP[0]->Write() ;

   heb_total[1]->Write() ; 
   heb_passE[1]->Write() ; 
   heb_passP[1]->Write() ; 

   heb_total[2]->Write() ; 
   heb_passE[2]->Write() ; 
   heb_passP[2]->Write() ; 

   heb_total[3]->Write() ;
   heb_passE[3]->Write() ;
   heb_passP[3]->Write() ;

   heb_total[4]->Write() ; 
   heb_passE[4]->Write() ; 
   heb_passP[4]->Write() ; 

   heb_total[5]->Write() ; 
   heb_passE[5]->Write() ; 
   heb_passP[5]->Write() ; 

   heb_total[6]->Write() ;
   heb_passE[6]->Write() ;
   heb_passP[6]->Write() ;

   heb_total[7]->Write() ; 
   heb_passE[7]->Write() ; 
   heb_passP[7]->Write() ; 

    outFile->Write();


        cout << "## " << fileName << " totalN " << TotalN <<endl;
        cout << " EB Total ele    = " << denomenator << "   " << int(heb_total[1]->Integral()) << endl; //eta
        cout << " EB Ele115 pass     = " << nomenator << "   " << int(heb_passE[1]->Integral()) << endl;
        cout << " EB Pho175 pass = " << nomenatorch << "   " <<int(heb_passP[1]->Integral()) << endl;
        cout << "  # Ele115/Tot = " <<  (double(nomenator)/double(denomenator)) << " %    " << ((double(heb_passE[1]->Integral())/double(heb_total[1]->Integral()))*100) <<  " %" << endl;
        cout << "  # Pho175/Tot = " << (double(nomenatorch)/double(denomenator)) << " %    " <<  ((double(heb_passP[1]->Integral())/double(heb_total[1]->Integral()))*100) <<  " %" << endl;

        outFile->Close();
        cout << "*** File Closed " << endl;

        return 0;
}
