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
#include "./../ProdNpKNU/src/NpKNU.hh"

   using namespace std;
   using namespace TMath;
   npknu::PUReweight* puWeight ;

   const char * eb_total[9] = {"eb_pt_tot", "eb_eta_tot", "eb_phi_tot", "eb_sigieie_tot", "eb_hoe_tot", "eb_deta_tot", "eb_dphi_tot", "eb_trkIso_tot", "eb_nPV_tot"} ; 
   const char * eb_pass1[9] = {"eb_pt_hlt0", "eb_eta_hlt0", "eb_phi_hlt0", "eb_sigieie_hlt0", "eb_hoe_hlt0", "eb_deta_hlt0", "eb_dphi_hlt0","eb_trkIso_hlt0", "eb_nPV_hlt0"} ; 
   const char * eb_pass2[9] = {"eb_pt_hlt2", "eb_eta_hlt2", "eb_phi_hlt2", "eb_sigieie_hlt2", "eb_hoe_hlt2", "eb_deta_hlt2", "eb_dphi_hlt2",  "eb_trkIso_hlt2", "eb_nPV_hlt2"} ; 

   //const char * ee_total[9] = {"ee_pt_tot", "ee_eta_tot", "ee_phi_tot", "ee_sigieie_tot", "ee_hoe_tot", "ee_deta_tot", "ee_dphi_tot","ee_trkIso_tot", "ee_nPV_tot"} ; 
   //const char * ee_pass1[9] = {"ee_pt_hlt0", "ee_eta_hlt0", "ee_phi_hlt0", "ee_sigieie_hlt0", "ee_hoe_hlt0", "ee_deta_hlt0", "ee_dphi_hlt0",  "ee_trkIso_hlt0", "ee_nPV_hlt0"} ; 
   //const char * ee_pass2[9] = {"ee_pt_hlt2", "ee_eta_hlt2", "ee_phi_hlt2", "ee_sigieie_hlt2", "ee_hoe_hlt2", "ee_deta_hlt2", "ee_dphi_hlt2",  "ee_trkIso_hlt2", "ee_nPV_hlt2"} ; 

  const Int_t ptbins = 9;
  Double_t edges[ptbins+1] = {80, 90, 100, 130, 150, 300, 400, 500, 800, 2000};

   Double_t  elebin[9]  = { 20,    16,  16,   20,  20,     20,    20,    50, 100};
   Double_t  elexmin[9]  = { 0.0, -2.5,-3.2,  0.0, 0.0, -0.008, -0.08,   0.0, 0.0};
   Double_t  elexmax[9]  = {2000,  2.5, 3.2, 0.05, 0.1,  0.008,  0.08,   100, 100};
   TString  elexTitle[9] = {"Electron p_{T} [GeV]", "Electron #eta_{SC}", "Electron #phi_{SC}", "#sigma_{I#etaI#eta}", "H/E", "|#delta#eta_{In}^{seed}|","|#delta#phi_{in}|", "TrackerPtIsol", "nPV"};

   TH1F* heb_total[9] ;  
   TH1F* heb_pass1[9] ;  
   TH1F* heb_pass2[9] ;  
   //TH1F* hee_total[9] ;  
   //TH1F* hee_pass1[9] ;  
   //TH1F* hee_pass2[9] ;  
        
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
     TClonesArray* genParticleTCA = new TClonesArray("npknu::GenParticle")    ;         inChain->SetBranchAddress("genParticle"  , &genParticleTCA );

     if(type == 1) {
        lumiWeight = 1.0;
     } else {
        lumiWeight = lumi *1000 * xsec / gen;
        puWeight = new npknu::PUReweight("../etc/pileup/pu_Run_Full.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
     }
     cout << "Type " << type << ": lumiWeight " << lumiWeight << "=" << lumi << "*" << xsec << "/" << gen << "  : " << fileName << endl;

     for(int j=0; j<9; j++) { //Draw Variables
        if(j==0){
           heb_total[j] = new TH1F(eb_total[j], eb_total[j], ptbins, edges )   ; heb_total[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_total[j]->Sumw2();    
           heb_pass1[j] = new TH1F(eb_pass1[j], eb_pass1[j], ptbins, edges )   ; heb_pass1[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_pass1[j]->Sumw2();    
           heb_pass2[j] = new TH1F(eb_pass2[j], eb_pass2[j], ptbins, edges )   ; heb_pass2[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_pass2[j]->Sumw2();    
         }if(j>0){
         heb_total[j] = new TH1F(eb_total[j], eb_total[j], elebin[j], elexmin[j], elexmax[j]) ; heb_total[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_total[j]->Sumw2();   
         heb_pass1[j] = new TH1F(eb_pass1[j], eb_pass1[j], elebin[j], elexmin[j], elexmax[j]) ; heb_pass1[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_pass1[j]->Sumw2();   
         heb_pass2[j] = new TH1F(eb_pass2[j], eb_pass2[j], elebin[j], elexmin[j], elexmax[j]) ; heb_pass2[j]->GetXaxis()->SetTitle(elexTitle[j]) ;   heb_pass2[j]->Sumw2();   
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

     std::vector<std::string> ECALpathNameVec;
     ECALpathNameVec.push_back("HLT_ECALHT800_v*");
     std::vector<std::string> phopathNameVec;
     phopathNameVec.push_back("HLT_Photon175_v*");

        TClonesArray* eleIdTCA         = new TClonesArray("npknu::Electron")       ;
        TClonesArray* referencePassTCA = new TClonesArray("npknu::Trigger");
        TClonesArray* signalPassTCA    = new TClonesArray("npknu::Trigger");
        TClonesArray* ecalPassTCA      = new TClonesArray("npknu::Trigger");
        TClonesArray* photonPassTCA    = new TClonesArray("npknu::Trigger");

     for(int eventLoop=0; eventLoop<TotalN ; eventLoop++) {
        inChain->GetEntry(eventLoop);
        if((eventLoop%per99) == 0) cout << "Running " << (per100++ * 10)<< " % " << eventLoop << " / " << TotalN << endl;

        int NsignalPass =0;
        int NecalPass =0;
        int NphotonPass=0;
        int nVtx = 0;

        eleIdTCA->Clear();
        signalPassTCA->Clear();   referencePassTCA->Clear();
        ecalPassTCA->Clear();     photonPassTCA->Clear();

        npknu::Evt* evtPtr = (npknu::Evt*)evtTCA->At(0);
        double pileup       = (evtPtr->isRealData) ? -999 : npknu::GetPileupTrue(pileupTCA);
        double pileupWeight = (evtPtr->isRealData) ?  1.0 : puWeight->getWeight(pileup);
        double tWeight      = pileupWeight * lumiWeight; // Notice:We need additional weight for gen failed event (Ngen/NGencut)


        npknu::GetPassTriggerTCA(triggerTCA, referencePassTCA, ReftrigNameVec);

        bool isRefHLT = (referencePassTCA->GetEntries() > 0) ;

        npknu::GetElectronIdTCA(electronTCA, eleIdTCA, 4, 80.0, 2.5); /*HEEP, ptCut, etaCut*/

        NsignalPass = npknu::GetPassTriggerTCA(triggerTCA, signalPassTCA, SignalpathNameVec);
        bool isHLT115 = (NsignalPass > 0) ;

        NecalPass   = npknu::GetPassTriggerTCA(triggerTCA, ecalPassTCA, ECALpathNameVec);
        bool isECALpass = (NecalPass > 0) ;

        NphotonPass = npknu::GetPassTriggerTCA(triggerTCA, photonPassTCA, phopathNameVec);
        bool isphopass  = (NphotonPass > 0) ;

        nVtx = vertexTCA->GetEntries();

        if(!isRefHLT) continue;
        if( eleIdTCA->GetEntries() != 1  ) continue; //HEEP 130Pt one Ele selection_2

        npknu::Electron* heepPtr = (npknu::Electron*)eleIdTCA->At(0) ;
        if(heepPtr->isEB() != 1) continue;

        heb_total[0]->Fill( heepPtr->pt ,tWeight  );  
        if(isHLT115 ) heb_pass1[0]->Fill( heepPtr->pt ,tWeight  );             
        if(isHLT115 || isECALpass || isphopass) heb_pass2[0]->Fill( heepPtr->pt ,tWeight  );                                

        if(heepPtr->pt < 130.0) continue;
        heb_total[1]->Fill( heepPtr->superCluster_eta  ,tWeight  );  
        heb_total[2]->Fill( heepPtr->superCluster_phi   ,tWeight );  
        heb_total[3]->Fill( heepPtr->full5x5_sigmaIetaIeta  ,tWeight  );  
        heb_total[4]->Fill( heepPtr->hadronicOverEm   ,tWeight );  
        heb_total[5]->Fill( heepPtr->dEtaSeedAtVtx   ,tWeight );  
        heb_total[6]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );  
        heb_total[7]->Fill( heepPtr->trkIsol   ,tWeight );  
        heb_total[8]->Fill( nVtx   ,tWeight );  
        if(isHLT115){
            heb_pass1[1]->Fill( heepPtr->superCluster_eta  ,tWeight  );                 
            heb_pass1[2]->Fill( heepPtr->superCluster_phi   ,tWeight );                 
            heb_pass1[3]->Fill( heepPtr->full5x5_sigmaIetaIeta  ,tWeight  );            
            heb_pass1[4]->Fill( heepPtr->hadronicOverEm   ,tWeight );                   
            heb_pass1[5]->Fill( heepPtr->dEtaSeedAtVtx   ,tWeight );                    
            heb_pass1[6]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );   
            heb_pass1[7]->Fill( heepPtr->trkIsol   ,tWeight );                          
            heb_pass1[8]->Fill( nVtx   ,tWeight );                                      
        }
        if(isphopass || isECALpass || isHLT115 ) {
             heb_pass2[1]->Fill( heepPtr->superCluster_eta  ,tWeight  );                 
             heb_pass2[2]->Fill( heepPtr->superCluster_phi   ,tWeight );                 
             heb_pass2[3]->Fill( heepPtr->full5x5_sigmaIetaIeta  ,tWeight  );            
             heb_pass2[4]->Fill( heepPtr->hadronicOverEm   ,tWeight );                   
             heb_pass2[5]->Fill( heepPtr->dEtaSeedAtVtx   ,tWeight );                    
             heb_pass2[6]->Fill( heepPtr->deltaPhiSuperClusterTrackAtVtx   ,tWeight );   
             heb_pass2[7]->Fill( heepPtr->trkIsol   ,tWeight );                          
             heb_pass2[8]->Fill( nVtx   ,tWeight );                                      
           }
        if(!isHLT115 && (isECALpass || isphopass) ){
             cout << "@ EB HLT_Ele115 Fail run:lumi:evt " << evtPtr->run << ":" << evtPtr->lumi << ":" << evtPtr->event << " #HLT_Ele115 " <<  isHLT115 << " #HLT_ECALHT " <<  isECALpass  << " #HLT_photon175 " <<isphopass << endl;
	     heepPtr->PrintObj();
         }
        if(!isHLT115 && !isECALpass && !isphopass){
             cout << "@@@@@ EB All path Fail run:lumi:evt " << evtPtr->run << ":" << evtPtr->lumi << ":" << evtPtr->event << " #HLT_Ele115 " <<  isHLT115 << " #HLT_ECALHT " <<  isECALpass  << " #HLT_photon175 " <<isphopass << endl;
	     heepPtr->PrintObj();
         }


  }

   heb_total[0]->Write() ;    heb_total[1]->Write() ;    heb_total[2]->Write() ;    heb_total[3]->Write() ;  
   heb_pass1[0]->Write() ;    heb_pass1[1]->Write() ;    heb_pass1[2]->Write() ;    heb_pass1[3]->Write() ;  
   heb_pass2[0]->Write() ;    heb_pass2[1]->Write() ;    heb_pass2[2]->Write() ;    heb_pass2[3]->Write() ;  

   heb_total[4]->Write() ;    heb_total[5]->Write() ;    heb_total[6]->Write() ;    heb_total[7]->Write() ; 
   heb_pass1[4]->Write() ;    heb_pass1[5]->Write() ;    heb_pass1[6]->Write() ;    heb_pass1[7]->Write() ; 
   heb_pass2[4]->Write() ;    heb_pass2[5]->Write() ;    heb_pass2[6]->Write() ;    heb_pass2[7]->Write() ; 

  heb_total[8]->Write() ; 
  heb_pass1[8]->Write() ; 
  heb_pass2[8]->Write() ; 

          outFile->Write();


        cout << "## " << fileName << endl;
        cout << " EB Total ele   = " << heb_total[1]->Integral() << endl;
        cout << " EB HLT pass    = " << heb_pass1[1]->Integral() << endl;
        cout << " EB HLT+ecalpho = " << heb_pass2[1]->Integral() << endl;
        cout << "  # hlt/Tot = " <<  ( (heb_pass1[1]->Integral())/(heb_total[1]->Integral())*100  ) <<  " %" << endl;
        cout << "  # hlt+ecalpho/Tot = " <<  ( (heb_pass2[1]->Integral())/(heb_total[1]->Integral())*100  ) <<  " %" << endl;

  delete heb_total ; 
  delete heb_pass1 ; 
  delete heb_pass2 ; 

        delete eleIdTCA ;
        delete referencePassTCA ;
        delete signalPassTCA ;
        delete ecalPassTCA;
        delete photonPassTCA ;
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
        delete genParticleTCA  ;

        outFile->Close();
        cout << "*** File Closed " << endl;

        return 0;
}
