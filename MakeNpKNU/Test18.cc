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

int main(int argc, char* argv[]) {
     WelcomeNpKNU();
     TChain* inChain = new TChain("MakeNpKNU/NpKNU");
     inChain->Add(argv[1]);
     TClonesArray* evtTCA         = new TClonesArray("npknu::Evt")            ;         inChain->SetBranchAddress("evt"          , &evtTCA         );
     TClonesArray* electronTCA    = new TClonesArray("npknu::Electron")       ;         inChain->SetBranchAddress("electron"     , &electronTCA    );
     TClonesArray* triggerTCA     = new TClonesArray("npknu::Trigger")        ;         inChain->SetBranchAddress("trigger"      , &triggerTCA     );
     TClonesArray* trigObjTCA     = new TClonesArray("npknu::TriggerObject")  ;         inChain->SetBranchAddress("triggerObject", &trigObjTCA     );
     TClonesArray* vertexTCA      = new TClonesArray("npknu::Vertex")         ;         inChain->SetBranchAddress("vertex"       , &vertexTCA      );

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

        npknu::GetPassTriggerTCA(triggerTCA, referencePassTCA, ReftrigNameVec);
        if(referencePassTCA->GetEntries() < 1) continue;
        bool isElePreselect = npknu::GetElectronPreselection(electronTCA, 800.0, -2.5, 2.5 ); /* ptCut, etaCut*/
        if(!isElePreselect) continue;
        npknu::GetElectronIdTCA(electronTCA, eleIdTCA, 4, 35.0, 2.5); /*HEEP, ptCut, etaCut*/
        nVtx = vertexTCA->GetEntries();
        if( eleIdTCA->GetEntries() != 1  ) continue; //HEEP 130Pt one Ele selection_2
        npknu::Electron* heepPtr = (npknu::Electron*)eleIdTCA->At(0) ;
        if(heepPtr->pt < 800.0) continue;
        if(heepPtr->isEB() != 1) cout << "EE ";
        if(heepPtr->isEB() == 1) cout << "EB ";

        npknu::GetPassTriggerTCA(triggerTCA, signalPassTCA, SignalpathNameVec);
        npknu::GetPassTriggerTCA(triggerTCA, PhotonPassTCA, PhotonNameVec);

        bool isHLT115     = ( signalPassTCA->GetEntries() > 0) ;
        bool isPhotonpass = ( PhotonPassTCA->GetEntries() > 0) ;

        if(isHLT115 ){ cout << "pass Ele115 " << endl; }
        if(!isHLT115 ){ cout << "NOTpass Ele115 " << endl; }
        if(isPhotonpass){ cout << "pass Pho175 " << endl; }                            
        if(!isPhotonpass){ cout << "@@ NOTpass Pho175 " << endl; }                            
	
        cout << "  scEta=" << heepPtr->superCluster_eta ;  
        cout << "  nVtx=" << nVtx   ;  
        cout << "  sigInIn=" << heepPtr->full5x5_sigmaIetaIeta  ;  
        cout << "  hoe=" << heepPtr->hadronicOverEm   ;  
        cout << "  deta=" << heepPtr->dEtaSeedAtVtx    ;  
        cout << "  dphi=" << heepPtr->deltaPhiSuperClusterTrackAtVtx  << endl;;  
  }
        delete eleIdTCA ;
        delete referencePassTCA ;
        delete signalPassTCA ;
        delete PhotonPassTCA ;
        delete evtTCA          ;
        delete electronTCA     ;
        delete triggerTCA      ;
        delete trigObjTCA      ;
        delete vertexTCA       ;
        cout << "*** File Closed " << endl;

        return 0;
}
