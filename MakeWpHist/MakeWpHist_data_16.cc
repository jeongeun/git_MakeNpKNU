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
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ProdNpKNU/src/NpKNU.hh"

   using namespace std;
   using namespace TMath;
   npknu::PUReweight* puWeight ;

   const char *  elestr_pre[9] = {"ele_EBPt_pre","ele_EEPt_pre","ele_Pt_pre", "ele_Eta_pre","ele_Phi_pre","ele_Sig_pre", "ele_HoE_pre", "ele_dEta_pre", "ele_dPhi_pre" };
   const char *  elestr_kin[9] = {"ele_EBPt_kin","ele_EEPt_kin","ele_Pt_kin", "ele_Eta_kin","ele_Phi_kin","ele_Sig_kin", "ele_HoE_kin", "ele_dEta_kin", "ele_dPhi_kin" };

   Int_t elebin[9]      = { 80,  80,  80,    16,  16,  50,  50,   100,   100 };
   Double_t elexmin[9]  = { 0.0, 0.0, 0.0, -2.5,-3.2, 0.0,  0.0, -0.01, -0.1 };
   Double_t elexmax[9]  = {2000,2000,2000,  2.5, 3.2, 0.05, 0.2,  0.01,  0.1 };
   TString elexTitle[9] = {"EB Electron p_{T} [GeV]", "EE Electron p_{T} [GeV]", "Electron p_{T} [GeV]", "Electron #eta_{SC}", "Electron #phi_{SC}","#sigma_{i#etai#eta}", "H/E", "#delta#eta_{in}^{seed}", "#delta#phi_{in}"};

   const char * metstr_pre[2] = {"met_Pt_pre","met_Phi_pre"}  ;
   const char * metstr_kin[2] = {"met_Pt_kin","met_Phi_kin"}  ;

   Int_t metbin[2]      = { 80,   50}  ;
   Double_t metxmin[2]  = { 0.0,-3.2}  ;
   Double_t metxmax[2]  = {2000, 3.2}  ;
   TString metxTitle[2] = {"MET [GeV]", "MET #phi"}  ;

   const char *  etcstr_pre[6] = {"invM_pre","nPV_pre","nPVw_pre","wp_Mt_pre","wp_Dphi_pre","wp_EtRatio_pre"}  ;
   const char *  etcstr_kin[6] = {"invM_kin","nPV_kin","nPVw_kin","wp_Mt_kin","wp_Dphi_kin","wp_EtRatio_kin"}  ;

   Int_t etcbin[6]      ={1000, 100,  100,   100,  50, 50 }   ;
   Double_t etcxmin[6]  ={500.0, 0.0,  0.0,  0.0, 0.0,0.0 }  ;
   Double_t etcxmax[6]  ={1500, 100,  100,  4000, 5.0, 10 }  ;
   TString etcxTitle[6] ={"inv M_{#ele#ele}","nPV(w/o PUreWeight)","nPV(w PUreweight)" ,"M_{T} [GeV]","#Delta #phi","E_{T}^{#ele}/E_{T}^{miss}"}  ;

   TH1F* hele_pre[9] ;
   TH1F* hele_kin[9] ;

   TH1F* hmet_pre[2] ;
   TH1F* hmet_kin[2] ;

   TH1F* hetc_pre[6] ;
   TH1F* hetc_kin[6] ;

double GetIntegral(THStack * st) {
        TList * histKeys = st->GetHists();
        histKeys->Print();
        TIter next(histKeys);
        TObject* object = 0;
        double sum = 0;
        while ((object = next())) sum += ((TH1*)object)->Integral();

        return sum;
}

/*Here we need to make scalefactor apply part*/
        
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
   
//     bool isEB =  true;
//     bool isEE =  false; //false;


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
     TClonesArray* eleIdTCA        = new TClonesArray("npknu::Electron")           ;

     if(type == 1) {
        lumiWeight = 1.0;
     } else {
        lumiWeight = lumi *1000 * xsec / gen;
        //puWeight = 1;//new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_RunBtoF_69p2mb.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
       // puWeight = new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_RunGtoH_69p2mb.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
     }
     cout << "Type " << type << ": lumiWeight " << lumiWeight << "=" << lumi << "*" << xsec << "/" << gen << "  : " << fileName << endl;


     for(int j=0; j<9; j++){ //Draw Variables
        hele_pre[j] = new TH1F(elestr_pre[j],elestr_pre[j], elebin[j], elexmin[j], elexmax[j]); hele_pre[j]->GetXaxis()->SetTitle(elexTitle[j]);  hele_pre[j]->Sumw2();
        hele_kin[j] = new TH1F(elestr_kin[j],elestr_kin[j], elebin[j], elexmin[j], elexmax[j]); hele_kin[j]->GetXaxis()->SetTitle(elexTitle[j]);  hele_kin[j]->Sumw2();
     }
     for(int k=0; k<2; k++){ //Draw Variables
        hmet_pre[k] = new TH1F(metstr_pre[k],metstr_pre[k], metbin[k], metxmin[k], metxmax[k]); hmet_pre[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_pre[k]->Sumw2();
        hmet_kin[k] = new TH1F(metstr_kin[k],metstr_kin[k], metbin[k], metxmin[k], metxmax[k]); hmet_kin[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_kin[k]->Sumw2();
     }
     for(int p=0; p<6; p++){ //Draw Variables
        hetc_pre[p] = new TH1F(etcstr_pre[p],etcstr_pre[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_pre[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_pre[p]->Sumw2();
        hetc_kin[p] = new TH1F(etcstr_kin[p],etcstr_kin[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_kin[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_kin[p]->Sumw2();
     }

     TFile* outFile = new TFile(argv[6],"recreate");
     int TotalN = (int)inChain->GetEntries();
     int per99 = (TotalN > 10) ? TotalN / 9 : 1;
     int per100 = 0;

     std::vector<std::string> trigNameVec;

     trigNameVec.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
     trigNameVec.push_back("HLT_Ele135_CaloIdVT_GsfTrkIdT_v*");
     //trigNameVec.push_back("HLT_Ele35_WPTight_Gsf_v*");
     //trigNameVec.push_back("HLT_Ele38_WPTight_Gsf_v*");
     trigNameVec.push_back("HLT_Photon175_v*");
     TClonesArray* trigNamePassTCA = new TClonesArray("npknu::Trigger");

//   bool isMETFilter = false;
   bool isMETAllFilter = false;
   bool isElePreselec = false;

   int Nmetfil            = 0 ; 
   int Nmetallfil         = 0 ; 
   int NpassMET           = 0 ; 
   int Ntrigg             = 0 ; 
   int NHEEP              = 0 ; 
   int NHEEP2             = 0 ; 
   int NHEEP1Veto25       = 0 ; 
   int Nhas2ndMuon        = 0 ;
   int NPreSelection      = 0 ; 
   int NKinCutPass        = 0 ; 
   int NHighMT            = 0 ; 
   int NHighMT1TeV        = 0 ; 
   int NHighMT2TeV        = 0 ; 
   int NHighMT3TeV        = 0 ; 

     for(int eventLoop=0; eventLoop<TotalN ; eventLoop++) {
        inChain->GetEntry(eventLoop);
        if((eventLoop%per99) == 0) cout << "Running " << (per100++ * 10)<< " % " << eventLoop << " / " << TotalN << endl;
        eleIdTCA->Clear();
        npknu::Evt* evtPtr = (npknu::Evt*)evtTCA->At(0);
        //if(evtPtr->run != 282037) continue;
        //if( 357 > evtPtr->lumi or 360 <evtPtr->lumi ) continue;
        //if((long(evtPtr->event) == 660799612 or long(evtPtr->event) == 662973378 or long(evtPtr->event) ==662847956)){ ///* (evtPtr->event != 663395404)*/) continue;
        //cout << "==> RunH " << evtPtr->run << " Event " << long(evtPtr->event) << " Lumi " << evtPtr->lumi << endl;
        double pileup       = (evtPtr->isRealData) ? -999 : npknu::GetPileupTrue(pileupTCA);
        double pileupWeight = (evtPtr->isRealData) ?  1.0 : puWeight->getWeight(pileup);
        double tWeight      = pileupWeight * lumiWeight; // Notice:We need additional weight for gen failed event (Ngen/NGencut)
        int nVtx = vertexTCA->GetEntries();

        /*************************************/
        /*Both Data and MC MET Filters passed*/
        /*************************************/
        npknu::MET* metPtr  = (npknu::MET*)metTCA->At(0);

        isMETAllFilter = ( metPtr->Flag_HBHENoiseFilter == 1 && 
                        metPtr->Flag_HBHENoiseIsoFilter == 1 && 
                        metPtr->Flag_goodVertices ==1 && 
                        /*metPtr->Flag_CSCTightHaloFilter ==1 && */
                        metPtr->Flag_EcalDeadCellTriggerPrimitiveFilter ==1 && 
                        metPtr->Flag_eeBadScFilter ==1 && 
                        metPtr->Flag_BadPFMuonFilter ==1 && 
                        metPtr->Flag_BadChargedCandidateFilter ==1 && 
                       (metPtr->Flag_globalSuperTightHalo2016Filter ==1  ||  metPtr->Flag_globalTightHalo2016Filter ==1) ) ;
       //isMETFilter = (metPtr->Flag_METFilters == 1) ;
       //if(!isMETFilter){ Nmetfil ++ ; cout << " METFilters unPASS" <<endl;}//cout BadMET
       //if(!isMETAllFilter){ Nmetallfil ++ ;//cout BadMET
       //    cout << " HBHENoise " << metPtr->Flag_HBHENoiseFilter <<
       //            " HBHENoiseIso " << metPtr->Flag_HBHENoiseIsoFilter<<
       //            " goodVertices " << metPtr->Flag_goodVertices <<
       //            " CSCTightHalo " << metPtr->Flag_CSCTightHaloFilter <<
       //            " EcalDeadCell " << metPtr->Flag_EcalDeadCellTriggerPrimitiveFilter <<
       //            " globalSuperTightHalo2016 "  << metPtr->Flag_globalSuperTightHalo2016Filter << " globalTightHalo2016 " << metPtr->Flag_globalTightHalo2016Filter << endl;
       //}

        if(!isMETAllFilter) continue;
        NpassMET ++ ; 

       bool Veto2ndheep = false;
       bool Veto2ndele = false;
       bool has2ndMuon = false;

        bool isElePreselec = npknu::GetElectronPreselection(electronTCA, 35.0, 0.0, 2.5); /* ptCut, etaCut*/
        if(!isElePreselec) continue;

        npknu::GetPassTriggerTCA(triggerTCA, trigNamePassTCA, trigNameVec);
        if(trigNamePassTCA->GetEntries() == 0) continue;
        Ntrigg++ ;

        npknu::GetElectronIdTCA(electronTCA, eleIdTCA, 4, 50.0, 2.5); /*HEEP, ptCut, etaCut*/
        if(eleIdTCA->GetEntries() ==0) continue;
        NHEEP++;
        npknu::Electron* heepPtr0 = (npknu::Electron*)eleIdTCA->At(0) ;
        if(eleIdTCA->GetEntries() > 1){
           NHEEP2++;
           npknu::Electron* heepPtr1 = (npknu::Electron*)eleIdTCA->At(1);
           if(heepPtr1->pt > 25.0){
		Veto2ndheep = true;
               if(heepPtr1->charge * heepPtr0->charge < 0){
                  hetc_pre[0]->Fill(heepPtr0->GetM(heepPtr1), tWeight    );
                  hetc_kin[0]->Fill(heepPtr0->GetM(heepPtr1), tWeight    );
                }
	  }
        }
	for(int i=0; i< electronTCA->GetEntries(); i++){
           npknu::Electron* elePtr = (npknu::Electron*)electronTCA->At(i);
           if(elePtr->isCutPassHEEP != 1 && elePtr->pt > 25.0 && (elePtr->isCutPassVeto == 1  or elePtr->isCutPassLoose == 1 or elePtr->isCutPassMedium == 1)){ 
              NHEEP1Veto25++;
              Veto2ndele = true;
           }
        }

       for(int i=0; i< muonTCA->GetEntries(); i++){
          npknu::Muon* muPtr = (npknu::Muon*)muonTCA->At(0);
          if( (muPtr->isLooseMuon == 1 or muPtr->isSoftMuonGV==1 or muPtr->isTightMuonGV ==1) && muPtr->pt > 25) {
             has2ndMuon = true;
	     Nhas2ndMuon++;
          }
       }

       if(Veto2ndheep==1 or Veto2ndele==1 or has2ndMuon==1 )continue; 
       NPreSelection++;

        /* Muon Plots */
       if(heepPtr0->isEB() == 1) hele_pre[0] ->Fill(heepPtr0->pt    ,tWeight  )  ;
       if(heepPtr0->isEE() == 1) hele_pre[1] ->Fill(heepPtr0->pt   ,tWeight  )  ;
        hele_pre[2] ->Fill(heepPtr0->pt   ,tWeight  )  ;
        hele_pre[3] ->Fill(heepPtr0->scEta() ,tWeight  )  ;
        hele_pre[4] ->Fill(heepPtr0->phi   ,tWeight  )  ;
        hele_pre[5] ->Fill(heepPtr0->full5x5_sigmaIetaIeta   ,tWeight  )  ;
        hele_pre[6] ->Fill(heepPtr0->hcalOverEcal   ,tWeight  )  ;
        hele_pre[7] ->Fill(heepPtr0->dEtaSeedAtVtx   ,tWeight  )  ;
        hele_pre[8] ->Fill(heepPtr0->deltaPhiSuperClusterTrackAtVtx   ,tWeight  )  ;
        /* MET Plots */
        hmet_pre[0] ->Fill(metPtr->pt     ,tWeight  ) ;
        hmet_pre[1] ->Fill(metPtr->phi    ,tWeight  ) ;
        /* kinematic Plots */
        hetc_pre[1] ->Fill(nVtx  ,lumiWeight ) ;
        hetc_pre[2] ->Fill(nVtx  ,tWeight )  ;
        hetc_pre[3] ->Fill(heepPtr0->GetMt(metPtr)  ,tWeight )  ;
        hetc_pre[4] ->Fill(heepPtr0->deltaPhi(metPtr)  ,tWeight )  ;
        hetc_pre[5] ->Fill(heepPtr0->etRatio(metPtr)  ,tWeight )  ;

       if(heepPtr0->etRatio(metPtr)<= 0.4 or heepPtr0->etRatio(metPtr)>= 1.5) continue ;
       if(heepPtr0->deltaPhi(metPtr->phi) <=2.5)  continue ;
       NKinCutPass++;
        /* Muon Plots */
       if(heepPtr0->isEB() == 1) hele_kin[0] ->Fill(heepPtr0->pt    ,tWeight  )  ;
       if(heepPtr0->isEE() == 1) hele_kin[1] ->Fill(heepPtr0->pt   ,tWeight  )  ;
        hele_kin[2] ->Fill(heepPtr0->pt   ,tWeight  )  ;
        hele_kin[3] ->Fill(heepPtr0->scEta() ,tWeight  )  ;
        hele_kin[4] ->Fill(heepPtr0->phi   ,tWeight  )  ;
        hele_kin[5] ->Fill(heepPtr0->full5x5_sigmaIetaIeta   ,tWeight  )  ;
        hele_kin[6] ->Fill(heepPtr0->hcalOverEcal   ,tWeight  )  ;
        hele_kin[7] ->Fill(heepPtr0->dEtaSeedAtVtx   ,tWeight  )  ;
        hele_kin[8] ->Fill(heepPtr0->deltaPhiSuperClusterTrackAtVtx   ,tWeight  )  ;
        /* MET Plots */
        hmet_kin[0] ->Fill(metPtr->pt     ,tWeight  ) ;
        hmet_kin[1] ->Fill(metPtr->phi    ,tWeight  ) ;
        /* kinematic Plots */
        hetc_kin[1] ->Fill(nVtx  ,lumiWeight ) ;
        hetc_kin[2] ->Fill(nVtx  ,tWeight )  ;
        hetc_kin[3] ->Fill(heepPtr0->GetMt(metPtr)  ,tWeight )  ;
        hetc_kin[4] ->Fill(heepPtr0->deltaPhi(metPtr)  ,tWeight )  ;
        hetc_kin[5] ->Fill(heepPtr0->etRatio(metPtr)  ,tWeight )  ;


       if(heepPtr0->GetMt(metPtr) >= 1000.0){
          NHighMT++ ; 
          if(heepPtr0->GetMt(metPtr) < 2000.0) { cout << "1TeV HighMT PASS ==> " << "MT " <<  heepPtr0->GetMt(metPtr) << " GeV "; NHighMT1TeV++ ;  }
          if(heepPtr0->GetMt(metPtr) >= 2000.0 and heepPtr0->GetMt(metPtr) < 3000.0){  cout << "2TeV HighMT PASS ==> " << "MT " <<  heepPtr0->GetMt(metPtr) << " GeV " ; NHighMT2TeV++ ;  }
          if(heepPtr0->GetMt(metPtr) >= 3000.0 and heepPtr0->GetMt(metPtr) < 4000.0){  cout << "3TeV HighMT PASS ==> " << "MT " <<  heepPtr0->GetMt(metPtr) << " GeV " ; NHighMT3TeV++ ;  }
          if(heepPtr0->GetMt(metPtr) >= 4000.0) { cout << "4TeV HighMT PASS ==> " << "MT " <<  heepPtr0->GetMt(metPtr) << " GeV ";}

          cout << "Event " << evtPtr->run << ":" << long(evtPtr->event) << ":" << evtPtr->lumi << " AllPass => MT " << heepPtr0->GetMt(metPtr) << " Et/MET " << heepPtr0->etRatio(metPtr) << " dPhi " << heepPtr0->deltaPhi(metPtr->phi) << " Info=> " ;
          heepPtr0->PrintObj(); metPtr->PrintMET();
      }
  }//event loop  

           hele_pre[0] ->Write() ;      hele_kin[0] ->Write() ;
           hele_pre[1] ->Write() ;      hele_kin[1] ->Write() ;
           hele_pre[2] ->Write() ;      hele_kin[2] ->Write() ;
           hele_pre[3] ->Write() ;      hele_kin[3] ->Write() ;
           hele_pre[4] ->Write() ;      hele_kin[4] ->Write() ;
           hele_pre[5] ->Write() ;      hele_kin[5] ->Write() ;
           hele_pre[6] ->Write() ;      hele_kin[6] ->Write() ;
           hele_pre[7] ->Write() ;      hele_kin[7] ->Write() ;
           hele_pre[8] ->Write() ;      hele_kin[8] ->Write() ;

          hmet_pre[0] ->Write() ;      hmet_kin[0] ->Write() ;
          hmet_pre[1] ->Write() ;      hmet_kin[1] ->Write() ;

          hetc_pre[0] ->Write() ;      hetc_kin[0] ->Write() ;
          hetc_pre[1] ->Write() ;      hetc_kin[1] ->Write() ;
          hetc_pre[2] ->Write() ;      hetc_kin[2] ->Write() ;
          hetc_pre[3] ->Write() ;      hetc_kin[3] ->Write() ;
          hetc_pre[4] ->Write() ;      hetc_kin[4] ->Write() ;
          hetc_pre[5] ->Write() ;      hetc_kin[5] ->Write() ;
          outFile->Write();


        cout << "## reMiniAOD Data " << fileName ;

        cout << "@Num of TotEvent  = " << setw(12) << TotalN             << " : " << setw(10) <<  double(TotalN )/double(TotalN )   <<  endl ;
        cout << " METAllFilter     = " << setw(12) << Nmetallfil            << " : " << setw(10) <<  double(Nmetallfil  )/double(TotalN ) << " Filtered in tot" << endl ;
        cout << " METFilter(BadMu) = " << setw(12) << Nmetfil           << " : " << setw(10) <<  double(Nmetfil )/double(TotalN ) << " Filtered in tot" << endl ;
        cout << "@Num of MET PASS  = " << setw(12) << NpassMET           << " : " << setw(10) <<  double(NpassMET)/double(TotalN ) << " Survived in tot" << endl ;
        cout << "@Num of HLT_Ele   = " << setw(12) << Ntrigg             << " : " << setw(10) <<  double(Ntrigg  )/double(NpassMET) << " Survived in NpassMET" << endl ;
        cout << " Num of HEEPexist = " << setw(12) << NHEEP       << " : " << setw(10) <<  double(NHEEP   )/double(Ntrigg) << " Survived in Ntrigg" << endl ;
        cout << " Num of HEEP 2    = " << setw(12) << NHEEP2       << " : " << setw(10) <<  double(NHEEP2   )/double(NHEEP) << " Survived" << endl ;
        cout << " Num of Veto25    = " << setw(12) << NHEEP1Veto25       << " : " << setw(10) <<  double(NHEEP1Veto25   )/double(NHEEP) << " Vetoed" << endl ;
        cout << " Num of 2ndMu Veto= " << setw(12) << Nhas2ndMuon       << " : " << setw(10) <<  double(Nhas2ndMuon)/double(NHEEP) << "Failed in HEEP " << endl;
        cout << "@Num of PreSelec  = " << setw(12) << NPreSelection     << " : " << setw(10) <<  double(NPreSelection)/double(Ntrigg) << " Survived in Ntrigg" << endl ;
        cout << "@Num of KinSelec  = " << setw(12) << NKinCutPass        << " : " << setw(10) <<  double(NKinCutPass)/double(NPreSelection) << " Survived in PreSelection" << endl ;
        cout << "@Num of HighMT >=1TeV = " << setw(12) << NHighMT            << " == " << NHighMT1TeV + NHighMT2TeV + NHighMT3TeV <<endl;
        cout << "# HighMT 1=< MT <2TeV = " << setw(12) << NHighMT1TeV        << endl;
        cout << "# HighMT 2=< MT <3TeV = " << setw(12) << NHighMT2TeV        << endl;
        cout << "# HighMT 3=< MT <4TeV = " << setw(12) << NHighMT3TeV        << endl;

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
        delete eleIdTCA         ;
        cout << "*** File Closed " << endl;

        return 0;
}
