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
//#include "./ProdNpKNU/src/NpKNU.hh"
#include "/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ProdNpKNU/src/NpKNU.hh"
   using namespace std;
   using namespace TMath;
   npknu::PUReweight* puWeight ;
//   const char *  mustr_ori[11] = {"mu_Pt_ori", "mu_Eta_ori","mu_Phi_ori", "mu_momQuality_ori","mu_TrkIsoR03_sumPt_ori" "mu_relTrkIsoR03_ori", "mu_nPixelHits_ori", "mu_nTrkLayers_ori", "mu_nMatchedStations_ori","mu_nValidMuonHits_ori", "mu_dz_ori", "mu_dxy_ori"};
//   const char *  mustr_pre[11] = {"mu_Pt_pre", "mu_Eta_pre","mu_Phi_pre", "mu_momQuality_pre","mu_TrkIsoR03_sumPt_pre" "mu_relTrkIsoR03_pre", "mu_nPixelHits_pre", "mu_nTrkLayers_pre", "mu_nMatchedStations_pre","mu_nValidMuonHits_pre", "mu_dz_pre", "mu_dxy_pre"};
//   const char *  mustr_kin[11] = {"mu_Pt_kin", "mu_Eta_kin","mu_Phi_kin", "mu_momQuality_kin","mu_TrkIsoR03_sumPt_kin" "mu_relTrkIsoR03_kin", "mu_nPixelHits_kin", "mu_nTrkLayers_kin", "mu_nMatchedStations_kin","mu_nValidMuonHits_kin", "mu_dz_kin", "mu_dxy_kin"};
//
//   Int_t mubin[11]      = { 300,  100, 100, 200, 200,  15, 30, 30, 100,100 ,60 };
//   Double_t muxmin[11]  = { 0.0, -2.5,-3.2, 0.0, 0.0, 0.0,0.0,0.0,0.0 , -0.6,-0.3 };
//   Double_t muxmax[11]  = {3000,  2.5, 3.2, 0.5, 0.2,  15, 30, 30, 100, 0.6, 0.3};
//   TString muxTitle[11] = {"Muon P_{T} [GeV]", "Muon #eta", "Muon #phi", "#sigma_{Pt}/P_{T}(<0.3)","relTrkIsoR03(<0.1)", "nPixelHits(>0)","nTrackerLayers(>5)","nMatchedStations(>1)","nValidMuonHits(>0)","dZ(<5mm)","dXY(<2mm)"};
//   
//   const char *  etcstr_ori[8] = {"nPV_ori","nPVw_ori","wp_MtT1_ori","wp_MtT1XY_ori","wp_DphiT1_ori","wp_DphiT1XY_ori","wp_EtRatioT1_ori","invM_ori"}  ;   
//   const char *  etcstr_pre[8] = {"nPV_pre","nPVw_pre","wp_MtT1_pre","wp_MtT1XY_pre","wp_DphiT1_pre","wp_DphiT1XY_pre","wp_EtRatioT1_pre","invM_pre"}  ;   
//   const char *  etcstr_kin[8] = {"nPV_kin","nPVw_kin","wp_MtT1_kin","wp_MtT1XY_kin","wp_DphiT1_kin","wp_DphiT1XY_kin","wp_EtRatioT1_kin","invM_kin"}  ;   
//
//
//   Int_t etcbin[8]      ={ 100,  100, 700 , 700,  200, 200, 200, 100 }  ;   
//   Double_t etcxmin[8]  ={ 0.0, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 60.0 }  ;   
//   Double_t etcxmax[8]  ={ 100,  100, 7000, 7000, 5.0, 5.0,  10, 1500 }  ;   
//   TString etcxTitle[8] ={"nPV(w/o PUreWeight)","nPV(w PUreweight)","M_{T} [GeV]","M_{T}(T1XY) [GeV]","#Delta #phi","#Delta #phi(T1XY)","E_{T}^{#mu}/E_{T}^{miss}", "inv M_{#mu#mu}" }  ;  
//   
//   TH1F* hgenMass_ori ;  
//   TH1F* hgenMass_pre ;  
//   TH1F* hgenMass_kin ;  
//   
//   TH1F* hmu_ori[11] ;  
//   TH1F* hmu_pre[11] ;  
//   TH1F* hmu_kin[11] ;  
//   
//   TH1F* hmet_ori[8] ;  
//   TH1F* hmet_pre[8] ;  
//   TH1F* hmet_kin[8] ; 
//    
//   TH1F* hetc_ori[8] ;  
//   TH1F* hetc_pre[8] ;  
//   TH1F* hetc_kin[8] ;  
//
//   int NTotal   = 0;
//   int NBadPFmu = 0;
//   int NBadch   = 0;
////   int Nmupt50  = 0;
//   int Ntrigg   = 0;
//   int Nveto2nd = 0;
//   int NhighPtID= 0;
//   int NmuIso   = 0;
//   int Nmupt53  = 0;
//   int Nkinet   = 0 ;
//   int Nkinphi  = 0;
//
////   int NGencut  = 0;
////   int NGenMasscut  = 0;
//
//   TH2F* h2_hlt ;
//   TH2F* h2_id  ;
//   TH2F* h2_iso ;

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
class getSFHist {
public:
   TH2F* sfHist;
   double HLTWeight;
   double IDWeight;
   double ISOWeight;
   double value;
   double error;

   getSFHist(TString sfFileName, TString H2Name){
        TFile* sfFile = TFile::Open(sfFileName, "READ");
        sfHist = (TH2F*)sfFile->Get(H2Name);
        sfHist->SetDirectory(0);
        sfFile->Close();
   };
  ~getSFHist() {};

   double getSFweight(double pt, double abseta){
        value =0;
        error=0;
        int nBins = 0;
        nBins = sfHist->GetNbinsY();
        double ptmax = -1;
        if( ptmax <= 0.){
            ptmax = sfHist->GetYaxis()->GetBinCenter(nBins);
        }
        double value = sfHist->GetBinContent(sfHist->FindBin(abseta, TMath::Min(pt, ptmax)));
        double error = sfHist->GetBinError  (sfHist->FindBin(abseta, TMath::Min(pt, ptmax)));
        /*cout << "pt " << pt << ", ptmax " << ptmax << ", AEta " << abseta << ", sf= "<< value << " +/- " << error << endl; */
        return value ;
           }
};
        
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
     TClonesArray* muIdTCA        = new TClonesArray("npknu::Muon")           ;

     if(type == 1) {
        lumiWeight = 1.0;
     } else {
        lumiWeight = lumi *1000 * xsec / gen;
        puWeight = new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_RunG-H.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
     }

//pu_RunB-F.root
//pu_RunG-H.root 

     cout << "Type " << type << ": lumiWeight " << lumiWeight << "=" << lumi << "*" << xsec << "/" << gen << "  : " << fileName << endl;
     //for(int j=0; j<11; j++){ //Draw Variables
     //   hmu_ori[j] = new TH1F(mustr_ori[j],mustr_ori[j], mubin[j], muxmin[j], muxmax[j]); hmu_ori[j]->GetXaxis()->SetTitle(muxTitle[j]);  hmu_ori[j]->Sumw2();
     //   hmu_pre[j] = new TH1F(mustr_pre[j],mustr_pre[j], mubin[j], muxmin[j], muxmax[j]); hmu_pre[j]->GetXaxis()->SetTitle(muxTitle[j]);  hmu_pre[j]->Sumw2();
     //   hmu_kin[j] = new TH1F(mustr_kin[j],mustr_kin[j], mubin[j], muxmin[j], muxmax[j]); hmu_kin[j]->GetXaxis()->SetTitle(muxTitle[j]);  hmu_kin[j]->Sumw2();
     //}
     //for(int k=0; k<8; k++){ //Draw Variables
     //   hmet_ori[k] = new TH1F(metstr_ori[k],metstr_ori[k], metbin[k], metxmin[k], metxmax[k]); hmet_ori[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_ori[k]->Sumw2();
     //   hmet_pre[k] = new TH1F(metstr_pre[k],metstr_pre[k], metbin[k], metxmin[k], metxmax[k]); hmet_pre[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_pre[k]->Sumw2();
     //   hmet_kin[k] = new TH1F(metstr_kin[k],metstr_kin[k], metbin[k], metxmin[k], metxmax[k]); hmet_kin[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_kin[k]->Sumw2();
     //}
     //for(int p=0; p<8; p++){ //Draw Variables
     //   hetc_ori[p] = new TH1F(etcstr_ori[p],etcstr_ori[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_ori[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_ori[p]->Sumw2();
     //   hetc_pre[p] = new TH1F(etcstr_pre[p],etcstr_pre[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_pre[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_pre[p]->Sumw2();
     //   hetc_kin[p] = new TH1F(etcstr_kin[p],etcstr_kin[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_kin[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_kin[p]->Sumw2();
     //}

     //TFile* outFile = new TFile(argv[6],"recreate");
     int TotalN = (int)inChain->GetEntries();
     int per99 = (TotalN > 10) ? TotalN / 9 : 1;
     int per100 = 0;

     std::vector<std::string> trigNameVec;
     trigNameVec.push_back("HLT_Mu50_v*");
     trigNameVec.push_back("HLT_TkMu50_v*");
     TClonesArray* trigNamePassTCA = new TClonesArray("npknu::Trigger");

//   int NVetoManyMu = 0;
////   int NVeto2nd = 0;
//   int NhighPtID1= 0;
//   int NhighPtID2= 0;
//   int NmuIso   = 0;
//   int NmuonTkcut =0;
//   int NBadmucut =0;
//   int Nmt60  = 0;
//   int Nmt80  = 0;
//   int Nkinet   = 0 ;
//   int Nkindphi  = 0;

   bool isMETFilter = false;
   bool isMETBadPFMuon = true;
   bool isMETBadChMuon = true;
   bool isMuPreselection = false;

   int Nmupt100  = 0;
   int NBadGBMuon = 0;
   int NBadPFMu = 0;
   int NnoBadch   = 0;
   int Nmetfil  = 0;
   int NpassMET =0;
   int Ntrigg   = 0;
   int NHighPt =0;
   int NHighPtFailTrk =0;
   int NHighPt2 =0;
   int NHighPt2FailTrk =0;
   int NnoHighPtButBig =0;
   int NnoHighPtnoVeto25 =0;
   int NClosedR =0;
   int NClosedR2 =0;
   int NGoodMu =0;
   int NBadMu =0;
   int NGoodMuTest =0;
   int NBadMuTest =0;
   int NBadFirstHighPt =0;

     for(int eventLoop=0; eventLoop<TotalN ; eventLoop++) {
        inChain->GetEntry(eventLoop);
        if((eventLoop%per99) == 0) cout << "Running " << (per100++ * 10)<< " % " << eventLoop << " / " << TotalN << endl;
        muIdTCA->Clear();
        npknu::Evt* evtPtr = (npknu::Evt*)evtTCA->At(0);

        /********************************************************/
        /*Both Data and MC MET Filters passed & MET > 100GeV cut*/
        /********************************************************/
        npknu::MET* metPtr  = (npknu::MET*)metTCA->At(0);
        isMETFilter = ( metPtr->Flag_HBHENoiseFilter == 1 &&  metPtr->Flag_HBHEIsoNoiseFilter == 1 && metPtr->Flag_goodVertices ==1 && metPtr->Flag_CSCTightHaloFilter ==1 && metPtr->Flag_EcalDeadCellTriggerPrimitiveFilter ==1 && metPtr->Flag_eeBadScFilter ==1 && (metPtr->Flag_globalSuperTightHalo2016Filter ==1 ||  metPtr->Flag_globalTightHalo2016Filter ==1) ) ;
        if (isMETFilter) Nmetfil ++ ; ////////////////////////////////

        isMETBadPFMuon = ( metPtr->filterbadPFMuon == 1 || metPtr->IsNotBadPFMuon == 1) ;
        if (!isMETBadPFMuon) NBadPFMu ++; ////////////////////////

        isMETBadChMuon = ( metPtr->filterbadChCandidate == 1 ||  metPtr->IsNotBadChargedCandidate == 1 ) ;
        if(!isMETBadChMuon ) NnoBadch ++; ///////////////////////
        if(metPtr->IsBadGlobalMuonTagger == 1){NBadGBMuon++;}//Count BadGBMu
        if(isMETFilter || !isMETBadPFMuon || !isMETBadChMuon  ) continue;
        NpassMET ++ ; ////////////////////////////////

        bool isMuPreselection = npknu::GetMuonPreselection(muonTCA, 100.0, 0.0, 2.4); /*highPtMuonGVId, ptCut, etaCut*/
        if(isMuPreselection) Nmupt100++; //////////////////////////
        if(!isMuPreselection) continue;
        ///*******************************************************/
        /*SingleMuon HLT Filter to Data and MC(17) (Mu50 OR TkMu50)     */
        /*******************************************************/

       double highPtmuonPt0 = 0.0;
       double highPtmuonPt1 = 0.0;
       bool highPtTrkPass = false;
       bool highPtTrkPass2 = false;
       bool VetoPass = false;

        npknu::GetPassTriggerTCA(triggerTCA, trigNamePassTCA, trigNameVec);
        if(trigNamePassTCA->GetEntries() == 0) continue;
        Ntrigg++ ;
        double pileup       = (evtPtr->isRealData) ? -999 : npknu::GetPileupTrue(pileupTCA);
        double pileupWeight = (evtPtr->isRealData) ?  1.0 : puWeight->getWeight(pileup);
        double tWeight      = pileupWeight * lumiWeight; // Notice:We need additional weight for gen failed event (Ngen/NGencut)
        int nVtx = vertexTCA->GetEntries();
        double dR = 0.0 ;
        npknu::GetMuonIdTCA(muonTCA, muIdTCA, 6, 53.0, 0.0, 2.4); /*highPtMuonGVId, ptCut, etaCut*/
        if(muIdTCA->GetEntries() ==0) continue;
        npknu::Muon* highmuPtr0 = (npknu::Muon*)muIdTCA->At(0);
        highPtmuonPt0 = highmuPtr0->pt;
        NHighPt++;
        cout << "Num " << muIdTCA->GetEntries() << endl;
        cout << "1st highPtId Muon " << highmuPtr0->pt_TPq << " gb " << highmuPtr0->pt_GBq << " bt " << highmuPtr0->pt_BTq << " in " << highmuPtr0->pt_INq <<  ", eta " << highmuPtr0->eta << " phi " << highmuPtr0->phi << " trkSumPt " << highmuPtr0->TrkIsoR03_sumPt <<" relTrkIso " << highmuPtr0->relTrkIsoR03_BT <<  endl;
        if(highmuPtr0->relTrkIsoR03_BT >0.1) {highPtTrkPass = false; cout << "Not Passed TrkIso Cut " << endl; NHighPtFailTrk++;}
        if(highmuPtr0->relTrkIsoR03_BT <0.1 ) {highPtTrkPass = true; }
        if(muIdTCA->GetEntries() > 1){  
           NHighPt2++;
           npknu::Muon* highmuPtr1 = (npknu::Muon*)muIdTCA->At(1);
           highPtmuonPt1 = highmuPtr1->pt;
        cout << " 2 highPtId Muon " << highmuPtr1->pt_TPq << " gb " << highmuPtr1->pt_GBq << " bt " << highmuPtr1->pt_BTq << " in " << highmuPtr1->pt_INq <<  ", eta " << highmuPtr1->eta << " phi " << highmuPtr1->phi << " trkSumPt " << highmuPtr1->TrkIsoR03_sumPt << " relTrkIso " << highmuPtr1->relTrkIsoR03_BT <<  " dR " << highmuPtr0->deltaR(highmuPtr1) << endl;
         if(highmuPtr1->relTrkIsoR03_BT >0.1) {cout << "Not Passed TrkIso Cut " << endl; NHighPt2FailTrk++;}
         if(highPtTrkPass == true and highmuPtr1->pt < 25.0 ) {cout << "Good Event Veto 2ndMu >25 Passed " << endl; VetoPass = true; }
         if(highPtTrkPass == false and highmuPtr1->relTrkIsoR03_BT <0.1 and  highmuPtr1->pt > 53.0 ) {cout << "Bad Reject First HighPt muon and Second need Check" << endl; highPtTrkPass2 = true; NBadFirstHighPt++;}
        }

        if((muIdTCA->GetEntries() == 1 and highPtTrkPass == 1) or (muIdTCA->GetEntries() > 1 and VetoPass == 1)){ 
            NGoodMu++;
        }else{
            NBadMu++;
        }

        bool foundbiggermuon = false ;
	for(int i=0; i< muonTCA->GetEntries(); i++){
           npknu::Muon* muPtr = (npknu::Muon*)muonTCA->At(i);
           if(muPtr->isHighPtMuonGV == 1) continue;
           if(muPtr->pt > highPtmuonPt0) { cout << "Bigger than HighptMuon :" << i << " Muon pt " << muPtr->pt_TPq << " gb " << muPtr->pt_GBq << " bt " << muPtr->pt_BTq << " in " << muPtr->pt_INq <<  ", eta " << muPtr->eta << " phi " << muPtr->phi << " trkSumPt " << muPtr->TrkIsoR03_sumPt <<" relTrkIso " << muPtr->relTrkIsoR03_BT << " ,ID Soft " << muPtr->isSoftMuonGV << ", Tight " <<  muPtr->isTightMuonGV << ", HighPt " << muPtr->isHighPtMuonGV <<  endl; NnoHighPtButBig++;
              if(highmuPtr0->deltaR(muPtr) < 0.3){ cout << "Find Close dR = " << highmuPtr0->deltaR(muPtr) << endl; NClosedR++;}
            }
           if(muPtr->pt > 25) foundbiggermuon = true;
 
	  if(muPtr->pt < 25.0) {cout << "Smaller than 25GeV :" << i <<  " Muon pt " << muPtr->pt_TPq << " gb " << muPtr->pt_GBq << " bt " << muPtr->pt_BTq << " in " << muPtr->pt_INq <<  ", eta " << muPtr->eta << " phi " << muPtr->phi << " trkSumPt " << muPtr->TrkIsoR03_sumPt <<" relTrkIso " << muPtr->relTrkIsoR03_BT << " ,ID Soft " << muPtr->isSoftMuonGV << ", Tight " <<  muPtr->isTightMuonGV << ", HighPt " << muPtr->isHighPtMuonGV <<  endl; NnoHighPtnoVeto25++;         
              if(highmuPtr0->deltaR(muPtr) < 0.3){ cout << "Find Close dR = " << highmuPtr0->deltaR(muPtr) << endl; NClosedR2++;}
          }
        }       
       
        if(((muIdTCA->GetEntries() == 1 and highPtTrkPass == 1) or (muIdTCA->GetEntries() > 1 and VetoPass == 1)) and (foundbiggermuon == 0)){ 
            NGoodMuTest++;
        }else{
            NBadMuTest++;
        }
        cout << "--------MET " << metPtr->pt << " vertex " << nVtx << " puWeight " << pileupWeight << endl;
}


   cout << "Num of TotEvent    = " << TotalN     << endl ;
   cout << "Mun of NBadGBMuon  = " << NBadGBMuon << endl ;
   cout << "Num of NBadPFMu  = " << NBadPFMu << endl ;
   cout << "Num of NnoBadch    = " << NnoBadch   << endl ;
   cout << "Num of Nmetfil     = " << Nmetfil    << endl ;
   cout << "Num of NpassMET    = " << NpassMET   << endl ;
   cout << "Num of Nmupt100    = " << Nmupt100   << endl ;
   cout << "Num of Ntrigg      = " << Ntrigg     << endl ;
   cout << "Num of NHighPt1    = " << NHighPt         << endl ;
   cout << "Num of 1FailTrk    = " << NHighPtFailTrk         << endl ;
   cout << "Num of NHighPt2    = " << NHighPt2         << endl ;
   cout << "Num of 2FailTrk    = " << NHighPt2FailTrk         << endl ;
   cout << "Num of noHighPt Big= " << NnoHighPtButBig         << endl ;
   cout << "Num of NClosedR    = " << NClosedR     << endl ;
   cout << "Num of noHighPt< 25= " << NnoHighPtnoVeto25         << endl ;
   cout << "Num of NClosedR2   = " << NClosedR2     << endl ;
   cout << "Num of Bad1stHighPt= " << NBadFirstHighPt     << endl ;
   cout << "Num of NGoodMu     = " << NGoodMu     << endl ;
   cout << "Num of NBadMu     = " << NBadMu     << endl ;
   cout << "Num of NGoodMuTest     = " << NGoodMuTest     << endl ;
   cout << "Num of NBadMuTest     = " << NBadMuTest     << endl ;

        /*****************************************************************************/
        /*Both Data and MC HighPtMuID & Mu Pt > 55 GeV& |eta|<2.4 &Loose TrkIso <0.1 */
        /*****************************************************************************/
//           /* Muon Plots */
//           hmu_ori[0] ->Fill(muPtr->pt    ,lumiWeight  )  ;  
//           hmu_ori[1] ->Fill(muPtr->eta   ,lumiWeight  )  ;  
//           hmu_ori[2] ->Fill(muPtr->phi   ,lumiWeight  )  ;  
//           hmu_ori[3] ->Fill( ((muPtr->ptErr_BT)/(muPtr->pt_BT)) ,lumiWeight  )  ;  
//           hmu_ori[4] ->Fill(muPtr->relTrkIsoR03_BT    ,lumiWeight  )  ;  
//           hmu_ori[5] ->Fill(muPtr->nValidPixelHits    ,lumiWeight  )  ;  
//           hmu_ori[6] ->Fill(muPtr->nTrackerLayers     ,lumiWeight  )  ;  
//           hmu_ori[7] ->Fill(muPtr->nMatchedStations   ,lumiWeight  )  ;  
//           hmu_ori[8] ->Fill(muPtr->nValidMuonHits     ,lumiWeight  )  ; 
//           hmu_ori[9] ->Fill(muPtr->dz_GV_BT           ,lumiWeight  )  ;
//           hmu_ori[10]->Fill(muPtr->dxy_GV_BT          ,lumiWeight  )  ;
//           /* MET Plots */
//           hmet_ori[0] ->Fill(metPtr->pt     ,lumiWeight  ) ;  
//           hmet_ori[1] ->Fill(metPtr->phi    ,lumiWeight  ) ;  
//           hmet_ori[2] ->Fill(metPtr->sumEt  ,lumiWeight  ) ;  
//           hmet_ori[3] ->Fill(metPtr->metSignificance    ,lumiWeight  ) ;  
//           hmet_ori[4] ->Fill(metPtr->shiftedPt_Type1XY  ,lumiWeight  ) ;  
//           hmet_ori[5] ->Fill(metPtr->shiftedPhi_Type1XY ,lumiWeight  ) ;  
//           hmet_ori[6] ->Fill(metPtr->METPuppi_pt  ,lumiWeight  ) ;  
//           hmet_ori[7] ->Fill(metPtr->METPuppi_phi ,lumiWeight  ) ;  
//           /* kinematic Plots */
//           hetc_ori[0] ->Fill(nVtx  ,lumiWeight ) ;  
//           hetc_ori[1] ->Fill(nVtx  ,tWeight )  ;  
//           hetc_ori[2] ->Fill(muPtr->GetMt(metPtr)  ,lumiWeight )  ;  
//           hetc_ori[3] ->Fill(muPtr->GetMt(metPtr->shiftedPt_Type1XY, metPtr->shiftedPhi_Type1XY) ,lumiWeight );
//           hetc_ori[4] ->Fill(muPtr->deltaPhi(metPtr)  ,lumiWeight )  ;  
//           hetc_ori[5] ->Fill(muPtr->deltaPhi(metPtr->shiftedPhi_Type1XY)  ,lumiWeight )  ;  
//           hetc_ori[6] ->Fill(muPtr->etRatio(metPtr)  ,lumiWeight )  ;  
//
//          if ((muPtr->relTrkIsoR03_BT) >  0.10 ) continue; /*MuIsoLoose(dR<0.3) TrkRel.Iso<0.1*/
//           NmuIso++ ;
//
//          if((muPtr->hastunePMuonBestTrack == 0) or (muPtr->hasglobalTrack == 0) or (muPtr->hasinnerTrack == 0) ) continue;
//           NmuonTkcut++;
//          if( (muPtr->isBadPFMuon == 1) or (muPtr->isBadChHadron == 1 )) NBadmucut++;
//
//        /*****************************************************************************/
//        /*New ScaleFactor HLT(only G,H availble so far ,ID & ISO BCDEF era1, GH era2 */
//        /*****************************************************************************/
//           double MuHLTWeight = 1; double MuIDWeight = 1; double MuISOWeight = 1;
//           double muonWeight = MuHLTWeight * MuIDWeight * MuISOWeight ;
//           /* Muon Plots */
//           hmu_pre[0] ->Fill(muPtr->pt    ,tWeight * muonWeight  )  ;  
//           hmu_pre[1] ->Fill(muPtr->eta   ,tWeight * muonWeight  )  ;  
//           hmu_pre[2] ->Fill(muPtr->phi   ,tWeight * muonWeight  )  ;  
//           hmu_pre[3] ->Fill( ((muPtr->ptErr_BT)/(muPtr->pt_BT)) ,tWeight * muonWeight  )  ;  
//           hmu_pre[4] ->Fill(muPtr->relTrkIsoR03_BT    ,tWeight * muonWeight  )  ;  
//           hmu_pre[5] ->Fill(muPtr->nValidPixelHits    ,tWeight * muonWeight  )  ;  
//           hmu_pre[6] ->Fill(muPtr->nTrackerLayers     ,tWeight * muonWeight  )  ;  
//           hmu_pre[7] ->Fill(muPtr->nMatchedStations   ,tWeight * muonWeight  )  ;  
//           hmu_pre[8] ->Fill(muPtr->nValidMuonHits     ,tWeight * muonWeight  )  ;  
//           hmu_pre[9] ->Fill(muPtr->dz_GV_BT           ,tWeight * muonWeight )  ;
//           hmu_pre[10]->Fill(muPtr->dxy_GV_BT          ,tWeight * muonWeight )  ;
//           /* MET Plots */
//           hmet_pre[0] ->Fill(metPtr->pt     ,tWeight * muonWeight  ) ;  
//           hmet_pre[1] ->Fill(metPtr->phi    ,tWeight * muonWeight  ) ;  
//           hmet_pre[2] ->Fill(metPtr->sumEt  ,tWeight * muonWeight  ) ;  
//           hmet_pre[3] ->Fill(metPtr->metSignificance    ,tWeight * muonWeight  ) ;  
//           hmet_pre[4] ->Fill(metPtr->shiftedPt_Type1XY  ,tWeight * muonWeight  ) ;  
//           hmet_pre[5] ->Fill(metPtr->shiftedPhi_Type1XY ,tWeight * muonWeight  ) ;  
//           hmet_pre[6] ->Fill(metPtr->METPuppi_pt  ,tWeight * muonWeight  ) ;  
//           hmet_pre[7] ->Fill(metPtr->METPuppi_phi ,tWeight * muonWeight  ) ;  
//           /* kinematic Plots */
//           hetc_pre[0] ->Fill(nVtx  ,lumiWeight * muonWeight ) ;  
//           hetc_pre[1] ->Fill(nVtx  ,tWeight * muonWeight )  ;  
//           hetc_pre[2] ->Fill(muPtr->GetMt(metPtr)  ,tWeight * muonWeight )  ;  
//           hetc_pre[3] ->Fill(muPtr->GetMt(metPtr->shiftedPt_Type1XY, metPtr->shiftedPhi_Type1XY) ,tWeight * muonWeight );
//           hetc_pre[4] ->Fill(muPtr->deltaPhi(metPtr)  ,tWeight * muonWeight )  ;  
//           hetc_pre[5] ->Fill(muPtr->deltaPhi(metPtr->shiftedPhi_Type1XY)  ,tWeight * muonWeight )  ;  
//           hetc_pre[6] ->Fill(muPtr->etRatio(metPtr)  ,tWeight * muonWeight )  ;  
//
//        /************************************************************************************/
//        /* Wp Kinematic Cut:(Back to back)dPhi < 2.5 and (Energy Balance)0.4 < Et/MET < 1.5 */
//        /************************************************************************************/
//           if(muPtr->GetMt(metPtr) < 60 ) continue;
//           Nmt60++;
//           if((muPtr->etRatio(metPtr)< 0.4 or muPtr->etRatio(metPtr)> 1.5) )  continue ;
//           Nkinet++;
//           if((muPtr->deltaPhi(metPtr->phi) <2.5) or (muPtr->deltaPhi(metPtr->shiftedPhi_Type1XY) <2.5))  continue ;
//            Nkindphi++;
//           if((muPtr->GetMt(metPtr) < 80)) continue;
//           Nmt80++;
//           /* Muon Plots */
//           hmu_kin[0] ->Fill(muPtr->pt    ,tWeight * muonWeight  )  ;  
//           hmu_kin[1] ->Fill(muPtr->eta   ,tWeight * muonWeight  )  ;  
//           hmu_kin[2] ->Fill(muPtr->phi   ,tWeight * muonWeight  )  ;  
//           hmu_kin[3] ->Fill( ((muPtr->ptErr_BT)/(muPtr->pt_BT)) ,tWeight * muonWeight  )  ;  
//           hmu_kin[4] ->Fill(muPtr->relTrkIsoR03_BT    ,tWeight * muonWeight  )  ;  
//           hmu_kin[5] ->Fill(muPtr->nValidPixelHits    ,tWeight * muonWeight  )  ;  
//           hmu_kin[6] ->Fill(muPtr->nTrackerLayers     ,tWeight * muonWeight  )  ;  
//           hmu_kin[7] ->Fill(muPtr->nMatchedStations   ,tWeight * muonWeight  )  ;  
//           hmu_kin[8] ->Fill(muPtr->nValidMuonHits     ,tWeight * muonWeight  )  ;  
//           hmu_kin[9] ->Fill(muPtr->dz_GV_BT           ,tWeight * muonWeight )  ;
//           hmu_kin[10]->Fill(muPtr->dxy_GV_BT          ,tWeight * muonWeight )  ;
//           /* MET Plots */
//           hmet_kin[0] ->Fill(metPtr->pt     ,tWeight * muonWeight  ) ;  
//           hmet_kin[1] ->Fill(metPtr->phi    ,tWeight * muonWeight  ) ;  
//           hmet_kin[2] ->Fill(metPtr->sumEt  ,tWeight * muonWeight  ) ;  
//           hmet_kin[3] ->Fill(metPtr->metSignificance    ,tWeight * muonWeight  ) ;  
//           hmet_kin[4] ->Fill(metPtr->shiftedPt_Type1XY  ,tWeight * muonWeight  ) ;  
//           hmet_kin[5] ->Fill(metPtr->shiftedPhi_Type1XY ,tWeight * muonWeight  ) ;  
//           hmet_kin[6] ->Fill(metPtr->METPuppi_pt  ,tWeight * muonWeight  ) ;  
//           hmet_kin[7] ->Fill(metPtr->METPuppi_phi ,tWeight * muonWeight  ) ;  
//           /* kinematic Plots */
//           hetc_kin[0] ->Fill(nVtx  ,lumiWeight * muonWeight ) ;  
//           hetc_kin[1] ->Fill(nVtx  ,tWeight * muonWeight )  ;  
//           hetc_kin[2] ->Fill(muPtr->GetMt(metPtr)  ,tWeight * muonWeight )  ;  
//           hetc_kin[3] ->Fill(muPtr->GetMt(metPtr->shiftedPt_Type1XY, metPtr->shiftedPhi_Type1XY) ,tWeight * muonWeight );
//           hetc_kin[4] ->Fill(muPtr->deltaPhi(metPtr)  ,tWeight * muonWeight )  ;  
//           hetc_kin[5] ->Fill(muPtr->deltaPhi(metPtr->shiftedPhi_Type1XY)  ,tWeight * muonWeight )  ;  
//           hetc_kin[6] ->Fill(muPtr->etRatio(metPtr)  ,tWeight * muonWeight )  ;  
//
//        }
//        if(muIdTCA->GetEntries() > 1 ){
//           NhighPtID2++;
//           npknu::Muon* muPtr0 = (npknu::Muon*)muIdTCA->At(0);
//           npknu::Muon* muPtr1 = (npknu::Muon*)muIdTCA->At(1);
//
//           double MuHLTWeight0 = 1; double MuHLTWeight1 = 1;
//           double MuIDWeight0 = 1 ; double MuIDWeight1 = 1;
//           double MuISOWeight0 = 1; double MuISOWeight1 = 1;
//
//           double muonWeight2 =  MuHLTWeight0 * MuHLTWeight1 * MuIDWeight0 * MuIDWeight1 * MuISOWeight0 * MuISOWeight1 ;
//           hetc_ori[7]->Fill(muPtr0->GetM(muPtr1), lumiWeight) ;
//           hetc_pre[7]->Fill(muPtr0->GetM(muPtr1),   tWeight) ;
//           hetc_kin[7]->Fill(muPtr0->GetM(muPtr1),    tWeight * muonWeight2) ;
//
//        }
//    }
//          hmu_ori[0] ->Write() ;      hmu_pre[0] ->Write() ;      hmu_kin[0] ->Write() ;
//          hmu_ori[1] ->Write() ;      hmu_pre[1] ->Write() ;      hmu_kin[1] ->Write() ;
//          hmu_ori[2] ->Write() ;      hmu_pre[2] ->Write() ;      hmu_kin[2] ->Write() ;
//          hmu_ori[3] ->Write() ;      hmu_pre[3] ->Write() ;      hmu_kin[3] ->Write() ;
//          hmu_ori[4] ->Write() ;      hmu_pre[4] ->Write() ;      hmu_kin[4] ->Write() ;
//          hmu_ori[5] ->Write() ;      hmu_pre[5] ->Write() ;      hmu_kin[5] ->Write() ;
//          hmu_ori[6] ->Write() ;      hmu_pre[6] ->Write() ;      hmu_kin[6] ->Write() ;
//          hmu_ori[7] ->Write() ;      hmu_pre[7] ->Write() ;      hmu_kin[7] ->Write() ;
//          hmu_ori[8] ->Write() ;      hmu_pre[8] ->Write() ;      hmu_kin[8] ->Write() ;
//          hmu_ori[9] ->Write() ;      hmu_pre[9] ->Write() ;      hmu_kin[9] ->Write() ;
//          hmu_ori[10] ->Write() ;      hmu_pre[10] ->Write() ;      hmu_kin[10] ->Write() ;
//
//          hmet_ori[0] ->Write() ;      hmet_pre[0] ->Write() ;      hmet_kin[0] ->Write() ;
//          hmet_ori[1] ->Write() ;      hmet_pre[1] ->Write() ;      hmet_kin[1] ->Write() ;
//          hmet_ori[2] ->Write() ;      hmet_pre[2] ->Write() ;      hmet_kin[2] ->Write() ;
//          hmet_ori[3] ->Write() ;      hmet_pre[3] ->Write() ;      hmet_kin[3] ->Write() ;
//          hmet_ori[4] ->Write() ;      hmet_pre[4] ->Write() ;      hmet_kin[4] ->Write() ;
//          hmet_ori[5] ->Write() ;      hmet_pre[5] ->Write() ;      hmet_kin[5] ->Write() ;
//          hmet_ori[6] ->Write() ;      hmet_pre[6] ->Write() ;      hmet_kin[6] ->Write() ;
//          hmet_ori[7] ->Write() ;      hmet_pre[7] ->Write() ;      hmet_kin[7] ->Write() ;
//
//          hetc_ori[0] ->Write() ;      hetc_pre[0] ->Write() ;      hetc_kin[0] ->Write() ;
//          hetc_ori[1] ->Write() ;      hetc_pre[1] ->Write() ;      hetc_kin[1] ->Write() ;
//          hetc_ori[2] ->Write() ;      hetc_pre[2] ->Write() ;      hetc_kin[2] ->Write() ;
//          hetc_ori[3] ->Write() ;      hetc_pre[3] ->Write() ;      hetc_kin[3] ->Write() ;
//          hetc_ori[4] ->Write() ;      hetc_pre[4] ->Write() ;      hetc_kin[4] ->Write() ;
//          hetc_ori[5] ->Write() ;      hetc_pre[5] ->Write() ;      hetc_kin[5] ->Write() ;
//          hetc_ori[6] ->Write() ;      hetc_pre[6] ->Write() ;      hetc_kin[6] ->Write() ;
//          hetc_ori[7] ->Write() ;      hetc_pre[7] ->Write() ;      hetc_kin[7] ->Write() ;
//          outFile->Write();
//
//        cout << "## " << fileName
//             << " Type1 oriMt Integral = " << hetc_ori[2]->Integral()
//             << "    preMt Integral = " << hetc_pre[2]->Integral()
//             << "    kinMt Integral = " << hetc_kin[2]->Integral()
//             << "    # kinMt/preMt = " <<  ( (hetc_kin[2]->Integral())/(hetc_pre[2]->Integral())*100  ) <<  " %" << endl;
//        cout << "Type1XY  oriMt Integral = " << hetc_ori[3]->Integral()
//             << "    preMt Integral = " << hetc_pre[3]->Integral()
//             << "    kinMt Integral = " << hetc_kin[3]->Integral()
//             << "    # kinMt/preMt = " <<  ( (hetc_kin[3]->Integral())/(hetc_pre[3]->Integral())*100  ) <<  " %" << endl;
//
//
//        cout << "@@@ NumEvt Total       = " << TotalN     <<  endl ;
//      //  cout << "@@@ NumEvt Mu pt>50    = " << Nmupt50    << " / " << TotalN     << " = " << (double(Nmupt50   )/double(TotalN  ))*100 <<" %"  << endl ;
//        cout << "@@@ NumEvt METnoBadPF  = " << NBadPFMu << " / " << TotalN    << " = " << 100 - ((double(NBadPFMu)/double(TotalN ))*100) <<" %"  << endl ;
//        cout << "@@@ NumEvt METnoBadCh  = " << NnoBadch   << " / " << TotalN    << " = " << 100 - ((double(NnoBadch  )/double(TotalN ))*100) <<" %"  << endl ;
//        cout << "@@@ NumEvt MET good    = " << Nmetfil    << " / " << TotalN    << " = " << 100 - ((double(Nmetfil   )/double(TotalN ))*100) <<" %"  << endl ;
//        cout << "@@@ NumEvt METallPass  = " << NpassMET   << " / " << TotalN    << " = " << (double(NpassMET  )/double(TotalN ))*100 <<" %"  << endl ;
//        cout << "@@@ NumEvt HLT(Tk)Mu50 = " << Ntrigg     << " / " << NpassMET   << " = " << (double(Ntrigg    )/double(NpassMET))*100 <<" %"  << endl ;
//        cout << "@@@ NumEvt VetoManyMu  = " << NVetoManyMu   << " / " << Ntrigg     << " = " << (double(NVetoManyMu  )/double(Ntrigg  ))*100 <<" %"  << endl ;
//    //    cout << "@@@ NumEvt Veto2mu>25  = " << NVeto2nd   << " / " << Ntrigg     << " = " << (double(NVeto2nd  )/double(Ntrigg  ))*100 <<" %"  << endl ;
//        cout << "@@@ NumEvt 2HighPt>53  = " << NhighPtID2 << endl ;
//    //    cout << "@@@ NumEvt 1HighPt>53  = " << NhighPtID1 << " / " << NVeto2nd   << " = " << (double(NhighPtID1)/double(NVeto2nd  ))*100 << " %" << endl ;
//        cout << "@@@ NumEvt LooseTrkIso = " << NmuIso     << " / " << NhighPtID1 << " = " << (double(NmuIso    )/double(NhighPtID1))*100 << " %" << endl ;
//        cout << "@@@ NumEvt MuNoTkcut   = " << NmuonTkcut << endl ;
//        cout << "@@@ NumEvt MuNBadmu    = " << NBadmucut  << endl ;
//        cout << "@@@ NumEvt Mt>60       = " << Nmt60      << " / " << NmuIso     << " = " << (double(Nmt60   )/double(NmuIso  ))*100 << " %" << endl;
//        cout << "@@@ NumEvt kin ET/MET  = " << Nkinet     << " / " << Nmt60      << " = " << (double(Nkinet  )/double(Nmt60   ))*100 << " %" << endl;
//        cout << "@@@ NumEvt kin dPhi    = " << Nkindphi   << " / " << Nkinet     << " = " << (double(Nkindphi)/double(Nkinet  ))*100 << " %" << endl;
//        cout << "@@@ NumEvt Mt>80       = " << Nmt80      << " / " << Nkindphi   << " = " << (double(Nmt80   )/double(Nkindphi))*100 << " %" << endl;

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
        delete muIdTCA         ;

        //outFile->Close();
        //cout << "*** File Closed " << endl;

        return 0;
}
