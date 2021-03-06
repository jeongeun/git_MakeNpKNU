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
#include "./ProdNpKNU/src/NpKNU.hh"

   using namespace std;
   using namespace TMath;
   npknu::PUReweight* puWeight ;

//   const char *  mustr_ori[12] = {"mu_Pt_ori", "mu_Eta_ori","mu_Phi_ori", "mu_momQuality_ori","mu_TrkIsoR03_sumPt_ori", "mu_relTrkIsoR03_ori", "mu_nPixelHits_ori", "mu_nTrkLayers_ori", "mu_nMatchedStations_ori","mu_nValidMuonHits_ori", "mu_dz_ori", "mu_dxy_ori"};
//   const char *  mustr_pre[12] = {"mu_Pt_pre", "mu_Eta_pre","mu_Phi_pre", "mu_momQuality_pre","mu_TrkIsoR03_sumPt_pre", "mu_relTrkIsoR03_pre", "mu_nPixelHits_pre", "mu_nTrkLayers_pre", "mu_nMatchedStations_pre","mu_nValidMuonHits_pre", "mu_dz_pre", "mu_dxy_pre"};
//   const char *  mustr_kin[12] = {"mu_Pt_kin", "mu_Eta_kin","mu_Phi_kin", "mu_momQuality_kin","mu_TrkIsoR03_sumPt_kin", "mu_relTrkIsoR03_kin", "mu_nPixelHits_kin", "mu_nTrkLayers_kin", "mu_nMatchedStations_kin","mu_nValidMuonHits_kin", "mu_dz_kin", "mu_dxy_kin"};
//
//   Int_t mubin[12]      = { 80,    50,  50,  50,  50,  50,  15, 30, 30, 100, 60,   60 };
//   Double_t muxmin[12]  = { 0.0, -2.5,-3.2, 0.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0,-0.6,-0.3 };
//   Double_t muxmax[12]  = {2000,  2.5, 3.2, 0.5, 100, 0.2,  15, 30, 30, 100, 0.6, 0.3 };
//   TString muxTitle[12] = {"Muon P_{T} [GeV]", "Muon #eta", "Muon #phi", "#sigma_{Pt}/P_{T}(<0.3)","#Sigmap_{T}^{Trk}(#DeltaR<0.3)","#Sigmap_{T}^{Trk}(#DeltaR<0.3)/p_{T}(#mu)(<0.1)", "nPixelHits(>0)","nTrackerLayers(>5)","nMatchedStations(>1)","nValidMuonHits(>0)","dZ(<5mm)","dXY(<2mm)"};
//
//   const char * metstr_ori[7] = {"met_Pt_ori","met_Phi_ori","met_sumEt_ori", "met_unCorr_ori","met_unCorrPhi_ori", "met_Puppi_ori","met_PuppiPhi_ori"}  ;
//   const char * metstr_pre[7] = {"met_Pt_pre","met_Phi_pre","met_sumEt_pre", "met_unCorr_pre","met_unCorrPhi_pre", "met_Puppi_pre","met_PuppiPhi_pre"}  ;
//   const char * metstr_kin[7] = {"met_Pt_kin","met_Phi_kin","met_sumEt_kin", "met_unCorr_kin","met_unCorrPhi_kin", "met_Puppi_kin","met_PuppiPhi_kin"}  ;
//   
//   Int_t metbin[7]      = { 80,   50, 100, 80,  50,   80,  50}  ;
//   Double_t metxmin[7]  = { 0.0,-3.2, 0.0, 0.0,-3.2,  0.0,-3.2}  ;
//   Double_t metxmax[7]  = {2000, 3.2,3000,2000, 3.2, 2000, 3.2}  ;
//   TString metxTitle[7] = {"MET(MuEGClean) [GeV]", "MET(MuEGClean) #phi", "sumEt(MuEGClean)", "MET(UnCorr) [GeV]", "MET (UnCorr) #phi","PUPPI MET [GeV]", "PUPPI MET #phi" }  ;
//
//   const char *  etcstr_ori[10] = {"nPVwTest_ori","invM_ori","nPV_ori","nPVw_ori","wp_Mt_ori","wp_MtUnCorr_ori","wp_Dphi_ori","wp_DphiUncorr_ori","wp_EtRatio_ori","wp_EtRatioUnCorr_ori"}  ;
//   const char *  etcstr_pre[10] = {"nPVwTest_pre","invM_pre","nPV_pre","nPVw_pre","wp_Mt_pre","wp_MtUnCorr_pre","wp_Dphi_pre","wp_DphiUncorr_pre","wp_EtRatio_pre","wp_EtRatioUnCorr_pre"}  ;
//   const char *  etcstr_kin[10] = {"nPVwTest_kin","invM_kin","nPV_kin","nPVw_kin","wp_Mt_kin","wp_MtUnCorr_kin","wp_Dphi_kin","wp_DphiUncorr_kin","wp_EtRatio_kin","wp_EtRatioUnCorr_kin"}  ;
//
//   Int_t etcbin[10]      ={100,   40, 100,  100,   100,  100,   50,  50,  50,  50 }   ;
//   Double_t etcxmin[10]  ={0.0, 60.0, 0.0,  0.0,   0.0,  0.0,  0.0, 0.0, 0.0, 0.0 }  ;
//   Double_t etcxmax[10]  ={100, 1400, 100,  100,  5000, 5000,  5.0, 5.0,  10,  10 }  ;
//   TString etcxTitle[10] ={"nPV_Test(w PUreweight)", "inv M_{#mu#mu}","nPV(w/o PUreWeight)","nPV(w PUreweight)" ,"M_{T} [GeV]","M_{T}(UnCorr) [GeV]","#Delta #phi","#Delta #phi(UnCorr)","E_{T}^{#mu}/E_{T}^{miss}","E_{T}^{#mu}/E_{T}^{miss}(UnCorr)"}  ;
//
//   TH1F* hmu_ori[12] ;
//   TH1F* hmu_pre[12] ;
//   TH1F* hmu_kin[12] ;
//
//   TH1F* hmet_ori[7] ;
//   TH1F* hmet_pre[7] ;
//   TH1F* hmet_kin[7] ;
//
//   TH1F* hetc_ori[10] ;
//   TH1F* hetc_pre[10] ;
//   TH1F* hetc_kin[10] ;

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
        puWeight = new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_RunBtoF_69p2mb.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
       // puWeight = new npknu::PUReweight("/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/etc/pileup/pu_RunGtoH_69p2mb.root" , argv[5] , "MakeNpKNUHist/pileupHist/h1_TrueNumInteractions100");
     }
     cout << "Type " << type << ": lumiWeight " << lumiWeight << "=" << lumi << "*" << xsec << "/" << gen << "  : " << fileName << endl;
//     for(int j=0; j<12; j++){ //Draw Variables
//        hmu_ori[j] = new TH1F(mustr_ori[j],mustr_ori[j], mubin[j], muxmin[j], muxmax[j]); hmu_ori[j]->GetXaxis()->SetTitle(muxTitle[j]);  hmu_ori[j]->Sumw2();
//        hmu_pre[j] = new TH1F(mustr_pre[j],mustr_pre[j], mubin[j], muxmin[j], muxmax[j]); hmu_pre[j]->GetXaxis()->SetTitle(muxTitle[j]);  hmu_pre[j]->Sumw2();
//        hmu_kin[j] = new TH1F(mustr_kin[j],mustr_kin[j], mubin[j], muxmin[j], muxmax[j]); hmu_kin[j]->GetXaxis()->SetTitle(muxTitle[j]);  hmu_kin[j]->Sumw2();
//     }
//     for(int k=0; k<7; k++){ //Draw Variables
//        hmet_ori[k] = new TH1F(metstr_ori[k],metstr_ori[k], metbin[k], metxmin[k], metxmax[k]); hmet_ori[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_ori[k]->Sumw2();
//        hmet_pre[k] = new TH1F(metstr_pre[k],metstr_pre[k], metbin[k], metxmin[k], metxmax[k]); hmet_pre[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_pre[k]->Sumw2();
//        hmet_kin[k] = new TH1F(metstr_kin[k],metstr_kin[k], metbin[k], metxmin[k], metxmax[k]); hmet_kin[k]->GetXaxis()->SetTitle(metxTitle[k]);  hmet_kin[k]->Sumw2();
//     }
//     for(int p=0; p<10; p++){ //Draw Variables
//        hetc_ori[p] = new TH1F(etcstr_ori[p],etcstr_ori[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_ori[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_ori[p]->Sumw2();
//        hetc_pre[p] = new TH1F(etcstr_pre[p],etcstr_pre[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_pre[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_pre[p]->Sumw2();
//        hetc_kin[p] = new TH1F(etcstr_kin[p],etcstr_kin[p], etcbin[p], etcxmin[p], etcxmax[p]); hetc_kin[p]->GetXaxis()->SetTitle(etcxTitle[p]);  hetc_kin[p]->Sumw2();
//     }
//     TFile* outFile = new TFile(argv[6],"recreate");
     int TotalN = (int)inChain->GetEntries();
     int per99 = (TotalN > 10) ? TotalN / 9 : 1;
     int per100 = 0;

     std::vector<std::string> trigNameVec;
     trigNameVec.push_back("HLT_Mu50_v*");
     trigNameVec.push_back("HLT_TkMu50_v*");
     TClonesArray* trigNamePassTCA = new TClonesArray("npknu::Trigger");

   bool isMETFilter = false;
   bool isMETBadPFMuon = false;
   bool isMETBadChMuon = false;
   bool isMuPreselection = false;

   bool isCorrectedGlobalMu = false;
   bool isCorrectedDupliMu = false;
   int NRunLumiPass       =0 ;
   int Nmetfil            = 0 ; 
   int NBadPFMu           = 0 ; 
   int NnoBadch           = 0 ; 
   int NCorrectedGBMu     = 0 ; 
   int NCorrectedDuMu     = 0 ; 
   int NpassMET           = 0 ; 
   int Nmupt50            = 0 ; 
   int Ntrigg             = 0 ; 
   int NHighPtExist       = 0 ; 
   int NHighPtPassTrk     = 0 ; 
   int NHighPt2           = 0 ; 
   int NDimuonZ           = 0 ;
   int NHighPt2ButBig     = 0 ; 
   int NCloseR            = 0 ; 
   int NSoftMuPassVeto25  = 0 ; 
   int NHighPt2PassVeto25 = 0 ; 
   int NMuonSelection     = 0 ; 
   int NMt60Cut           = 0 ; 
   int MKinCutPass        = 0 ; 

     for(int eventLoop=0; eventLoop<TotalN ; eventLoop++) {
        inChain->GetEntry(eventLoop);
        if((eventLoop%per99) == 0) cout << "Running " << (per100++ * 10)<< " % " << eventLoop << " / " << TotalN << endl;
        muIdTCA->Clear();
        npknu::Evt* evtPtr = (npknu::Evt*)evtTCA->At(0);
        if((evtPtr->run != 282037) or (evtPtr->lumi != 359) or (evtPtr->event != 663395404)) continue;
        NRunLumiPass++;

        double pileup       = (evtPtr->isRealData) ? -999 : npknu::GetPileupTrue(pileupTCA);
        double pileupWeight = (evtPtr->isRealData) ?  1.0 : puWeight->getWeight(pileup);
        double tWeight      = pileupWeight * lumiWeight; // Notice:We need additional weight for gen failed event (Ngen/NGencut)
        int nVtx = vertexTCA->GetEntries();

//        hetc_ori[0] ->Fill(nVtx  ,tWeight ) ;

        /*************************************/
        /*Both Data and MC MET Filters passed*/
        /*************************************/
        npknu::MET* metPtr  = (npknu::MET*)metTCA->At(0);

        isMETFilter = ( metPtr->Flag_HBHENoiseFilter == 1 && 
                        metPtr->Flag_HBHENoiseIsoFilter == 1 && 
                        metPtr->Flag_goodVertices ==1 && 
                        metPtr->Flag_CSCTightHaloFilter ==1 && 
                        metPtr->Flag_EcalDeadCellTriggerPrimitiveFilter ==1 && 
                        metPtr->Flag_eeBadScFilter ==1 && 
                       (metPtr->Flag_globalSuperTightHalo2016Filter ==1  ||  metPtr->Flag_globalTightHalo2016Filter ==1) ) ;

        if (!isMETFilter) Nmetfil ++ ; ////Count BadMET

        isMETBadPFMuon = ( metPtr->IsNoBadPFMuon == 0) ;
        if (isMETBadPFMuon) NBadPFMu ++; //Count BadMET

        isMETBadChMuon = ( metPtr->IsNoBadChargedCandidate == 0 ) ;
        if(isMETBadChMuon) NnoBadch ++; 
        //MC17 case if(metPtr->IsBadGlobalMuonTagger == 1){NBadGlobalMuon++;}//Count BadGBMu

        isCorrectedGlobalMu = (  metPtr->IsCorrectedBadGlobalMuon == 1 ) ;//goodEvent
        if(isCorrectedGlobalMu ) NCorrectedGBMu ++;

        isCorrectedDupliMu = (  metPtr->IsCorrectedDuplicateMuon == 1 ) ;//goodEvent
        if(isCorrectedDupliMu ) NCorrectedDuMu ++;

        cout << "==> RunH " << evtPtr->run << " Event " << long(evtPtr->event) << " Lumi " << evtPtr->lumi << endl;
        cout << "isGoodMET isBadPFMuon isBadChMuon isCorrectedGlobalMu isCorrectedDupliMu "  << isMETFilter << " " <<isMETBadPFMuon <<" " << isMETBadChMuon << " " << isCorrectedGlobalMu << isCorrectedDupliMu << endl;

        if(!isMETFilter || isMETBadPFMuon || isMETBadChMuon ) continue;
        NpassMET ++ ; 
        cout << " PASS METFilter " << endl;
        metPtr->PrintMET();

        bool isMuPreselection = npknu::GetMuonPreselection(muonTCA, 53.0, 0.0, 2.4); /*highPtMuonGVId, ptCut, etaCut*/
        if(isMuPreselection) Nmupt50++; //////////////////////////
        if(!isMuPreselection) continue;
//        hetc_pre[0] ->Fill(nVtx  ,tWeight ) ;
        cout << " PASS Muon Pt > 53GeV " << endl;

        /****************************************************************/
        /*SingleMuon HLT Filter to Data and MC(17) (Mu50 OR TkMu50)     */
        /****************************************************************/
        npknu::GetPassTriggerTCA(triggerTCA, trigNamePassTCA, trigNameVec);
        if(trigNamePassTCA->GetEntries() == 0) continue;
        cout << " PASS HLT_(Tk)Mu50_v* " << endl;
        Ntrigg++ ;
//        hetc_kin[0] ->Fill(nVtx  ,tWeight ) ;
        /******************************************/
        /*Both Data HighPtMuID &Loose TrkIso <0.1 */
        /******************************************/
       bool highPtTrkPass1 = false;
       bool highPtTrkPass2 = false;
       bool Veto2ndhighPtPass = false;
       bool Veto2ndsoftPass = false;
       bool noSoftbiggermuon = true ;

        double dR = 0.0 ;
        npknu::GetMuonIdTCA(muonTCA, muIdTCA, 6, 53.0, 0.0, 2.4); /*highPtMuonGVId, ptCut, etaCut*/
        if(muIdTCA->GetEntries() ==0) continue;
        
        NHighPtExist++;
        npknu::Muon* highmuPtr0 = (npknu::Muon*)muIdTCA->At(0);
        cout << " PASS Muon HighPtID Exist " << endl;
        highmuPtr0->PrintMuon();
        if(highmuPtr0->relTrkIsoR03_BT <0.1 ){highPtTrkPass1 = true; cout << " PASS Muon HighPtID+Isolation " << highmuPtr0->pt << "GeV " << endl;}
        if(muIdTCA->GetEntries() > 1){
           NHighPt2++;
           npknu::Muon* highmuPtr1 = (npknu::Muon*)muIdTCA->At(1);
         cout << " SecondMuon " << highmuPtr1->pt << " GeV" ; highmuPtr1->PrintMuon();
           if(highmuPtr1->relTrkIsoR03_BT <0.1 and highPtTrkPass1 and highmuPtr1->pt > 25.0 and highmuPtr0->charge * highmuPtr1->charge < 1 ){
               NDimuonZ++;
//              hetc_ori[1]->Fill(highmuPtr0->GetM(highmuPtr1), lumiWeight );
//              hetc_pre[1]->Fill(highmuPtr0->GetM(highmuPtr1), tWeight    );
//              hetc_kin[1]->Fill(highmuPtr0->GetM(highmuPtr1), tWeight    );
           }
           if(highPtTrkPass1 and highmuPtr1->pt < 25.0){Veto2ndhighPtPass = true ; cout << "PASS Veto2nd Muon " << highmuPtr1->pt << "GeV " <<endl;}
        }
	for(int i=0; i< muonTCA->GetEntries(); i++){
           npknu::Muon* muPtr = (npknu::Muon*)muonTCA->At(i);
           cout << " all muon info " << i << " HighPtID "<<muPtr->isHighPtMuonGV << "LooseId " << muPtr->isLooseMuon << "SoftId " << muPtr->isSoftMuonGV  ;  muPtr->PrintMuon();
           if(muPtr->isHighPtMuonGV == 1) continue;
           if(muPtr->pt > 25.0){ 
              NHighPt2ButBig++;
              noSoftbiggermuon = false;
              if(highmuPtr0->deltaR(muPtr) < 0.3){ NCloseR++;}
           }
        }       
       bool OnehighPtMu = false;
       bool Pass2ndHighPtsmallMu = false;
       if(noSoftbiggermuon == true ){ NSoftMuPassVeto25++; }   
       if(muIdTCA->GetEntries() == 1 and highPtTrkPass1 == 1 and noSoftbiggermuon == true){OnehighPtMu = true ; NHighPtPassTrk++; } 
       if(muIdTCA->GetEntries() > 1  and Veto2ndhighPtPass == 1 and noSoftbiggermuon == true){Pass2ndHighPtsmallMu = true ; NHighPt2PassVeto25++; } 
  
       if(!OnehighPtMu and !Pass2ndHighPtsmallMu)continue; 
       NMuonSelection++;
       cout << " PASS Muon Selection HighPt, Isolation, 2ndMu Veto " << endl; 
        cout << "Mt " << highmuPtr0->GetMt(metPtr) << "GeV, mu " << highmuPtr0->pt << ", met " << metPtr->pt << ", et/MET " << highmuPtr0->etRatio(metPtr) << ", dPhi " << highmuPtr0->deltaPhi(metPtr) << endl;
       // /* Muon Plots */
       // hmu_ori[0] ->Fill(highmuPtr0->pt    ,lumiWeight  )  ;
       // hmu_ori[1] ->Fill(highmuPtr0->eta   ,lumiWeight  )  ;
       // hmu_ori[2] ->Fill(highmuPtr0->phi   ,lumiWeight  )  ;
       // hmu_ori[3] ->Fill( ((highmuPtr0->ptErr_BT)/(highmuPtr0->pt_BT)) ,lumiWeight  )  ;
       // hmu_ori[4] ->Fill(highmuPtr0->TrkIsoR03_sumPt    ,lumiWeight  )  ;
       // hmu_ori[5] ->Fill(highmuPtr0->relTrkIsoR03_BT    ,lumiWeight  )  ;
       // hmu_ori[6] ->Fill(highmuPtr0->nValidPixelHits    ,lumiWeight  )  ;
       // hmu_ori[7] ->Fill(highmuPtr0->nTrackerLayers     ,lumiWeight  )  ;
       // hmu_ori[8] ->Fill(highmuPtr0->nMatchedStations   ,lumiWeight  )  ;
       // hmu_ori[9] ->Fill(highmuPtr0->nValidMuonHits     ,lumiWeight  )  ;
       // hmu_ori[10] ->Fill(highmuPtr0->dz_GV_BT           ,lumiWeight  )  ;
       // hmu_ori[11]->Fill(highmuPtr0->dxy_GV_BT          ,lumiWeight  )  ;
       // /* MET Plots */
       // hmet_ori[0] ->Fill(metPtr->pt     ,lumiWeight  ) ;
       // hmet_ori[1] ->Fill(metPtr->phi    ,lumiWeight  ) ;
       // hmet_ori[2] ->Fill(metPtr->sumEt  ,lumiWeight  ) ;
       // hmet_ori[3] ->Fill(metPtr->Uncorr_pt ,lumiWeight  ) ;
       // hmet_ori[4] ->Fill(metPtr->Uncorr_phi ,lumiWeight  ) ;
       // hmet_ori[5] ->Fill(metPtr->METPuppi_pt  ,lumiWeight  ) ;
       // hmet_ori[6] ->Fill(metPtr->METPuppi_phi ,lumiWeight  ) ;
       // /* kinematic Plots */
       // hetc_ori[2] ->Fill(nVtx  ,lumiWeight ) ;
       // hetc_ori[3] ->Fill(nVtx  ,tWeight )  ;
       // hetc_ori[4] ->Fill(highmuPtr0->GetMt(metPtr)  ,lumiWeight )  ;
       // hetc_ori[5] ->Fill(highmuPtr0->GetMt(metPtr->Uncorr_pt, metPtr->Uncorr_phi) ,lumiWeight );
       // hetc_ori[6] ->Fill(highmuPtr0->deltaPhi(metPtr)  ,lumiWeight )  ;
       // hetc_ori[7] ->Fill(highmuPtr0->deltaPhi(metPtr->Uncorr_phi)  ,lumiWeight )  ;
       // hetc_ori[8] ->Fill(highmuPtr0->etRatio(metPtr)  ,lumiWeight )  ;
       // hetc_ori[9] ->Fill(highmuPtr0->etRatio(metPtr->Uncorr_pt)  ,lumiWeight )  ;

        if(highmuPtr0->GetMt(metPtr) < 60.0) continue;
        cout << " PASS Mt >= 60GeV cut " << endl;
        cout << "Mt " << highmuPtr0->GetMt(metPtr) << "GeV, mu " << highmuPtr0->pt << ", met " << metPtr->pt << ", et/MET " << highmuPtr0->etRatio(metPtr) << ", dPhi " << highmuPtr0->deltaPhi(metPtr) << endl;
        NMt60Cut++;
        /* Muon Plots */
       // hmu_pre[0] ->Fill(highmuPtr0->pt    ,tWeight  )  ;
       // hmu_pre[1] ->Fill(highmuPtr0->eta   ,tWeight  )  ;
       // hmu_pre[2] ->Fill(highmuPtr0->phi   ,tWeight  )  ;
       // hmu_pre[3] ->Fill( ((highmuPtr0->ptErr_BT)/(highmuPtr0->pt_BT)) ,tWeight  )  ;
       // hmu_pre[4] ->Fill(highmuPtr0->TrkIsoR03_sumPt    ,tWeight  )  ;
       // hmu_pre[5] ->Fill(highmuPtr0->relTrkIsoR03_BT    ,tWeight  )  ;
       // hmu_pre[6] ->Fill(highmuPtr0->nValidPixelHits    ,tWeight  )  ;
       // hmu_pre[7] ->Fill(highmuPtr0->nTrackerLayers     ,tWeight  )  ;
       // hmu_pre[8] ->Fill(highmuPtr0->nMatchedStations   ,tWeight  )  ;
       // hmu_pre[9] ->Fill(highmuPtr0->nValidMuonHits     ,tWeight  )  ;
       // hmu_pre[10] ->Fill(highmuPtr0->dz_GV_BT           ,tWeight  )  ;
       // hmu_pre[11]->Fill(highmuPtr0->dxy_GV_BT          ,tWeight  )  ;
       // /* MET Plots */
       // hmet_pre[0] ->Fill(metPtr->pt     ,tWeight  ) ;
       // hmet_pre[1] ->Fill(metPtr->phi    ,tWeight  ) ;
       // hmet_pre[2] ->Fill(metPtr->sumEt  ,tWeight  ) ;
       // hmet_pre[3] ->Fill(metPtr->Uncorr_pt ,tWeight  ) ;
       // hmet_pre[4] ->Fill(metPtr->Uncorr_phi ,tWeight  ) ;
       // hmet_pre[5] ->Fill(metPtr->METPuppi_pt  ,tWeight  ) ;
       // hmet_pre[6] ->Fill(metPtr->METPuppi_phi ,tWeight  ) ;
       // /* kinematic Plots */
       // hetc_pre[2] ->Fill(nVtx  ,lumiWeight ) ;
       // hetc_pre[3] ->Fill(nVtx  ,tWeight )  ;
       // hetc_pre[4] ->Fill(highmuPtr0->GetMt(metPtr)  ,tWeight )  ;
       // hetc_pre[5] ->Fill(highmuPtr0->GetMt(metPtr->Uncorr_pt, metPtr->Uncorr_phi) ,tWeight );
       // hetc_pre[6] ->Fill(highmuPtr0->deltaPhi(metPtr)  ,tWeight )  ;
       // hetc_pre[7] ->Fill(highmuPtr0->deltaPhi(metPtr->Uncorr_phi)  ,tWeight )  ;
       // hetc_pre[8] ->Fill(highmuPtr0->etRatio(metPtr)  ,tWeight )  ;
       // hetc_pre[9] ->Fill(highmuPtr0->etRatio(metPtr->Uncorr_pt)  ,tWeight )  ;

        if(highmuPtr0->GetMt(metPtr) < 80.0) continue;
        cout << " PASS Mt >= 80GeV cut " << endl;
        if(highmuPtr0->etRatio(metPtr)<= 0.4 or highmuPtr0->etRatio(metPtr)>= 1.5) continue ;
        cout << " PASS 0.4< et/MET <1.5 cut " << endl;
        if(highmuPtr0->deltaPhi(metPtr->phi) <=2.5)  continue ;
        cout << " PASS dPhi > 2.5 cut " << endl;
        MKinCutPass++;

        cout << "==> RunH " << evtPtr->run << " Event " << long(evtPtr->event) << " Lumi " << evtPtr->lumi << endl;
        cout << "Mt " << highmuPtr0->GetMt(metPtr) << "GeV, mu " << highmuPtr0->pt << ", met " << metPtr->pt << ", et/MET " << highmuPtr0->etRatio(metPtr) << ", dPhi " << highmuPtr0->deltaPhi(metPtr) << endl;
        highmuPtr0->PrintMuon();
        metPtr->PrintMET();

       // /* Muon Plots */
       // hmu_kin[0] ->Fill(highmuPtr0->pt    ,tWeight  )  ;
       // hmu_kin[1] ->Fill(highmuPtr0->eta   ,tWeight  )  ;
       // hmu_kin[2] ->Fill(highmuPtr0->phi   ,tWeight  )  ;
       // hmu_kin[3] ->Fill( ((highmuPtr0->ptErr_BT)/(highmuPtr0->pt_BT)) ,tWeight  )  ;
       // hmu_kin[4] ->Fill(highmuPtr0->TrkIsoR03_sumPt    ,tWeight  )  ;
       // hmu_kin[5] ->Fill(highmuPtr0->relTrkIsoR03_BT    ,tWeight  )  ;
       // hmu_kin[6] ->Fill(highmuPtr0->nValidPixelHits    ,tWeight  )  ;
       // hmu_kin[7] ->Fill(highmuPtr0->nTrackerLayers     ,tWeight  )  ;
       // hmu_kin[8] ->Fill(highmuPtr0->nMatchedStations   ,tWeight  )  ;
       // hmu_kin[9] ->Fill(highmuPtr0->nValidMuonHits     ,tWeight  )  ;
       // hmu_kin[10] ->Fill(highmuPtr0->dz_GV_BT           ,tWeight  )  ;
       // hmu_kin[11]->Fill(highmuPtr0->dxy_GV_BT          ,tWeight  )  ;
       // /* MET Plots */
       // hmet_kin[0] ->Fill(metPtr->pt     ,tWeight  ) ;
       // hmet_kin[1] ->Fill(metPtr->phi    ,tWeight  ) ;
       // hmet_kin[2] ->Fill(metPtr->sumEt  ,tWeight  ) ;
       // hmet_kin[3] ->Fill(metPtr->Uncorr_pt ,tWeight  ) ;
       // hmet_kin[4] ->Fill(metPtr->Uncorr_phi ,tWeight  ) ;
       // hmet_kin[5] ->Fill(metPtr->METPuppi_pt  ,tWeight  ) ;
       // hmet_kin[6] ->Fill(metPtr->METPuppi_phi ,tWeight  ) ;
       // /* kinematic Plots */
       // hetc_kin[2] ->Fill(nVtx  ,lumiWeight ) ;
       // hetc_kin[3] ->Fill(nVtx  ,tWeight )  ;
       // hetc_kin[4] ->Fill(highmuPtr0->GetMt(metPtr)  ,tWeight )  ;
       // hetc_kin[5] ->Fill(highmuPtr0->GetMt(metPtr->Uncorr_pt, metPtr->Uncorr_phi) ,tWeight );
       // hetc_kin[6] ->Fill(highmuPtr0->deltaPhi(metPtr)  ,tWeight )  ;
       // hetc_kin[7] ->Fill(highmuPtr0->deltaPhi(metPtr->Uncorr_phi)  ,tWeight )  ;
       // hetc_kin[8] ->Fill(highmuPtr0->etRatio(metPtr)  ,tWeight )  ;
       // hetc_kin[9] ->Fill(highmuPtr0->etRatio(metPtr->Uncorr_pt)  ,tWeight )  ;
  }//event loop


          //hmu_ori[0] ->Write() ;      hmu_pre[0] ->Write() ;      hmu_kin[0] ->Write() ;
          //hmu_ori[1] ->Write() ;      hmu_pre[1] ->Write() ;      hmu_kin[1] ->Write() ;
          //hmu_ori[2] ->Write() ;      hmu_pre[2] ->Write() ;      hmu_kin[2] ->Write() ;
          //hmu_ori[3] ->Write() ;      hmu_pre[3] ->Write() ;      hmu_kin[3] ->Write() ;
          //hmu_ori[4] ->Write() ;      hmu_pre[4] ->Write() ;      hmu_kin[4] ->Write() ;
          //hmu_ori[5] ->Write() ;      hmu_pre[5] ->Write() ;      hmu_kin[5] ->Write() ;
          //hmu_ori[6] ->Write() ;      hmu_pre[6] ->Write() ;      hmu_kin[6] ->Write() ;
          //hmu_ori[7] ->Write() ;      hmu_pre[7] ->Write() ;      hmu_kin[7] ->Write() ;
          //hmu_ori[8] ->Write() ;      hmu_pre[8] ->Write() ;      hmu_kin[8] ->Write() ;
          //hmu_ori[9] ->Write() ;      hmu_pre[9] ->Write() ;      hmu_kin[9] ->Write() ;
          //hmu_ori[10] ->Write() ;     hmu_pre[10] ->Write() ;     hmu_kin[10] ->Write() ;
          //hmu_ori[11] ->Write() ;     hmu_pre[11] ->Write() ;     hmu_kin[11] ->Write() ;

          //hmet_ori[0] ->Write() ;     hmet_pre[0] ->Write() ;      hmet_kin[0] ->Write() ;
          //hmet_ori[1] ->Write() ;     hmet_pre[1] ->Write() ;      hmet_kin[1] ->Write() ;
          //hmet_ori[2] ->Write() ;     hmet_pre[2] ->Write() ;      hmet_kin[2] ->Write() ;
          //hmet_ori[3] ->Write() ;     hmet_pre[3] ->Write() ;      hmet_kin[3] ->Write() ;
          //hmet_ori[4] ->Write() ;     hmet_pre[4] ->Write() ;      hmet_kin[4] ->Write() ;
          //hmet_ori[5] ->Write() ;     hmet_pre[5] ->Write() ;      hmet_kin[5] ->Write() ;
          //hmet_ori[6] ->Write() ;     hmet_pre[6] ->Write() ;      hmet_kin[6] ->Write() ;

          //hetc_ori[0] ->Write() ;     hetc_pre[0] ->Write() ;      hetc_kin[0] ->Write() ;
          //hetc_ori[1] ->Write() ;     hetc_pre[1] ->Write() ;      hetc_kin[1] ->Write() ;
          //hetc_ori[2] ->Write() ;     hetc_pre[2] ->Write() ;      hetc_kin[2] ->Write() ;
          //hetc_ori[3] ->Write() ;     hetc_pre[3] ->Write() ;      hetc_kin[3] ->Write() ;
          //hetc_ori[4] ->Write() ;     hetc_pre[4] ->Write() ;      hetc_kin[4] ->Write() ;
          //hetc_ori[5] ->Write() ;     hetc_pre[5] ->Write() ;      hetc_kin[5] ->Write() ;
          //hetc_ori[6] ->Write() ;     hetc_pre[6] ->Write() ;      hetc_kin[6] ->Write() ;
          //hetc_ori[7] ->Write() ;     hetc_pre[7] ->Write() ;      hetc_kin[7] ->Write() ;
          //hetc_ori[8] ->Write() ;     hetc_pre[8] ->Write() ;      hetc_kin[8] ->Write() ;
          //hetc_ori[9] ->Write() ;     hetc_pre[9] ->Write() ;      hetc_kin[9] ->Write() ;
          //outFile->Write();

        //cout << "## Muon Barrel 1.2 < |eta| < 2.4" << fileName
        cout << "## Muon Total 0.0 < |eta| < 2.4" << fileName << endl;
             //<< " Type1 oriMt Integral = " << hetc_ori[4]->Integral()
             //<< "    preMt Integral = " << hetc_pre[4]->Integral()
             //<< "    kinMt Integral = " << hetc_kin[4]->Integral()
             //<< "    # kinMt/preMt = " <<  ( (hetc_kin[4]->Integral())/(hetc_pre[4]->Integral())*100  ) <<  " %" << endl;

        cout << "@Num of TotEvent      = " << setw(12) << TotalN             << " : " << setw(10) <<  double(TotalN )/double(TotalN )   <<  endl ;
        cout << "@Num of RunLumiPass   = " << setw(12) << NRunLumiPass       << " : " << setw(10) <<  double(NRunLumiPass )/double(TotalN )   <<  endl ;
        cout << " BadNum of BadMet     = " << setw(12) << Nmetfil            << " : " << setw(10) <<  double(Nmetfil  )/double(TotalN ) << " Filtered in tot" << endl ;
        cout << " BadNum of BadPFMu    = " << setw(12) << NBadPFMu           << " : " << setw(10) <<  double(NBadPFMu )/double(TotalN ) << " Filtered in tot" << endl ;
        cout << " BadNum of BadChCan   = " << setw(12) << NnoBadch           << " : " << setw(10) <<  double(NnoBadch )/double(TotalN ) << " Filtered in tot" << endl ;
        cout << " Num of CorrGBMu      = " << setw(12) << NCorrectedGBMu     << " : " << setw(10) <<  double(NCorrectedGBMu)/double(TotalN) << " Corrected in tot" << endl;
        cout << " Num of CorrDuMu      = " << setw(12) << NCorrectedDuMu     << " : " << setw(10) <<  double(NCorrectedDuMu)/double(TotalN) << " Corrected in tot" << endl;
        cout << "@Num of passGoodMET   = " << setw(12) << NpassMET           << " : " << setw(10) <<  double(NpassMET)/double(TotalN ) << " Survived in tot" << endl ;
        cout << "@Num of Mupt53        = " << setw(12) << Nmupt50            << " : " << setw(10) <<  double(Nmupt50 )/double(NpassMET ) << " Survived in NpassMET" << endl ;
        cout << "@Num of HLT           = " << setw(12) << Ntrigg             << " : " << setw(10) <<  double(Ntrigg  )/double(Nmupt50) << " Survived in NMupt53" << endl ;
        cout << " Num of HighPtExist   = " << setw(12) << NHighPtExist       << " : " << setw(10) <<  double(NHighPtExist   )/double(Ntrigg) << " Survived in Ntrigg" << endl ;
        cout << "@Num of HighPt1PassTrk= " << setw(12) << NHighPtPassTrk     << " : " << setw(10) <<  double(NHighPtPassTrk )/double(Ntrigg) << " Survived in Ntrigg" << endl ;
        cout << " Num of HighPt2       = " << setw(12) << NHighPt2           << " : " << setw(10) <<  double(NHighPt2)/double(NHighPtExist) << " Exist in NHithPtExist" << endl ;
        cout << " Num of HighPt2DimuZ  = " << setw(12) << NDimuonZ           << " : " << setw(10) <<  double(NDimuonZ)/double(NHighPtExist) << " Exist in NHithPtExist" << endl ;
        cout << " BadNum of HighPt2Big = " << setw(12) << NHighPt2ButBig     << " : " << setw(10) <<  double(NHighPt2ButBig)/double(NHighPt2) << " Exist in HighPt2" << endl ;
        cout << " BadNum of Close dR03 = " << setw(12) << NCloseR            << " : " << setw(10) <<  double(NCloseR)/double(NHighPt2ButBig) << " Exist in HighPt2ButBig" << endl ;
        cout << " Num of SoftMuPass25  = " << setw(12) << NSoftMuPassVeto25  << " : " << setw(10) <<  double(NSoftMuPassVeto25)/double(NHighPtPassTrk) << " Passed in HighPtExist" << endl ;
        cout << " Num of HighPt2Pass25 = " << setw(12) << NHighPt2PassVeto25 << " : " << setw(10) <<  double(NHighPt2PassVeto25)/double(NHighPt2) << " Passed in HighPt2" << endl;
        cout << "@Num of MuIDISOVeto2mu= " << setw(12) << NMuonSelection     << " : " << setw(10) <<  double(NMuonSelection)/double(Ntrigg) << " Survived in Ntrigg" << endl ;
        cout << " Num of Mt>60GeV      = " << setw(12) << NMt60Cut           << " : " << setw(10) <<  double(NMt60Cut)/double(NMuonSelection) << " Survived in MuonIDISOVeto2" << endl ;
        cout << "@Num of Mt>80 KinCut  = " << setw(12) << MKinCutPass        << " : " << setw(10) <<  double(MKinCutPass)/double(NMt60Cut) << " Survived in Mt60" << endl ;

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

//        outFile->Close();
//        cout << "*** File Closed " << endl;

        return 0;
}
