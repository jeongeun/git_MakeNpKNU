using namespace std;
using namespace TMath;

void Draw_17_0301(){
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  const int NumFiles = 13 ;
  const int NumBGs= 10 ;
  const int NumSigs= 2;

  bool isBF = false;
  bool isGH = true;
  bool isFull = false;

  TFile* openFile[NumFiles];
  TString inFileName[NumFiles];
  inFileName[0] = "../condorOut/final/hData_RunGH_reMiniaod.root" ;
  inFileName[1] = "../condorOut/final/hMC_W_M100_GH.root"         ;
  inFileName[2] = "../condorOut/final/hMC_W_M200_GH.root"         ;
  inFileName[3] = "../condorOut/final/hMC_W_M500_GH.root"         ;
  inFileName[4] = "../condorOut/final/hMC_W_M1000_GH.root"        ;
  inFileName[5] = "../condorOut/final/hMC_W_M2000_GH.root"        ;
  inFileName[6] = "../condorOut/final/hMC_W_M3000_GH.root"        ;
  inFileName[7] = "../condorOut/final/hMC_WJets_amc_GH.root"      ;
  inFileName[8] = "../condorOut/final/hMCtot_DY_GH.root"          ;
  inFileName[9] = "../condorOut/final/hMCtot_TOP_GH.root"         ;
  inFileName[10]= "../condorOut/final/hMCtot_VV_GH.root"          ;
  inFileName[11]= "../condorOut/final/hMC_Wprime_1800_GH.root"    ;
  inFileName[12]= "../condorOut/final/hMC_Wprime_3800_GH.root"    ;

TString xTitle[29]={"Muon P_{T} [GeV]", "Muon #eta", "Muon #phi", "#sigma_{Pt}/P_{T}(<0.3)","#Sigmap_{T}^{Trk}(#DeltaR<0.3)","#Sigmap_{T}^{Trk}(#DeltaR<0.3)/p_{T}(#mu)(<0.1)", "nPixelHits(>0)","nTrackerLayers(>5)","nMatchedStations(>1)","nValidMuonHits(>0)","dZ(<5mm)","dXY(<2mm)","MET(MuEGClean) [GeV]", "MET(MuEGClean) #phi", "sumEt(MuEGClean)", "MET(UnCorr) [GeV]", "MET (UnCorr) #phi","PUPPI MET [GeV]", "PUPPI MET #phi" ,"nPV_Test(w PUreweight)", "inv M_{#mu#mu}","nPV(w/o PUreWeight)","nPV(w PUreweight)" ,"M_{T} [GeV]","M_{T}(UnCorr) [GeV]","#Delta #phi","#Delta #phi(UnCorr)","E_{T}^{#mu}/E_{T}^{miss}","E_{T}^{#mu}/E_{T}^{miss}(UnCorr)"} ;

TString save[29] = { "mu_Pt", "mu_Eta","mu_Phi", "mu_momQuality","mu_TrkIsoR03_sumPt", "mu_relTrkIsoR03", "mu_nPixelHits", "mu_nTrkLayers", "mu_nMatchedStations","mu_nValidMuonHits", "mu_dz", "mu_dxy","met_Pt","met_Phi","met_sumEt", "met_unCorr","met_unCorrPhi", "met_Puppi","met_PuppiPhi","nPVwTest","invM","nPV","nPVw","wp_Mt","wp_MtUnCorr","wp_Dphi","wp_DphiUncorr","wp_EtRatio","wp_EtRatioUnCorr" } ;

const char * hname_ori[29] = { "mu_Pt_ori", "mu_Eta_ori","mu_Phi_ori", "mu_momQuality_ori","mu_TrkIsoR03_sumPt_ori", "mu_relTrkIsoR03_ori", "mu_nPixelHits_ori", "mu_nTrkLayers_ori", "mu_nMatchedStations_ori","mu_nValidMuonHits_ori", "mu_dz_ori", "mu_dxy_ori","met_Pt_ori","met_Phi_ori","met_sumEt_ori", "met_unCorr_ori","met_unCorrPhi_ori", "met_Puppi_ori","met_PuppiPhi_ori","nPVwTest_ori","invM_ori","nPV_ori","nPVw_ori","wp_Mt_ori","wp_MtUnCorr_ori","wp_Dphi_ori","wp_DphiUncorr_ori","wp_EtRatio_ori","wp_EtRatioUnCorr_ori" } ;


TH1F* hist_ori[NumFiles][29] ;
TH1F* hist_Data_ori[29] ;
TH1F* hist_BG_ori[NumBGs][29] ;
TH1F* hist_BGall_ori[5][29] ;//W WJet DY Top VV
TH1F* hist_AllBG_ori[29] ;
TH1F* hist_Sig_ori[NumSigs][29] ;
THStack* hs_BG_ori[29];
TH1F* hist_Ratio_ori[29] ;
TString save_ori[29] ;


 // TFile* outfile = new TFile("can_20170323.root","recreate");

  for(int i=0; i< NumFiles; i++){
        openFile[i] = TFile::Open(inFileName[i]);
        cout << "# OpenFile[" << i << "] : "<< inFileName[i] << endl;

        for(int var = 0; var < 29 ; var++){
            hist_ori[i][var] = (TH1F*)openFile[i]->Get(hname_ori[var]) ;
            //if(!((var>4 and var<8) || var == 16 || var ==17 || var ==18)) {hist_ori[i][var]->Rebin(5); }
            hs_BG_ori[var] = new THStack(hname_ori[var], hname_ori[var] ) ;

            if( i == 0 ) {
                cout << "## i = 0 Data " << var << " " << hname_ori[var] << endl ;
                hist_Data_ori[var] = hist_ori[i][var] ;
            }if( i > 0 && i < NumBGs+1){
                hist_BG_ori[i-1][var] = hist_ori[i][var] ;

                if(i==1) { hist_BG_ori[0][var]->Scale(1.04319 * 1.0) ;//W100 (NTotGen/NGenMassCut * NLO-kfactor)
                }if(i==2){ hist_BG_ori[1][var]->Scale(1.03528 * 1.1) ;//W200
                }if(i==3){ hist_BG_ori[2][var]->Scale(1.06388 * 1.3) ;//W500
                }if(i==4){ hist_BG_ori[3][var]->Scale(1.03258 * 1.1) ;//W1000
                }if(i==5){ hist_BG_ori[4][var]->Scale(1.07822 * 1.2) ;//W2000
                }if(i==6){ hist_BG_ori[5][var]->Scale(1.0 * 0.85)    ;//W3000
                }if(i==7){ hist_BG_ori[6][var]->Scale(1.01939 )      ;//WJetamc
                }
            }/*Backgrounds*/
            if( i > NumBGs ){
                int sigN = i - (NumBGs+1) ;
                hist_Sig_ori[sigN][var] = hist_ori[i][var] ;
                //if(!((var>4 and var<8) || var == 16 || var ==17 || var ==18)) {hist_Sig_ori[sigN][var]->Rebin(5); }
                hist_Sig_ori[sigN][var]->SetLineWidth(2) ; hist_Sig_ori[sigN][var]->SetLineColor(i-5) ;
            }/*Signals*/
       }/*For 29vars*/
   }/*For numFiles*/

   for(int j=0; j<29; j++){
       /* Merge W SM Background */
       hist_BGall_ori[0][j] = new TH1F(*hist_BG_ori[0][j]);
       hist_BGall_ori[0][j]->Add(hist_BG_ori[1][j]);
       hist_BGall_ori[0][j]->Add(hist_BG_ori[2][j]);
       hist_BGall_ori[0][j]->Add(hist_BG_ori[3][j]);
       hist_BGall_ori[0][j]->Add(hist_BG_ori[4][j]);
       hist_BGall_ori[0][j]->Add(hist_BG_ori[5][j]);
       hist_BGall_ori[0][j]->SetFillColor(kAzure+1);
       /* WJet Background */
       hist_BGall_ori[1][j] = new TH1F(*hist_BG_ori[6][j]);
       hist_BGall_ori[1][j]->SetFillColor(kAzure+2);
       /* DY Background */
       hist_BGall_ori[2][j] = new TH1F(*hist_BG_ori[7][j]);
       hist_BGall_ori[2][j]->SetFillColor(kRed);
       /* TOP Background */
       hist_BGall_ori[3][j] = new TH1F(*hist_BG_ori[8][j]);
       hist_BGall_ori[3][j]->SetFillColor(6);//pink
       /* VV Background */
       hist_BGall_ori[4][j] = new TH1F(*hist_BG_ori[9][j]);
       hist_BGall_ori[4][j]->SetFillColor(kGreen);

       hs_BG_ori[j]->Add(hist_BGall_ori[4][j]) ;//VV
       hs_BG_ori[j]->Add(hist_BGall_ori[2][j]) ;//DY
       hs_BG_ori[j]->Add(hist_BGall_ori[3][j]) ;//Top
       hs_BG_ori[j]->Add(hist_BGall_ori[1][j]) ;//WJet
       hs_BG_ori[j]->Add(hist_BGall_ori[0][j]) ;//W

       hist_AllBG_ori[j] = new TH1F(*hist_BGall_ori[0][j]);
       hist_AllBG_ori[j]->Add(hist_BGall_ori[1][j]);
       hist_AllBG_ori[j]->Add(hist_BGall_ori[2][j]);
       hist_AllBG_ori[j]->Add(hist_BGall_ori[3][j]);
       hist_AllBG_ori[j]->Add(hist_BGall_ori[4][j]);

       hist_Data_ori[j]->SetMarkerSize(0.9); hist_Data_ori[j]->SetMarkerStyle(20); hist_Data_ori[j]->SetMarkerColor(1); //hist_Data_ori[j]->SetLineColor(1) ;
       hist_Ratio_ori[j] = (TH1F*)hist_Data_ori[j]->Clone();
       hist_Ratio_ori[j]->Divide(hist_AllBG_ori[j]);

       hist_Ratio_ori[j]->SetMarkerSize(0.9); hist_Ratio_ori[j]->SetMarkerStyle(20); hist_Ratio_ori[j]->SetMarkerColor(1); //hist_Ratio_ori[j]->SetLineColor(1) ;
    }/*For 29vars*/

  TCanvas*   c_ori[25];
  TPad*     p1_ori[29];
  TPad*     p2_ori[29];
  TH2F*     n1_ori[29];
  TH2F*     n2_ori[29];
  TLegend* leg_ori[29];

  double xmin[29] ; 
  double xmax[29] ;
  double ymin_ori[29] ;
  double ymax_ori[29] ;
  double YMax_ori[29] ;

/* "0mu_Pt", "1mu_Eta","2mu_Phi", "3mu_momQuality","4mu_TrkIsoR03_sumPt", "5mu_relTrkIsoR03", "6mu_nPixelHits", "7mu_nTrkLayers", "8mu_nMatchedStations","9mu_nValidMuonHits", "10mu_dz", "11mu_dxy",
 * "12met_Pt","13met_Phi","14met_sumEt", "15met_unCorr","16met_unCorrPhi", "17met_Puppi","18met_PuppiPhi",
 * "19nPVwTest","20invM","21nPV","22nPVw",
 * "23wp_Mt","24wp_MtUnCorr","25wp_Dphi","26wp_DphiUncorr","27wp_EtRatio","28wp_EtRatioUnCorr" } ;*/


 for(int k=0; k<29; k++){
        bool isLogy = true ;
        if(k==19 || k==21 || k==22 ) isLogy = false; //nPV

        xmin[k] = hist_Data_ori[k]->GetXaxis()->GetXmin();
        xmax[k] = hist_Data_ori[k]->GetXaxis()->GetXmax();

       // if(k==0) xmax[k] = 3000. ; //Muon Pt distribution
       // if(k==19 || k==20 ) xmin[k] = 150. ; //Mt distribution
       // if(k==19 || k==20 ) xmax[k] = 3000. ; //Mt distribution
        ymin_ori[k] = (!isLogy) ? 0.0 : 10e-2;
        YMax_ori[k] = hs_BG_ori[k]->GetMaximum(); cout << YMax_ori[k] << endl;
        ymax_ori[k] = (!isLogy) ? hs_BG_ori[k]->GetMaximum() * 1.2 : 200 * YMax_ori[k];

        c_ori[k] = new TCanvas(Form("can%d",k), "canvas", 600, 600) ;
        p1_ori[k] = new TPad(Form("p1_%d",k), "", 0, 0.31, 1, 1.0);
        p1_ori[k]->SetBottomMargin(0.005)  ;
        p1_ori[k]->SetGrid()  ;
        p1_ori[k]->SetTicks(1,1)  ;
        if(isLogy) p1_ori[k]->SetLogy() ;
        p1_ori[k]->Draw()  ;
        p1_ori[k]->cd()  ;

        n1_ori[k] = new TH2F(Form("n1_%d",k), "", 2, xmin[k], xmax[k], 2, ymin_ori[k], ymax_ori[k]);
        n1_ori[k]->SetStats(kFALSE);
        n1_ori[k]->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (oricision 3)
        n1_ori[k]->GetYaxis()->SetLabelSize(20);
        n1_ori[k]->GetXaxis()->SetTitleSize(25);
        n1_ori[k]->GetXaxis()->SetTitleFont(45);
        n1_ori[k]->GetXaxis()->SetTitleOffset(3.5);
        n1_ori[k]->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (oricision 3)
        n1_ori[k]->GetXaxis()->SetLabelSize(0);
        n1_ori[k]->Draw();
        hs_BG_ori[k]->Draw("hist same");
        hist_Sig_ori[0][k]->Draw("same");
        hist_Sig_ori[1][k]->Draw("same");
        hist_Data_ori[k]->Draw("ep same");

        text = new TLatex();
        text->SetNDC();
        text->SetTextColor(1);
        text->SetTextSize(0.030);
        text->DrawLatex(0.11,0.84,"#color[1]{#scale[1.5]{CMS}}");
        text->DrawLatex(0.11,0.81,"#color[1]{#scale[1.0]{#it{Preliminary} } }");
        text->DrawLatex(0.25,0.82,"#color[1]{#scale[1.2]{#mu + #slash{E}_{T}^{miss} channel}}");
        if(isBF)text->DrawLatex(0.58,0.92,"#color[1]{#scale[1.2]{2016 RunB-F_reMiniaod  19.72 fb^{-1} (13 TeV)} }");
        if(isGH)text->DrawLatex(0.58,0.92,"#color[1]{#scale[1.2]{2016 RunG-H_reMiniaod  16.15 fb^{-1} (13 TeV)} }");
        if(isFull)text->DrawLatex(0.58,0.92,"#color[1]{#scale[1.2]{2016 Full Run 35.9 fb^{-1} (13 TeV)} }");

        leg_ori[k] = new TLegend(0.63, 0.50, 0.97, 0.90,"");
        leg_ori[k]->SetBorderSize(0);
        leg_ori[k]->SetTextSize(0.036);
        leg_ori[k]->SetFillStyle(0);
        leg_ori[k]->AddEntry(hist_Data_ori[k]     , "Data"  , "pl");
        leg_ori[k]->AddEntry(hist_BGall_ori[0][k] , "W#rightarrowl#nu" ,"f");
        leg_ori[k]->AddEntry(hist_BGall_ori[1][k] , "WJet#rightarrowl#nu" ,"f");
        leg_ori[k]->AddEntry(hist_BGall_ori[3][k] , "t#bar{t}"   ,"f");
        leg_ori[k]->AddEntry(hist_BGall_ori[2][k] , "Z/#gamma*#rightarrowll"  ,"f");
        leg_ori[k]->AddEntry(hist_BGall_ori[4][k] , "DiBoson"    ,"f");
        leg_ori[k]->AddEntry(hist_Sig_ori[0][k]   , "SSM W' M=1.8TeV" ,"l");
        leg_ori[k]->AddEntry(hist_Sig_ori[1][k]   , "SSM W' M=3.8TeV" ,"l");
        leg_ori[k]->Draw();

        c_ori[k]->cd();
        p2_ori[k] = new TPad(Form("p2_%d",k), "", 0, 0.020, 1, 0.3);
        p2_ori[k]->SetTopMargin(0.01);
        p2_ori[k]->SetBottomMargin(0.25) ; //230);
        p2_ori[k]->SetGrid(); // vertical grid
        p2_ori[k]->SetTicks(1,1); 
        p2_ori[k]->Draw();
        p2_ori[k]->cd();       // pad2 becomes the current pad

        n2_ori[k] = new TH2F(Form("n2_%d",k),"", 2, xmin[k], xmax[k], 2, 0.5, 1.5);
        n2_ori[k]->SetStats(kFALSE);
        n2_ori[k]->GetYaxis()->SetTitle("Data/MC ");
        n2_ori[k]->GetXaxis()->SetTitle(xTitle[k]);
        n2_ori[k]->GetYaxis()->SetNdivisions(505);
        n2_ori[k]->GetYaxis()->SetTitleSize(25);
        n2_ori[k]->GetYaxis()->SetTitleFont(45);
        n2_ori[k]->GetYaxis()->SetTitleOffset(1.);
        n2_ori[k]->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (oricision 3)
        n2_ori[k]->GetYaxis()->SetLabelSize(20);
        n2_ori[k]->GetXaxis()->SetTitleSize(25);
        n2_ori[k]->GetXaxis()->SetTitleFont(45);
        n2_ori[k]->GetXaxis()->SetTitleOffset(3.);
        n2_ori[k]->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (oricision 3)
        n2_ori[k]->GetXaxis()->SetLabelSize(20);
        n2_ori[k]->Draw();
        hist_Ratio_ori[k]->Draw("ep same");

        TLine* line = new TLine(xmin[k], 1.0 , xmax[k], 1.0);
        line->SetLineStyle(1);  line->SetLineWidth(1);  line->SetLineColor(kBlue);
        line->Draw();

        //c_ori[k]->Write();
        
        if(isBF)save_ori[k] = save[k]+"_ori_0322_BF.png" ;
        if(isGH)save_ori[k] = save[k]+"_ori_0322_GH.png" ;
        if(isFull)save_ori[k] = save[k]+"_ori_0322_Full.png" ;
        c_ori[k]->Print(save_ori[k]);
   }
   //outfile->Close();


//  return 0;

}

