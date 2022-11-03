#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TPad.h>
#include "./../ProdNpKNU/src/NpKNU.hh"


void DrawEff_new(){
   gROOT->Reset();
   gROOT->SetStyle("Default");
   gStyle->SetOptStat(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetFrameBorderMode(0);

   using namespace std ;
   using namespace TMath ;

   TString EBdataname = "./EB/hEff_EB_Data_RunFull_v3.root"; //"./ReminiaodEB/hEff_EB_Data_Full_reminiaod.root";   //
   TString EEdataname = "./EE/hEff_EE_Data_RunFull_v3.root"; //"./ReminiaodEE/hEff_EE_Data_Full_reminiaod.root";   //
   TString EEmcname =   "./EE/hEff_EE_MC_Full_v3.root" ; //hEff_EE_MC_WJet_Full.root";
   TString EBmcname =   "./EB/hEff_EB_MC_Full_v3.root" ;// hEff_EB_MC_WJet_Full_v3.root";
   TFile* EBmcfile = TFile::Open(EBmcname);
   TFile* EBdafile = TFile::Open(EBdataname);
   TFile* EEmcfile = TFile::Open(EEmcname);
   TFile* EEdafile = TFile::Open(EEdataname);
   //delete gROOT->GetListOfFiles()->FindObject(rootname); // clear memory of file name
   //file->GetListOfKeys()->Print();
   //file->ls();
   cout << "InputFile Okay " << endl;

//   const char * eb_total[3] = {"eb_pt_tot", "eb_eta_tot",   "eb_nPV_tot", /* "eb_sig_tot", "eb_hoe_tot"  */ } ;
//   const char * eb_pass1[3] = {"eb_pt_hlt0", "eb_eta_hlt0", "eb_nPV_hlt0",/* "eb_sig_hlt0", "eb_hoe_hlt0"*/ } ;
//   const char * eb_pass2[3] = {"eb_pt_hlt2", "eb_eta_hlt2", "eb_nPV_hlt2",/* "eb_sig_hlt2", "eb_hoe_hlt2"*/ } ;
//   const char * ee_total[3] = {"ee_pt_tot", "ee_eta_tot",    "ee_nPV_tot"} ;
//   const char * ee_pass1[3] = {"ee_pt_hlt0", "ee_eta_hlt0", "ee_nPV_hlt0"} ;
//   const char * ee_pass2[3] = {"ee_pt_hlt2", "ee_eta_hlt2", "ee_nPV_hlt2"} ;
   Double_t  elebin[3]  =  { 20,    16, 100 /*, 50,   50*/};
   Double_t  elexmin[3]  = {90.0, -2.5, 0.0 /*, 0.0, 0.0*/};
   Double_t  elexmax[3]  = {2000,  2.5, 100 /*, 0.03,0.1*/};
   TString  elexTitle[3] = {"Electron p_{T} [GeV]", "Electron #eta_{SC}",  "nPV" /*, "#sigma_{i#etai#eta}", "H/E"*/ };

  TH1F* hdata_ee_hlt0 ;   TH1F*   hmc_ee_hlt0 ;
  TH1F* hdata_ee_hlt2 ;   TH1F*   hmc_ee_hlt2 ;
  TH1F* hdata_eb_tot   ;   TH1F*   hmc_eb_tot  ;
  TH1F* hdata_eb_hlt0 ;   TH1F*   hmc_eb_hlt0 ;
  TH1F* hdata_eb_hlt2 ;   TH1F*   hmc_eb_hlt2 ;

  TH1F*   hdata_ee_eff0 ;   TH1F*   hmc_ee_eff0 ;
  TH1F*   hdata_ee_eff1 ;   TH1F*   hmc_ee_eff1 ;
  TH1F*   hdata_eb_eff0 ;   TH1F*   hmc_eb_eff0 ;
  TH1F*   hdata_eb_eff1 ;   TH1F*   hmc_eb_eff1 ;

  TH1F*   hratio_eb0;        TH1F*   hratio_ee0   ;
  TH1F*   hratio_eb1;        TH1F*   hratio_ee1   ;

   hdata_ee_tot  = (TH1F*)EEdafile->Get("ee_pt_tot");
   hdata_ee_hlt0 = (TH1F*)EEdafile->Get("ee_pt_hlt0");
   hdata_ee_hlt2 = (TH1F*)EEdafile->Get("ee_pt_hlt2");
   hdata_eb_tot  = (TH1F*)EBdafile->Get("eb_pt_tot");
   hdata_eb_hlt0 = (TH1F*)EBdafile->Get("eb_pt_hlt0");
   hdata_eb_hlt2 = (TH1F*)EBdafile->Get("eb_pt_hlt2");

    hmc_ee_tot = (TH1F*)EEmcfile->Get("ee_pt_tot");
   hmc_ee_hlt0 = (TH1F*)EEmcfile->Get("ee_pt_hlt0");
   hmc_ee_hlt2 = (TH1F*)EEmcfile->Get("ee_pt_hlt2");
    hmc_eb_tot = (TH1F*)EBmcfile->Get("eb_pt_tot");
   hmc_eb_hlt0 = (TH1F*)EBmcfile->Get("eb_pt_hlt0");
   hmc_eb_hlt2 = (TH1F*)EBmcfile->Get("eb_pt_hlt2");

     hdata_ee_eff0 = (TH1F*) hdata_ee_hlt0->Clone();
   hdata_ee_eff1 = (TH1F*) hdata_ee_hlt2->Clone();
   hdata_eb_eff0 = (TH1F*) hdata_eb_hlt0->Clone();
   hdata_eb_eff1 = (TH1F*) hdata_eb_hlt2->Clone();
   hmc_ee_eff0 = (TH1F*) hmc_ee_hlt0->Clone();
   hmc_ee_eff1 = (TH1F*) hmc_ee_hlt2->Clone();
   hmc_eb_eff0 = (TH1F*) hmc_eb_hlt0->Clone();
   hmc_eb_eff1 = (TH1F*) hmc_eb_hlt2->Clone();


   hdata_ee_eff0 ->Divide(hdata_ee_eff0/* hdata_ee_hlt0*/, hdata_ee_tot , 1, 1, "B");
   hdata_ee_eff1 ->Divide(hdata_ee_eff1/* hdata_ee_hlt2*/, hdata_ee_tot , 1, 1, "B");
   hmc_ee_eff0   ->Divide(hmc_ee_eff0  /* hmc_ee_hlt0  */, hmc_ee_tot , 1, 1, "B");
   hmc_ee_eff1   ->Divide(hmc_ee_eff1  /* hmc_ee_hlt2  */, hmc_ee_tot , 1, 1, "B");

   hdata_eb_eff0 ->Divide(hdata_eb_eff0/* hdata_eb_hlt0*/, hdata_eb_tot , 1, 1, "B");
   hdata_eb_eff1 ->Divide(hdata_eb_eff1/* hdata_eb_hlt2*/, hdata_eb_tot , 1, 1, "B");
   hmc_eb_eff0   ->Divide(hmc_eb_eff0  /* hmc_eb_hlt0  */, hmc_eb_tot , 1, 1, "B");
   hmc_eb_eff1   ->Divide(hmc_eb_eff1  /* hmc_eb_hlt2  */, hmc_eb_tot , 1, 1, "B");

      hratio_eb0 =  (TH1F*) hdata_eb_eff0->Clone()  ;
      hratio_eb1 =  (TH1F*) hdata_eb_eff1->Clone()  ;
      hratio_ee0 =  (TH1F*) hdata_ee_eff0->Clone()  ;
      hratio_ee1 =  (TH1F*) hdata_ee_eff1->Clone()  ;


   hratio_eb0->Divide(hratio_eb0, hmc_eb_eff0, 1, 1, "B");
   hratio_eb1->Divide(hratio_eb1, hmc_eb_eff1, 1, 1, "B");
   hratio_ee0->Divide(hratio_ee0, hmc_ee_eff0, 1, 1, "B");
   hratio_ee1->Divide(hratio_ee1, hmc_ee_eff1, 1, 1, "B");

//      cout << "---------------------------" << endl;
//        cout << "EE Data Total N = " << hdata_ee_tot->GetEntries() << "  Pass0N " << hdata_ee_hlt0->GetEntries()<< "  Pass1N " << hdata_ee_hlt2->GetEntries() << endl;
//        cout << "EE Data Eff     = " << hdata_ee_hlt0->GetEntries()/hdata_ee_tot->GetEntries() << "+/-" <<data_eeErr  << " " << hdata_ee_hlt2->GetEntries()/hdata_ee_tot->GetEntries() << endl;
//        cout << "EE MC   Total N = " << hmc_ee_tot->GetEntries()   << "  Pass0N " << hmc_ee_hlt0->GetEntries()  << "  Pass1N " << hmc_ee_hlt2->GetEntries()  << endl;
//        cout << "EE MC Eff     = " << hmc_ee_hlt0->GetEntries()/hmc_ee_tot->GetEntries() << "+/-" <<mc_eeErr  << "  " << hmc_ee_hlt2->GetEntries()/hmc_ee_tot->GetEntries() << endl;
//        cout << "---------------------------" << endl;
  
  
     int eebins = hdata_ee_eff0->GetNbinsX();
     int ebbins = hdata_eb_eff0->GetNbinsX();
        /////////////////////////////////
        //////////EE Canvas//////////////
        /////////////////////////////////
       TCanvas* cane = new TCanvas("cane", "canvas", 800, 800) ;
       TPad* p1e = new TPad("p1e", "", 0, 0.3, 1, 1.0);
       p1e->SetBottomMargin(0.0)  ;
       p1e->SetGrid()  ;
       p1e->Draw()  ;
       p1e->cd()  ;
       TH2F* n1e = new TH2F("n1e", "", 2, 98, 1000, 2,0.0, 1.20);
       n1e->SetStats(kFALSE);
       n1e->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
       n1e->GetYaxis()->SetLabelSize(20);
       n1e->GetYaxis()->SetTitleSize(25);
       n1e->GetXaxis()->SetTitleFont(45);
       n1e->GetYaxis()->SetTitle("HLT efficiency for HEEPV7.0 e");
       n1e->GetYaxis()->SetTitleFont(45);
       n1e->GetXaxis()->SetTitleOffset(3.5);
       n1e->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
       n1e->GetXaxis()->SetLabelSize(20);
       n1e->Draw();

       hmc_ee_eff0->SetMarkerColor(kSpring-1);
       hmc_ee_eff0->SetMarkerStyle(1);
       hmc_ee_eff0->SetLineColor(kSpring-1);
       hmc_ee_eff0->SetFillColor(kSpring-1);
     hdata_ee_eff0->SetMarkerColor(kBlue);
     hdata_ee_eff0->SetMarkerStyle(22);
     hdata_ee_eff0->SetLineColor(kBlue);

       hmc_ee_eff1->SetMarkerColor(kSpring-1);
       hmc_ee_eff1->SetMarkerStyle(1);
       hmc_ee_eff1->SetLineColor(kSpring-1);
       hmc_ee_eff1->SetFillColor(kSpring-1);
     hdata_ee_eff1->SetMarkerColor(kBlue);
     hdata_ee_eff1->SetMarkerStyle(22);
     hdata_ee_eff1->SetLineColor(kBlue);
/////////////FixMe
       hmc_ee_eff0->Draw("e2same");
     hdata_ee_eff0->Draw("same");
       TLatex* text = new TLatex();
       text->SetNDC();
       text->SetTextColor(1);
       text->SetTextSize(0.030);
       text->DrawLatex(0.11,0.86,"#color[1]{#scale[1.8]{CMS}}");
       text->DrawLatex(0.11,0.80,"#color[1]{#scale[1.5]{#it{Preliminary} } }");
       text->DrawLatex(0.55,0.92,"#color[1]{#scale[1.3]{2016 full Run 35.9 fb^{-1} (13 TeV)} }");
       //text->DrawLatex(0.35,0.25,"#color[1.1]{#scale[1.7]{HLT_Ele115 || ECALHT800 }}");
       text->DrawLatex(0.35,0.25,"#color[1.1]{#scale[1.7]{HLT_Ele115_CaloIdVT_GsfTrkIdT }}");

       text->DrawLatex(0.35,0.20,"#color[1.1]{#scale[1.7]{ReReco B-G PromptReco H}}");
       //text->DrawLatex(0.35,0.20,"#color[1.1]{#scale[1.7]{03Feb2017 ReMiniAOD}}");

       TLegend* lege = new TLegend(0.7, 0.05, 0.94, 0.20,"");
       lege->SetBorderSize(0);
       lege->SetTextSize(0.05);
       lege->SetFillStyle(0);
       lege->AddEntry(hmc_ee_eff1     ,"EE MC "  , "f");
       lege->AddEntry(hdata_ee_eff1   ,"EE Data  "  , "pl");
       lege->Draw();
       cane->cd();
       TPad* p2e = new TPad("p2e", "", 0, 0.022, 1,0.3);
       p2e->SetTopMargin(0.0);
       p2e->SetBottomMargin(0.230);
       p2e->SetGrid(); // vertical grid
       p2e->Draw();
       p2e->cd();       // pad2 becomes the current pad

       TH2F* n2e = new TH2F("n2e","", 2, 98, 1000, 2, 0.8, 1.2);
       n2e->SetStats(kFALSE);
       n2e->SetStats(kFALSE);
       n2e->GetYaxis()->SetTitle("Data/MC");
       n2e->GetXaxis()->SetTitle("Electron P_{T} [GeV]");
       n2e->GetYaxis()->SetNdivisions(505);
       n2e->GetYaxis()->SetTitleSize(23);
       n2e->GetYaxis()->SetTitleFont(45);
       n2e->GetYaxis()->SetTitleOffset(1.);
       n2e->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
       n2e->GetYaxis()->SetLabelSize(20);
       n2e->GetXaxis()->SetTitleSize(25);
       n2e->GetXaxis()->SetTitleFont(45);
       n2e->GetXaxis()->SetTitleOffset(3.);
       n2e->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
       n2e->GetXaxis()->SetLabelSize(20);
       n2e->Draw();
       hratio_ee0->SetMarkerStyle(20);
       hratio_ee0->SetMarkerColor(kBlue);
       hratio_ee0->SetLineColor(kBlue);
       hratio_ee1->SetMarkerStyle(20);
       hratio_ee1->SetMarkerColor(kBlue);
       hratio_ee1->SetLineColor(kBlue);
/////////////FixMe
       hratio_ee0->Draw("epsame");
//       hratio_ee1->Draw("epsame");

       TLine* line = new TLine(98, 1.0 , 1000, 1.0);
       line->SetLineStyle(2);  line->SetLineWidth(2);  line->SetLineColor(kGreen);
       line->Draw();
       cane->Print("heffvspt_ee_rereco_Type1.png");
