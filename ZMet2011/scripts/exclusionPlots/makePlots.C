#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"

#include "fedorContours.C"
#include "Pieter_observed/LeftSlepton_Combo_25.C"
#include "Pieter_observed/LeftSlepton_Combo_50.C"
#include "Pieter_observed/LeftSlepton_Combo_75.C"
#include "Pieter_observed/TauEnriched_Combo_25.C"
#include "Pieter_observed/TauEnriched_Combo_50.C"
#include "Pieter_observed/TauEnriched_Combo_75.C"
#include "Pieter_expected/LeftSlepton_Combo_25_expected.C"
#include "Pieter_expected/LeftSlepton_Combo_50_expected.C"
#include "Pieter_expected/LeftSlepton_Combo_75_expected.C"
#include "Pieter_expected/TauEnriched_Combo_25_expected.C"
#include "Pieter_expected/TauEnriched_Combo_50_expected.C"
#include "Pieter_expected/TauEnriched_Combo_75_expected.C"

using namespace std;

bool isPreliminary = false;
char* isPrelimChar = (char*) "";

int width1 = 8;
int width2 = 4;

// int width1 = 4;
// int width2 = 2;

void makeFloridaPlot(char* sample, int x , bool printplot);

void removeDiagonal( TH2D* h , float deltaM ){

  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

      float mg = h->GetXaxis()->GetBinCenter(ibin);
      float ml = h->GetYaxis()->GetBinCenter(jbin);

      cout << "mg ml " << mg << " " << ml << endl;

      if( mg - ml < deltaM ) h->SetBinContent(ibin,jbin,0);
      
    }
  }
}

TGraph* extendGraph(TGraph* gin){
  
  const unsigned int nold = gin->GetN();
  //cout << "npoints " << nold << endl;

  const unsigned int nnew = nold + 2;
  
  float xnew[nnew];
  float ynew[nnew];

  float xold[nold];
  float yold[nold];

  Double_t x;
  Double_t y;

  //cout << endl << "old" << endl;
  for(int i = 0 ; i < nold ; ++i){
    gin->GetPoint(i,x,y);
    //cout << i << " " << x << " " << y << endl;
    xold[i] = x;
    yold[i] = y;

    xnew[i+1] = x;
    ynew[i+1] = y;
  }

  xnew[0]      = xold[0];      ynew[0] = 0.0;
  xnew[nnew-1] = xold[nold-1]; ynew[0] = 0.0;

  // cout << endl << "new" << endl;
  // for(int i = 0 ; i < nnew ; ++i){
  //   cout << i << " " << xnew[i] << " " << ynew[i] << endl;
  // }

  TGraph* gout = new TGraph(nnew,xnew,ynew);

  gout->SetLineColor(gin->GetLineColor());
  gout->SetLineWidth(gin->GetLineWidth());
  gout->SetLineStyle(gin->GetLineStyle());

  return gout;
}

TGraph* fixupGraph(TGraph* gin){
  
  const unsigned int nold = gin->GetN();
  
  float xgr[nold];
  float ygr[nold];

  Double_t x;
  Double_t y;

  //cout << endl << "old" << endl;
  for(int i = 0 ; i < nold ; ++i){
  
    gin->GetPoint(i,x,y);
    cout << i << " " << x << " " << y << endl;

    if( fabs(x-250)<0.1 && fabs(y-25)<0.1 ){
      y = -10.0;
      cout << "SWITCHING Y FROM 25 TO -10" << endl;
    }

    xgr[i] = x;
    ygr[i] = y;

  }

  TGraph* gout = new TGraph(nold,xgr,ygr);

  gout->SetLineColor(gin->GetLineColor());
  gout->SetLineWidth(gin->GetLineWidth());
  gout->SetLineStyle(gin->GetLineStyle());

  return gout;
}

void formatHist( TH2D* hist ){

  //hist->GetXaxis()->SetTitle("m(#tilde{#chi}_{2}^{0}) = m(#tilde{#chi}_{1}^{#pm}) [GeV]");
  //hist->GetYaxis()->SetTitle("m(#tilde{#chi}_{1}^{0}) [GeV]");
  hist->GetXaxis()->SetTitle("m_{#tilde{#chi}_{2}^{0}} = m_{#tilde{#chi}_{1}^{#pm}} [GeV]");
  hist->GetYaxis()->SetTitle("m_{#tilde{#chi}_{1}^{0}} [GeV]");
  hist->GetZaxis()->SetTitle("95% CL UL #sigma#timesBF [fb]");
  hist->GetXaxis()->SetTitleOffset(1.12);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->SetMinimum(1);
  hist->SetMaximum(1000);
  
}

// void cmsPrelim( double intLumi )
// {

//   TLatex latex;
//   latex.SetNDC();
//   latex.SetTextSize(0.04);
//   latex.DrawLatex(0.18,0.93,"CMS Preliminary,  #sqrt{s}=7 TeV,  L_{int}=4.98 fb^{-1}");
// }

void cmsPrelim(double intLumi, bool prelim)
{
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62);
        if(prelim) latex.SetTextSize(0.035);
        else       latex.SetTextSize(0.045);

        latex.SetTextAlign(11); // align left
        if(prelim) latex.DrawLatex(0.17,0.92,"CMS Preliminary");
        else       latex.DrawLatex(0.17,0.92,"CMS");

        latex.SetTextAlign(31); // align right
        latex.DrawLatex(0.85, 0.92, Form("#sqrt{s} = 7 TeV, L_{int} = %4.2f fb^{-1}", intLumi));
}

void makePlots( bool printPlots = false){

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.05);


  if( isPreliminary ) isPrelimChar = (char*) "_prelim";

  //-------------------------------------------
  // Rutgers/KIT
  //-------------------------------------------

  TGraph* mod2i_observed   = model2i_observed();
  TGraph* mod2i_expected   = model2i_expected();
  TGraph* mod2i_expectedP1 = model2i_expectedP1();
  TGraph* mod2i_expectedM1 = model2i_expectedM1();
  TGraph* mod2i_observedP  = model2i_observedp();
  TGraph* mod2i_observedM  = model2i_observedm();

  mod2i_observed->SetLineWidth(width1);

  mod2i_expected->SetLineWidth(width1);
  mod2i_expected->SetLineStyle(2);

  mod2i_expectedP1->SetLineColor(2);
  mod2i_expectedP1->SetLineWidth(width2);
  mod2i_expectedP1->SetLineStyle(3);

  mod2i_expectedM1->SetLineColor(2);
  mod2i_expectedM1->SetLineWidth(width2);
  mod2i_expectedM1->SetLineStyle(3);

  mod2i_observedP->SetLineColor(4);
  mod2i_observedP->SetLineWidth(width2);
  mod2i_observedP->SetLineStyle(4);

  mod2i_observedM->SetLineColor(4);
  mod2i_observedM->SetLineWidth(width2);
  mod2i_observedM->SetLineStyle(4);

  //-----------------
  // model 2i
  //-----------------

  TFile *fkit_2i  = TFile::Open("KIT_2i.root");

  TH2D*   h2i     = (TH2D*)   fkit_2i->Get("xSecObserved");
  TH2D*   h2i1    = (TH2D*)   fkit_2i->Get("xSecObserved1");
  TH2D*   h2i0    = (TH2D*)   fkit_2i->Get("xSecObserved0");
  TGraph* gr2i_1  = (TGraph*) fkit_2i->Get("Graph1");
  TGraph* gr2i_2  = (TGraph*) fkit_2i->Get("Graph2");
  // TGraph* gr2i_3  = (TGraph*) fkit_2i->Get("Graph3");
  // TGraph* gr2i_4  = (TGraph*) fkit_2i->Get("Graph4");
  
  TCanvas *can_2i = new TCanvas("can_2i","can_2i",600,600);
  can_2i->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(h2i0);
  //removeDiagonal(h2i,50);
  //removeDiagonal(h2i1,50);
  h2i0->Draw("colz");
  h2i->Draw("colsame");
  h2i1->Draw("colsame");
  //gr2i_1->Draw("l");
  //gr2i_2->Draw("l");

  mod2i_observed->Draw("l");
  mod2i_expected->Draw("l");
  mod2i_expectedP1->Draw("l");
  mod2i_expectedM1->Draw("l");
  mod2i_observedP->Draw("l");
  mod2i_observedM->Draw("l");

  h2i0->Draw("axissame");
  //gr2i_3->Draw("l");
  //gr2i_4->Draw("l");
  cmsPrelim(4.98,isPreliminary);

  TLegend *leg = new TLegend(0.2,0.72,0.6,0.88);
  //leg->AddEntry(gr2i_1,"observed","l");
  //leg->AddEntry(gr2i_2,"median expected","l");
  //leg->AddEntry(gr2i_3,"expected #pm1#sigma","l");
  leg->AddEntry(mod2i_observed  ,"Observed","l");
  leg->AddEntry(mod2i_observedP ,"Observed (#pm1#sigma^{theory})","l");
  leg->AddEntry(mod2i_expected  ,"Median expected","l");
  leg->AddEntry(mod2i_expectedP1,"Expected (#pm1#sigma)","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();


  tex->SetTextSize(0.04);
  tex->DrawLatex(0.19,0.65,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  //tex->DrawLatex(0.18,0.65,"m(#tilde{#font[12]{l}}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  //tex->DrawLatex(0.19,0.65,"m_{#tilde{#font[12]{l}}} = 0.5(^{}m_{#tilde{#chi}_{2}^{0}} = m_{#tilde{#chi}_{1}^{#pm}}) + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");
  //tex->DrawLatex(0.19,0.61,"x_{#tilde{#font[12]{l}}} = 0.5");
  tex->DrawLatex(0.19,0.59,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{#font[12]{l}}#font[12]{l} (BF=0.5)");
  tex->DrawLatex(0.19,0.53,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#font[12]{l}}#nu_{#font[12]{l}} , #font[12]{l}#tilde{#nu}_{#font[12]{l}}");

  //tex->DrawLatex(0.01,0.03,"m_{#tilde{#font[12]{l}}} = 0.5(^{}m_{#tilde{#chi}_{2}^{0}}=m_{#tilde{#chi}_{1}^{#pm}}) + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");
  tex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.5^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");

  can_2i->Modified();
  can_2i->Update();
  if( printPlots) can_2i->Print(Form("multilepton_flavordemocratic_Fig7%s.pdf",isPrelimChar));
  if( printPlots){
    can_2i->Print("Figure7a.pdf");
    can_2i->Print("Figure7a.png");
  }

  //-----------------
  // model 2a
  //-----------------

  TGraph* mod2a_observed   = model2a_observed();
  TGraph* mod2a_expected   = model2a_expected();
  TGraph* mod2a_expectedP1 = model2a_expectedP1();
  TGraph* mod2a_expectedM1 = model2a_expectedM1();
  TGraph* mod2a_observedP  = model2a_observedp();
  TGraph* mod2a_observedM  = model2a_observedm();

  mod2a_observed->SetLineWidth(width1);

  mod2a_expected->SetLineWidth(width1);
  mod2a_expected->SetLineStyle(2);

  mod2a_expectedP1->SetLineColor(2);
  mod2a_expectedP1->SetLineWidth(width2);
  mod2a_expectedP1->SetLineStyle(3);

  mod2a_expectedM1->SetLineColor(2);
  mod2a_expectedM1->SetLineWidth(width2);
  mod2a_expectedM1->SetLineStyle(3);

  mod2a_observedP->SetLineColor(4);
  mod2a_observedP->SetLineWidth(width2);
  mod2a_observedP->SetLineStyle(4);

  mod2a_observedM->SetLineColor(4);
  mod2a_observedM->SetLineWidth(width2);
  mod2a_observedM->SetLineStyle(4);

  TFile *fkit_2a  = TFile::Open("KIT_2a.root");

  TH2D*   h2a     = (TH2D*)   fkit_2a->Get("xSecObserved");
  TH2D*   h2a0    = (TH2D*)   fkit_2a->Get("xSecObserved0");
  TGraph* gr2a_1  = (TGraph*) fkit_2a->Get("Graph1");
  TGraph* gr2a_2  = (TGraph*) fkit_2a->Get("Graph2");
  
  TCanvas *can_2a = new TCanvas();
  can_2a->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(h2a0);
  h2a0->Draw("colz");
  //removeDiagonal(h2a,50);
  h2a->Draw("colsame");
  //gr2a_1->Draw("l");
  //gr2a_2->Draw("l");

  mod2a_observed->Draw("l");
  mod2a_expected->Draw("l");
  mod2a_expectedP1->Draw("l");
  mod2a_expectedM1->Draw("l");
  mod2a_observedP->Draw("l");
  mod2a_observedM->Draw("l");

  h2a0->Draw("axissame");
  cmsPrelim(4.98,isPreliminary);

  leg->Draw();

  // tex->SetTextSize(0.03);
  // tex->DrawLatex(0.18,0.70,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  // tex->DrawLatex(0.18,0.65,"m(#tilde{l}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  // tex->DrawLatex(0.18,0.60,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l (BF=100%)");
  // tex->DrawLatex(0.18,0.55,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#tau}#nu_{#tau}");

  tex->SetTextSize(0.04);
  tex->DrawLatex(0.19,0.65,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  //tex->DrawLatex(0.18,0.65,"m(#tilde{#font[12]{l}}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  //tex->DrawLatex(0.19,0.61,"m_{#tilde{#font[12]{l}}} = 0.5(^{}m_{#tilde{#chi}_{2}^{0}} = m_{#tilde{#chi}_{1}^{#pm}}) + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");
  //tex->DrawLatex(0.19,0.61,"x_{#tilde{#font[12]{l}}} = 0.5");
  tex->DrawLatex(0.19,0.59,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{#font[12]{l}}#font[12]{l} (BF=1)");
  tex->DrawLatex(0.19,0.53,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#tau}#nu_{#tau}");

  //tex->DrawLatex(0.01,0.03,"m_{#tilde{#font[12]{l}}} = 0.5(^{}m_{#tilde{#chi}_{2}^{0}}=m_{#tilde{#chi}_{1}^{#pm}}) + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");
  tex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.5^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");

  can_2a->Modified();
  can_2a->Update();
  if( printPlots) can_2a->Print(Form("multilepton_tauenriched_Fig8%s.pdf",isPrelimChar));
  if( printPlots){
    can_2a->Print("Figure7b.pdf");
    can_2a->Print("Figure7b.png");
  }

  //-----------------
  // TChiWZ
  //-----------------

  TFile *fwz  = TFile::Open("combinePlots_VZ_Trilepton_AllUncertainties.root");

  TH2F*   hwz             = (TH2F*)   fwz->Get("hexcl");
  TGraph* grwz_combo      = (TGraph*) fwz->Get("gr_combo");
  TGraph* grwz_combo_exp  = (TGraph*) fwz->Get("gr_combo_exp");
  TGraph* grwz_tri        = (TGraph*) fwz->Get("gr_tri");
  TGraph* grwz_vzmet      = (TGraph*) fwz->Get("gr_vzmet");
  TGraph* grwz_combo_expp1  = (TGraph*) fwz->Get("gr_combo_expp1");
  TGraph* grwz_combo_expm1  = (TGraph*) fwz->Get("gr_combo_expm1");
  TGraph* grwz_comboUp      = (TGraph*) fwz->Get("gr_comboTheoryUp");
  TGraph* grwz_comboDn      = (TGraph*) fwz->Get("gr_comboTheoryDown");

  grwz_combo->SetLineWidth(width1);
  grwz_combo_exp->SetLineWidth(width1);

  grwz_tri->SetLineWidth(width2);
  grwz_tri->SetLineColor(6);
  grwz_tri->SetLineStyle(7);

  grwz_vzmet->SetLineWidth(width2);
  grwz_vzmet->SetLineColor(kViolet-5);
  grwz_vzmet->SetLineStyle(9);

  grwz_combo_expp1->SetLineColor(2);
  grwz_combo_expp1->SetLineWidth(width2);
  grwz_combo_expp1->SetLineStyle(3);

  grwz_combo_expm1->SetLineColor(2);
  grwz_combo_expm1->SetLineWidth(width2);
  grwz_combo_expm1->SetLineStyle(3);

  grwz_comboUp->SetLineColor(4);
  grwz_comboUp->SetLineWidth(width2);
  grwz_comboUp->SetLineStyle(4);

  grwz_comboDn->SetLineColor(4);
  grwz_comboDn->SetLineWidth(width2);
  grwz_comboDn->SetLineStyle(4);

  TCanvas *can_wz = new TCanvas();
  can_wz->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist((TH2D*)hwz);
  hwz->GetXaxis()->SetRangeUser(87.5,512.5);
  hwz->GetYaxis()->SetRangeUser(-12.5,512.5);
  //hwz->GetXaxis()->SetRangeUser(100,300);
  //hwz->GetYaxis()->SetRangeUser(  0,300);
  hwz->SetMinimum(50);
  hwz->SetMaximum(10000);
  hwz->Draw("colz");
  grwz_combo->Draw("l");
  grwz_combo_exp->Draw("l");
  grwz_combo_expp1->Draw("l");
  grwz_combo_expm1->Draw("l");
  grwz_comboUp->Draw("l");
  grwz_comboDn->Draw("l");
  grwz_tri->Draw("l");
  grwz_vzmet->Draw("l");
  cmsPrelim(4.98,isPreliminary);

  TLegend *legwz = new TLegend(0.2,0.6,0.65,0.88);
  legwz->AddEntry(grwz_combo       ,"Combined observed","l");
  legwz->AddEntry(grwz_comboUp     ,"Combined observed (#pm1#sigma^{theory})","l");
  legwz->AddEntry(grwz_combo_exp   ,"Combined median expected","l");
  legwz->AddEntry(grwz_combo_expp1 ,"Combined expected (#pm1#sigma)","l");
  legwz->AddEntry(grwz_vzmet       ,"2#font[12]{l}2j observed","l");
  legwz->AddEntry(grwz_tri         ,"Trilepton (M_{T}) observed","l");
  legwz->SetBorderSize(0);
  legwz->SetFillColor(0);
  legwz->Draw();
  legwz->Draw();
  
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.18,0.55,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm} #rightarrow WZ+E_{T}^{miss}");

  can_wz->Modified();
  can_wz->Update();
  if( printPlots){
    can_wz->Print(Form("WZ_Fig11%s.pdf",isPrelimChar));
    can_wz->Print(Form("WZ_Fig11%s.png",isPrelimChar));
    can_wz->Print("Figure10.pdf");
    can_wz->Print("Figure10.png");
  }

  hwz->GetXaxis()->SetRangeUser(100,300);
  hwz->GetYaxis()->SetRangeUser(  0,300);

  can_wz->Modified();
  can_wz->Update();
  if( printPlots){
    can_wz->Print(Form("WZ_zoom_Fig11%s.pdf",isPrelimChar));
    can_wz->Print(Form("WZ_zoom_Fig11%s.png",isPrelimChar));
    can_wz->Print("Figure14_Add.pdf");
    can_wz->Print("Figure14_Add.png");
  }

  //-----------------------------
  // Florida/ETH plots
  //-----------------------------

  makeFloridaPlot("LeftSlepton",25,printPlots);
  makeFloridaPlot("LeftSlepton",50,printPlots);
  makeFloridaPlot("LeftSlepton",75,printPlots);

  makeFloridaPlot("TauEnriched",25,printPlots);
  makeFloridaPlot("TauEnriched",50,printPlots);
  makeFloridaPlot("TauEnriched",75,printPlots);

}


TH2D* cloneHist( TH2D* hin ){

  TH2D* hout = new TH2D(Form("%s_clone",hin->GetName()),Form("%s_clone",hin->GetName()),
			hin->GetXaxis()->GetNbins(),hin->GetXaxis()->GetXmin(),hin->GetXaxis()->GetXmax(),
			hin->GetYaxis()->GetNbins(),hin->GetYaxis()->GetXmin(),hin->GetYaxis()->GetXmax());

  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ++ibin ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; ++jbin ){
      hout->SetBinContent(ibin,jbin,hin->GetBinContent(ibin,jbin));
    }
  }

  return hout;

}

void makeFloridaPlot(char* sample, int x, bool printPlots ){


  bool plotss = true;
  if( TString(sample).Contains("Tau") && x==50 ) plotss = false;

  TFile *fcombo   = TFile::Open(Form("%s_Combo_%i.root",sample,x));
  TFile *fflorida = TFile::Open(Form("%s_%i.root",sample,x));

  //TFile *fcombo_band = TFile::Open(Form("%s_Combo_%i_UNCBANDS.root",sample,x));
  TFile *fcombo_band = TFile::Open(Form("ExpectedPlusMinusLimits/%s_Combo_%i_UNCBANDS.root",sample,x));

  TFile *fcombo_theoryUp = TFile::Open(Form("TheoryPlusMinusLimits/%s_Combo_%i_THEORYUP.root",sample,x));
  TFile *fcombo_theoryDn = TFile::Open(Form("TheoryPlusMinusLimits/%s_Combo_%i_THEORYDOWN.root",sample,x));

  //TFile *fcomboNew   = TFile::Open(Form("NormalLimits/%s_Combo_%i.root",sample,x));

  TH2D*    hobs_temp       = (TH2D*)   fcombo->Get("BestObsSxBR");
  TH2D*    hobs            = cloneHist(hobs_temp);
  //TGraph*  gr_combo_obs    = (TGraph*) fcombo->Get("ObservedExclusion");
  //TGraph*  gr_combo_exp    = (TGraph*) fcombo->Get("ExpectedExclusion");
  TGraph*  gr_florida      = (TGraph*) fflorida->Get("ObservedExclusion");
  TGraph*  gr_ss           = new TGraph();


  // TGraph*  gr_combo_expp1  = (TGraph*) fcombo_band->Get("ExpectedP1SigmaExclusion");
  // TGraph*  gr_combo_expm1  = (TGraph*) fcombo_band->Get("ExpectedM1SigmaExclusion");

  // TGraph*  gr_combo_theoryUp = (TGraph*) fcombo_theoryUp->Get("ObservedExclusion_THEORYUP");
  // TGraph*  gr_combo_theoryDn = (TGraph*) fcombo_theoryDn->Get("ObservedExclusion_THEORYDOWN");

  TGraph*  gr_combo_expp1;
  TGraph*  gr_combo_expm1;

  TGraph*  gr_combo_theoryUp;
  TGraph*  gr_combo_theoryDn;

  TGraph*  gr_combo_obs;
  TGraph*  gr_combo_exp;

  if( TString(sample).Contains("LeftSlepton") ){
    if( x==25){
      gr_combo_expp1     = Left25_expectedup();
      gr_combo_expm1     = Left25_expecteddown();
      gr_combo_theoryUp  = Left25_observedup();
      gr_combo_theoryDn  = Left25_observeddown();
      gr_combo_obs       = Left25_observed();
      gr_combo_exp       = Left25_expected();
    }
    else if( x==50){
      gr_combo_expp1     = Left50_expectedup();
      gr_combo_expm1     = Left50_expecteddown();
      gr_combo_theoryUp  = Left50_observedup();
      gr_combo_theoryDn  = Left50_observeddown();
      gr_combo_obs       = Left50_observed();
      gr_combo_exp       = Left50_expected();
    }
    else if( x==75){
      gr_combo_expp1     = Left75_expectedup();
      gr_combo_expm1     = Left75_expecteddown();
      gr_combo_theoryUp  = Left75_observedup();
      gr_combo_theoryDn  = Left75_observeddown();
      gr_combo_obs       = Left75_observed();
      gr_combo_exp       = Left75_expected();
    }
  }
  else{
    if( x==25){
      gr_combo_expp1     = Tau25_expectedup();
      gr_combo_expm1     = Tau25_expecteddown();
      gr_combo_theoryUp  = Tau25_observedup();
      gr_combo_theoryDn  = Tau25_observeddown();
      gr_combo_obs       = Tau25_observed();
      gr_combo_exp       = Tau25_expected();
    }
    else if( x==50){
      gr_combo_expp1     = Tau50_expectedup();
      gr_combo_expm1     = Tau50_expecteddown();
      gr_combo_theoryUp  = Tau50_observedup();
      gr_combo_theoryDn  = Tau50_observeddown();
      gr_combo_obs       = Tau50_observed();
      gr_combo_exp       = Tau50_expected();
    }
    else if( x==75){
      gr_combo_expp1     = Tau75_expectedup();
      gr_combo_expm1     = Tau75_expecteddown();
      gr_combo_theoryUp  = Tau75_observedup();
      gr_combo_theoryDn  = Tau75_observeddown();
      gr_combo_obs       = Tau75_observed();
      gr_combo_exp       = Tau75_expected();
    }
  }

  //TGraph*  gr_combo_obs_new    = (TGraph*) fcomboNew->Get("ObservedExclusion");
  //TGraph*  gr_combo_exp_new    = (TGraph*) fcomboNew->Get("ExpectedExclusion");

  //gr_combo_obs_new->SetLineColor(4);
  //gr_combo_exp_new->SetLineColor(4);

  if( plotss ) gr_ss = (TGraph*) fcombo->Get("SSObservedExclusion");

  // 3l observed
  gr_florida->SetLineWidth(width2);
  gr_florida->SetLineStyle(7);
  gr_florida->SetLineColor(6);

  // SS observed
  gr_ss->SetLineWidth(width2);
  gr_ss->SetLineStyle(9);
  gr_ss->SetLineColor(kViolet-5);

  // combined expected +/-1 sigma
  gr_combo_expp1->SetLineColor(2);
  gr_combo_expp1->SetLineWidth(width2);
  gr_combo_expp1->SetLineStyle(3);

  gr_combo_expm1->SetLineColor(2);
  gr_combo_expm1->SetLineWidth(width2);
  gr_combo_expm1->SetLineStyle(3);

  // combined observed +/-1 sigma theory
  gr_combo_theoryUp->SetLineColor(4);
  gr_combo_theoryUp->SetLineWidth(width2);
  gr_combo_theoryUp->SetLineStyle(4);

  gr_combo_theoryDn->SetLineColor(4);
  gr_combo_theoryDn->SetLineWidth(width2);
  gr_combo_theoryDn->SetLineStyle(4);

  // combined observed
  gr_combo_obs->SetLineWidth(width1);

  // combined expected
  gr_combo_exp->SetLineWidth(width1);
  gr_combo_exp->SetLineStyle(2);

  //gr_ss_new->SetLineColor(2);
  //gr_ss_new->SetLineWidth(2);

  if( TString(sample).Contains("Tau")  && (x==25||x==75) ) gr_ss = extendGraph(gr_ss);
  if( TString(sample).Contains("Left") && (x==50) )        gr_ss = fixupGraph(gr_ss);

  //fixupGraph(gr_ss);

  TH2D *hdummy = new TH2D("hdummy","",65,100,750,72,0,725);
  
  TCanvas *can = new TCanvas(Form("%s_%i",sample,x),Form("%s_%i",sample,x),600,600);
  can->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(hdummy);
  formatHist(hobs);

  if( TString(sample).Contains("Tau") ){
    hdummy->GetXaxis()->SetRangeUser(100,475);
    hdummy->GetYaxis()->SetRangeUser(  0,450);
  }

  hdummy->Draw();
  hobs->Draw("colzsame");

  gr_combo_obs->Draw("l");
  gr_combo_exp->Draw("l");
  gr_florida->Draw("l");

  //gr_combo_obs_new->Draw("l");
  //gr_combo_exp_new->Draw("l");

  gr_combo_theoryUp->Draw("l");
  gr_combo_theoryDn->Draw("l");

  gr_combo_expp1->Draw("l");
  gr_combo_expm1->Draw("l");

  if( plotss ) gr_ss->Draw("l");

  //gr_ss_new->Draw("l");

  hobs->Draw("axissame");
  cmsPrelim(4.98,isPreliminary);

  float ymin = 0.67;
  if( !plotss ) ymin = 0.88 - (5.0/6.0) * ( 0.88-0.67 );
  
  //TLegend *leg = new TLegend(0.2,0.67,0.6,0.88);
  TLegend *leg = new TLegend(0.2,ymin,0.6,0.88);
  leg->AddEntry(gr_combo_obs      ,"Comb. observed","l");
  leg->AddEntry(gr_combo_theoryUp ,"Comb. observed (#pm1#sigma^{theory})","l");
  leg->AddEntry(gr_combo_exp      ,"Comb. median expected","l");
  leg->AddEntry(gr_combo_expp1    ,"Comb. expected (#pm1#sigma)","l");
  leg->AddEntry(gr_florida        ,"Trilepton (M_{T}) observed","l");
  if( plotss ) leg->AddEntry(gr_ss             ,"SS observed","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();
  leg->Draw();

  TLatex *thistex = new TLatex();
  thistex->SetNDC();
  thistex->SetTextSize(0.04);
  thistex->DrawLatex(0.19,0.6,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  // if     (x==25) thistex->DrawLatex(0.01,0.03,"m(#tilde{l}) = 0.25m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.75m(#tilde{#chi}_{1}^{0})");
  // else if(x==50) thistex->DrawLatex(0.01,0.03,"m(#tilde{l}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  // else if(x==75) thistex->DrawLatex(0.01,0.03,"m(#tilde{l}) = 0.75m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.25m(#tilde{#chi}_{1}^{0})");

  // if     (x==25) thistex->DrawLatex(0.01,0.03,"m_{#tilde{#font[12]{l}}} = 0.25(^{}m_{#tilde{#chi}_{2}^{0}}=m_{#tilde{#chi}_{1}^{#pm}}) + 0.75^{}m_{#tilde{#chi}_{1}^{0}}");
  // else if(x==50) thistex->DrawLatex(0.01,0.03,"m_{#tilde{#font[12]{l}}} = 0.5(^{}m_{#tilde{#chi}_{2}^{0}}=m_{#tilde{#chi}_{1}^{#pm}}) + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");
  // else if(x==75) thistex->DrawLatex(0.01,0.03,"m_{#tilde{#font[12]{l}}} = 0.75(^{}m_{#tilde{#chi}_{2}^{0}}=m_{#tilde{#chi}_{1}^{#pm}}) + 0.25^{}m_{#tilde{#chi}_{1}^{0}}");

  if     (x==25) thistex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.25^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.75^{}m_{#tilde{#chi}_{1}^{0}}");
  else if(x==50) thistex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.5^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");
  else if(x==75) thistex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.75^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.25^{}m_{#tilde{#chi}_{1}^{0}}");

  if( TString(sample).Contains("Left") ){
    thistex->DrawLatex(0.19,0.54,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{#font[12]{l}}#font[12]{l} (BF=0.5)");
    thistex->DrawLatex(0.19,0.48,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#font[12]{l}}#nu_{#font[12]{l}} , #font[12]{l}#tilde{#nu}_{#font[12]{l}}");
  }
  else{
    thistex->DrawLatex(0.19,0.54,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{#font[12]{l}}#font[12]{l} (BF=1)");
    thistex->DrawLatex(0.19,0.48,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#tau}#nu_{#tau}");
  }

  can->Modified();
  can->Update();

  if( printPlots ){
    if     ( TString(sample).Contains("Left") ) can->Print(Form("%s_%i_Fig9%s.pdf" ,sample,x,isPrelimChar));
    else if( TString(sample).Contains("Tau") )  can->Print(Form("%s_%i_Fig10%s.pdf",sample,x,isPrelimChar));

    char* fignum = "8";
    char* figlet = "a";

    if( TString(sample).Contains("Tau") ) fignum = "9";

    if( x==50 ) figlet = "b";
    if( x==75 ) figlet = "c";

    can->Print(Form("Figure%s%s.pdf",fignum,figlet));
    can->Print(Form("Figure%s%s.png",fignum,figlet));
  }

}




