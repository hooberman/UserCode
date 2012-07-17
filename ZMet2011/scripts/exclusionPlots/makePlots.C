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


using namespace std;

void makeFloridaPlot(char* sample, int x );

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

void formatHist( TH2D* hist ){

  hist->GetXaxis()->SetTitle("m(#tilde{#chi}_{2}^{0}) = m(#tilde{#chi}_{1}^{#pm}) [GeV]");
  hist->GetYaxis()->SetTitle("m(#tilde{#chi}_{1}^{0}) [GeV]");
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

void cmsPrelim( double intLumi )
{

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.18,0.93,"CMS Preliminary,  #sqrt{s}=7 TeV,  L_{int}=4.98 fb^{-1}");
}

void makePlots(){

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.03);

  //-------------------------------------------
  // Rutgers/KIT
  //-------------------------------------------

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
  
  TCanvas *can_2i = new TCanvas();
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
  gr2i_1->Draw("l");
  gr2i_2->Draw("l");
  h2i0->Draw("axissame");
  //gr2i_3->Draw("l");
  //gr2i_4->Draw("l");
  cmsPrelim(4.98);

  TLegend *leg = new TLegend(0.2,0.78,0.6,0.88);
  leg->AddEntry(gr2i_1,"observed","l");
  leg->AddEntry(gr2i_2,"median expected","l");
  //leg->AddEntry(gr2i_3,"expected #pm1#sigma","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();


  tex->SetTextSize(0.03);
  tex->DrawLatex(0.18,0.70,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  tex->DrawLatex(0.18,0.65,"m(#tilde{l}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  tex->DrawLatex(0.18,0.60,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l (BF=50%)");
  tex->DrawLatex(0.18,0.55,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{l}#nu_{l}");


  can_2i->Modified();
  can_2i->Update();
  can_2i->Print("multilepton_flavordemocratic_Fig7.pdf");


  //-----------------
  // model 2a
  //-----------------

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
  gr2a_1->Draw("l");
  gr2a_2->Draw("l");
  h2a0->Draw("axissame");
  cmsPrelim(4.98);

  leg->Draw();

  tex->SetTextSize(0.03);
  tex->DrawLatex(0.18,0.70,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  tex->DrawLatex(0.18,0.65,"m(#tilde{l}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  tex->DrawLatex(0.18,0.60,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l (BF=100%)");
  tex->DrawLatex(0.18,0.55,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#tau}#nu_{#tau}");

  can_2a->Modified();
  can_2a->Update();
  can_2a->Print("multilepton_tauenriched_Fig8.pdf");

  //-----------------
  // TChiWZ
  //-----------------

  TFile *fwz  = TFile::Open("combinePlots_VZ_Trilepton.root");

  TH2F*   hwz             = (TH2F*)   fwz->Get("hexcl");
  TGraph* grwz_combo      = (TGraph*) fwz->Get("gr_combo");
  TGraph* grwz_combo_exp  = (TGraph*) fwz->Get("gr_combo_exp");
  TGraph* grwz_tri        = (TGraph*) fwz->Get("gr_tri");
  TGraph* grwz_vzmet      = (TGraph*) fwz->Get("gr_vzmet");
  
  TCanvas *can_wz = new TCanvas();
  can_wz->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist((TH2D*)hwz);
  hwz->GetXaxis()->SetRangeUser(87.5,512.5);
  hwz->SetMinimum(50);
  hwz->SetMaximum(5000);
  hwz->Draw("colz");
  grwz_combo->Draw("l");
  grwz_combo_exp->Draw("l");
  grwz_tri->Draw("l");
  grwz_vzmet->Draw("l");
  cmsPrelim(4.98);

  TLegend *legwz = new TLegend(0.2,0.68,0.7,0.88);
  legwz->AddEntry(grwz_combo    ,"combined observed","l");
  legwz->AddEntry(grwz_combo_exp,"combined median expected","l");
  legwz->AddEntry(grwz_vzmet    ,"2l2j observed","l");
  legwz->AddEntry(grwz_tri      ,"trilepton (M_{T}) observed","l");
  legwz->SetBorderSize(0);
  legwz->SetFillColor(0);
  legwz->Draw();
  legwz->Draw();
  
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.18,0.6,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm} #rightarrow WZ+E_{T}^{miss}");

  can_wz->Modified();
  can_wz->Update();
  can_wz->Print("WZ_Fig11.pdf");

  //-----------------------------
  // Florida/ETH plots
  //-----------------------------

  makeFloridaPlot("LeftSlepton",25);
  makeFloridaPlot("LeftSlepton",50);
  makeFloridaPlot("LeftSlepton",75);

  makeFloridaPlot("TauEnriched",25);
  makeFloridaPlot("TauEnriched",50);
  makeFloridaPlot("TauEnriched",75);

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

void makeFloridaPlot(char* sample, int x ){


  bool plotss = true;
  if( TString(sample).Contains("Tau") && x==50 ) plotss = false;

  TFile *fcombo   = TFile::Open(Form("%s_Combo_%i.root",sample,x));
  TFile *fflorida = TFile::Open(Form("%s_%i.root",sample,x));

  TH2D*    hobs_temp       = (TH2D*)   fcombo->Get("BestObsSxBR");
  TH2D*    hobs            = cloneHist(hobs_temp);
  TGraph*  gr_combo_obs    = (TGraph*) fcombo->Get("ObservedExclusion");
  TGraph*  gr_combo_exp    = (TGraph*) fcombo->Get("ExpectedExclusion");
  TGraph*  gr_florida      = (TGraph*) fflorida->Get("ObservedExclusion");
  TGraph*  gr_ss           = new TGraph();

  if( plotss ) gr_ss = (TGraph*) fcombo->Get("SSObservedExclusion");

  gr_ss->SetLineWidth(4);
  gr_ss->SetLineStyle(3);
  gr_ss->SetLineColor(6);

  gr_combo_obs->SetLineWidth(4);

  gr_combo_exp->SetLineWidth(4);
  gr_combo_exp->SetLineStyle(2);

  gr_florida->SetLineWidth(4);
  gr_florida->SetLineStyle(4);
  gr_florida->SetLineColor(2);

  TH2D *hdummy = new TH2D("hdummy","",65,100,750,72,0,725);
  
  TCanvas *can = new TCanvas(Form("%s_%i",sample,x),Form("%s_%i",sample,x),600,600);
  can->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(hdummy);
  formatHist(hobs);
  hdummy->Draw();
  hobs->Draw("colzsame");

  gr_combo_obs->Draw("l");
  gr_combo_exp->Draw("l");
  gr_florida->Draw("l");
  if( plotss ) gr_ss->Draw("l");
  hobs->Draw("axissame");
  cmsPrelim(4.98);

  
  TLegend *leg = new TLegend(0.2,0.73,0.65,0.88);
  leg->AddEntry(gr_combo_obs    ,"combined observed","l");
  leg->AddEntry(gr_combo_exp    ,"combined median expected","l");
  leg->AddEntry(gr_florida      ,"trilepton (M_{T}) observed","l");
  leg->AddEntry(gr_ss           ,"SS observed","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();
  leg->Draw();

  TLatex *thistex = new TLatex();
  thistex->SetNDC();
  thistex->SetTextSize(0.03);
  thistex->DrawLatex(0.18,0.7,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  if     ( x == 25 ) thistex->DrawLatex(0.18,0.65,"m(#tilde{l}) = 0.25m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.75m(#tilde{#chi}_{1}^{0})");
  else if( x == 50 ) thistex->DrawLatex(0.18,0.65,"m(#tilde{l}) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");
  else if( x == 75 ) thistex->DrawLatex(0.18,0.65,"m(#tilde{l}) = 0.75m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.25m(#tilde{#chi}_{1}^{0})");
  if( TString(sample).Contains("Left") ){
    thistex->DrawLatex(0.18,0.60,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l (BF=50%)");
    thistex->DrawLatex(0.18,0.55,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{l}#nu_{l}");
  }
  else{
    thistex->DrawLatex(0.18,0.60,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l (BF=100%)");
    thistex->DrawLatex(0.18,0.55,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#tau}#nu_{#tau}");
  }

  can->Modified();
  can->Update();
  if     ( TString(sample).Contains("Left") ) can->Print(Form("%s_%i_Fig9.pdf" ,sample,x));
  else if( TString(sample).Contains("Tau") )  can->Print(Form("%s_%i_Fig10.pdf",sample,x));
  
}
