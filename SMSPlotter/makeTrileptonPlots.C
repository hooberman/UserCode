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

#include "contours.C"
#include "tdrstyle_SUSY.C"

using namespace std;

bool isPreliminary = false;
char* isPrelimChar = (char*) "";

int width1 = 8;
int width2 = 4;

void formatHist( TH2D* hist ){

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

void makeTrileptonPlots( bool printPlots = false){

  //gROOT->ProcessLine(".L tdrstyle_SUSY.C");
  setTDRStyle();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.05);


  if( isPreliminary ) isPrelimChar = (char*) "_prelim";

  //-------------------------------------------
  // model 2i (flavor-democratic)
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

  TFile *file_2i  = TFile::Open("trilepton_2i.root");

  TH2D*   h2i     = (TH2D*)   file_2i->Get("xSecObserved");
  TH2D*   h2i1    = (TH2D*)   file_2i->Get("xSecObserved1");
  TH2D*   h2i0    = (TH2D*)   file_2i->Get("xSecObserved0");
  
  TCanvas *can_2i = new TCanvas("can_2i","can_2i",600,600);
  can_2i->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(h2i0);
  h2i0->Draw("colz");
  h2i->Draw("colsame");
  h2i1->Draw("colsame");

  mod2i_observed->Draw("l");
  mod2i_expected->Draw("l");
  mod2i_expectedP1->Draw("l");
  mod2i_expectedM1->Draw("l");
  mod2i_observedP->Draw("l");
  mod2i_observedM->Draw("l");

  h2i0->Draw("axissame");
  cmsPrelim(4.98,isPreliminary);

  TLegend *leg = new TLegend(0.2,0.72,0.6,0.88);
  leg->AddEntry(mod2i_observed  ,"observed","l");
  leg->AddEntry(mod2i_observedP ,"observed (#pm1#sigma^{theory})","l");
  leg->AddEntry(mod2i_expected  ,"median expected","l");
  leg->AddEntry(mod2i_expectedP1,"expected (#pm1#sigma)","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  tex->SetTextSize(0.04);
  tex->DrawLatex(0.19,0.65,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  tex->DrawLatex(0.19,0.59,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{#font[12]{l}}#font[12]{l} (BF=0.5)");
  tex->DrawLatex(0.19,0.53,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#font[12]{l}}#nu_{#font[12]{l}} , #font[12]{l}#tilde{#nu}_{#font[12]{l}}");
  tex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.5^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");

  can_2i->Modified();
  can_2i->Update();
  if( printPlots) can_2i->Print(Form("multilepton_flavordemocratic_Fig7%s.pdf",isPrelimChar));
  if( printPlots){
    can_2i->Print("Figure7a.pdf");
    can_2i->Print("Figure7a.png");
  }

  //-------------------------------
  // model 2a (tau-enriched)
  //-------------------------------

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

  TFile *file_2a  = TFile::Open("trilepton_2a.root");

  TH2D*   h2a     = (TH2D*)   file_2a->Get("xSecObserved");
  TH2D*   h2a0    = (TH2D*)   file_2a->Get("xSecObserved0");
  
  TCanvas *can_2a = new TCanvas();
  can_2a->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(h2a0);
  h2a0->Draw("colz");
  h2a->Draw("colsame");

  mod2a_observed->Draw("l");
  mod2a_expected->Draw("l");
  mod2a_expectedP1->Draw("l");
  mod2a_expectedM1->Draw("l");
  mod2a_observedP->Draw("l");
  mod2a_observedM->Draw("l");

  h2a0->Draw("axissame");
  cmsPrelim(4.98,isPreliminary);

  leg->Draw();

  tex->SetTextSize(0.04);
  tex->DrawLatex(0.19,0.65,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  tex->DrawLatex(0.19,0.59,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{#font[12]{l}}#font[12]{l} (BF=1)");
  tex->DrawLatex(0.19,0.53,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{#tau}#nu_{#tau}");
  tex->DrawLatex(0.10,0.03,"m_{#tilde{#font[12]{l}}} = 0.5^{}m_{#tilde{#chi}_{1}^{#pm}} + 0.5^{}m_{#tilde{#chi}_{1}^{0}}");

  can_2a->Modified();
  can_2a->Update();
  if( printPlots) can_2a->Print(Form("multilepton_tauenriched_Fig8%s.pdf",isPrelimChar));
  if( printPlots){
    can_2a->Print("Figure7b.pdf");
    can_2a->Print("Figure7b.png");
  }
}
