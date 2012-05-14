#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
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
#include "TBox.h"
#include <sstream>
#include <iomanip>

using namespace std;


void makeWZPlot( bool printplots = false ){

  TFile *f = TFile::Open("cards/V00-02-09/observed_limit_line.root");

  TGraph* g0  = (TGraph*) f->Get("grobs0");
  TGraph* g50 = (TGraph*) f->Get("grobs50");

  //TGraph* g0  = (TGraph*) f->Get("grexp0");
  //TGraph* g50 = (TGraph*) f->Get("grexp50");

  const unsigned int n = 12;
  float x[n];
  float yw[n];
  float yh[n];

  x[0]  = 125;   yw[0]=3960; yh[0]=1980;
  x[1]  = 150;   yw[1]=1949; yh[1]=975;
  x[2]  = 175;   yw[2]=1052; yh[2]=526;
  x[3]  = 200;   yw[3]=614;  yh[3]=307;
  x[4]  = 225;   yw[4]=377;  yh[4]=189;
  x[5]  = 250;   yw[5]=239;  yh[5]=120;
  x[6]  = 275;   yw[6]=157;  yh[6]=79;
  x[7]  = 300;   yw[7]=106;  yh[7]=53;    
  x[8]  = 325;   yw[8]=73;   yh[8]=36;
  x[9]  = 350;   yw[9]=51;   yh[9]=26;
  x[10] = 375;   yw[10]=37;  yh[10]=19;
  x[11] = 400;   yw[11]=26;  yh[11]=13;    

  TGraph* gw  = new TGraph(n,x,yw);
  TGraph* gh  = new TGraph(n,x,yh);

  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  gPad->SetGridx();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",100,125,400,100,0,5000);
  hdummy->Draw();

  c1->cd();

  g0->SetLineColor(2);
  g50->SetLineColor(2);
  gw->SetLineColor(4);
  gh->SetLineColor(4);

  g0->SetLineWidth(4);
  g50->SetLineWidth(4);
  gw->SetLineWidth(2);
  gh->SetLineWidth(2);

  g50->SetLineStyle(2);
  gh->SetLineStyle(2);

  hdummy->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma #times BR [fb]");
  hdummy->GetYaxis()->SetLabelSize(0.04);
  hdummy->GetXaxis()->SetLabelSize(0.04);

  /*
  TBox* box = new TBox();
  //box->SetBorderStyle(2);
  //box->SetBorderSize(1);
  //box->SetFillColor(5);
  //box->SetFillStyle(3002);
  box->DrawBox(169,0,230,5000);
  TLine line;
  line.DrawLine(169,0,169,5000);
  line.DrawLine(230,0,230,5000);
  hdummy->Draw("axissame");
  */

  g0->Draw("samel");
  g50->Draw("samel");
  gw->Draw("samec");
  gh->Draw("samec");

  // g1->SetMinimum(0);
  // g1->SetMaximum(5000);
  // g1->Draw("samel");
  // g2->Draw("samel");

  //g1->Draw("Al");
  //g2->Draw("samel");
  
  TLegend *leg = new TLegend(0.4,0.6,0.9,0.8);
  leg->AddEntry(g0 ,"observed UL m_{LSP} = 0 GeV","l");
  leg->AddEntry(g50,"observed UL m_{LSP} = 50 GeV","l");
  leg->AddEntry(gw,"theory wino-like: g#gamma^{#mu}","l");
  leg->AddEntry(gh,"theory higgsino-like: g#gamma^{#mu} / #sqrt{2}","l");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();								
  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");
  t->SetTextSize(0.04);
  //t->DrawLatex(0.47,0.45,"");
  t->DrawLatex(0.57,0.5,"#chi^{#pm}#chi^{0} #rightarrow WZ + E_{T}^{miss}");

  if( printplots ){
    c1->Print("WZnew.pdf");
    c1->Print("WZnew.png");
    c1->Print("WZnew.eps");
    gROOT->ProcessLine(".! ps2pdf WZnew.eps WZnew_ppt.pdf");
  }	   
}
