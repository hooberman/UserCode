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


void makeWZPlot(){

  const unsigned int n = 5;
  float x[n];
  float y0[n];
  float y50[n];
  float yw[n];
  float yh[n];

  x[0] = 200;  y0[0] =336;  y50[0] =399;  yw[0]=614; yh[0]=307;
  x[1] = 225;  y0[1] =231;  y50[1] =269;  yw[1]=377; yh[1]=189;
  x[2] = 250;  y0[2] =211;  y50[2] =231;  yw[2]=239; yh[2]=120;
  x[3] = 275;  y0[3] =178;  y50[3] =192;  yw[3]=157; yh[3]=79;
  x[4] = 300;  y0[4] =152;  y50[4] =167;  yw[4]=106; yh[4]=53;    
	   
  TGraph* g0  = new TGraph(n,x,y0);
  TGraph* g50 = new TGraph(n,x,y50);
  TGraph* gw  = new TGraph(n,x,yw);
  TGraph* gh  = new TGraph(n,x,yh);

  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  gPad->SetGridx();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",100,200,300,100,0,700);
  hdummy->Draw();

  c1->cd();

  g0->SetLineColor(2);
  g50->SetLineColor(2);
  gw->SetLineColor(4);
  gh->SetLineColor(4);

  g0->SetLineWidth(3);
  g50->SetLineWidth(3);
  gw->SetLineWidth(3);
  gh->SetLineWidth(3);

  g50->SetLineStyle(2);
  gh->SetLineStyle(2);

  hdummy->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma #times BR [GeV]");
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
  gw->Draw("samel");
  gh->Draw("samel");

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
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
  t->SetTextSize(0.04);
  //t->DrawLatex(0.47,0.45,"");
  t->DrawLatex(0.57,0.5,"#chi^{#pm}#chi^{0} #rightarrow WZ + E_{T}^{miss}");

  c1->Print("WZ.pdf");
  c1->Print("WZ.png");
  c1->Print("WZ.eps");
  gROOT->ProcessLine(".! ps2pdf WZ.eps WZ_ppt.pdf");
}	   
