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


void makeGMSBPlot( bool printplots = false ){

  TFile *f = TFile::Open("cards/V00-02-04/observed_limit.root");

  TGraph* gul  = (TGraph*) f->Get("grobs");
  //TGraph* g = (TGraph*) f->Get("grexp");

  const unsigned int n = 15;
  float x[n];
  float y[n];

  x[0]  = 130;   y[0]  = 3057;
  x[1]  = 150;   y[1]  = 1719; 
  x[2]  = 170;   y[2]  = 1035; 
  x[3]  = 190;   y[3]  =  656;  
  x[4]  = 210;   y[4]  =  433;  
  x[5]  = 230;   y[5]  =  293;  
  x[6]  = 250;   y[6]  =  205;  
  x[7]  = 270;   y[7]  =  146;     
  x[8]  = 290;   y[8]  =  105;   
  x[9]  = 310;   y[9]  =   77;   
  x[10] = 330;   y[10] =   57;  
  x[11] = 350;   y[11] =   43;      
  x[12] = 370;   y[12] =   33;      
  x[13] = 390;   y[13] =   25;      
  x[14] = 410;   y[14] =   20;      

  TGraph* g  = new TGraph(n,x,y);

  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  gPad->SetGridx();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",100,130,410,100,0,5000);
  hdummy->Draw();

  c1->cd();

  gul->SetLineColor(2);
  g->SetLineColor(4);

  gul->SetLineWidth(4);
  g->SetLineWidth(2);

  hdummy->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma [fb]");
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

  gul->Draw("samel");
  g->Draw("samec");

  // g1->SetMinimum(0);
  // g1->SetMaximum(5000);
  // g1->Draw("samel");
  // g2->Draw("samel");

  //g1->Draw("Al");
  //g2->Draw("samel");
  
  float xmin = 165;
  float xmax = 239;

  TBox* box = new TBox();
  //box->SetBorderStyle(2);
  //box->SetBorderSize(1);
  //box->SetFillColor(5);
  //box->SetFillStyle(3002);
  box->DrawBox(xmin,0,xmax,5000);

  TLine line;
  line.DrawLine(xmin,0,xmin,5000);
  line.DrawLine(xmax,0,xmax,5000);

  hdummy->Draw("axissame");
  g->SetMinimum(0);
  g->SetMaximum(5000);
  g->Draw("samec");
  gul->Draw("samel");

  TLegend *leg = new TLegend(0.55,0.6,0.9,0.8);
  leg->AddEntry(gul ,"observed UL","l");
  leg->AddEntry(g,  "theory","l");
  leg->AddEntry(box,"excluded region","f");
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
  t->DrawLatex(0.57,0.5,"GMSB  ZZ + E_{T}^{miss}");

  if( printplots ){
    c1->Print("GMSB.pdf");
    c1->Print("GMSB.png");
    c1->Print("GMSB.eps");
    gROOT->ProcessLine(".! ps2pdf GMSB.eps GMSB_ppt.pdf");
  }	   
}
