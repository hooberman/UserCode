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


void makeZZPlot(){

  const unsigned int n = 15;
  float x[n];
  float y1[n];
  float y2[n];

  x[0] =130;  y1[0] =27952;  y2[0] =3057;                 
  x[1] =150;  y1[1] = 2951;  y2[1] =1702;              
  x[2] =170;  y1[2] =  993;  y2[2] = 994;           
  x[3] =190;  y1[3] =  501;  y2[3] = 617;        
  x[4] =210;  y1[4] =  348;  y2[4] = 398;     
  x[5] =230;  y1[5] =  287;  y2[5] = 264;  
  x[6] =250;  y1[6] =  240;  y2[6] = 182;
  x[7] =270;  y1[7] =  204;  y2[7] = 128;
  x[8] =290;  y1[8] =  197;  y2[8] =  91;
  x[9] =310;  y1[9] =  179;  y2[9] =  67;
  x[10]=330;  y1[10]=  164;  y2[10]=  49;
  x[11]=350;  y1[11]=  159;  y2[11]=  37;
  x[12]=370;  y1[12]=  152;  y2[12]=  28;
  x[13]=390;  y1[13]=  149;  y2[13]=  21;
  x[14]=410;  y1[14]=  142;  y2[14]=  17;
	   
  TGraph* g1 = new TGraph(n,x,y2);
  TGraph* g2 = new TGraph(n,x,y1);

  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  //gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",100,130,410,100,0,5000);
  hdummy->Draw();

  c1->cd();
  g1->SetLineColor(4);
  g2->SetLineColor(2);
  g1->SetLineWidth(2);
  g2->SetLineWidth(5);
  g2->SetLineStyle(2);

  hdummy->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma #times BR [fb]");
  hdummy->GetYaxis()->SetLabelSize(0.04);
  hdummy->GetXaxis()->SetLabelSize(0.04);

  TBox* box = new TBox();
  //box->SetBorderStyle(2);
  //box->SetBorderSize(1);
  //box->SetFillColor(5);
  //box->SetFillStyle(3002);
  box->DrawBox(171,0,221,5000);
  TLine line;
  line.DrawLine(171,0,171,5000);
  line.DrawLine(221,0,221,5000);

  hdummy->Draw("axissame");
  g1->SetMinimum(0);
  g1->SetMaximum(5000);
  g1->Draw("samec");
  g2->Draw("samel");

  //g1->Draw("Al");
  //g2->Draw("samel");
  
  
  TLegend *leg = new TLegend(0.5,0.6,0.9,0.8);
  leg->AddEntry(g2,"observed UL","l");
  leg->AddEntry(g1,"theory","l");
  leg->AddEntry(box,"excluded region","f");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();								
  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
  t->SetTextSize(0.04);
  t->DrawLatex(0.46,0.5,"#chi^{0}#chi^{0} #rightarrow ZZ + E_{T}^{miss}");
  t->SetTextSize(0.038);
  t->DrawLatex(0.46,0.45,"GGMSB Z-enriched higgsino");

  c1->Print("GMSB.pdf");
  c1->Print("GMSB.png");
  c1->Print("GMSB.eps");
  gROOT->ProcessLine(".! ps2pdf GMSB.eps GMSB_ppt.pdf");
}	   
