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

  x[0] =130;  y1[0] =26580;  y2[0] =3057;                 
  x[1] =150;  y1[1] =2734 ;  y2[1] =1702;              
  x[2] =170;  y1[2] =916  ;  y2[2] = 994;           
  x[3] =190;  y1[3] =461  ;  y2[3] = 617;        
  x[4] =210;  y1[4] =321  ;  y2[4] = 398;     
  x[5] =230;  y1[5] =264  ;  y2[5] = 264;  
  x[6] =250;  y1[6] =221  ;  y2[6] = 182;
  x[7] =270;  y1[7] =188  ;  y2[7] = 128;
  x[8] =290;  y1[8] =182  ;  y2[8] =  91;
  x[9] =310;  y1[9] =166  ;  y2[9] =  67;
  x[10]=330;  y1[10]=152  ;  y2[10]=  49;
  x[11]=350;  y1[11]=147  ;  y2[11]=  37;
  x[12]=370;  y1[12]=140  ;  y2[12]=  28;
  x[13]=390;  y1[13]=138  ;  y2[13]=  21;
  x[14]=410;  y1[14]=132  ;  y2[14]=  17;
	   
  TGraph* g1 = new TGraph(n,x,y1);
  TGraph* g2 = new TGraph(n,x,y2);

  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.1);
  //gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",100,130,410,100,0,5000);
  hdummy->Draw();

  c1->cd();
  g1->SetLineColor(4);
  g2->SetLineColor(2);
  g1->SetLineWidth(3);
  g2->SetLineWidth(3);

  hdummy->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma #times BR [GeV]");
  hdummy->GetYaxis()->SetLabelSize(0.04);

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
  g1->SetMinimum(0);
  g1->SetMaximum(5000);
  g1->Draw("samel");
  g2->Draw("samel");

  //g1->Draw("Al");
  //g2->Draw("samel");
  
  
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->AddEntry(g1,"theory","l");
  leg->AddEntry(g2,"observed UL","l");
  leg->AddEntry(box,"excluded region","f");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();								
  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
  t->SetTextSize(0.04);
  t->DrawLatex(0.47,0.45,"GGMSB Z-enriched higgsino");
  t->DrawLatex(0.47,0.4,"#chi^{0}#chi^{0} #rightarrow ZZ + E_{T}^{miss}");

  c1->Print("GMSB.pdf");
  c1->Print("GMSB.png");
  c1->Print("GMSB.eps");
  gROOT->ProcessLine(".! ps2pdf GMSB.eps GMSB_ppt.pdf");
}	   
