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

bool plotExpected = false;
bool plotObserved = true;
bool logplot      = false;

void makeGMSBPlot( bool printplots = false ){

  // VZ+MET exclusion
  TFile *f       = TFile::Open("cards/V00-02-08/observed_limit.root");
  TGraph* gul    = (TGraph*) f->Get("grobs");
  TGraph* gulexp = (TGraph*) f->Get("grexp");

  // Rutgers exclusion
  //TFile *frutgers = TFile::Open("20120411_UCSD_GMSB_datacard/observed_limit.root ");
  //TFile *frutgers = TFile::Open("20120419_UCSD_GMSB_datacard/observed_limit.root ");
  TFile *frutgers = TFile::Open("cards/20120420_UCSD_GMSB_datacard/observed_limit.root ");
  TGraph* gul2    = (TGraph*) frutgers->Get("grobs");
  TGraph* gul2exp = (TGraph*) frutgers->Get("grexp");

  // VZ+MET exclusion
  TFile *fc       = TFile::Open("cards/V00-02-08/observed_limit_combined.root");
  TGraph* gulc    = (TGraph*) fc->Get("grobs");
  TGraph* gulcexp = (TGraph*) fc->Get("grexp");


  Double_t xp;
  Double_t yp;

  Double_t xp2;
  Double_t yp2;

  Double_t xpc;
  Double_t ypc;

  cout << setw(15) << "mass"        << setw(4) << "&"
       << setw(15) << "\\wzzmet"    << setw(4) << "&"
       << setw(15) << "mult-lepton" << setw(4) << "&"
       << setw(15) << "combined"    << setw(4) << "&"
       << setw(15) << "asdf"        << setw(4) << "\\\\" << endl;

  for( int i = 0 ; i < 15 ; ++i ){

    // gulexp->GetPoint ((Int_t) i,xp,yp);
    // gul2exp->GetPoint((Int_t) i,xp2,yp2);
    // gulcexp->GetPoint((Int_t) i,xpc,ypc);
    // float exp = 1.0 / sqrt( 1.0/(yp*yp) + 1.0/(yp2*yp2) ); 

    gul->GetPoint ((Int_t) i,xp,yp);
    gul2->GetPoint((Int_t) i,xp2,yp2);
    gulc->GetPoint((Int_t) i,xpc,ypc);
    float exp = 1.0 / sqrt( 1.0/(yp*yp) + 1.0/(yp2*yp2) ); 


    // gul->GetPoint ((Int_t) i,xp,yp);
    // gul2->GetPoint((Int_t) i,xp2,yp2);
    // gulc->GetPoint((Int_t) i,xpc,ypc);

    cout << setw(15) << xp               << setw(4) << "&"
	 << setw(15) << Form("%.0f",yp)  << setw(4) << "&"
	 << setw(15) << Form("%.0f",yp2) << setw(4) << "&"
	 << setw(15) << Form("%.0f",ypc) << setw(4) << "&"
	 << setw(15) << Form("%.0f",exp) << setw(4) << "\\\\" << endl;
    

    // cout << "mass    " << Form("%.0f",xp) << endl;
    // cout << "VZ+MET  " << Form("%.0f",yp) << endl;
    // cout << "4l      " << Form("%.0f",yp2) << endl;
    // cout << "combo   " << Form("%.0f",ypc) << endl;
    // cout << "exp     " << Form("%.0f",exp) << endl << endl;    



  }


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

  float ymin = 0;
  if( logplot ) ymin = 80;

  TH2F* hdummy = new TH2F("hdummy","",100,130,300,100,ymin,3000);
  hdummy->Draw();

  c1->cd();
  gPad->SetLogy();

  g->SetLineColor(4);
  g->SetLineWidth(2);

  gul->SetLineColor(2);
  gul->SetLineWidth(4);
  gulexp->SetLineColor(2);
  gulexp->SetLineWidth(4);

  gul2->SetLineWidth(4);
  gul2exp->SetLineWidth(4);

  gulc->SetLineWidth(4);
  gulc->SetLineColor(kGreen+2);
  gulcexp->SetLineWidth(4);
  gulcexp->SetLineColor(kGreen+2);

  gul->SetLineStyle(2);
  gulexp->SetLineStyle(2);
  gul2->SetLineStyle(2);
  gul2exp->SetLineStyle(2);

  // gulexp->SetLineWidth(2);
  // gulexp->SetLineStyle(2);
  // gulexp->SetLineColor(2);  
  // gul2exp->SetLineWidth(2);
  // gul2exp->SetLineStyle(2);


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

  if( plotObserved ){
    gul->Draw("samel");
    gul2->Draw("samel");
    gulc->Draw("samel");
  }

  if( plotExpected ){
    gulexp->Draw("samel");
    gul2exp->Draw("samel");
    gulcexp->Draw("samel");
  }

  g->Draw("samec");

  //gulexp->Draw("samel");
  //gul2exp->Draw("samel");

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
  //box->DrawBox(xmin,0,xmax,5000);

  TLine line;
  //line.DrawLine(xmin,0,xmin,5000);
  //line.DrawLine(xmax,0,xmax,5000);

  hdummy->Draw("axissame");
  // g->SetMinimum(0);
  // g->SetMaximum(3000);
  // g->Draw("samec");
  // gul->Draw("samel");
  // gul2->Draw("samel");

  TLegend *leg = new TLegend(0.45,0.7,0.95,0.9);
  if( plotObserved ){
    leg->AddEntry(gul  ,"observed UL (VZ+E_{T}^{miss})","l");
    leg->AddEntry(gul2 ,"observed UL (multi-lepton)","l");
    leg->AddEntry(gulc ,"observed UL (combined)","l");
  }
  if( plotExpected ){
    leg->AddEntry(gulexp  ,"expected UL (VZ+E_{T}^{miss})","l");
    leg->AddEntry(gul2exp ,"expected UL (multi-lepton)","l");
    leg->AddEntry(gulcexp ,"expected UL (combined)","l");
  }

  leg->AddEntry(g,  "theory","l");
  //leg->AddEntry(box,"excluded region","f");
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
  t->DrawLatex(0.57,0.63,"GMSB  ZZ + E_{T}^{miss}");

  if( printplots ){
    c1->Print("GMSB.pdf");
    c1->Print("GMSB.png");
    c1->Print("GMSB.eps");
    gROOT->ProcessLine(".! ps2pdf GMSB.eps GMSB_ppt.pdf");
  }	   
}
