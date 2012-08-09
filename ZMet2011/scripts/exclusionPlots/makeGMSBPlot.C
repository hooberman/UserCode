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

bool plotExpected  = false;
bool plotObserved  = true;
bool logplot       = true;
bool isPreliminary = true;

void cmsPrelim(double intLumi, bool prelim)
{
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62);
        if(prelim) latex.SetTextSize(0.04);
        else       latex.SetTextSize(0.045);

        latex.SetTextAlign(11); // align left
        //if(prelim) latex.DrawLatex(0.13,0.92,"CMS Preliminary");
        //else       latex.DrawLatex(0.13,0.92,"CMS");

        if(prelim) latex.DrawLatex(0.19,0.92,"CMS Preliminary");
        else       latex.DrawLatex(0.19,0.92,"CMS");

        latex.SetTextAlign(31); // align right
        //latex.DrawLatex(0.89, 0.92, Form("#sqrt{s} = 7 TeV, L_{int} = %4.2f fb^{-1}", intLumi));
        latex.DrawLatex(0.92, 0.92, Form("#sqrt{s} = 7 TeV, L_{int} = %4.2f fb^{-1}", intLumi));
}


void getUncertainties(){

  ifstream *ip = new ifstream();
  ifstream *im = new ifstream();

  ip->open("7TeVc1pn2_finer_less.dat");
  im->open("7TeVc1mn2_finer_less.dat");

  float m;
  float xsec;
  float unc;
  
  while( *ip >> m >> xsec >> unc ){
    cout << m << " " << xsec << " " << unc << " " << Form("%.3f",unc/xsec) << endl;
  }

  while( *im >> m >> xsec >> unc ){
    cout << m << " " << xsec << " " << unc << " " << Form("%.3f",unc/xsec) << endl;
  }

}

TGraph* uncertaintyBand( TGraph* gup , TGraph* gdn ){

  float x[30];
  float y[30];

  Double_t thisx;
  Double_t thisy;

  for( int i = 0 ; i < 15; ++i ){
    gup->GetPoint(i,thisx,thisy);
    x[i] = thisx;
    y[i] = thisy;
    //cout << x[i] << " " << y[i] << endl;
  }

  for( int i = 0 ; i < 15; ++i ){
    gdn->GetPoint(14-i,thisx,thisy);
    x[i+15] = thisx;
    y[i+15] = thisy;
    //cout << x[i+15] << " " << y[i+15] << endl;
  }

  TGraph* gr = new TGraph(30,x,y);
  gr->SetFillColor(7);

  return gr;

}

void makeGMSBPlot( bool printplots = false ){

  //getUncertainties();

  // VZ+MET exclusion
  TFile *f       = TFile::Open("/tas/benhoob/home/LandS/VZMet_LandS/fullShapeAnalysis/cards/V00-02-08/observed_limit.root");
  TGraph* gul    = (TGraph*) f->Get("grobs");
  TGraph* gulexp = (TGraph*) f->Get("grexp");

  // Rutgers exclusion
  //TFile *frutgers = TFile::Open("20120411_UCSD_GMSB_datacard/observed_limit.root ");
  //TFile *frutgers = TFile::Open("20120419_UCSD_GMSB_datacard/observed_limit.root ");
  TFile *frutgers = TFile::Open("/tas/benhoob/home/LandS/VZMet_LandS/fullShapeAnalysis/cards/20120420_UCSD_GMSB_datacard/observed_limit.root ");
  TGraph* gul2    = (TGraph*) frutgers->Get("grobs");
  TGraph* gul2exp = (TGraph*) frutgers->Get("grexp");

  // VZ+MET exclusion
  TFile *fc       = TFile::Open("/tas/benhoob/home/LandS/VZMet_LandS/fullShapeAnalysis/cards/V00-02-08/observed_limit_combined_band.root");
  TGraph* gulc      = (TGraph*) fc->Get("grobs");
  TGraph* gulcexp   = (TGraph*) fc->Get("grexp");
  TGraph* gulcexpp1 = (TGraph*) fc->Get("grexpp1");
  TGraph* gulcexpm1 = (TGraph*) fc->Get("grexpm1");
  TGraph* gulcband  = uncertaintyBand( gulcexpp1 , gulcexpm1 );

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

    gulexp->GetPoint ((Int_t) i,xp,yp);
    gul2exp->GetPoint((Int_t) i,xp2,yp2);
    gulcexp->GetPoint((Int_t) i,xpc,ypc);
    float exp = 1.0 / sqrt( 1.0/(yp*yp) + 1.0/(yp2*yp2) ); 

    // gul->GetPoint ((Int_t) i,xp,yp);
    // gul2->GetPoint((Int_t) i,xp2,yp2);
    // gulc->GetPoint((Int_t) i,xpc,ypc);
    // float exp = 1.0 / sqrt( 1.0/(yp*yp) + 1.0/(yp2*yp2) ); 


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
  float yup[n];
  float ydn[n];

  float xerr[n];
  float yerr[n];

  float xband[30];
  float yband[30];

  x[0]  = 130;   y[0]  = 3057;   yerr[0]  = 0.055 * y[0];
  x[1]  = 150;   y[1]  = 1719;   yerr[1]  = 0.054 * y[1];
  x[2]  = 170;   y[2]  = 1035;   yerr[2]  = 0.050 * y[2];
  x[3]  = 190;   y[3]  =  656;   yerr[3]  = 0.047 * y[3];
  x[4]  = 210;   y[4]  =  433;   yerr[4]  = 0.048 * y[4];
  x[5]  = 230;   y[5]  =  293;   yerr[5]  = 0.051 * y[5];
  x[6]  = 250;   y[6]  =  205;   yerr[6]  = 0.047 * y[6];
  x[7]  = 270;   y[7]  =  146;   yerr[7]  = 0.048 * y[7];
  x[8]  = 290;   y[8]  =  105;   yerr[8]  = 0.049 * y[8];
  x[9]  = 310;   y[9]  =   77;   yerr[9]  = 0.048 * y[9];
  x[10] = 330;   y[10] =   57;   yerr[10] = 0.053 * y[10];
  x[11] = 350;   y[11] =   43;   yerr[11] = 0.055 * y[11];   
  x[12] = 370;   y[12] =   33;   yerr[12] = 0.057 * y[12];   
  x[13] = 390;   y[13] =   25;   yerr[13] = 0.057 * y[13];   
  x[14] = 410;   y[14] =   20;   yerr[14] = 0.060 * y[14];   

  for( int i = 0 ; i < 15; ++i ){
    xerr[i] = 0.0;
    yup[i]  = y[i] + yerr[i];
    ydn[i]  = y[i] - yerr[i];
  }

  for( int i = 0 ; i < 15; ++i ){
    xband[i] = x[i];
    yband[i] = y[i] + yerr[i];
  }

  for( int i = 0 ; i < 15; ++i ){
    xband[i+15] = x[14-i];
    yband[i+15] = y[14-i] - yerr[14-i];
  }
  
  // cout << endl << endl;
  // for( int i = 0 ; i < 30 ; ++i ){
  //   cout << xband[i] << " " << yband[i] << endl;
  // }
  // cout << endl << endl;

  TGraph* g     = new TGraph(n,x,y);
  TGraph* gup   = new TGraph(n,x,yup);
  TGraph* gdn   = new TGraph(n,x,ydn);
  TGraph* gband = new TGraph(30,xband,yband);


  // UP:   248
  // DOWN: 148

  //TGraphErrors* g  = new TGraphErrors(n,x,y,xerr,yerr);

  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  //gPad->SetGridx();
  //gPad->SetGridy();

  float ymin = 0;
  if( logplot ) ymin = 99.9;

  TH2F* hdummy = new TH2F("hdummy","",100,130,300,100,ymin,3000);
  hdummy->Draw();

  c1->cd();
  if( logplot ) gPad->SetLogy();

  g->SetLineColor(2);
  g->SetLineWidth(3);
  g->SetFillColor(5);
  gup->SetLineColor(2);
  gdn->SetLineColor(2);
  gup->SetLineStyle(2);
  gdn->SetLineStyle(2);
  gband->SetFillColor(5);

  //2l2j observed
  gul->SetLineColor(6);
  gul->SetLineWidth(3);
  gul->SetLineStyle(4);

  //2l2j expected
  gulexp->SetLineColor(2);
  gulexp->SetLineWidth(3);
  gulexp->SetLineStyle(2);

  //4l observed
  gul2->SetLineWidth(3);
  gul2->SetLineStyle(4);
  gul2->SetLineColor(kGreen+2);

  //4l expected
  gul2exp->SetLineWidth(3);
  gul2exp->SetLineStyle(2);

  //combined observed
  gulc->SetLineWidth(5);
  gulc->SetLineColor(1);

  //combined expected
  gulcexp->SetLineWidth(5);
  gulcexp->SetLineColor(4);
  gulcexp->SetLineStyle(2);






  // gulexp->SetLineWidth(2);
  // gulexp->SetLineStyle(2);
  // gulexp->SetLineColor(2);  
  // gul2exp->SetLineWidth(2);
  // gul2exp->SetLineStyle(2);


  hdummy->GetXaxis()->SetTitle("#mu [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma [fb]");
  hdummy->GetYaxis()->SetLabelSize(0.04);
  hdummy->GetXaxis()->SetLabelSize(0.04);
  hdummy->GetYaxis()->SetTitleSize(0.05);
  hdummy->GetXaxis()->SetTitleSize(0.05);
  hdummy->GetXaxis()->SetTitleOffset(1.12);
  hdummy->GetYaxis()->SetTitleOffset(1.5);

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

  gband->Draw("samef");
  g->Draw("samel");

  gulcband->SetFillStyle(3002);
  gulcband->Draw("samef");


  if( plotObserved ){
    gul->Draw("samel");
    gul2->Draw("samel");
    gulc->Draw("samel");
    gulcexp->Draw("samel");
    //gulcexpp1->Draw("samel");
    //gulcexpm1->Draw("samel");
  }

  if( plotExpected ){
    gulexp->Draw("samel");
    gul2exp->Draw("samel");
    gulcexp->Draw("samel");
  }


  //gband->Draw("samef");
  //g->Draw("samel");
  //gup->Draw("samel");
  //gdn->Draw("samel");

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

  //TLegend *leg = new TLegend(0.4,0.6,0.9,0.8);

  TH1F* hgexp = new TH1F("hgexp","",1,0,1);
  hgexp->SetLineColor(4);
  hgexp->SetLineWidth(5);
  hgexp->SetLineStyle(2);
  hgexp->SetFillColor(7);
  hgexp->SetFillStyle(3002);



  TLegend *leg = new TLegend(0.33,0.7,0.9,0.88);
  if( plotObserved ){
    leg->AddEntry(gulc    ,"combined observed UL","l");
    //leg->AddEntry(gulcexp ,"combined median expected UL","l");
    leg->AddEntry(hgexp   ,"combined median expected UL (#pm1#sigma)","lf");
    leg->AddEntry(gul     ,"2#font[12]{l}2j observed UL","l");
    leg->AddEntry(gul2    ,"4#font[12]{l} observed UL","l");

  }
  if( plotExpected ){
    leg->AddEntry(gulexp  ,"expected UL (VZ+E_{T}^{miss})","l");
    leg->AddEntry(gul2exp ,"expected UL (multi-lepton)","l");
    leg->AddEntry(gulcexp ,"expected UL (combined)","l");
  }

  TH1F* hg = new TH1F("h","",1,0,1);
  hg->SetLineColor(2);
  hg->SetLineWidth(3);
  hg->SetFillColor(5);

  //leg->AddEntry(g,  "theory","l");
  leg->AddEntry(hg,  "#sigma^{NLO} theory (#pm1#sigma)","lf");

  //leg->AddEntry(box,"excluded region","f");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();								
  t->SetTextSize(0.04);
  //t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");
  //t->DrawLatex(0.18,0.93,"CMS Preliminary,  #sqrt{s}=7 TeV,  L_{int}=4.98 fb^{-1}");
  cmsPrelim(4.98,isPreliminary);
  t->SetTextSize(0.04);
  //t->DrawLatex(0.47,0.45,"");
  t->DrawLatex(0.57,0.63,"GMSB  ZZ + E_{T}^{miss}");

  if( printplots ){
    if( isPreliminary) c1->Print("GMSB_Fig12_prelim.pdf");
    else               c1->Print("GMSB_Fig12.pdf");

    // c1->Print("GMSB.png");
    // c1->Print("GMSB.eps");
    // gROOT->ProcessLine(".! ps2pdf GMSB.eps GMSB_ppt.pdf");
  }	   
}
