#include <algorithm>
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
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>


using namespace std;

double mfitf (double* x, double* par);
double fitf (double* x, double* par);

void makeMETPlots( bool printplot = false ){

  gStyle->SetOptFit(0);

  TChain *ch = new TChain("T1");
  ch->Add("../../output/V00-02-21/wz_summer11_madgraph_gen_baby.root");

  vector<TCut> metcuts;
  vector<float> metcutvals;

  metcuts.push_back(TCut("pfmet>100")); metcutvals.push_back(100);
  metcuts.push_back(TCut("pfmet>200")); metcutvals.push_back(200);
  metcuts.push_back(TCut("pfmet>300")); metcutvals.push_back(300);

  TCut sel("dilmass>81&&dilmass<101&&njets>=2");

  const unsigned int n = metcuts.size();

  TH1F* hpass[n];
  TH1F* hall[n];

  for( unsigned int i = 0 ; i < metcuts.size() ; ++i){

    hpass[i]   = new TH1F(Form("hpass_%i",i),Form("hpass_%i",i),30,0,600);
    hall[i]    = new TH1F(Form("hall_%i",i), Form("hall_%i",i) ,30,0,600);

    ch->Draw(Form("genmet>>hpass_%i",i),sel+metcuts.at(i));
    ch->Draw(Form("genmet>>hall_%i",i)  ,sel);

  }

  TCanvas *can = new TCanvas();
  can->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTopMargin(0.08);

  TGraphAsymmErrors* gr[n];  
  TLegend *leg = new TLegend(0.6,0.2,0.95,0.4);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.035);

  TF1* erf[n];

  for( unsigned int i = 0 ; i < metcuts.size() ; ++i){

    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),600,600);
    //can[i]->cd();
    
    TF1* efunc = new TF1("efitf", fitf, 0, 600, 3);
    efunc->SetParameters(1, 100, 10);
    efunc->SetParNames("norm", "offset", "width");

    erf[i] = new TF1("efitf", fitf, 0, 600, 3);
    erf[i]->SetParameters(1, 100, 10);
    erf[i]->SetParNames("norm", "offset", "width");
    erf[i]->SetLineWidth(2);
    //erf[i]->FixParameter(0,1);

    //erf[i] = new TF1(Form("erf_%i",i),mfitf,0,400);

    //erf[i]->SetParameter(0,100*(i+1));
    //erf[i]->SetParameter(1,10);

    gr[i] = new TGraphAsymmErrors();
    if( i==0 ){
      erf[i]->SetLineColor(1);
    }

    if( i==1 ){
      gr[i]->SetLineColor(2);
      gr[i]->SetMarkerColor(2);
      gr[i]->SetMarkerStyle(21);
      erf[i]->SetLineColor(2);
    }
    if( i==2 ){
      gr[i]->SetLineColor(4);
      gr[i]->SetMarkerColor(4);
      gr[i]->SetMarkerStyle(25);
      erf[i]->SetLineColor(4);
    }

    leg->AddEntry(gr[i],Form("E_{T}^{miss}>%.0f GeV",metcutvals.at(i)),"p");

    gr[i]->GetXaxis()->SetTitle("generator-level E_{T}^{miss} (GeV)");
    gr[i]->GetYaxis()->SetTitle("efficiency");
    gr[i]->SetMaximum(1.05);
    gr[i]->BayesDivide(hpass[i],hall[i]);

    //gr[i]->Fit(efunc,"R");
    gr[i]->Fit(erf[i],"R");


    if( i==0 ) gr[i]->Draw("AP");
    else       gr[i]->Draw("sameP");

    gr[i]->GetXaxis()->SetTitle("generator E_{T}^{miss} [GeV]");
    gr[i]->GetYaxis()->SetTitle("efficiency");

    //erf[i]->Draw("same");
  }

  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);
  t->DrawLatex(0.28,0.95,"CMS Simulation,  #sqrt{s} = 7 TeV");


  if (printplot) can->Print("../plots/met_turnon_LM4.pdf");





}

double fitf (double* x, double* par) {
double arg = 0;
if (par[2] != 0)
arg = (x[0] - par[1])/par[2];

double fitval = 0.5 * par[0] * (TMath::Erf(arg) + 1);
return fitval;
}

double mfitf (double* x, double* par) {
double arg = 0;

 double ptcutoff = 10;
 //double ptcutoff = 20;

if (par[2] != 0)
arg = (x[0] - ptcutoff)/par[2];  

double fitval = par[0]*TMath::Erf(arg)+par[1]*(1.-TMath::Erf(arg));
return fitval;
}
