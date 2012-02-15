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

void makeMETPlots(){

  gStyle->SetOptFit(0);

  TChain *ch = new TChain("t");
  //ch->Add("../output/V00-02-09/highpt/LM6v2_smallTree.root");
  //ch->Add("../output/V00-02-00/LM4_baby.root");
  //ch->Add("../output/V00-02-16/highpt/LM6v2_smallTree_gen_TEMP.root");
  ch->Add("../output/V00-02-18/highpt/LM6v2_smallTree_gen.root");

  vector<TCut> metcuts;
  vector<float> metcutvals;

  metcuts.push_back(TCut("pfmet>200")); metcutvals.push_back(200);
  metcuts.push_back(TCut("pfmet>275")); metcutvals.push_back(275);

  //TCut sel("njets>=2 && ht>100 && !passz");
  TCut sel("foundPair==1 && reco1==1 && reco2==1 && ht>100 && htgen2>100 && ngenjets>=2");

  const unsigned int n = metcuts.size();

  TH1F* hpass[n];
  TH1F* hall[n];

  for( unsigned int i = 0 ; i < metcuts.size() ; ++i){

    hpass[i]   = new TH1F(Form("hpass_%i",i),Form("hpass_%i",i),25,0,500);
    hall[i]    = new TH1F(Form("hall_%i",i), Form("hall_%i",i) ,25,0,500);

    ch->Draw(Form("genmet>>hpass_%i",i),sel+metcuts.at(i));
    ch->Draw(Form("genmet>>hall_%i",i)  ,sel);

  }




  TCanvas *can = new TCanvas();
  can->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetTopMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();  
  gStyle->SetOptFit(0);

  TGraphAsymmErrors* gr[n];  
  TLegend *leg = new TLegend(0.6,0.2,0.87,0.4);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);

  TF1* erf[n];

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);

  for( unsigned int i = 0 ; i < metcuts.size() ; ++i){

    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),500,500);
    //can[i]->cd();
    
    // TF1* efunc = new TF1("efitf", fitf, 0, 500, 3);
    // efunc->SetParameters(1, 100, 10);
    // efunc->SetParNames("norm", "offset", "width");

    erf[i] = new TF1("efitf", fitf, 0, 500, 3);
    erf[i]->SetParameters(1, 100, 30);
    erf[i]->SetParNames("norm", "offset", "width");
    erf[i]->SetLineWidth(2);
    //erf[i]->FixParameter(0,1);

    //erf[i] = new TF1(Form("erf_%i",i),mfitf,0,400);

    //erf[i]->SetParameter(0,100*(i+1));
    //erf[i]->SetParameter(1,10);

    gr[i] = new TGraphAsymmErrors();
    if( i==0 ){
      erf[i]->SetLineColor(1);
      line.SetLineColor(1);
    }

    if( i==1 ){
      line.SetLineColor(2);
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

    gr[i]->GetXaxis()->SetRangeUser(0,500);
    gr[i]->GetXaxis()->SetTitle("generator-level E_{T}^{miss} (GeV)");
    gr[i]->GetYaxis()->SetTitle("efficiency");
    gr[i]->SetMaximum(1);
    gr[i]->BayesDivide(hpass[i],hall[i]);

    //gr[i]->Fit(efunc,"R");
    gr[i]->Fit(erf[i],"R");


    if( i==0 ) gr[i]->Draw("AP");
    else       gr[i]->Draw("sameP");

    gr[i]->GetXaxis()->SetRangeUser(0,500);
    gr[i]->GetXaxis()->SetTitle("generator-level E_{T}^{miss} (GeV)");
    gr[i]->GetYaxis()->SetTitle("efficiency");

    line.DrawLine(metcutvals.at(i),0,metcutvals.at(i),1);
    //erf[i]->Draw("same");
  }

  leg->Draw();



  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.05);
  t->DrawLatex(0.25,0.92,"CMS Simulation, #sqrt{s} = 7 TeV");


  can->Print("../plots/met_turnon_LM6.pdf");





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
if (par[2] != 0)
arg = (x[0] - 100.)/par[2];

double fitval = par[0]*TMath::Erf(arg)+par[1]*(1.-TMath::Erf(arg));
return fitval;
}
