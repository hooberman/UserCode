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




double mfitf (double* x, double* par) {

  double arg = 0;
  if (par[2] != 0)
    arg = (x[0] - 10.)/par[2];
  
  double fitval = par[0]*TMath::Erf(arg)+par[1]*(1.-TMath::Erf(arg));
  return fitval;
}

void makeLeptonPlots(){

    gStyle->SetOptFit(0);

    TChain *ch = new TChain("t");
    //ch->Add("../output/V00-02-09/highpt/LM6v2_smallTree.root");
    //ch->Add("../output/V00-02-10/highpt/LM6v2_smallTree.root");
    //ch->Add("../output/V00-02-15/highpt/LM6v2_smallTree.root");
    ch->Add("../output/V00-02-18/highpt/LM6v2_smallTree_gen.root");
    
    char* var1 = "genlep1.pt()";
    char* var2 = "genlep2.pt()";
    int nbins =  28;
    int xmin  =  10;
    int xmax  = 150;

    TF1* fel = new TF1("fel", mfitf, xmin, xmax, 3);
    fel->SetParameters(0.9, 0.5, 20);
    fel->SetParNames("epsinf", "epsC", "width");

    TF1* fmu = new TF1("fmu", mfitf, xmin, xmax, 3);
    fmu->SetParameters(0.9, 0.5, 20);
    fmu->SetParNames("epsinf", "epsC", "width");


    // char* var1 = "htgen2";
    // char* var2 = "htgen2";
    // int nbins = 20;
    // int xmin  = 0;
    // int xmax  = 1000;

    // TF1* fel = new TF1("fel", "pol1", xmin, xmax);
    // TF1* fmu = new TF1("fmu", "pol1", xmin, xmax);


    //TCut sel("njets>=2 && pfmet>50 && !passz");
    TCut sel("foundPair==1 && genmet>200 && htgen2>300 ");
    //TCut sel("");
    TCut el1("abs(genid1)==11");
    TCut el2("abs(genid2)==11");
    TCut mu1("abs(genid1)==13");
    TCut mu2("abs(genid2)==13");
    TCut pass1("reco1==1");
    TCut pass2("reco2==1");
    
    const unsigned int n = 2;

    // TH1F* hpass1[2];
    // TH1F* hpass2[2];
    // TH1F* hpass[2];
    // TH1F* hall1[2];
    // TH1F* hall2[2];
    // TH1F* hall[2];
    
    //--------------------
    // electrons
    //--------------------
    
    TH1F* hpassel1 = new TH1F("hpassel1","",nbins,xmin,xmax);
    TH1F* hallel1  = new TH1F("hallel1" ,"",nbins,xmin,xmax);
    TH1F* hpassel2 = new TH1F("hpassel2","",nbins,xmin,xmax);
    TH1F* hallel2  = new TH1F("hallel2" ,"",nbins,xmin,xmax);

    ch->Draw(Form("%s>>hallel1" ,var1) , sel+el1);
    ch->Draw(Form("%s>>hpassel1",var1) , sel+el1+pass1);
    ch->Draw(Form("%s>>hallel2" ,var2) , sel+el2);
    ch->Draw(Form("%s>>hpassel2",var2) , sel+el2+pass2);
    
    TH1F* hpassel = (TH1F*) hpassel1->Clone("hpassel");
    TH1F* hallel  = (TH1F*) hallel1->Clone("hallel");
    hpassel->Add(hpassel2);
    hallel->Add(hallel2);

    //--------------------
    // muons
    //--------------------
    
    TH1F* hpassmu1 = new TH1F("hpassmu1","",nbins,xmin,xmax);
    TH1F* hallmu1  = new TH1F("hallmu1" ,"",nbins,xmin,xmax);
    TH1F* hpassmu2 = new TH1F("hpassmu2","",nbins,xmin,xmax);
    TH1F* hallmu2  = new TH1F("hallmu2" ,"",nbins,xmin,xmax);

    ch->Draw(Form("%s>>hallmu1" ,var1) , sel+mu1);
    ch->Draw(Form("%s>>hpassmu1",var1) , sel+mu1+pass1);
    ch->Draw(Form("%s>>hallmu2" ,var2) , sel+mu2);
    ch->Draw(Form("%s>>hpassmu2",var2) , sel+mu2+pass2);
    
    TH1F* hpassmu = (TH1F*) hpassmu1->Clone("hpassmu");
    TH1F* hallmu  = (TH1F*) hallmu1->Clone("hallmu");
    hpassmu->Add(hpassmu2);
    hallmu->Add(hallmu2);


    TCanvas *c1 = new TCanvas();
    c1->cd();

    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
    gPad->SetGridx();
    gPad->SetGridy();

    fel->SetLineWidth(2);
    fmu->SetLineWidth(2);
    fmu->SetLineColor(4);

    TGraphAsymmErrors *grel = new TGraphAsymmErrors();
    grel->BayesDivide(hpassel,hallel);    
    grel->Fit(fel,"R");

    TGraphAsymmErrors *grmu = new TGraphAsymmErrors();
    grmu->BayesDivide(hpassmu,hallmu);
    grmu->Fit(fmu,"R");

    grel->SetMarkerColor(2);
    grel->SetLineColor(2);
    grel->SetMarkerStyle(21);
    grmu->SetMarkerColor(4);
    grmu->SetLineColor(4);
    grmu->SetMarkerStyle(25);

    grel->GetXaxis()->SetTitle("generated lepton p_{T} (GeV)");
    grel->GetYaxis()->SetTitle("efficiency");
    grel->SetMinimum(0.4);
    grel->SetMaximum(1.0);

    grel->GetXaxis()->SetRangeUser(10,150);
    grel->Draw("AP");
    grmu->Draw("sameP");

    TLegend *leg = new TLegend(0.5,0.2,0.8,0.4);
    leg->AddEntry(grel,"electrons","lp");
    leg->AddEntry(grmu,"muons","lp");
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextSize(0.05);
    t->DrawLatex(0.25,0.92,"CMS Simulation, #sqrt{s} = 7 TeV");

    c1->Print("makeLeptonPlots.pdf");

}
