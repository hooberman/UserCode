#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCut.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include <iomanip>

using namespace std;


void format(TH1F *h1, bool data)
{

    h1->SetLineWidth(2);
    if (data) {
        h1->SetLineColor(kRed);
        h1->SetMarkerColor(kRed);
    }
    else {
        h1->SetLineColor(kBlue);
        h1->SetMarkerColor(kBlue);
        h1->SetFillStyle(3003);
        h1->SetFillColor(kBlue);
    }

}

void format(TGraph *h1, bool data)
{

    h1->SetLineWidth(2);
    if (data) {
        h1->SetLineColor(kRed);
        h1->SetMarkerColor(kRed);
    }
    else {
        h1->SetLineColor(kBlue);
        h1->SetMarkerColor(kBlue);
        h1->SetMarkerStyle(4);
    }

}

void plotSimple_OS(TString data, TString mc, TCut var1, TCut var2, Int_t nBinsVar, 
		   Float_t minVar, Float_t maxVar, Int_t type = 0, bool iso = false, bool fastjet = false ){

  char* leptype = "";
  float pt2cut  = -1;
  if     ( type == 0 ){
    leptype = "ee";
    pt2cut = 10.;
  }
  else if( type == 1 ){
    leptype = "mm";
    pt2cut = 5.;
  }
  else{
    cout << "Unrecognized leptype " << type << endl;
    exit(0);
  }

  char* efftype = "";
  if( iso ) efftype = "iso";
  else      efftype = "ID";

  char* fastjetchar = "";
  if( fastjet ) fastjetchar = "F";

  cout << "leptype  " << leptype     << endl;
  cout << "efftype  " << efftype     << endl;
  cout << "fastjet  " << fastjetchar << endl;

  TFile f_mc(mc, "READ");
  TTree *t_mc = (TTree*)f_mc.Get("tree");

  TFile f_data(data, "READ");
  TTree *t_data = (TTree*)f_data.Get("tree");

  //gROOT->cd();
  gStyle->SetOptTitle(1);

  TH1F *h1_mc_var_2TT = new TH1F("h1_mc_var_2TT", Form("%s (2TT); %s", var1.GetTitle(), var1.GetTitle()), nBinsVar, minVar, maxVar);
  TH1F *h1_mc_var_TP  = new TH1F("h1_mc_var_TP",  Form("%s (TP); %s", var1.GetTitle(), var1.GetTitle()),  nBinsVar, minVar, maxVar);
  TH1F *h1_mc_var_TF  = new TH1F("h1_mc_var_TF",  Form("%s (TF); %s", var1.GetTitle(), var1.GetTitle()),  nBinsVar, minVar, maxVar);
  TH1F *h1_data_var_2TT = new TH1F("h1_data_var_2TT", Form("%s (2TT); %s", var1.GetTitle(), var1.GetTitle()), nBinsVar, minVar, maxVar);
  TH1F *h1_data_var_TP  = new TH1F("h1_data_var_TP",  Form("%s (TP); %s", var1.GetTitle(), var1.GetTitle()),  nBinsVar, minVar, maxVar);
  TH1F *h1_data_var_TF  = new TH1F("h1_data_var_TF",  Form("%s (TF); %s", var1.GetTitle(), var1.GetTitle()),  nBinsVar, minVar, maxVar);
  format(h1_mc_var_2TT, false);
  format(h1_mc_var_TP, false);
  format(h1_mc_var_TF, false);
  format(h1_data_var_2TT, true);
  format(h1_data_var_TP, true);
  format(h1_data_var_TF, true);

  Int_t nBinsMass = 34;
  Float_t minMass = 30.0;
  Float_t maxMass = 200.0;
  TH1F *h1_mc_mass_2TT = new TH1F("h1_mc_mass_2TT", "Mass (2TT); Mass", nBinsMass, minMass, maxMass);
  TH1F *h1_mc_mass_TP  = new TH1F("h1_mc_mass_TP",  "Mass (TP); Mass",  nBinsMass, minMass, maxMass);
  TH1F *h1_mc_mass_TF  = new TH1F("h1_mc_mass_TF",  "Mass (TF); Mass",  nBinsMass, minMass, maxMass);
  TH1F *h1_data_mass_2TT = new TH1F("h1_data_mass_2TT", "Mass (2TT); Mass", nBinsMass, minMass, maxMass);
  TH1F *h1_data_mass_TP  = new TH1F("h1_data_mass_TP",  "Mass (TP); Mass",  nBinsMass, minMass, maxMass);
  TH1F *h1_data_mass_TF  = new TH1F("h1_data_mass_TF",  "Mass (TF); Mass",  nBinsMass, minMass, maxMass);
  format(h1_mc_mass_2TT, false);
  format(h1_mc_mass_TP, false);
  format(h1_mc_mass_TF, false);
  format(h1_data_mass_2TT, true);
  format(h1_data_mass_TP, true);
  format(h1_data_mass_TF, true);

  // di-electron event type
  // probe leg must pass di-lepton trigger
  // - in effect this means the "loose" leg, which is SC8
  TCut cut_base("cut_base", Form("type == %i && passTrigger && pt1>20" , type ) );
  //TCut cut_mass("cut_mass", "mass > 76 && mass < 106");
  TCut cut_mass("cut_mass", "mass > 86 && mass < 96");
  
  // what cut to apply on the probe,
  // remember - it MUST be a subset of the tag
  //TCut cut_probe("cut_probe", "passId == 1 && pt2>10");

  //TCut cut_pass("cut_pass", "passIso == 1");
  //TCut cut_pass("cut_pass", "passIsoF == 1");
  //TCut cut_pass("cut_pass", Form("reliso2%s < %f",fastjetchar,relisocut));

  TCut cut_probe;
  TCut cut_pass;

  TCut cut_tightiso("pt2>10 || relisont2 <0.15");
    
  if( iso ){
    cut_probe = TCut("cut_probe" , Form("passId  == 1 && pt2>%.0f",pt2cut));
    if( fastjet ) cut_pass  = TCut("cut_pass"  , "passIsoF == 1");
    else          cut_pass  = TCut("cut_pass"  , "passIso == 1")  + cut_tightiso;
  }else{
    if( fastjet ) cut_probe = TCut("cut_probe" , Form("passIsoF == 1 && pt2>%.0f" , pt2cut));
    else          cut_probe = TCut("cut_probe" , Form("passIso == 1 && pt2>%.0f"  , pt2cut));
    cut_pass  = TCut("cut_pass"  , "passId  == 1");
  }

  cout << "Base selection  " << cut_base.GetTitle()  << endl;
  cout << "Probe selection " << cut_probe.GetTitle() << endl;
  cout << "Pass selection  " << cut_pass.GetTitle()  << endl;

  //
  // now do drawing
  //

  // remember, if either leg could be a tag
  // then both legs are unbiassed and can be used
  char *drawcommand1 = Form("%s>>", var1.GetTitle());

  t_mc->Draw(Form("%s>>h1_mc_var_2TT", var1.GetTitle()), cut_probe+cut_base+cut_mass+"tt==1", "goff");
  t_mc->Draw(Form("%s>>h1_mc_var_2TT", var2.GetTitle()), cut_probe+cut_base+cut_mass+"tt==1", "goff");
  t_mc->Draw("mass >> h1_mc_mass_2TT", cut_probe+cut_base+"tt==1", "goff");
  t_mc->Draw("mass >>+ h1_mc_mass_2TT", cut_probe+cut_base+"tt==1", "goff");

  t_data->Draw(Form("%s>>h1_data_var_2TT", var1.GetTitle()), cut_probe+cut_base+cut_mass+"tt==1", "goff");
  t_data->Draw(Form("%s>>h1_data_var_2TT", var2.GetTitle()), cut_probe+cut_base+cut_mass+"tt==1", "goff");
  t_data->Draw("mass >> h1_data_mass_2TT", cut_probe+cut_base+"tt==1", "goff");
  t_data->Draw("mass >>+ h1_data_mass_2TT", cut_probe+cut_base+"tt==1", "goff");

  // tag and passing probe
  t_mc->Draw(Form("%s>>h1_mc_var_TP", var2.GetTitle()), cut_probe+cut_base+cut_mass+"tt!=1"+cut_pass, "goff");
  t_mc->Draw("mass >> h1_mc_mass_TP", cut_probe+cut_base+"tt!=1"+cut_pass, "goff");

  t_data->Draw(Form("%s>>h1_data_var_TP", var2.GetTitle()), cut_probe+cut_base+cut_mass+"tt!=1"+cut_pass, "goff");
  t_data->Draw("mass >> h1_data_mass_TP", cut_probe+cut_base+"tt!=1"+cut_pass, "goff");

  // tag and failing probe
  t_mc->Draw(Form("%s>>h1_mc_var_TF", var2.GetTitle()), cut_probe+cut_base+cut_mass+"tt!=1"+!cut_pass, "goff");
  t_mc->Draw("mass >> h1_mc_mass_TF", cut_probe+cut_base+"tt!=1"+!cut_pass, "goff");

  t_data->Draw(Form("%s>>h1_data_var_TF", var2.GetTitle()), cut_probe+cut_base+cut_mass+"tt!=1"+!cut_pass, "goff");
  t_data->Draw("mass >> h1_data_mass_TF", cut_probe+cut_base+"tt!=1"+!cut_pass, "goff");

  //
  // now compute the efficiency
  //

  TH1F *h1_mc_numerator = (TH1F*)h1_mc_var_2TT->Clone("h1_mc_numerator");
  h1_mc_numerator->Add(h1_mc_var_TP);
  TH1F *h1_mc_denominator = (TH1F*)h1_mc_var_2TT->Clone("h1_mc_denominator");
  h1_mc_denominator->Add(h1_mc_var_TP);
  h1_mc_denominator->Add(h1_mc_var_TF);
  TGraphAsymmErrors *gr_mc_eff = new TGraphAsymmErrors();
  gr_mc_eff->Divide(h1_mc_numerator, h1_mc_denominator);
  gr_mc_eff->GetXaxis()->SetTitle(h1_mc_var_2TT->GetXaxis()->GetTitle());
  gr_mc_eff->SetTitle("Efficiency");

  TH1F *h1_data_numerator = (TH1F*)h1_data_var_2TT->Clone("h1_data_numerator");
  h1_data_numerator->Add(h1_data_var_TP);
  TH1F *h1_data_denominator = (TH1F*)h1_data_var_2TT->Clone("h1_data_denominator");
  h1_data_denominator->Add(h1_data_var_TP);
  h1_data_denominator->Add(h1_data_var_TF);
  TGraphAsymmErrors *gr_data_eff = new TGraphAsymmErrors();
  gr_data_eff->Divide(h1_data_numerator, h1_data_denominator);
  gr_data_eff->GetXaxis()->SetTitle(h1_data_var_2TT->GetXaxis()->GetTitle());
  gr_data_eff->SetTitle("Efficiency");

  format(gr_data_eff, true);
  format(gr_mc_eff, false);

  //
  // scale MC to data
  // where scale factor from number of 2TT 
  // in Z window in data and MC
  // NOTE: this is just cosmetic and
  // --MUST-- be done after computing the efficiency
  //

  Float_t nGoodTPMC   = h1_mc_var_2TT->Integral(0, nBinsVar+1);
  Float_t nGoodTPData = h1_data_var_2TT->Integral(0, nBinsVar+1);
  h1_mc_var_2TT->Scale(nGoodTPData/nGoodTPMC);
  h1_mc_var_TP->Scale(nGoodTPData/nGoodTPMC);
  h1_mc_var_TF->Scale(nGoodTPData/nGoodTPMC);
  h1_mc_mass_2TT->Scale(nGoodTPData/nGoodTPMC);
  h1_mc_mass_TP->Scale(nGoodTPData/nGoodTPMC);
  h1_mc_mass_TF->Scale(nGoodTPData/nGoodTPMC);


  //
  // and finally display it
  //

  TLegend *l1 = new TLegend(0.4, 0.2, 0.7, 0.4);
  l1->SetFillColor(kWhite);
  //l1->SetLineColor(kWhite);
  l1->SetShadowColor(kWhite);
  if( type == 0 ) l1->AddEntry(h1_mc_mass_2TT, "DY#rightarrow ee (Spring11)", "f");
  if( type == 1 ) l1->AddEntry(h1_mc_mass_2TT, "DY#rightarrow #mu#mu (Spring11)", "f");
  l1->AddEntry(h1_data_mass_2TT, "DATA (2011A)", "lp");

  /*
  // the variable
  TCanvas *c1 = new TCanvas();
  c1->Divide(2, 2);
  c1->cd(1);
  h1_data_var_2TT->Draw("E1");
  h1_mc_var_2TT->Draw("SAME HIST");
  l1->Draw();
  c1->cd(2);
  h1_data_var_TP->Draw("E1");
  h1_mc_var_TP->Draw("SAME HIST");
  l1->Draw();
  c1->cd(3);
  h1_data_var_TF->Draw("E1");
  h1_mc_var_TF->Draw("SAME HIST");
  l1->Draw();
  c1->cd(4);
  gr_data_eff->Draw("AP");
  gr_mc_eff->Draw("P");
  c1->RedrawAxis();

  // mass, as diagnostic
  TCanvas *c2 = new TCanvas();
  c2->Divide(2, 2);
  c2->cd(1);
  h1_data_mass_2TT->Draw("E1");
  h1_mc_mass_2TT->Draw("SAME HIST");
  l1->Draw();
  c2->cd(2);
  h1_data_mass_TP->Draw("E1");
  h1_mc_mass_TP->Draw("SAME HIST");
  l1->Draw();
  c2->cd(3);
  h1_data_mass_TF->Draw("E1");
  h1_mc_mass_TF->Draw("SAME HIST");
  l1->Draw();
  c2->RedrawAxis();
  */

  //float ymin = 0.5;
  //if( fabs( relisocut - 0.05 ) < 1e-5 ) ymin = 0.4;
  //if( fabs( relisocut - 0.10 ) < 1e-5 ) ymin = 0.6;
  //if( fabs( relisocut - 0.15 ) < 1e-5 ) ymin = 0.8;

  TCanvas *c3 = new TCanvas();
  c3->cd();
  gPad->SetGridy();
  gr_data_eff->GetYaxis()->SetRangeUser(0.9,1.);
  if( type == 0 ) gr_data_eff->GetXaxis()->SetTitle("electron p_{T} (GeV)");
  if( type == 1 ) gr_data_eff->GetXaxis()->SetTitle("muon p_{T} (GeV)");
  //gr_data_eff->GetXaxis()->SetTitle("N_{PV}");
  gr_data_eff->SetTitle("");
  if( iso ) gr_data_eff->GetYaxis()->SetTitle("Iso Eff");
  else      gr_data_eff->GetYaxis()->SetTitle("ID Eff");
  gr_data_eff->Draw("AP");
  gr_mc_eff->Draw("P");
  l1->Draw();
  c3->RedrawAxis();

  c3->Print(Form("plots/%s_%s%s.gif",efftype,leptype,fastjetchar));

}

