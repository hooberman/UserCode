#include <algorithm>
#include <iostream>
#include <iomanip>
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
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include <sstream>

double fitf (double* x, double* par) {
  double arg = x[0];

  //if( par[1] == 0 ) cout << "Error!
  double fitval = par[6] + par[4] * ( TMath::Erf((arg-par[0])/par[1]) + 1 ) + par[5] * ( TMath::Erf((arg-par[2])/par[3]) + 1 );
  return fitval;
}



void jetEfficiencies( bool printplots = false ){

  gStyle->SetOptFit(0);

  TChain *ch = new TChain("T1");
  //ch->Add("../output/V00-02-20/wzsms_gen_baby.root");
  //TFile* f = TFile::Open("../output/V00-02-20/wzsms.root");

  //ch->Add("../output/V00-02-20/wz_summer11_madgraph_gen_baby.root");
  //TFile* f = TFile::Open("../output/V00-02-20/wz_summer11_madgraph.root");

  ch->Add("../output/V00-02-20/ttbar_gen_baby.root");
  TFile* f = TFile::Open("../output/V00-02-20/ttbar.root");

  //---------------------------------------
  // jets
  //---------------------------------------

  TH1F* hjetpt_all      = (TH1F*) f->Get("hjetpt_all");
  TH1F* hjetpt_q_all    = (TH1F*) f->Get("hjetpt_q_all");
  TH1F* hjetpt_c_all    = (TH1F*) f->Get("hjetpt_c_all");
  TH1F* hjetpt_b_all    = (TH1F*) f->Get("hjetpt_b_all");
  TH1F* hjetpt_g_all    = (TH1F*) f->Get("hjetpt_g_all");
  TH1F* hjetpt_pass04   = (TH1F*) f->Get("hjetpt_pass04");
  TH1F* hbtag_q_all     = (TH1F*) f->Get("hbtag_q_all");
  TH1F* hbtag_q_pass    = (TH1F*) f->Get("hbtag_q_pass");
  TH1F* hbtag_c_all     = (TH1F*) f->Get("hbtag_c_all");
  TH1F* hbtag_c_pass    = (TH1F*) f->Get("hbtag_c_pass");
  TH1F* hbtag_b_all     = (TH1F*) f->Get("hbtag_b_all");
  TH1F* hbtag_b_pass    = (TH1F*) f->Get("hbtag_b_pass");
  TH1F* hbtag_g_all     = (TH1F*) f->Get("hbtag_g_all");
  TH1F* hbtag_g_pass    = (TH1F*) f->Get("hbtag_g_pass");

  int rebin = 10;

  hjetpt_all->Rebin(rebin);
  hjetpt_q_all->Rebin(rebin);
  hjetpt_b_all->Rebin(rebin);
  hjetpt_c_all->Rebin(rebin);
  hjetpt_g_all->Rebin(rebin);
  hjetpt_pass04->Rebin(rebin);
  hbtag_q_all->Rebin(rebin);
  hbtag_q_pass->Rebin(rebin);
  hbtag_c_all->Rebin(rebin);
  hbtag_c_pass->Rebin(rebin);
  hbtag_b_all->Rebin(rebin);
  hbtag_b_pass->Rebin(rebin);
  hbtag_g_all->Rebin(rebin);
  hbtag_g_pass->Rebin(rebin);

  cout << "Total entries   " << hjetpt_all->GetEntries()   << endl;
  cout << "Total q entries " << hjetpt_q_all->GetEntries() << endl;
  cout << "Total c entries " << hjetpt_b_all->GetEntries() << endl;
  cout << "Total b entries " << hjetpt_c_all->GetEntries() << endl;
  cout << "Total g entries " << hjetpt_g_all->GetEntries() << endl;

  cout << "Total jet entries   " << hjetpt_pass04->GetEntries()   << endl;
  cout << "Total jet q entries " << hbtag_q_all->GetEntries() << endl;
  cout << "Total jet c entries " << hbtag_b_all->GetEntries() << endl;
  cout << "Total jet b entries " << hbtag_c_all->GetEntries() << endl;
  cout << "Total jet g entries " << hbtag_g_all->GetEntries() << endl;

  TGraphAsymmErrors* grjet = new TGraphAsymmErrors();
  grjet->BayesDivide(hjetpt_pass04,hjetpt_all);

  TGraphAsymmErrors* grjet_q = new TGraphAsymmErrors();
  grjet_q->BayesDivide(hbtag_q_all,hjetpt_q_all);

  TGraphAsymmErrors* grjet_c = new TGraphAsymmErrors();
  grjet_c->BayesDivide(hbtag_c_all,hjetpt_c_all);

  TGraphAsymmErrors* grjet_b = new TGraphAsymmErrors();
  grjet_b->BayesDivide(hbtag_b_all,hjetpt_b_all);

  TGraphAsymmErrors* grjet_g = new TGraphAsymmErrors();
  grjet_g->BayesDivide(hbtag_g_all,hjetpt_g_all);

  TGraphAsymmErrors* grbtag_q = new TGraphAsymmErrors();
  grbtag_q->BayesDivide(hbtag_q_pass,hbtag_q_all);

  TGraphAsymmErrors* grbtag_c = new TGraphAsymmErrors();
  grbtag_c->BayesDivide(hbtag_c_pass,hbtag_c_all);

  TGraphAsymmErrors* grbtag_b = new TGraphAsymmErrors();
  grbtag_b->BayesDivide(hbtag_b_pass,hbtag_b_all);

  TGraphAsymmErrors* grbtag_g = new TGraphAsymmErrors();
  grbtag_g->BayesDivide(hbtag_g_pass,hbtag_g_all);

  TGraphAsymmErrors* grbtagjet_q = new TGraphAsymmErrors();
  grbtagjet_q->BayesDivide(hbtag_q_pass,hjetpt_q_all);

  TGraphAsymmErrors* grbtagjet_c = new TGraphAsymmErrors();
  grbtagjet_c->BayesDivide(hbtag_c_pass,hjetpt_c_all);

  TGraphAsymmErrors* grbtagjet_b = new TGraphAsymmErrors();
  grbtagjet_b->BayesDivide(hbtag_b_pass,hjetpt_b_all);

  TGraphAsymmErrors* grbtagjet_g = new TGraphAsymmErrors();
  grbtagjet_g->BayesDivide(hbtag_g_pass,hjetpt_g_all);

  //---------------------------------------------
  // Nik canvas
  //---------------------------------------------

  TCanvas *c1 = new TCanvas();
  c1->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grjet->SetMarkerColor(4);
  grjet->SetLineColor(4);
  grbtagjet_b->SetMarkerColor(2);
  grbtagjet_b->SetLineColor(2);
  grbtagjet_c->SetMarkerColor(1);
  grbtagjet_c->SetLineColor(1);
  grbtagjet_q->SetMarkerColor(kGreen+2);
  grbtagjet_q->SetLineColor(kGreen+2);
  grbtagjet_g->SetMarkerColor(6);
  grbtagjet_g->SetLineColor(6);

  grjet->GetXaxis()->SetNdivisions(7);
  grjet->GetXaxis()->SetTitle("parton p_{T} [GeV]");
  grjet->GetYaxis()->SetTitle("efficiency");
  grjet->Draw("AP");
  grbtagjet_b->Draw("sameP");
  grbtagjet_c->Draw("sameP");

  if( printplots ) c1->Print("../plots/nik.pdf");

  //---------------------------------------------
  // jets
  //---------------------------------------------

  TCanvas *c2 = new TCanvas();
  c2->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grjet->SetMinimum(0);
  grjet->SetMaximum(1);
  grjet->GetXaxis()->SetTitle("generator parton p_{T} [GeV]");
  grjet->GetXaxis()->SetLabelSize(0.04);
  grjet->GetXaxis()->SetNdivisions(11);
  grjet->GetYaxis()->SetTitle("jet reco efficiency");

  grjet->SetMarkerColor(4);
  grjet->SetMarkerSize(1.5);
  grjet->SetLineWidth(2);
  grjet->SetLineColor(4);
  grjet_q->SetMarkerColor(2);
  grjet_q->SetLineColor(2);
  grjet_c->SetMarkerColor(1);
  grjet_c->SetLineColor(1);
  grjet_b->SetMarkerColor(kGreen+2);
  grjet_b->SetLineColor(kGreen+2);
  grjet_g->SetMarkerColor(6);
  grjet_g->SetLineColor(6);

  grjet->Draw("AP");
  grjet_q->Draw("sameP");
  grjet_c->Draw("sameP");
  grjet_b->Draw("sameP");
  grjet_g->Draw("sameP");

  TLegend *leg = new TLegend(0.5,0.3,0.7,0.6);
  leg->AddEntry(grjet,"all","lp");
  leg->AddEntry(grjet_q,"uds","lp");
  leg->AddEntry(grjet_c,"c","lp");
  leg->AddEntry(grjet_b,"b","lp");
  leg->AddEntry(grjet_g,"g","lp");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.DrawLine(20,0,20,1);

  if( printplots ) c2->Print("../plots/jetreco.pdf");


  //---------------------------------
  // b-tagging X reco jet
  //---------------------------------

  TCanvas *c3 = new TCanvas();
  c3->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grbtagjet_b->SetMinimum(0);
  grbtagjet_b->SetMaximum(1);

  grbtagjet_b->GetXaxis()->SetTitle("parton p_{T} [GeV]");
  grbtagjet_b->GetYaxis()->SetTitle("eff(reco jet #times b-tagging)");

  grbtagjet_b->Draw("AP");
  grbtagjet_c->Draw("sameP");
  grbtagjet_q->Draw("sameP");
  grbtagjet_g->Draw("sameP");

  TLegend *leg3 = new TLegend(0.4,0.4,0.6,0.6);
  leg3->AddEntry(grbtagjet_b,"b","lp");
  leg3->AddEntry(grbtagjet_c,"c","lp");
  leg3->AddEntry(grbtagjet_g,"g","lp");
  leg3->AddEntry(grbtagjet_q,"uds","lp");
  leg3->SetFillColor(0);
  leg3->SetBorderSize(1);
  leg3->Draw();

  if( printplots ) c3->Print("../plots/btag.pdf");

  //---------------------------------
  // b-tagging
  //---------------------------------

  TCanvas *c4 = new TCanvas();
  c4->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grbtag_q->SetMarkerColor(kGreen+2);
  grbtag_q->SetLineColor(kGreen+2);
  grbtag_c->SetMarkerColor(1);
  grbtag_c->SetLineColor(1);
  grbtag_b->SetMarkerColor(2);
  grbtag_b->SetLineColor(2);
  grbtag_g->SetMarkerColor(6);
  grbtag_g->SetLineColor(6);

  grbtag_q->GetXaxis()->SetRangeUser(0,200);
  grbtag_q->GetXaxis()->SetTitle("parton p_{T} [GeV]");
  grbtag_q->GetYaxis()->SetTitle("b-tagging efficiency");
  grbtag_q->SetMinimum(0.001);
  grbtag_q->SetMaximum(1.0);

  // TF1*  fit_c = new TF1("fit_c", fitf, 20, 200, 7);
  // fit_c->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  // fit_c->FixParameter(6,0);
  // fit_c->SetLineWidth(2);

  // TF1*  fit_b = new TF1("fit_c", fitf, 20, 200, 7);
  // fit_b->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  // fit_b->FixParameter(6,0);
  // fit_b->SetLineWidth(2);
  // fit_b->SetLineColor(1);

  // TF1*  fit_q = new TF1("fit_q", fitf, 20, 200, 7);
  // fit_q->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  // fit_q->FixParameter(6,0);
  // fit_q->SetLineWidth(2);
  // fit_q->SetLineColor(4);

  // grbtag_c->Fit(fit_c);
  // grbtag_b->Fit(fit_b);
  // grbtag_q->Fit(fit_q);

  grbtag_q->Draw("AP");
  grbtag_c->Draw("sameP");
  grbtag_b->Draw("sameP");
  grbtag_g->Draw("sameP");

  TLegend *leg4 = new TLegend(0.4,0.4,0.6,0.6);
  leg4->AddEntry(grbtag_b,"b","lp");
  leg4->AddEntry(grbtag_c,"c","lp");
  leg4->AddEntry(grbtag_g,"g","lp");
  leg4->AddEntry(grbtag_q,"uds","lp");
  leg4->SetBorderSize(1);
  leg4->SetFillColor(0);
  leg4->Draw();
  
  if( printplots ) c4->Print("btag2.pdf");
}
