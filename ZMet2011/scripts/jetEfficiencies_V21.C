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



void jetEfficiencies_V21( bool printplots = false ){

  gStyle->SetOptFit(0);

  TFile* f = TFile::Open("../output/V00-02-21/wz_summer11_madgraph.root");

  //TChain *ch = new TChain("T1");
  //ch->Add("../output/V00-02-20/wzsms_gen_baby.root");
  //TFile* f = TFile::Open("../output/V00-02-20/wzsms.root");
  //ch->Add("../output/V00-02-21/wz_summer11_madgraph_gen_baby.root");
  //ch->Add("../output/V00-02-20/ttbar_gen_baby.root");
  //TFile* f = TFile::Open("../output/V00-02-20/ttbar.root");

  //---------------------------------------
  // jets
  //---------------------------------------

  TH1F* hjetpt_all      = (TH1F*) f->Get("hjetpt_all");
  TH1F* hjetpt_q_all    = (TH1F*) f->Get("hjetpt_q_all");
  TH1F* hjetpt_c_all    = (TH1F*) f->Get("hjetpt_c_all");
  TH1F* hjetpt_b_all    = (TH1F*) f->Get("hjetpt_b_all");
  TH1F* hjetpt_g_all    = (TH1F*) f->Get("hjetpt_g_all");

  TH1F* hjetpt_pass     = (TH1F*) f->Get("hjetpt_pass");
  TH1F* hjetpt_q_pass   = (TH1F*) f->Get("hjetpt_q_pass");
  TH1F* hjetpt_c_pass   = (TH1F*) f->Get("hjetpt_c_pass");
  TH1F* hjetpt_b_pass   = (TH1F*) f->Get("hjetpt_b_pass");
  TH1F* hjetpt_g_pass   = (TH1F*) f->Get("hjetpt_g_pass");

  TH1F* hbtag_q_all     = (TH1F*) f->Get("hbtag_q_all");
  TH1F* hbtag_c_all     = (TH1F*) f->Get("hbtag_c_all");
  TH1F* hbtag_b_all     = (TH1F*) f->Get("hbtag_b_all");
  TH1F* hbtag_g_all     = (TH1F*) f->Get("hbtag_g_all");

  TH1F* hbtag_q_passL   = (TH1F*) f->Get("hbtag_q_passL");
  TH1F* hbtag_c_passL   = (TH1F*) f->Get("hbtag_c_passL");
  TH1F* hbtag_b_passL   = (TH1F*) f->Get("hbtag_b_passL");
  TH1F* hbtag_g_passL   = (TH1F*) f->Get("hbtag_g_passL");

  TH1F* hbtag_q_passM   = (TH1F*) f->Get("hbtag_q_passM");
  TH1F* hbtag_c_passM   = (TH1F*) f->Get("hbtag_c_passM");
  TH1F* hbtag_b_passM   = (TH1F*) f->Get("hbtag_b_passM");
  TH1F* hbtag_g_passM   = (TH1F*) f->Get("hbtag_g_passM");

  int rebin = 5;

  hjetpt_all->Rebin(rebin);
  hjetpt_q_all->Rebin(rebin);
  hjetpt_b_all->Rebin(rebin);
  hjetpt_c_all->Rebin(rebin);
  hjetpt_g_all->Rebin(rebin);

  hjetpt_pass->Rebin(rebin);
  hjetpt_q_pass->Rebin(rebin);
  hjetpt_b_pass->Rebin(rebin);
  hjetpt_c_pass->Rebin(rebin);
  hjetpt_g_pass->Rebin(rebin);

  int rebinb = 5;

  hbtag_q_all->Rebin(rebinb);
  hbtag_c_all->Rebin(rebinb);
  hbtag_b_all->Rebin(rebinb);
  hbtag_g_all->Rebin(rebinb);

  hbtag_q_passL->Rebin(rebinb);
  hbtag_c_passL->Rebin(rebinb);
  hbtag_b_passL->Rebin(rebinb);
  hbtag_g_passL->Rebin(rebinb);

  hbtag_q_passM->Rebin(rebinb);
  hbtag_c_passM->Rebin(rebinb);
  hbtag_b_passM->Rebin(rebinb);
  hbtag_g_passM->Rebin(rebinb);

  // cout << "Total entries   " << hjetpt_all->GetEntries()   << endl;
  // cout << "Total q entries " << hjetpt_q_all->GetEntries() << endl;
  // cout << "Total c entries " << hjetpt_b_all->GetEntries() << endl;
  // cout << "Total b entries " << hjetpt_c_all->GetEntries() << endl;
  // cout << "Total g entries " << hjetpt_g_all->GetEntries() << endl;

  // cout << "Total jet entries   " << hjetpt_pass04->GetEntries()   << endl;
  // cout << "Total jet q entries " << hbtag_q_all->GetEntries() << endl;
  // cout << "Total jet c entries " << hbtag_b_all->GetEntries() << endl;
  // cout << "Total jet b entries " << hbtag_c_all->GetEntries() << endl;
  // cout << "Total jet g entries " << hbtag_g_all->GetEntries() << endl;

  //---------------------------------------
  // jet reco
  //---------------------------------------

  TGraphAsymmErrors* grjet = new TGraphAsymmErrors();
  grjet->BayesDivide(hjetpt_pass,hjetpt_all);

  TGraphAsymmErrors* grjet_q = new TGraphAsymmErrors();
  grjet_q->BayesDivide(hjetpt_q_pass,hjetpt_q_all);

  TGraphAsymmErrors* grjet_c = new TGraphAsymmErrors();
  grjet_c->BayesDivide(hjetpt_c_pass,hjetpt_c_all);

  TGraphAsymmErrors* grjet_b = new TGraphAsymmErrors();
  grjet_b->BayesDivide(hjetpt_b_pass,hjetpt_b_all);

  TGraphAsymmErrors* grjet_g = new TGraphAsymmErrors();
  grjet_g->BayesDivide(hjetpt_g_pass,hjetpt_g_all);

  //---------------------------------------
  // b-jet loose
  //---------------------------------------

  TGraphAsymmErrors* grbtag_qL = new TGraphAsymmErrors();
  grbtag_qL->BayesDivide(hbtag_q_passL,hbtag_q_all);

  TGraphAsymmErrors* grbtag_cL = new TGraphAsymmErrors();
  grbtag_cL->BayesDivide(hbtag_c_passL,hbtag_c_all);

  TGraphAsymmErrors* grbtag_bL = new TGraphAsymmErrors();
  grbtag_bL->BayesDivide(hbtag_b_passL,hbtag_b_all);

  TGraphAsymmErrors* grbtag_gL = new TGraphAsymmErrors();
  grbtag_gL->BayesDivide(hbtag_g_passL,hbtag_g_all);

  //---------------------------------------
  // b-jet medium
  //---------------------------------------

  TGraphAsymmErrors* grbtag_qM = new TGraphAsymmErrors();
  grbtag_qM->BayesDivide(hbtag_q_passM,hbtag_q_all);

  TGraphAsymmErrors* grbtag_cM = new TGraphAsymmErrors();
  grbtag_cM->BayesDivide(hbtag_c_passM,hbtag_c_all);

  TGraphAsymmErrors* grbtag_bM = new TGraphAsymmErrors();
  grbtag_bM->BayesDivide(hbtag_b_passM,hbtag_b_all);

  TGraphAsymmErrors* grbtag_gM = new TGraphAsymmErrors();
  grbtag_gM->BayesDivide(hbtag_g_passM,hbtag_g_all);

  //---------------------------------------------
  // jet reco
  //---------------------------------------------
  /*
  TCanvas *c1 = new TCanvas();
  c1->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grjet->SetMaximum(1.05);
  grjet->SetMarkerColor(1);
  grjet->SetMarkerSize(1);
  grjet->SetLineColor(1);
  grjet_b->SetMarkerColor(2);
  grjet_b->SetLineColor(2);
  grjet_c->SetMarkerColor(4);
  grjet_c->SetLineColor(4);
  grjet_q->SetMarkerColor(kGreen+2);
  grjet_q->SetLineColor(kGreen+2);
  grjet_g->SetMarkerColor(6);
  grjet_g->SetLineColor(6);

  TF1*  fit_jet = new TF1("fit_jet", fitf, 20, 200, 7);
  fit_jet->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_jet->FixParameter(6,0);
  fit_jet->SetLineWidth(3);
  fit_jet->SetLineColor(1);

  grjet->Fit(fit_jet);

  grjet->GetXaxis()->SetNdivisions(7);
  grjet->GetXaxis()->SetTitle("parton p_{T} [GeV]");
  grjet->GetYaxis()->SetTitle("jet reco efficiency");
  grjet->Draw("AP");
  // grjet_b->Draw("sameP");
  // grjet_c->Draw("sameP");
  // grjet_q->Draw("sameP");
  // grjet_g->Draw("sameP");

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
  line.DrawLine(30,0,30,1.05);

  if( printplots ) c1->Print("../plots/jetreco.pdf");
  */

  //---------------------------------
  // b-tagging (loose)
  //---------------------------------

  TCanvas *c2 = new TCanvas();
  c2->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grbtag_qL->SetMarkerColor(kGreen+2);
  grbtag_qL->SetLineColor(kGreen+2);
  grbtag_cL->SetMarkerColor(4);
  grbtag_cL->SetLineColor(4);
  grbtag_bL->SetMarkerColor(2);
  grbtag_bL->SetLineColor(2);
  grbtag_gL->SetMarkerColor(6);
  grbtag_gL->SetLineColor(6);

  grbtag_qL->GetXaxis()->SetRangeUser(0,200);
  grbtag_qL->GetXaxis()->SetTitle("parton p_{T} [GeV]");
  grbtag_qL->GetYaxis()->SetTitle("loose b-tagging efficiency");
  grbtag_qL->SetMinimum(0.001);
  grbtag_qL->SetMaximum(1.0);

  TF1*  fit_cL = new TF1("fit_cL", fitf, 30, 200, 7);
  fit_cL->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_cL->FixParameter(6,0);
  fit_cL->SetLineColor(4);
  fit_cL->SetLineWidth(2);

  TF1*  fit_bL = new TF1("fit_bL", fitf, 30, 200, 7);
  fit_bL->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_bL->FixParameter(6,0);
  fit_bL->SetLineColor(2);
  fit_bL->SetLineWidth(2);

  TF1*  fit_qL = new TF1("fit_qL", fitf, 30, 200, 7);
  fit_qL->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_qL->FixParameter(6,0);
  fit_qL->SetLineColor(kGreen+2);
  fit_qL->SetLineWidth(2);

  TF1*  fit_gL = new TF1("fit_gL", fitf, 30, 200, 7);
  fit_gL->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_gL->FixParameter(6,0);
  fit_gL->SetLineColor(6);
  fit_gL->SetLineWidth(2);

  grbtag_qL->Fit(fit_qL,"R");
  grbtag_bL->Fit(fit_bL,"R");
  grbtag_cL->Fit(fit_cL,"R");
  grbtag_gL->Fit(fit_gL,"R");

  grbtag_qL->Draw("AP");
  grbtag_cL->Draw("sameP");
  grbtag_bL->Draw("sameP");
  grbtag_gL->Draw("sameP");

  TLegend *leg2 = new TLegend(0.5,0.3,0.7,0.6);
  leg2->AddEntry(grbtag_qL,"uds","lp");
  leg2->AddEntry(grbtag_cL,"c","lp");
  leg2->AddEntry(grbtag_bL,"b","lp");
  leg2->AddEntry(grbtag_gL,"g","lp");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(1);
  //leg2->Draw();
  
  if( printplots ) c2->Print("../plots/btagL.pdf");


  //---------------------------------
  // b-tagging (medium)
  //---------------------------------
  /*
  TCanvas *c3 = new TCanvas();
  c3->cd();

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.05);  

  grbtag_qM->SetMarkerColor(kGreen+2);
  grbtag_qM->SetLineColor(kGreen+2);
  grbtag_cM->SetMarkerColor(4);
  grbtag_cM->SetLineColor(4);
  grbtag_bM->SetMarkerColor(2);
  grbtag_bM->SetLineColor(2);
  grbtag_gM->SetMarkerColor(6);
  grbtag_gM->SetLineColor(6);

  grbtag_qM->GetXaxis()->SetRangeUser(0,200);
  grbtag_qM->GetXaxis()->SetTitle("parton p_{T} [GeV]");
  grbtag_qM->GetYaxis()->SetTitle("medium b-tagging efficiency");
  grbtag_qM->SetMinimum(0.001);
  grbtag_qM->SetMaximum(1.0);

  TF1*  fit_cM = new TF1("fit_cM", fitf, 30, 200, 7);
  fit_cM->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_cM->FixParameter(6,0);
  fit_cM->SetLineColor(4);
  fit_cM->SetLineWidth(2);

  TF1*  fit_bM = new TF1("fit_bM", fitf, 30, 200, 7);
  fit_bM->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_bM->FixParameter(6,0);
  fit_bM->SetLineColor(2);
  fit_bM->SetLineWidth(2);

  TF1*  fit_qM = new TF1("fit_qM", fitf, 30, 200, 7);
  fit_qM->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_qM->FixParameter(6,0);
  fit_qM->SetLineColor(kGreen+2);
  fit_qM->SetLineWidth(2);

  TF1*  fit_gM = new TF1("fit_gM", fitf, 30, 200, 7);
  fit_gM->SetParameters( 60 , 60 , 20 , 20 , 0.1 , 0.1 , 0 );
  fit_gM->FixParameter(6,0);
  fit_gM->SetLineColor(6);
  fit_gM->SetLineWidth(2);

  grbtag_qM->Fit(fit_qM,"R");
  grbtag_bM->Fit(fit_bM,"R");
  grbtag_cM->Fit(fit_cM,"R");
  grbtag_gM->Fit(fit_gM,"R");

  grbtag_qM->Draw("AP");
  grbtag_cM->Draw("sameP");
  grbtag_bM->Draw("sameP");
  grbtag_gM->Draw("sameP");
  
  if( printplots ) c3->Print("../plots/btagM.pdf");
*/




  // grbtag_c->Fit(fit_c);
  // grbtag_b->Fit(fit_b);
  // grbtag_q->Fit(fit_q);

}
