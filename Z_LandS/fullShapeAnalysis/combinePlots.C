#include "Utils/SMS_utils.C"
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
#include <sstream>
#include <iomanip>

using namespace std;

void combinePlots(bool print = false){

  char* version = "V00-00-00";

  TFile *file = TFile::Open(Form("cards/%s/observed_limit.root",version));
  TH2F* hexcl = (TH2F*) file->Get("hexcl");

  // TCanvas *can = new TCanvas("can","",1200,600);
  // can->cd();
  // can->Divide(2,1);

  TCanvas *can = new TCanvas("can","",600,600);
  can->cd();

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // efficiency
  //-------------------------------
  /*  
  can->cd(1);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  TH2F* heff = (TH2F*) file->Get("heff_1");
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  heff->GetXaxis()->SetTitle("gluino mass (GeV)");
  heff->GetZaxis()->SetTitle("efficiency");
  heff->Draw("colz");
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
  t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
  t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 200 GeV");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
  */

  //-------------------------------
  // cross section limit
  //-------------------------------

  //can->cd(2);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  hexcl->GetXaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetLabelSize(0.035);
  hexcl->GetZaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  hexcl->GetXaxis()->SetTitle("gluino mass (GeV)");
  hexcl->GetZaxis()->SetTitle("#sigma upper limit");
  hexcl->Draw("colz");
  hexcl->SetMinimum(0.01);
  hexcl->SetMaximum(10);
  
  TGraph* gr_excl      = getRefXsecGraph(hexcl, "T5zz", 1.0);
  TGraph* gr_excl_down = getRefXsecGraph(hexcl, "T5zz", 1./3.);
  TGraph* gr_excl_up   = getRefXsecGraph(hexcl, "T5zz", 3.);
  
  gr_excl->SetLineWidth(1);
  gr_excl_up->SetLineWidth(1);
  gr_excl_down->SetLineWidth(1);
  gr_excl_up->SetLineStyle(2);
  gr_excl_down->SetLineStyle(3);
  gr_excl->Draw("same");
  gr_excl_up->Draw("same");
  gr_excl_down->Draw("same");
  
  TLegend *leg = new TLegend(0.2,0.53,0.53,0.67);
  leg->AddEntry(gr_excl,     "#sigma^{NLO-QCD}","l");
  leg->AddEntry(gr_excl_up,  "3 #times #sigma^{NLO-QCD}","l");
  leg->AddEntry(gr_excl_down,"1/3 #times #sigma^{NLO-QCD}","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
  t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
  t->DrawLatex(0.2,0.71,"best limit");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");


  if( print ){
    can->Print(Form("cards/%s/plots/SMS.eps",version));
    can->Print(Form("cards/%s/plots/SMS.png",version));
    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS.eps cards/%s/plots/SMS.pdf",version,version));
  }

}
