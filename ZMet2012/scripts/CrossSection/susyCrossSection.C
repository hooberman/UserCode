#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

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
#include "C1C1.C"
#include "C1N2.C"
#include "slepton.C"
#include "stop.C"
#include "gluino.C"



using namespace std;

void susyCrossSection(){

  //---------------------------------
  // get total cross section
  //---------------------------------

  TGraph *C1C1    = gr_C1C1();
  TGraph *C1N2    = gr_C1N2();
  TGraph *slepton = gr_slepton();
  TGraph *stop    = gr_stop();
  TGraph *gluino  = gr_gluino();


  //---------------------------------
  // create and scale histograms
  //---------------------------------

  C1N2->SetLineWidth(2);
  C1N2->SetLineColor(7);

  C1C1->SetLineWidth(2);
  C1C1->SetLineColor(7);

  slepton->SetLineWidth(2);
  slepton->SetLineColor(7);

  stop->SetLineWidth(2);
  stop->SetLineColor(4);

  stop->SetLineWidth(2);
  stop->SetLineColor(2);


  //---------------------------------
  // make plots
  //---------------------------------

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetLogy();

  C1N2->Draw("Al");
  C1C1->Draw("samel");
  slepton->Draw("samel");
  stop->Draw("samel");

  // h_C1N2->GetXaxis()->SetTitle("SUSY particle mass [GeV]");
  // h_C1N2->GetYaxis()->SetTitle("#sigma(pp#rightarrow SUSY) [pb]");
  // h_C1N2->Draw("c");
  // h_C1N2->SetMinimum(0.0001);
  // h_C1N2->SetMaximum(100);
  // h_C1N2->GetXaxis()->SetRangeUser(125,400);

  TLatex *t = new TLatex();
  t->SetNDC();

  t->SetTextSize(0.05);
  t->DrawLatex(0.2,0.92,"8 TeV NLO cross sections");

  // TLegend *leg = new TLegend(0.6,0.57,0.85,0.858);
  // leg->SetFillColor(0);
  // leg->SetBorderSize(0);
  // leg->AddEntry(h,"total","l");
  // leg->AddEntry(h_Wlv_Hbb,"W(#font[12]{l}#nu)H(b#bar{b})","l");
  // leg->AddEntry(h_Wlv_HWW,"W(#font[12]{l}#nu)H(WW)","l");
  // leg->AddEntry(h_Wlv_Htt,"W(#font[12]{l}#nu)H(#tau#tau)","l");
  // leg->AddEntry(h_Wlv_HZZ,"W(#font[12]{l}#nu)H(ZZ)","l");
  // leg->AddEntry(h_Wjj_Hgg,"W(jj)H(#gamma#gamma)","l");
  // leg->AddEntry(h_Wlv_Hgg,"W(#font[12]{l}#nu)H(#gamma#gamma)","l");

  // leg->Draw();

  c1->Print("../plots/susyCrossSection.pdf");


}
