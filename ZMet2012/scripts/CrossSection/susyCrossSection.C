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
#include "squark.C"



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
  TGraph *squark  = gr_squark();


  //---------------------------------
  // create and scale histograms
  //---------------------------------

  C1N2->SetLineWidth(4);
  C1N2->SetLineColor(6);

  C1C1->SetLineWidth(4);
  C1C1->SetLineColor(6);
  C1C1->SetLineStyle(9);

  slepton->SetLineWidth(4);
  slepton->SetLineColor(6);
  slepton->SetLineStyle(2);

  stop->SetLineWidth(4);
  stop->SetLineColor(4);

  gluino->SetLineWidth(4);
  gluino->SetLineColor(2);

  squark->SetLineWidth(4);
  squark->SetLineColor(2);
  squark->SetLineStyle(9);

  //---------------------------------
  // make plots
  //---------------------------------

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetLogy();

  TH2F* hdummy = new TH2F("hdummy","hdummy",100,100,1500,100,0.0005,100);

  hdummy->Draw();

  C1N2->Draw("samel");
  C1C1->Draw("samel");
  slepton->Draw("samel");
  stop->Draw("samel");
  gluino->Draw("samel");
  squark->Draw("samel");

  hdummy->GetXaxis()->SetTitle("SUSY particle mass [GeV]");
  hdummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow SUSY) [pb]");

  TLatex *t = new TLatex();
  t->SetNDC();

  t->SetTextSize(0.05);
  t->DrawLatex(0.2,0.92,"8 TeV NLO cross sections");

  TLegend *leg = new TLegend(0.65,0.4,0.85,0.88);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(gluino  ,"#tilde{g}#tilde{g}","l");
  leg->AddEntry(squark  ,"#tilde{q}#bar{#tilde{q}}","l");
  leg->AddEntry(stop    ,"#tilde{t_{1}}#bar{#tilde{t_{1}}}","l");
  leg->AddEntry(C1N2    ,"#tilde{#chi}_{1}^{#pm}#tilde{#chi}_{2}^{0} (wino)","l");
  leg->AddEntry(C1C1    ,"#tilde{#chi}_{1}^{+}#tilde{#chi}_{1}^{-} (wino)","l");
  leg->AddEntry(slepton ,"#tilde{#font[12]{l}_{e}^{+}}#tilde{#font[12]{l}_{e}^{-}}","l");
  leg->SetTextSize(0.04);
  leg->Draw();

  c1->Print("susyCrossSection.pdf");


}
