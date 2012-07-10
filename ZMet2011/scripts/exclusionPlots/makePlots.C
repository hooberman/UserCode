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
#include "TPaveText.h"


using namespace std;

void formatHist( TH2D* hist ){

  hist->GetXaxis()->SetTitle("m(#tilde{#chi}_{2}^{0}) = m(#tilde{#chi}_{1}^{#pm}) [GeV]");
  hist->GetYaxis()->SetTitle("m(#tilde{#chi}_{1}^{0}) [GeV]");
  hist->GetZaxis()->SetTitle("95% CL UL #sigma [fb]");
  hist->GetXaxis()->SetTitleOffset(1.12);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->SetMinimum(1);
  hist->SetMaximum(1000);
  
}

void cmsPrelim( double intLumi )
{

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.18,0.93,"CMS Preliminary,  #sqrt{s}=7 TeV,  L_{int}=4.98 fb^{-1}");
}

void makePlots(){

  //-------------------------------------------
  // Rutgers/KIT
  //-------------------------------------------

  TFile *fkit = TFile::Open("KIT_2i.root");

  TH2D*   h2i     = (TH2D*) fkit->Get("xSecObserved");
  TH2D*   h2i1    = (TH2D*) fkit->Get("xSecObserved1");
  TH2D*   h2i0    = (TH2D*) fkit->Get("xSecObserved0");
  TGraph* gr2i_1  = (TGraph*) fkit->Get("Graph1");
  TGraph* gr2i_2  = (TGraph*) fkit->Get("Graph2");
  TGraph* gr2i_3  = (TGraph*) fkit->Get("Graph3");
  TGraph* gr2i_4  = (TGraph*) fkit->Get("Graph4");
  
  TCanvas *can = new TCanvas();
  can->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(h2i0);
  h2i0->Draw("colz");
  h2i->Draw("colsame");
  h2i1->Draw("colsame");
  gr2i_1->Draw("l");
  gr2i_2->Draw("l");
  gr2i_3->Draw("l");
  gr2i_4->Draw("l");
  cmsPrelim(4.98);

  TLegend *leg = new TLegend(0.2,0.78,0.6,0.88);
  leg->AddEntry(gr2i_1,"observed","l");
  leg->AddEntry(gr2i_2,"median expected","l");
  leg->AddEntry(gr2i_3,"expected #pm1#sigma","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  tex->DrawLatex(0.18,0.73,"pp #rightarrow #tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm}");
  tex->DrawLatex(0.18,0.68,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{e}e, #tilde{#mu}#mu, #tilde{#tau}#tau (BF=50%)");
  tex->DrawLatex(0.18,0.63,"#tilde{#chi}_{1}^{#pm} #rightarrow #tilde{e}#nu_{e}, #tilde{#mu}#nu_{#mu}, #tilde{#tau}#nu_{#tau}");
  tex->DrawLatex(0.18,0.58,"m_{ #tilde{l}} = (m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + m(#tilde{#chi}_{1}^{0})) / 2");









}
