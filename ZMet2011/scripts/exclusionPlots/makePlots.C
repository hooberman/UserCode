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

  TH2D* h2i  = (TH2D*) fkit->Get("xSecObserved");
  TH2D* h2i1 = (TH2D*) fkit->Get("xSecObserved1");
  TH2D* h2i0 = (TH2D*) fkit->Get("xSecObserved0");

  

  TCanvas *can = new TCanvas();
  can->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  gPad->SetLogz();
  formatHist(h2i0);
  h2i0->Draw("colz");
  h2i->Draw("colsame");
  h2i1->Draw("colsame");
  cmsPrelim(4.98);

  //h2i->Draw("colsame");
















}
