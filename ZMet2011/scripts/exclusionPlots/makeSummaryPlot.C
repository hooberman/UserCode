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

bool drawExpected  = false;
bool isPreliminary = true;
char* isPrelimChar = (char*) "";


void removeDiagonal( TH2D* h , float deltaM ){

  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

      float mg = h->GetXaxis()->GetBinCenter(ibin);
      float ml = h->GetYaxis()->GetBinCenter(jbin);

      cout << "mg ml " << mg << " " << ml << endl;

      if( mg - ml < deltaM ) h->SetBinContent(ibin,jbin,0);
      
    }
  }
}

void formatHist( TH2D* hist ){

  hist->GetXaxis()->SetTitle("m(#tilde{#chi}_{2}^{0}) = m(#tilde{#chi}_{1}^{#pm}) [GeV]");
  hist->GetYaxis()->SetTitle("m(#tilde{#chi}_{1}^{0}) [GeV]");
  hist->GetZaxis()->SetTitle("95% CL UL #sigma#timesBF [fb]");
  hist->GetXaxis()->SetTitleOffset(1.12);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->SetMinimum(1);
  hist->SetMaximum(1000);
  
}

// void cmsPrelim( double intLumi )
// {

//   TLatex latex;
//   latex.SetNDC();
//   latex.SetTextSize(0.04);
//   latex.DrawLatex(0.18,0.93,"CMS Preliminary,  #sqrt{s}=7 TeV,  L_{int}=4.98 fb^{-1}");
// }

void cmsPrelim(double intLumi, bool prelim)
{
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62);
        if(prelim) latex.SetTextSize(0.04);
        else       latex.SetTextSize(0.045);

        latex.SetTextAlign(11); // align left
        if(prelim) latex.DrawLatex(0.17,0.92,"CMS Preliminary");
        else       latex.DrawLatex(0.17,0.92,"CMS");

        latex.SetTextAlign(31); // align right
        latex.DrawLatex(0.89, 0.92, Form("#sqrt{s} = 7 TeV, L_{int} = %4.2f fb^{-1}", intLumi));
}

TGraph *getSleptonGraph(){

  float x[4];
  float y[4];

  //http://pdglive.lbl.gov/Rsummary.brl?nodein=S046&exp=Y
  float exc = 82.0;

  x[0] = 100.0; y[0] = 0;
  x[1] = 2*exc; y[1] = 0;
  x[2] = 100;   y[2] = 2*exc-100;
  x[3] = 100.0; y[3] = 0;

  TGraph *gr = new TGraph(4,x,y);
  gr->SetFillColor(kGray);
  gr->SetLineColor(0);

  return gr;
}

TGraph *getCharginoGraph(){

  float x[4];
  float y[4];

  //http://lepsusy.web.cern.ch/lepsusy/www/inos_moriond01/charginos_pub.html
  float exc = 103.5;

  x[0] = 100.0; y[0] =     0;
  x[1] =   exc; y[1] =     0;
  x[2] =   exc; y[2] =   exc;
  x[3] = 100.0; y[3] = 100.0;
  x[4] = 100.0; y[4] =     0;

  TGraph *gr = new TGraph(5,x,y);
  gr->SetFillColor(kRed);
  gr->SetLineColor(0);

  return gr;
}


void makeSummaryPlot(){

  if( isPreliminary ) isPrelimChar = (char*) "_prelim";

  //-------------------------------------------
  // Rutgers/KIT
  //-------------------------------------------

  //-----------------
  // model 2i
  //-----------------

  TFile *f2i      = TFile::Open("KIT_2i.root");
  TFile *f2a      = TFile::Open("KIT_2a.root");
  TFile *fwz      = TFile::Open("combinePlots_VZ_Trilepton.root");
  TFile *f2imt     = TFile::Open("LeftSlepton_50.root");
  TFile *f2amt     = TFile::Open("TauEnriched_50.root");

  TGraph* gr2i          = (TGraph*) f2i->Get("Graph1");
  TGraph* gr2a          = (TGraph*) f2a->Get("Graph1");
  TGraph* grwz          = (TGraph*) fwz->Get("gr_combo");
  TGraph* gr2imt        = (TGraph*) f2imt->Get("ObservedExclusion");
  TGraph* gr2amt        = (TGraph*) f2amt->Get("ObservedExclusion");

  TGraph* gr2i_exp      = (TGraph*) f2i->Get("Graph2");
  TGraph* gr2a_exp      = (TGraph*) f2a->Get("Graph2");
  TGraph* grwz_exp      = (TGraph*) fwz->Get("gr_combo_exp");
  TGraph* gr2imt_exp    = (TGraph*) f2imt->Get("ExpectedExclusion");
  TGraph* gr2amt_exp    = (TGraph*) f2amt->Get("ExpectedExclusion");

  TGraph* grslepton     = getSleptonGraph();
  TGraph* grchargino    = getCharginoGraph();

  //------------------
  // observed
  //------------------

  gr2i->SetLineColor(1);
  gr2i->SetFillColor(kGray);
  gr2i->SetLineWidth(3);
  if( !drawExpected ) gr2i->SetLineStyle(1);

  gr2a->SetLineColor(2);
  gr2a->SetFillColor(2);
  gr2a->SetLineWidth(3);
  if( !drawExpected ) gr2a->SetLineStyle(2);

  grwz->SetLineColor(4);
  grwz->SetFillColor(4);
  grwz->SetLineWidth(3);
  if( !drawExpected ) grwz->SetLineStyle(3);

  gr2imt->SetLineColor(kGreen+2);
  gr2imt->SetFillColor(kGreen+2);
  gr2imt->SetLineWidth(3);
  if( !drawExpected ) gr2imt->SetLineStyle(4);

  gr2amt->SetLineColor(6);
  gr2amt->SetFillColor(6);
  gr2amt->SetLineWidth(3);
  if( !drawExpected ) gr2amt->SetLineStyle(5);

  //------------------
  // expected
  //------------------

  gr2i_exp->SetLineColor(1);
  gr2i_exp->SetLineWidth(3);
  gr2i_exp->SetLineStyle(2);

  gr2a_exp->SetLineColor(2);
  gr2a_exp->SetLineWidth(3);
  gr2a_exp->SetLineStyle(2);

  grwz_exp->SetLineColor(4);
  grwz_exp->SetLineWidth(3);
  grwz_exp->SetLineStyle(2);

  gr2imt_exp->SetLineColor(kGreen+2);
  gr2imt_exp->SetLineWidth(3);
  gr2imt_exp->SetLineStyle(2);

  gr2amt_exp->SetLineColor(6);
  gr2amt_exp->SetLineWidth(3);
  gr2amt_exp->SetLineStyle(2);



  grslepton->SetFillColor(kGreen+2);
  grslepton->SetFillStyle(3002);

  TCanvas *can = new TCanvas();
  can->cd();
  gPad->SetTopMargin(0.1);

  TH2D* hdummy = new TH2D("hdummy","",500,100,600,800,0,800);

  formatHist(hdummy);
  hdummy->Draw();

  grslepton->Draw("f");
  grchargino->Draw("samef");

  gr2imt->Draw("l");
  gr2i->Draw("l");
  gr2a->Draw("l");
  gr2amt->Draw("l");
  grwz->Draw("l");

  if( drawExpected ){
    gr2i_exp->Draw("l");
    gr2a_exp->Draw("l");
    grwz_exp->Draw("l");
    gr2imt_exp->Draw("l");
    gr2amt_exp->Draw("l");
  }

  hdummy->Draw("sameaxis");

  TLegend *leg = new TLegend(0.18,0.6,0.78,0.88);
  leg->AddEntry(grslepton  ,"LEP2 slepton limit m( #tilde{#font[12]{l}} ) > 82 GeV","f");
  leg->AddEntry(grchargino ,"LEP2 chargino limit m(#chi^{#pm}) > 103.5 GeV","f");
  leg->AddEntry(gr2i  ,"3#font[12]{l}+E_{T}^{miss} ( #tilde{#font[12]{l}}_{L} , BF(3#font[12]{l})=0.5)","l");
  leg->AddEntry(gr2a  ,"3#font[12]{l}+E_{T}^{miss} ( #tilde{#font[12]{l}}_{R} , BF(#font[12]{l^{+}l^{-}}#tau)=1)","l");
  leg->AddEntry(gr2imt,"3#font[12]{l}+M_{#font[12]{ll}}+M_{T} & 2#font[12]{l}(SS) ( #tilde{#font[12]{l}}_{L} , BF(3#font[12]{l})=0.5)","l");
  leg->AddEntry(gr2amt,"3#font[12]{l}+M_{#font[12]{ll}}+M_{T} & 2#font[12]{l}(SS) ( #tilde{#font[12]{l}}_{R} , BF(#font[12]{l^{+}l^{-}}#tau)=1)","l");
  leg->AddEntry(grwz  ,"2#font[12]{l}2j & 3#font[12]{l}+M_{#font[12]{ll}}+M_{T} (no #tilde{#font[12]{l}} , BF(WZ)=1)","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TH1F *hobs = new TH1F();
  TH1F *hexp = new TH1F();

  hobs->SetLineWidth(3);
  hexp->SetLineWidth(3);
  hexp->SetLineStyle(2);

  TLegend *leg2 = new TLegend(0.65,0.45,0.85,0.55);
  leg2->AddEntry(hobs  ,"observed","l");
  leg2->AddEntry(hexp  ,"expected","l");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  if( drawExpected ) leg2->Draw();



  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.03);  
  tex->DrawLatex(0.2,0.55,"m( #tilde{#font[12]{l}} ) = 0.5m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) + 0.5m(#tilde{#chi}_{1}^{0})");

  tex->SetTextAngle(32);
  tex->DrawLatex(0.25,0.31,"m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{#pm}) > m(#tilde{#chi}_{1}^{0})");

  TLine line;

  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.DrawLine(100.0,100.0,600.0,600.0);

  cmsPrelim(4.98,isPreliminary);

  if( isPreliminary ){
    if( drawExpected ) can->Print("ewkino_summaryPlot_expected_prelim.pdf");
    else               can->Print("ewkino_summaryPlot_prelim.pdf");
  }
  else{
    if( drawExpected ) can->Print("ewkino_summaryPlot_expected.pdf");
    else               can->Print("ewkino_summaryPlot.pdf");
  }

}


TH2D* cloneHist( TH2D* hin ){

  TH2D* hout = new TH2D(Form("%s_clone",hin->GetName()),Form("%s_clone",hin->GetName()),
			hin->GetXaxis()->GetNbins(),hin->GetXaxis()->GetXmin(),hin->GetXaxis()->GetXmax(),
			hin->GetYaxis()->GetNbins(),hin->GetYaxis()->GetXmin(),hin->GetYaxis()->GetXmax());

  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ++ibin ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; ++jbin ){
      hout->SetBinContent(ibin,jbin,hin->GetBinContent(ibin,jbin));
    }
  }

  return hout;

}
