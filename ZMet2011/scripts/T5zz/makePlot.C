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
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector.h"
#include <sstream>
#include <iomanip>



using namespace std;

TH1F* smoothHist( TH1F* hin , int n ){

  if( n%2==0 ){
    cout << "ERROR! n " << n << " must be odd!" << endl;
    exit(0);
  }

  //int firstbin = (n+1)/2;
  //int lastbin  = hin->GetXaxis()->GetNbins() - firstbin + 1;
  int diff     = (n-1)/2;

  TH1F* hout = (TH1F*) hin->Clone(Form("%s_clone",hin->GetName()));

  //for( int ibin = firstbin ; ibin <= lastbin ; ++ibin ){
  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ++ibin ){
    float ave = 0;

    float skip = false;

    for( int jbin = ibin-diff ; jbin <= ibin+diff ; ++jbin ){
      if( hin->GetBinContent(jbin) < 1e-10   ) skip = true;
      if( jbin > hin->GetXaxis()->GetNbins() ) skip = true;
      ave += hin->GetBinContent(jbin);
    }

    if( skip ){
      hout->SetBinContent(ibin,0);
      continue;
    }

    ave = ave / (float) n;

    hout->SetBinContent(ibin,ave);

  }
  
  return hout;
}

TGraph* getGraph(TH1F* hist){


  vector<float> x;
  vector<float> y;

  for(int ibin = 1 ; ibin <= hist->GetXaxis()->GetNbins() ; ibin++ ){
    if( hist->GetBinContent(ibin) > 1.e-10 ){
      x.push_back( hist->GetBinCenter(ibin) );
      y.push_back( hist->GetBinContent(ibin) );
    }
  }

  const unsigned int n = x.size();
  float x_[n];
  float y_[n];

  for( unsigned int i = 0 ; i < n ; ++i ){
    x_[i] = x[i];
    y_[i] = y[i];
  }

  TGraph *gr = new TGraph(n,x_,y_);

  return gr;
  
}

void makePlot(){

  gStyle->SetPaintTextFormat(".1f");

  TFile *f = TFile::Open("histos.root");
  TH2F* hist = (TH2F*) f->Get("h_T5zzgmsb");

  TFile *fxsec = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* hxsec = (TH1F*) fxsec->Get("gluino");


  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.25);

  hist->Draw("colz");
  hist->Draw("sametext");


  TH1F* h100 = (TH1F*) hist->ProjectionX("h100",5,5);
  TH1F* h200 = (TH1F*) hist->ProjectionX("h200",9,9);
  TH1F* h300 = (TH1F*) hist->ProjectionX("h300",13,13);
  TH1F* h400 = (TH1F*) hist->ProjectionX("h400",17,17);
  TH1F* h500 = (TH1F*) hist->ProjectionX("h500",21,21);

  // cout << "m(chi10) = 100 GeV" << endl;
  // for(unsigned int ibin = 1 ; ibin <= h100->GetXaxis()->GetNbins() ; ibin++ ){
  //   cout << ibin << " " << h100->GetBinCenter(ibin) << " " << h100->GetBinContent(ibin) << endl;
  // }
  // cout << "m(chi10) = 200 GeV" << endl;
  // for(unsigned int ibin = 1 ; ibin <= h100->GetXaxis()->GetNbins() ; ibin++ ){
  //   cout << ibin << " " << h100->GetBinCenter(ibin) << " " << h200->GetBinContent(ibin) << endl;
  // }
  // cout << "m(chi10) = 300 GeV" << endl;
  // for(unsigned int ibin = 1 ; ibin <= h100->GetXaxis()->GetNbins() ; ibin++ ){
  //   cout << ibin << " " << h100->GetBinCenter(ibin) << " " << h300->GetBinContent(ibin) << endl;
  // }



  TCanvas *c2 = new TCanvas();
  c2->cd();
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.1);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();

  //--------------------------------------
  // smooth histograms
  //--------------------------------------

  TH1F* h100_smooth = smoothHist(h100,7);
  TH1F* h200_smooth = smoothHist(h200,7);
  TH1F* h300_smooth = smoothHist(h300,7);
  TH1F* h400_smooth = smoothHist(h400,7);
  TH1F* h500_smooth = smoothHist(h500,7);

  //--------------------------------------
  // convert: histograms to graphs
  //--------------------------------------

  TGraph* gxsec = getGraph(hxsec);
  // TGraph* g100 = getGraph(h100);
  // TGraph* g200 = getGraph(h200);
  // TGraph* g300 = getGraph(h300);
  // TGraph* g400 = getGraph(h400);
  // TGraph* g500 = getGraph(h500);
  TGraph* g100 = getGraph(h100_smooth);
  TGraph* g200 = getGraph(h200_smooth);
  TGraph* g300 = getGraph(h300_smooth);
  TGraph* g400 = getGraph(h400_smooth);
  TGraph* g500 = getGraph(h500_smooth);

  hxsec->SetLineWidth(3);
  gxsec->SetLineWidth(3);

  g100->SetLineWidth(3);
  g200->SetLineWidth(3);
  g300->SetLineWidth(3);
  g400->SetLineWidth(3);
  g500->SetLineWidth(3);

  g100->SetLineColor(kGreen+2);
  g200->SetLineColor(2);
  g300->SetLineColor(4);
  g400->SetLineColor(4);
  g500->SetLineColor(2);

  g100->SetLineStyle(2);
  g200->SetLineStyle(4);
  g300->SetLineStyle(3);
  g400->SetLineStyle(3);
  g500->SetLineStyle(4);

  h100->SetLineWidth(2);
  h200->SetLineWidth(2);
  h300->SetLineWidth(2);
  h400->SetLineWidth(2);
  h500->SetLineWidth(2);

  h100->SetLineColor(kGreen+2);
  h200->SetLineColor(2);
  h300->SetLineColor(4);
  h400->SetLineColor(4);
  h500->SetLineColor(2);

  gxsec->SetMinimum(0.001);
  gxsec->SetMaximum(10000);
  gxsec->GetYaxis()->SetLabelSize(0.04);
  gxsec->GetXaxis()->SetLabelSize(0.04);
  gxsec->GetYaxis()->SetTitle("#sigma #times BF [pb]");
  gxsec->GetXaxis()->SetTitle("gluino mass [GeV]");
  gxsec->GetXaxis()->SetRangeUser(200,1100);
  //gxsec->GetXaxis()->SetLimits(200,1100);

  gxsec->Draw("Al");
  //gxsec->SetLineColor(2);
  gxsec->Draw("samel");
  //h100->Draw("same");
  //h200->Draw("same");
  //h300->Draw("same");
  g100->Draw("samec");
  //g200->Draw("samel");
  g300->Draw("samec");
  //g400->Draw("samel");
  g500->Draw("samec");
  
  // TF1* f100 = new TF1("f100","[0]*TMath::Exp(x/[1])",200,1100);
  // f100->SetParameter(0,10);
  // f100->SetParameter(1,200);

  // h100->SetLineWidth(4);
  // h300->SetLineWidth(4);
  // h500->SetLineWidth(4);
  //h100->Fit(f100,"R");

  //h100->Draw("same");
  //h300->Draw("same");
  //h500->Draw("same");


  TLegend *leg = new TLegend(0.45,0.68,0.9,0.88);
  leg->AddEntry(hxsec  ,"#sigma(pp#rightarrow#tilde{g}#tilde{g})   NLO-QCD","l");
  leg->AddEntry(g100,"#sigma_{UL} M(#chi^{0}_{1}) = 100 GeV","l");
  //leg->AddEntry(g200,"#sigma_{UL} M(#chi^{0}_{1}) = 200 GeV","l");
  leg->AddEntry(g300,"#sigma_{UL} M(#chi^{0}_{1}) = 300 GeV","l");
  //leg->AddEntry(g400,"#sigma_{UL} M(#chi^{0}_{1}) = 400 GeV","l");
  leg->AddEntry(g500,"#sigma_{UL} M(#chi^{0}_{1}) = 500 GeV","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.04);
  t->SetNDC();
  t->DrawLatex(0.2,0.92,"CMS preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}L dt = 4.7 fb^{-1}");

  t->DrawLatex(0.45,0.62,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{1}^{0}, #chi_{1}^{0} #rightarrow Z G");
  t->DrawLatex(0.45,0.57,"m(#tilde{q})>>m(#tilde{g})");
  t->DrawLatex(0.45,0.52,"n_{jets} #geq 2");

  c2->Print("makePlot.pdf");
  c2->Print("makePlot.eps");
  c2->Print("makePlot.png");
  gROOT->ProcessLine(".! ps2pdf makePlot.eps makePlot_ppt.pdf");

}
