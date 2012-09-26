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
#include "TF1.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>

using namespace std;

TGraph* getGraph_2011(){

  float x[6];
  float y[6];
  int npoints = -1;

  x[0] =  212.5;  y[0] = -12.5;
  x[1] =  212.5;  y[1] =  25.0;
  x[2] =  200.0;  y[2] =  37.5;
  x[3] =  150.0;  y[3] =  37.5;
  x[4] =  137.5;  y[4] =  25.0;
  x[5] =  137.5;  y[5] = -12.5;
  npoints = 6;

  TGraph *gr = new TGraph(npoints,x,y);

  gr->SetLineWidth(3);
  gr->SetLineStyle(2);
  return gr;

}

TGraph* getGraph_expected(){

  float x[7];
  float y[7];
  int npoints = -1;

  x[0] =  255.0;  y[0] =  -5.0;
  x[1] =  245.0;  y[1] =  30.0;
  x[2] =  215.0;  y[2] =  35.0;
  x[3] =  195.0;  y[3] =  35.0;
  x[4] =  180.0;  y[4] =  25.0;
  x[5] =  165.0;  y[5] =   5.0;
  x[6] =  165.0;  y[6] =  -5.0;
  npoints = 7;

  TGraph *gr = new TGraph(npoints,x,y);
  gr->SetLineWidth(3);

  return gr;

}

TGraph* getGraph_observed(){

  float x[7];
  float y[7];
  int npoints = -1;

  x[0] =  285.0;  y[0] =  -5.0;
  x[1] =  275.0;  y[1] =  30.0;
  x[2] =  235.0;  y[2] =  35.0;
  x[3] =  215.0;  y[3] =  35.0;
  x[4] =  200.0;  y[4] =  25.0;
  x[5] =  185.0;  y[5] =   5.0;
  x[6] =  185.0;  y[6] =  -5.0;
  npoints = 7;

  TGraph *gr = new TGraph(npoints,x,y);
  gr->SetLineWidth(3);

  return gr;

}

void smoothHist( TH2F* h ){

  vector<int> binx;
  vector<int> biny;
  binx.push_back(19);    biny.push_back(8);
  binx.push_back(33);    biny.push_back(7);
  binx.push_back(35);    biny.push_back(7);

  binx.push_back(28);    biny.push_back(26);
  binx.push_back(34);    biny.push_back(30);
  binx.push_back(34);    biny.push_back(31);

  binx.push_back(42);    biny.push_back(37);
  binx.push_back(45);    biny.push_back(37);
  binx.push_back(45);    biny.push_back(36);

  const unsigned int npoints = binx.size();

  for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

      float val = h->GetBinContent(ibin,jbin);

      if( val < 1e-10 ) continue;
      //cout << ibin << " " << jbin << " " << val << endl;

      for(unsigned int ipoint = 0 ; ipoint < npoints ; ipoint++ ){
	if( ibin == binx.at(ipoint) && jbin == biny.at(ipoint) ){

	  float valup  = h->GetBinContent(ibin+1,jbin);
	  float valdn  = h->GetBinContent(ibin-1,jbin);
	  float valavg = 0.5 * (valup+valdn);

	  h->SetBinContent(ibin,jbin,valavg);
	}
	
      }
    }
  }  
}

void plotProjections( TH2F* h , TH1F* hxsec ){

  const unsigned int n = 1;

  TH1F*    hproj[n];
  TCanvas* can[n];

  for( int i = 0 ; i < n ; i++ ){

    cout << endl;
    cout << "LSP mass " << 10*i << endl;
    hproj[i] = (TH1F*) h->ProjectionX(Form("hproj_%i",i),i+1,i+1);
    hproj[i]->SetLineColor(2);

    for( int ibin = 1 ; ibin <= hproj[i]->GetXaxis()->GetNbins() ; ibin++ ){
      cout << ibin << " " << hproj[i]->GetBinCenter(ibin) << " " << hproj[i]->GetBinContent(ibin) << endl;
    }

    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),600,600);
    can[i]->cd();
    gPad->SetLogy();

    TF1* f = new TF1("f","[0]*TMath::Exp(-1*x/[1])",150,300);

    hproj[i]->Fit(f,"R");
    hproj[i]->SetMinimum(0.001);
    hproj[i]->Draw();
    hxsec->Draw("same");
  }



}

void combinePlots(string version = "V00-00-02" , bool print = false){

  char* xtitle  = "m_{#chi_{2}^{0}} = m_{#chi_{1}^{#pm}} [GeV]";
  char* ytitle  = "m_{#chi_{1}^{0}} [GeV]";

  bool smooth = false;

  char* sample;
  char* title;
  char* xsectype;
  char* xsechist;
  float denom;

  if( version == "V00-00-02" ){
    sample   = (char*) "wzsms";
    title    = (char*) "pp#rightarrow #chi^{#pm}#chi^{0} #rightarrow WZ + E_{T}^{miss}";
    xsectype = (char*) "C1N2_8TeV_finer";
    xsechist = (char*) "C1N2_8TeV_NLO";
    denom    = 100000.0;
  } 
  
  //float ymin = 0.;

  TFile *file = TFile::Open(Form("cards/%s/observed_limit.root",version.c_str()));

  TH2F* hexcl      = (TH2F*) file->Get("hexcl");
  TH2F* hexp       = (TH2F*) file->Get("hexp");
  
  // hexcl->Scale(0.9);
  //hexp->Scale(sqrt(9.2/15.0));

  // hexcl->RebinX(2);
  // hexcl->RebinY(2);

  // hexp->RebinX(2);
  // hexp->RebinY(2);

  // hexcl->Scale(0.25);
  // hexp->Scale(0.25);

  //hexcl->Scale(1000.0); // pb --> fb

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // calculate efficiency
  //-------------------------------

  char* babyversion = (char*) "V00-01-05";

  TChain *ch = new TChain("T1");
  ch->Add(Form("output/%s/%s_baby_oldIso.root",babyversion,sample));

  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets >= 2");
  TCut bveto("nbcsvm==0");
  TCut mjj("mjj>70 && mjj<110");
  TCut nlep2("nlep==2");
  TCut rho("rho>0 && rho<40");
  TCut sf("leptype==0||leptype==1");
  TCut met80("pfmet>80");
  TCut weight("vtxweight * trgeff");

  TCut sel  = pt2020 + zmass + njets2 + bveto + mjj + nlep2 + sf + met80;

  cout << "Using selection: " << sel.GetTitle() << endl;

  int   nx   =    31;
  float xmin =  -5.0;
  float xmax = 305.0;

  TH2F* heff = new TH2F("heff","heff", nx , xmin , xmax , nx , xmin , xmax );

  TCanvas *ctemp = new TCanvas();
  ch->Draw("ml:mg>>heff",sel*weight);
  delete ctemp;
  heff->Scale(1.0/denom);
  heff->Scale(1000);

  int effbin = heff->FindBin(200,0);
  cout << "Efficiency for 200,0  " << heff->GetBinContent(effbin) << endl;

  effbin = heff->FindBin(200,50);
  cout << "Efficiency for 200,50 " << heff->GetBinContent(effbin) << endl;

  effbin = heff->FindBin(150,0);
  cout << "Efficiency for 150,0  " << heff->GetBinContent(effbin) << endl;

  //-------------------------------
  // find excluded points
  //-------------------------------

  //TFile *xsecfile = TFile::Open("reference_xSec.root");
  TFile *xsecfile = TFile::Open(Form("%s.root",xsectype));
  TH1F* refxsec   = (TH1F*) xsecfile->Get(xsechist);

  //plotProjections(hexcl,refxsec);

  TH2F* hexcluded     = new TH2F("hexcluded"     ,"hexcluded"     , nx , xmin , xmax , nx , xmin , xmax );
  TH2F* hexcluded_exp = new TH2F("hexcluded_exp" ,"hexcluded_exp" , nx , xmin , xmax , nx , xmin , xmax );

  for( unsigned int ibin = 1 ; ibin <= nx ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= nx ; jbin++ ){

      float xsecul     = hexcl->GetBinContent(ibin,jbin);
      float xsecul_exp = hexp->GetBinContent(ibin,jbin);

      if( xsecul < 1.e-10 ) continue;

      float mg = hexcluded->GetXaxis()->GetBinCenter(ibin);
      float ml = hexcluded->GetYaxis()->GetBinCenter(jbin);

      int   bin  = refxsec->FindBin(mg);
      float xsec = refxsec->GetBinContent(bin);
      
      //cout << ibin << " " << jbin << " " << mg << " " << ml << " " << xsecul << " " << xsec << " " << (xsec > xsecul) << endl;

      hexcluded->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul )     hexcluded->SetBinContent(ibin,jbin,1);

      hexcluded_exp->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_exp ) hexcluded_exp->SetBinContent(ibin,jbin,1);

      //cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
    }
  }

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);

  //-------------------------------
  // draw efficiency
  //-------------------------------

  TCanvas *can = new TCanvas("can","",600,600);
  can->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);

  if( TString(sample).Contains("gmsb") && smooth ) smoothHist( heff );

  heff->GetYaxis()->SetRangeUser(-5,305);
  heff->GetXaxis()->SetRangeUser(95,305);
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetTitle(ytitle);
  heff->GetXaxis()->SetTitle(xtitle);
  heff->GetZaxis()->SetTitle("efficiency #times acceptance (10^{-3})");
  heff->Draw("colz");
  //heff->GetYaxis()->SetRangeUser(ymin,1200);

  t->SetTextSize(0.04);

  //t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.85,title);
  t->DrawLatex(0.2,0.78,"E_{T}^{miss} > 80 GeV");

  t->SetTextSize(0.04);
  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 8 TeV, L_{int} = 9.2 fb^{-1}");

  //-------------------------------
  // cross section limit
  //-------------------------------

  TCanvas *can2 = new TCanvas("can2","",600,600);
  can2->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);

  if( TString(sample).Contains("gmsb") && smooth ) smoothHist( hexcl );

  hexcl->GetYaxis()->SetRangeUser(-5,305);
  hexcl->GetXaxis()->SetRangeUser(95,305);
  gPad->SetLogz();
  hexcl->GetXaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetLabelSize(0.035);
  hexcl->GetZaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetTitle(ytitle);
  hexcl->GetXaxis()->SetTitle(xtitle);
  hexcl->GetZaxis()->SetTitle("95% CL UL #sigma #times BR [pb]");
  gStyle->SetPaintTextFormat(".2f");
  hexcl->Draw("colz");

  TGraph* grobs  = getGraph_observed();
  TGraph* grexp  = getGraph_expected();
  TGraph* gr2011 = getGraph_2011();

  grobs->SetLineWidth(5);
  grexp->SetLineWidth(5);
  gr2011->SetLineWidth(3);

  grexp->SetLineStyle(9);

  grobs->SetLineColor(1);
  grexp->SetLineColor(kOrange+3);
  gr2011->SetLineColor(2);

  gr2011->Draw("l");
  grexp->Draw("l");
  grobs->Draw("l");
  
  TLegend *leg = new TLegend(0.2,0.6,0.55,0.75);
  leg->AddEntry(grobs,  "observed","l");
  leg->AddEntry(grexp,  "expected","l");
  leg->AddEntry(gr2011, "2011 observed","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->Draw();

  t->SetTextSize(0.04);
  //t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.85,title);

  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 8 TeV, L_{int} = 9.2 fb^{-1}");

  if( print ){
    can->Print(Form("cards/%s/plots/%s_eff.pdf"   ,version.c_str(),sample));
    can->Print(Form("cards/%s/plots/%s_eff.C"     ,version.c_str(),sample));
    can2->Print(Form("cards/%s/plots/%s_xsec.pdf" ,version.c_str(),sample));
    can2->Print(Form("cards/%s/plots/%s_xsec.C"   ,version.c_str(),sample));
  }

  TFile* fout = TFile::Open(Form("cards/%s/limit.root",version.c_str()),"RECREATE");
  fout->cd();
  hexcl->Write();
  grobs->Write();
  grexp->Write();
  fout->Close();

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  t->SetTextSize(0.07);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded->GetXaxis()->SetTitle(xtitle);
  hexcluded->GetYaxis()->SetTitle(ytitle);
  hexcluded->GetXaxis()->SetRangeUser(150,350);
  hexcluded->GetYaxis()->SetRangeUser(0,100);
  hexcluded->Draw("colz");
  grobs->Draw("lp");
  t->DrawLatex(0.3,0.8,"observed");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded_exp->GetXaxis()->SetTitle(xtitle);
  hexcluded_exp->GetYaxis()->SetTitle(ytitle);
  hexcluded_exp->GetXaxis()->SetRangeUser(150,350);
  hexcluded_exp->GetYaxis()->SetRangeUser(0,100);
  hexcluded_exp->Draw("colz");
  grexp->Draw("lp");
  t->DrawLatex(0.3,0.8,"expected");

  // TCanvas *c3 = new TCanvas("c3","c3",600,600);
  // c3->cd();
  // hexcluded_exp->Draw("colz");

  if( print ){
    c2->Print(Form("cards/%s/plots/%s_points.pdf",version.c_str(),sample));
    c2->Print(Form("cards/%s/plots/%s_points.C",version.c_str(),sample));
  }


}
