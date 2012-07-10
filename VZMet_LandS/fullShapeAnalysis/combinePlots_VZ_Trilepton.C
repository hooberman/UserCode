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

/*
TGraph* getGraph_WZ(string type){

  float x[6];
  float y[6];
  int npoints = -1;

  if( type == "nom" ){

    //-------------------
    // combination
    //-------------------

    // x[0] =  212.5;  y[0] = -12.5;
    // x[1] =  212.5;  y[1] =  62.5;
    // x[2] =  162.5;  y[2] =  62.5;
    // x[3] =  162.5;  y[3] =  37.5;
    // x[4] =  112.5;  y[4] =  37.5;
    // x[5] =  112.5;  y[5] = -12.5;
    // npoints = 6;

    //-------------------
    // trileptons only
    //-------------------

    x[0] =  162.5;  y[0] = -12.5;
    x[1] =  162.5;  y[1] =  12.5;
    x[2] =  187.5;  y[2] =  12.5;
    x[3] =  187.5;  y[3] =  37.5;
    x[4] =  112.5;  y[4] =  37.5;
    x[5] =  112.5;  y[5] = -12.5;
    npoints = 6;
  }
  else if( type == "down" ){
    x[0] = 600;   y[0] =  50;
    x[1] = 600;   y[1] = 150;
    x[2] = 525;   y[2] = 200;
    x[3] = 475;   y[3] = 300;
    x[4] = 525;   y[4] = 450;
    npoints = 5;
  }
  else if( type == "up" ){
    x[0] = 1000;  y[0] =   50;
    x[1] = 1000;  y[1] =  300;
    x[2] =  950;  y[2] =  475;
    x[3] = 712.5; y[3] =  475;
    x[4] = 712.5; y[4] =  637.5;
    npoints = 5;
  }

  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}
*/

TGraph* getGraph_Trilepton(){

  float x[6];
  float y[6];
  int npoints = -1;

  x[0] =  187.5;  y[0] = -12.5;
  x[1] =  187.5;  y[1] =  25.0;
  x[2] =  175.0;  y[2] =  37.5;
  x[3] =  127.5;  y[3] =  37.5;
  x[4] =   87.5;  y[4] =  12.5;
  x[5] =   87.5;  y[5] = -12.5;
  npoints = 6;

  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}

TGraph* getGraph_Combo(){

  float x[6];
  float y[6];
  int npoints = -1;

  x[0] =  237.5;  y[0] = -12.5;
  x[1] =  212.5;  y[1] =   75;
  x[2] =  187.5;  y[2] =   75;
  x[3] =   87.5;  y[3] =  12.5;
  x[4] =   87.5;  y[4] = -12.5;

  npoints = 5;

  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}

TGraph* getGraph_VZMet(){

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

void combinePlots_VZ_Trilepton(bool print = false){

  ifstream ifile("T1ChiWZ_Combo_UL.txt");
  //ifstream ifile("T1ChiWZ_Trileptons_UL.txt");
  //ifstream ifile("T1ChiWZ_LLJJ_UL.txt");

  char* rootfilename = "T1ChiWZ_LLJJ_UL.txt";

  bool smooth = false;

  char* sample   = (char*) "wzsms";
  char* title    = (char*) "pp#rightarrow #chi^{#pm}#chi^{0} #rightarrow WZ + E_{T}^{miss}";
  char* xsectype = (char*) "C1N2";

  TLatex *t = new TLatex();
  t->SetNDC();

  TH2F* hexp   = new TH2F("hexp","hexp"    , 21 , -12.5 , 512.5 , 21 , -12.5 , 512.5 );
  TH2F* hexcl   = new TH2F("hexcl","hhexcl"    , 21 , -12.5 , 512.5 , 21 , -12.5 , 512.5 );

  int m1;
  int m2;
  float obs;
  float exp;

  while( !ifile.eof() ){
    ifile >> m1 >> m2 >> exp >> obs;
    cout << m1 << " " << m2 << " " << exp << " " << obs << endl;

    hexcl->Fill(m1,m2,obs);
    hexp->Fill(m1,m2,exp);
    
  }


  gStyle->SetPaintTextFormat(".0f");

  // TCanvas *c1 = new TCanvas("c1","",2400,1200);
  // c1->Divide(2,1);

  // c1->cd(1);
  // hexp->Draw("colz");
  // hexp->Draw("sametext");

  // c1->cd(2);
  // hexcl->Draw("colz");
  // hexcl->Draw("sametext");


  //-------------------------------
  // find excluded points
  //-------------------------------

  //TFile *xsecfile = TFile::Open("reference_xSec.root");
  TFile *xsecfile = TFile::Open(Form("%s_referencexSec.root",xsectype));
  TH1F* refxsec   = (TH1F*) xsecfile->Get(xsectype);

  TH2F* hexcluded   = new TH2F("hexcluded","hexcluded"    , 21 , -12.5 , 512.5 , 21 , -12.5 , 512.5 );
  TH2F* hexcluded13 = new TH2F("hexcluded13","hexcluded13", 21 , -12.5 , 512.5 , 21 , -12.5 , 512.5 );
  TH2F* hexcluded3  = new TH2F("hexcluded3","hexcluded3"  , 21 , -12.5 , 512.5 , 21 , -12.5 , 512.5 );
  
  for( unsigned int ibin = 1 ; ibin <= 21 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 21 ; jbin++ ){

      float xsecul = hexcl->GetBinContent(ibin,jbin);

      if( xsecul < 1.e-10 ) continue;

      float mg = hexcluded->GetXaxis()->GetBinCenter(ibin)-12.5;
      float ml = hexcluded->GetYaxis()->GetBinCenter(jbin)-12.5;

      int   bin  = refxsec->FindBin(mg);
      float xsec = refxsec->GetBinContent(bin);
      
      //cout << ibin << " " << jbin << " " << mg << " " << ml << " " << xsecul << " " << xsec << " " << (xsec > xsecul) << endl;

      hexcluded->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul )   hexcluded->SetBinContent(ibin,jbin,1);

      hexcluded3->SetBinContent(ibin,jbin,0);
      if( 3 * xsec > xsecul )   hexcluded3->SetBinContent(ibin,jbin,1);

      hexcluded13->SetBinContent(ibin,jbin,0);
      if( (1./2.) * xsec > xsecul )   hexcluded13->SetBinContent(ibin,jbin,1);

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

  /*
  can->Divide(2,1);

  can->cd(1);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);

  if( TString(sample).Contains("gmsb") && smooth ) smoothHist( heff );

  //heff->GetYaxis()->SetRangeUser(ymin,1200);
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  //heff->SetMaximum(0.35);
  heff->GetYaxis()->SetTitle("LSP mass [GeV]");
  heff->GetXaxis()->SetTitle("#chi mass [GeV]");
  heff->GetZaxis()->SetTitle("A #times #varepsilon (%)");
  heff->Draw("colz");
  //heff->GetYaxis()->SetRangeUser(ymin,1200);

  t->SetTextSize(0.04);

  t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.77,title);
  t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 60 GeV");

  t->SetTextSize(0.04);
  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}");

  //-------------------------------
  // cross section limit
  //-------------------------------

  can->cd(2);
  */
  can->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);

  if( TString(sample).Contains("gmsb") && smooth ) smoothHist( hexcl );

  //hexcl->GetYaxis()->SetRangeUser(ymin,1200);
  //hexcl->GetXaxis()->SetRangeUser(0,950);
  gPad->SetLogz();
  hexcl->GetXaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetLabelSize(0.035);
  hexcl->GetZaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetTitle("LSP mass [GeV]");
  hexcl->GetXaxis()->SetTitle("#chi mass [GeV]");
  hexcl->GetZaxis()->SetTitle("95% CL upper limit on #sigma [fb]");
  hexcl->Draw("colz");
  hexcl->SetMinimum(50);
  hexcl->SetMaximum(5000);
  //hexcl->GetYaxis()->SetRangeUser(ymin,1200);

  // TGraph* gr_excl      = getGraph_WZ("nom");
  // TGraph* gr_excl_down = getGraph_WZ("down");
  // TGraph* gr_excl_up   = getGraph_WZ("up");

  TGraph* gr_vzmet   = getGraph_VZMet();
  TGraph* gr_tri     = getGraph_Trilepton();
  TGraph* gr_combo   = getGraph_Combo();
  
  if( TString(sample).Contains("wzsms") ) {
    
    gr_vzmet->SetLineWidth(3);
    gr_combo->SetLineWidth(6);
    gr_tri->SetLineWidth(3);
    gr_combo->SetLineStyle(1);
    gr_combo->SetLineColor(1);
    gr_tri->SetLineStyle(2);
    gr_vzmet->SetLineStyle(2);
    gr_tri->SetLineColor(2);
    gr_vzmet->SetLineColor(4);
    gr_combo->Draw("same");

    gr_vzmet->Draw("same");
    gr_tri->Draw("same");

    TLegend *leg = new TLegend(0.2,0.53,0.55,0.75);
    //leg->AddEntry(gr_vzmet,     "#sigma^{wino-like}","l");
    leg->AddEntry(gr_vzmet  , "VZ+E_{T}^{miss}","l");
    leg->AddEntry(gr_tri    , "trilepton","l");
    leg->AddEntry(gr_combo  , "combination","l");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->Draw();
  }
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.77,title);

  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}");

  if( print ){
    // can->Print(Form("cards/%s/plots/SMS.eps",version.c_str()));
    // can->Print(Form("cards/%s/plots/SMS.pdf",version.c_str()));
    // can->Print(Form("cards/%s/plots/SMS.png",version.c_str()));
    // can->Print(Form("cards/%s/plots/SMS.C",version.c_str()));

    // gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS.eps cards/%s/plots/SMS_ppt.pdf",version.c_str(),version.c_str()));
  }

  // TH2F* hexcluded_shifted   = shiftHist( hexcluded   );
  // TH2F* hexcluded13_shifted = shiftHist( hexcluded13 );
  // TH2F* hexcluded3_shifted  = shiftHist( hexcluded3  );

  TH2F* hexcluded_shifted   = (TH2F*) hexcluded->Clone("hexcluded_shifted");
  TH2F* hexcluded13_shifted = (TH2F*) hexcluded13->Clone("hexcluded13_shifted");
  TH2F* hexcluded3_shifted  = (TH2F*) hexcluded3->Clone("hexcluded3_shifted");

  // TFile* fout = TFile::Open(Form("cards/%s/limit.root",version.c_str()),"RECREATE");
  // fout->cd();
  // hexcl->Write();
  // gr_vzmet->Write();
  // fout->Close();

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  t->SetTextSize(0.07);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded13->Draw("colz");
  hexcluded13_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded13_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded13_shifted->Draw("colz");
  //gr_tri->Draw();
  t->DrawLatex(0.3,0.8,"#sigma^{higgsino}");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded->Draw("colz");
  hexcluded_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded_shifted->Draw("colz");
  gr_tri->Draw();
  gr_vzmet->Draw();
  gr_combo->Draw();
  t->DrawLatex(0.3,0.8,"#sigma^{wino}");

  if( print ){
    // c2->Print(Form("cards/%s/plots/SMS_points.eps",version.c_str()));
    // c2->Print(Form("cards/%s/plots/SMS_points.pdf",version.c_str()));
    // c2->Print(Form("cards/%s/plots/SMS_points.png",version.c_str()));
    // c2->Print(Form("cards/%s/plots/SMS_points.C",version.c_str()));

    // gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS_points.eps cards/%s/plots/SMS_points_ppt.pdf",version.c_str(),version.c_str()));
  }


  TFile *fout = TFile::Open(rootfilename,"RECREATE");
  fout->cd();
  hexcl->Write();
  hexp->Write();
  fout->Close();

}
