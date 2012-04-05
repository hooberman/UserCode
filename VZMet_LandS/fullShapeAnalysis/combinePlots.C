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

// TH2F* shiftHist(TH2F* hin){

//   TH2F* hout = new TH2F(Form("%s_out",hin->GetName()),Form("%s_out",hin->GetName()), 48,0-12.5,1200-12.5,48,0-12.5,1200-12.5);

//   for(int ibin = 1 ; ibin <= 48 ; ibin++ ){
//     for(int jbin = 1 ; jbin <= 48 ; jbin++ ){
//       hout->SetBinContent(ibin,jbin,hin->GetBinContent(ibin,jbin));
//     }
//   }

//   return hout;
// }


TGraph* getGraph_WZ(string type){

  float x[5];
  float y[5];
  int npoints = -1;

  if( type == "nom" ){
    x[0] =  212.5;  y[0] = -12.5;
    x[1] =  212.5;  y[1] =  37.5;
    x[2] =  137.5;  y[2] =  37.5;
    x[3] =  137.5;  y[3] = -12.5;
    npoints = 4;
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

void combinePlots(string version , bool print = false){

  bool smooth = false;

  char* sample;
  char* title;
  char* xsectype;
  float denom;

  if( version == "V00-02-02" ){
    sample   = (char*) "wzsms";
    title    = (char*) "pp#rightarrow #chi^{#pm}#chi^{0} #rightarrow WZ + E_{T}^{miss}";
    xsectype = (char*) "C1N2";
    denom    = 100000.0;
  } 

  else if( version == "V00-02-03" ){
    sample   = (char*) "zzsms";
    title    = (char*) "pp#rightarrow #chi^{0}#chi^{0} #rightarrow ZZ + E_{T}^{miss}";
    xsectype = (char*) "N1N2";
    denom    = 52600.0;
  } 
  
  //float ymin = 0.;

  TFile *file = TFile::Open(Form("cards/%s/observed_limit.root",version.c_str()));

  TH2F* hexcl      = (TH2F*) file->Get("hexcl");
  hexcl->Scale(1000.0); // pb --> fb

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // calculate efficiency
  //-------------------------------

  char* babyversion = (char*) "V00-02-14";

  TChain *ch = new TChain("T1");
  ch->Add(Form("output/%s/%s_baby.root",babyversion,sample));

  TCut zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets >= 2");
  TCut bveto("nbvz==0");
  TCut mjj("mjj>70 && mjj<110");
  TCut nlep2("nlep==2");
  TCut rho("rho>0 && rho<40");
  TCut sf("leptype==0||leptype==1");
  TCut met60("pfmet>60");
  TCut weight("btagweight * davtxweight * trgeff");

  TCut sel  = zmass + njets2 + bveto + mjj + nlep2 + sf + met60;

  cout << "Using selection: " << sel.GetTitle() << endl;

  TH2F* heff = new TH2F("heff","heff", 21 , -12.5 , 512.5 , 21 , -12.5 , 512.5);

  TCanvas *ctemp = new TCanvas();
  ch->Draw("ml:mg>>heff",sel);
  delete ctemp;
  heff->Scale(1.0/denom);
  heff->Scale(100);

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
      if( (1./3.) * xsec > xsecul )   hexcluded13->SetBinContent(ibin,jbin,1);

      //cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
    }
  }

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);

  //-------------------------------
  // draw efficiency
  //-------------------------------

  TCanvas *can = new TCanvas("can","",1200,600);
  can->cd();
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
  hexcl->SetMinimum(99);
  hexcl->SetMaximum(5000);
  //hexcl->GetYaxis()->SetRangeUser(ymin,1200);

  TGraph* gr_excl;      
  TGraph* gr_excl_down;
  TGraph* gr_excl_up;   
  
  if( TString(sample).Contains("wzsms") ) {
    gr_excl      = getGraph_WZ("nom");
    gr_excl_down = getGraph_WZ("down");
    gr_excl_up   = getGraph_WZ("up");
  }
  else if( TString(sample).Contains("zzsms") ) {
    // gr_excl      = getGraph_T5zzl("nom");
    // gr_excl_down = getGraph_T5zzl("down");
    // gr_excl_up   = getGraph_T5zzl("up");
  }

  gr_excl->SetLineWidth(2.5);
  gr_excl_up->SetLineWidth(2.5);
  gr_excl_down->SetLineWidth(2.5);
  gr_excl_up->SetLineStyle(2);
  gr_excl_down->SetLineStyle(3);
  gr_excl->Draw("same");
  //gr_excl_up->Draw("same");
  //gr_excl_down->Draw("same");

  TLegend *leg = new TLegend(0.2,0.53,0.55,0.67);
  leg->AddEntry(gr_excl,     "#sigma^{wino-like}","l");
  //leg->AddEntry(gr_excl_up,  "3 #times #sigma^{NLO-QCD}","l");
  //leg->AddEntry(gr_excl_down,"1/3 #times #sigma^{NLO-QCD}","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
  leg->Draw();
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.77,title);

  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}");

  if( print ){
    can->Print(Form("cards/%s/plots/SMS.eps",version.c_str()));
    can->Print(Form("cards/%s/plots/SMS.pdf",version.c_str()));
    can->Print(Form("cards/%s/plots/SMS.png",version.c_str()));
    can->Print(Form("cards/%s/plots/SMS.C",version.c_str()));

    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS.eps cards/%s/plots/SMS_ppt.pdf",version.c_str(),version.c_str()));
  }

  // TH2F* hexcluded_shifted   = shiftHist( hexcluded   );
  // TH2F* hexcluded13_shifted = shiftHist( hexcluded13 );
  // TH2F* hexcluded3_shifted  = shiftHist( hexcluded3  );

  TH2F* hexcluded_shifted   = (TH2F*) hexcluded->Clone("hexcluded_shifted");
  TH2F* hexcluded13_shifted = (TH2F*) hexcluded13->Clone("hexcluded13_shifted");
  TH2F* hexcluded3_shifted  = (TH2F*) hexcluded3->Clone("hexcluded3_shifted");

  TFile* fout = TFile::Open(Form("cards/%s/limit.root",version.c_str()),"RECREATE");
  fout->cd();
  hexcl->Write();
  gr_excl->Write();
  fout->Close();

  TCanvas *c2 = new TCanvas("c2","c2",1500,500);
  c2->Divide(3,1);

  t->SetTextSize(0.07);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded13->Draw("colz");
  hexcluded13_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded13_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded13_shifted->Draw("colz");
  gr_excl_down->Draw();
  t->DrawLatex(0.3,0.8,"1/3 #times #sigma^{NLO-QCD}");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded->Draw("colz");
  hexcluded_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded_shifted->Draw("colz");
  gr_excl->Draw();
  t->DrawLatex(0.3,0.8,"#sigma^{NLO-QCD}");

  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded3->Draw("colz");
  hexcluded3_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded3_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded3_shifted->Draw("colz");
  gr_excl_up->Draw();
  t->DrawLatex(0.3,0.8,"3 #times #sigma^{NLO-QCD}");

  if( print ){
    c2->Print(Form("cards/%s/plots/SMS_points.eps",version.c_str()));
    c2->Print(Form("cards/%s/plots/SMS_points.pdf",version.c_str()));
    c2->Print(Form("cards/%s/plots/SMS_points.png",version.c_str()));
    c2->Print(Form("cards/%s/plots/SMS_points.C",version.c_str()));

    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS_points.eps cards/%s/plots/SMS_points_ppt.pdf",version.c_str(),version.c_str()));
  }


}
