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

void plotProjections( TH2F* h ){

  TH1F* hproj = (TH1F*) h->ProjectionX("hproj",1,1);


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
  TCut met100("pfmet>100");
  TCut weight("vtxweight * trgeff");

  TCut sel  = pt2020 + zmass + njets2 + bveto + mjj + nlep2 + sf + met100;

  cout << "Using selection: " << sel.GetTitle() << endl;

  TH2F* heff = new TH2F("heff","heff", 31 , -5.0 , 305.0 , 31 , -5.0 , 305.0 );

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

  int   nx   =    31;
  float xmin =  -5.0;
  float xmax = 305.0;

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
  heff->GetYaxis()->SetTitle(ytitle);
  heff->GetXaxis()->SetTitle(xtitle);
  heff->GetZaxis()->SetTitle("efficiency (%)");
  heff->Draw("colz");
  //heff->GetYaxis()->SetRangeUser(ymin,1200);

  t->SetTextSize(0.04);

  //t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.77,title);
  t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 100 GeV");

  t->SetTextSize(0.04);
  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 8 TeV, L_{int} = 9.2 fb^{-1}");

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
  hexcl->GetYaxis()->SetTitle(ytitle);
  hexcl->GetXaxis()->SetTitle(xtitle);
  hexcl->GetZaxis()->SetTitle("95% CL upper limit on #sigma [pb]");
  hexcl->Draw("colz");
  //hexcl->SetMinimum(50);
  //hexcl->SetMaximum(5000);
  //hexcl->GetYaxis()->SetRangeUser(ymin,1200);

  /*
  TGraph* gr_excl      = getGraph_WZ("nom");
  TGraph* gr_excl_down = getGraph_WZ("down");
  TGraph* gr_excl_up   = getGraph_WZ("up");
  
  if( TString(sample).Contains("wzsms") ) {
    
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
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.06);
    leg->Draw();
  }
  */

  t->SetTextSize(0.04);
  //t->DrawLatex(0.2,0.83,"E_{T}^{miss} templates");
  t->DrawLatex(0.2,0.77,title);

  t->DrawLatex(0.15,0.93,"CMS Preliminary  #sqrt{s} = 8 TeV, L_{int} = 9.2 fb^{-1}");

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

  // TH2F* hexcluded_shifted   = (TH2F*) hexcluded->Clone("hexcluded_shifted");
  // TH2F* hexcluded13_shifted = (TH2F*) hexcluded13->Clone("hexcluded13_shifted");
  // TH2F* hexcluded3_shifted  = (TH2F*) hexcluded3->Clone("hexcluded3_shifted");

  TFile* fout = TFile::Open(Form("cards/%s/limit.root",version.c_str()),"RECREATE");
  fout->cd();
  hexcl->Write();
  //gr_excl->Write();
  fout->Close();

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  t->SetTextSize(0.07);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded->GetXaxis()->SetTitle(xtitle);
  hexcluded->GetYaxis()->SetTitle(ytitle);
  hexcluded->Draw("colz");
  //gr_excl_down->Draw();
  t->DrawLatex(0.3,0.8,"observed");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded_exp->GetXaxis()->SetTitle(xtitle);
  hexcluded_exp->GetYaxis()->SetTitle(ytitle);
  hexcluded_exp->Draw("colz");
  //gr_excl_down->Draw();
  t->DrawLatex(0.3,0.8,"expected");


  if( print ){
    c2->Print(Form("cards/%s/plots/SMS_points.eps",version.c_str()));
    c2->Print(Form("cards/%s/plots/SMS_points.pdf",version.c_str()));
    c2->Print(Form("cards/%s/plots/SMS_points.png",version.c_str()));
    c2->Print(Form("cards/%s/plots/SMS_points.C",version.c_str()));

    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS_points.eps cards/%s/plots/SMS_points_ppt.pdf",version.c_str(),version.c_str()));
  }


}
