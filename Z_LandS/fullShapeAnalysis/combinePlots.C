#include "Utils/SMS_utils.C"
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

TGraph* getGraph(bool do3jets,string type){

  float x[5];
  float y[5];
  int npoints = -1;

  if( !do3jets && type == "nom" ){
    x[0] =  700;  y[0] = 100;
    x[1] =  850;  y[1] = 300;
    x[2] =  925;  y[2] = 500;
    x[3] =  975;  y[3] = 850;
    x[4] =  925;  y[4] = 925;
    x[6] =  100;  y[6] = 100;
    npoints = 6;
  }
  else if( !do3jets && type == "down" ){
    x[0] = 500;   y[0] = 100;
    x[1] = 720;   y[1] = 300;
    x[2] = 825;   y[2] = 550;
    x[3] = 825;   y[3] = 825;
    x[4] = 100;   y[4] = 100;
    npoints = 5;
  }
  else if( !do3jets && type == "up" ){
    x[0] = 800;   y[0] =  100;
    x[1] = 950;   y[1] =  200;
    x[2] = 1025;  y[2] =  400;
    x[3] = 1100;  y[3] = 1000;
    x[4] = 1000;  y[4] = 1000;
    x[5] =  100;  y[5] =  100;
    npoints = 6;
  }
  else if( do3jets && type == "nom" ){
    x[0] = 700;   y[0] = 100;
    x[1] = 850;   y[1] = 300;
    x[2] = 925;   y[2] = 500;
    x[3] = 925;   y[3] = 900;
    npoints = 4;
  }
  else if( do3jets && type == "down" ){
    x[0] = 525;   y[0] = 100;
    x[1] = 720;   y[1] = 300;
    x[2] = 800;   y[2] = 500;
    x[3] = 800;   y[3] = 800;
    npoints = 4;
  }

  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}

TGraph* getGraph_T5zzh(string type){

  float x[5];
  float y[5];
  int npoints = -1;

  if( type == "nom" ){
    x[0] =  837;  y[0] =  50;
    x[1] =  837;  y[1] = 300;
    x[2] =  800;  y[2] = 350;
    x[3] =  625;  y[3] = 350;
    x[4] =  625;  y[4] = 550;
    npoints = 5;
  }
  else if( type == "down" ){
    x[0] = 500;   y[0] = 100;
    x[1] = 720;   y[1] = 300;
    x[2] = 825;   y[2] = 550;
    x[3] = 825;   y[3] = 800;
    npoints = 4;
  }
  else if( type == "up" ){
    x[0] = 800;   y[0] =  100;
    x[1] = 950;   y[1] =  200;
    x[2] = 1025;  y[2] =  400;
    x[3] = 1100;  y[3] = 1000;
    x[4] = 1000;  y[4] = 1000;
    npoints = 5;
  }

  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}

TGraph* getGraph_T5zz(string type){

  float x[6];
  float y[6];
  int npoints = -1;

  if( type == "nom" ){
    x[0] =  925;  y[0] =  50;
    x[1] =  900;  y[1] = 325;
    x[2] =  850;  y[2] = 400;
    x[3] =  825;  y[3] = 425;
    x[4] =  625;  y[4] = 425;
    x[5] =  625;  y[5] = 537.5;
    x[6] =  125;  y[6] =  37.5;
    npoints = 7;
  }
  else if( type == "up" ){
    x[0] =  1050;  y[0] =  50;
    x[1] =  1050;  y[1] = 375;
    x[2] =  1000;  y[2] = 475;
    x[3] =   950;  y[3] = 525;
    x[4] =   725;  y[4] = 525;
    x[5] =   725;  y[5] = 637.5;
    x[6] =  125;   y[6] =  37.5;
    npoints = 7;
  }
  else if( type == "down" ){
    x[0] =  775;  y[0] =  50;
    x[1] =  775;  y[1] = 150;
    x[2] =  725;  y[2] = 300;
    x[3] =  525;  y[3] = 300;
    x[4] =  525;  y[4] = 437.5;
    x[5] =  125;  y[5] =  37.5;
    npoints = 6;
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

void combinePlots(bool print = false){
  
  char* version        = "V00-01-00";
  char* sample         = "T5zz";
  bool  do3jets        = false;
  char* title          = "m(#tilde{q}) >> m(#tilde{g}), x = 0.5";
  float dm             = 182.0;

  // char* version        = "V00-01-01";
  // char* sample         = "T5zz";
  // bool  do3jets        = true;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g}), x = 0.5";

  // char* version        = "V00-01-02";
  // char* sample         = "T5zzl";
  // bool  do3jets        = false;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g}), x = 0.75";

  // char* version        = "V00-01-03";
  // char* sample         = "T5zzl";
  // bool  do3jets        = true;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g}), x = 0.75";
  
  // char* version        = "V00-01-04";
  // char* sample         = "T5zzgmsb";
  // bool  do3jets        = false;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g})";
  // float dm             = 0.0;

  // char* version        = "V00-01-05";
  // char* sample         = "T5zzgmsb";
  // bool  do3jets        = true;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g})";

  // char* version        = "V00-01-07";
  // char* sample         = "T5zzgmsb_hadoop";
  // bool  do3jets        = false;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g})";

  // char* version        = "V00-01-08";
  // char* sample         = "T5zzh";
  // bool  do3jets        = false;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g}), x = 0.25";

  char* njets          = "n_{jets} #geq 2";
  if( do3jets )  njets = "n_{jets} #geq 3";

  float ymin = 50.;
  if( TString(sample).Contains("gmsb") ) ymin = 100.;

  TFile *file = TFile::Open(Form("cards/%s/observed_limit.root",version));
  TH2F* hexcl = (TH2F*) file->Get("hexcl");
  //TH2F* hexcl = (TH2F*) file->Get("hexp");
  //cout << "USING EXPECTED LIMIT!!!!!!!!!!!!" << endl;

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // calculate efficiency
  //-------------------------------

  char* babyversion = "V00-02-04";
  if( TString(sample).Contains("hadoop") ) babyversion = "V00-02-05";
  if( TString(sample).Contains("T5zzh") )  babyversion = "V00-02-05";

  TChain *ch = new TChain("T1");
  //ch->Add(Form("output/V00-02-04/%s_baby.root",sample));
  ch->Add(Form("output/%s/%s_baby.root",babyversion,sample));

  TCut zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets >= 2");
  TCut njets3("njets >= 3");
  TCut sf("leptype==0||leptype==1");
  TCut met100("pfmet>100");

  TCut sel  = zmass + njets2 + sf + met100;
  if( do3jets ) sel  = zmass + njets3 + sf + met100;

  cout << "Using selection: " << sel.GetTitle() << endl;

  TH2F* heff = new TH2F("heff","heff", 48,0,1200,48,0,1200);
  TCanvas *ctemp = new TCanvas();
  ch->Draw("ml:mg>>heff",sel);
  delete ctemp;
  heff->Scale(1./20000.);

  int bin = heff->FindBin(600,200);
  cout << "Efficiency for 600,200 " << heff->GetBinContent(bin) << endl;

  //-------------------------------
  // find excluded points
  //-------------------------------

  //TFile *xsecfile = TFile::Open("reference_xSec.root");
  TFile *xsecfile = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("gluino");

  TH2F* hexcluded   = new TH2F("hexcluded","hexcluded", 48,0,1200,48,0,1200);
  TH2F* hexcluded13 = new TH2F("hexcluded13","hexcluded13", 48,0,1200,48,0,1200);
  TH2F* hexcluded3  = new TH2F("hexcluded3","hexcluded3", 48,0,1200,48,0,1200);
  
  for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

      float xsecul = hexcl->GetBinContent(ibin,jbin);

      if( xsecul < 1.e-10 ) continue;

      float mg = hexcluded->GetXaxis()->GetBinCenter(ibin)-12.5;
      float ml = hexcluded->GetYaxis()->GetBinCenter(jbin)-12.5;

      int   bin = refxsec->FindBin(mg);
      float xsec = refxsec->GetBinContent(bin);

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

  if( TString(sample).Contains("gmsb") ) smoothHist( heff );

  heff->GetYaxis()->SetRangeUser(ymin,1200);
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetTitle("#chi^{0}_{1} mass [GeV]");
  heff->GetXaxis()->SetTitle("gluino mass [GeV]");
  heff->GetZaxis()->SetTitle("A #times #varepsilon (#geq1 Z(ll))");
  heff->Draw("colz");
  heff->GetYaxis()->SetRangeUser(ymin,1200);

  t->SetTextSize(0.04);

  if(TString(sample).Contains("gmsb") )  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{1}^{0}, #chi_{1}^{0} #rightarrow Z G");
  else                                   t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");

  t->DrawLatex(0.2,0.77,title);
  t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 100 GeV");
  t->DrawLatex(0.2,0.65,njets);
  t->DrawLatex(0.2,0.55,"E_{T}^{miss} templates");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

  if(TString(sample).Contains("gmsb") )   line.DrawLine(100,100,1200,1200);
  else                                    line.DrawLine(50+dm,50,1200,1200-dm);

  //-------------------------------
  // cross section limit
  //-------------------------------

  can->cd(2);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);

  if( TString(sample).Contains("gmsb") ) smoothHist( hexcl );

  hexcl->GetYaxis()->SetRangeUser(ymin,1200);
  //hexcl->GetXaxis()->SetRangeUser(0,950);
  gPad->SetLogz();
  hexcl->GetXaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetLabelSize(0.035);
  hexcl->GetZaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetTitle("#chi^{0}_{1} mass [GeV]");
  hexcl->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcl->GetZaxis()->SetTitle("95% CL upper limit on #sigma [pb]");
  hexcl->Draw("colz");
  hexcl->SetMinimum(0.001);
  hexcl->SetMaximum(10);
  hexcl->GetYaxis()->SetRangeUser(ymin,1200);

  TGraph* gr_excl;      
  TGraph* gr_excl_down;
  TGraph* gr_excl_up;   
  
  if( TString(sample).Contains("gmsb") ) {
    gr_excl      = getGraph(do3jets,"nom");
    gr_excl_down = getGraph(do3jets,"down");
    gr_excl_up   = getGraph(do3jets,"up");
  }
  else if( TString(sample).Contains("T5zz") ) {
    gr_excl      = getGraph_T5zz("nom");
    gr_excl_down = getGraph_T5zz("down");
    gr_excl_up   = getGraph_T5zz("up");
  }
  else{
    gr_excl      = getRefXsecGraph(hexcl, "T5zz", 1.0);
    gr_excl_down = getRefXsecGraph(hexcl, "T5zz", 1./3.);
    gr_excl_up   = getRefXsecGraph(hexcl, "T5zz", 3.);
  }

  gr_excl->SetLineWidth(2.5);
  gr_excl_up->SetLineWidth(2.5);
  gr_excl_down->SetLineWidth(2.5);
  gr_excl_up->SetLineStyle(2);
  gr_excl_down->SetLineStyle(3);
  gr_excl->Draw("same");
  gr_excl_up->Draw("same");
  gr_excl_down->Draw("same");

  if(TString(sample).Contains("gmsb") )   line.DrawLine(100,100,1200,1200);
  else                                    line.DrawLine(50+dm,50,1200,1200-dm);

  TLegend *leg = new TLegend(0.2,0.53,0.45,0.67);
  leg->AddEntry(gr_excl,     "#sigma^{NLO-QCD}","l");
  leg->AddEntry(gr_excl_up,  "3 #times #sigma^{NLO-QCD}","l");
  leg->AddEntry(gr_excl_down,"1/3 #times #sigma^{NLO-QCD}","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw();
  
  t->SetTextSize(0.04);
  if(TString(sample).Contains("gmsb") )  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{1}^{0}, #chi_{1}^{0} #rightarrow Z G");
  else                                   t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
  t->DrawLatex(0.2,0.77,title);
  t->DrawLatex(0.2,0.71,njets);
  t->DrawLatex(0.2,0.47,"E_{T}^{miss} templates");
  t->SetTextSize(0.035);
  //t->DrawLatex(0.18,0.92,"CMS                     #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

  if( print ){
    can->Print(Form("cards/%s/plots/SMS.eps",version));
    can->Print(Form("cards/%s/plots/SMS.pdf",version));
    can->Print(Form("cards/%s/plots/SMS.png",version));
    can->Print(Form("cards/%s/plots/SMS.C",version));

    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS.eps cards/%s/plots/SMS_ppt.pdf",version,version));
  }


  TFile* fout = TFile::Open(Form("cards/%s/limit.root",version),"RECREATE");
  fout->cd();
  hexcl->Write();
  gr_excl->Write();
  fout->Close();

  TCanvas *c2 = new TCanvas("c2","c2",1500,500);
  c2->Divide(3,1);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded13->Draw("colz");
  gr_excl_down->Draw();

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded->Draw("colz");
  gr_excl->Draw();

  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hexcluded3->Draw("colz");
  gr_excl_up->Draw();

}
