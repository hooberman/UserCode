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

  float x[4];
  float y[4];
  int npoints = -1;

  if( !do3jets && type == "nom" ){
    x[0] = 700;   y[0] = 100;
    x[1] = 850;   y[1] = 300;
    x[2] = 925;   y[2] = 500;
    x[3] = 925;   y[3] = 900;
    npoints = 4;
  }
  else if( !do3jets && type == "down" ){
    x[0] = 500;   y[0] = 100;
    x[1] = 720;   y[1] = 300;
    x[2] = 800;   y[2] = 500;
    x[3] = 800;   y[3] = 800;
    npoints = 4;
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

void combinePlots(bool print = false){

  // char* version        = "V00-01-00";
  // char* sample         = "T5zz";
  // bool  do3jets        = false;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g}), x = 0.5";

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
  
  char* version        = "V00-01-04";
  char* sample         = "T5zzgmsb";
  bool  do3jets        = false;
  char* title          = "m(#tilde{q}) >> m(#tilde{g})";

  // char* version        = "V00-01-05";
  // char* sample         = "T5zzgmsb";
  // bool  do3jets        = true;
  // char* title          = "m(#tilde{q}) >> m(#tilde{g})";

  char* njets          = "n_{jets} #geq 2";
  if( do3jets )  njets = "n_{jets} #geq 3";

  TFile *file = TFile::Open(Form("cards/%s/observed_limit.root",version));
  TH2F* hexcl = (TH2F*) file->Get("hexcl");
  

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // calculate efficiency
  //-------------------------------

  TChain *ch = new TChain("T1");
  ch->Add(Form("output/V00-02-04/%s_baby.root",sample));

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

  TFile *xsecfile = TFile::Open("reference_xSec.root");
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



  //-------------------------------
  // draw efficiency
  //-------------------------------

  TCanvas *can = new TCanvas("can","",1200,600);
  can->cd();
  can->Divide(2,1);

  can->cd(1);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  //heff->GetXaxis()->SetRangeUser(0,950);
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  heff->GetXaxis()->SetTitle("gluino mass (GeV)");
  heff->GetZaxis()->SetTitle("efficiency");
  heff->Draw("colz");
  
  t->SetTextSize(0.04);

  if(TString(sample).Contains("gmsb") )  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{1}^{0}, #chi_{1}^{0} #rightarrow Z G");
  else                                   t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");

  t->DrawLatex(0.2,0.77,title);
  t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 100 GeV");
  t->DrawLatex(0.2,0.65,njets);
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");


  //-------------------------------
  // cross section limit
  //-------------------------------

  can->cd(2);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  //hexcl->GetXaxis()->SetRangeUser(0,950);
  gPad->SetLogz();
  hexcl->GetXaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetLabelSize(0.035);
  hexcl->GetZaxis()->SetLabelSize(0.035);
  hexcl->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  hexcl->GetXaxis()->SetTitle("gluino mass (GeV)");
  hexcl->GetZaxis()->SetTitle("#sigma upper limit");
  hexcl->Draw("colz");
  hexcl->SetMinimum(0.001);
  hexcl->SetMaximum(10);
  
  TGraph* gr_excl;      
  TGraph* gr_excl_down;
  TGraph* gr_excl_up;   
  
  if( TString(sample).Contains("gmsb") ) {
    gr_excl      = getGraph(do3jets,"nom");
    gr_excl_down = getGraph(do3jets,"down");
    gr_excl_up   = getGraph(do3jets,"nom");
  }else{
    gr_excl      = getRefXsecGraph(hexcl, "T5zz", 1.0);
    gr_excl_down = getRefXsecGraph(hexcl, "T5zz", 1./3.);
    gr_excl_up   = getRefXsecGraph(hexcl, "T5zz", 3.);
  }

  gr_excl->SetLineWidth(2);
  gr_excl_up->SetLineWidth(2);
  gr_excl_down->SetLineWidth(2);
  gr_excl_up->SetLineStyle(2);
  gr_excl_down->SetLineStyle(3);
  gr_excl->Draw("same");
  //gr_excl_up->Draw("same");
  gr_excl_down->Draw("same");


  TLegend *leg = new TLegend(0.2,0.53,0.45,0.67);
  leg->AddEntry(gr_excl,     "#sigma^{NLO-QCD}","l");
  //leg->AddEntry(gr_excl_up,  "3 #times #sigma^{NLO-QCD}","l");
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
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");


  if( print ){
    can->Print(Form("cards/%s/plots/SMS.eps",version));
    can->Print(Form("cards/%s/plots/SMS.png",version));
    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS.eps cards/%s/plots/SMS.pdf",version,version));
  }




  TCanvas *c2 = new TCanvas("c2","c2",1500,500);
  c2->Divide(3,1);

  c2->cd(1);
  hexcluded13->Draw("colz");
  gr_excl_down->Draw();

  c2->cd(2);
  hexcluded->Draw("colz");
  gr_excl->Draw();

  c2->cd(3);
  hexcluded3->Draw("colz");


}
