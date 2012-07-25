//#include "Utils/SMS_utils.C"
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

TH2F* shiftHist(TH2F* hin){

  TH2F* hout = new TH2F(Form("%s_shift",hin->GetName()),Form("%s_shift",hin->GetTitle()), 48,0-12.5,1200-12.5,48,0-12.5,1200-12.5);

  for(int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for(int jbin = 1 ; jbin <= 48 ; jbin++ ){
      hout->SetBinContent(ibin,jbin,hin->GetBinContent(ibin,jbin));
    }
  }

  return hout;
}

TGraph* getGraph_T1lh(string type){

  float x[15];
  float y[15];
  int npoints = -1;

  if( type == "nom" ){
    x[0]  =  800;    y[0]  = 37.5;
    x[1]  =  900;    y[1]  = 187.5;
    x[2]  =  900;    y[2]  = 250.0;
    x[3]  =  875;    y[3]  = 350;
    x[4]  =  837.5;  y[4]  = 425;
    x[5]  =  787.5;  y[5]  = 437.5;
    x[6]  =  700;    y[6]  = 437.5;
    x[7]  =  687.5;  y[7]  = 500;
    x[8]  =  650;    y[8]  = 512.5;
    x[9]  =  575;    y[9]  = 512.5;
    x[10] =  187.5;  y[10] = 162.5;
    x[11] =  187.4;  y[11] =  37.5;
    npoints = 12;
  }
  else if( type == "down" ){
    x[0]  =  775;    y[0]  = 37.5;
    x[1]  =  850;    y[1]  = 187.5;
    x[2]  =  850;    y[2]  = 337.5;
    x[3]  =  800;    y[3]  = 400;
    x[4]  =  775;    y[4]  = 412.5;
    x[5]  =  700;    y[5]  = 412.5;
    x[6]  =  650;    y[6]  = 437.5;
    x[7]  =  637.5;  y[7]  = 500;
    x[8]  =  575;    y[8]  = 500;
    x[9]  =  500;    y[9]  = 445;
    x[10] =  187.5;  y[10] = 162.5;
    x[11] =  187.4;  y[11] =  37.5;
    npoints = 11;
  }
  else if( type == "up" ){
    x[0]  =  850;    y[0]  = 37.5;
    x[1]  =  925;    y[1]  = 212.5;
    x[2]  =  925;    y[2]  = 325.0;
    x[3]  =  912.5;  y[3]  = 412.5;
    x[4]  =  875.0;  y[4]  = 450.0;
    x[5]  =  787.5;  y[5]  = 462.5;
    x[6]  =  712.5;  y[6]  = 462.5;
    x[7]  =  687.5;  y[7]  = 537.5;
    x[8]  =  650;    y[8]  = 537.5;
    x[9]  =  575;    y[9]  = 512.5;
    x[10] =  187.5;  y[10] = 162.5;
    x[11] =  187.4;  y[11] =  37.5;
    npoints = 12;
  }


  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}

TGraph* getGraph_T1lh_expected(string type){

  float x[15];
  float y[15];
  int npoints = -1;

  if( type == "nom" ){
    x[0]  =  800;    y[0]  = 37.5;
    x[1]  =  887.5;  y[1]  = 187.5;
    x[2]  =  887.5;  y[2]  = 387.5;
    x[3]  =  837.5;  y[3]  = 425;
    x[4]  =  787.5;  y[4]  = 437.5;
    x[5]  =  687.5;  y[5]  = 437.5;
    x[6]  =  675;    y[6]  = 487.5;
    x[7]  =  650;    y[7]  = 512.5;
    x[8]  =  575;    y[8]  = 512.5;
    x[9]  =  187.5;  y[9]  = 162.5;
    x[10] =  187.4;  y[10] =  37.5;
    npoints = 11;
  }
  else if( type == "down" ){
    x[0]  =  762.5;  y[0]  = 37.5;
    x[1]  =  837.5;  y[1]  = 187.5;
    x[2]  =  837.5;  y[2]  = 337.5;
    x[3]  =  800;    y[3]  = 400;
    x[4]  =  775;    y[4]  = 412.5;
    x[5]  =  700;    y[5]  = 387.5;
    x[6]  =  637.5;  y[6]  = 387.5;
    x[7]  =  612.5;  y[7]  = 487.5;
    x[8]  =  575;    y[8]  = 487.5;
    x[9]  =  500;    y[9]  = 445;
    x[10] =  187.5;  y[10] = 162.5;
    x[11] =  187.4;  y[11] =  37.5;
    npoints = 11;
  }
  else if( type == "up" ){
    x[0]  =  862.5;  y[0]  = 37.5;
    x[1]  =  937.5;  y[1]  = 212.5;
    x[2]  =  937.5;  y[2]  = 325.0;
    x[3]  =  912.5;  y[3]  = 412.5;
    x[4]  =  875.0;  y[4]  = 512.5;
    x[5]  =  825;    y[5]  = 512.5;
    x[6]  =  712.5;  y[6]  = 475.0;
    x[7]  =  687.5;  y[7]  = 562.5;
    x[8]  =  650;    y[8]  = 562.5;
    x[9]  =  575;    y[9]  = 512.5;
    x[10] =  187.5;  y[10] = 162.5;
    x[11] =  187.4;  y[11] =  37.5;
    npoints = 12;
  }


  TGraph *gr = new TGraph(npoints,x,y);
  return gr;

}

void smoothHist( TH2F* h ){

  vector<int> binx;
  vector<int> biny;
  binx.push_back(8);    biny.push_back(4);
  // binx.push_back(33);    biny.push_back(7);
  // binx.push_back(35);    biny.push_back(7);

  // binx.push_back(28);    biny.push_back(26);
  // binx.push_back(34);    biny.push_back(30);
  // binx.push_back(34);    biny.push_back(31);

  // binx.push_back(42);    biny.push_back(37);
  // binx.push_back(45);    biny.push_back(37);
  // binx.push_back(45);    biny.push_back(36);

  const unsigned int npoints = binx.size();

  for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

      float val = h->GetBinContent(ibin,jbin);

      //if( val < 1e-10 ) continue;
      //cout << ibin << " " << jbin << " " << val << endl;

      for(unsigned int ipoint = 0 ; ipoint < npoints ; ipoint++ ){
	if( ibin == binx.at(ipoint) && jbin == biny.at(ipoint) ){

	  float valup  = h->GetBinContent(ibin+1,jbin+1);
	  float valdn  = h->GetBinContent(ibin-1,jbin);
	  float valavg = 0.5 * (valup+valdn);

	  h->SetBinContent(ibin,jbin,valavg);
	  //h->SetBinContent(ibin,jbin,100000);
	}
	
      }
    }
  }  
}

TH2F* fixupHist( TH2F* hin ){

  TH2F* hclone = (TH2F*) hin->Clone("hfixup");

  for( int i = 1 ; i <= 46 ; ++i ){
    int ibin = i + 2;
    int jbin = i;

    float val = hclone->GetBinContent(ibin,jbin);
    if( val < 1e-10){

      float npoints = 0;
      float tot     = 0;

      if( hclone->GetBinContent(ibin+1,jbin+1) > 1e-10 ){
	tot += hclone->GetBinContent(ibin+1,jbin+1);
	npoints += 1.0;
      }

      if( hclone->GetBinContent(ibin-1,jbin-1) > 1e-10 ){
	tot += hclone->GetBinContent(ibin-1,jbin-1);
	npoints += 1.0;
      }

      if( npoints > 1e-10) hin->SetBinContent(ibin,jbin,tot/npoints);
      //if( npoints > 1e-10 ) hin->SetBinContent(ibin,jbin,100000);

    }
  }

  return hin;
}

TH2F* smoothAllHist( TH2F* hin ){

  TH2F* hclone = (TH2F*) hin->Clone("hclone_smooth");

  for( unsigned int ibin = 6 ; ibin <= 48 ; ++ibin ){
    for( unsigned int jbin = 4 ; jbin <= 48 ; ++jbin ){
      
      float npoints = 0;
      float tot     = 0;
      
      if( hclone->GetBinContent(ibin,jbin) > 1e-10 ){
	tot += hclone->GetBinContent(ibin,jbin);
	npoints += 1.0;
      }
      if( hclone->GetBinContent(ibin+1,jbin+1) > 1e-10 ){
	tot += hclone->GetBinContent(ibin+1,jbin+1);
	npoints += 1.0;
      }
      if( hclone->GetBinContent(ibin-1,jbin-1) > 1e-10 ){
	tot += hclone->GetBinContent(ibin-1,jbin-1);
	npoints += 1.0;
      }
      if( hclone->GetBinContent(ibin+1,jbin-1) > 1e-10 ){
	tot += hclone->GetBinContent(ibin+1,jbin-1);
	npoints += 1.0;
      }
      if( hclone->GetBinContent(ibin-1,jbin+1) > 1e-10 ){
	tot += hclone->GetBinContent(ibin-1,jbin+1);
	npoints += 1.0;
      }
      
      if( npoints > 1e-10 ) hin->SetBinContent( ibin , jbin , tot/npoints );

    }
  }

  return hin;

}

void combineSMSPlots(bool print = false){
  
  char* version        = "V00-00-14";
  char* sample         = "T1lh";
  char* title          = "m(#tilde{q}) >> m(#tilde{g})";
  float dm             = 0.0;
  bool  fixup          = true;
  //bool  smooth         = true;

  float ymin = 50.;
  if( TString(sample).Contains("gmsb") ) ymin = 100.;

  TFile *file = TFile::Open(Form("cards/%s/observed_limit.root",version));
  //TH2F* hexcl = (TH2F*) file->Get("hexcl");

  TH2F* hexcl_temp   = (TH2F*) file->Get("hexcl");
  TH2F* hexp_temp    = (TH2F*) file->Get("hexp");
  TH2F* hexpp1_temp  = (TH2F*) file->Get("hexpp1");
  TH2F* hexpm1_temp  = (TH2F*) file->Get("hexpm1");

  //TH2F* hexcl_temp = (TH2F*) file->Get("hexp");
  TH2F* hexcl   = shiftHist(hexcl_temp);
  TH2F* hexp    = shiftHist(hexp_temp);
  TH2F* hexpp1  = shiftHist(hexpp1_temp);
  TH2F* hexpm1  = shiftHist(hexpm1_temp);

  if( fixup  ) fixupHist    ( hexcl );
  //if( smooth ) smoothAllHist( hexcl );
  smoothHist( hexcl );

  // TH2F* hexcl_temp  = (TH2F*) file->Get("hexcl");
  // TH2F* hexcl_temp2 = shiftHist(hexcl_temp);
  // TH2F* hexcl       = fixupHist(hexcl_temp2);

  //TH2F* hexcl = (TH2F*) file->Get("hexp");
  //cout << "USING EXPECTED LIMIT!!!!!!!!!!!!" << endl;

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // calculate efficiency
  //-------------------------------

  char* babyversion = "V00-02-24/highpt";
  //if( TString(sample).Contains("hadoop") ) babyversion = "V00-02-05";
  //if( TString(sample).Contains("T5zzh") )  babyversion = "V00-02-05";

  TChain *ch = new TChain("t");
  //ch->Add(Form("output/V00-02-04/%s_baby.root",sample));
  ch->Add(Form("output/%s/%s_smallTree.root",babyversion,sample));

  TCut presel("pass && !passz");
  TCut sig("(pfmet>275 && ht>300) || (pfmet>200 && ht>600)");

  TCut sel  = presel + sig;
  TCut weight("trgeff * lepscale * ndavtxweight * (1.0/1.07)");

  cout << "Using selection: " << sel.GetTitle() << endl;

  TH2F* heff = new TH2F("heff","heff", 48,0-12.5,1200-12.5,48,0-12.5,1200-12.5);
  //TH2F* heff = new TH2F("heff","heff", 48,0,1200,48,0,1200);
  TCanvas *ctemp = new TCanvas();
  ch->Draw("ml:mg>>heff",sel*weight);
  //ch->Draw("ml:mg>>heff",sel,"trgeff * lepscale & ndavtxweight");

  delete ctemp;
  heff->Scale(1./10000.);

  int bin = heff->FindBin(600,200);
  cout << "Efficiency for 600,200 " << heff->GetBinContent(bin) << endl;

  //-------------------------------
  // find excluded points
  //-------------------------------

  TFile *xsecfile = TFile::Open("referenceXSecs.root");
  //TFile *xsecfile = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* refxsec       = (TH1F*) xsecfile->Get("gluino_NLONLL");
  TH1F* refxsec_unc   = (TH1F*) xsecfile->Get("gluino_NLONLL_unc");

  TH2F* hexcluded   = new TH2F("hexcluded","hexcluded", 48,0,1200,48,0,1200);
  TH2F* hexcluded13 = new TH2F("hexcluded13","hexcluded13", 48,0,1200,48,0,1200);
  TH2F* hexcluded3  = new TH2F("hexcluded3","hexcluded3", 48,0,1200,48,0,1200);

  TH2F* hexcl_exp     = new TH2F("hexcl_exp"  ,"hexcl_exp"  , 48,0,1200,48,0,1200);
  TH2F* hexcl_expp1   = new TH2F("hexcl_expp1","hexcl_expp1", 48,0,1200,48,0,1200);
  TH2F* hexcl_expm1   = new TH2F("hexcl_expm1","hexcl_expm1", 48,0,1200,48,0,1200);
  
  for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

      float xsecul     = hexcl->GetBinContent(ibin,jbin);
      float xsec_exp   = hexp->GetBinContent(ibin,jbin);
      float xsec_expp1 = hexpp1->GetBinContent(ibin,jbin);
      float xsec_expm1 = hexpm1->GetBinContent(ibin,jbin);

      if( xsecul < 1.e-10 ) continue;

      float mg = hexcluded->GetXaxis()->GetBinCenter(ibin)-12.5;
      float ml = hexcluded->GetYaxis()->GetBinCenter(jbin)-12.5;

      int   bin      = refxsec->FindBin(mg);
      float xsec     = refxsec->GetBinContent(bin);
      float xsec_unc = refxsec_unc->GetBinContent(bin);

      hexcluded->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul )   hexcluded->SetBinContent(ibin,jbin,1);

      hexcluded3->SetBinContent(ibin,jbin,0);
      //if( 3 * xsec > xsecul )   hexcluded3->SetBinContent(ibin,jbin,1);
      if( (xsec+xsec_unc) > xsecul )   hexcluded3->SetBinContent(ibin,jbin,1);

      hexcluded13->SetBinContent(ibin,jbin,0);
      //if( (1./3.) * xsec > xsecul )   hexcluded13->SetBinContent(ibin,jbin,1);
      if( (xsec-xsec_unc) > xsecul )   hexcluded13->SetBinContent(ibin,jbin,1);

      hexcl_exp->SetBinContent(ibin,jbin,0);
      if( xsec > xsec_exp )  hexcl_exp->SetBinContent(ibin,jbin,1);

      hexcl_expp1->SetBinContent(ibin,jbin,0);
      if( xsec > xsec_expp1 )  hexcl_expp1->SetBinContent(ibin,jbin,1);

      hexcl_expm1->SetBinContent(ibin,jbin,0);
      if( xsec > xsec_expm1 )  hexcl_expm1->SetBinContent(ibin,jbin,1);

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
  if( fixup ) fixupHist(heff);

  heff->GetYaxis()->SetRangeUser(ymin,1200);
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetTitle("#chi^{0}_{1} mass [GeV]");
  heff->GetXaxis()->SetTitle("gluino mass [GeV]");
  heff->GetZaxis()->SetTitle("efficiency #times acceptance");
  heff->Draw("colz");
  heff->GetYaxis()->SetRangeUser(ymin,1200);

  t->SetTextSize(0.04);

  //t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow l^{+}l^{-} #chi_{1}^{0}");
  t->DrawLatex(0.2,0.85,"pp #rightarrow #tilde{g}#tilde{g}");
  t->DrawLatex(0.2,0.79,"#tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow l^{+}l^{-} #chi_{1}^{0}");
  t->DrawLatex(0.2,0.73,"#tilde{g} #rightarrow 2j+#chi_{1}^{0}");
  t->DrawLatex(0.2,0.67,title);


  //t->DrawLatex(0.2,0.77,title);
  //t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 100 GeV");
  //t->DrawLatex(0.2,0.65,njets);
  t->DrawLatex(0.2,0.55,"l^{+}l^{-} + E_{T}^{miss} + H_{T}");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}");

  if(TString(sample).Contains("gmsb") )   line.DrawLine(100-12.5,100-12.5,1200-12.5,1200-12.5);
  else                                    line.DrawLine(50-12.5+dm,50-12.5,1200-12.5,1200-12.5-dm);

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

  TGraph* gr_excl_exp;      
  TGraph* gr_excl_exp_down;
  TGraph* gr_excl_exp_up;   
  
  if( TString(sample).Contains("T1lh") ) {
    gr_excl      = getGraph_T1lh("nom");
    gr_excl_down = getGraph_T1lh("down");
    gr_excl_up   = getGraph_T1lh("up");

    gr_excl_exp      = getGraph_T1lh_expected("nom");
    gr_excl_exp_down = getGraph_T1lh_expected("down");
    gr_excl_exp_up   = getGraph_T1lh_expected("up");
  }
  else{
    // gr_excl      = getRefXsecGraph(hexcl, "T1lh", 1.0);
    // gr_excl_down = getRefXsecGraph(hexcl, "T1lh", 1./3.);
    // gr_excl_up   = getRefXsecGraph(hexcl, "T1lh", 3.);
  }

  gr_excl->SetLineWidth(3);
  gr_excl_up->SetLineWidth(3);
  gr_excl_down->SetLineWidth(3);

  gr_excl->SetMarkerSize(1);
  gr_excl_up->SetMarkerSize(1);
  gr_excl_down->SetMarkerSize(1);

  gr_excl->SetMarkerStyle(20);
  gr_excl_up->SetMarkerStyle(20);
  gr_excl_down->SetMarkerStyle(20);

  gr_excl_exp->SetLineWidth(3);
  gr_excl_exp_up->SetLineWidth(3);
  gr_excl_exp_down->SetLineWidth(3);

  gr_excl_exp->SetMarkerSize(1);
  gr_excl_exp_up->SetMarkerSize(1);
  gr_excl_exp_down->SetMarkerSize(1);

  gr_excl_exp->SetMarkerStyle(20);
  gr_excl_exp_up->SetMarkerStyle(20);
  gr_excl_exp_down->SetMarkerStyle(20);

  gr_excl_down->SetLineStyle(2);
  gr_excl_up->SetLineStyle(2);

  gr_excl_exp->SetLineColor(2);
  gr_excl_exp_down->SetLineColor(2);
  gr_excl_exp_up->SetLineColor(2);

  gr_excl_exp_up->SetLineStyle(2);
  gr_excl_exp_down->SetLineStyle(2);

  gr_excl->Draw("same");
  gr_excl_up->Draw("same");
  gr_excl_down->Draw("same");

  gr_excl_exp->Draw("same");
  gr_excl_exp_up->Draw("same");
  gr_excl_exp_down->Draw("same");

  if(TString(sample).Contains("gmsb") )   line.DrawLine(100-12.5,100-12.5,1200-12.5,1200-12.5);
  else                                    line.DrawLine(50-12.5+dm,50-12.5,1200-12.5,1200-12.5-dm);

  // if(TString(sample).Contains("gmsb") )   line.DrawLine(100,100,1200,1200);
  // else                                    line.DrawLine(50+dm,50,1200,1200-dm);

  TLegend *leg = new TLegend(0.2,0.45,0.35,0.65);
  leg->AddEntry(gr_excl,         "observed","l");
  leg->AddEntry(gr_excl_up,      "observed (#pm1#sigma theory)","l");
  leg->AddEntry(gr_excl_exp,     "expected","l");
  leg->AddEntry(gr_excl_exp_up,  "expected (#pm1#sigma)","l");
  //leg->AddEntry(gr_excl_down,"#sigma^{NLO-QCD}_{DOWN}","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->Draw();
  
  t->SetTextSize(0.025);
  t->DrawLatex(0.2,0.85,"pp #rightarrow #tilde{g}#tilde{g}");
  t->DrawLatex(0.2,0.79,"#tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow l^{+}l^{-} #chi_{1}^{0}");
  t->DrawLatex(0.2,0.73,"#tilde{g} #rightarrow 2j+#chi_{1}^{0}");
  t->DrawLatex(0.2,0.67,title);
  //t->DrawLatex(0.2,0.71,njets);
  t->DrawLatex(0.2,0.38,"l^{+}l^{-} + E_{T}^{miss} + H_{T}");
  t->SetTextSize(0.035);
  //t->DrawLatex(0.18,0.92,"CMS                     #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
  t->DrawLatex(0.18,0.92,"CMS Preliminary       #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}");

  if( print ){
    //can->Print(Form("cards/%s/plots/SMS.eps",version));
    can->Print(Form("cards/%s/plots/SMS.pdf",version));
    //can->Print(Form("cards/%s/plots/SMS.png",version));
    //can->Print(Form("cards/%s/plots/SMS.C",version));

    //gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS.eps cards/%s/plots/SMS_ppt.pdf",version,version));
  }

  TH2F* hexcluded_shifted   = shiftHist( hexcluded   );
  TH2F* hexcluded13_shifted = shiftHist( hexcluded13 );
  TH2F* hexcluded3_shifted  = shiftHist( hexcluded3  );

  TH2F* hexcl_exp_shifted   = shiftHist( hexcl_exp   );
  TH2F* hexcl_expp1_shifted = shiftHist( hexcl_expp1 );
  TH2F* hexcl_expm1_shifted = shiftHist( hexcl_expm1 );

  // TH2F* hexcluded_shifted   = (TH2F*) hexcluded->Clone("hexcluded_shifted");
  // TH2F* hexcluded13_shifted = (TH2F*) hexcluded13->Clone("hexcluded13_shifted");
  // TH2F* hexcluded3_shifted  = (TH2F*) hexcluded3->Clone("hexcluded3_shifted");
  gr_excl->SetName("obs");
  gr_excl->SetTitle("obs");

  gr_excl_up->SetName("obs_theoryup");
  gr_excl_up->SetTitle("obs_theoryup");

  gr_excl_down->SetName("obs_theorydown");
  gr_excl_down->SetTitle("obs_theorydown");

  gr_excl_exp->SetName("expected");
  gr_excl_exp->SetTitle("expected");

  gr_excl_exp_up->SetName("exp_up");
  gr_excl_exp_up->SetTitle("exp_up");

  gr_excl_exp_down->SetName("exp_down");
  gr_excl_exp_down->SetTitle("exp_down");

  TFile* fout = TFile::Open(Form("cards/%s/OS_shape_T1lh_limit.root",version),"RECREATE");
  fout->cd();
  heff->Write();
  hexcl->Write();
  gr_excl->Write();
  gr_excl_up->Write();
  gr_excl_down->Write();
  gr_excl_exp->Write();
  gr_excl_exp_up->Write();
  gr_excl_exp_down->Write();
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
  gr_excl_down->Draw("sameP");
  t->DrawLatex(0.3,0.8,"#sigma^{NLO-QCD}_{DOWN}");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded->Draw("colz");
  hexcluded_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded_shifted->Draw("colz");
  gr_excl->Draw();
  gr_excl->Draw("sameP");
  t->DrawLatex(0.3,0.8,"#sigma^{NLO-QCD}");

  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded3->Draw("colz");
  hexcluded3_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded3_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded3_shifted->Draw("colz");
  gr_excl_up->Draw();
  gr_excl_up->Draw("sameP");
  t->DrawLatex(0.3,0.8,"#sigma^{NLO-QCD}_{UP}");


  TCanvas *c3 = new TCanvas("c3","c3",1500,500);
  c3->Divide(3,1);

  t->SetTextSize(0.07);

  c3->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded13->Draw("colz");
  hexcl_expp1_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcl_expp1_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcl_expp1_shifted->Draw("colz");
  gr_excl_exp_down->Draw();
  gr_excl_exp_down->Draw("sameP");
  t->DrawLatex(0.3,0.8,"expected (-1#sigma)");

  c3->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded->Draw("colz");
  hexcl_exp_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcl_exp_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcl_exp_shifted->Draw("colz");
  gr_excl_exp->Draw();
  gr_excl_exp->Draw("sameP");
  t->DrawLatex(0.3,0.8,"expected");

  c3->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded3->Draw("colz");
  hexcl_expm1_shifted->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcl_expm1_shifted->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcl_expm1_shifted->Draw("colz");
  gr_excl_exp_up->Draw();
  gr_excl_exp_up->Draw("sameP");
  t->DrawLatex(0.3,0.8,"expected (+1#sigma)");

  if( print ){

    c2->Print(Form("cards/%s/plots/SMS_points.pdf",version));
    c3->Print(Form("cards/%s/plots/SMS_points_expected.pdf",version));

    // c2->Print(Form("cards/%s/plots/SMS_points.eps",version));
    // c2->Print(Form("cards/%s/plots/SMS_points.pdf",version));
    // c2->Print(Form("cards/%s/plots/SMS_points.png",version));
    // c2->Print(Form("cards/%s/plots/SMS_points.C",version));

    // gROOT->ProcessLine(Form(".! ps2pdf cards/%s/plots/SMS_points.eps cards/%s/plots/SMS_points_ppt.pdf",version,version));
  }


}
