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

float getObservedLimit( int metcut , float jes );

using namespace std;

void SMS(bool print = false){

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------
  
  const float denom    = 20000;
  const float lumi     = 3500;
  const char* filename = "../../output/V00-02-01/T5zz_baby.root";

  cout << "Using file        " << filename << endl;
  cout << "Using denominator " << denom    << " events" << endl;
  cout << "Using lumi        " << lumi     << " pb-1" << endl;

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("T1");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("reference_xSec.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("gluino");

  //--------------------------------------------------
  // preselection
  //--------------------------------------------------

  TCut zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets >= 2");
  TCut sf("leptype==0||leptype==1");
  TCut met30 ("pfmet>30");
  TCut met60 ("pfmet>60");
  TCut met100("pfmet>100");
  TCut met200("pfmet>200");
  TCut met300("pfmet>300");

  TCut presel  = zmass + njets2 + sf;

  cout << "Using selection   " << presel.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<float>   ul;
  vector<int>     cuts;

  //sigcuts.push_back(TCut(presel+met30));      signames.push_back("E_{T}^{miss} > 30 GeV");     ul.push_back(2518.);   labels.push_back("met030"); cuts.push_back(30);
  //sigcuts.push_back(TCut(presel+met60));      signames.push_back("E_{T}^{miss} > 60 GeV");     ul.push_back(134.);    labels.push_back("met060"); cuts.push_back(60);
  //sigcuts.push_back(TCut(presel+met100));     signames.push_back("E_{T}^{miss} > 100 GeV");    ul.push_back(35.0);    labels.push_back("met100"); cuts.push_back(100);
  sigcuts.push_back(TCut(presel+met200));     signames.push_back("E_{T}^{miss} > 200 GeV");    ul.push_back(7.2);     labels.push_back("met200"); cuts.push_back(200);
  //sigcuts.push_back(TCut(presel+met300));     signames.push_back("E_{T}^{miss} > 300 GeV");    ul.push_back(3.0);     labels.push_back("met300"); cuts.push_back(300);

  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* heff[nsig];
  TH2F* heffup[nsig];
  TH2F* heffdn[nsig];
  TH2F* hxsec[nsig];
  TH2F* hexcl[nsig];
  TH2F* hjes[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    TString jesup(sigcuts.at(i));
    jesup.ReplaceAll("njets","njetsup");
    jesup.ReplaceAll("pfmet","pfmetup");

    TString jesdn(sigcuts.at(i));
    jesdn.ReplaceAll("njets","njetsdn");
    jesdn.ReplaceAll("pfmet","pfmetdn");

    TCut jesupcut(jesup);
    TCut jesdncut(jesdn);

    cout << endl << endl;
    cout << "Signal region : " << labels.at(i)  << endl;
    cout << "Selection     : " << sigcuts.at(i) << endl;
    cout << "Selection up  : " << jesupcut      << endl;
    cout << "Selection dn  : " << jesdncut      << endl;

    heff[i]     = new TH2F(Form("heff_%i",i)    , Form("heff_%i",i)    , 48,0,1200,48,0,1200);
    heffup[i]   = new TH2F(Form("heffup_%i",i)  , Form("heffup_%i",i)  , 48,0,1200,48,0,1200);
    heffdn[i]   = new TH2F(Form("heffdn_%i",i)  , Form("heffdn_%i",i)  , 48,0,1200,48,0,1200);
    hxsec[i]    = new TH2F(Form("hxsec_%i",i)   , Form("hxsec_%i",i)   , 48,0,1200,48,0,1200);
    hexcl[i]    = new TH2F(Form("hexcl_%i",i)   , Form("hexcl_%i",i)   , 48,0,1200,48,0,1200);
    hjes[i]     = new TH2F(Form("hjes_%i",i)    , Form("hjes_%i",i)    , 48,0,1200,48,0,1200);

    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i));
    heff[i]->Scale(1./denom);

    ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut);
    heffup[i]->Scale(1./denom);

    ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut);
    heffdn[i]->Scale(1./denom);

    for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin)-12.5;
	float ml = heff[i]->GetYaxis()->GetBinCenter(jbin)-12.5;

	float eff    = heff[i]->GetBinContent(ibin,jbin);
	float effup  = heffup[i]->GetBinContent(ibin,jbin);
	float effdn  = heffdn[i]->GetBinContent(ibin,jbin);

	if( eff < 1e-20 ) continue;
	if( effdn < 1e-20 ) cout << "Error eff JES down " << effdn << endl;

	// cout << endl;
	// cout << "mg ml     " << mg << " " << ml << endl;
	// cout << "eff       " << eff   << endl;
	// cout << "effup     " << effup << endl;
	// cout << "effdn     " << effdn << endl;

	//if( effdn > 0 ){
	float dup    = effup/eff-1;
	float ddn    = 1-effdn/eff;
	float djes   = 0.5 * (dup+ddn);
	hjes[i]->SetBinContent(ibin,jbin,djes);

	// cout << "dup       " << dup  << endl;
	// cout << "ddn       " << ddn  << endl;
	// cout << "djes      " << djes << endl;
	//}

	float this_ul = getObservedLimit( cuts.at(i) , djes );

	float xsecul = this_ul / ( lumi * eff * 0.19 );
	//float xsecul = ul.at(i) / ( lumi * eff * 0.19 );
	if( eff > 0 ) hxsec[i]->SetBinContent(ibin,jbin, xsecul );
	
	int   bin = refxsec->FindBin(mg);
	float xsec = refxsec->GetBinContent(bin);

	hexcl[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul )   hexcl[i]->SetBinContent(ibin,jbin,1);
	//cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
      }
    }
  }

  delete ctemp;

  //--------------------------------------------------
  // make pretty pictures
  //--------------------------------------------------
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);

  TCanvas* can[nsig];

  for( unsigned int i = 0 ; i < nsig ; ++i ){
  
    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    //can[i]->Divide(2,1);
    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1800,600);
    //can[i]->Divide(3,1);
    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,1200);
    can[i]->Divide(2,2);

    //-------------------------------
    // efficiency
    //-------------------------------

    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    heff[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    heff[i]->GetZaxis()->SetTitle("efficiency");
    heff[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 3.5 fb^{-1}");

    //-------------------------------
    // cross section
    //-------------------------------
  
    can[i]->cd(2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    gPad->SetLogz();
    hxsec[i]->GetXaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hxsec[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hxsec[i]->GetZaxis()->SetTitle("#sigma upper limit");
    hxsec[i]->Draw("colz");
    hxsec[i]->SetMinimum(0.01);
    hxsec[i]->SetMaximum(10);

    TGraph* gr_excl      = getRefXsecGraph(hxsec[i], "T5zz", 1.0);
    TGraph* gr_excl_down = getRefXsecGraph(hxsec[i], "T5zz", 1./3.);
    TGraph* gr_excl_up   = getRefXsecGraph(hxsec[i], "T5zz", 3.);

    gr_excl->SetLineWidth(2);
    gr_excl_up->SetLineWidth(2);
    gr_excl_down->SetLineWidth(2);
    gr_excl_up->SetLineStyle(2);
    gr_excl_down->SetLineStyle(3);
    gr_excl->Draw("same");
    gr_excl_up->Draw("same");
    gr_excl_down->Draw("same");

    TLegend *leg = new TLegend(0.2,0.53,0.53,0.67);
    leg->AddEntry(gr_excl,"#sigma^{prod} = #sigma^{NLO-QCD}","l");
    leg->AddEntry(gr_excl_up,"#sigma^{prod} = 3 #times #sigma^{NLO-QCD}","l");
    leg->AddEntry(gr_excl_down,"#sigma^{prod} = 1/3 #times #sigma^{NLO-QCD}","l");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 3.5 fb^{-1}");

    //-------------------------------
    // excluded points
    //-------------------------------

    can[i]->cd(3);    
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hexcl[i]->GetZaxis()->SetTitle("excluded points");
    hexcl[i]->Draw("colz");
    gr_excl->Draw("same");

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 3.5 fb^{-1}");

    //-------------------------------
    // JES uncertainty
    //-------------------------------

    can[i]->cd(4);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hjes[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hjes[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hjes[i]->GetZaxis()->SetTitle("JES uncertainty");
    hjes[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 3.5 fb^{-1}");




    if( print ){
      can[i]->Print(Form("../plots/%s.eps",labels.at(i).c_str()));
      gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.eps  ../plots/%s.pdf",labels.at(i).c_str(),labels.at(i).c_str()));
    }

    int bin = heff[i]->FindBin(600,200);
    cout << "efficiency (600,200) " << heff[i]->GetBinContent(bin) << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin) << endl;

  }

}



float getObservedLimit( int metcut , float jes ){

  float ul = 999;

  if( jes > 0.5 ){
    cout << "ERROR! JES too high " << jes << ", quitting" << endl;
    exit(0);
  }

  if( metcut == 100 ){
    if     ( jes > 0.00 && jes < 0.10 ) ul = 39.2;
    else if( jes > 0.10 && jes < 0.20 ) ul = 42.5;
    else if( jes > 0.20 && jes < 0.30 ) ul = 45.4;
    else if( jes > 0.30 && jes < 0.40 ) ul = 49.7;
    else if( jes > 0.40 && jes < 0.50 ) ul = 54.3;
  }
  else if( metcut == 200 ){
    if     ( jes > 0.00 && jes < 0.10 ) ul = 6.6;
    else if( jes > 0.10 && jes < 0.20 ) ul = 6.9;
    else if( jes > 0.20 && jes < 0.30 ) ul = 7.9;
    else if( jes > 0.30 && jes < 0.40 ) ul = 8.7;
    else if( jes > 0.40 && jes < 0.50 ) ul = 9.5;
  }
  else if( metcut == 300 ){
    if     ( jes > 0.00 && jes < 0.10 ) ul = 3.1;
    else if( jes > 0.10 && jes < 0.20 ) ul = 3.2;
    else if( jes > 0.20 && jes < 0.30 ) ul = 3.3;
    else if( jes > 0.30 && jes < 0.40 ) ul = 3.5;
    else if( jes > 0.40 && jes < 0.50 ) ul = 3.8;
  }  
  else{
    cout << "ERROR! unrecognized met cut " << metcut << ", quitting" << endl;
    exit(0);
  }

  return ul;


}



float getExpectedLimit( int metcut , float jes ){

  float ul = 999;

  if( jes > 0.5 ){
    cout << "ERROR! JES too high " << jes << ", quitting" << endl;
    exit(0);
  }

  if( metcut == 100 ){
    if     ( jes > 0.00 && jes < 0.10 ) ul = 46.5;
    else if( jes > 0.10 && jes < 0.20 ) ul = 48.3;
    else if( jes > 0.20 && jes < 0.30 ) ul = 52.2;
    else if( jes > 0.30 && jes < 0.40 ) ul = 57.0;
    else if( jes > 0.40 && jes < 0.50 ) ul = 61.9;
  }
  else if( metcut == 200 ){
    if     ( jes > 0.00 && jes < 0.10 ) ul = 9.2;
    else if( jes > 0.10 && jes < 0.20 ) ul = 9.8;
    else if( jes > 0.20 && jes < 0.30 ) ul = 10.4;
    else if( jes > 0.30 && jes < 0.40 ) ul = 11.3;
    else if( jes > 0.40 && jes < 0.50 ) ul = 12.0;
  }
  else if( metcut == 300 ){
    if     ( jes > 0.00 && jes < 0.10 ) ul = 4.3;
    else if( jes > 0.10 && jes < 0.20 ) ul = 4.5;
    else if( jes > 0.20 && jes < 0.30 ) ul = 4.8;
    else if( jes > 0.30 && jes < 0.40 ) ul = 5.2;
    else if( jes > 0.40 && jes < 0.50 ) ul = 5.2;
  }  
  else{
    cout << "ERROR! unrecognized met cut " << metcut << ", quitting" << endl;
    exit(0);
  }

  return ul;


}




