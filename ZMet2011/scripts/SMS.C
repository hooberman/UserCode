
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

void SMS(bool print = false){

  float denom = 20000;
  float lumi  = 3500;
  cout << "Using denominator " << denom << endl;
  cout << "Using lumi        " << lumi  << " pb-1" << endl;

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("T1");
  ch->Add("../output/V00-02-00/T5zz_baby.root");

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

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<float>   ul;

  sigcuts.push_back(TCut(presel+met30));      signames.push_back("E_{T}^{miss} > 30 GeV");     ul.push_back(2518.);   labels.push_back("met30");
  sigcuts.push_back(TCut(presel+met60));      signames.push_back("E_{T}^{miss} > 60 GeV");     ul.push_back(134.);    labels.push_back("met60");
  sigcuts.push_back(TCut(presel+met100));     signames.push_back("E_{T}^{miss} > 100 GeV");    ul.push_back(35.0);    labels.push_back("met100");
  sigcuts.push_back(TCut(presel+met200));     signames.push_back("E_{T}^{miss} > 200 GeV");    ul.push_back(7.2);     labels.push_back("met200");
  sigcuts.push_back(TCut(presel+met300));     signames.push_back("E_{T}^{miss} > 300 GeV");    ul.push_back(3.0);     labels.push_back("met300");

  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* heff[nsig];
  TH2F* hxsec[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    cout << endl << endl;
    cout << "Signal region" << endl;
    cout << signames.at(i)  << endl;
    cout << sigcuts.at(i)   << endl;

    heff[i]   = new TH2F(Form("heff_%i",i)  , Form("heff_%i",i)  , 48,0,1200,48,0,1200);
    hxsec[i]  = new TH2F(Form("hxsec_%i",i) , Form("hxsec_%i",i) , 48,0,1200,48,0,1200);

    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i));
    heff[i]->Scale(1./denom);

    for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

	float eff  = heff[i]->GetBinContent(ibin,jbin);
	float xsec = ul.at(i) / ( lumi * eff * 0.19 );
	if( eff > 0 ) hxsec[i]->SetBinContent(ibin,jbin, xsec );
	
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
  
    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    can[i]->Divide(2,1);

    can[i]->cd(1);
    gPad->SetRightMargin(0.2);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    heff[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    heff[i]->Draw("colz");
    t->DrawLatex(0.2,0.70,signames.at(i).c_str());
  
    can[i]->cd(2);
    gPad->SetRightMargin(0.2);
    gPad->SetLogz();
    hxsec[i]->GetXaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hxsec[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hxsec[i]->Draw("colz");
    t->DrawLatex(0.2,0.70,signames.at(i).c_str());

    if( print ){
      can[i]->Print(Form("../plots/%s.eps",labels.at(i).c_str()));
      gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.eps  ../plots/%s.pdf",labels.at(i).c_str(),labels.at(i).c_str()));
    }

    int bin = heff[i]->FindBin(600,200);
    cout << "efficiency (600,200) " << heff[i]->GetBinContent(bin) << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin) << endl;

  }

}
