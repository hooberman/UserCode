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

void plotYields( int mgbin , int mlbinbin , char* version ){

  TFile *f = TFile::Open(Form("rootfiles/%s/SMS_%i_%i.root",version,mgbin,mlbinbin));

  TH1F* hsig       = (TH1F*) f->Get(Form("histo_SMS_%i_%i",mgbin,mlbinbin));
  TH1F* hsigup     = (TH1F*) f->Get(Form("histo_SMS_%i_%i_JES_shapeUp",mgbin,mlbinbin));
  TH1F* hsigdn     = (TH1F*) f->Get(Form("histo_SMS_%i_%i_JES_shapeDown",mgbin,mlbinbin));
  TH1F* hbkg       = (TH1F*) f->Get("histo_bkg");
  TH1F* hbkgerrup  = (TH1F*) f->Get("histo_bkg_errUp");
  TH1F* hbkgerrdn  = (TH1F*) f->Get("histo_bkg_errDown");
  TH1F* hdata      = (TH1F*) f->Get("histo_Data");

  // TCanvas* can = new TCanvas();

  // hsig->GetXaxis()->SetBinLabel(1,"SR1 SF");
  // hsig->GetXaxis()->SetBinLabel(2,"SR1 OF");
  // hsig->GetXaxis()->SetBinLabel(3,"SR2 SF");
  // hsig->GetXaxis()->SetBinLabel(4,"SR2 OF");
  // hsig->GetXaxis()->SetBinLabel(5,"SR3 SF");
  // hsig->GetXaxis()->SetBinLabel(6,"SR3 OF");

  // hsig->Draw();
  // hbkg->SetLineColor(2);
  // hbkg->Draw("same");


  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "yield      " << ibin << " " << Form("%.1f",hsig->GetBinContent(ibin)) << endl;
  }

  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "yieldup    " << ibin << " " << Form("%.1f",hsigup->GetBinContent(ibin)) << endl;
  }

  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "yielddn    " << ibin << " " << Form("%.1f",hsigdn->GetBinContent(ibin)) << endl;
  }

  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "bkg        " << ibin << " " << Form("%.1f",hbkg->GetBinContent(ibin)) << endl;
  }

  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "bkgerrup  " << ibin << " " << Form("%.1f",hbkgerrup->GetBinContent(ibin)) << endl;
  }

  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "bkgerrdn  " << ibin << " " << Form("%.1f",hbkgerrdn->GetBinContent(ibin)) << endl;
  }

  cout << endl;
  for( unsigned int ibin = 1 ; ibin <= 3 ; ibin++ ){
    cout << "data       " << ibin << " " << Form("%.1f",hdata->GetBinContent(ibin)) << endl;
  }





}
