#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TStyle.h"
#include <sstream>
#include "/tas03/home/benhoob/.root/benstyle.C"

using namespace std;

enum metType   { e_tcmet = 0, e_tcmetNew = 1, e_pfmet = 2};

//----------------------------------------------------------------------------------------
//user params
//----------------------------------------------------------------------------------------
const int      nTrigBins        = 5;
const int      nJetBins         = 3;
const int      nSumJetPtBins    = 7;
const metType  myMetType        = e_pfmet;
const float    maxmet           = 100;
const int      rebin            = 10;
const bool     twoTemplates     = true;
const bool     vertical         = false;

//const char*    filename1        = "../photon_output/V00-00-07/Photon_templates.root";
//const char*    filename2        = "../photon_output/V00-00-07/DoubleElectron_templates.root";

//const char*    filename1        = "../photon_output/V00-00-05/DoubleElectron_templates_1btag.root";
//const char*    filename2        = "../photon_output/V00-00-05/DoubleElectron_templates.root";

//const char*    filename1        = "../photon_output/V00-00-07/DoubleElectron_templates.root";
//const char*    filename2        = "../photon_output/V00-00-07/DoubleElectron_templates_vtxreweight.root";

const char*    filename1        = "../photon_output/V00-00-12/DoubleElectron_templates_vtxreweight.root";
const char*    filename2        = "../photon_output/V00-00-13/DoubleElectron_templates_vtxreweight_pt40.root";

//----------------------------------------------------------------------------------------

void displayPhotonTemplates( bool printgif = false ){
  setBenStyle();

  gStyle->SetOptStat(0);

  TFile *file1 = TFile::Open(filename1);
  TFile *file2;
  if( twoTemplates ) file2 = TFile::Open(filename2);

  TH1F* hmet1[nTrigBins][nJetBins][nSumJetPtBins];
  TH1F* hmet2[nTrigBins][nJetBins][nSumJetPtBins];

  TCanvas *can[nTrigBins];

  Double_t nentries1 = 0;
  Double_t nentries2 = 0;
  TLatex *t = new TLatex();
  t->SetNDC();
  
  for( int iT = 0 ; iT < nTrigBins ; ++iT ){

    if( vertical ){
      can[iT] = new TCanvas(Form("can_%i",iT),Form("can_%i",iT),600,1500);
      can[iT]->Divide(2,5);
    }
    else{
      can[iT] = new TCanvas(Form("can_%i",iT),Form("can_%i",iT),1500,600);
      can[iT]->Divide(5,2);
    }
  
    for( int iJ = 1 ; iJ < nJetBins ; iJ++ ){
      
      for( int iS = 2 ; iS < nSumJetPtBins ; iS++ ){
        
        //can[iJ]->cd(iS-1);
        if( vertical ) can[iT]->cd( 2 * (iS-1) + iJ - 2);
        else           can[iT]->cd( (iJ-1) * 5 + iS-1 );
        
        string mettype = "";
        if( myMetType == e_pfmet ) mettype = "pfmet";
        if( myMetType == e_tcmet ) mettype = "tcmet";
        
        hmet1[iT][iJ][iS] = (TH1F*) file1->Get(Form("%sTemplate_photon_%i_%i_%i",mettype.c_str(),iT,iJ,iS));
        
        if( twoTemplates )
          hmet2[iT][iJ][iS] = (TH1F*) file2->Get(Form("%sTemplate_photon_%i_%i_%i",mettype.c_str(),iT,iJ,iS));
        
        if( hmet1[iT][iJ][iS]->Integral() > 0 )
          gPad->SetLogy(1);
        
        hmet1[iT][iJ][iS]->Rebin(rebin);
        if( maxmet > 0 ) hmet1[iT][iJ][iS]->GetXaxis()->SetRangeUser(0,maxmet);
        hmet1[iT][iJ][iS]->SetLineColor(2);
        hmet1[iT][iJ][iS]->SetMarkerColor(2);
        hmet1[iT][iJ][iS]->SetMarkerSize(0.5);
        hmet1[iT][iJ][iS]->Draw("E1");
        
        hmet1[iT][iJ][iS]->SetMinimum(1e-4);
        hmet1[iT][iJ][iS]->SetMaximum(1.0);
        
        stringstream s1;
        s1 << hmet1[iT][iJ][iS]->GetEntries();
        //s1 << hmet1[iT][iJ][iS]->Integral(  hmet1[iT][iJ][iS]->FindBin(180) , 100000) ;
        nentries1 += hmet1[iT][iJ][iS]->GetEntries();
        t->SetTextColor(2);
        t->SetTextSize(0.1);
        t->DrawLatex(0.6,0.8,s1.str().c_str());
        
        if( twoTemplates ){
          stringstream s2;
          s2 << hmet2[iT][iJ][iS]->GetEntries();
          nentries2 += hmet2[iT][iJ][iS]->GetEntries();
          t->SetTextColor(4);
          t->DrawLatex(0.6,0.7,s2.str().c_str());
          hmet2[iT][iJ][iS]->Rebin(rebin);
          hmet2[iT][iJ][iS]->SetLineColor(1);
          hmet2[iT][iJ][iS]->SetMarkerColor(38);
          hmet2[iT][iJ][iS]->SetFillColor(38);
          hmet2[iT][iJ][iS]->Draw("samehist");
          hmet1[iT][iJ][iS]->Draw("sameE1");
          hmet1[iT][iJ][iS]->Draw("axissame");
        }
      }
    }
    if( printgif) can[iT]->Print(Form("../plots/template_targeted_%i.pdf",iT));
  }

  cout << "nentries1 " << nentries1 << endl;
  cout << "nentries2 " << nentries2 << endl;
}
