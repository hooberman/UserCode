#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCut.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include <iomanip>

using namespace std;

void plotMVAOutput(){

  //char* path    = "Trainings/H130_WWTo2L2Nu/output";
  //char* path    = "Trainings/H130_WWTo2L2Nu_WJetsToLNu/output";
  char* path    = "Trainings/H130_allbkg/output";

  //char* mvaname = "MVA_PDERS";
  //char* mvaname = "MVA_MLPBNN";
  char* mvaname = "MVA_BDT";
  //char* mvaname = "LikelihoodPCA";

  int rebin     = 10;
  int colors[]  = { 5 , 2 , 4 , 3 , 7 , 8 , 6 , 9 , 10};

  vector<char*> samples;
  samples.push_back("WWTo2L2Nu");
  samples.push_back("GluGluToWWTo4L");
  samples.push_back("WZ");
  samples.push_back("ZZ");
  samples.push_back("TTJets");
  samples.push_back("tW");
  samples.push_back("WJetsToLNu");
  samples.push_back("DY");
  const unsigned int nsamples = samples.size();

  char* higgssample = "Higgs130";
  //char* higgssample = "Higgs160";
  //char* higgssample = "Higgs200";



  TFile*   file     = new TFile();
  TH1F*    hist     = new TH1F();
  TH1F*    bkghist  = new TH1F();
  THStack* bkgstack = new THStack("bkgstack","bkgstack");
  TLegend *leg      = new TLegend(0.3,0.7,0.5,0.9);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);

  //loop over backgrounds
  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    file = TFile::Open(Form("%s/%s.root",path,samples.at(i)));
    hist = (TH1F*) file->Get( mvaname );
    hist->SetFillColor(colors[i]);
    
    leg->AddEntry(hist,samples.at(i),"f");
  
    if( i == 0 ) bkghist = (TH1F*) hist->Clone();
    else         bkghist -> Add(hist);

    hist->Rebin( rebin );
    bkgstack->Add(hist);

  }

  //higgs sample
  file = TFile::Open(Form("%s/%s.root",path,higgssample));
  TH1F* higgshist = (TH1F*) file->Get( mvaname );
  higgshist->SetLineWidth(2);
  leg->AddEntry(higgshist,higgssample,"l");



 

  float bkg = 0;
  float sig = 0;
  float minbkg = 1.48;

  for( int ibin = 1 ; ibin < bkghist->GetNbinsX() ; ibin++ ){

    bkg = bkghist->Integral( ibin , 10000 );
    sig = higgshist->Integral( ibin , 10000 );

    if( bkg < minbkg ){
      cout << "cut value " << ibin << " " << bkghist->GetBinCenter(ibin) << endl;
      cout << "sig " << sig << " bkg " << bkg << endl;   
      break;
    }

  }



  bkghist->Rebin( rebin );
  higgshist->Rebin( rebin );

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  gPad->SetLogy();
  bkghist->GetXaxis()->SetTitle(Form("%s output",mvaname));
  bkghist->Draw();
  bkgstack->Draw("same");
  higgshist->Draw("same");
  bkghist->Draw("axissame");
  leg->Draw();

}
