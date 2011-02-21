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
#include "TLatex.h"
#include "TLine.h"
#include <iomanip>


using namespace std;

TH1F* cloneHist( TH1F* hin ){
  
  TH1F *hout = new TH1F(hin->GetTitle(),hin->GetName(),
                        hin->GetNbinsX(), hin->GetXaxis()->GetXmin() , hin->GetXaxis()->GetXmax() );

  for( unsigned int ibin = 1 ; ibin < hin->GetNbinsX() ; ibin++ ){
    float val = hin->GetBinContent( ibin );
    hout->SetBinContent( ibin , val );
  }

  return hout;

}

void plotMVAOutput( bool printgif = false ){

  //gROOT->ProcessLine(".x selection.h");

  //char* path    = "Trainings/H130_WWTo2L2Nu/output";
  //char* path    = "Trainings/H130_WWTo2L2Nu_WJetsToLNu/output";
  char* path    = "Trainings/H130_allbkg_4vars/output";

  //char* mvaname = "MVA_PDERS";
  //char* mvaname = "MVA_MLPBNN";
  //char* mvaname = "MVA_BDT";
  //char* mvaname = "LikelihoodPCA";

  vector<char*> mvanames;
  mvanames.push_back("BDT");
  mvanames.push_back("MLPBNN");
  const unsigned int nmva = mvanames.size();

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

  TCanvas *can[nmva];


  for( unsigned int imva = 0 ; imva < nmva ; ++imva ){


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
      hist = cloneHist( (TH1F*) file->Get( Form("MVA_%s",mvanames.at(imva) ) ) );
   
      hist->SetFillColor(colors[i]);
    
      leg->AddEntry(hist,samples.at(i),"f");
  
      if( i == 0 ) bkghist = (TH1F*) hist->Clone();
      else         bkghist -> Add(hist);

      hist->Rebin( rebin );
      bkgstack->Add(hist);

    }

    //higgs sample
    file = TFile::Open(Form("%s/%s.root",path,higgssample));
    TH1F* higgshist = cloneHist( (TH1F*) file->Get( Form("MVA_%s",mvanames.at(imva) ) ) );
    higgshist->SetLineWidth(2);
    leg->AddEntry(higgshist,higgssample,"l");



    float bkg    = 0;
    float sig    = 0;
    float minbkg = 1.48;
    //float minbkg = 1.10;
    float cut    = 0.;

    for( int ibin = 1 ; ibin < bkghist->GetNbinsX() ; ibin++ ){

      bkg = bkghist->Integral( ibin , 10000 );
      sig = higgshist->Integral( ibin , 10000 );

      if( bkg < minbkg ){
        cut = bkghist->GetBinCenter(ibin);
        cout << endl;
        cout << "S/B       " << sig/bkg    << endl;
        cout << "Sig       " << sig        << endl;
        cout << "Bkg       " << bkg        << endl;
        cout << "cut value " << cut        << endl;
        break;
      }

    }

    float cutsig = sig;
    float cutbkg = bkg;
 
    float maxfom      = -1;
    float maxfom_sig  = -1;
    float maxfom_bkg  = -1;
    float cutval      = -1;
  

    for( int ibin = 1 ; ibin < bkghist->GetNbinsX() ; ibin++ ){

      bkg = bkghist->Integral( ibin , 10000 );
      sig = higgshist->Integral( ibin , 10000 );

      float fom = sig / sqrt( sig + bkg + pow( 0.35 * bkg , 2 ) );
    
      if( fom > maxfom ){
        maxfom  = fom;
        maxfom_sig = sig;
        maxfom_bkg = bkg;
        cutval  = bkghist->GetBinCenter(ibin);
      }

    }

    cout << endl;
    cout << "Max FOM   " << maxfom        << endl;
    cout << "Sig       " << maxfom_sig    << endl;
    cout << "Bkg       " << maxfom_bkg    << endl;
    cout << "cut value " << cutval        << endl;


    bkghist->Rebin( rebin );
    higgshist->Rebin( rebin );

    can[imva] = new TCanvas(Form("can_%i",imva),Form("can_%i",imva),800,600);
    can[imva]->cd();

    //gPad->SetLogy();
    bkghist->GetXaxis()->SetTitle(Form("%s output",mvanames.at(imva)));
    bkghist->Draw();
    bkgstack->Draw("same");
    higgshist->Scale(10.);
    higgshist->Draw("same");
    bkghist->Draw("axissame");
    //leg->Draw();

    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextColor(2);
    t->DrawLatex(0.2,0.85,Form("FOM: %.2f",maxfom));
    t->SetTextColor(1);
    t->DrawLatex(0.2,0.80,Form("Sig: %.2f",maxfom_sig));
    t->DrawLatex(0.2,0.75,Form("Bkg: %.2f",maxfom_bkg));

    t->SetTextColor(4);
    t->DrawLatex(0.2,0.55,Form("S/B: %.2f",cutsig/cutbkg));
    t->SetTextColor(1);
    t->DrawLatex(0.2,0.50,Form("Sig: %.2f",cutsig));
    t->DrawLatex(0.2,0.45,Form("Bkg: %.2f",cutbkg));
  
    TLine  line;
    line.SetLineColor(2);
    line.DrawLine( cutval , bkghist->GetMinimum() , cutval , 1.05 * bkghist->GetMaximum() );
    line.SetLineColor(4);
    line.DrawLine( cut , bkghist->GetMinimum() , cut , 1.05 * bkghist->GetMaximum() );

    if( printgif ) can[imva]->Print(Form("plots/%s.gif",mvanames.at(imva)));
  }
}
