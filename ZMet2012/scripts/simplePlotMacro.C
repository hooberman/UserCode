#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

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
#include "TColor.h"
#include "TMath.h"

#include "mycolors.h"

//#include "histtools.h"
using namespace std;

bool  doKscaling =   true;
float K          =   0.14;
float Rem        =   1.22;   // OR of all triggers
int   rebin      =     10;
bool  bveto      =   false;

float histError( TH1F* hist , int lowbin ){

  float err2 = 0;

  for( int ibin = lowbin ; ibin <= hist->GetXaxis()->GetNbins() ; ibin++ ){
    err2 += pow( hist->GetBinError(ibin) , 2 );
  }

  return sqrt(err2);
}


void simplePlotMacro( bool printplots = false ){

  //-----------------------------------
  // data/MC files
  //-----------------------------------
  
  char* iter = "V00-00-17";

  TFile *f   = new TFile();
  TFile *fwz = new TFile();
  TFile *fzz = new TFile();

  TChain *chwz = new TChain("T1");
  TChain *chzz = new TChain("T1");

  chwz->Add(Form("../output/%s/wz_baby.root",iter));
  chzz->Add(Form("../output/%s/zz_baby.root",iter));

  if( bveto ){
    f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));
    fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));
    fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));

    K = 0.13;
  }

  else{
    // f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplatenjetsgeq2.root",iter));
    // fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplatenjetsgeq2.root",iter));
    // fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplatenjetsgeq2.root",iter));

    f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplate_pfmet.root",iter));
    fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplate_pfmet.root",iter));
    fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplate_pfmet.root",iter));

    // f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplate_t1pfmet.root",iter));
    // fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplate_t1pfmet.root",iter));
    // fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplate_t1pfmet.root",iter));
  }

  cout << "B-veto?   " << bveto << endl;
  cout << "K         " << K     << endl;

  //-----------------------------------
  // OF prediction
  //-----------------------------------

  TH1F* h_of = new TH1F();
  
  if( doKscaling ){
    h_of = (TH1F*) f->Get("metObserved_df_nozveto");
    h_of->Sumw2();
    h_of->Scale(K);
  }

  else{
    h_of    = (TH1F*) f->Get("metObserved_df");
    h_of->Sumw2();
  }

  //-----------------------------------
  // define plots to make
  //-----------------------------------

  vector<char*> observedHisto;
  vector<char*> predictedHisto;

  observedHisto.push_back((char*)"metObserved");        predictedHisto.push_back((char*)"metPredicted");
  observedHisto.push_back((char*)"metObserved_ee");     predictedHisto.push_back((char*)"metPredicted_ee");
  observedHisto.push_back((char*)"metObserved_mm");     predictedHisto.push_back((char*)"metPredicted_mm");

  //-----------------------------------
  // define histos, canvas, pads
  //-----------------------------------

  const unsigned int nplots = observedHisto.size();

  TCanvas* can[nplots];
  TPad*    mainpad[nplots];
  TPad*    respad[nplots];
  THStack* pred[nplots];
  TH1F*    h_sf[nplots];
  TH1F*    h_gjets[nplots];
  TH1F*    h_ofpred[nplots];
  TH1F*    hratio[nplots];
  TH1F*    htotpred[nplots];
  TH1F*    h_wz[nplots];
  TH1F*    h_zz[nplots];
  TH1F*    h_vz[nplots];
  TH1F*    hsysterr[nplots];

  TLatex *t = new TLatex();
  t->SetNDC();

  for( unsigned int i = 0 ; i < nplots ; ++i ){

    //------------------------------------------
    // get data observed and predicted histos
    //------------------------------------------
    
    h_sf[i]     = (TH1F*) f->Get(observedHisto.at(i));           
    h_gjets[i]  = (TH1F*) f->Get(predictedHisto.at(i));          
    h_ofpred[i] = (TH1F*) h_of->Clone(Form("hofpred_%i",i));     

    h_wz[i]     = (TH1F*) fwz->Get(observedHisto.at(i));
    h_zz[i]     = (TH1F*) fzz->Get(observedHisto.at(i));

    h_vz[i]     = (TH1F*) h_wz[i]->Clone(Form("h_vz_%i",i));
    h_vz[i]->Add(h_zz[i]);
    
    //------------------------------------------
    // scale OF prediction for ee vs. mm
    //------------------------------------------

    char* title     = (char*) "ee/#mu#mu events";
    bool  ee_and_mm = true;

    if( TString(observedHisto.at(i)).Contains("ee") ){
      h_ofpred[i]->Scale(0.5/Rem);
      title     = (char*) "ee events";
      ee_and_mm = false;
    }

    else if( TString(observedHisto.at(i)).Contains("mm") ){
      h_ofpred[i]->Scale(0.5*Rem);
      title     = (char*) "#mu#mu events";
      ee_and_mm = false;
    }

    //------------------------------------------
    // rebin, style, sum of predictions
    //------------------------------------------

    h_ofpred[i]->Rebin(rebin);
    h_sf[i]->Rebin(rebin);
    h_gjets[i]->Rebin(rebin);
    h_wz[i]->Rebin(rebin);
    h_zz[i]->Rebin(rebin);
    h_vz[i]->Rebin(rebin);

    // h_gjets[i]->SetFillColor(kAzure-9);
    // h_ofpred[i]->SetFillColor(kGreen+2);
    // h_zz[i]->SetFillColor(kRed);
    // h_wz[i]->SetFillColor(kYellow);

    // h_gjets[i]->SetFillColor(kYellow-9);
    // h_ofpred[i]->SetFillColor(kAzure-9);
    // h_zz[i]->SetFillColor(kYellow-9);
    // h_vz[i]->SetFillColor(kOrange);
    // h_vz[i]->SetFillColor(kGreen-6);
    // h_zz[i]->SetFillColor(kYellow-9);
    // h_wz[i]->SetFillColor(kGray+1);

    // h_gjets[i]->SetFillColor(myBurntOrange);
    // h_ofpred[i]->SetFillColor(kYellow);
    // h_zz[i]->SetFillColor(kRed);
    // h_wz[i]->SetFillColor(kGreen+2);

    //h_gjets[i]->SetFillColor(kYellow);
    //h_gjets[i]->SetFillColor(kCyan-10);
    //h_wz[i]->SetFillColor(kRed);
    // h_gjets[i]->SetFillStyle(3003);
    // h_ofpred[i]->SetFillStyle(3004);
    // h_wz[i]->SetFillStyle(3005);
    // h_zz[i]->SetFillStyle(3002);

    h_gjets[i]->SetFillColor(50);
    h_ofpred[i]->SetFillColor(42);
    h_vz[i]->SetFillColor(31);

    cout << "SF events " << h_sf[i]->GetEntries() << endl;
    cout << "OF events " << h_ofpred[i]->GetEntries() << endl;

    pred[i] = new THStack();
    // pred[i]->Add(h_wz[i]);
    // pred[i]->Add(h_zz[i]);
    pred[i]->Add(h_vz[i]);
    pred[i]->Add(h_ofpred[i]);
    pred[i]->Add(h_gjets[i]);

    htotpred[i] = (TH1F*) h_ofpred[i]->Clone(Form("htotpred_%i",i));
    htotpred[i]->Add(h_gjets[i]);
    htotpred[i]->Add(h_wz[i]);
    htotpred[i]->Add(h_zz[i]);
    for( int ibin = 1 ; ibin <= htotpred[i]->GetXaxis()->GetNbins() ; ibin++ ){
      htotpred[i]->SetBinError(ibin,0);
    }

    //-----------------------------------------------
    // make tables
    //-----------------------------------------------

    const unsigned int nbins = 6;

    //int bins[nbins]={0,30,60,100,150,200,300};
    int      bins[nbins]  = {0,30,60,100,200,300};
    Double_t xbins[nbins] = {0.0,30.0,60.0,100.0,200.0,300.0};

    int width1 = 16;
    int width2 =  2;

    cout << endl;
    cout << title << endl << endl;
    
    cout << "|" << setw(width1) << "" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << "|" << setw(width1) << Form("MET>%i GeV",bins[ibin]) << setw(width2);
    }
    cout << "|" << endl;
      
    int   ndata[nbins];
    float ngjets[nbins];
    float nof[nbins];
    float nwz[nbins];
    float nzz[nbins];
    float ntot[nbins];

    float ngjets_syst[nbins];
    float nof_syst[nbins];
    float nwz_syst[nbins];
    float nzz_syst[nbins];
    float ntot_syst[nbins];

    float ngjets_stat[nbins];
    float nof_stat[nbins];
    float nwz_stat[nbins];
    float nzz_stat[nbins];
    float ntot_stat[nbins];

    float ngjets_toterr[nbins];
    float nof_toterr[nbins];
    float nwz_toterr[nbins];
    float nzz_toterr[nbins];
    float ntot_toterr[nbins];

    float excess[nbins];

    hsysterr[i] = new TH1F(Form("hsysterr_%i",i),Form("hsysterr_%i",i),nbins-1,xbins);

    for( unsigned int ibin = 0 ; ibin < nbins ; ++ibin ){
      int bin      = h_sf[i]->FindBin(bins[ibin]);

      if( bin > h_sf[i]->GetNbinsX() ) bin = h_sf[i]->GetNbinsX();

      // values
      ndata[ibin]  = h_sf[i]->Integral(bin,1000);
      ngjets[ibin] = h_gjets[i]->Integral(bin,1000);
      nof[ibin]    = h_ofpred[i]->Integral(bin,1000);
      nwz[ibin]    = h_wz[i]->Integral(bin,1000);
      nzz[ibin]    = h_zz[i]->Integral(bin,1000);
      ntot[ibin]   = htotpred[i]->Integral(bin,1000);

      // syst uncertainties
      ngjets_syst[ibin] = 0.3 * h_gjets[i]->Integral(bin,1000);

      //float ofsyst = 0.1;
      //if( doKscaling && bins[ibin] >= 200 ) ofsyst = 0.25;

      float ofsyst = 1;

      if( bveto ){
	ofsyst = 0.23;
	if( !ee_and_mm ) ofsyst = 0.25;
      }
      else{
	ofsyst = 0.15;
	if( !ee_and_mm ) ofsyst = 0.2;
      }

      nof_syst[ibin]    = ofsyst * h_ofpred[i]->Integral(bin,1000);
      nwz_syst[ibin]    = 0.5 * h_wz[i]->Integral(bin,1000);
      nzz_syst[ibin]    = 0.5 * h_zz[i]->Integral(bin,1000);
      ntot_syst[ibin]   = sqrt( pow(ngjets_syst[ibin],2) + pow(nof_syst[ibin],2) + pow(nwz_syst[ibin],2) + pow(nzz_syst[ibin],2));
      
      // stat uncertainties
      ngjets_stat[ibin] = histError( h_gjets[i]  , bin );
      nof_stat[ibin]    = histError( h_ofpred[i] , bin );
      nwz_stat[ibin]    = histError( h_wz[i]     , bin );
      nzz_stat[ibin]    = histError( h_zz[i]     , bin );
      ntot_stat[ibin]   = sqrt( pow(ngjets_stat[ibin],2) + pow(nof_stat[ibin],2) + pow(nwz_stat[ibin],2) + pow(nzz_stat[ibin],2));

      // tot uncertainties
      ngjets_toterr[ibin] = sqrt( pow(ngjets_syst[ibin],2) + pow(ngjets_stat[ibin],2) );
      nof_toterr[ibin]    = sqrt( pow(nof_syst[ibin],2) + pow(nof_stat[ibin],2) );
      nwz_toterr[ibin]    = sqrt( pow(nwz_syst[ibin],2) + pow(nwz_stat[ibin],2) );
      nzz_toterr[ibin]    = sqrt( pow(nzz_syst[ibin],2) + pow(nzz_stat[ibin],2) );
      ntot_toterr[ibin]   = sqrt( pow(ngjets_toterr[ibin],2) + pow(nof_toterr[ibin],2) + pow(nwz_toterr[ibin],2) + pow(nzz_toterr[ibin],2));

      excess[ibin]        = (ndata[ibin]-ntot[ibin])/sqrt( pow(ntot_toterr[ibin],2) + ndata[ibin]);

      if( ibin+1 < nbins ){
	hsysterr[i]->SetBinContent(ibin+1,1);
	hsysterr[i]->SetBinError(ibin+1,ntot_toterr[ibin]/ntot[ibin]);
      }
    }


    //-----------------------------
    // g+jets
    //-----------------------------

    cout << "|" << setw(width1) << "Z bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( ngjets[ibin] > 100 ) cout << "|" << setw(width1) << Form("%.0f +/- %.0f",ngjets[ibin],ngjets_toterr[ibin]) << setw(width2);
      else                     cout << "|" << setw(width1) << Form("%.1f +/- %.1f",ngjets[ibin],ngjets_toterr[ibin]) << setw(width2);
    }
    cout << "|" << endl;

    //-----------------------------
    // OF
    //-----------------------------

    cout << "|" << setw(width1) << "OF bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( nof[ibin] > 100 )    cout << "|" << setw(width1) << Form("%.0f +/- %.0f",nof[ibin],nof_toterr[ibin]) << setw(width2);
      else                     cout << "|" << setw(width1) << Form("%.1f +/- %.1f",nof[ibin],nof_toterr[ibin]) << setw(width2);
    }
    cout << "|" << endl;
    
    //-----------------------------
    // WZ
    //-----------------------------

    cout << "|" << setw(width1) << "WZ bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << "|" << setw(width1) << Form("%.1f +/- %.1f",nwz[ibin],nwz_toterr[ibin]) << setw(width2);
    }
    cout << "|" << endl;
    
    //-----------------------------
    // ZZ
    //-----------------------------

    cout << "|" << setw(width1) << "ZZ bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << "|" << setw(width1) << Form("%.1f +/- %.1f",nzz[ibin],nzz_toterr[ibin]) << setw(width2);
      //cout << "|" << setw(width1) << Form("%.1f",nzz[ibin]) << setw(width2);
    }
    cout << "|" << endl;

    //-----------------------------
    // total bkg
    //-----------------------------

    cout << "|" << setw(width1) << "total bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( ntot[ibin] > 100 )   cout << "|" << setw(width1) << Form("%.0f +/- %.0f",ntot[ibin],ntot_toterr[ibin]) << setw(width2);
      else                     cout << "|" << setw(width1) << Form("%.1f +/- %.1f",ntot[ibin],ntot_toterr[ibin]) << setw(width2);

      //      cout << "|" << setw(width1) << Form("%.1f",ntot[ibin]) << setw(width2);
    }
    cout << "|" << endl;

    //-----------------------------
    // data
    //-----------------------------

    cout << "|" << setw(width1) << "data" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << "|" << setw(width1) << ndata[ibin] << setw(width2);
    }
    cout << "|" << endl;
    
    //-----------------------------
    // significance
    //-----------------------------

    cout << "|" << setw(width1) << "significance" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << "|" << setw(width1) << Form("%.1f",excess[ibin]) << setw(width2);
    }
    cout << "|" << endl;

    //------------------------------------------
    // draw plots
    //------------------------------------------
  
    can[i] = new TCanvas(Form("can_%i",i),"",800,800);
    can[i]->cd();

    mainpad[i] = new TPad(Form("mainpad_%i",i),Form("mainpad_%i",i),0.0,0.0,1.0,0.8);
    mainpad[i]->Draw();
    mainpad[i]->cd();
    mainpad[i]->SetLeftMargin(0.15);
    mainpad[i]->SetRightMargin(0.05);    
    mainpad[i]->SetLogy();

    h_sf[i]->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    h_sf[i]->GetYaxis()->SetTitle("entries / 10 GeV");
    h_sf[i]->GetYaxis()->SetTitleOffset(1.0);

    h_sf[i]->SetMinimum(0.1);
    h_sf[i]->Draw("E1");
    pred[i]->Draw("histsame");
    h_sf[i]->Draw("sameE1");
    h_sf[i]->Draw("sameaxis");

    TLegend* leg = new TLegend(0.75,0.6,0.9,0.9);
    leg->AddEntry(h_sf[i],"data","lp");
    leg->AddEntry(h_gjets[i],"Z+jets","f");
    leg->AddEntry(h_ofpred[i],"OF","f");
    //leg->AddEntry(h_zz[i],"ZZ","f");
    //leg->AddEntry(h_wz[i],"WZ","f");
    leg->AddEntry(h_vz[i],"VZ","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    t->SetTextSize(0.04);
    t->DrawLatex(0.4,0.85,"CMS Preliminary");
    t->DrawLatex(0.4,0.79,"#sqrt{s} = 8 TeV, L_{int} = 5.1 fb^{-1}");
    t->DrawLatex(0.4,0.73,title);

    can[i]->cd();

    respad[i] = new TPad(Form("respad_%i",i),Form("respad_%i",i),0.0,0.8,1.0,1.0);
    respad[i]->Draw();
    respad[i]->cd();
    respad[i]->SetTopMargin(0.1);
    respad[i]->SetGridy();
    respad[i]->SetLeftMargin(0.15);
    respad[i]->SetRightMargin(0.05);
    
    gStyle->SetErrorX(0.5);
    //hsysterr[i]->SetFillColor(5);
    hsysterr[i]->SetFillColor(kBlue+1);
    hsysterr[i]->SetFillStyle(3004);
    hsysterr[i]->SetMarkerSize(0);
   
    hratio[i]   = (TH1F*) h_sf[i]->Clone(Form("hratio_%i",i));
    hratio[i]->Divide(htotpred[i]);
    
    hratio[i]->GetXaxis()->SetLabelSize(0.0);
    hratio[i]->GetYaxis()->SetLabelSize(0.2);
    hratio[i]->GetYaxis()->SetTitleSize(0.25);
    hratio[i]->GetYaxis()->SetTitleOffset(0.25);
    hratio[i]->GetYaxis()->SetNdivisions(3);
    hratio[i]->GetYaxis()->SetTitle("ratio");
    hratio[i]->GetXaxis()->SetTitle("");
    
    hratio[i]->SetMinimum(0.5);
    hratio[i]->SetMaximum(1.5);
    
    hratio[i]->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio[i]->Draw();
    hsysterr[i]->Draw("sameE2");
    hratio[i]->Draw("same");
    hratio[i]->Draw("axissame");
	
    TLine line;
    line.DrawLine(0.0,1.0,300,1.0);

    if( printplots ){
      can[i]->Print(Form("../plots/met_%i.pdf",i));
      can[i]->Print(Form("../plots/met_%i.root",i));
      can[i]->Print(Form("../plots/met_%i.C",i));
    }



    

  }
  
}
