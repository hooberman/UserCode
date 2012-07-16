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

bool  doKscaling =    true;
float K          =    0.14;
float Rem        =    1.22;   // OR of all triggers
int   rebin      =      10;
bool  bveto      =   false;
bool  latex      =    true;

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

  //-----------------------------------
  // selection
  //-----------------------------------

  TCut pfleptons("pflep1.pt() > 20 && pflep2.pt() > 20");
  TCut Zmass("dilmasspf>81 && dilmasspf<101");
  TCut njets2("njets>=2");
  TCut ee("leptype==0 && (ee==1 || isdata==0)");
  TCut mm("leptype==1 && (mm==1 || isdata==0)");
  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut sf = ee||mm;
  TCut nu("ngennu>0");
  TCut bvetocut("nbm==0");
  TCut nlep2("nlep==2");
  TCut mjj("mjj>70.0 && mjj<110.0");

  TCut sel;
  sel += njets2;
  sel += pfleptons;
  sel += Zmass;
  sel += (ee||mm);
  sel += nu;
  if( bveto ){
    sel += bvetocut;
    sel += nlep2;
    sel += mjj;
  }

  TCut weight("weight * 5.10 * vtxweight");


  if( bveto ){
    // f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));
    // fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));
    // fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));

    f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplate_pfmet_bveto.root",iter));
    fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplate_pfmet_bveto.root",iter));
    fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplate_pfmet_bveto.root",iter));

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

    //h_wz[i]     = (TH1F*) fwz->Get(observedHisto.at(i));
    //h_zz[i]     = (TH1F*) fzz->Get(observedHisto.at(i));

    //------------------------------------------
    // scale OF prediction for ee vs. mm
    //------------------------------------------

    char* title     = (char*) "ee/#mu#mu events";
    bool  ee_and_mm = true;
    TCut  mysel = sel;

    if( TString(observedHisto.at(i)).Contains("ee") ){
      //h_ofpred[i]->Scale(0.5/Rem);
      h_ofpred[i]->Scale(0.41);
      cout << "ee channel: scale em yield by 0.41" << endl;
      title     = (char*) "ee events";
      ee_and_mm = false;
      mysel = sel + ee;
    }

    else if( TString(observedHisto.at(i)).Contains("mm") ){
      //h_ofpred[i]->Scale(0.5*Rem);
      h_ofpred[i]->Scale(0.58);
      cout << "mm channel: scale em yield by 0.58" << endl;
      title     = (char*) "#mu#mu events";
      ee_and_mm = false;
      mysel = sel + mm;
    }

    else{
      h_ofpred[i]->Scale(0.99);
      cout << "ee+mm channels: scale em yield by 0.99" << endl;
      mysel = sel;
    }

    //------------------------------------------
    // WZ/ZZ histos
    //------------------------------------------

    int   nmetbins = 350;
    float metmin   = 0.0;
    float metmax   = 350.0;

    if( bveto ){
      nmetbins = 250;
      metmin   = 0.0;
      metmax   = 250.0;
    }

    h_wz[i]     = new TH1F(Form("h_wz_%i",i),Form("h_wz_%i",i),nmetbins,metmin,metmax);
    h_zz[i]     = new TH1F(Form("h_zz_%i",i),Form("h_wz_%i",i),nmetbins,metmin,metmax);

    h_wz[i]->Sumw2();
    h_zz[i]->Sumw2();

    TCanvas *ctemp = new TCanvas();
    ctemp->cd();
    chwz->Draw(Form("pfmet>>h_wz_%i",i),mysel*weight);
    chzz->Draw(Form("pfmet>>h_zz_%i",i),mysel*weight);
    delete ctemp;

    h_vz[i]     = (TH1F*) h_wz[i]->Clone(Form("h_vz_%i",i));
    h_vz[i]->Add(h_zz[i]);
    

    //------------------------------------------
    // rebin, style, sum of predictions
    //------------------------------------------

    h_ofpred[i]->Rebin(rebin);
    h_sf[i]->Rebin(rebin);
    h_gjets[i]->Rebin(rebin);
    h_wz[i]->Rebin(rebin);
    h_zz[i]->Rebin(rebin);
    h_vz[i]->Rebin(rebin);

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

    int mynbins = 6;
    if( bveto ) mynbins = 7;

    const unsigned int nbins = mynbins;

    //int bins[nbins]={0,30,60,100,150,200,300};
    // int      bins[nbins]  = {0,30,60,100,200,300};
    // Double_t xbins[nbins] = {0.0,30.0,60.0,100.0,200.0,300.0};

    int      bins[nbins];
    Double_t xbins[nbins+1];

    if( bveto ){
      bins[0] =   0;   xbins[0] =   0.0;
      bins[1] =  30;   xbins[1] =  30.0;
      bins[2] =  60;   xbins[2] =  60.0;
      bins[3] =  80;   xbins[3] =  80.0;
      bins[4] = 100;   xbins[4] = 100.0;
      bins[5] = 150;   xbins[5] = 150.0;
      bins[6] = 200;   xbins[6] = 200.0;
      xbins[7] = 250.0;
    }

    else{
      bins[0] =   0;   xbins[0] =   0.0;
      bins[1] =  30;   xbins[1] =  30.0;
      bins[2] =  60;   xbins[2] =  60.0;
      bins[3] = 100;   xbins[3] = 100.0;
      bins[4] = 200;   xbins[4] = 200.0;
      bins[5] = 300;   xbins[5] = 300.0;
      xbins[6] = 350.0;
    }

    int width1 = 18;
    int width2 =  4;

    cout << endl;
    cout << title << endl << endl;
    
    char* delim_start =   "|";
    char* delim_end   =   "|";
    char* delim       =   "|";
    char* pm          = "+/-";

    if( latex ){
      delim_start = " ";
      delim_end   = "  \\\\";
      delim       = "&";
      pm          = "$\\pm$";
    }

    // cout << "|" << setw(width1) << "" << setw(width2);
    // for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
    //   cout << "|" << setw(width1) << Form("MET>%i GeV",bins[ibin]) << setw(width2);
    // }
    // cout << "|" << endl;

    cout << delim_start << setw(width1) << "" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( latex ) cout << delim << setw(width1) << Form("\\MET\\ $>$ %i GeV",bins[ibin]) << setw(width2);
      else        cout << delim << setw(width1) << Form("MET>%i GeV",bins[ibin]) << setw(width2);
    }
    cout << delim_end << endl;
      
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

    //hsysterr[i] = new TH1F(Form("hsysterr_%i",i),Form("hsysterr_%i",i),nbins-1,xbins);
    hsysterr[i] = new TH1F(Form("hsysterr_%i",i),Form("hsysterr_%i",i),nbins,xbins);

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


      //float ofsyst = 0.1;
      //if( doKscaling && bins[ibin] >= 200 ) ofsyst = 0.25;

      float ofsyst = 1.0;
      float Ksyst  = 1.0;
      float Rsyst  = 1.0;

      if( bveto ){
	Ksyst = 0.02/0.13;
	if( bins[ibin] >= 150 ) Ksyst = 0.07/0.13;

	Rsyst = 0.06;
	if( !ee_and_mm ) Rsyst = 0.12;

	ofsyst = sqrt( Ksyst*Ksyst + Rsyst*Rsyst );
      }
      else{
	Ksyst = 0.02/0.14;
	if( bins[ibin] == 300 ) Ksyst = 0.08/0.14;

	Rsyst = 0.06;
	if( !ee_and_mm ) Rsyst = 0.12;

	ofsyst = sqrt( Ksyst*Ksyst + Rsyst*Rsyst );
      }

      // syst uncertainties
      nof_syst[ibin]    = ofsyst * h_ofpred[i]->Integral(bin,1000);   
      ngjets_syst[ibin] = 0.3 * h_gjets[i]->Integral(bin,1000);       // 30% uncertainty on Z+jets
      nwz_syst[ibin]    = 0.8 * h_wz[i]->Integral(bin,1000);          // 80% uncertainty on WZ
      nzz_syst[ibin]    = 0.8 * h_zz[i]->Integral(bin,1000);          // 80% uncertainty on ZZ
      ntot_syst[ibin]   = sqrt( pow(ngjets_syst[ibin],2) + pow(nof_syst[ibin],2) + pow(nwz_syst[ibin],2) + pow(nzz_syst[ibin],2));
      
      // stat uncertainties
      ngjets_stat[ibin] = histError( h_gjets[i]  , bin );
      nof_stat[ibin]    = histError( h_ofpred[i] , bin );
      nwz_stat[ibin]    = histError( h_wz[i]     , bin );
      nzz_stat[ibin]    = histError( h_zz[i]     , bin );
      ntot_stat[ibin]   = sqrt( pow(ngjets_stat[ibin],2) + pow(nof_stat[ibin],2) + pow(nwz_stat[ibin],2) + pow(nzz_stat[ibin],2));

      // ngjets_stat[ibin] = 0.0;
      // nof_stat[ibin]    = 0.0;
      // nwz_stat[ibin]    = 0.0;
      // nzz_stat[ibin]    = 0.0;
      // ntot_stat[ibin]   = 0.0;

      // tot uncertainties
      ngjets_toterr[ibin] = sqrt( pow(ngjets_syst[ibin],2) + pow(ngjets_stat[ibin],2) );
      nof_toterr[ibin]    = sqrt( pow(nof_syst[ibin],2) + pow(nof_stat[ibin],2) );
      nwz_toterr[ibin]    = sqrt( pow(nwz_syst[ibin],2) + pow(nwz_stat[ibin],2) );
      nzz_toterr[ibin]    = sqrt( pow(nzz_syst[ibin],2) + pow(nzz_stat[ibin],2) );
      ntot_toterr[ibin]   = sqrt( pow(ngjets_toterr[ibin],2) + pow(nof_toterr[ibin],2) + pow(nwz_toterr[ibin],2) + pow(nzz_toterr[ibin],2));

      excess[ibin]        = (ndata[ibin]-ntot[ibin])/sqrt( pow(ntot_toterr[ibin],2) + ndata[ibin]);

      //if( ibin+1 < nbins ){
      hsysterr[i]->SetBinContent(ibin+1,1);
      hsysterr[i]->SetBinError(ibin+1,ntot_toterr[ibin]/ntot[ibin]);
      //}
    }


    //-----------------------------
    // g+jets
    //-----------------------------

    cout << delim_start << setw(width1) << "\\zjets\\ bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( ngjets[ibin] > 100 ) cout << delim << setw(width1) << Form("%.0f %s %.0f",ngjets[ibin],pm,ngjets_toterr[ibin]) << setw(width2);
      else                     cout << delim << setw(width1) << Form("%.1f %s %.1f",ngjets[ibin],pm,ngjets_toterr[ibin]) << setw(width2);
    }
    cout << delim_end << endl;

    //-----------------------------
    // OF
    //-----------------------------

    cout << delim_start << setw(width1) << "FS bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( nof[ibin] > 100 )    cout << delim << setw(width1) << Form("%.0f %s %.0f",nof[ibin],pm,nof_toterr[ibin]) << setw(width2);
      else                     cout << delim << setw(width1) << Form("%.1f %s %.1f",nof[ibin],pm,nof_toterr[ibin]) << setw(width2);
    }
    cout << delim_end << endl;
    
    //-----------------------------
    // WZ
    //-----------------------------

    cout << delim_start << setw(width1) << "WZ bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << delim << setw(width1) << Form("%.1f %s %.1f",nwz[ibin],pm,nwz_toterr[ibin]) << setw(width2);
    }
    cout << delim_end << endl;
    
    //-----------------------------
    // ZZ
    //-----------------------------

    cout << delim_start << setw(width1) << "ZZ bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << delim << setw(width1) << Form("%.1f %s %.1f",nzz[ibin],pm,nzz_toterr[ibin]) << setw(width2);
      //cout << "|" << setw(width1) << Form("%.1f",nzz[ibin]) << setw(width2);
    }
    cout << delim_end << endl;

    //-----------------------------
    // total bkg
    //-----------------------------

    cout << delim_start << setw(width1) << "total bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      if( ntot[ibin] > 100 )   cout << delim << setw(width1) << Form("%.0f %s %.0f",ntot[ibin],pm,ntot_toterr[ibin]) << setw(width2);
      else                     cout << delim << setw(width1) << Form("%.1f %s %.1f",ntot[ibin],pm,ntot_toterr[ibin]) << setw(width2);

      //      cout << "|" << setw(width1) << Form("%.1f",ntot[ibin]) << setw(width2);
    }
    cout << delim_end << endl;

    //-----------------------------
    // data
    //-----------------------------

    cout << delim_start << setw(width1) << "data" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << delim << setw(width1) << ndata[ibin] << setw(width2);
    }
    cout << delim_end << endl;
    
    //-----------------------------
    // significance
    //-----------------------------

    cout << delim_start << setw(width1) << "significance" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << delim << setw(width1) << Form("%.1f",excess[ibin]) << setw(width2);
    }
    cout << delim_end << endl;

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

    if( bveto ) h_sf[i]->SetMinimum(0.01);
    else        h_sf[i]->SetMinimum(0.1);
    h_sf[i]->Draw("E1");
    pred[i]->Draw("histsame");
    h_sf[i]->Draw("sameE1");
    h_sf[i]->Draw("sameaxis");

    TLegend* leg = new TLegend(0.75,0.6,0.9,0.9);
    leg->AddEntry(h_sf[i],"data","lp");
    leg->AddEntry(h_gjets[i],"Z+jets","f");
    leg->AddEntry(h_ofpred[i],"FS","f");
    //leg->AddEntry(h_zz[i],"ZZ","f");
    //leg->AddEntry(h_wz[i],"WZ","f");
    leg->AddEntry(h_vz[i],"WZ+ZZ","f");
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
    
    if( bveto ){
      hratio[i]->SetMinimum(0.0);
      hratio[i]->SetMaximum(2.0);
      hratio[i]->GetYaxis()->SetRangeUser(0.0,2.0);
    }

    else{
      hratio[i]->SetMinimum(0.5);
      hratio[i]->SetMaximum(1.5);
      hratio[i]->GetYaxis()->SetRangeUser(0.5,1.5);
    }


    hratio[i]->Draw();
    hsysterr[i]->Draw("sameE2");
    hratio[i]->Draw("same");
    hratio[i]->Draw("axissame");
	
    TLine line;
    if( bveto ) line.DrawLine(0.0,1.0,250,1.0);
    else        line.DrawLine(0.0,1.0,350,1.0);

    if( printplots ){
      if( bveto ){
	can[i]->Print(Form("../plots/met_bveto_%i.pdf",i));
	can[i]->Print(Form("../plots/met_bveto_%i.root",i));
	can[i]->Print(Form("../plots/met_bveto_%i.C",i));
      }
      else{
	can[i]->Print(Form("../plots/met_%i.pdf",i));
	can[i]->Print(Form("../plots/met_%i.root",i));
	can[i]->Print(Form("../plots/met_%i.C",i));
      }
    }



    

  }
  
}
