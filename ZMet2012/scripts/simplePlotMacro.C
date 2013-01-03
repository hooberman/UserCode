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

//#include "mycolors.h"

enum metType        { e_tcmet    = 0 , e_tcmetNew      = 1 , e_pfmet             = 2 , e_t1pfmet = 3 , e_t1newpfmet = 4 };

//#include "histtools.h"
using namespace std;

bool     doKscaling   =    true;
float    K            =    0.14;
int      rebin        =      10;
bool     bveto        =    true;
bool     mjjTemplates =    true;
char*    mybvetochar  = "_bveto";
bool     pt40         =   false;
char*    signalRegion = "highMet";
float    xmin         =      -1;
bool     latex        =    true;
metType  myMetType    = e_pfmet; 
bool     normToLowMet =    true;
bool     exclusive    =    true;
bool     blind        =    false;
float    lumi         =    19.3;
bool     printCards   =    false;

//metType  myMetType  = e_t1newpfmet; 
//metType  myMetType  = e_t1pfmet; 
//float    Rem        =    1.18; 

void printCard( char* filename , int nobs , float nbkg , float errbkg );

float histError( TH1F* hist , int lowbin , int binhigh ){

  float err2 = 0;

  for( int ibin = lowbin ; ibin <= binhigh ; ibin++ ){
    err2 += pow( hist->GetBinError(ibin) , 2 );
  }

  return sqrt(err2);
}


void simplePlotMacro( bool printplots = false ){

  char* metvar = "";
  
  if( myMetType == e_pfmet ){
    metvar = "pfmet";
    cout << "Using pfmet out-of-the-box" << endl;
  }

  else if( myMetType == e_t1pfmet ){
    metvar = "pfmett1";
    cout << "Using pfmet type1 42X" << endl;
  }

  else if( myMetType == e_t1newpfmet ){
    metvar = "pfmett1new";
    cout << "Using pfmet type1 52X" << endl;
  }

  char* bvetochar = "";
  if( bveto ){
    bvetochar = mybvetochar;
    K = 0.13;
  }

  char* mjjTemplatesChar = "";
  if( mjjTemplates ){
    mjjTemplatesChar = "_mjjcut";
    cout << "Apply mjj cut in templates" << endl;
  }
  else{
    cout << "DON'T Apply mjj cut in templates" << endl;
  }

  char* pt40char = "";
  if( pt40 ){
    if( TString(signalRegion).Contains("lowMet") ){
      cout << "Using pT > 40 GeV jets, low MET signal region" << endl;
      //pt40char = "_pt40_lowMet";
      //pt40char = "_pt40_2012AB_lowMet";
      pt40char = "_pt40_2012C_lowMet";
      K = 0.14;
    }
    else if( TString(signalRegion).Contains("highMet") ){
      cout << "Using pT > 40 GeV jets, high MET signal region" << endl;
      //pt40char = "_pt40_highMet";
      //pt40char = "_pt40_2012AB_highMet";
      pt40char = "_pt40_2012C_highMet";
      K = 0.13;
    }
  }

  //-----------------------------------
  // data/MC files
  //-----------------------------------
  
  char* iter = "V00-02-00";

  TFile *f   = new TFile();
  //TFile *fwz = new TFile();
  //TFile *fzz = new TFile();

  TChain *chwz   = new TChain("T1");
  TChain *chzz   = new TChain("T1");
  TChain *chrare = new TChain("T1");

  chwz->Add(Form("../output/%s/wz_53X_baby.root",iter));
  chzz->Add(Form("../output/%s/zz_53X_baby.root",iter));
  chrare->Add(Form("../output/V00-01-05/ttZ_53X_baby.root",iter));
  chrare->Add(Form("../output/V00-01-05/VVV_53X_baby.root",iter));

  //-----------------------------------
  // selection
  //-----------------------------------

  TCut pfleptons("pflep1.pt() > 20 && pflep2.pt() > 20 && abs(lep1.pt()-pflep1.pt())<5.0 && abs(lep2.pt()-pflep2.pt())<5.0 ");
  TCut Zmass("dilmass>81 && dilmass<101");
  TCut Zmasspf("dilmasspf>81 && dilmasspf<101");
  TCut njets2("njets>=2");
  TCut ee("leptype==0 && (ee==1 || isdata==0)");
  TCut mm("leptype==1 && (mm==1 || isdata==0)");
  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut sf = ee||mm;
  TCut nu("ngennu>0");
  //TCut bvetocut("nbm==0");
  TCut bvetocut;
  if( TString(bvetochar).Contains("Loose")  ) bvetocut = TCut("nbcsvl==0");
  if( TString(bvetochar).Contains("Medium") ) bvetocut = TCut("nbcsvm==0");
  
  TCut njets3_40("njets40>=3");
  TCut njets2_40("njets40>=2");
  TCut ht100_40("ht40>=100.0");
  TCut nlep2("nlep==2");
  TCut mjj("mjj>70.0 && mjj<110.0");
  TCut pt40cuts("njets40>=2 && ht40>=100.0");
  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut pt2010("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut filters("csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1");
  TCut zdilep("zdilep==1");
  //TCut pfleptons2010("pflep1.pt() > 20 && pflep2.pt() > 10 && abs(lep1.pt()-pflep1.pt())<5.0 && abs(lep2.pt()-pflep2.pt())<5.0");

  TCut sel;
  sel += (ee||mm);
  sel += nu;
  sel += filters;
  sel += Zmass;

  if( pt40 ){
    if( TString(signalRegion).Contains("lowMet") ){
      sel += pt2020;
      sel += njets3_40;
    }
    else if( TString(signalRegion).Contains("highMet") ){
      sel += pt2010;
      sel += njets2_40;
      sel += ht100_40;
    }
  }

  else{
    sel += njets2;
    sel += pt2020;
  }

  if( bveto ){
    sel += bvetocut;
    sel += nlep2;
    sel += mjj;
  }

  TCut weight(Form("weight * %.1f * vtxweight * trgeff",lumi));
  //TCut weight("1");

  cout << "WZ/ZZ selection : " << sel.GetTitle() << endl;
  cout << "WZ/ZZ weight    : " << weight.GetTitle() << endl;

  //char* datafilename = (char*) Form("../output/%s/babylooper_dataskim2010_PhotonStitchedTemplate_%s%s%s_HT100.root",iter,metvar,bvetochar,pt40char);
  char* datafilename = (char*) Form("../output/%s/babylooper_data_53X_2012ALL_PhotonStitchedTemplate_%s%s%s.root",iter,metvar,bvetochar,mjjTemplatesChar,pt40char);

  cout << "Opening " << datafilename << endl;
  f   = TFile::Open(datafilename);

/*
  if( bveto ){
    // f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));
    // fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));
    // fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplatenjetsgeq2_bveto.root",iter));

    f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplate_%s_bveto.root",iter,metvar));
    cout << "Opening file " << Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplate_%s_bveto.root",iter,metvar) << endl;
    //fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplate_pfmet_bveto.root",iter));
    //fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplate_pfmet_bveto.root",iter));

    K = 0.13;
  }

  else{
    // f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplatenjetsgeq2.root",iter));
    // fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplatenjetsgeq2.root",iter));
    // fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplatenjetsgeq2.root",iter));

    f   = TFile::Open(Form("../output/%s/babylooper_dataskim2010_PhotonStitchedTemplate_%s.root",iter,metvar));
    cout << "Opening file " << Form("../output/%s/babylooper_dataskim2010_PhotonStitchedTemplate_%s.root",iter,metvar) << endl;

    //fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplate_t1pfmet.root",iter));
    //fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplate_t1pfmet.root",iter));

    // f   = TFile::Open(Form("../output/%s/babylooper_dataskim_PhotonStitchedTemplate_t1pfmet.root",iter));
    // fwz = TFile::Open(Form("../output/%s/babylooper_wz_PhotonStitchedTemplate_t1pfmet.root",iter));
    // fzz = TFile::Open(Form("../output/%s/babylooper_zz_PhotonStitchedTemplate_t1pfmet.root",iter));
  }
*/

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
  TH1F*    h_rare[nplots];
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

    if( blind ){
      int bin100 = h_sf[i]->FindBin(100);
      for( int ibin = bin100 ; ibin <= h_sf[i]->GetNbinsX() ; ++ibin ){
	h_sf[i]->SetBinContent(ibin,0);
      }
    }


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
      //h_ofpred[i]->Scale(0.41);
      //cout << "ee channel: scale em yield by 0.41" << endl;
      h_ofpred[i]->Scale(0.44);
      cout << "ee channel: scale em yield by 0.44" << endl;
      title     = (char*) "ee events";
      ee_and_mm = false;
      mysel = sel + ee;
    }

    else if( TString(observedHisto.at(i)).Contains("mm") ){
      //h_ofpred[i]->Scale(0.5*Rem);
      //h_ofpred[i]->Scale(0.58);
      //cout << "mm channel: scale em yield by 0.58" << endl;
      h_ofpred[i]->Scale(0.54);
      cout << "mm channel: scale em yield by 0.54" << endl;
      title     = (char*) "#mu#mu events";
      ee_and_mm = false;
      mysel = sel + mm;
    }

    else{
      h_ofpred[i]->Scale(0.98);
      cout << "ee+mm channels: scale em yield by 0.98" << endl;
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

    h_wz[i]     = new TH1F(Form("h_wz_%i",i)   , Form("h_wz_%i",i)   , nmetbins,metmin,metmax);
    h_zz[i]     = new TH1F(Form("h_zz_%i",i)   , Form("h_wz_%i",i)   , nmetbins,metmin,metmax);
    h_rare[i]   = new TH1F(Form("h_rare_%i",i) , Form("h_rare_%i",i) , nmetbins,metmin,metmax);

    h_wz[i]->Sumw2();
    h_zz[i]->Sumw2();
    h_rare[i]->Sumw2();

    TCut zzxsecweight("1.0");
    //TCut zzxsecweight("1.96");
    //cout << "Scaling ZZ yield by 1.96" << endl;

    TCanvas *ctemp = new TCanvas();
    ctemp->cd();
    chwz->  Draw(Form("%s>>h_wz_%i"   , metvar,i),mysel*weight);
    chzz->  Draw(Form("%s>>h_zz_%i"   , metvar,i),mysel*weight*zzxsecweight);
    chrare->Draw(Form("%s>>h_rare_%i" , metvar,i),(mysel+zdilep)*weight);
    delete ctemp;

    h_vz[i]     = (TH1F*) h_wz[i]->Clone(Form("h_vz_%i",i));
    h_vz[i]->Add(h_zz[i]);


    //-----------------------------------------------------------
    // normalize templates prediction to data in low MET region
    //-----------------------------------------------------------

    if( normToLowMet ){
      int bin60 = h_sf[i]->FindBin(60)-1;

      float ndata60     = h_sf[i]->Integral(1,bin60);
      float ngjets60    = h_gjets[i]->Integral(1,bin60);
      float nof60       = h_ofpred[i]->Integral(1,bin60);
      float nwz60       = h_wz[i]->Integral(1,bin60);
      float nzz60       = h_zz[i]->Integral(1,bin60);
      float nrare60     = h_rare[i]->Integral(1,bin60);

      cout << "Yields in 0-60 GeV region" << endl;
      cout << "data   : " << ndata60  << endl;
      cout << "gjets  : " << ngjets60 << endl;
      cout << "OF     : " << nof60    << endl;
      cout << "WZ     : " << nwz60    << endl;
      cout << "ZZ     : " << nzz60    << endl;
      cout << "Rare   : " << nrare60  << endl;

      float SF = ( ndata60 - nof60 - nwz60 - nzz60 - nrare60 ) / ngjets60;
      cout << "Scaling gjets by : " << SF << endl;
      h_gjets[i]->Scale(SF);
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
    h_rare[i]->Rebin(rebin);

    h_rare[i]->SetFillColor(kMagenta-5);
    h_gjets[i]->SetFillColor(50);
    h_ofpred[i]->SetFillColor(42);
    h_vz[i]->SetFillColor(31);

    cout << "SF events " << h_sf[i]->GetEntries() << endl;
    cout << "OF events " << h_ofpred[i]->GetEntries() << endl;

    pred[i] = new THStack();
    // pred[i]->Add(h_wz[i]);
    // pred[i]->Add(h_zz[i]);
    pred[i]->Add(h_rare[i]);
    pred[i]->Add(h_vz[i]);
    pred[i]->Add(h_ofpred[i]);
    pred[i]->Add(h_gjets[i]);

    htotpred[i] = (TH1F*) h_ofpred[i]->Clone(Form("htotpred_%i",i));
    htotpred[i]->Add(h_rare[i]);
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
    if( bveto ) mynbins = 8;

    const unsigned int nbins = mynbins;

    //int bins[nbins]={0,30,60,100,150,200,300};
    // int      bins[nbins]  = {0,30,60,100,200,300};
    // Double_t xbins[nbins] = {0.0,30.0,60.0,100.0,200.0,300.0};

    int      bins[nbins];
    Double_t xbins[nbins+1];

    if( bveto ){
      // bins[0] =   0;   xbins[0] =   0.0;
      // bins[1] =  30;   xbins[1] =  30.0;
      // bins[2] =  60;   xbins[2] =  60.0;
      // bins[3] =  80;   xbins[3] =  80.0;
      // bins[4] = 100;   xbins[4] = 100.0;
      // bins[5] = 150;   xbins[5] = 150.0;
      // bins[6] = 200;   xbins[6] = 200.0;
      // xbins[7] = 250.0;

      /*
      bins[0] =   0;   xbins[0] =   0.0;
      bins[1] =  30;   xbins[1] =  30.0;
      bins[2] =  60;   xbins[2] =  60.0;
      bins[3] =  80;   xbins[3] =  80.0;
      bins[4] = 100;   xbins[4] = 100.0;
      bins[5] = 120;   xbins[5] = 120.0;
      bins[6] = 140;   xbins[6] = 140.0;
      bins[7] = 160;   xbins[7] = 160.0;
      bins[8] = 180;   xbins[8] = 180.0;
      bins[9] = 200;   xbins[9] = 200.0;
      xbins[10] = 250.0;
      */

      bins[0] =   0;   xbins[0] =   0.0;
      bins[1] =  30;   xbins[1] =  30.0;
      bins[2] =  60;   xbins[2] =  60.0;
      bins[3] =  80;   xbins[3] =  80.0;
      bins[4] = 100;   xbins[4] = 100.0;
      bins[5] = 120;   xbins[5] = 120.0;
      bins[6] = 150;   xbins[6] = 150.0;
      bins[7] = 200;   xbins[7] = 200.0;
      xbins[8] = 250.0;


      /*
      bins[0] =  60;   xbins[0] =  60.0;
      bins[1] =  80;   xbins[1] =  80.0;
      bins[2] = 100;   xbins[2] = 100.0;
      bins[3] = 120;   xbins[3] = 120.0;
      bins[4] = 140;   xbins[4] = 140.0;
      bins[5] = 160;   xbins[5] = 160.0;
      bins[6] = 180;   xbins[6] = 180.0;
      bins[7] = 200;   xbins[7] = 200.0;
      xbins[8] = 250.0;
      */
    }

    else{
      bins[0] =   0;   xbins[0] =   0.0;
      bins[1] =  30;   xbins[1] =  30.0;
      bins[2] =  60;   xbins[2] =  60.0;
      bins[3] = 100;   xbins[3] = 100.0;
      if( pt40 && TString(signalRegion).Contains("highMet") ) { bins[4] = 150;   xbins[4] = 150.0; }
      else                                                    { bins[4] = 200;   xbins[4] = 200.0; }
      //bins[4] = 200;   xbins[4] = 200.0;
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
      if( exclusive && ibin < nbins - 1 ){
	if( latex ) cout << delim << setw(width1) << Form("\\MET\\ %i--%i GeV",bins[ibin],bins[ibin+1]) << setw(width2);
	else        cout << delim << setw(width1) << Form("MET %i-%i GeV",bins[ibin],bins[ibin+1]) << setw(width2);
      }

      else{
	if( latex ) cout << delim << setw(width1) << Form("\\MET\\ $>$ %i GeV",bins[ibin]) << setw(width2);
	else        cout << delim << setw(width1) << Form("MET>%i GeV",bins[ibin]) << setw(width2);
      }
    }
    cout << delim_end << endl;
      
    int   ndata[nbins];
    float ngjets[nbins];
    float nof[nbins];
    float nwz[nbins];
    float nzz[nbins];
    float nrare[nbins];
    float ntot[nbins];

    float ngjets_syst[nbins];
    float nof_syst[nbins];
    float nwz_syst[nbins];
    float nzz_syst[nbins];
    float nrare_syst[nbins];
    float ntot_syst[nbins];

    float ngjets_stat[nbins];
    float nof_stat[nbins];
    float nwz_stat[nbins];
    float nzz_stat[nbins];
    float nrare_stat[nbins];
    float ntot_stat[nbins];

    float ngjets_toterr[nbins];
    float nof_toterr[nbins];
    float nwz_toterr[nbins];
    float nzz_toterr[nbins];
    float nrare_toterr[nbins];
    float ntot_toterr[nbins];

    float excess[nbins];

    //hsysterr[i] = new TH1F(Form("hsysterr_%i",i),Form("hsysterr_%i",i),nbins-1,xbins);
    hsysterr[i] = new TH1F(Form("hsysterr_%i",i),Form("hsysterr_%i",i),nbins,xbins);

    for( unsigned int ibin = 0 ; ibin < nbins ; ++ibin ){
      int bin      = h_sf[i]->FindBin(bins[ibin]);
      int binhigh  = 1000;

      if( exclusive ){
	if( ibin == nbins - 1 ) binhigh = 1000;
	else                    binhigh = h_sf[i]->FindBin(bins[ibin+1]) - 1;
      }

      if( bin > h_sf[i]->GetNbinsX() ) bin = h_sf[i]->GetNbinsX();

      // values
      ndata[ibin]  = h_sf[i]->Integral(bin,binhigh);
      ngjets[ibin] = h_gjets[i]->Integral(bin,binhigh);
      nof[ibin]    = h_ofpred[i]->Integral(bin,binhigh);
      nwz[ibin]    = h_wz[i]->Integral(bin,binhigh);
      nzz[ibin]    = h_zz[i]->Integral(bin,binhigh);
      nrare[ibin]  = h_rare[i]->Integral(bin,binhigh);
      ntot[ibin]   = htotpred[i]->Integral(bin,binhigh);


      //float ofsyst = 0.1;
      //if( doKscaling && bins[ibin] >= 200 ) ofsyst = 0.25;

      float ofsyst = 1.0;
      float Ksyst  = 1.0;
      float Rsyst  = 1.0;

      if( !pt40 ){

	// nominal analysis with b-veto
	if( bveto ){
	  Ksyst = 0.02/0.13;
	  if( bins[ibin] >= 150 ) Ksyst = 0.05/0.13;
	  if( !doKscaling) Ksyst = 0.0;

	  Rsyst = 0.06;
	  if( !ee_and_mm ) Rsyst = 0.12;

	  ofsyst = sqrt( Ksyst*Ksyst + Rsyst*Rsyst );
	}

	// nominal analysis without b-veto
	else{
	  Ksyst = 0.02/0.14;
	  if( bins[ibin] == 200 ) Ksyst = 0.04/0.14;
	  if( bins[ibin] == 300 ) Ksyst = 0.08/0.14;
	  if( !doKscaling) Ksyst = 0.0;

	  Rsyst = 0.06;
	  if( !ee_and_mm ) Rsyst = 0.12;

	  ofsyst = sqrt( Ksyst*Ksyst + Rsyst*Rsyst );
	}

      }

      else{

	if( bveto ){
	  cout << "ERROR! not set up for bveto pt40 analysis" << endl;
	  exit(0);
	}

	// pt40 analysis without bveto
	else{

	  if( TString(signalRegion).Contains("lowMet") ){
	    Ksyst = 0.02/0.14;	  
	    if( bins[ibin] == 200 ) Ksyst = 0.03/0.14;
	    if( bins[ibin] == 300 ) Ksyst = 0.07/0.14;
	  }

	  else if( TString(signalRegion).Contains("highMet") ){
	    Ksyst = 0.02/0.13;	  
	    if( bins[ibin] == 200 ) Ksyst = 0.04/0.13;
	    if( bins[ibin] == 300 ) Ksyst = 0.05/0.13;
	  }

	  if( !doKscaling) Ksyst = 0.0;

	  Rsyst = 0.06;
	  if( !ee_and_mm ) Rsyst = 0.12;

	  ofsyst = sqrt( Ksyst*Ksyst + Rsyst*Rsyst );
	}
      }

      //cout << bins[ibin] << ": ofsyst " << ofsyst << endl;

      // syst uncertainties
      nof_syst[ibin]    = ofsyst * h_ofpred[i]->Integral(bin,binhigh);   
      ngjets_syst[ibin] = 0.3 * h_gjets[i]->Integral(bin,binhigh);       // 30% uncertainty on Z+jets
      nwz_syst[ibin]    = 0.7 * h_wz[i]->Integral(bin,binhigh);          // 80% uncertainty on WZ
      nzz_syst[ibin]    = 0.5 * h_zz[i]->Integral(bin,binhigh);          // 50% uncertainty on ZZ
      nrare_syst[ibin]  = 0.5 * h_rare[i]->Integral(bin,binhigh);        // 50% uncertainty on rare
      ntot_syst[ibin]   = sqrt( pow(ngjets_syst[ibin],2) + pow(nof_syst[ibin],2) + pow(nwz_syst[ibin],2) + pow(nzz_syst[ibin],2) + pow(nrare_syst[ibin],2));
      
      // stat uncertainties
      ngjets_stat[ibin] = histError( h_gjets[i]  , bin , binhigh );
      nof_stat[ibin]    = histError( h_ofpred[i] , bin , binhigh );
      nwz_stat[ibin]    = histError( h_wz[i]     , bin , binhigh );
      nzz_stat[ibin]    = histError( h_zz[i]     , bin , binhigh );
      nrare_stat[ibin]  = histError( h_rare[i]   , bin , binhigh );
      ntot_stat[ibin]   = sqrt( pow(ngjets_stat[ibin],2) + pow(nof_stat[ibin],2) + pow(nwz_stat[ibin],2) + pow(nzz_stat[ibin],2) + pow(nrare_stat[ibin],2));

      // ngjets_stat[ibin] = 0.0;
      // nof_stat[ibin]    = 0.0;
      // nwz_stat[ibin]    = 0.0;
      // nzz_stat[ibin]    = 0.0;
      // ntot_stat[ibin]   = 0.0;

      // tot uncertainties
      ngjets_toterr[ibin] = TMath::Min( sqrt( pow(ngjets_syst[ibin],2) + pow(ngjets_stat[ibin],2) ), ngjets[ibin] );
      nof_toterr[ibin]    = TMath::Min( sqrt( pow(nof_syst[ibin],2)   + pow(nof_stat[ibin],2) )   , nof[ibin] );
      nwz_toterr[ibin]    = TMath::Min( sqrt( pow(nwz_syst[ibin],2)   + pow(nwz_stat[ibin],2) )   , nwz[ibin] );
      nzz_toterr[ibin]    = TMath::Min( sqrt( pow(nzz_syst[ibin],2)   + pow(nzz_stat[ibin],2) )   , nzz[ibin] );
      nrare_toterr[ibin]  = TMath::Min( sqrt( pow(nrare_syst[ibin],2) + pow(nrare_stat[ibin],2) ) , nrare[ibin] );
      ntot_toterr[ibin]   = sqrt( pow(ngjets_toterr[ibin],2) + pow(nof_toterr[ibin],2) + pow(nwz_toterr[ibin],2) + pow(nzz_toterr[ibin],2)  + pow(nrare_toterr[ibin],2));

      excess[ibin]        = (ndata[ibin]-ntot[ibin])/sqrt( pow(ntot_toterr[ibin],2) + ndata[ibin]);

      //if( ibin+1 < nbins ){
      hsysterr[i]->SetBinContent(ibin+1,1);
      hsysterr[i]->SetBinError(ibin+1,ntot_toterr[ibin]/ntot[ibin]);
      //}

      //if( printCards ) printCard( Form("met%i%s.txt",bins[ibin],mybvetochar) , (int)ntot[ibin] , ntot[ibin] , ntot_toterr[ibin]/ntot[ibin] );
      if( printCards ) printCard( Form("met%i%s.txt",bins[ibin],mybvetochar) , ndata[ibin] , ntot[ibin] , ntot_toterr[ibin]/ntot[ibin] );
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
    // rare MC
    //-----------------------------

    cout << delim_start << setw(width1) << "rare SM bkg" << setw(width2);
    for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
      cout << delim << setw(width1) << Form("%.1f %s %.1f",nrare[ibin],pm,nrare_toterr[ibin]) << setw(width2);
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
      cout << delim << setw(width1) << Form("%.1f$\\sigma$",excess[ibin]) << setw(width2);
    }
    cout << delim_end << endl;




    /*
  //signal regions                          60-80      80-100    100-150    150-200  >200
  int     data_yield[nbins]           = {   47       , 7       , 6        , 2       , 0    };

  float   Zbkg_yield[nbins]           = {   32.9     , 5.2     , 1.7      , 0.4     , 0.20 };
  float   Zbkg_err[nbins]             = {   11.1     , 1.8     , 0.6      , 0.2     , 0.09 };

  float   OFbkg_yield[nbins]          = {   6.6      , 4.6     , 4.6      , 0.8     , 0.06 };
  float   OFbkg_err[nbins]            = {   1.6      , 1.2     , 1.2      , 0.3     , 0.07 };     

  float   VZbkg_yield[nbins]          = {   3.9      , 2.2     , 2.5      , 0.7     , 0.4  };
  float   VZbkg_err[nbins]            = {   2.0      , 1.1     , 1.3      , 0.4     , 0.2  };     
    */

    

    //-----------------------------
    // g+jets
    //-----------------------------
    
    width1 = 7;

    cout << endl << endl << endl;

    cout << "float Zbkg_yield[nbins]    = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",ngjets[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,ngjets[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    cout << "float Zbkg_err[nbins]      = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",ngjets_toterr[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,ngjets_toterr[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    //-----------------------------
    // FS
    //-----------------------------
    
    cout << "float OFbkg_yield[nbins]   = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nof[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nof[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    cout << "float OFbkg_err[nbins]     = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nof_toterr[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nof_toterr[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    //-----------------------------
    // WZ
    //-----------------------------
    
    cout << "float WZbkg_yield[nbins]   = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nwz[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nwz[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    cout << "float WZbkg_err[nbins]     = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nwz_toterr[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nwz_toterr[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    //-----------------------------
    // ZZ
    //-----------------------------
    
    cout << "float ZZbkg_yield[nbins]   = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nzz[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nzz[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    cout << "float ZZbkg_err[nbins]     = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nzz_toterr[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nzz_toterr[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    //-----------------------------
    // rare
    //-----------------------------
    
    cout << "float rarebkg_yield[nbins] = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nrare[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nrare[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    cout << "float rarebkg_err[nbins]   = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.1f , ",nrare_toterr[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.1f"   ,nrare_toterr[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    //-----------------------------
    // data
    //-----------------------------
    
    cout << "int   data_yield[nbins]    = { " << setw(width2);
    for( unsigned int ibin = 3 ; ibin < nbins ; ibin++ ){
      if( ibin < nbins - 1 ) cout << setw(width1) << Form("%.i , ",ndata[ibin]) << setw(width2);
      else                   cout << setw(width1) << Form("%.i"   ,ndata[ibin]) << setw(width2);
    }
    cout << "};" << endl;

    cout << endl << endl << endl;






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
    if( xmin < 0 ) mainpad[i]->SetLogy();

    h_sf[i]->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    h_sf[i]->GetYaxis()->SetTitle("entries / 10 GeV");
    h_sf[i]->GetYaxis()->SetTitleOffset(1.0);
    if( xmin > 0 ) h_sf[i]->GetXaxis()->SetRangeUser(xmin,350);

    if( bveto ) h_sf[i]->SetMinimum(0.01);
    else        h_sf[i]->SetMinimum(0.1);

    h_sf[i]->Draw("E1");
    pred[i]->Draw("histsame");
    h_sf[i]->Draw("sameE1");
    h_sf[i]->Draw("sameaxis");

    TLegend* leg = new TLegend(0.65,0.6,0.93,0.9);
    leg->AddEntry(h_sf[i],"data","lp");
    leg->AddEntry(h_gjets[i],"Z+jets","f");
    //leg->AddEntry(h_ofpred[i],"FS","f");
    leg->AddEntry(h_ofpred[i],"Flavor Symmetric","f");
    //leg->AddEntry(h_zz[i],"ZZ","f");
    //leg->AddEntry(h_wz[i],"WZ","f");
    leg->AddEntry(h_vz[i],"WZ+ZZ","f");
    leg->AddEntry(h_rare[i],"Rare SM","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    t->SetTextSize(0.035);
    t->DrawLatex(0.35,0.88,"CMS Preliminary");
    //t->DrawLatex(0.4,0.79,"#sqrt{s} = 8 TeV, L_{int} = 9.2 fb^{-1}");
    t->DrawLatex(0.35,0.83,Form("#sqrt{s} = 8 TeV, L_{int} = %.1f fb^{-1}",lumi));
    t->DrawLatex(0.35,0.78,title);

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
    if( xmin > 0 ) hratio[i]->GetXaxis()->SetRangeUser(xmin,350);
    
    if( bveto ){
      hratio[i]->SetMinimum(0.0);
      hratio[i]->SetMaximum(2.0);
      hratio[i]->GetYaxis()->SetRangeUser(0.0,2.0);
    }

    else{
      hratio[i]->SetMinimum(0.0);
      hratio[i]->SetMaximum(3.0);
      hratio[i]->GetYaxis()->SetRangeUser(0.0,2.0);
    }

    if( blind ){
      int bin100 = h_sf[i]->FindBin(100);
      for( int ibin = bin100 ; ibin <= h_sf[i]->GetNbinsX() ; ++ibin ){
	hratio[i]->SetBinContent(ibin,100);
	hratio[i]->SetBinError(ibin,0);
      }
    }


    hratio[i]->Draw();
    hsysterr[i]->Draw("sameE2");
    hratio[i]->Draw("same");
    hratio[i]->Draw("axissame");
	
    TLine line;

    if( xmin > 0 ){
      if( bveto ) line.DrawLine(xmin,1.0,250,1.0);
      else        line.DrawLine(xmin,1.0,350,1.0);
    }
    else{
      if( bveto ) line.DrawLine(0.0,1.0,250,1.0);
      else        line.DrawLine(0.0,1.0,350,1.0);
    }

    char* lep[3] = {"_all","_ee" ,"_mm"};

    if( printplots ){
      can[i]->Print(Form("../plots/%s%s%s%s.pdf" ,metvar,bvetochar,pt40char,lep[i]));
      can[i]->Print(Form("../plots/%s%s%s%s.C"   ,metvar,bvetochar,pt40char,lep[i]));
      can[i]->Print(Form("../plots/%s%s%s%s.root",metvar,bvetochar,pt40char,lep[i]));
    }



    

  }
  
}


void printCard( char* filename , int nobs , float nbkg , float errbkg ){

  ofstream* ofile = new ofstream(Form("cards/%s",filename),ios::trunc);

  *ofile << "imax 1  number of channels" << endl;
  *ofile << "jmax 1  number of backgrounds" << endl;
  *ofile << "kmax 2  number of nuisance parameters (sources of systematical uncertainties)" << endl;
  *ofile << "------------" << endl;
  *ofile << "# we have just one channel, in which we observe 0 events" << endl;
  *ofile << "bin         1" << endl;
  *ofile << Form("observation %i",nobs) << endl;
  *ofile << "------------" << endl;
  *ofile << "bin             1      1" << endl;
  *ofile << "process       signal  background" << endl;
  *ofile << "process         0      1" << endl;
  *ofile << Form("rate            1   %.2f",nbkg) << endl;
  *ofile << "------------" << endl;
  *ofile << "deltaS  lnN   1.10    -     uncertainty on signal" << endl;
  *ofile << Form("deltaB  lnN     -   %.2f    uncertainty on background",1.0+errbkg) << endl;

  ofile->close();

}
