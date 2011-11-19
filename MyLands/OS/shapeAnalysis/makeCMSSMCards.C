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

void printCard( char* name , float sigtot , float kerr ){

  ofstream* ofile = new ofstream();

  ofile->open(Form("cards/%s.txt",name));

  *ofile <<      "imax 1 number of channels"                                                            << endl;
  *ofile <<      "jmax 1 number of background"                                                          << endl;
  *ofile <<      "kmax * number of nuisance parameters"                                                 << endl;
  *ofile <<      "Observation 35                                                            "           << endl;
  *ofile << Form("shapes      *   * rootfiles/%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC",name) << endl;
  *ofile << Form("shapes data_obs * rootfiles/%s.root  histo_Data" , name )                             << endl;
  *ofile <<      "bin                                  1       1"                                       << endl;
  *ofile << Form("process                      %s     bkg" , name )                                     << endl;
  *ofile <<      "process                              0       1"                                       << endl;
  *ofile << Form("rate                              %.1f    25.5" , sigtot)                             << endl;
  *ofile <<      "lumi                       lnN   1.060       -"                                       << endl;
  *ofile <<      "eff_leptons                lnN   1.050       -"                                       << endl;
  *ofile << Form("k                          lnN    %.2f       -"  , kerr)                              << endl;
  *ofile <<      "JES_shape                shape     1.0       -"                                       << endl;
  *ofile <<      "stat                 shapeStat       -     1.0"                                       << endl;
  *ofile <<      "syst                     shape       -     1.0"                                       << endl;
  
  ofile->close();

}


float getYield( TChain *ch , TCut sel , TCut weight ){

  TH1F* hee = new TH1F("hee","hee",1,0,1);
  TH1F* hmm = new TH1F("hmm","hmm",1,0,1);
  TH1F* hem = new TH1F("hem","hem",1,0,1);

  hee->Sumw2();
  hmm->Sumw2();
  hem->Sumw2();

  TCut ee("leptype==0");
  TCut mm("leptype==1");
  TCut em("leptype==2");

  ch->Draw("0.5>>hee",(sel+ee)*weight);
  ch->Draw("0.5>>hmm",(sel+mm)*weight);
  ch->Draw("0.5>>hem",(sel+em)*weight);

  float nee = hee->GetBinContent(1);
  float nmm = hmm->GetBinContent(1);
  float nem = hem->GetBinContent(1);

  float eee = hee->GetBinError(1);
  float emm = hmm->GetBinError(1);
  float eem = hem->GetBinError(1);

  float tot    = 1.05 * nee + 1.08 * nem + 1.12 * nmm;
  float toterr = 1.05 * eee + 1.08 * eem + 1.12 * emm;

  cout << "Yield " << Form("%.1f +/- %.1f",tot,toterr) << endl;
  return tot;

}

void makeCMSSMCards(){

  //---------------------------------------
  // load TChain
  //---------------------------------------
  
  TChain *ch = new TChain("t");
  ch->Add("output/V00-02-07/highpt/LMscan_smallTree.root");

  //---------------------------------------
  // selection
  //---------------------------------------

  TCut weight   ("weight * 3.5 * ndavtxweight * trgeff * lepscale");
  TCut weightkup("weight * 3.5 * ndavtxweight * trgeff * lepscale * ksusyup/ksusy");
  TCut weightkdn("weight * 3.5 * ndavtxweight * trgeff * lepscale * ksusydn/ksusy");
  TCut presel("pfmet>50 && njets>=2 && ht>100 && !passz");
  TCut preseljup("pfmetUp>50   && njetsUp>=2   && htUp>100   && !passz");
  TCut preseljdn("pfmetDown>50 && njetsDown>=2 && htDown>100 && !passz");
  TCut SF("leptype==0 || leptype==1");
  TCut OF("leptype==2");

  TCut SR1("ht>300 && ht<600 && pfmet>275");
  TCut SR2("ht>600 && pfmet>275");
  TCut SR3("ht>600 && pfmet>200 && pfmet<275");
  TCut sig("(ht>300 && pfmet>275) || (ht>600 && pfmet>200)");

  TCut SR1jup("htUp>300 && htUp<600 && pfmetUp>275");
  TCut SR2jup("htUp>600 && pfmetUp>275");
  TCut SR3jup("htUp>600 && pfmetUp>200 && pfmetUp<275");

  TCut SR1jdn("htDown>300 && htDown<600 && pfmetDown>275");
  TCut SR2jdn("htDown>600 && pfmetDown>275");
  TCut SR3jdn("htDown>600 && pfmetDown>200 && pfmetDown<275");

  //---------------------------------------
  // preselection and SR1,SR2,SR3 yields
  //---------------------------------------

  const int   nm0points    = 100;
  const float m0min        = 20.;
  const float m0max        = 2020.;
  const int   nm12points   = 38;
  const float m12min       = 20.;
  const float m12max       = 780.;

  const unsigned int nbins = 6;

  TH2F* h[nbins];
  TH2F* hjup[nbins];
  TH2F* hjdn[nbins];
  TH2F* hkup[nbins];
  TH2F* hkdn[nbins];

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]      = new TH2F( Form("h_%i",i)    , Form("h_%i",i)    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hjup[i]   = new TH2F( Form("hjup_%i",i) , Form("hjup_%i",i) , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hjdn[i]   = new TH2F( Form("hjdn_%i",i) , Form("hjdn_%i",i) , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hkup[i]   = new TH2F( Form("hkup_%i",i) , Form("hkup_%i",i) , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hkdn[i]   = new TH2F( Form("hkdn_%i",i) , Form("hkdn_%i",i) , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

    h[i]    ->Sumw2();
    hjup[i] ->Sumw2();
    hjdn[i] ->Sumw2();
    hkup[i] ->Sumw2();
    hkdn[i] ->Sumw2();
  }

  TH2F* hall    = new TH2F( "hall"    , "hall"    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hkupall = new TH2F( "hkupall" , "hkupall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hkdnall = new TH2F( "hkdnall" , "hkdnall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  
  hall->Sumw2();
  hkupall->Sumw2();
  hkdnall->Sumw2();

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

 
  //nominal
  cout << "Filling nominal histos" << endl;
  ch->Draw("m12:m0>>h_0"        , (presel + SR1 + SF) * weight );
  ch->Draw("m12:m0>>h_1"        , (presel + SR1 + OF) * weight );
  ch->Draw("m12:m0>>h_2"        , (presel + SR2 + SF) * weight );
  ch->Draw("m12:m0>>h_3"        , (presel + SR2 + OF) * weight );
  ch->Draw("m12:m0>>h_4"        , (presel + SR3 + SF) * weight );
  ch->Draw("m12:m0>>h_5"        , (presel + SR3 + OF) * weight );
  ch->Draw("m12:m0>>hall"       , (presel + sig     ) * weight );
  
  // //JES up
  // cout << "Filling JES up histos" << endl;
  // ch->Draw("m12:m0>>hjup_0"     , (preseljup + SR1jup + SF) * weight );
  // ch->Draw("m12:m0>>hjup_1"     , (preseljup + SR1jup + OF) * weight );
  // ch->Draw("m12:m0>>hjup_2"     , (preseljup + SR2jup + SF) * weight );
  // ch->Draw("m12:m0>>hjup_3"     , (preseljup + SR2jup + OF) * weight );
  // ch->Draw("m12:m0>>hjup_4"     , (preseljup + SR3jup + SF) * weight );
  // ch->Draw("m12:m0>>hjup_5"     , (preseljup + SR3jup + OF) * weight );

  // //JES down
  // cout << "Filling JES dn histos" << endl;
  // ch->Draw("m12:m0>>hjdn_0"     , (preseljdn + SR1jdn + SF) * weight );
  // ch->Draw("m12:m0>>hjdn_1"     , (preseljdn + SR1jdn + OF) * weight );
  // ch->Draw("m12:m0>>hjdn_2"     , (preseljdn + SR2jdn + SF) * weight );
  // ch->Draw("m12:m0>>hjdn_3"     , (preseljdn + SR2jdn + OF) * weight );
  // ch->Draw("m12:m0>>hjdn_4"     , (preseljdn + SR3jdn + SF) * weight );
  // ch->Draw("m12:m0>>hjdn_5"     , (preseljdn + SR3jdn + OF) * weight );

  // //k-factor up
  // cout << "Filling k up histos" << endl;
  // ch->Draw("m12:m0>>hkup_0"     , (presel + SR1 + SF) * weightkup );
  // ch->Draw("m12:m0>>hkup_1"     , (presel + SR1 + OF) * weightkup );
  // ch->Draw("m12:m0>>hkup_2"     , (presel + SR2 + SF) * weightkup );
  // ch->Draw("m12:m0>>hkup_3"     , (presel + SR2 + OF) * weightkup );
  // ch->Draw("m12:m0>>hkup_4"     , (presel + SR3 + SF) * weightkup );
  // ch->Draw("m12:m0>>hkup_5"     , (presel + SR3 + OF) * weightkup );
  ch->Draw("m12:m0>>hkupall"    , (presel + sig     ) * weightkup );

  // //k-factor down
  // cout << "Filling k down histos" << endl;
  // ch->Draw("m12:m0>>hkdn_0"     , (presel + SR1 + SF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_1"     , (presel + SR1 + OF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_2"     , (presel + SR2 + SF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_3"     , (presel + SR2 + OF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_4"     , (presel + SR3 + SF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_5"     , (presel + SR3 + OF) * weightkdn );
  ch->Draw("m12:m0>>hkdnall"    , (presel + sig     ) * weightkdn );

  delete ctemp;


  for( int m0bin = 1 ; m0bin <= hall->GetXaxis()->GetNbins() ; m0bin++ ){
    for( int m12bin = 1 ; m12bin <= hall->GetYaxis()->GetNbins() ; m12bin++ ){

      if( !( m0bin == 10 && m12bin == 10 ) ) continue;

      float nom = hall->GetBinContent(m0bin,m12bin);
      float kdn = hkupall->GetBinContent(m0bin,m12bin);
      float kup = hkdnall->GetBinContent(m0bin,m12bin);
	     
      float dup = kup / nom - 1;
      float ddn = 1 - kdn / nom;
      float kerr = 0.5 * ( dup + ddn );

      float sigtot = hall->GetBinContent(m0bin,m12bin);
      printCard( Form("CMSSM_%i_%i",m0bin,m12bin) , sigtot , 1+kerr );

    }
  }


  // TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  // c1->cd();
  // gStyle->SetPaintTextFormat(".1f");
  // hall->Draw("colz");


  //int bin = h[0]->FindBin(80,400);
  //cout << "yield " << h[2]->GetBinContent(bin)+h[3]->GetBinContent(bin)+h[4]->GetBinContent(bin)+h[5]->GetBinContent(bin) << endl;




  //for( unsigned int m0

  //cout << "yield " << yield_highHT_SF << endl;

  /*
  TH1F* h = new TH1F("h","h",1,0,1);
  h->Sumw2();

f
  cout << endl;
  LM->Draw("0.5>>h",presel*myweight);
  cout << "presel   : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  cout << endl;
  LM->Draw("0.5>>h",(presel+SR1)*myweight);
  cout << "SR1      : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  LM->Draw("0.5>>h",(presel+SR2)*myweight);
  cout << "SR2      : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 
  float SR2_nominal = h->GetBinContent(1);

  LM->Draw("0.5>>h",(presel+SR3)*myweight);
  cout << "SR3      : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  LM->Draw("0.5>>h",(presel+(SR1||SR2||SR3))*myweight);
  cout << "allSR    : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  //---------------------------------------
  // yields in 6 bins for shape analysis
  //---------------------------------------

  TH1F* h1SF = new TH1F("h1SF","h1SF",1,0,1);
  TH1F* h1OF = new TH1F("h1OF","h1OF",1,0,1);
  TH1F* h2SF = new TH1F("h2SF","h2SF",1,0,1);
  TH1F* h2OF = new TH1F("h2OF","h2OF",1,0,1);
  TH1F* h3SF = new TH1F("h3SF","h3SF",1,0,1);
  TH1F* h3OF = new TH1F("h3OF","h3OF",1,0,1);

  h1SF->Sumw2();
  h1OF->Sumw2();
  h2SF->Sumw2();
  h2OF->Sumw2();
  h3SF->Sumw2();
  h3OF->Sumw2();

  //------nominal-------//
  LM->Draw("0.5>>h1SF",(presel+SR1+SF)*myweight);
  LM->Draw("0.5>>h1OF",(presel+SR1+OF)*myweight);
  LM->Draw("0.5>>h2SF",(presel+SR2+SF)*myweight);
  LM->Draw("0.5>>h2OF",(presel+SR2+OF)*myweight);
  LM->Draw("0.5>>h3SF",(presel+SR3+SF)*myweight);
  LM->Draw("0.5>>h3OF",(presel+SR3+OF)*myweight);
  
  cout << endl << endl;
  cout << "float   LM_yield[nbins]           = {  " 
       << Form("%.1f",h1SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h1OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3OF->GetBinContent(1)) << "  }" << endl;
  //cout << endl << endl;

  TH1F* hnominal = new TH1F("hnominal","hnominal",6,0,6);
  hnominal->Sumw2();
  hnominal->SetBinContent(1,h1SF->GetBinContent(1));  hnominal->SetBinError(1,h1SF->GetBinError(1));
  hnominal->SetBinContent(2,h1OF->GetBinContent(1));  hnominal->SetBinError(2,h1OF->GetBinError(1));
  hnominal->SetBinContent(3,h2SF->GetBinContent(1));  hnominal->SetBinError(3,h2SF->GetBinError(1));
  hnominal->SetBinContent(4,h2OF->GetBinContent(1));  hnominal->SetBinError(4,h2OF->GetBinError(1));
  hnominal->SetBinContent(5,h3SF->GetBinContent(1));  hnominal->SetBinError(5,h3SF->GetBinError(1));
  hnominal->SetBinContent(6,h3OF->GetBinContent(1));  hnominal->SetBinError(6,h3OF->GetBinError(1));

  //------JESup-------//
  TString preselupstring(presel);
  TString SR1upstring(SR1);
  TString SR2upstring(SR2);
  TString SR3upstring(SR3);

  preselupstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  preselupstring.ReplaceAll("ht"     ,  "htUp");
  preselupstring.ReplaceAll("njets"  ,  "njetsUp");

  SR1upstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  SR1upstring.ReplaceAll("ht"     ,  "htUp");
  SR1upstring.ReplaceAll("njets"  ,  "njetsUp");

  SR2upstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  SR2upstring.ReplaceAll("ht"     ,  "htUp");
  SR2upstring.ReplaceAll("njets"  ,  "njetsUp");

  SR3upstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  SR3upstring.ReplaceAll("ht"     ,  "htUp");
  SR3upstring.ReplaceAll("njets"  ,  "njetsUp");

  TCut preselup(preselupstring);
  TCut SR1up(SR1upstring);
  TCut SR2up(SR2upstring);
  TCut SR3up(SR3upstring);

  LM->Draw("0.5>>h1SF",(preselup+SR1up+SF)*myweight);
  LM->Draw("0.5>>h1OF",(preselup+SR1up+OF)*myweight);
  LM->Draw("0.5>>h2SF",(preselup+SR2up+SF)*myweight);
  LM->Draw("0.5>>h2OF",(preselup+SR2up+OF)*myweight);
  LM->Draw("0.5>>h3SF",(preselup+SR3up+SF)*myweight);
  LM->Draw("0.5>>h3OF",(preselup+SR3up+OF)*myweight);
  
  float SR2_up = h2SF->GetBinContent(1) + h2OF->GetBinContent(1);

  //cout << endl << endl;
  cout << "float   LM_yield_JESup[nbins]      = {  " 
       << Form("%.1f",h1SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h1OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3OF->GetBinContent(1)) << "  }" << endl;
  //cout << endl << endl;

  TH1F* hup = new TH1F("hup","hup",6,0,6);
  hup->Sumw2();
  hup->SetBinContent(1,h1SF->GetBinContent(1));  hup->SetBinError(1,h1SF->GetBinError(1));
  hup->SetBinContent(2,h1OF->GetBinContent(1));  hup->SetBinError(2,h1OF->GetBinError(1));
  hup->SetBinContent(3,h2SF->GetBinContent(1));  hup->SetBinError(3,h2SF->GetBinError(1));
  hup->SetBinContent(4,h2OF->GetBinContent(1));  hup->SetBinError(4,h2OF->GetBinError(1));
  hup->SetBinContent(5,h3SF->GetBinContent(1));  hup->SetBinError(5,h3SF->GetBinError(1));
  hup->SetBinContent(6,h3OF->GetBinContent(1));  hup->SetBinError(6,h3OF->GetBinError(1));

  //------JESdn-------//
  TString preseldnstring(presel);
  TString SR1dnstring(SR1);
  TString SR2dnstring(SR2);
  TString SR3dnstring(SR3);

  preseldnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  preseldnstring.ReplaceAll("ht"     ,  "htDown");
  preseldnstring.ReplaceAll("njets"  ,  "njetsDown");

  SR1dnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  SR1dnstring.ReplaceAll("ht"     ,  "htDown");
  SR1dnstring.ReplaceAll("njets"  ,  "njetsDown");

  SR2dnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  SR2dnstring.ReplaceAll("ht"     ,  "htDown");
  SR2dnstring.ReplaceAll("njets"  ,  "njetsDown");

  SR3dnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  SR3dnstring.ReplaceAll("ht"     ,  "htDown");
  SR3dnstring.ReplaceAll("njets"  ,  "njetsDown");

  TCut preseldn(preseldnstring);
  TCut SR1dn(SR1dnstring);
  TCut SR2dn(SR2dnstring);
  TCut SR3dn(SR3dnstring);

  LM->Draw("0.5>>h1SF",(preseldn+SR1dn+SF)*myweight);
  LM->Draw("0.5>>h1OF",(preseldn+SR1dn+OF)*myweight);
  LM->Draw("0.5>>h2SF",(preseldn+SR2dn+SF)*myweight);
  LM->Draw("0.5>>h2OF",(preseldn+SR2dn+OF)*myweight);
  LM->Draw("0.5>>h3SF",(preseldn+SR3dn+SF)*myweight);
  LM->Draw("0.5>>h3OF",(preseldn+SR3dn+OF)*myweight);

  float SR2_down = h2SF->GetBinContent(1) + h2OF->GetBinContent(1);

  TH1F* hdn = new TH1F("hdn","hdn",6,0,6);
  hdn->Sumw2();  
  hdn->SetBinContent(1,h1SF->GetBinContent(1));  hdn->SetBinError(1,h1SF->GetBinError(1));
  hdn->SetBinContent(2,h1OF->GetBinContent(1));  hdn->SetBinError(2,h1OF->GetBinError(1));
  hdn->SetBinContent(3,h2SF->GetBinContent(1));  hdn->SetBinError(3,h2SF->GetBinError(1));
  hdn->SetBinContent(4,h2OF->GetBinContent(1));  hdn->SetBinError(4,h2OF->GetBinError(1));
  hdn->SetBinContent(5,h3SF->GetBinContent(1));  hdn->SetBinError(5,h3SF->GetBinError(1));
  hdn->SetBinContent(6,h3OF->GetBinContent(1));  hdn->SetBinError(6,h3OF->GetBinError(1));

  //cout << endl << endl;
  cout << "float   LM_yield_JESdn[nbins]      = {  " 
       << Form("%.1f",h1SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h1OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3OF->GetBinContent(1)) << "  }" << endl;
  cout << endl << endl;


  cout << endl;
  cout << "SR2 nominal " << SR2_nominal << endl;
  cout << "SR2 up      " << SR2_up      << endl;
  cout << "SR2 down    " << SR2_down    << endl;

  float dup = SR2_up/SR2_nominal-1;
  float ddn = 1-SR2_down/SR2_nominal;
  float dave = 0.5 * ( dup + ddn );

  cout << "dup                " << dup << endl;
  cout << "ddn                " << ddn << endl;
  cout << "JES uncertainty:   " << dave << endl;

  float tot = sqrt( 0.06*0.06 + 0.05*0.05 + dave * dave );
  cout << "tot uncertainty:   " << tot << endl;

  TCanvas *c1 = new TCanvas();
  c1->cd();

  hnominal->SetLineWidth(2);
  hup->SetLineColor(2);
  hup->SetLineStyle(2);
  hdn->SetLineColor(4);
  hdn->SetLineStyle(2);

  hup->GetXaxis()->SetBinLabel(1,"SR1 SF");
  hup->GetXaxis()->SetBinLabel(2,"SR1 OF");
  hup->GetXaxis()->SetBinLabel(3,"SR2 SF");
  hup->GetXaxis()->SetBinLabel(4,"SR2 OF");
  hup->GetXaxis()->SetBinLabel(5,"SR3 SF");
  hup->GetXaxis()->SetBinLabel(6,"SR3 OF");
  hup->GetYaxis()->SetTitle("Signal Region Yield");

  hup->Draw("hist");
  hdn->Draw("samehist");
  hnominal->Draw("samehist");

  TLegend *leg = new TLegend(0.65,0.7,0.95,0.9);
  leg->AddEntry(hnominal,"nominal","l");
  leg->AddEntry(hup,"JES +7.5%","l");
  leg->AddEntry(hdn,"JES -7.5%","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  c1->Print("shapeYields.pdf");
  c1->Print("shapeYields.png");
  
  */
}




  // cout << endl;
  // cout << "SR1 SF   : " << Form("%.1f +/- %.1f",h1SF->GetBinContent(1),h1SF->GetBinError(1)) << endl;
  // cout << "SR1 OF   : " << Form("%.1f +/- %.1f",h1OF->GetBinContent(1),h1OF->GetBinError(1)) << endl;
  // cout << "SR2 SF   : " << Form("%.1f +/- %.1f",h2SF->GetBinContent(1),h2SF->GetBinError(1)) << endl;
  // cout << "SR2 OF   : " << Form("%.1f +/- %.1f",h2OF->GetBinContent(1),h2OF->GetBinError(1)) << endl;
  // cout << "SR3 SF   : " << Form("%.1f +/- %.1f",h3SF->GetBinContent(1),h3SF->GetBinError(1)) << endl;
  // cout << "SR3 OF   : " << Form("%.1f +/- %.1f",h3OF->GetBinContent(1),h3OF->GetBinError(1)) << endl;







/*




  TH2F* h_SR1_SF = new TH2F("h_SR1_SF"      , "h_SR1_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR1_OF = new TH2F("h_SR1_OF"      , "h_SR1_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2_SF = new TH2F("h_SR2_SF"      , "h_SR2_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2_OF = new TH2F("h_SR2_OF"      , "h_SR2_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3_SF = new TH2F("h_SR3_SF"      , "h_SR3_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3_OF = new TH2F("h_SR3_OF"      , "h_SR3_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  TH2F* h_SR1jup_SF = new TH2F("h_SR1jup_SF"      , "h_SR1jup_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR1jup_OF = new TH2F("h_SR1jup_OF"      , "h_SR1jup_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2jup_SF = new TH2F("h_SR2jup_SF"      , "h_SR2jup_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2jup_OF = new TH2F("h_SR2jup_OF"      , "h_SR2jup_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3jup_SF = new TH2F("h_SR3jup_SF"      , "h_SR3jup_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3jup_OF = new TH2F("h_SR3jup_OF"      , "h_SR3jup_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  TH2F* h_SR1jdn_SF = new TH2F("h_SR1jdn_SF"      , "h_SR1jdn_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR1jdn_OF = new TH2F("h_SR1jdn_OF"      , "h_SR1jdn_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2jdn_SF = new TH2F("h_SR2jdn_SF"      , "h_SR2jdn_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2jdn_OF = new TH2F("h_SR2jdn_OF"      , "h_SR2jdn_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3jdn_SF = new TH2F("h_SR3jdn_SF"      , "h_SR3jdn_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3jdn_OF = new TH2F("h_SR3jdn_OF"      , "h_SR3jdn_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  TH2F* h_SR1kup_SF = new TH2F("h_SR1kup_SF"      , "h_SR1kup_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR1kup_OF = new TH2F("h_SR1kup_OF"      , "h_SR1kup_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2kup_SF = new TH2F("h_SR2kup_SF"      , "h_SR2kup_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2kup_OF = new TH2F("h_SR2kup_OF"      , "h_SR2kup_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3kup_SF = new TH2F("h_SR3kup_SF"      , "h_SR3kup_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3kup_OF = new TH2F("h_SR3kup_OF"      , "h_SR3kup_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  TH2F* h_SR1kdn_SF = new TH2F("h_SR1kdn_SF"      , "h_SR1kdn_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR1kdn_OF = new TH2F("h_SR1kdn_OF"      , "h_SR1kdn_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2kdn_SF = new TH2F("h_SR2kdn_SF"      , "h_SR2kdn_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR2kdn_OF = new TH2F("h_SR2kdn_OF"      , "h_SR2kdn_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3kdn_SF = new TH2F("h_SR3kdn_SF"      , "h_SR3kdn_SF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* h_SR3kdn_OF = new TH2F("h_SR3kdn_OF"      , "h_SR3kdn_OF" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  h_SR1_SF->Sumw2();
  h_SR1_OF->Sumw2();
  h_SR2_SF->Sumw2();
  h_SR2_OF->Sumw2();
  h_SR3_SF->Sumw2();
  h_SR3_OF->Sumw2();

  h_SR1jup_SF->Sumw2();
  h_SR1jup_OF->Sumw2();
  h_SR2jup_SF->Sumw2();
  h_SR2jup_OF->Sumw2();
  h_SR3jup_SF->Sumw2();
  h_SR3jup_OF->Sumw2();

  h_SR1jdn_SF->Sumw2();
  h_SR1jdn_OF->Sumw2();
  h_SR2jdn_SF->Sumw2();
  h_SR2jdn_OF->Sumw2();
  h_SR3jdn_SF->Sumw2();
  h_SR3jdn_OF->Sumw2();

  h_SR1kup_SF->Sumw2();
  h_SR1kup_OF->Sumw2();
  h_SR2kup_SF->Sumw2();
  h_SR2kup_OF->Sumw2();
  h_SR3kup_SF->Sumw2();
  h_SR3kup_OF->Sumw2();

  h_SR1kdn_SF->Sumw2();
  h_SR1kdn_OF->Sumw2();
  h_SR2kdn_SF->Sumw2();
  h_SR2kdn_OF->Sumw2();
  h_SR3kdn_SF->Sumw2();
  h_SR3kdn_OF->Sumw2();
*/
