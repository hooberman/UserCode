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

void printCard( char* name , float sigtot , float Ztot , float OFtot , float WZtot , float ZZtot , float raretot , int datatot , char* version ){

  ofstream* ofile = new ofstream();

  ofile->open(Form("cards/%s/%s.txt",version,name));

  *ofile <<      "imax 1 number of channels"                                                                  << endl;
  *ofile <<      "jmax 5 number of background"                                                                << endl;
  *ofile <<      "kmax * number of nuisance parameters"                                                       << endl;
  *ofile << Form("Observation %i                                                           ",datatot)         << endl;
  *ofile << Form("shapes      *   * ../../rootfiles/%s/%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC" , version , name) << endl;
  *ofile << Form("shapes data_obs * ../../rootfiles/%s/%s.root  histo_Data" , version , name )                << endl;
  //*ofile << Form("shapes      *   * %s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC" ,  name) << endl;
  //*ofile << Form("shapes data_obs * %s.root  histo_Data" ,  name )                << endl;
  *ofile <<      "bin                                  1        1      1      1     1     1"                  << endl;
  *ofile << Form("process                        %s     Zbkg  OFbkg  WZbkg  ZZbkg  rarebkg" , name )          << endl;
  *ofile <<      "process                              0        1      2      3     4     5"                  << endl;
  *ofile << Form("rate                              %.1f    %.1f    %.1f   %.1f   %.1f   %.1f" , sigtot,Ztot,OFtot,WZtot,ZZtot,raretot) << endl;
  *ofile <<      "lumi                       lnN   1.040       -       -      -     -     -"                  << endl;
  *ofile <<      "eff_leptons                lnN   1.050       -       -      -     -     -"                  << endl;
  *ofile <<      "btagerr                    lnN   1.100       -       -      -     -     -"                  << endl;
  *ofile <<      "JES_shape                shape     1.0       -       -      -     -     -"                  << endl;
  *ofile <<      "errZ                     shape       -     1.0       -      -     -     -"                  << endl;
  *ofile <<      "errOF                    shape       -       -     1.0      -     -     -"                  << endl;
  *ofile <<      "errWZ                    shape       -       -       -    1.0     -     -"                  << endl;
  *ofile <<      "errZZ                    shape       -       -       -      -   1.0     -"                  << endl;
  *ofile <<      "errRARE                  shape       -       -       -      -     -   1.0"                  << endl;
  
  ofile->close();

}


void makeSMSCards(){

  //---------------------------------------
  // load TChain
  //---------------------------------------
  
  TChain *ch = new TChain("T1");
  ch->Add("output/V00-02-13/wzsms_baby_oldIso.root");
  char* version = (char*) "V00-00-11";

  //---------------------------------------
  // load denominator histogram
  //---------------------------------------

  TFile* fdenom = TFile::Open("output/V00-02-13/wzsms_ngen.root");
  TH2F*  hdenom = (TH2F*) fdenom->Get("hmass");
  hdenom->Scale(10.0);

  //---------------------------------------
  // selection
  //---------------------------------------

  TCut weight   ("19500.0 * trgeff * vtxweight");
  //TCut weight   ("19500.0 * trgeff * vtxweight * (1./100000.)");
  //TCut weight   ("9.2 * trgeff * vtxweight * weight");

  // MEDIUM WP
  TCut presel   ("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvm==0 && mjj>70   && mjj<110   && nlep==2 && njets>=2     && leptype<2");
  TCut preseljup("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvm==0 && mjjup>70 && mjjup<110 && nlep==2 && njetsup>=2   && leptype<2");
  TCut preseljdn("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvm==0 && mjjdn>70 && mjjdn<110 && nlep==2 && njetsdn>=2   && leptype<2");

  /*
  // LOOSE WP
  TCut presel   ("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvl==0 && mjj>70   && mjj<110   && nlep==2 && njets>=2     && leptype<2");
  TCut preseljup("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvl==0 && mjjup>70 && mjjup<110 && nlep==2 && njetsup>=2   && leptype<2");
  TCut preseljdn("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvl==0 && mjjdn>70 && mjjdn<110 && nlep==2 && njetsdn>=2   && leptype<2");
  */

  const unsigned int nbins = 5;
  float metcuts[nbins+1] = {80,100,120,150,200,9999999};

  TCut sigall(Form("pfmet>%.0f"  ,metcuts[0]));
  TCut sigalldn(Form("pfmetdn>%.0f",metcuts[0]));

  cout << "Full signal region            " << sigall.GetTitle()   << endl;
  cout << "Full signal region (JES down) " << sigalldn.GetTitle() << endl;

  TCut sigcut[nbins];
  TCut sigcutup[nbins];
  TCut sigcutdn[nbins];

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    sigcut[i]   = TCut(Form("pfmet  >%.0f && pfmet<%.0f"   , metcuts[i] , metcuts[i+1]));
    sigcutup[i] = TCut(Form("pfmetup>%.0f && pfmetup<%.0f" , metcuts[i] , metcuts[i+1]));
    sigcutdn[i] = TCut(Form("pfmetdn>%.0f && pfmetdn<%.0f" , metcuts[i] , metcuts[i+1]));

    cout << "Region  : " << i << endl;
    cout << "Nominal : " << sigcut[i].GetTitle()   << endl;
    cout << "JES up  : " << sigcutup[i].GetTitle() << endl;
    cout << "JES dn  : " << sigcutdn[i].GetTitle() << endl << endl;
  }
  
  //---------------------------------------
  // preselection and signal region yields
  //---------------------------------------

  TH2F* h[nbins];
  TH2F* hjup[nbins];
  TH2F* hjdn[nbins];

  int   nx   =    41;
  float xmin =  -5.0;
  float xmax = 405.0;

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]        = new TH2F( Form("h_%i",i)        , Form("h_%i",i)           , nx,xmin,xmax,nx,xmin,xmax);
    hjup[i]     = new TH2F( Form("hjup_%i",i)     , Form("hjup_%i",i)        , nx,xmin,xmax,nx,xmin,xmax);
    hjdn[i]     = new TH2F( Form("hjdn_%i",i)     , Form("hjdn_%i",i)        , nx,xmin,xmax,nx,xmin,xmax);
    
    h[i]    ->Sumw2();
    hjup[i] ->Sumw2();
    hjdn[i] ->Sumw2();
  }
  
  TH2F* hall    = new TH2F( "hall"    , "hall"    , nx,xmin,xmax,nx,xmin,xmax);
  TH2F* hjdnall = new TH2F( "hjdnall" , "hjdnall" , nx,xmin,xmax,nx,xmin,xmax);
  
  hall->Sumw2();
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  cout << "Filling histos..." << endl;
  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
    ch->Draw(Form("ml:mg>>h_%i"    , ibin) , (presel    + sigcut[ibin]   ) * weight);
    ch->Draw(Form("ml:mg>>hjup_%i" , ibin) , (preseljup + sigcutup[ibin] ) * weight);
    ch->Draw(Form("ml:mg>>hjdn_%i" , ibin) , (preseljdn + sigcutdn[ibin] ) * weight);

    h[ibin]->Divide(hdenom);
    hjup[ibin]->Divide(hdenom);
    hjdn[ibin]->Divide(hdenom);
  }

  ch->Draw("ml:mg>>hall"        , (presel    + sigall     ) * weight);
  ch->Draw("ml:mg>>hjdnall"     , (preseljdn + sigalldn   ) * weight);

  hall->Divide(hdenom);
  hjdnall->Divide(hdenom);

  delete ctemp;

  ofstream* doScript = new ofstream();
  doScript->open(Form("cards/%s/doLimits.sh",version));

  ofstream* doScript_CLs = new ofstream();
  doScript_CLs->open(Form("cards/%s/doLimits_CLs.sh",version));

  ofstream* doScript_CLs1 = new ofstream();
  doScript_CLs1->open(Form("cards/%s/doLimits_CLs1.sh",version));

  ofstream* doScript_CLs2 = new ofstream();
  doScript_CLs2->open(Form("cards/%s/doLimits_CLs2.sh",version));

  ofstream* doScript_CLs3 = new ofstream();
  doScript_CLs3->open(Form("cards/%s/doLimits_CLs3.sh",version));

  ofstream* doScript_CLs4 = new ofstream();
  doScript_CLs4->open(Form("cards/%s/doLimits_CLs4.sh",version));

  ofstream* doScript_CLs_mL0 = new ofstream();
  doScript_CLs_mL0->open(Form("cards/%s/doLimits_mL0.sh",version));

  ofstream* doScript_CLs_mL50 = new ofstream();
  doScript_CLs_mL50->open(Form("cards/%s/doLimits_mL50.sh",version));

  ofstream* doScript_CLs_mL100 = new ofstream();
  doScript_CLs_mL100->open(Form("cards/%s/doLimits_mL100.sh",version));

  ofstream* filelist = new ofstream();
  filelist->open(Form("cards/%s/file_list.txt",version));

  ofstream* filelist_CLs = new ofstream();
  filelist_CLs->open(Form("cards/%s/file_list_CLs.txt",version));


  //---------------------------------------
  // make and fill data and bkg histos
  //---------------------------------------

  /*
  // MEDIUM WP
  //signal regions             80-100 100-120 120-150 150-200    >200
  float Zbkg_yield[nbins]    = { 40.9 ,  7.0 ,  3.1 ,  1.6 ,     0.8  };
  float Zbkg_err[nbins]      = { 12.4 ,  2.2 ,  0.9 ,  0.5 ,     0.3  };
  float OFbkg_yield[nbins]   = { 17.9 , 11.3 ,  6.9 ,  2.4 ,     0.4  };
  float OFbkg_err[nbins]     = {  3.3 ,  2.2 ,  1.5 ,  1.1 ,     0.3  };
  float WZbkg_yield[nbins]   = {  3.9 ,  2.1 ,  1.6 ,  1.0 ,     0.5  };
  float WZbkg_err[nbins]     = {  2.7 ,  1.5 ,  1.1 ,  0.7 ,     0.5  };
  float ZZbkg_yield[nbins]   = {  1.8 ,  1.0 ,  1.1 ,  0.8 ,     0.7  };
  float ZZbkg_err[nbins]     = {  0.9 ,  0.5 ,  0.6 ,  0.4 ,     0.7  };
  float rarebkg_yield[nbins] = {  0.3 ,  0.2 ,  0.3 ,  0.2 ,     0.2  };
  float rarebkg_err[nbins]   = {  0.2 ,  0.1 ,  0.1 ,  0.1 ,     0.2  };
  int   data_yield[nbins]    = {   56 ,   24 ,   16 ,    3 ,       1  };
  */

  /*
  // LOOSE WP
  //signal regions             80-100 100-120 120-150 150-200    >200
  float Zbkg_yield[nbins]    = { 29.7 ,  3.8 ,  2.2 ,  1.4 ,     0.5  };
  float Zbkg_err[nbins]      = {  9.1 ,  1.2 ,  0.7 ,  0.4 ,     0.2  };
  float OFbkg_yield[nbins]   = {  6.3 ,  5.0 ,  2.7 ,  1.4 ,     0.1  };
  float OFbkg_err[nbins]     = {  1.4 ,  1.2 ,  0.7 ,  0.7 ,     0.1  };
  float WZbkg_yield[nbins]   = {  2.6 ,  1.5 ,  1.0 ,  0.7 ,     0.3  };
  float WZbkg_err[nbins]     = {  1.8 ,  1.0 ,  0.7 ,  0.5 ,     0.3  };
  float ZZbkg_yield[nbins]   = {  1.4 ,  0.8 ,  0.8 ,  0.6 ,     0.5  };
  float ZZbkg_err[nbins]     = {  0.7 ,  0.4 ,  0.4 ,  0.3 ,     0.5  };
  float rarebkg_yield[nbins] = {  0.2 ,  0.1 ,  0.2 ,  0.2 ,     0.1  };
  float rarebkg_err[nbins]   = {  0.1 ,  0.1 ,  0.1 ,  0.1 ,     0.1  };
  int   data_yield[nbins]    = {   40 ,   10 ,   10 ,    2 ,       1  };
  */

  /*
  // MEDIUM WP, 19.3/fb RESULTS
  float Zbkg_yield[nbins]    = { 68.9 ,  7.8 ,  4.8 ,  2.1 ,     0.5  };
  float Zbkg_err[nbins]      = { 21.2 ,  2.5 ,  1.5 ,  0.7 ,     0.1  };
  float OFbkg_yield[nbins]   = { 35.2 , 21.9 , 13.2 ,  5.7 ,     0.8  };
  float OFbkg_err[nbins]     = {  6.2 ,  4.0 ,  2.5 ,  1.6 ,     0.4  };
  float WZbkg_yield[nbins]   = {  7.4 ,  4.0 ,  3.3 ,  2.0 ,     0.9  };
  float WZbkg_err[nbins]     = {  3.7 ,  2.0 ,  1.6 ,  1.0 ,     0.9  };
  float ZZbkg_yield[nbins]   = {  3.2 ,  1.9 ,  2.1 ,  1.5 ,     1.4  };
  float ZZbkg_err[nbins]     = {  1.6 ,  1.0 ,  1.1 ,  0.8 ,     1.4  };
  float rarebkg_yield[nbins] = {  0.9 ,  0.4 ,  0.9 ,  0.6 ,     0.4  };
  float rarebkg_err[nbins]   = {  0.5 ,  0.2 ,  0.5 ,  0.3 ,     0.4  };
  int   data_yield[nbins]    = {  115 ,   36 ,   25 ,   13 ,       4  };
  */

  // MEDIUM WP, 19.5/fb RESULTS
  float Zbkg_yield[nbins]    = { 64.5 ,  7.8 ,  3.7 ,  2.0 ,     0.4  };
  float Zbkg_err[nbins]      = { 22.2 ,  3.1 ,  1.6 ,  1.0 ,     0.3  };
  float OFbkg_yield[nbins]   = { 35.2 , 21.9 , 13.2 ,  5.7 ,     0.8  };
  float OFbkg_err[nbins]     = {  6.2 ,  4.0 ,  2.5 ,  1.6 ,     0.4  };
  float WZbkg_yield[nbins]   = {  7.4 ,  4.0 ,  3.3 ,  2.0 ,     0.9  };
  float WZbkg_err[nbins]     = {  3.7 ,  2.0 ,  1.6 ,  1.0 ,     0.9  };
  float ZZbkg_yield[nbins]   = {  3.2 ,  1.9 ,  2.1 ,  1.5 ,     1.4  };
  float ZZbkg_err[nbins]     = {  1.6 ,  1.0 ,  1.1 ,  0.8 ,     1.4  };
  float rarebkg_yield[nbins] = {  0.9 ,  0.4 ,  0.9 ,  0.6 ,     0.4  };
  float rarebkg_err[nbins]   = {  0.5 ,  0.2 ,  0.5 ,  0.3 ,     0.4  };
  int   data_yield[nbins]    = {  115 ,   36 ,   25 ,   13 ,       4  };

  int   data_tot  = 0;
  float Zbkg_tot  = 0;
  float OFbkg_tot = 0;
  float WZbkg_tot = 0;
  float ZZbkg_tot = 0;
  float rarebkg_tot = 0;

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    cout << "Bin      " << i << endl;
    cout << "Data     " << data_yield[i] << endl;
    cout << "Z  bkg   " << Form("%.1f +/- %.1f" ,Zbkg_yield[i]    , Zbkg_err[i])    << endl;
    cout << "OF bkg   " << Form("%.1f +/- %.1f" ,OFbkg_yield[i]   , OFbkg_err[i])   << endl;
    cout << "WZ bkg   " << Form("%.1f +/- %.1f" ,WZbkg_yield[i]   , WZbkg_err[i])   << endl;
    cout << "ZZ bkg   " << Form("%.1f +/- %.1f" ,ZZbkg_yield[i]   , ZZbkg_err[i])   << endl;
    cout << "rare bkg " << Form("%.1f +/- %.1f" ,rarebkg_yield[i] , rarebkg_err[i]) << endl;

    data_tot     += data_yield[i];
    Zbkg_tot     += Zbkg_yield[i];
    OFbkg_tot    += OFbkg_yield[i];
    WZbkg_tot    += WZbkg_yield[i];
    ZZbkg_tot    += ZZbkg_yield[i];
    rarebkg_tot  += rarebkg_yield[i];
  }

  cout << "Total data     " << data_tot    << endl;
  cout << "Total Z  bkg   " << Zbkg_tot    << endl;
  cout << "Total OF bkg   " << OFbkg_tot   << endl;
  cout << "Total WZ bkg   " << WZbkg_tot   << endl;
  cout << "Total ZZ bkg   " << ZZbkg_tot   << endl;
  cout << "Total rare bkg " << rarebkg_tot << endl;

  TH1F* histo_Data = new TH1F("histo_Data","histo_Data",nbins,0,nbins);

  TH1F* histo_Zbkg               = new TH1F("histo_Zbkg"                ,"histo_Zbkg"           ,nbins,0,nbins);
  TH1F* histo_Zbkg_errUp         = new TH1F("histo_Zbkg_errZUp"         ,"histo_Zbkg_errZUp"    ,nbins,0,nbins);
  TH1F* histo_Zbkg_errDown       = new TH1F("histo_Zbkg_errZDown"       ,"histo_Zbkg_errZDown"  ,nbins,0,nbins);

  TH1F* histo_OFbkg              = new TH1F("histo_OFbkg"               ,"histo_OFbkg"          ,nbins,0,nbins);
  TH1F* histo_OFbkg_errUp        = new TH1F("histo_OFbkg_errOFUp"       ,"histo_OFbkg_errOFUp"  ,nbins,0,nbins);
  TH1F* histo_OFbkg_errDown      = new TH1F("histo_OFbkg_errOFDown"     ,"histo_OFbkg_errOFDown",nbins,0,nbins);

  TH1F* histo_WZbkg              = new TH1F("histo_WZbkg"               ,"histo_WZbkg"          ,nbins,0,nbins);
  TH1F* histo_WZbkg_errUp        = new TH1F("histo_WZbkg_errWZUp"       ,"histo_WZbkg_errWZUp"  ,nbins,0,nbins);
  TH1F* histo_WZbkg_errDown      = new TH1F("histo_WZbkg_errWZDown"     ,"histo_WZbkg_errWZDown",nbins,0,nbins);

  TH1F* histo_ZZbkg              = new TH1F("histo_ZZbkg"               ,"histo_ZZbkg"          ,nbins,0,nbins);
  TH1F* histo_ZZbkg_errUp        = new TH1F("histo_ZZbkg_errZZUp"       ,"histo_ZZbkg_errZZUp"  ,nbins,0,nbins);
  TH1F* histo_ZZbkg_errDown      = new TH1F("histo_ZZbkg_errZZDown"     ,"histo_ZZbkg_errZZDown",nbins,0,nbins);
            
  TH1F* histo_rarebkg            = new TH1F("histo_rarebkg"             ,"histo_rarebkg"             ,nbins,0,nbins);
  TH1F* histo_rarebkg_errUp      = new TH1F("histo_rarebkg_errRAREUp"   ,"histo_rarebkg_errRAREUp"   ,nbins,0,nbins);
  TH1F* histo_rarebkg_errDown    = new TH1F("histo_rarebkg_errRAREDown" ,"histo_rarebkg_errRAREDown" ,nbins,0,nbins);
      
  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){

    histo_Data                -> SetBinContent(ibin+1, data_yield[ibin] );

    histo_Zbkg                -> SetBinContent(ibin+1, Zbkg_yield[ibin]);
    histo_Zbkg_errUp          -> SetBinContent(ibin+1, Zbkg_yield[ibin] + Zbkg_err[ibin] );
    histo_Zbkg_errDown        -> SetBinContent(ibin+1, TMath::Max(Zbkg_yield[ibin] - Zbkg_err[ibin],(float)0.0) );

    histo_OFbkg               -> SetBinContent(ibin+1, OFbkg_yield[ibin]);
    histo_OFbkg_errUp         -> SetBinContent(ibin+1, OFbkg_yield[ibin] + OFbkg_err[ibin] );
    histo_OFbkg_errDown       -> SetBinContent(ibin+1, TMath::Max(OFbkg_yield[ibin] - OFbkg_err[ibin],(float)0.0) );

    histo_WZbkg               -> SetBinContent(ibin+1, WZbkg_yield[ibin]);
    histo_WZbkg_errUp         -> SetBinContent(ibin+1, WZbkg_yield[ibin] + WZbkg_err[ibin] );
    histo_WZbkg_errDown       -> SetBinContent(ibin+1, TMath::Max(WZbkg_yield[ibin] - WZbkg_err[ibin],(float)0.0) );

    histo_ZZbkg               -> SetBinContent(ibin+1, ZZbkg_yield[ibin]);
    histo_ZZbkg_errUp         -> SetBinContent(ibin+1, ZZbkg_yield[ibin] + ZZbkg_err[ibin] );
    histo_ZZbkg_errDown       -> SetBinContent(ibin+1, TMath::Max(ZZbkg_yield[ibin] - ZZbkg_err[ibin],(float)0.0) );

    histo_rarebkg             -> SetBinContent(ibin+1, rarebkg_yield[ibin]);
    histo_rarebkg_errUp       -> SetBinContent(ibin+1, rarebkg_yield[ibin] + rarebkg_err[ibin] );
    histo_rarebkg_errDown     -> SetBinContent(ibin+1, TMath::Max(rarebkg_yield[ibin] - rarebkg_err[ibin],(float)0.0) );
  }


  //------------------------------------------
  // loop over SMS points
  //------------------------------------------

  int counter = 0;

  for( int mgbin = 1 ; mgbin <= hall->GetXaxis()->GetNbins() ; mgbin++ ){
    for( int mlbin = 1 ; mlbin <= hall->GetYaxis()->GetNbins() ; mlbin++ ){

      //if( !( mgbin == 17 && mlbin == 5 ) ) continue;

      int mg  = hall->GetXaxis()->GetBinCenter(mgbin);
      int ml  = hall->GetXaxis()->GetBinCenter(mlbin);

      // bool pass = false;

      // if( mg==150 && ml==0  ) pass = true;
      // if( mg==200 && ml==0  ) pass = true;
      // if( mg==250 && ml==0  ) pass = true;
      // if( mg==150 && ml==25 ) pass = true;
      // if( mg==200 && ml==50 ) pass = true;
      // if( mg==250 && ml==50 ) pass = true;
      // if( mg==200 && ml==80 ) pass = true;
      // if( mg==250 && ml==80 ) pass = true;
      
      // if( !pass ) continue;

      cout << endl;
      cout << "----------------------------------" << endl;
      cout << "mg    " << mg    << " ml    " << ml    << endl;
      cout << "mgbin " << mgbin << " mlbin " << mlbin << endl;
      cout << "----------------------------------" << endl;

      if( hjdnall->GetBinContent(mgbin,mlbin) < 1e-10 ) continue;

      float nom = hall->GetBinContent(mgbin,mlbin);

      //---------------------------------------
      // make signal histos
      //---------------------------------------

      TH1F* histo_SMS               = new TH1F( Form("histo_SMS_%i_%i"              ,mgbin,mlbin) , Form("histo_SMS_%i_%i"              ,mgbin,mlbin) , nbins,0,nbins);
      TH1F* histo_SMS_JES_shapeUp   = new TH1F( Form("histo_SMS_%i_%i_JES_shapeUp"  ,mgbin,mlbin) , Form("histo_SMS_%i_%i_JES_shapeUp"  ,mgbin,mlbin) , nbins,0,nbins);
      TH1F* histo_SMS_JES_shapeDown = new TH1F( Form("histo_SMS_%i_%i_JES_shapeDown",mgbin,mlbin) , Form("histo_SMS_%i_%i_JES_shapeDown",mgbin,mlbin) , nbins,0,nbins);

      float sigtot    = 0;
      float sigtotjdn = 0;

      cout << "Signal yield" << endl;
      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
	float yieldnom = h[ibin]->GetBinContent(mgbin,mlbin);
	float yieldjup = hjup[ibin]->GetBinContent(mgbin,mlbin);
	float yieldjdn = hjdn[ibin]->GetBinContent(mgbin,mlbin);

	sigtot    += yieldnom;
	sigtotjdn += yieldjdn;

	histo_SMS->SetBinContent              ( ibin + 1 , yieldnom );
	histo_SMS_JES_shapeUp->SetBinContent  ( ibin + 1 , yieldjup );
	histo_SMS_JES_shapeDown->SetBinContent( ibin + 1 , yieldjdn );
      
	cout << "Bin " << ibin << " " << yieldnom << endl;
      }

      cout << "Total signal yield " << sigtot << endl;

      if( sigtotjdn < 1e-10 ) continue;
      //if( sigtot    < 2     ) continue;

      counter++;

      *doScript << Form("../../../../test/lands.exe -M Bayesian -d SMS_%i_%i.txt",mgbin,mlbin)         << endl;

      *doScript_CLs << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      *filelist << Form("cards/%s/SMS_%i_%i.txt_Bayesian_bysObsLimit.root",version,mgbin,mlbin)        << endl;

      *filelist_CLs << Form("cards/%s/SMS_%i_%i_m2lnQ.root",version,mgbin,mlbin)        << endl;

      if( counter%4 == 0 ) *doScript_CLs1 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 1 ) *doScript_CLs2 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 2 ) *doScript_CLs3 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 3 ) *doScript_CLs4 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( mlbin==1 ) *doScript_CLs_mL0 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( mlbin==3 ) *doScript_CLs_mL50 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( mlbin==5 ) *doScript_CLs_mL100 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      printCard( Form("SMS_%i_%i",mgbin,mlbin) , sigtot , Zbkg_tot , OFbkg_tot , WZbkg_tot , ZZbkg_tot , rarebkg_tot , data_tot , version );


      TFile *f = TFile::Open( Form("rootfiles/%s/SMS_%i_%i.root",version,mgbin,mlbin) , "RECREATE");
      f->cd();
      histo_Data->Write();
      histo_Zbkg->Write();
      histo_Zbkg_errUp->Write();
      histo_Zbkg_errDown->Write();
      histo_OFbkg->Write();
      histo_OFbkg_errUp->Write();
      histo_OFbkg_errDown->Write();
      histo_WZbkg->Write();
      histo_WZbkg_errUp->Write();
      histo_WZbkg_errDown->Write();
      histo_ZZbkg->Write();
      histo_ZZbkg_errUp->Write();
      histo_ZZbkg_errDown->Write();
      histo_rarebkg->Write();
      histo_rarebkg_errUp->Write();
      histo_rarebkg_errDown->Write();
      histo_SMS->Write();
      histo_SMS_JES_shapeUp->Write();
      histo_SMS_JES_shapeDown->Write();
      f->Close();

      delete histo_SMS;
      delete histo_SMS_JES_shapeUp;
      delete histo_SMS_JES_shapeDown;
    }
  }

  doScript->close();
  doScript_CLs->close();
  doScript_CLs1->close();
  doScript_CLs2->close();
  doScript_CLs3->close();
  doScript_CLs4->close();
  doScript_CLs_mL0->close();
  doScript_CLs_mL50->close();
  filelist->close();
  filelist_CLs->close();

}
