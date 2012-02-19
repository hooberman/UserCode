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

void printCard( char* name , float sigtot , float kerr , char* version ){

  ofstream* ofile = new ofstream();

  ofile->open(Form("cards/%s/%s.txt",version,name));

  *ofile <<      "imax 1 number of channels"                                                            << endl;
  *ofile <<      "jmax 1 number of background"                                                          << endl;
  *ofile <<      "kmax * number of nuisance parameters"                                                 << endl;
  *ofile <<      "Observation 48                                                            "           << endl;
  *ofile << Form("shapes      *   * ../../rootfiles/%s/%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC" , version , name) << endl;
  *ofile << Form("shapes data_obs * ../../rootfiles/%s/%s.root  histo_Data" , version , name )                               << endl;
  *ofile <<      "bin                                  1       1"                                       << endl;
  *ofile << Form("process                      %s     bkg" , name )                                     << endl;
  *ofile <<      "process                              0       1"                                       << endl;
  *ofile << Form("rate                              %.1f    33.29" , sigtot)                            << endl;
  *ofile <<      "lumi                       lnN   1.045       -"                                       << endl;
  *ofile <<      "eff_leptons                lnN   1.050       -"                                       << endl;
  //*ofile << Form("k                          lnN    %.2f       -"  , kerr)                              << endl;
  //*ofile <<      "pdf                        lnN     1.2       -"                                       << endl;
  *ofile <<      "JES_shape                shape     1.0       -"                                       << endl;
  *ofile <<      "stat                 shapeStat       -     1.0"                                       << endl;
  *ofile <<      "syst                     shape       -     1.0"                                       << endl;
  
  ofile->close();

}


void makeCMSSMCards(){

  //---------------------------------------
  // load TChain
  //---------------------------------------
  
  TChain *ch = new TChain("t");
  ch->Add("output/V00-02-21/highpt/LMscanFall11dil_combined_smallTree.root");
  char* version = "V00-00-11";
  bool doSigCont = true;

  //---------------------------------------
  // selection
  //---------------------------------------

  //TCut weight   ("weight * 4.7 * ndavtxweight * trgeff * lepscale");
  //TCut weight   ("weight * 4.7 * ndavtxweight * trgeff * lepscale * ( 1 - sqrt(pow(ksusyup/ksusy-1,2)+0.2*0.2) )");
  //TCut weight("weight * 4.7 * ndavtxweight * trgeff * lepscale * ksusyup/ksusy");
  TCut weight("weight * 4.7 * ndavtxweight * trgeff * lepscale * ksusydn/ksusy");
  TCut presel("pfmet>50 && njets>=2 && ht>100 && !passz");
  TCut preselptll("pfmet>50 && njets>=2 && ht>100 && !passz && ( (leptype==2) || (leptype<2 && pfmet>75) )");
  TCut preseljup("pfmetUp>50   && njetsUp>=2   && htUp>100   && !passz");
  TCut preseljdn("pfmetDown>50 && njetsDown>=2 && htDown>100 && !passz");
  TCut SF("leptype==0 || leptype==1");
  TCut OF("leptype==2");

  TCut sig("(ht>300 && pfmet>275) || (ht>600 && pfmet>200)");  
  TCut sigdn("(htDown>300 && pfmetDown>275) || (htDown>600 && pfmetDown>200)");

  TCut SR1("ht>300 && ht<600 && pfmet>275");
  TCut SR2("ht>600 && pfmet>275");
  TCut SR3("ht>600 && pfmet>200 && pfmet<275");

  TCut SR1ptll("ht>300 && ht<600 && dilpt>275");
  TCut SR2ptll("ht>600 && dilpt>275");
  TCut SR3ptll("ht>600 && dilpt>200 && dilpt<275");

  TCut SR1jup("htUp>300 && htUp<600 && pfmetUp>275");
  TCut SR2jup("htUp>600 && pfmetUp>275");
  TCut SR3jup("htUp>600 && pfmetUp>200 && pfmetUp<275");

  TCut SR1ptlljup("htUp>300 && htUp<600 && dilpt>275");
  TCut SR2ptlljup("htUp>600 && dilpt>275");
  TCut SR3ptlljup("htUp>600 && dilpt>200 && dilpt<275");

  TCut SR1jdn("htDown>300 && htDown<600 && pfmetDown>275");
  TCut SR2jdn("htDown>600 && pfmetDown>275");
  TCut SR3jdn("htDown>600 && pfmetDown>200 && pfmetDown<275");

  TCut SR1ptlljdn("htDown>300 && htDown<600 && dilpt>275");
  TCut SR2ptlljdn("htDown>600 && dilpt>275");
  TCut SR3ptlljdn("htDown>600 && dilpt>200 && dilpt<275");

  //---------------------------------------
  // preselection and SR1,SR2,SR3 yields
  //---------------------------------------

  const int   nm0points    = 150;
  const float m0min        = 20.;
  const float m0max        = 3020.;
  const int   nm12points   = 38;
  const float m12min       = 20.;
  const float m12max       = 780.;

  const unsigned int nbins = 6;

  float KKC[nbins] = { 2.74 , 2.74 , 2.59 , 2.59 , 1.81 , 1.81 };

  TH2F* h[nbins];
  TH2F* hptll[nbins];
  TH2F* hjup[nbins];
  TH2F* hjdn[nbins];
  TH2F* hptlljup[nbins];
  TH2F* hptlljdn[nbins];
  TH2F* hkup[nbins];
  TH2F* hkdn[nbins];

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]        = new TH2F( Form("h_%i",i)        , Form("h_%i",i)           , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hjup[i]     = new TH2F( Form("hjup_%i",i)     , Form("hjup_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hjdn[i]     = new TH2F( Form("hjdn_%i",i)     , Form("hjdn_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hkup[i]     = new TH2F( Form("hkup_%i",i)     , Form("hkup_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hkdn[i]     = new TH2F( Form("hkdn_%i",i)     , Form("hkdn_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hptll[i]    = new TH2F( Form("hptll_%i",i)    , Form("hptll_%i",i)       , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hptlljup[i] = new TH2F( Form("hptlljup_%i",i) , Form("hptlljup_%i",i)    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hptlljdn[i] = new TH2F( Form("hptlljdn_%i",i) , Form("hptlljdn_%i",i)    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

    h[i]    ->Sumw2();
    hjup[i] ->Sumw2();
    hjdn[i] ->Sumw2();
    hkup[i] ->Sumw2();
    hkdn[i] ->Sumw2();

    hptll[i]->Sumw2();
    hptlljup[i]->Sumw2();
    hptlljdn[i]->Sumw2();
  }

  TH2F* hall    = new TH2F( "hall"    , "hall"    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hkupall = new TH2F( "hkupall" , "hkupall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hjdnall = new TH2F( "hjdnall" , "hjdnall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
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
  
  ch->Draw("m12:m0>>hptll_0"    , (preselptll + SR1ptll) * weight );
  ch->Draw("m12:m0>>hptll_1"    , (preselptll + SR1ptll) * weight );
  ch->Draw("m12:m0>>hptll_2"    , (preselptll + SR2ptll) * weight );
  ch->Draw("m12:m0>>hptll_3"    , (preselptll + SR2ptll) * weight );
  ch->Draw("m12:m0>>hptll_4"    , (preselptll + SR3ptll) * weight );
  ch->Draw("m12:m0>>hptll_5"    , (preselptll + SR3ptll) * weight );

  //JES up
  cout << "Filling JES up histos" << endl;
  ch->Draw("m12:m0>>hjup_0"     , (preseljup + SR1jup + SF) * weight );
  ch->Draw("m12:m0>>hjup_1"     , (preseljup + SR1jup + OF) * weight );
  ch->Draw("m12:m0>>hjup_2"     , (preseljup + SR2jup + SF) * weight );
  ch->Draw("m12:m0>>hjup_3"     , (preseljup + SR2jup + OF) * weight );
  ch->Draw("m12:m0>>hjup_4"     , (preseljup + SR3jup + SF) * weight );
  ch->Draw("m12:m0>>hjup_5"     , (preseljup + SR3jup + OF) * weight );
  
  //JES down
  cout << "Filling JES dn histos" << endl;
  ch->Draw("m12:m0>>hjdn_0"     , (preseljdn + SR1jdn + SF) * weight );
  ch->Draw("m12:m0>>hjdn_1"     , (preseljdn + SR1jdn + OF) * weight );
  ch->Draw("m12:m0>>hjdn_2"     , (preseljdn + SR2jdn + SF) * weight );
  ch->Draw("m12:m0>>hjdn_3"     , (preseljdn + SR2jdn + OF) * weight );
  ch->Draw("m12:m0>>hjdn_4"     , (preseljdn + SR3jdn + SF) * weight );
  ch->Draw("m12:m0>>hjdn_5"     , (preseljdn + SR3jdn + OF) * weight );
  ch->Draw("m12:m0>>hjdnall"    , (preseljdn + sigdn      ) * weight );
  
  //k-factor up
  // cout << "Filling k up histos" << endl;
  // ch->Draw("m12:m0>>hkup_0"     , (presel + SR1 + SF) * weightkup );
  // ch->Draw("m12:m0>>hkup_1"     , (presel + SR1 + OF) * weightkup );
  // ch->Draw("m12:m0>>hkup_2"     , (presel + SR2 + SF) * weightkup );
  // ch->Draw("m12:m0>>hkup_3"     , (presel + SR2 + OF) * weightkup );
  // ch->Draw("m12:m0>>hkup_4"     , (presel + SR3 + SF) * weightkup );
  // ch->Draw("m12:m0>>hkup_5"     , (presel + SR3 + OF) * weightkup );
  // ch->Draw("m12:m0>>hkupall"    , (presel + sig     ) * weightkup );
  
  //k-factor down
  // cout << "Filling k down histos" << endl;
  // ch->Draw("m12:m0>>hkdn_0"     , (presel + SR1 + SF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_1"     , (presel + SR1 + OF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_2"     , (presel + SR2 + SF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_3"     , (presel + SR2 + OF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_4"     , (presel + SR3 + SF) * weightkdn );
  // ch->Draw("m12:m0>>hkdn_5"     , (presel + SR3 + OF) * weightkdn );
  // ch->Draw("m12:m0>>hkdnall"    , (presel + sig     ) * weightkdn );

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

  ofstream* filelist = new ofstream();
  filelist->open(Form("cards/%s/file_list.txt",version));

  ofstream* filelist_CLs = new ofstream();
  filelist_CLs->open(Form("cards/%s/file_list_CLs.txt",version));

  //------------------------------------------
  // open histogram of #entries/CMSSM point
  //------------------------------------------

  // TFile* corfile = TFile::Open("mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0.root");
  // TH2F*  hscan   = (TH2F*) corfile->Get("hscan");
  
  // if( hscan == 0 ){
  //   cout << "Can't find TH2 hscan!!!" << endl;
  //   exit(0);
  // }

  //------------------------------------------
  // loop over CMSSM points
  //------------------------------------------

  int counter = 0;

  for( int m0bin = 1 ; m0bin <= hall->GetXaxis()->GetNbins() ; m0bin++ ){
    for( int m12bin = 1 ; m12bin <= hall->GetYaxis()->GetNbins() ; m12bin++ ){

      //------------------------------------------
      // require nentries = 10
      //------------------------------------------
      
      //int ngen = hscan->GetBinContent(m0bin,m12bin);
      //if( ngen != 10000 )  continue;
      
      int m0  = hall->GetXaxis()->GetBinCenter(m0bin);
      int m12 = hall->GetXaxis()->GetBinCenter(m12bin);

      //if( m0 > 1000 )                                  continue;
      if( hjdnall->GetBinContent(m0bin,m12bin) < 1e-10 ) continue;


      float nom = hall->GetBinContent(m0bin,m12bin);
      float kdn = hkupall->GetBinContent(m0bin,m12bin);
      float kup = hkdnall->GetBinContent(m0bin,m12bin);
	     
      float dup = kup / nom - 1;
      float ddn = 1 - kdn / nom;
      float kerr = 0.5 * ( dup + ddn );

      //---------------------------------------
      // make root file
      //---------------------------------------

      TH1F* histo_CMSSM               = new TH1F( Form("histo_CMSSM_%i_%i"              ,m0bin,m12bin) , Form("histo_CMSSM_%i_%i"              ,m0bin,m12bin) , nbins,0,nbins);
      TH1F* histo_CMSSM_JES_shapeUp   = new TH1F( Form("histo_CMSSM_%i_%i_JES_shapeUp"  ,m0bin,m12bin) , Form("histo_CMSSM_%i_%i_JES_shapeUp"  ,m0bin,m12bin) , nbins,0,nbins);
      TH1F* histo_CMSSM_JES_shapeDown = new TH1F( Form("histo_CMSSM_%i_%i_JES_shapeDown",m0bin,m12bin) , Form("histo_CMSSM_%i_%i_JES_shapeDown",m0bin,m12bin) , nbins,0,nbins);

      float sigtot    = 0;
      float sigtotjdn = 0;

      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
	float yieldnom = h[ibin]->GetBinContent(m0bin,m12bin);
	float yieldjup = hjup[ibin]->GetBinContent(m0bin,m12bin);
	float yieldjdn = hjdn[ibin]->GetBinContent(m0bin,m12bin);

	if( doSigCont ){
	  float ptll = hptll[ibin]->GetBinContent(m0bin,m12bin) * 0.5 * KKC[ibin];

	  if( ptll < yieldnom ){
	    float scale = (yieldnom-ptll)/yieldnom;
	    yieldnom -= ptll;
	    yieldjup *= scale;
	    yieldjdn *= scale;
	  }
	  
	  else{
	    yieldnom = 0.0;
	    yieldjup = 0.0;
	    yieldjdn = 0.0;
	  }

	}

	sigtot    += yieldnom;
	sigtotjdn += yieldjdn;

	histo_CMSSM->SetBinContent              ( ibin + 1 , yieldnom );
	histo_CMSSM_JES_shapeUp->SetBinContent  ( ibin + 1 , yieldjup );
	histo_CMSSM_JES_shapeDown->SetBinContent( ibin + 1 , yieldjdn );
      }

      if( sigtotjdn < 1e-10 ) continue;
      if( sigtot    < 2     ) continue;
      if( sigtot    > 80    ) continue;

      //float sigtot = hall->GetBinContent(m0bin,m12bin);

      counter++;

      char* fitoptions = "-M Hybrid --freq  --nToysForCLsb 1500 --nToysForCLb 500  --scanRs 1 -vR [0.2,5,x1.1]";
      //char* fitoptions = "-M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -rMin 0 -rMax 5";

      *doScript << Form("../../../../test/lands.exe -M Bayesian -d CMSSM_%i_%i.txt",m0bin,m12bin)         << endl;

      *doScript_CLs << Form("../../../../test/lands.exe -d CMSSM_%i_%i.txt %s  -n CMSSM_%i_%i",m0bin,m12bin,fitoptions,m0bin,m12bin) << endl;

       if( counter%4 == 0 ) *doScript_CLs1 << Form("../../../../test/lands.exe -d CMSSM_%i_%i.txt %s  -n CMSSM_%i_%i",m0bin,m12bin,fitoptions,m0bin,m12bin) << endl;

       if( counter%4 == 1 ) *doScript_CLs2 << Form("../../../../test/lands.exe -d CMSSM_%i_%i.txt %s  -n CMSSM_%i_%i",m0bin,m12bin,fitoptions,m0bin,m12bin) << endl;

       if( counter%4 == 2 ) *doScript_CLs3 << Form("../../../../test/lands.exe -d CMSSM_%i_%i.txt %s  -n CMSSM_%i_%i",m0bin,m12bin,fitoptions,m0bin,m12bin) << endl;

       if( counter%4 == 3 ) *doScript_CLs4 << Form("../../../../test/lands.exe -d CMSSM_%i_%i.txt %s  -n CMSSM_%i_%i",m0bin,m12bin,fitoptions,m0bin,m12bin) << endl;

      *filelist << Form("cards/%s/CMSSM_%i_%i.txt_Bayesian_bysObsLimit.root",version,m0bin,m12bin)        << endl;

      *filelist_CLs << Form("cards/%s/CMSSM_%i_%i_m2lnQ2.root",version,m0bin,m12bin)        << endl;

      printCard( Form("CMSSM_%i_%i",m0bin,m12bin) , sigtot , 1+kerr , version );
      
      //signal regions                         R1(SF)      R1(OF)    R2(SF)     R2(OF)      R3(SF)   R3(OF)
      int     data_yield[nbins]           = {    9     ,    10    ,    6     ,    5     ,     5    ,    13  };
      float   bkg_yield[nbins]            = {   5.7    ,    5.7   ,   5.2    ,   5.2    ,    5.6   ,   5.6  };
      float   bkg_syst[nbins]             = {   2.8    ,    2.8   ,   1.9    ,   1.9    ,    2.1   ,   2.1  };
      float   bkg_stat[nbins]             = {   5.1    ,    5.1   ,   4.1    ,   4.1    ,    3.4   ,   3.4  };

      TH1F* histo_Data = new TH1F("histo_Data","histo_Data",nbins,0,nbins);

      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
	histo_Data->SetBinContent( ibin+1 , data_yield[ibin] );
      }

      TH1F* histo_bkg               = new TH1F("histo_bkg","histo_bkg",nbins,0,nbins);
      TH1F* histo_bkg_statUp        = new TH1F("histo_bkg_statUp"  ,"histo_bkg_statUp"  ,nbins,0,nbins);
      TH1F* histo_bkg_statDown      = new TH1F("histo_bkg_statDown","histo_bkg_statDown",nbins,0,nbins);
      TH1F* histo_bkg_systUp        = new TH1F("histo_bkg_systUp"  ,"histo_bkg_systUp"  ,nbins,0,nbins);
      TH1F* histo_bkg_systDown      = new TH1F("histo_bkg_systDown","histo_bkg_systDown",nbins,0,nbins);
      
      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
	histo_bkg               -> SetBinContent(ibin+1,bkg_yield[ibin]);
	histo_bkg_statUp        -> SetBinContent(ibin+1, bkg_yield[ibin] + bkg_stat[ibin] );
	histo_bkg_statDown      -> SetBinContent(ibin+1, TMath::Max(bkg_yield[ibin] - bkg_stat[ibin],(float)0.0) );
	histo_bkg_systUp        -> SetBinContent(ibin+1, bkg_yield[ibin] + bkg_syst[ibin] );
	histo_bkg_systDown      -> SetBinContent(ibin+1, TMath::Max(bkg_yield[ibin] - bkg_syst[ibin],(float)0.0) );
      }


      TFile *f = TFile::Open( Form("rootfiles/%s/CMSSM_%i_%i.root",version,m0bin,m12bin) , "RECREATE");
      f->cd();
      histo_Data->Write();
      histo_bkg->Write();
      histo_bkg_statUp->Write();
      histo_bkg_statDown->Write();
      histo_bkg_systUp->Write();
      histo_bkg_systDown->Write();
      histo_CMSSM->Write();
      histo_CMSSM_JES_shapeUp->Write();
      histo_CMSSM_JES_shapeDown->Write();
      f->Close();

      delete histo_Data;
      delete histo_bkg;
      delete histo_bkg_statUp;
      delete histo_bkg_statDown;
      delete histo_bkg_systUp;
      delete histo_bkg_systDown;
      delete histo_CMSSM;
      delete histo_CMSSM_JES_shapeUp;
      delete histo_CMSSM_JES_shapeDown;
    }
  }

  doScript->close();
  filelist->close();

}
