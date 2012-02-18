#include "babylooper.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TChain.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"


using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}


void printCard( char* name , float sigtot , char* version ){

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


void babylooper::ScanChain (TChain* chain, const char* version, const char* prefix, bool isData){

  char* outversion = (char*) "V00-00-09";
  bool doSigCont   = true;

  TObjArray *listOfFiles = chain->GetListOfFiles();
  unsigned int nEventsTotal = 0;

  //---------------------------------------
  // declare histograms
  //---------------------------------------

  const int   nm0points    = 100;
  const float m0min        = 20.;
  const float m0max        = 2020.;
  const int   nm12points   = 38;
  const float m12min       = 20.;
  const float m12max       = 780.;

  const unsigned int nbins = 6;

  float KKC[nbins] = { 2.74 , 2.74 , 2.59 , 2.59 , 1.81 , 1.81 };

  TH2F* h[nbins];
  TH2F* hptll[nbins];
  TH2F* hjup[nbins];
  TH2F* hjdn[nbins];
  //TH2F* hptlljup[nbins];
  //TH2F* hptlljdn[nbins];
  TH2F* hkup[nbins];
  TH2F* hkdn[nbins];

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]        = new TH2F( Form("h_%i",i)        , Form("h_%i",i)           , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hjup[i]     = new TH2F( Form("hjup_%i",i)     , Form("hjup_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hjdn[i]     = new TH2F( Form("hjdn_%i",i)     , Form("hjdn_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hkup[i]     = new TH2F( Form("hkup_%i",i)     , Form("hkup_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hkdn[i]     = new TH2F( Form("hkdn_%i",i)     , Form("hkdn_%i",i)        , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    hptll[i]    = new TH2F( Form("hptll_%i",i)    , Form("hptll_%i",i)       , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    //hptlljup[i] = new TH2F( Form("hptlljup_%i",i) , Form("hptlljup_%i",i)    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
    //hptlljdn[i] = new TH2F( Form("hptlljdn_%i",i) , Form("hptlljdn_%i",i)    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

    h[i]    ->Sumw2();
    hjup[i] ->Sumw2();
    hjdn[i] ->Sumw2();
    hkup[i] ->Sumw2();
    hkdn[i] ->Sumw2();

    hptll[i]->Sumw2();
    //hptlljup[i]->Sumw2();
    //hptlljdn[i]->Sumw2();
  }

  TH2F* hall    = new TH2F( "hall"    , "hall"    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hkupall = new TH2F( "hkupall" , "hkupall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hjdnall = new TH2F( "hjdnall" , "hjdnall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hkdnall = new TH2F( "hkdnall" , "hkdnall" , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  
  hall->Sumw2();
  hkupall->Sumw2();
  hkdnall->Sumw2();


  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())){
    
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("t");

    setBranches(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();
 
    for (unsigned int event = 0 ; event < nEvents; ++event){
        
      tree->GetEntry(event);
      ++nEventsTotal;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
            
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEvents*0.01));
          fflush(stdout);
        }
      }
      
      float weight      = weight_ * 4.7 * ndavtxweight_ * trgeff_ * lepscale_;

      bool presel	= pfmet_>50 && njets_>=2 && ht_>100 && passz_==0;
      bool preselptll	= pfmet_>50 && njets_>=2 && ht_>100 && passz_==0 && ( (leptype_==2) || (leptype_<2 && pfmet_>75) );
      bool preseljup	= pfmetUp_>50   && njetsUp_>=2   && htUp_>100   && passz_==0;
      bool preseljdn	= pfmetDown_>50 && njetsDown_>=2 && htDown_>100 && passz_==0;
      bool SF		= leptype_==0 || leptype_==1;
      bool OF		= leptype_==2;

      bool sig		= (ht_>300 && pfmet_>275) || (ht_>600 && pfmet_>200);  
      bool sigdn	= (htDown_>300 && pfmetDown_>275) || (htDown_>600 && pfmetDown_>200);

      bool SR1		= ht_>300 && ht_<600 && pfmet_>275;
      bool SR2		= ht_>600 && pfmet_>275;
      bool SR3		= ht_>600 && pfmet_>200 && pfmet_<275;

      bool SR1ptll	= ht_>300 && ht_<600 && dilpt_>275;
      bool SR2ptll	= ht_>600 && dilpt_>275;
      bool SR3ptll	= ht_>600 && dilpt_>200 && dilpt_<275;

      bool SR1jup	= htUp_>300 && htUp_<600 && pfmetUp_>275;
      bool SR2jup	= htUp_>600 && pfmetUp_>275;
      bool SR3jup	= htUp_>600 && pfmetUp_>200 && pfmetUp_<275;

      // bool SR1ptlljup	= htUp_>300 && htUp_<600 && dilpt_>275;
      // bool SR2ptlljup	= htUp_>600 && dilpt_>275;
      // bool SR3ptlljup	= htUp_>600 && dilpt_>200 && dilpt_<275;

      bool SR1jdn	= htDown_>300 && htDown_<600 && pfmetDown_>275;
      bool SR2jdn	= htDown_>600 && pfmetDown_>275;
      bool SR3jdn	= htDown_>600 && pfmetDown_>200 && pfmetDown_<275;

      // bool SR1ptlljdn	= htDown_>300 && htDown_<600 && dilpt_>275;
      // bool SR2ptlljdn	= htDown_>600 && dilpt_>275;
      // bool SR3ptlljdn	= htDown_>600 && dilpt_>200 && dilpt_<275;

      //----------------------------
      // nominal
      //----------------------------

      if( presel ){
	if( SR1 && SF ) h[0]->Fill( m0_ , m12_ , weight );
	if( SR1 && OF ) h[1]->Fill( m0_ , m12_ , weight );
	if( SR2 && SF ) h[2]->Fill( m0_ , m12_ , weight );
	if( SR2 && OF ) h[3]->Fill( m0_ , m12_ , weight );
	if( SR3 && SF ) h[4]->Fill( m0_ , m12_ , weight );
	if( SR3 && OF ) h[5]->Fill( m0_ , m12_ , weight );
	if( sig       ) hall->Fill( m0_ , m12_ , weight );
      }

      // ch->Draw("m12:m0>>h_0"        , (presel + SR1 + SF) * weight );
      // ch->Draw("m12:m0>>h_1"        , (presel + SR1 + OF) * weight );
      // ch->Draw("m12:m0>>h_2"        , (presel + SR2 + SF) * weight );
      // ch->Draw("m12:m0>>h_3"        , (presel + SR2 + OF) * weight );
      // ch->Draw("m12:m0>>h_4"        , (presel + SR3 + SF) * weight );
      // ch->Draw("m12:m0>>h_5"        , (presel + SR3 + OF) * weight );
      // ch->Draw("m12:m0>>hall"       , (presel + sig     ) * weight );

      //----------------------------
      // ptll control region
      //----------------------------

      if( preselptll ){
	if( SR1ptll && SF ) hptll[0]->Fill( m0_ , m12_ , weight );
	if( SR1ptll && OF ) hptll[1]->Fill( m0_ , m12_ , weight );
	if( SR2ptll && SF ) hptll[2]->Fill( m0_ , m12_ , weight );
	if( SR2ptll && OF ) hptll[3]->Fill( m0_ , m12_ , weight );
	if( SR3ptll && SF ) hptll[4]->Fill( m0_ , m12_ , weight );
	if( SR3ptll && OF ) hptll[5]->Fill( m0_ , m12_ , weight );
      }

      // ch->Draw("m12:m0>>hptll_0"    , (preselptll + SR1ptll) * weight );
      // ch->Draw("m12:m0>>hptll_1"    , (preselptll + SR1ptll) * weight );
      // ch->Draw("m12:m0>>hptll_2"    , (preselptll + SR2ptll) * weight );
      // ch->Draw("m12:m0>>hptll_3"    , (preselptll + SR2ptll) * weight );
      // ch->Draw("m12:m0>>hptll_4"    , (preselptll + SR3ptll) * weight );
      // ch->Draw("m12:m0>>hptll_5"    , (preselptll + SR3ptll) * weight );

      //----------------------------
      // JES up
      //----------------------------

      if( preseljup ){
	if( SR1jup && SF ) hjup[0]->Fill( m0_ , m12_ , weight );
	if( SR1jup && OF ) hjup[1]->Fill( m0_ , m12_ , weight );
	if( SR2jup && SF ) hjup[2]->Fill( m0_ , m12_ , weight );
	if( SR2jup && OF ) hjup[3]->Fill( m0_ , m12_ , weight );
	if( SR3jup && SF ) hjup[4]->Fill( m0_ , m12_ , weight );
	if( SR3jup && OF ) hjup[5]->Fill( m0_ , m12_ , weight );
      }

      //cout << "Filling JES up histos" << endl;
      // ch->Draw("m12:m0>>hjup_0"     , (preseljup + SR1jup + SF) * weight );
      // ch->Draw("m12:m0>>hjup_1"     , (preseljup + SR1jup + OF) * weight );
      // ch->Draw("m12:m0>>hjup_2"     , (preseljup + SR2jup + SF) * weight );
      // ch->Draw("m12:m0>>hjup_3"     , (preseljup + SR2jup + OF) * weight );
      // ch->Draw("m12:m0>>hjup_4"     , (preseljup + SR3jup + SF) * weight );
      // ch->Draw("m12:m0>>hjup_5"     , (preseljup + SR3jup + OF) * weight );

      //----------------------------
      // JES down
      //----------------------------

      if( preseljdn ){
	if( SR1jdn && SF ) hjdn[0]->Fill( m0_ , m12_ , weight );
	if( SR1jdn && OF ) hjdn[1]->Fill( m0_ , m12_ , weight );
	if( SR2jdn && SF ) hjdn[2]->Fill( m0_ , m12_ , weight );
	if( SR2jdn && OF ) hjdn[3]->Fill( m0_ , m12_ , weight );
	if( SR3jdn && SF ) hjdn[4]->Fill( m0_ , m12_ , weight );
	if( SR3jdn && OF ) hjdn[5]->Fill( m0_ , m12_ , weight );
	if( sigdn        ) hjdnall->Fill( m0_ , m12_ , weight );
      }
  
      //JES down
      // cout << "Filling JES dn histos" << endl;
      // ch->Draw("m12:m0>>hjdn_0"     , (preseljdn + SR1jdn + SF) * weight );
      // ch->Draw("m12:m0>>hjdn_1"     , (preseljdn + SR1jdn + OF) * weight );
      // ch->Draw("m12:m0>>hjdn_2"     , (preseljdn + SR2jdn + SF) * weight );
      // ch->Draw("m12:m0>>hjdn_3"     , (preseljdn + SR2jdn + OF) * weight );
      // ch->Draw("m12:m0>>hjdn_4"     , (preseljdn + SR3jdn + SF) * weight );
      // ch->Draw("m12:m0>>hjdn_5"     , (preseljdn + SR3jdn + OF) * weight );
      // ch->Draw("m12:m0>>hjdnall"    , (preseljdn + sigdn      ) * weight );

      
    }// end loop over events
  } // end loop over files
  

  ofstream* doScript = new ofstream();
  doScript->open(Form("cards/%s/doLimits.sh",outversion));

  ofstream* doScript_CLs = new ofstream();
  doScript_CLs->open(Form("cards/%s/doLimits_CLs.sh",outversion));

  ofstream* doScript_CLs1 = new ofstream();
  doScript_CLs1->open(Form("cards/%s/doLimits_CLs1.sh",outversion));

  ofstream* doScript_CLs2 = new ofstream();
  doScript_CLs2->open(Form("cards/%s/doLimits_CLs2.sh",outversion));

  ofstream* doScript_CLs3 = new ofstream();
  doScript_CLs3->open(Form("cards/%s/doLimits_CLs3.sh",outversion));

  ofstream* doScript_CLs4 = new ofstream();
  doScript_CLs4->open(Form("cards/%s/doLimits_CLs4.sh",outversion));

  ofstream* filelist = new ofstream();
  filelist->open(Form("cards/%s/file_list.txt",outversion));

  ofstream* filelist_CLs = new ofstream();
  filelist_CLs->open(Form("cards/%s/file_list_CLs.txt",outversion));


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
      
      //int m0  = hall->GetXaxis()->GetBinCenter(m0bin);
      //int m12 = hall->GetXaxis()->GetBinCenter(m12bin);

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

      *filelist << Form("cards/%s/CMSSM_%i_%i.txt_Bayesian_bysObsLimit.root",outversion,m0bin,m12bin)        << endl;

      *filelist_CLs << Form("cards/%s/CMSSM_%i_%i_m2lnQ2.root",outversion,m0bin,m12bin)        << endl;

      printCard( Form("CMSSM_%i_%i",m0bin,m12bin) , sigtot ,  outversion );
      
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


      TFile *f = TFile::Open( Form("rootfiles/%s/CMSSM_%i_%i.root",outversion,m0bin,m12bin) , "RECREATE");
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

  cout << "Looped over " << nEventsTotal << " events" << endl;
  
  TFile* f = TFile::Open("histos.root","RECREATE");
  f->cd();
  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]        ->Write();
    hjup[i]     ->Write();
    hjdn[i]     ->Write();
    hkup[i]     ->Write();
    hkdn[i]     ->Write();
    hptll[i]    ->Write();
    //hptlljup[i] ->Write();
    //hptlljdn[i] ->Write();
  }
  f->Close();

} // end ScanChain


void babylooper::setBranches (TTree* tree){

  tree->SetBranchAddress("run"		,       &run_		);
  tree->SetBranchAddress("lumi"		,       &lumi_          );
  tree->SetBranchAddress("event"	,       &event_         );
  tree->SetBranchAddress("njets"	,       &njets_         );
  tree->SetBranchAddress("njetsUp"	,       &njetsUp_       );
  tree->SetBranchAddress("njetsDown"	,       &njetsDown_     );
  tree->SetBranchAddress("passz"	,       &passz_         );
  tree->SetBranchAddress("m0"	        ,       &m0_            );
  tree->SetBranchAddress("m12"	        ,       &m12_           );
  tree->SetBranchAddress("leptype"	,       &leptype_       );
  tree->SetBranchAddress("dilpt"	,       &dilpt_         );
  tree->SetBranchAddress("pfmet"	,       &pfmet_         );
  tree->SetBranchAddress("pfmetUp"	,       &pfmetUp_       );
  tree->SetBranchAddress("pfmetDown"	,       &pfmetDown_     );
  tree->SetBranchAddress("ht"		,       &ht_		);
  tree->SetBranchAddress("htUp"		,       &htUp_		);
  tree->SetBranchAddress("htDown"	,       &htDown_	);
  tree->SetBranchAddress("weight"	,       &weight_	);
  tree->SetBranchAddress("ndavtxweight"	,       &ndavtxweight_	);
  tree->SetBranchAddress("trgeff"	,       &trgeff_	);
  tree->SetBranchAddress("lepscale"	,       &lepscale_	);


}
