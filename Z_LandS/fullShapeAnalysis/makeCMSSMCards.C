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

void printCard( char* name , float sigtot , char* version , bool do3jets ){

  ofstream* ofile = new ofstream();

  ofile->open(Form("cards/%s/%s.txt",version,name));

  *ofile <<      "imax 1 number of channels"                                                            << endl;
  *ofile <<      "jmax 1 number of background"                                                          << endl;
  *ofile <<      "kmax * number of nuisance parameters"                                                 << endl;
  if( !do3jets ){
    *ofile <<      "Observation 290                                                           "           << endl;
  }else{
    *ofile <<      "Observation 137                                                           "           << endl;
  }
  *ofile << Form("shapes      *   * ../../rootfiles/%s/%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC" , version , name) << endl;
  *ofile << Form("shapes data_obs * ../../rootfiles/%s/%s.root  histo_Data" , version , name )                               << endl;
  *ofile <<      "bin                                  1       1"                                       << endl;
  *ofile << Form("process                      %s     bkg" , name )                                     << endl;
  *ofile <<      "process                              0       1"                                       << endl;
  if( !do3jets ){
    *ofile << Form("rate                              %.1f   295" , sigtot)                             << endl;
  }else{
    *ofile << Form("rate                              %.1f   129" , sigtot)                             << endl;
  }
  *ofile <<      "lumi                       lnN   1.060       -"                                       << endl;
  *ofile <<      "eff_leptons                lnN   1.050       -"                                       << endl;
  *ofile <<      "JES_shape                shape     1.0       -"                                       << endl;
  *ofile <<      "err                      shape       -     1.0"                                       << endl;
  
  ofile->close();

}


void makeCMSSMCards(){

  //---------------------------------------
  // load TChain
  //---------------------------------------
  
  TChain *ch = new TChain("T1");
  ch->Add("output/V00-02-05/T5zzgmsb_hadoop_baby.root");
  char* version = "V00-01-07";

  bool do3jets = false;

  //---------------------------------------
  // selection
  //---------------------------------------

  TCut weight   ("4.7 * davtxweight * 0.95 * (1000./105000.)");

  TCut presel   ("dilmass>81 && dilmass<101 && njets>=2     && leptype<2");
  TCut preseljup("dilmass>81 && dilmass<101 && njetsup>=2   && leptype<2");
  TCut preseljdn("dilmass>81 && dilmass<101 && njetsdn>=2   && leptype<2");

  if( do3jets ){
    presel    = TCut("dilmass>81 && dilmass<101 && njets>=3     && leptype<2");
    preseljup = TCut("dilmass>81 && dilmass<101 && njetsup>=3   && leptype<2");
    preseljdn = TCut("dilmass>81 && dilmass<101 && njetsdn>=3   && leptype<2");
  }

  TCut met100_200("pfmet>100 && pfmet<200");
  TCut met200_300("pfmet>200 && pfmet<300");
  TCut met300    ("pfmet>300");
  TCut met30     ("pfmet>30");

  TCut met100_200up("pfmetup>100 && pfmetup<200");
  TCut met200_300up("pfmetup>200 && pfmetup<300");
  TCut met300up    ("pfmetup>300");

  TCut met100_200dn("pfmetdn>100 && pfmetdn<200");
  TCut met200_300dn("pfmetdn>200 && pfmetdn<300");
  TCut met300dn    ("pfmetdn>300");
  TCut met30dn     ("pfmetdn>30");

  //---------------------------------------
  // preselection and signal region yields
  //---------------------------------------

  const unsigned int nbins = 3;

  TH2F* h[nbins];
  TH2F* hjup[nbins];
  TH2F* hjdn[nbins];

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]        = new TH2F( Form("h_%i",i)        , Form("h_%i",i)           , 48,0,1200,48,0,1200);
    hjup[i]     = new TH2F( Form("hjup_%i",i)     , Form("hjup_%i",i)        , 48,0,1200,48,0,1200);
    hjdn[i]     = new TH2F( Form("hjdn_%i",i)     , Form("hjdn_%i",i)        , 48,0,1200,48,0,1200);
    
    h[i]    ->Sumw2();
    hjup[i] ->Sumw2();
    hjdn[i] ->Sumw2();
  }
  
  TH2F* hall    = new TH2F( "hall"    , "hall"    , 48,0,1200,48,0,1200);
  TH2F* hjdnall = new TH2F( "hjdnall" , "hjdnall" , 48,0,1200,48,0,1200);
  
  hall->Sumw2();
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  cout << "Filling nominal histos" << endl;
  ch->Draw("ml:mg>>h_0"         , (presel+met100_200) * weight);
  ch->Draw("ml:mg>>h_1"         , (presel+met200_300) * weight);
  ch->Draw("ml:mg>>h_2"         , (presel+met300    ) * weight);
  ch->Draw("ml:mg>>hall"        , (presel+met30     ) * weight);

  cout << "Filling JES up histos" << endl;
  ch->Draw("ml:mg>>hjup_0"      , (preseljup+met100_200up) * weight);
  ch->Draw("ml:mg>>hjup_1"      , (preseljup+met200_300up) * weight);
  ch->Draw("ml:mg>>hjup_2"      , (preseljup+met300up    ) * weight);

  cout << "Filling JES down histos" << endl;
  ch->Draw("ml:mg>>hjdn_0"      , (preseljdn+met100_200dn) * weight);
  ch->Draw("ml:mg>>hjdn_1"      , (preseljdn+met200_300dn) * weight);
  ch->Draw("ml:mg>>hjdn_2"      , (preseljdn+met300dn    ) * weight);
  ch->Draw("ml:mg>>hjdnall"     , (preseljdn+met30dn     ) * weight);

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
  // loop over CMSSM points
  //------------------------------------------

  int counter = 0;

  for( int mgbin = 1 ; mgbin <= hall->GetXaxis()->GetNbins() ; mgbin++ ){
    for( int mlbin = 1 ; mlbin <= hall->GetYaxis()->GetNbins() ; mlbin++ ){

      //if( !( mgbin == 40 && mlbin == 5 ) && !(mgbin == 20 && mlbin == 10) ) continue;

      int mg  = hall->GetXaxis()->GetBinCenter(mgbin)-12.5;
      int ml  = hall->GetXaxis()->GetBinCenter(mlbin)-12.5;

      //cout << "mg " << mg << " ml " << ml << endl;

      if( hjdnall->GetBinContent(mgbin,mlbin) < 1e-10 ) continue;

      float nom = hall->GetBinContent(mgbin,mlbin);

      //---------------------------------------
      // make root file
      //---------------------------------------

      TH1F* histo_SMS               = new TH1F( Form("histo_SMS_%i_%i"              ,mgbin,mlbin) , Form("histo_SMS_%i_%i"              ,mgbin,mlbin) , nbins,0,nbins);
      TH1F* histo_SMS_JES_shapeUp   = new TH1F( Form("histo_SMS_%i_%i_JES_shapeUp"  ,mgbin,mlbin) , Form("histo_SMS_%i_%i_JES_shapeUp"  ,mgbin,mlbin) , nbins,0,nbins);
      TH1F* histo_SMS_JES_shapeDown = new TH1F( Form("histo_SMS_%i_%i_JES_shapeDown",mgbin,mlbin) , Form("histo_SMS_%i_%i_JES_shapeDown",mgbin,mlbin) , nbins,0,nbins);

      float sigtot    = 0;
      float sigtotjdn = 0;

      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
	float yieldnom = h[ibin]->GetBinContent(mgbin,mlbin);
	float yieldjup = hjup[ibin]->GetBinContent(mgbin,mlbin);
	float yieldjdn = hjdn[ibin]->GetBinContent(mgbin,mlbin);

	sigtot    += yieldnom;
	sigtotjdn += yieldjdn;

	histo_SMS->SetBinContent              ( ibin + 1 , yieldnom );
	histo_SMS_JES_shapeUp->SetBinContent  ( ibin + 1 , yieldjup );
	histo_SMS_JES_shapeDown->SetBinContent( ibin + 1 , yieldjdn );
      }

      if( sigtotjdn < 1e-10 ) continue;
      //if( sigtot    < 2     ) continue;

      counter++;

      *doScript << Form("../../../../test/lands.exe -M Bayesian -d SMS_%i_%i.txt",mgbin,mlbin)         << endl;

      *doScript_CLs << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      *filelist << Form("cards/%s/SMS_%i_%i.txt_Bayesian_bysObsLimit.root",version,mgbin,mlbin)        << endl;

      *filelist_CLs << Form("cards/%s/SMS_%i_%i_m2lnQ2.root",version,mgbin,mlbin)        << endl;

      if( counter%4 == 0 ) *doScript_CLs1 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 1 ) *doScript_CLs2 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 2 ) *doScript_CLs3 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 3 ) *doScript_CLs4 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      printCard( Form("SMS_%i_%i",mgbin,mlbin) , sigtot , version , do3jets);
      
      //signal regions                          met100    met200    met300
      int     data_yield[nbins]           = {     276   ,    14   ,    0 };
      float   bkg_yield[nbins]            = {     276   ,  15.7   , 3.09 };
      float   bkg_err[nbins]              = {      27   ,  2.60   ,  1.0 };

      if( do3jets ){

	data_yield[0]   =  129;
	data_yield[1]   =    8;
	data_yield[2]   =    0;

	bkg_yield[0]    =  119.;
	bkg_yield[1]    =   8.7;
	bkg_yield[2]    =   1.8;

	bkg_err[0]      =   12.;
	bkg_err[1]      =   1.7;
	bkg_err[2]      =   0.6;

      }

      TH1F* histo_Data = new TH1F("histo_Data","histo_Data",nbins,0,nbins);

      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
	histo_Data->SetBinContent( ibin+1 , data_yield[ibin] );
      }

      TH1F* histo_bkg               = new TH1F("histo_bkg","histo_bkg"                ,nbins,0,nbins);
      TH1F* histo_bkg_errUp         = new TH1F("histo_bkg_errUp"  ,"histo_bkg_errUp"  ,nbins,0,nbins);
      TH1F* histo_bkg_errDown       = new TH1F("histo_bkg_errDown","histo_bkg_errDown",nbins,0,nbins);
      
      for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
	histo_bkg               -> SetBinContent(ibin+1, bkg_yield[ibin]);
	histo_bkg_errUp         -> SetBinContent(ibin+1, bkg_yield[ibin] + bkg_err[ibin] );
	histo_bkg_errDown       -> SetBinContent(ibin+1, TMath::Max(bkg_yield[ibin] - bkg_err[ibin],(float)0.0) );
      }


      TFile *f = TFile::Open( Form("rootfiles/%s/SMS_%i_%i.root",version,mgbin,mlbin) , "RECREATE");
      f->cd();
      histo_Data->Write();
      histo_bkg->Write();
      histo_bkg_errUp->Write();
      histo_bkg_errDown->Write();
      histo_SMS->Write();
      histo_SMS_JES_shapeUp->Write();
      histo_SMS_JES_shapeDown->Write();
      f->Close();

      delete histo_Data;
      delete histo_bkg;
      delete histo_bkg_errUp;
      delete histo_bkg_errDown;
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
  filelist->close();
  filelist_CLs->close();

}
