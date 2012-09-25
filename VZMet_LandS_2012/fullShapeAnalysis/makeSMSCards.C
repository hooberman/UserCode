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

void printCard( char* name , float sigtot , float Ztot , float OFtot , float VZtot , int datatot , char* version ){

  ofstream* ofile = new ofstream();

  ofile->open(Form("cards/%s/%s.txt",version,name));

  *ofile <<      "imax 1 number of channels"                                                                  << endl;
  *ofile <<      "jmax 3 number of background"                                                                << endl;
  *ofile <<      "kmax * number of nuisance parameters"                                                       << endl;
  *ofile << Form("Observation %i                                                           ",datatot)         << endl;
  //*ofile << Form("shapes      *   * ../../rootfiles/%s/%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC" , version , name) << endl;
  //*ofile << Form("shapes data_obs * ../../rootfiles/%s/%s.root  histo_Data" , version , name )                << endl;
  *ofile << Form("shapes      *   * %s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC" ,  name) << endl;
  *ofile << Form("shapes data_obs * %s.root  histo_Data" ,  name )                << endl;
  *ofile <<      "bin                                  1        1      1      1"                              << endl;
  *ofile << Form("process                        %s     Zbkg  OFbkg  VZbkg" , name )                          << endl;
  *ofile <<      "process                              0        1      2      3"                              << endl;
  *ofile << Form("rate                              %.1f    %.1f    %.1f   %.1f" , sigtot,Ztot,OFtot,VZtot)   << endl;
  *ofile <<      "lumi                       lnN   1.022       -       -      -"                              << endl;
  *ofile <<      "eff_leptons                lnN   1.050       -       -      -"                              << endl;
  *ofile <<      "btagerr                    lnN   1.040       -       -      -"                              << endl;
  *ofile <<      "JES_shape                shape     1.0       -       -      -"                              << endl;
  *ofile <<      "errZ                     shape       -     1.0       -      -"                              << endl;
  *ofile <<      "errOF                    shape       -       -     1.0      -"                              << endl;
//*ofile <<      "errVZ                    shape       -       -       -    1.0"                              << endl;
  *ofile <<      "wz                       shape       -       -       -    1.0"                              << endl;
  
  ofile->close();

}


void makeSMSCards(){

  //---------------------------------------
  // load TChain
  //---------------------------------------
  
  TChain *ch = new TChain("T1");
  ch->Add("output/V00-01-04/wzsms_baby_pt2020_oldIso.root ");
  char* version = (char*) "V00-00-00";

  //---------------------------------------
  // selection
  //---------------------------------------

  TCut weight   ("9.2 * trgeff * vtxweight * (1./100000.)");

  TCut presel   ("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvm==0 && mjj>70   && mjj<110   && nlep==2 && njets>=2     && leptype<2");
  TCut preseljup("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvm==0 && mjjup>70 && mjjup<110 && nlep==2 && njetsup>=2   && leptype<2");
  TCut preseljdn("lep2.pt()>20.0 && dilmass>81 && dilmass<101 && nbcsvm==0 && mjjdn>70 && mjjdn<110 && nlep==2 && njetsdn>=2   && leptype<2");

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

  for( unsigned int i = 0 ; i < nbins ; ++i ){
    h[i]        = new TH2F( Form("h_%i",i)        , Form("h_%i",i)           , 24,0,600,24,0,600);
    hjup[i]     = new TH2F( Form("hjup_%i",i)     , Form("hjup_%i",i)        , 24,0,600,24,0,600);
    hjdn[i]     = new TH2F( Form("hjdn_%i",i)     , Form("hjdn_%i",i)        , 24,0,600,24,0,600);
    
    h[i]    ->Sumw2();
    hjup[i] ->Sumw2();
    hjdn[i] ->Sumw2();
  }
  
  TH2F* hall    = new TH2F( "hall"    , "hall"    , 24,0,600,24,0,600);
  TH2F* hjdnall = new TH2F( "hjdnall" , "hjdnall" , 24,0,600,24,0,600);
  
  hall->Sumw2();
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  cout << "Filling histos..." << endl;
  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++ ){
    ch->Draw(Form("ml:mg>>h_%i"    , ibin) , (presel    + sigcut[ibin]   ) * weight);
    ch->Draw(Form("ml:mg>>hjup_%i" , ibin) , (preseljup + sigcutup[ibin] ) * weight);
    ch->Draw(Form("ml:mg>>hjdn_%i" , ibin) , (preseljdn + sigcutdn[ibin] ) * weight);
  }

  ch->Draw("ml:mg>>hall"        , (presel    + sigall     ) * weight);
  ch->Draw("ml:mg>>hjdnall"     , (preseljdn + sigalldn   ) * weight);

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

  TH1F* histo_Zbkg               = new TH1F("histo_Zbkg"           ,"histo_Zbkg"           ,nbins,0,nbins);
  TH1F* histo_Zbkg_errUp         = new TH1F("histo_Zbkg_errZUp"    ,"histo_Zbkg_errZUp"    ,nbins,0,nbins);
  TH1F* histo_Zbkg_errDown       = new TH1F("histo_Zbkg_errZDown"  ,"histo_Zbkg_errZDown"  ,nbins,0,nbins);

  TH1F* histo_OFbkg              = new TH1F("histo_OFbkg"          ,"histo_OFbkg"          ,nbins,0,nbins);
  TH1F* histo_OFbkg_errUp        = new TH1F("histo_OFbkg_errOFUp"  ,"histo_OFbkg_errOFUp"  ,nbins,0,nbins);
  TH1F* histo_OFbkg_errDown      = new TH1F("histo_OFbkg_errOFDown","histo_OFbkg_errOFDown",nbins,0,nbins);

  TH1F* histo_WZbkg              = new TH1F("histo_WZbkg"          ,"histo_WZbkg"          ,nbins,0,nbins);
  TH1F* histo_WZbkg_errUp        = new TH1F("histo_WZbkg_wzUp"     ,"histo_WZbkg_wzUp"     ,nbins,0,nbins);
  TH1F* histo_WZbkg_errDown      = new TH1F("histo_WZbkg_wzDown"   ,"histo_WZbkg_wzDown"   ,nbins,0,nbins);
      
  TH1F* histo_ZZbkg              = new TH1F("histo_ZZbkg"          ,"histo_ZZbkg"          ,nbins,0,nbins);
  TH1F* histo_ZZbkg_errUp        = new TH1F("histo_ZZbkg_wzUp"     ,"histo_ZZbkg_wzUp"     ,nbins,0,nbins);
  TH1F* histo_ZZbkg_errDown      = new TH1F("histo_ZZbkg_wzDown"   ,"histo_ZZbkg_wzDown"   ,nbins,0,nbins);
      
  TH1F* histo_rarebkg            = new TH1F("histo_rarebkg"        ,"histo_rarebkg"        ,nbins,0,nbins);
  TH1F* histo_rarebkg_errUp      = new TH1F("histo_rarebkg_wzUp"   ,"histo_rarebkg_wzUp"   ,nbins,0,nbins);
  TH1F* histo_rarebkg_errDown    = new TH1F("histo_rarebkg_wzDown" ,"histo_rarebkg_wzDown" ,nbins,0,nbins);
      
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

      int mg  = hall->GetXaxis()->GetBinCenter(mgbin)-12.5;
      int ml  = hall->GetXaxis()->GetBinCenter(mlbin)-12.5;

      cout << "mg " << mg << " ml " << ml << endl;

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

      *filelist_CLs << Form("cards/%s/SMS_%i_%i_m2lnQ2.root",version,mgbin,mlbin)        << endl;

      if( counter%4 == 0 ) *doScript_CLs1 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 1 ) *doScript_CLs2 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 2 ) *doScript_CLs3 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( counter%4 == 3 ) *doScript_CLs4 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( mlbin==1 ) *doScript_CLs_mL0 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( mlbin==3 ) *doScript_CLs_mL50 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      if( mlbin==5 ) *doScript_CLs_mL100 << Form("../../../../test/lands.exe -d SMS_%i_%i.txt  -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 3000 --nToysForCLb 1500 --seed 1234 -n SMS_%i_%i -rMin 0 -rMax 100",mgbin,mlbin,mgbin,mlbin) << endl;

      printCard( Form("SMS_%i_%i",mgbin,mlbin) , sigtot , Zbkg_tot , OFbkg_tot , VZbkg_tot , data_tot , version );


      TFile *f = TFile::Open( Form("rootfiles/%s/SMS_%i_%i.root",version,mgbin,mlbin) , "RECREATE");
      f->cd();
      histo_Data->Write();
      histo_Zbkg->Write();
      histo_Zbkg_errUp->Write();
      histo_Zbkg_errDown->Write();
      histo_OFbkg->Write();
      histo_OFbkg_errUp->Write();
      histo_OFbkg_errDown->Write();
      histo_VZbkg->Write();
      histo_VZbkg_errUp->Write();
      histo_VZbkg_errDown->Write();
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
