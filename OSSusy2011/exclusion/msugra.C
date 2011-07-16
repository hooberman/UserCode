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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"

//#include "cl95cms_landsberg.c"

using namespace std;


//-------------------------------------------
//parameters for CL95 function
//-------------------------------------------

const Double_t mylumi        = 0.976;
const bool     doCorrection  = false;
const bool     addHighm0     = false;

/*
const Double_t ilum            = 976.0; // lumi
const Double_t slum            = 0.;    // lumi uncertainty (=0 b/c uncertainty is included in sig acceptance)
const Double_t eff             = 1.;    // sig efficiency
const Double_t bck             = 5.1;   // expected background
const Double_t sbck            = 1.7;   // background error
const int      n               = 4;     // observed yield
const int      nuissanceModel  = 1;     // nuissance model (0 - Gaussian, 1 - lognormal, 2 - gamma)
*/

//-------------------------------------------
// uncertainties
//-------------------------------------------

const float lumierr = 0.06;  // 11% lumi error
const float leperr  = 0.04;  // lepton ID/iso + trigger uncertainty
const float pdferr  = 0.20;  // PDF uncertainty

//----------------------------------------------
// ABCD yields (for signal contamination study)
//----------------------------------------------

const bool  doABCD  = false; //turned off by default
const float yieldA  =  12.;
const float yieldB  =  37.;
const float yieldC  =   4.;

//---------------------------------------------------------------------------
// calculate the expected UL? this is time-consuming so turn off if not needed
//---------------------------------------------------------------------------

const bool calculateExpectedUL = false;

//-----------------------------------------------
//rebin the TH2 yield histos (NOT RECOMMENDED)
//-----------------------------------------------

const int   rebin   = 1;     // rebin yield histos


float getUpperLimit( float seff ){

  float ul = 1e10;

  //-----------
  // CLs
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 5.35;
  if(seff > 0.1 && seff < 0.2) ul = 5.83;
  if(seff > 0.2 && seff < 0.3) ul = 6.09;
  if(seff > 0.3 && seff < 0.4) ul = 6.57;
  if(seff > 0.4 && seff < 0.5) ul = 7.04;
  if(seff > 0.5 && seff < 0.6) ul = 7.57;
  if(seff > 0.6 && seff < 0.7) ul = 8.23;
  if(seff > 0.7 && seff < 0.8) ul = 8.56;
  if(seff > 0.8 && seff < 0.9) ul = 9.17;
  
  //-----------
  // bayesian
  //----------
  /*  
  if(seff > 0.00 && seff < 0.05) ul = 5.79;
  if(seff > 0.05 && seff < 0.10) ul = 5.87;
  if(seff > 0.10 && seff < 0.15) ul = 6.01;
  if(seff > 0.15 && seff < 0.20) ul = 6.19;
  if(seff > 0.20 && seff < 0.25) ul = 6.42;
  if(seff > 0.25 && seff < 0.30) ul = 6.68;
  if(seff > 0.30 && seff < 0.35) ul = 6.99;
  if(seff > 0.35 && seff < 0.40) ul = 7.33;
  if(seff > 0.40 && seff < 0.45) ul = 7.71;
  if(seff > 0.45 && seff < 0.50) ul = 8.14;
  if(seff > 0.50 && seff < 0.55) ul = 8.61;
  if(seff > 0.55 && seff < 0.60) ul = 9.11;
  if(seff > 0.60 && seff < 0.65) ul = 9.65;
  if(seff > 0.65 && seff < 0.70) ul = 10.24;
  if(seff > 0.70 && seff < 0.75) ul = 10.88;
  if(seff > 0.75 && seff < 0.80) ul = 11.58;
  if(seff > 0.80 && seff < 0.85) ul = 12.35;
  if(seff > 0.85 && seff < 0.90) ul = 13.16;
  */

  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;
  return ul;
}

float getExpectedUpperLimit( float seff ){

  float ul = 1e10;

  //-----------
  // expected
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 6.09;
  if(seff > 0.1 && seff < 0.2) ul = 6.38;
  if(seff > 0.2 && seff < 0.3) ul = 6.77;
  if(seff > 0.3 && seff < 0.4) ul = 7.60;
  if(seff > 0.4 && seff < 0.5) ul = 7.93;
  if(seff > 0.5 && seff < 0.6) ul = 8.37;
  if(seff > 0.6 && seff < 0.7) ul = 9.07;
  if(seff > 0.7 && seff < 0.8) ul = 9.64;
  if(seff > 0.8 && seff < 0.9) ul = 10.11;

  //-----------
  // bayesian
  //----------
  /*  
  if(seff > 0.00 && seff < 0.05) ul = 6.99;
  if(seff > 0.05 && seff < 0.10) ul = 7.11;
  if(seff > 0.10 && seff < 0.15) ul = 7.29;
  if(seff > 0.15 && seff < 0.20) ul = 7.53;
  if(seff > 0.20 && seff < 0.25) ul = 7.82;
  if(seff > 0.25 && seff < 0.30) ul = 8.17;
  if(seff > 0.30 && seff < 0.35) ul = 8.57;
  if(seff > 0.35 && seff < 0.40) ul = 9.02;
  if(seff > 0.40 && seff < 0.45) ul = 9.52;
  if(seff > 0.45 && seff < 0.50) ul = 10.07;
  if(seff > 0.50 && seff < 0.55) ul = 10.67;
  if(seff > 0.55 && seff < 0.60) ul = 11.33;
  if(seff > 0.60 && seff < 0.65) ul = 12.03;
  if(seff > 0.65 && seff < 0.70) ul = 12.79;
  if(seff > 0.70 && seff < 0.75) ul = 13.62;
  if(seff > 0.75 && seff < 0.80) ul = 14.54;
  if(seff > 0.80 && seff < 0.85) ul = 15.54;
  if(seff > 0.85 && seff < 0.90) ul = 16.61;
  */

  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;
}

float getExpectedP1UpperLimit( float seff ){

  float ul = 1e10;

  //-----------
  // expected
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 8.70;
  if(seff > 0.1 && seff < 0.2) ul = 9.17;
  if(seff > 0.2 && seff < 0.3) ul = 9.81;
  if(seff > 0.3 && seff < 0.4) ul = 11.30;
  if(seff > 0.4 && seff < 0.5) ul = 12.70;
  if(seff > 0.5 && seff < 0.6) ul = 14.02;
  if(seff > 0.6 && seff < 0.7) ul = 15.08;
  if(seff > 0.7 && seff < 0.8) ul = 16.32;
  if(seff > 0.8 && seff < 0.9) ul = 17.26;

  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;
}

float getExpectedM1UpperLimit( float seff ){

  float ul = 1e10;

  //-----------
  // expected
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 4.23;
  if(seff > 0.1 && seff < 0.2) ul = 4.38;
  if(seff > 0.2 && seff < 0.3) ul = 4.93;
  if(seff > 0.3 && seff < 0.4) ul = 5.29;
  if(seff > 0.4 && seff < 0.5) ul = 5.54;
  if(seff > 0.5 && seff < 0.6) ul = 6.02;
  if(seff > 0.6 && seff < 0.7) ul = 6.07;
  if(seff > 0.7 && seff < 0.8) ul = 6.64;
  if(seff > 0.8 && seff < 0.9) ul = 6.84;
  
  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;

}

TH1F* getCurve               ( TH2I *hist , char* name );
TGraphErrors* getCurve_TGraph( TH2I *hist , char* name );

void msugra( char* filename ){
  
  TFile *corfile = new TFile();
  TH2F* hscan    = new TH2F();
  
  if( doCorrection ){
    
    corfile = TFile::Open("mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0.root");
    hscan = (TH2F*) corfile->Get("hscan");

    if( hscan == 0 ){
      cout << "Can't find TH2 hscan!!!" << endl;
      exit(0);
    }
  }

  TFile *outfile = new TFile("exclusion.root","RECREATE");
  //TFile *outfile = new TFile("exclusion_Fall10_pfmet_pfjets.root","RECREATE");

  //-----------------------------------
  // Here we load the yield histograms
  //-----------------------------------
  TFile *f = TFile::Open( filename );

  /*
  //TH2F* hyield     = (TH2F*) f->Get("msugra",prefix));
  TH2F* hyield_k   = (TH2F*) f->Get("msugra_highht");
  TH2F* hyield_kup = (TH2F*) f->Get("msugra_highht_kup");
  TH2F* hyield_kdn = (TH2F*) f->Get("msugra_highht_kdn");
  TH2F* hyield_jup = (TH2F*) f->Get("msugra_highht_jup");
  TH2F* hyield_jdn = (TH2F*) f->Get("msugra_highht_jdn");

  if( addHighm0 ){

    TFile *fm0 = TFile::Open( "../output/V00-01-07/highpt/ossusy_pfjet_pfmet_LMscan_m0_1000.root" );

    //TH2F* hyield     = (TH2F*) f->Get("msugra",prefix));
    TH2F* hyield_k_m0   = (TH2F*) fm0->Get("msugra_highht");
    TH2F* hyield_kup_m0 = (TH2F*) fm0->Get("msugra_highht_kup");
    TH2F* hyield_kdn_m0 = (TH2F*) fm0->Get("msugra_highht_kdn");
    TH2F* hyield_jup_m0 = (TH2F*) fm0->Get("msugra_highht_jup");
    TH2F* hyield_jdn_m0 = (TH2F*) fm0->Get("msugra_highht_jdn");

    float kfactor = 1.5;

    hyield_k_m0->Scale(kfactor);
    hyield_kup_m0->Scale(kfactor);
    hyield_kdn_m0->Scale(kfactor);
    hyield_jup_m0->Scale(kfactor);
    hyield_jdn_m0->Scale(kfactor);

    hyield_k  ->Add( hyield_k_m0   );
    hyield_kup->Add( hyield_kup_m0 );
    hyield_kdn->Add( hyield_kdn_m0 );
    hyield_jup->Add( hyield_jup_m0 );
    hyield_jdn->Add( hyield_jdn_m0 );
  }
  */



  TH2F* hyield_k   = (TH2F*) f->Get("LMscan10_lmgridyield_k");
  TH2F* hyield_kup = (TH2F*) f->Get("LMscan10_lmgridyield_kup");
  TH2F* hyield_kdn = (TH2F*) f->Get("LMscan10_lmgridyield_kdn");
  TH2F* hyield_jup = (TH2F*) f->Get("LMscan10_lmgridyield_jup");
  TH2F* hyield_jdn = (TH2F*) f->Get("LMscan10_lmgridyield_jdn");







  hyield_k->Scale(mylumi);
  hyield_kup->Scale(mylumi);
  hyield_kdn->Scale(mylumi);
  hyield_jup->Scale(mylumi);
  hyield_jdn->Scale(mylumi);

  TH2F* hyield_A; 
  TH2F* hyield_B; 
  TH2F* hyield_C; 
  TH2F* hyield_D; 

  TH2F* hUL_NLO       = (TH2F*) hyield_k->Clone("hUL_NLO");
  TH2F* hUL_NLO_exp   = (TH2F*) hyield_k->Clone("hUL_NLO_exp");
  TH2F* hUL_NLO_expp1 = (TH2F*) hyield_k->Clone("hUL_NLO_expp1");
  TH2F* hUL_NLO_expm1 = (TH2F*) hyield_k->Clone("hUL_NLO_expm1");
  TH2F* hUL_LO        = (TH2F*) hyield_k->Clone("hUL_LO");

  TH1F* htotuncertainty = new TH1F("htotuncertainty","",200,0,2);

  hUL_NLO->Reset();
  hUL_LO->Reset();
  hUL_NLO_exp->Reset();

  const unsigned int nm0bins  = hyield_k->GetXaxis()->GetNbins();
  const unsigned int nm12bins = hyield_k->GetYaxis()->GetNbins();

  //---------------------------------------------------------
  //calculate acceptance error and UL at each scan point
  //---------------------------------------------------------

  float    accerr_NLO[nm0bins][nm12bins];
  Double_t ul_NLO[nm0bins][nm12bins];
  Double_t ul_NLO_SC[nm0bins][nm12bins];
  Double_t ul_NLO_nobkg[nm0bins][nm12bins];
  Double_t ul_NLO_exp[nm0bins][nm12bins];
  Double_t ul_NLO_expp1[nm0bins][nm12bins];
  Double_t ul_NLO_expm1[nm0bins][nm12bins];

  float accerr_LO[nm0bins][nm12bins];
  Double_t ul_LO[nm0bins][nm12bins];

  for( unsigned int m0bin = 1 ; m0bin <= nm0bins ; ++m0bin ){
    for( unsigned int m12bin = 1 ; m12bin <= nm12bins ; ++m12bin ){
      accerr_NLO[m0bin-1][m12bin-1]	= 0.;
      accerr_LO[m0bin-1][m12bin-1]	= 0.;
      ul_NLO[m0bin-1][m12bin-1]		= 9999.;
      ul_NLO_SC[m0bin-1][m12bin-1]	= 9999.;
      ul_NLO_nobkg[m0bin-1][m12bin-1]	= 9999.;
      ul_NLO_exp[m0bin-1][m12bin-1]	= 9999.;
      ul_NLO_expp1[m0bin-1][m12bin-1]   = 9999.;
      ul_NLO_expm1[m0bin-1][m12bin-1]   = 9999.;
      ul_LO[m0bin-1][m12bin-1]		= 9999.;
    }
  }

  //-----------------------------------------------
  // loop over scan points
  //-----------------------------------------------
      
  for( unsigned int m0bin = 1  ; m0bin  <= nm0bins  ; ++m0bin ){
    for( unsigned int m12bin = 1 ; m12bin <= nm12bins ; ++m12bin ){
     
      //---------------------------
      //get yields from TH2's
      //---------------------------
     
      //float yield     = hyield->GetBinContent( m0bin , m12bin );
      float yield_k   = hyield_k->GetBinContent( m0bin , m12bin );
      float yield_kup = hyield_kup->GetBinContent( m0bin , m12bin );
      float yield_kdn = hyield_kdn->GetBinContent( m0bin , m12bin );
      float yield_jup = hyield_jup->GetBinContent( m0bin , m12bin );
      float yield_jdn = hyield_jdn->GetBinContent( m0bin , m12bin );

      int   ngen  = 1;
      float scale = 1.0;

      if( doCorrection ){
	ngen = hscan->GetBinContent(m0bin,m12bin);

	if( ngen != 10000 ){
	  cout << "Skipping point with " << ngen << " entries" << endl;
	  continue;
	}
      }

      //---------------------------
      // find LM6
      //---------------------------

      // if( fabs( hyield_k->GetXaxis()->GetBinCenter(m0bin) - 80 ) < 0.1 && fabs ( hyield_k->GetYaxis()->GetBinCenter(m12bin) - 400 ) < 0.1 ){
      //   cout << "Found LM6" << endl;
      // }
      // else{ continue; }
     
      cout << endl << endl << "------------------------------------------------------------------------" << endl;
      cout << endl << "m0 " << m0bin-1 << " m12 " << m12bin-1 << endl;
      cout << hyield_k->GetXaxis()->GetBinCenter(m0bin) << " " << hyield_k->GetYaxis()->GetBinCenter(m12bin) << endl;

      if( doCorrection ){
	cout << "ngen " << ngen << endl;
	cout << "Rescaling yields by " << scale << endl;
      }

      //------------------------
      // skip bins with 0 yield
      //------------------------
     
      if( fabs( yield_k ) < 1.e-10 ){
	cout << "zero yield, skipping!" << endl;
	continue; 
      }

      //------------------------
      // skip bins m0 > 1000 GeV or m12 > 500 GeV
      //------------------------

      // if( hyield_k->GetXaxis()->GetBinCenter(m0bin) > 1000 ){
      //   cout << "m0 > 1000 GeV, skipping" << endl;
      //   continue;
      // }

      // if( hyield_k->GetYaxis()->GetBinCenter(m12bin) > 500 ){
      //   cout << "m12 > 500 GeV, skipping" << endl;
      //   continue;
      // }

      //-------------------------------------------------
      // this can save a lot of time, turned off here
      //-------------------------------------------------
     
      //      //a point with LO yield > 10 is definitely excluded
      //      if( yield > 10. ){
      //        ul_NLO[m0bin-1][m12bin-1]        = -1;
      //        ul_NLO_SC[m0bin-1][m12bin-1]     = -1;
      //        ul_NLO_nobkg[m0bin-1][m12bin-1]  = -1;
      //        ul_NLO_exp[m0bin-1][m12bin-1]    = -1;
      //        ul_LO[m0bin-1][m12bin-1]         = -1;
      //        hUL_NLO->SetBinContent     ( m0bin , m12bin , -1 );
      //        hUL_NLO_exp->SetBinContent ( m0bin , m12bin , -1 );
      //        hUL_LO->SetBinContent      ( m0bin , m12bin , -1 );
      //        cout << "yield " << yield << " point is excluded, skipping" << endl;
      //        continue;
      //      }

      //      //a point with NLO yield < 3 is definitely not excluded
      //      if( yield_k < 3. ){
      //        ul_NLO[m0bin-1][m12bin-1]        = 100;
      //        ul_NLO_SC[m0bin-1][m12bin-1]     = 100;
      //        ul_NLO_nobkg[m0bin-1][m12bin-1]  = 100;
      //        ul_NLO_exp[m0bin-1][m12bin-1]    = 100;
      //        ul_LO[m0bin-1][m12bin-1]         = 100;
      //        cout << "yield " << yield << " point is NOT excluded, skipping" << endl;
      //        hUL_NLO->SetBinContent     ( m0bin , m12bin , 100 );
      //        hUL_NLO_exp->SetBinContent ( m0bin , m12bin , 100 );
      //        hUL_LO->SetBinContent      ( m0bin , m12bin , 100 );
      //        continue;
      //      }
     

      //--------------------------
      //uncertainty from k-factor
      //--------------------------

      float kerr_up   = fabs( ( yield_kup - yield_k ) / yield_k );
      float kerr_dn   = fabs( ( yield_kdn - yield_k ) / yield_k );
      float kerr      = TMath::Max( kerr_up , kerr_dn );

      //--------------------------     
      //uncertainty from JES
      //--------------------------

      float jerr_up   = fabs( ( yield_jup - yield_k ) / yield_k );
      float jerr_dn   = fabs( ( yield_jdn - yield_k ) / yield_k );
      float jerr      = TMath::Max( jerr_up , jerr_dn );
     
      //-----------------------------------------------------------
      //add up NLO uncertainties (including k-factor uncertainty)
      //-----------------------------------------------------------

      float err2_NLO = 0.;
      err2_NLO += kerr * kerr;                                              //k-factor
      err2_NLO += jerr * jerr;                                              //JES
      err2_NLO += lumierr * lumierr;                                        //lumi
      err2_NLO += leperr * leperr;                                          //lep efficiency
      err2_NLO += pdferr * pdferr;                                          //PDF uncertainty
      accerr_NLO[m0bin-1][m12bin-1] = err2_NLO > 0 ? sqrt( err2_NLO ) : 0.; //total

      if( yield_k > 5 && yield_k < 8 )
	htotuncertainty->Fill( accerr_NLO[m0bin-1][m12bin-1] );

      //ul_NLO[m0bin-1][m12bin-1] = ilum * CL95( ilum, slum, eff, accerr_NLO[m0bin-1][m12bin-1], bck, sbck, n, false, nuissanceModel );
      ul_NLO[m0bin-1][m12bin-1]       = getUpperLimit          ( accerr_NLO[m0bin-1][m12bin-1] );
      ul_NLO_exp[m0bin-1][m12bin-1]   = getExpectedUpperLimit  ( accerr_NLO[m0bin-1][m12bin-1] );
      ul_NLO_expp1[m0bin-1][m12bin-1] = getExpectedP1UpperLimit( accerr_NLO[m0bin-1][m12bin-1] );
      ul_NLO_expm1[m0bin-1][m12bin-1] = getExpectedM1UpperLimit( accerr_NLO[m0bin-1][m12bin-1] );

      //----------------------------------------------------------------
      //add up LO uncertainties (NOT including k-factor uncertainty)
      //----------------------------------------------------------------

      float err2_LO = 0.;
      err2_LO += jerr * jerr;                                               //JES
      err2_LO += lumierr * lumierr;                                         //lumi
      err2_LO += leperr * leperr;                                           //lep efficiency
      err2_LO += pdferr * pdferr;                                           //PDF uncertainty
      accerr_LO[m0bin-1][m12bin-1] = err2_LO > 0 ? sqrt( err2_LO ) : 0.;    //total

      //--------------------------
      //printout errors and UL
      //--------------------------
     
      //cout << "yield            " << yield   << endl;
      cout << "yield * K        " << yield_k << endl;
      cout << "yield * Kup      " << yield_kup << endl;
      cout << "yield * Kdn      " << yield_kdn << endl;
      cout << "yield * JESup    " << yield_jup << endl;
      cout << "yield * JESdn    " << yield_jdn << endl;
      cout << "K error          " << kerr << endl;
      cout << "JES error        " << jerr << endl;
      cout << "total error NLO  " << accerr_NLO[m0bin-1][m12bin-1] << endl;
      cout << "NLO UL           " << ul_NLO[m0bin-1][m12bin-1] << endl;
      cout << "NLO UL sig cont  " << ul_NLO_SC[m0bin-1][m12bin-1] << endl;
      cout << "NLO UL no bkg    " << ul_NLO_nobkg[m0bin-1][m12bin-1] << endl;
      cout << "NLO UL exp       " << ul_NLO_exp[m0bin-1][m12bin-1] << endl;
      cout << "total error LO   " << accerr_LO[m0bin-1][m12bin-1] << endl;
      cout << "LO UL            " << ul_LO[m0bin-1][m12bin-1] << endl;
     
      hUL_NLO->SetBinContent( m0bin , m12bin , ul_NLO[m0bin-1][m12bin-1] );
      hUL_NLO_exp->SetBinContent( m0bin , m12bin , ul_NLO_exp[m0bin-1][m12bin-1] );
      hUL_NLO_expp1->SetBinContent( m0bin , m12bin , ul_NLO_expp1[m0bin-1][m12bin-1] );
      hUL_NLO_expm1->SetBinContent( m0bin , m12bin , ul_NLO_expm1[m0bin-1][m12bin-1] );
      hUL_LO->SetBinContent( m0bin , m12bin , ul_LO[m0bin-1][m12bin-1] );
    }
  }

  //--------------------------------------------------------------------------------------------------------
  // at this point we have calculated the yields and upper limits at each point
  // next step is to build TH2 histos (ie. hexcl_NLO_obs) which are 1 if point is excluded, otherwise 0
  //--------------------------------------------------------------------------------------------------------
 
  float xmin = hyield_k->GetXaxis()->GetXmin();
  float xmax = hyield_k->GetXaxis()->GetXmax();
  float ymin = hyield_k->GetYaxis()->GetXmin();
  float ymax = hyield_k->GetYaxis()->GetXmax();
  int nx     = hyield_k->GetXaxis()->GetNbins();
  int ny     = hyield_k->GetYaxis()->GetNbins();

  TH2I* hexcl_NLO_obs       = new TH2I("hexcl_NLO_obs",       "Observed NLO Exclusion",nx,xmin,xmax,ny,ymin,ymax);
  TH2I* hexcl_NLO_obs_SC    = new TH2I("hexcl_NLO_obs_SC",    "Observed NLO Exclusion (Sig Cont)",nx,xmin,xmax,ny,ymin,ymax);
  TH2I* hexcl_NLO_obs_nobkg = new TH2I("hexcl_NLO_obs_nobkg", "Observed NLO Exclusion (No Bkg)",nx,xmin,xmax,ny,ymin,ymax);
  TH2I* hexcl_LO_obs        = new TH2I("hexcl_LO_obs",        "Observed LO Exclusion", nx,xmin,xmax,ny,ymin,ymax);
  TH2I* hexcl_NLO_exp       = new TH2I("hexcl_NLO_exp",       "Expected NLO Exclusion",nx,xmin,xmax,ny,ymin,ymax);
  TH2I* hexcl_NLO_expp1     = new TH2I("hexcl_NLO_expp1",     "Expected(+1) NLO Exclusion",nx,xmin,xmax,ny,ymin,ymax);
  TH2I* hexcl_NLO_expm1     = new TH2I("hexcl_NLO_expm1",     "Expected(-1) NLO Exclusion",nx,xmin,xmax,ny,ymin,ymax);

  for( unsigned int m0bin  = 1 ; m0bin  <= nm0bins  ; ++m0bin  ){
    for( unsigned int m12bin = 1 ; m12bin <= nm12bins ; ++m12bin ){

      hexcl_NLO_obs->SetBinContent( m0bin , m12bin , 0 );
      hexcl_NLO_obs_SC->SetBinContent( m0bin , m12bin , 0 );
      hexcl_NLO_obs_nobkg->SetBinContent( m0bin , m12bin , 0 );
      hexcl_LO_obs->SetBinContent( m0bin , m12bin , 0 );
      hexcl_NLO_exp->SetBinContent( m0bin , m12bin , 0 );
      hexcl_NLO_expp1->SetBinContent( m0bin , m12bin , 0 );
      hexcl_NLO_expm1->SetBinContent( m0bin , m12bin , 0 );

      if( doCorrection ){

	int ngen = hscan->GetBinContent(m0bin,m12bin);

	if( ngen != 10000 ){
	  //cout << "Skipping point with " << ngen << " entries" << endl;
	  hexcl_NLO_obs->SetBinContent( m0bin , m12bin , 2 );
	  hexcl_NLO_exp->SetBinContent( m0bin , m12bin , 2 );
	  continue;
	}
      }

      //---------------------------
      //NLO observed
      //---------------------------

      int excluded_NLO_obs = 0;
      if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO[m0bin-1][m12bin-1] ) excluded_NLO_obs = 1;
      hexcl_NLO_obs->SetBinContent( m0bin , m12bin , excluded_NLO_obs );

      //---------------------------
      //NLO observed (sig cont)
      //---------------------------

      int excluded_NLO_obs_SC = 0;
      if( doABCD ){
	if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_SC[m0bin-1][m12bin-1] ) excluded_NLO_obs_SC = 1;
	hexcl_NLO_obs_SC->SetBinContent( m0bin , m12bin , excluded_NLO_obs_SC );
      }

      //---------------------------
      //NLO observed (nobkg)
      //---------------------------

      int excluded_NLO_obs_nobkg = 0;
      if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_nobkg[m0bin-1][m12bin-1] ) excluded_NLO_obs_nobkg = 1;
      hexcl_NLO_obs_nobkg->SetBinContent( m0bin , m12bin , excluded_NLO_obs_nobkg );

      //---------------------------
      //NLO expected
      //---------------------------

      int excluded_NLO_exp = 0;
      if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_exp[m0bin-1][m12bin-1] ) excluded_NLO_exp = 1;
      hexcl_NLO_exp->SetBinContent( m0bin , m12bin , excluded_NLO_exp );

      //---------------------------
      //NLO expected(+1)
      //---------------------------

      int excluded_NLO_expp1 = 0;
      if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_expp1[m0bin-1][m12bin-1] ) excluded_NLO_expp1 = 1;
      hexcl_NLO_expp1->SetBinContent( m0bin , m12bin , excluded_NLO_expp1 );

      //---------------------------
      //NLO expected(-1)
      //---------------------------

      int excluded_NLO_expm1 = 0;
      if( hyield_k->GetBinContent( m0bin , m12bin ) > ul_NLO_expm1[m0bin-1][m12bin-1] ) excluded_NLO_expm1 = 1;
      hexcl_NLO_expm1->SetBinContent( m0bin , m12bin , excluded_NLO_expm1 );

      //---------------------------
      //print out stuff
      //---------------------------

      if( hyield_k->GetBinContent( m0bin , m12bin ) > 0. ){
     
	cout << endl << "m0 " << m0bin-1 << " m12 " << m12bin-1 << endl;
	cout << "NLO Yield     " << hyield_k->GetBinContent( m0bin , m12bin ) << endl;
       
	cout << "NLO UL        " << ul_NLO[m0bin-1][m12bin-1] << endl;
	cout << "Excluded?     " << excluded_NLO_obs << endl;
       
	cout << "NLO UL SC     " << ul_NLO_SC[m0bin-1][m12bin-1] << endl;
	cout << "Excluded?     " << excluded_NLO_obs_SC << endl;

	cout << "NLO UL no bkg " << ul_NLO_nobkg[m0bin-1][m12bin-1] << endl;
	cout << "Excluded?     " << excluded_NLO_obs_nobkg << endl;

	cout << "NLO UL exp    " << ul_NLO_exp[m0bin-1][m12bin-1] << endl;
	cout << "Excluded?     " << excluded_NLO_exp << endl;
       
	// cout << "LO UL         " << ul_LO[m0bin-1][m12bin-1] << endl;
	// cout << "Excluded?     " << excluded_LO_obs << endl;
       
      }
    }
  }


  //----------------------------------------------------------------------------------------
  // these objects specify the exclusion contours, they are built from the hexcl_* histos
  //----------------------------------------------------------------------------------------

  TH1F*         limit_NLO_obs			= getCurve(        hexcl_NLO_obs	, "limit_NLO_obs");
  TGraphErrors* limitgraph_NLO_obs		= getCurve_TGraph( hexcl_NLO_obs	, "limitgraph_NLO_obs");

  TH1F*         limit_NLO_obs_SC		= getCurve(        hexcl_NLO_obs_SC	, "limit_NLO_obs_SC");
  TGraphErrors* limitgraph_NLO_obs_SC		= getCurve_TGraph( hexcl_NLO_obs_SC	, "limitgraph_NLO_obs_SC");

  TH1F*         limit_NLO_obs_nobkg		= getCurve(        hexcl_NLO_obs_nobkg	, "limit_NLO_obs_nobkg");
  TGraphErrors* limitgraph_NLO_obs_nobkg	= getCurve_TGraph( hexcl_NLO_obs_nobkg	, "limitgraph_NLO_obs_nobkg");

  TH1F*         limit_NLO_exp			= getCurve(        hexcl_NLO_exp	, "limit_NLO_exp");
  TGraphErrors* limitgraph_NLO_exp		= getCurve_TGraph( hexcl_NLO_exp	, "limitgraph_NLO_exp");

  TH1F*         limit_NLO_expp1			= getCurve(        hexcl_NLO_expp1	, "limit_NLO_expp1");
  TGraphErrors* limitgraph_NLO_expp1		= getCurve_TGraph( hexcl_NLO_expp1	, "limitgraph_NLO_expp1");

  TH1F*         limit_NLO_expm1			= getCurve(        hexcl_NLO_expm1	, "limit_NLO_expm1");
  TGraphErrors* limitgraph_NLO_expm1		= getCurve_TGraph( hexcl_NLO_expm1	, "limitgraph_NLO_expm1");

  TH1F*         limit_LO_obs			= getCurve(        hexcl_LO_obs		, "limit_LO_obs");
  TGraphErrors* limitgraph_LO_obs		= getCurve_TGraph( hexcl_LO_obs		, "limitgraph_LO_obs");

  //----------------------------------
  // now draw a bunch of plots
  //----------------------------------

  //excluded points
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->Divide(2,2);

  TLatex *t = new TLatex();
  t->SetNDC();

  c1->cd(1);
  hexcl_NLO_obs->Draw("colz");
  t->DrawLatex(0.5,0.8,"observed NLO");

  c1->cd(2);
  hexcl_NLO_exp->Draw("colz");
  t->DrawLatex(0.5,0.8,"expected NLO");
 
  c1->cd(3);
  hexcl_NLO_expp1->Draw("colz");
  t->DrawLatex(0.5,0.8,"expected(+1) NLO");

  c1->cd(4);
  hexcl_NLO_expm1->Draw("colz");
  t->DrawLatex(0.5,0.8,"expected(-1) NLO");
 
  c1->Print("exclusion/exclusion.png");

  /*
  //exclusion TH1
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->cd();

  limit_NLO_obs->SetMinimum(90.);
  limit_NLO_obs->SetTitle("Exclusion Curve in mSUGRA Space");
  limit_NLO_obs->GetXaxis()->SetTitle("m_{0} (GeV)");
  limit_NLO_obs->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  //limit_NLO_obs->SetLineWidth(2);
  //limit_NLO_obs->SetLineColor(2);
  limit_NLO_obs->GetXaxis()->SetRangeUser(0,500);
  limit_NLO_obs->GetYaxis()->SetRangeUser(100,400);
  //limit_NLO_obs->Draw("c");
  limit_NLO_obs->Draw();
  limit_LO_obs->SetLineColor(2);
  limit_NLO_exp->SetLineColor(4);
  limit_LO_obs->Draw("same");
  limit_NLO_exp->Draw("same");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(limit_NLO_obs, "NLO obs","l");
  leg->AddEntry(limit_LO_obs,  "LO obs","l");
  leg->AddEntry(limit_NLO_exp, "NLO exp","l");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  c2->Print("exclusion/limit_hist.png");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->cd();

  limitgraph_LO_obs->SetMarkerColor(2);
  limitgraph_NLO_exp->SetMarkerColor(4);
  limitgraph_NLO_obs->Draw("AP");
  limitgraph_LO_obs->Draw("sameP");
  limitgraph_NLO_exp->Draw("sameP");
  leg->Draw();

  c3->Print("exclusion/limit_graph.png");
  */

  outfile->cd();

  hexcl_NLO_obs->Write();
  limit_NLO_obs->Write();
  limitgraph_NLO_obs->Write();

  hexcl_NLO_obs_SC->Write();
  limit_NLO_obs_SC->Write();
  limitgraph_NLO_obs_SC->Write();

  hexcl_NLO_obs_nobkg->Write();
  limit_NLO_obs_nobkg->Write();
  limitgraph_NLO_obs_nobkg->Write();

  hexcl_LO_obs->Write();
  limit_LO_obs->Write();
  limitgraph_LO_obs->Write();

  hexcl_NLO_exp->Write();
  limit_NLO_exp->Write();
  limitgraph_NLO_exp->Write();

  hUL_NLO->Write();
  hUL_NLO_exp->Write();
  hUL_LO->Write();
  hyield_k->Write();

  htotuncertainty->Write();
  outfile->Close();
 
}

TH1F* getCurve( TH2I *hist , char* name){

  //--------------------------------------------
  // Get the parameters of the yield histogram
  //--------------------------------------------
  float xmin = hist->GetXaxis()->GetXmin();
  float xmax = hist->GetXaxis()->GetXmax();
  float ymin = hist->GetYaxis()->GetXmin();
  float ymax = hist->GetYaxis()->GetXmax();
  int nx     = hist->GetXaxis()->GetNbins();
  int ny     = hist->GetYaxis()->GetNbins();

  cout << "xmin = " << xmin <<endl;
  cout << "xmax = " << xmax <<endl;
  cout << "ymin = " << ymin <<endl;
  cout << "ymax = " << ymax <<endl;
  cout << "nx  = " << nx << endl;
  cout << "ny  = " << ny << endl;

  //---------------------------------------------------------------
  // We book a 1D histogram to keep the results in
  //---------------------------------------------------------------
  TH1F* limit = new TH1F(name,name,nx,xmin,xmax);
  //---------------------------------------------------------------
  // Use Sanjay's method, scan from the top and hope for the best
  //---------------------------------------------------------------
  float ybinsize = (ymax-ymin)/ny;
  float xbinsize = (xmax-xmin)/nx;
  for (int ix=1; ix<=nx; ix++) {
    float x = xmin + ix*xbinsize - 0.5*xbinsize;
    bool foundOne = false;

    for (int iy=ny; iy>0; iy--) {
      float this_ = hist->GetBinContent(ix,iy);
      //this_ = kfact*fudge*this_;
      //if (this_ > nev) {
      if (this_ > 0.5) {
	float yupperedge = ymin + iy*ybinsize;
	limit->Fill(x,yupperedge);
	// cout << yupperedge << " " << iy << endl;
	foundOne=true;
	break;
      }
    } //close iy loop
    //    if (!foundOne) limit->Fill(x,ymin);

  }   //close ix loop
 
  return limit;
}

TGraphErrors* getCurve_TGraph( TH2I *hist , char* name ){

  //--------------------------------------------
  // Get the parameters of the yield histogram
  //--------------------------------------------
  float xmin = hist->GetXaxis()->GetXmin();
  float xmax = hist->GetXaxis()->GetXmax();
  float ymin = hist->GetYaxis()->GetXmin();
  float ymax = hist->GetYaxis()->GetXmax();
  int nx     = hist->GetXaxis()->GetNbins();
  int ny     = hist->GetYaxis()->GetNbins();

  cout << "xmin = " << xmin <<endl;
  cout << "xmax = " << xmax <<endl;
  cout << "ymin = " << ymin <<endl;
  cout << "ymax = " << ymax <<endl;
  cout << "nx  = " << nx << endl;
  cout << "ny  = " << ny << endl;

  const unsigned int npoints = nx;
  float xpoint[npoints];
  float ypoint[npoints];
  float ex[npoints];
  float ey[npoints];

  for( unsigned int i = 0 ; i < npoints ; ++i ){
    xpoint[i] = 0.;
    ypoint[i] = 0.;
    ex[i]     = 0.;
    ey[i]     = 0.;
  }


  //---------------------------------------------------------------
  // Use Sanjay's method, scan from the top and hope for the best
  //---------------------------------------------------------------
  float ybinsize = (ymax-ymin)/ny;
  float xbinsize = (xmax-xmin)/nx;
  for (int ix=1; ix<=nx; ix++) {
    float x = xmin + ix*xbinsize - 0.5*xbinsize;
    bool foundOne = false;

    for (int iy=ny; iy>0; iy--) {
      float this_ = hist->GetBinContent(ix,iy);
      //this_ = kfact*fudge*this_;
      //if (this_ > nev) {
      if (this_ > 0.5) {
	float yupperedge = ymin + iy*ybinsize;
	// limit->Fill(x,yupperedge);
	xpoint[ix-1] = x;
	ypoint[ix-1] = yupperedge;
	// cout << yupperedge << " " << iy << endl;
	foundOne=true;
	break;
      }
    } //close iy loop
    if (!foundOne){
      xpoint[ix-1] = x;
      ypoint[ix-1] = 0;
    }
    //limit->Fill(x,ymin);
  }   //close ix loop
 
  //---------------------------------------------------------------
  // We book a TGraph to store the results
  //---------------------------------------------------------------
  TGraphErrors *gr = new TGraphErrors(npoints,xpoint,ypoint,ex,ey);
  gr->SetName(name);
  gr->SetTitle(name);

  return gr;
 
}
