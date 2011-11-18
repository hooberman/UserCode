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
#include "TROOT.h"

//#include "cl95cms_landsberg.c"

using namespace std;


//-------------------------------------------
//parameters for CL95 function
//-------------------------------------------

const Double_t mylumi        = 3.5;
const bool     doCorrection  = true;

//-------------------------------------------
// uncertainties
//-------------------------------------------

const float lumierr = 0.06;  // 6% lumi error
const float leperr  = 0.05;  // lepton ID/iso + trigger uncertainty
const float pdferr  = 0.20;  // PDF uncertainty

//----------------------------------------------
// ABCD yields (for signal contamination study)
//----------------------------------------------

const bool  doABCD  = false; //turned off by default
const float yieldA  =  12.;
const float yieldB  =  37.;
const float yieldC  =   4.;

//-----------------------------------------------
//rebin the TH2 yield histos (NOT RECOMMENDED)
//-----------------------------------------------

const int   rebin   = 1;     // rebin yield histos


float getUpperLimit( float seff ){

  float ul = 9999.;

  //-----------
  // CLs
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 11.8;
  if(seff > 0.1 && seff < 0.2) ul = 12.4;
  if(seff > 0.2 && seff < 0.3) ul = 13.6;
  if(seff > 0.3 && seff < 0.4) ul = 14.5;
  if(seff > 0.4 && seff < 0.5) ul = 16.0;
  if(seff > 0.5 && seff < 0.6) ul = 17.6;
  if(seff > 0.6 && seff < 0.7) ul = 19.4;
  if(seff > 0.7 && seff < 0.8) ul = 21.0;
  if(seff > 0.8 && seff < 0.9) ul = 22.7;

  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;
}

float getExpectedUpperLimit( float seff ){

  float ul = 9999.;

  //-----------
  // expected
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 10.3;
  if(seff > 0.1 && seff < 0.2) ul = 10.9;
  if(seff > 0.2 && seff < 0.3) ul = 11.9;
  if(seff > 0.3 && seff < 0.4) ul = 12.7;
  if(seff > 0.4 && seff < 0.5) ul = 13.7;
  if(seff > 0.5 && seff < 0.6) ul = 14.9;
  if(seff > 0.6 && seff < 0.7) ul = 16.8;
  if(seff > 0.7 && seff < 0.8) ul = 17.7;
  if(seff > 0.8 && seff < 0.9) ul = 18.6;

  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;
}

float getExpectedP1UpperLimit( float seff ){

  float ul = 9999.;

  //-----------
  // expected
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 8.1;
  if(seff > 0.1 && seff < 0.2) ul = 8.1;
  if(seff > 0.2 && seff < 0.3) ul = 8.3;
  if(seff > 0.3 && seff < 0.4) ul = 9.0;
  if(seff > 0.4 && seff < 0.5) ul = 9.8;
  if(seff > 0.5 && seff < 0.6) ul = 10.5;
  if(seff > 0.6 && seff < 0.7) ul = 10.5;
  if(seff > 0.7 && seff < 0.8) ul = 11.2;
  if(seff > 0.8 && seff < 0.9) ul = 13.3;

  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;
}

float getExpectedM1UpperLimit( float seff ){

  float ul = 9999.;

  //-----------
  // expected
  //----------

  if(seff > 0.0 && seff < 0.1) ul = 14.4;
  if(seff > 0.1 && seff < 0.2) ul = 15.7;
  if(seff > 0.2 && seff < 0.3) ul = 17.0;
  if(seff > 0.3 && seff < 0.4) ul = 19.8;
  if(seff > 0.4 && seff < 0.5) ul = 22.4;
  if(seff > 0.5 && seff < 0.6) ul = 25.1;
  if(seff > 0.6 && seff < 0.7) ul = 28.0;
  if(seff > 0.7 && seff < 0.8) ul = 31.5;
  if(seff > 0.8 && seff < 0.9) ul = 35.9;
  
  if( seff > 0.9 ) cout << "Error, signal efficiency uncertainty too large! " << seff << endl;

  return ul;

}

TH1F* getCurve               ( TH2I *hist , char* name );
TGraphErrors* getCurve_TGraph( TH2I *hist , char* name );

void msugra( char* filename , bool print = false ){
  
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

  //-----------------------------------
  // Here we load the yield histograms
  //-----------------------------------
  TFile *f = TFile::Open( filename );
  
  //TH2F* hyield     = (TH2F*) f->Get("msugra",prefix));
  TH2F* hyield_k   = (TH2F*) f->Get("msugra");
  TH2F* hyield_kup = (TH2F*) f->Get("msugra_kup");
  TH2F* hyield_kdn = (TH2F*) f->Get("msugra_kdn");
  TH2F* hyield_jup = (TH2F*) f->Get("msugra_jup");
  TH2F* hyield_jdn = (TH2F*) f->Get("msugra_jdn");

  hyield_k->  Scale(mylumi);
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

      int m0  = hyield_k->GetXaxis()->GetBinCenter(m0bin);
      int m12 = hyield_k->GetXaxis()->GetBinCenter(m12bin);

      if( doCorrection ){
	ngen = hscan->GetBinContent(m0bin,m12bin);

	if( ngen != 10000 ){
	  cout << "Skipping point with " << ngen << " entries" << endl;
	  continue;
	}
      }
     
      cout << endl << endl << "------------------------------------------------------------------------" << endl;
      cout << endl << "m0 " << m0bin-1 << " m12 " << m12bin-1 << endl;
      cout << hyield_k->GetXaxis()->GetBinCenter(m0bin) << " " << hyield_k->GetYaxis()->GetBinCenter(m12bin) << endl;

      //------------------------
      // skip bins with 0 yield
      //------------------------
     
      if( fabs( yield_k ) < 1.e-10 ){
	cout << "zero yield, skipping!" << endl;
	continue; 
      }

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

      ul_NLO[m0bin-1][m12bin-1]       = getUpperLimit          ( accerr_NLO[m0bin-1][m12bin-1] );
      ul_NLO_exp[m0bin-1][m12bin-1]   = getExpectedUpperLimit  ( accerr_NLO[m0bin-1][m12bin-1] );
      ul_NLO_expp1[m0bin-1][m12bin-1] = getExpectedP1UpperLimit( accerr_NLO[m0bin-1][m12bin-1] );
      ul_NLO_expm1[m0bin-1][m12bin-1] = getExpectedM1UpperLimit( accerr_NLO[m0bin-1][m12bin-1] );

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
     
      hUL_NLO->SetBinContent( m0bin , m12bin , ul_NLO[m0bin-1][m12bin-1] );
      hUL_NLO_exp->SetBinContent( m0bin , m12bin , ul_NLO_exp[m0bin-1][m12bin-1] );
      hUL_NLO_expp1->SetBinContent( m0bin , m12bin , ul_NLO_expp1[m0bin-1][m12bin-1] );
      hUL_NLO_expm1->SetBinContent( m0bin , m12bin , ul_NLO_expm1[m0bin-1][m12bin-1] );
    }
  }

  cout << __LINE__ << endl;
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
	  hexcl_NLO_obs->SetBinContent( m0bin , m12bin , 2 );
	  hexcl_NLO_exp->SetBinContent( m0bin , m12bin , 2 );
	  hexcl_NLO_expp1->SetBinContent( m0bin , m12bin , 2 );
	  hexcl_NLO_expm1->SetBinContent( m0bin , m12bin , 2 );
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

  TLatex *t = new TLatex();
  t->SetNDC();

  bool overlayTGraph = false;

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  hexcl_NLO_obs->GetXaxis()->SetTitle("m_{0} (GeV)");
  hexcl_NLO_obs->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexcl_NLO_obs->GetXaxis()->SetNdivisions(5);
  hexcl_NLO_obs->Draw("colz");
  t->DrawLatex(0.2,0.92,"observed NLO");
  if( overlayTGraph ){
    TGraphErrors* gr_NLO_obs = getCurve_TGraph( hexcl_NLO_obs , "gr_NLO_obs" );
    gr_NLO_obs->Draw("same");
  }

  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  hexcl_NLO_exp->GetXaxis()->SetTitle("m_{0} (GeV)");
  hexcl_NLO_exp->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexcl_NLO_exp->GetXaxis()->SetNdivisions(5);
  hexcl_NLO_exp->Draw("colz");
  t->DrawLatex(0.2,0.92,"expected NLO");
 
  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  c3->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  hexcl_NLO_expp1->GetXaxis()->SetTitle("m_{0} (GeV)");
  hexcl_NLO_expp1->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexcl_NLO_expp1->GetXaxis()->SetNdivisions(5);
  hexcl_NLO_expp1->Draw("colz");
  t->DrawLatex(0.2,0.92,"expected(+1) NLO");

  TCanvas *c4 = new TCanvas("c4","c4",1000,800);
  c4->cd(1);
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.1);
  hexcl_NLO_expm1->GetXaxis()->SetTitle("m_{0} (GeV)");
  hexcl_NLO_expm1->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexcl_NLO_expm1->GetXaxis()->SetNdivisions(5);
  hexcl_NLO_expm1->Draw("colz");
  t->DrawLatex(0.2,0.92,"expected(-1) NLO");
 
  if( print ){
    c1->Print("../../plots/exclusion_obs.png");
    c2->Print("../../plots/exclusion_exp.png");
    c3->Print("../../plots/exclusion_expp1.png");
    c4->Print("../../plots/exclusion_expm1.png");

    c1->Print("../../plots/exclusion_obs.eps");
    c2->Print("../../plots/exclusion_exp.eps");
    c3->Print("../../plots/exclusion_expp1.eps");
    c4->Print("../../plots/exclusion_expm1.eps");

    gROOT->ProcessLine(".! ps2pdf ../../plots/exclusion_obs.eps   ../../plots/exclusion_obs.pdf");
    gROOT->ProcessLine(".! ps2pdf ../../plots/exclusion_exp.eps   ../../plots/exclusion_exp.pdf");
    gROOT->ProcessLine(".! ps2pdf ../../plots/exclusion_expp1.eps ../../plots/exclusion_expp1.pdf");
    gROOT->ProcessLine(".! ps2pdf ../../plots/exclusion_expm1.eps ../../plots/exclusion_expm1.pdf");
  }


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
  hexcl_NLO_expp1->Write();
  hexcl_NLO_expm1->Write();
  limit_NLO_exp->Write();
  limitgraph_NLO_exp->Write();

  hUL_NLO->Write();
  hUL_NLO_exp->Write();
  hUL_NLO_expp1->Write();
  hUL_NLO_expm1->Write();
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
  gr->SetLineWidth(3);
  gr->SetLineColor(4);

  return gr;
 
}
