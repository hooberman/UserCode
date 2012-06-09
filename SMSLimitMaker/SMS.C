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

float getObservedLimit( int metcut , float seff );
float getExpectedLimit( int metcut , float seff );

using namespace std;

void SMS( bool print = false){

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------
  
  const float denom    = 20000;   // number of generated events per point
  const float lumi     = 4980;    // lumi, in pb
  const char* sample   = "T5zz";  // name of sample
  const char* path     = "./";    // path to babies
  const char* filename = Form("%s/%s_baby.root",path,sample); // input filename

  cout << "Using file        " << filename << endl;
  cout << "Using denominator " << denom    << " events" << endl;
  cout << "Using lumi        " << lumi     << " pb-1" << endl;

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("T1");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("gluino");

  //--------------------------------------------------
  // preselection: store the preselection in TCut presel
  //--------------------------------------------------

  TCut zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets >= 2");
  TCut njets3("njets >= 3");
  TCut sf("leptype==0||leptype==1");
  TCut met30 ("pfmet>30");
  TCut met60 ("pfmet>60");
  TCut met100("pfmet>100");
  TCut met200("pfmet>200");
  TCut met300("pfmet>300");

  TCut presel  = zmass + njets2 + sf;

  cout << "Using selection   " << presel.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions: store the signal region cuts in vector<TCut> sigcuts
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<int>     cuts;

  sigcuts.push_back(TCut(presel+met100));     signames.push_back("E_{T}^{miss} > 100 GeV");     labels.push_back("met100"); cuts.push_back(100);
  sigcuts.push_back(TCut(presel+met200));     signames.push_back("E_{T}^{miss} > 200 GeV");     labels.push_back("met200"); cuts.push_back(200);
  sigcuts.push_back(TCut(presel+met300));     signames.push_back("E_{T}^{miss} > 300 GeV");     labels.push_back("met300"); cuts.push_back(300);

  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* heff[nsig];
  TH2F* heffup[nsig];
  TH2F* heffdn[nsig];
  TH2F* hxsec[nsig];
  TH2F* hxsec_exp[nsig];
  TH2F* hexcl[nsig];
  TH2F* hjes[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    //----------------------------------------------------------
    // calculate JES uncertainty by varying njets and pfmet
    //----------------------------------------------------------

    TString jesup(sigcuts.at(i));
    jesup.ReplaceAll("njets","njetsup");
    jesup.ReplaceAll("pfmet","pfmetup");

    TString jesdn(sigcuts.at(i));
    jesdn.ReplaceAll("njets","njetsdn");
    jesdn.ReplaceAll("pfmet","pfmetdn");

    TCut jesupcut(jesup);
    TCut jesdncut(jesdn);

    cout << endl << endl;
    cout << "Signal region : " << labels.at(i)  << endl;
    cout << "Selection     : " << sigcuts.at(i) << endl;
    cout << "Selection up  : " << jesupcut      << endl;
    cout << "Selection dn  : " << jesdncut      << endl;

    heff[i]      = new TH2F(Form("heff_%i",i)        , Form("heff_%i",i)       , 48,0,1200,48,0,1200);
    heffup[i]    = new TH2F(Form("heffup_%i",i)      , Form("heffup_%i",i)     , 48,0,1200,48,0,1200);
    heffdn[i]    = new TH2F(Form("heffdn_%i",i)      , Form("heffdn_%i",i)     , 48,0,1200,48,0,1200);
    hxsec[i]     = new TH2F(Form("hxsec_%i",i)       , Form("hxsec_%i",i)      , 48,0,1200,48,0,1200);
    hxsec_exp[i] = new TH2F(Form("hxsec_exp_%i",i)   , Form("hxsec_exp_%i",i)  , 48,0,1200,48,0,1200);
    hexcl[i]     = new TH2F(Form("hexcl_%i",i)       , Form("hexcl_%i",i)      , 48,0,1200,48,0,1200);
    hjes[i]      = new TH2F(Form("hjes_%i",i)        , Form("hjes_%i",i)       , 48,0,1200,48,0,1200);

    //--------------------------------------------------------------------------
    // here ml is the LSP mass and mg is the gluino mass
    // store the efficiency across 2D SMS space, nominal, JES up, JES down
    //--------------------------------------------------------------------------
    
    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i));
    heff[i]->Scale(0.95/denom);

    ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut);
    heffup[i]->Scale(0.95/denom);

    ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut);
    heffdn[i]->Scale(0.95/denom);

    //----------------------------------------------------------------------
    // at this point we have the efficiency for each bin the SMS space
    // now we loop over the 2D points and calculate the limit at each point
    //----------------------------------------------------------------------
    
    for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin)-12.5;  // gluino mass
	float ml = heff[i]->GetYaxis()->GetBinCenter(jbin)-12.5;  // LSP mass

	float eff    = heff[i]->GetBinContent(ibin,jbin);
	float effup  = heffup[i]->GetBinContent(ibin,jbin);
	float effdn  = heffdn[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	//calculate the JES uncertainty, store in float djes
	float dup    = effup/eff-1;
	float ddn    = 1-effdn/eff;
	float djes   = 0.5 * (dup+ddn);
	hjes[i]->SetBinContent(ibin,jbin,djes);

	//total uncertainty: JES, lumi, and leptons
	float toterr = sqrt( 0.022*0.022 + 0.05*0.05 + djes*djes );

	// get observed yield UL, 3 inputes:
	// cuts.at(i) : this is just the MET cut defining the signal region
	// toterr     : total signal efficiency uncertainty
	float this_ul = getObservedLimit( cuts.at(i) , toterr );

	// UL(xsec) = UL(yield)/( lumi X eff X 0.19 )
	// 0.19 accounts for the fact that efficiency is defined w.r.t. 
	// #events with >=1 Z->ll, but xsec limit is for inclusive Z decays
	float xsecul  = this_ul / ( lumi * eff * 0.19 );

	// ditto for expected UL
	float this_ul_exp = getExpectedLimit( cuts.at(i) , toterr );
	float xsecul_exp  = this_ul_exp / ( lumi * eff * 0.19 );

	// store observed and expected xsec UL's
	if( eff > 0 ){
	  hxsec[i]->SetBinContent(ibin,jbin, xsecul );
	  hxsec_exp[i]->SetBinContent(ibin,jbin, xsecul_exp );
	}

	// get cross section, check whether point is excluded
	int   bin  = refxsec->FindBin(mg);
	float xsec = refxsec->GetBinContent(bin);

	hexcl[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul )   hexcl[i]->SetBinContent(ibin,jbin,1); // point is excluded

	//cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
      }
    }
  }

  delete ctemp;

  cout << endl << endl;
  
  //-----------------------------------------------------------------------
  // make root file
  // store observed and expected limits and signal efficiency TH2's
  //-----------------------------------------------------

  TFile *outfile = TFile::Open(Form("%s_histos.root",sample),"RECREATE");
  outfile->cd();
  for( unsigned int i = 0 ; i < nsig ; ++i ){
    hxsec[i]->Write();
    heff[i]->Write();
    hxsec_exp[i]->Write();
  }
  outfile->Close();

}

//-------------------------------------------------------------------
// these functions store the observed and expected limits
// metcut is just the MET requirement defining the signal region
// seff is the total signal efficiency uncertainty
//-------------------------------------------------------------------

float getObservedLimit( int metcut , float seff ){

  float ul = 999;

  if( metcut == 100 ){
    if(seff >= 0.0 && seff < 0.1) ul = 57.6;
    if(seff >= 0.1 && seff < 0.2) ul = 59.9;
    if(seff >= 0.2 && seff < 0.3) ul = 64.5;
    if(seff >= 0.3 && seff < 0.4) ul = 69.0;
    if(seff >= 0.4 && seff < 0.5) ul = 74.5;
    if(seff >= 0.5 && seff < 0.6) ul = 79.9;
    if(seff >= 0.6 && seff < 0.7) ul = 86.1;
    if(seff >= 0.7 && seff < 0.8) ul = 92.1;
    if(seff >= 0.8 && seff < 0.9) ul = 101.2;
  }
  else if( metcut == 200 ){
    if(seff >= 0.0 && seff < 0.1) ul = 8.3;
    if(seff >= 0.1 && seff < 0.2) ul = 8.9;
    if(seff >= 0.2 && seff < 0.3) ul = 8.8;
    if(seff >= 0.3 && seff < 0.4) ul = 9.5;
    if(seff >= 0.4 && seff < 0.5) ul = 10.1;
    if(seff >= 0.5 && seff < 0.6) ul = 10.7;
    if(seff >= 0.6 && seff < 0.7) ul = 11.2;
    if(seff >= 0.7 && seff < 0.8) ul = 11.6;
    if(seff >= 0.8 && seff < 0.9) ul = 12.3;
  }
  else if( metcut == 300 ){
    if(seff >= 0.0 && seff < 0.1) ul = 2.8;
    if(seff >= 0.1 && seff < 0.2) ul = 2.8;
    if(seff >= 0.2 && seff < 0.3) ul = 2.8;
    if(seff >= 0.3 && seff < 0.4) ul = 3.0;
    if(seff >= 0.4 && seff < 0.5) ul = 3.1;
    if(seff >= 0.5 && seff < 0.6) ul = 3.3;
    if(seff >= 0.6 && seff < 0.7) ul = 3.5;
    if(seff >= 0.7 && seff < 0.8) ul = 3.6;
    if(seff >= 0.8 && seff < 0.9) ul = 3.6;
  }  
  

  if( ul > 998 ){
    cout << "Error ul " << ul << " metcut " << metcut << " SEFF " << seff << endl;
  }

  return ul;
}



float getExpectedLimit( int metcut , float seff ){

  float ul = 999;

  if( metcut == 100 ){
    if(seff >= 0.0 && seff < 0.1) ul = 61.3;
    if(seff >= 0.1 && seff < 0.2) ul = 63.8;
    if(seff >= 0.2 && seff < 0.3) ul = 68.1;
    if(seff >= 0.3 && seff < 0.4) ul = 73.6;
    if(seff >= 0.4 && seff < 0.5) ul = 79.5;
    if(seff >= 0.5 && seff < 0.6) ul = 86.3;
    if(seff >= 0.6 && seff < 0.7) ul = 93.3;
    if(seff >= 0.7 && seff < 0.8) ul = 100.3;
    if(seff >= 0.8 && seff < 0.9) ul = 109.5;
  }
  else if( metcut == 200 ){
    if(seff >= 0.0 && seff < 0.1) ul = 11.4;
    if(seff >= 0.1 && seff < 0.2) ul = 11.8;
    if(seff >= 0.2 && seff < 0.3) ul = 12.1;
    if(seff >= 0.3 && seff < 0.4) ul = 13.5;
    if(seff >= 0.4 && seff < 0.5) ul = 13.7;
    if(seff >= 0.5 && seff < 0.6) ul = 14.8;
    if(seff >= 0.6 && seff < 0.7) ul = 16.1;
    if(seff >= 0.7 && seff < 0.8) ul = 17.2;
    if(seff >= 0.8 && seff < 0.9) ul = 18.2;
  }
  else if( metcut == 300 ){
    if(seff >= 0.0 && seff < 0.1) ul = 4.8;
    if(seff >= 0.1 && seff < 0.2) ul = 4.9;
    if(seff >= 0.2 && seff < 0.3) ul = 5.1;
    if(seff >= 0.3 && seff < 0.4) ul = 5.3;
    if(seff >= 0.4 && seff < 0.5) ul = 5.7;
    if(seff >= 0.5 && seff < 0.6) ul = 5.9;
    if(seff >= 0.6 && seff < 0.7) ul = 6.3;
    if(seff >= 0.7 && seff < 0.8) ul = 6.8;
    if(seff >= 0.8 && seff < 0.9) ul = 7.0;
  }  
  else{
    cout << "ERROR! unrecognized met cut " << metcut << ", quitting" << endl;
    exit(0);
  }

  if( ul > 998 ){
    cout << "Error ul " << ul << " metcut " << metcut << " SEFF " << seff << endl;
  }

  return ul;

}




