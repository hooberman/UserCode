#include "../../test/fitRvsCLs.C"
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

char* version             = "V00-00-04";
//bool  writeExpectedLimits = false;

bool fileInList(string thisfilename);

void extractLimits( bool print = false ){

  //------------------------------------------
  // create exclusion histogram
  //------------------------------------------

  const int   nm0points    = 100;
  const float m0min        = 20.;
  const float m0max        = 2020.;
  const int   nm12points   = 38;
  const float m12min       = 20.;
  const float m12max       = 780.;
  TH2F* hexcl    = new TH2F( "hexcl"    , "hexcl"    , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hexp     = new TH2F( "hexp"     , "hexp"     , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hexpp1   = new TH2F( "hexpp1"   , "hexpp1"   , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* hexpm1   = new TH2F( "hexpm1"   , "hexpm1"   , nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

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

  TFile* f = new TFile();

  // ofstream* ofile    = new ofstream();
  // ofstream* filelist = new ofstream();

  // if( writeExpectedLimits ){
  //   ofile->open(Form("cards/%s/doLimits_expected.sh",version));
  //   filelist->open(Form("cards/%s/file_list_expected.txt",version));
  // }

  ofstream* doScript_failed = new ofstream();
  doScript_failed->open(Form("cards/%s/doLimits_failed.sh",version));

  for( int m0bin = 1 ; m0bin <= hexcl->GetXaxis()->GetNbins() ; m0bin++ ){
    for( int m12bin = 1 ; m12bin <= hexcl->GetYaxis()->GetNbins() ; m12bin++ ){

      //------------------------------------------
      // require nentries = 10
      //------------------------------------------
      
      // int ngen = hscan->GetBinContent(m0bin,m12bin);

      // if( ngen != 10000 ){
      // 	//cout << "Skipping point with " << ngen << " entries" << endl;
      // 	hexcl->SetBinContent(m0bin,m12bin,2);
      // 	continue;
      // }

      //------------------------------------------
      // restrict range
      //------------------------------------------

      int m0  = hexcl->GetXaxis()->GetBinCenter(m0bin);
      int m12 = hexcl->GetXaxis()->GetBinCenter(m12bin);

      //if( m0 == 80 && m12 == 400 ) cout << "Found LM6 " << endl;
      //else continue;
      //if( m0bin >= 50 ) continue;

      hexcl->SetBinContent(m0bin,m12bin,0);

      //------------------------------------------
      // open file, if available
      //------------------------------------------

      char* filename = Form("cards/%s/CMSSM_%i_%i_m2lnQ2.root",version,m0bin,m12bin);
      //char* filename = Form("cards/%s/CMSSM_%i_%i.txt_Bayesian_bysObsLimit.root",version,m0bin,m12bin);

      bool found = fileInList( filename );
      if( !found ) continue;

      f = TFile::Open(filename);

      //------------------------------------------
      // check if point is excluded
      //------------------------------------------

      limitResult mylimit = run(filename,"plot");

      if( mylimit.obs < 1.e-10 ){
	*doScript_failed << Form("../../../../test/lands.exe -d CMSSM_%i_%i.txt -M Hybrid --freq  --nToysForCLsb 1500 --nToysForCLb 500  --scanRs 1 -vR [0.1,10,x1.1] -n CMSSM_%i_%i",m0bin,m12bin,m0bin,m12bin) << endl;
      }
      
      else{
	cout << "Writing limits to histos" << endl;
	cout << "Observed      " << mylimit.obs << endl;
	cout << "Expected      " << mylimit.exp << endl;
	cout << "Expected(+1)  " << mylimit.expp1 << endl;
	cout << "Expected(-1)  " << mylimit.expm1 << endl;

	hexcl-> SetBinContent(m0bin,m12bin,mylimit.obs);
	hexp->  SetBinContent(m0bin,m12bin,mylimit.exp);
	hexpp1->SetBinContent(m0bin,m12bin,mylimit.expp1);
	hexpm1->SetBinContent(m0bin,m12bin,mylimit.expm1);
      }


      // TTree* t = (TTree*) f->Get("T");
      // Double_t limit;
      // t->SetBranchAddress( "limit" , &limit ); 

      // t->GetEntry(0);
      // //cout << "limit " << limit << endl;
      
      // int excluded = 0;
      // if( limit < 1 ) excluded = 1;
      // hexcl->SetBinContent(m0bin,m12bin,excluded);

      // f->Close();

      // if( writeExpectedLimits ){
      // 	if( limit > 0.5 && limit < 3 ){
      // 	  *ofile    << Form("../../../../test/lands.exe -M Bayesian -d CMSSM_%i_%i.txt --doExpectation 1 -t 100",m0bin,m12bin) << endl;
      // 	  *filelist << Form("cards/%s/CMSSM_%i_%i.txt_BayesianBayesian_limitbands.root",version,m0bin,m12bin)                  << endl;
      // 	}
      // }
    }
  }

  doScript_failed->close();

  //------------------------------------------
  // draw exclusion histogram
  //------------------------------------------

  TCanvas *can = new TCanvas("can","can",1000,800);
  can->cd();
  gPad->SetRightMargin(0.2);
  hexcl->GetXaxis()->SetRangeUser(0,2000);
  hexcl->GetYaxis()->SetRangeUser(100,700);
  hexcl->GetXaxis()->SetLabelSize(0.03);
  hexcl->GetYaxis()->SetLabelSize(0.03);
  hexcl->GetXaxis()->SetTitle("m_{0} (GeV)");
  hexcl->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexcl->Draw("colz");

  TFile* outfile = TFile::Open(Form("cards/%s/observed_limit.root",version),"RECREATE");
  hexcl->Write();
  hexp->Write();
  hexpp1->Write();
  hexpm1->Write();
  outfile->Close();

  if( print ){
    can->Print(Form("cards/%s/CMSSM.eps",version));
    can->Print(Form("cards/%s/CMSSM.png",version));
    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/CMSSM.eps cards/%s/CMSSM.pdf",version,version));
  }

}

//------------------------------------------
// check if this file appears in file list
//------------------------------------------

bool fileInList(string thisfilename){

  ifstream* ifile = new ifstream();
  ifile->open(Form("cards/%s/file_list_CLs.txt",version));

  string filename;

  bool found = false;

  while( ifile->good() ){
    *ifile >> filename;
    if( filename == thisfilename ){
      found = true;
      break;
    }
  }

  ifile->close();
  delete ifile;

  return found;
}
