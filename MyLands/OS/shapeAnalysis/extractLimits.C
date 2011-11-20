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

bool fileInList(string thisfilename);

void extractLimits(){

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

  //------------------------------------------
  // open histogram of #entries/CMSSM point
  //------------------------------------------

  TFile* corfile = corfile = TFile::Open("mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0.root");
  TH2F*  hscan   = (TH2F*) corfile->Get("hscan");
  
  if( hscan == 0 ){
    cout << "Can't find TH2 hscan!!!" << endl;
    exit(0);
  }

  //------------------------------------------
  // loop over CMSSM points
  //------------------------------------------

  TFile* f = new TFile();

  for( int m0bin = 1 ; m0bin <= hexcl->GetXaxis()->GetNbins() ; m0bin++ ){
    for( int m12bin = 1 ; m12bin <= hexcl->GetYaxis()->GetNbins() ; m12bin++ ){

      //------------------------------------------
      // require nentries = 10
      //------------------------------------------
      
      int ngen = hscan->GetBinContent(m0bin,m12bin);

      if( ngen != 10000 ){
	cout << "Skipping point with " << ngen << " entries" << endl;
	hexcl->SetBinContent(m0bin,m12bin,2);
	continue;
      }

      //------------------------------------------
      // restrict range
      //------------------------------------------

      int m0  = hexcl->GetXaxis()->GetBinCenter(m0bin);
      int m12 = hexcl->GetXaxis()->GetBinCenter(m12bin);

      if( m0bin > 33 ) continue;
      //if( !( m0bin == 17 && m12bin == 31 ) ) continue;

      hexcl->SetBinContent(m0bin,m12bin,0);

      //------------------------------------------
      // open file, if available
      //------------------------------------------

      char* filename = Form("cards/CMSSM_%i_%i.txt_Bayesian_bysObsLimit.root",m0bin,m12bin);
      cout << "Opening file " << filename << endl;

      bool found = fileInList( filename );
      if( !found ) continue;

      f = TFile::Open(Form("cards/CMSSM_%i_%i.txt_Bayesian_bysObsLimit.root",m0bin,m12bin));

      //------------------------------------------
      // check if point is excluded
      //------------------------------------------

      TTree* t = (TTree*) f->Get("T");
      Double_t limit;
      t->SetBranchAddress( "limit" , &limit ); 

      t->GetEntry(0);
      cout << "limit " << limit << endl;
      
      int excluded = 0;
      if( limit < 1 ) excluded = 1;
      hexcl->SetBinContent(m0bin,m12bin,excluded);

      f->Close();
    }
  }

  //------------------------------------------
  // draw exclusion histogram
  //------------------------------------------

  TCanvas *can = new TCanvas("can","can",1000,800);
  can->cd();
  gPad->SetRightMargin(0.2);
  hexcl->GetXaxis()->SetLabelSize(0.03);
  hexcl->GetYaxis()->SetLabelSize(0.03);
  hexcl->GetXaxis()->SetTitle("m_{0} (GeV)");
  hexcl->GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hexcl->Draw("colz");



}


//------------------------------------------
// check if this file appears in file list
//------------------------------------------

bool fileInList(string thisfilename){

  ifstream* ifile = new ifstream();
  ifile->open("cards/file_list.txt");

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

  if( found )  cout << "FOUND       " << thisfilename << endl;
  else         cout << "Didn't find " << thisfilename << endl;

  return found;
}
