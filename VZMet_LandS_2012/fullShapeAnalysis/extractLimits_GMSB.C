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

char* version             = (char*) "V00-00-03";
bool combined             = false;

bool fileInList(string thisfilename);

// 8 TeV cross section (in pb)
float crossSectionGMSB( int mu ){

  float xsec = -1;

  if( mu ==110 ) xsec = 7.288;
  if( mu ==130 ) xsec = 3.764;
  if( mu ==150 ) xsec = 2.141;
  if( mu ==170 ) xsec = 1.304;
  if( mu ==190 ) xsec = 0.837;
  if( mu ==210 ) xsec = 0.558;
  if( mu ==230 ) xsec = 0.382;
  if( mu ==250 ) xsec = 0.271;
  if( mu ==270 ) xsec = 0.195;
  if( mu ==290 ) xsec = 0.142;
  if( mu ==310 ) xsec = 0.106;
  if( mu ==330 ) xsec = 0.0798;
  if( mu ==350 ) xsec = 0.0608;
  if( mu ==370 ) xsec = 0.0468;
  if( mu ==390 ) xsec = 0.0366;
  if( mu ==410 ) xsec = 0.0287;
      
  if( xsec < 0 ){
    cout << __FILE__ << " " << __LINE << " ERROR! couldn't find GMSB cross section for mu " << mu << endl;
  }

  return xsec;

}

/*
// 7 TeV cross section (in fb)
float crossSectionGMSB( int m ){

  float xsec = -1;

  if     ( m == 130 ) xsec = 3057;
  else if( m == 150 ) xsec = 1719;
  else if( m == 170 ) xsec = 1035;
  else if( m == 190 ) xsec =  656;
  else if( m == 210 ) xsec =  433;
  else if( m == 230 ) xsec =  293;
  else if( m == 250 ) xsec =  205;
  else if( m == 270 ) xsec =  146;
  else if( m == 290 ) xsec =  105;
  else if( m == 310 ) xsec =   77;
  else if( m == 330 ) xsec =   57;
  else if( m == 350 ) xsec =   43;
  else if( m == 370 ) xsec =   33;
  else if( m == 390 ) xsec =   25;
  else if( m == 410 ) xsec =   20;
  else{
    cout << "ERROR! unrecognized GMSB mass " << m << endl;
  }

  return xsec;
}
*/

void extractLimits_GMSB(){

  //------------------------------------------
  // create exclusion histogram
  //------------------------------------------

  ofstream* doScript_failed = new ofstream();
  doScript_failed->open(Form("cards/%s/doLimits_failed.sh",version));

  //------------------------------------------
  // loop over CMSSM points
  //------------------------------------------

  const unsigned int nbins = 7;
  
  vector<float> mgvec;
  vector<float> expvec;
  vector<float> obsvec;
  vector<float> expp1vec;
  vector<float> expm1vec;

  for( unsigned int mgbin = 1 ; mgbin <= nbins ; mgbin++ ){

    float mg  = ( mgbin - 1 ) * 20 + 130;

    //------------------------------------------
    // open file, if available
    //------------------------------------------
      
    char* filename = "";
    
    if( combined ) filename = Form("cards/%s/SMS_combined_%i_m2lnQ2.root" ,version,mgbin);
    else           filename = Form("cards/%s/SMS_%i_m2lnQ.root"           ,version,mgbin);

    //bool found = fileInList( filename );
    //if( !found ) continue;
      
    //------------------------------------------
    // check if point is excluded
    //------------------------------------------
      
    limitResult mylimit = run(filename,"plot");
      
    if( mylimit.obs < 1.e-10 ){
      *doScript_failed << Form("../../../../test/lands.exe -d SMS_%i.txt -M Hybrid --freq  --nToysForCLsb 1500 --nToysForCLb 500  --scanRs 1 -vR [0.005,0.5,x1.1] -n SMS_%i",mgbin,mgbin) << endl;
    }
      
    else{

      //----------------------
      // normalized to xsec
      //----------------------
      
      //float xsec = crossSectionGMSB(mg);
      float xsec = 1.0;
      mgvec.push_back(mg);
      expvec.push_back(xsec * mylimit.exp);
      obsvec.push_back(xsec * mylimit.obs);
      expp1vec.push_back(xsec * mylimit.expp1);
      expm1vec.push_back(xsec * mylimit.expm1);

      cout << "---------------------------------------------------" << endl;
      cout << "filename       " << filename             << endl;
      cout << "mgbin          " << mgbin                << endl;
      cout << "xsec           " << xsec                 << endl;
      cout << "mg             " << mg                   << endl;
      cout << "Observed:      " << xsec * mylimit.obs   << endl; 
      cout << "Expected:      " << xsec * mylimit.exp   << endl; 
      cout << "Expected(+1):  " << xsec * mylimit.expp1 << endl; 
      cout << "Expected(-1):  " << xsec * mylimit.expm1 << endl; 
      cout << "---------------------------------------------------" << endl;

      //----------------------
      // normalized to 1/pb
      //----------------------

      // mgvec.push_back(mg);
      // expvec.push_back(1000 * mylimit.exp);
      // obsvec.push_back(1000 * mylimit.obs);
      
      // cout << "---------------------------------------------------" << endl;
      // cout << "mgbin       " << mgbin       << endl;
      // cout << "mg          " << mg          << endl;
      // cout << "Observed:   " << 1000 * mylimit.obs << endl; 
      // cout << "Expected:   " << 1000 * mylimit.exp << endl; 
      // cout << "---------------------------------------------------" << endl;

    }
  }

  doScript_failed->close();
  
  //------------------------------------------
  // save TGraphs
  //------------------------------------------

  const unsigned int n  = mgvec.size();

  float mg[n];
  float obs[n];
  float exp[n];
  float expp1[n];
  float expm1[n];

  for( unsigned int i0 = 0 ; i0 < n ; ++i0 ){
    mg[i0]      = mgvec.at(i0);
    exp[i0]     = expvec.at(i0);
    obs[i0]     = obsvec.at(i0);
    expp1[i0]   = expp1vec.at(i0);
    expm1[i0]   = expm1vec.at(i0);
  }

  TGraph grobs(n,mg,obs);    
  TGraph grexp(n,mg,exp);
  TGraph grexpp1(n,mg,expp1);
  TGraph grexpm1(n,mg,expm1);

  grobs.SetName("grobs");
  grexp.SetName("grexp");
  grexpp1.SetName("grexpp1");
  grexpm1.SetName("grexpm1");

  grobs.SetTitle("grobs");
  grexp.SetTitle("grexp");
  grexpp1.SetTitle("grexpp1");
  grexpm1.SetTitle("grexpm1");

  char* outfilename         = Form("cards/%s/observed_limit.root",version);
  if( combined) outfilename = Form("cards/%s/observed_limit_combined_band.root",version);

  TFile* outfile = TFile::Open(outfilename,"RECREATE");

  grobs.Write();
  grexp.Write();
  grexpp1.Write();
  grexpm1.Write();
  outfile->Close();

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
