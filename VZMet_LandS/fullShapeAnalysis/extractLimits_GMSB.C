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

char* version             = (char*) "V00-02-08";
bool combined             = true;

bool fileInList(string thisfilename);

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

void extractLimits_GMSB(){

  //------------------------------------------
  // create exclusion histogram
  //------------------------------------------

  ofstream* doScript_failed = new ofstream();
  doScript_failed->open(Form("cards/%s/doLimits_failed.sh",version));

  //------------------------------------------
  // loop over CMSSM points
  //------------------------------------------

  const unsigned int nbins = 15;
  
  vector<float> mgvec;
  vector<float> expvec;
  vector<float> obsvec;

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
      
      float xsec = crossSectionGMSB(mg);
      mgvec.push_back(mg);
      expvec.push_back(xsec * mylimit.exp);
      obsvec.push_back(xsec * mylimit.obs);
      
      cout << "---------------------------------------------------" << endl;
      cout << "filename    " << filename    << endl;
      cout << "mgbin       " << mgbin       << endl;
      cout << "xsec        " << xsec        << endl;
      cout << "mg          " << mg          << endl;
      cout << "Observed:   " << xsec * mylimit.obs << endl; 
      cout << "Expected:   " << xsec * mylimit.exp << endl; 
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

  for( unsigned int i0 = 0 ; i0 < n ; ++i0 ){
    mg[i0]    = mgvec.at(i0);
    exp[i0]   = expvec.at(i0);
    obs[i0]   = obsvec.at(i0);
  }

  TGraph grobs(n,mg,obs);    
  TGraph grexp(n,mg,exp);

  grobs.SetName("grobs");
  grexp.SetName("grexp");

  grobs.SetTitle("grobs");
  grexp.SetTitle("grexp");

  char* outfilename         = Form("cards/%s/observed_limit.root",version);
  if( combined) outfilename = Form("cards/%s/observed_limit_combined.root",version);

  TFile* outfile = TFile::Open(outfilename,"RECREATE");

  grobs.Write();
  grexp.Write();
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
