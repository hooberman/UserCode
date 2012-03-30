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

char* version             = (char*) "V00-02-04";

bool fileInList(string thisfilename);

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
      
    char* filename = Form("cards/%s/SMS_%i_m2lnQ2.root",version,mgbin);
      
    bool found = fileInList( filename );
    if( !found ) continue;
      
    //------------------------------------------
    // check if point is excluded
    //------------------------------------------
      
    limitResult mylimit = run(filename,"plot");
      
    if( mylimit.obs < 1.e-10 ){
      *doScript_failed << Form("../../../../test/lands.exe -d SMS_%i.txt -M Hybrid --freq  --nToysForCLsb 1500 --nToysForCLb 500  --scanRs 1 -vR [0.005,0.5,x1.1] -n SMS_%i",mgbin,mgbin) << endl;
    }
      
    else{

      mgvec.push_back(mg);
      expvec.push_back(1000 * mylimit.exp);
      obsvec.push_back(1000 * mylimit.obs);
      
      cout << "---------------------------------------------------" << endl;
      cout << "mgbin       " << mgbin       << endl;
      cout << "mg          " << mg          << endl;
      cout << "Observed:   " << mylimit.obs << endl; 
      cout << "Expected:   " << mylimit.exp << endl; 
      cout << "---------------------------------------------------" << endl;

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

  TFile* outfile = TFile::Open(Form("cards/%s/observed_limit.root",version),"RECREATE");
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
