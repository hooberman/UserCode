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

char* version             = (char*) "V00-00-00";

bool fileInList(string thisfilename);

void extractLimits_line( bool print = false ){

  //------------------------------------------
  // create exclusion histogram
  //------------------------------------------

  ofstream* doScript_failed = new ofstream();
  doScript_failed->open(Form("cards/%s/doLimits_failed.sh",version));

  //------------------------------------------
  // loop over CMSSM points
  //------------------------------------------

  const unsigned int nbins = 12;
  
  vector<float> mg0vec;
  vector<float> mg50vec;
  vector<float> exp0vec;
  vector<float> exp50vec;
  vector<float> obs0vec;
  vector<float> obs50vec;

  int mlbin = 1;

  for( int ml = 0 ; ml < 2 ; ml++ ){
    
    if( ml == 0 ) mlbin = 1;
    if( ml == 1 ) mlbin = 3;

    for( int mgbin = 6 ; mgbin < 18 ; mgbin++ ){

      float mg  = ( mgbin - 1 ) * 25;
      //mgval[mgbin-6] = mg;

      //------------------------------------------
      // open file, if available
      //------------------------------------------
      
      char* filename = Form("cards/%s/SMS_%i_%i_m2lnQ2.root",version,mgbin,mlbin);
      
      bool found = fileInList( filename );
      if( !found ) continue;
      
      //------------------------------------------
      // check if point is excluded
      //------------------------------------------
      
      limitResult mylimit = run(filename,"plot");
      
      if( mylimit.obs < 1.e-10 ){
	*doScript_failed << Form("../../../../test/lands.exe -d SMS_%i_%i.txt -M Hybrid --freq  --nToysForCLsb 1500 --nToysForCLb 500  --scanRs 1 -vR [0.005,0.5,x1.1] -n SMS_%i_%i",mgbin,mlbin,mgbin,mlbin) << endl;
      }
      
      else{

	if( mlbin == 1 ){
	  //obs0[mgbin-6]  = 1000 * mylimit.obs;
	  //exp0[mgbin-6]  = 1000 * mylimit.exp;

	  mg0vec.push_back(mg);
	  exp0vec.push_back(1000 * mylimit.exp);
	  obs0vec.push_back(1000 * mylimit.obs);
	}
	else if( mlbin == 3 ){
	  //obs50[mgbin-6] = 1000 * mylimit.obs;
	  //exp50[mgbin-6] = 1000 * mylimit.exp;

	  mg50vec.push_back(mg);
	  exp50vec.push_back(1000 * mylimit.exp);
	  obs50vec.push_back(1000 * mylimit.obs);
	}

	cout << "---------------------------------------------------" << endl;
	cout << "mgbin mlbin " << mgbin << " " << mlbin << endl;
	cout << "mg          " << mg          << endl;
	cout << "Observed:   " << mylimit.obs << endl; 
	cout << "Expected:   " << mylimit.exp << endl; 
	cout << "---------------------------------------------------" << endl;

      }
    }
  }
  

  doScript_failed->close();
  
  //------------------------------------------
  // save TGraphs
  //------------------------------------------

  const unsigned int n0  = mg0vec.size();
  const unsigned int n50 = mg50vec.size();

  float mg0[n0];
  float obs0[n0];
  float exp0[n0];

  float mg50[n50];
  float obs50[n50];
  float exp50[n50];

  // float mgval[nbins];
  // float exp0[nbins];
  // float obs0[nbins];
  // float exp50[nbins];
  // float obs50[nbins];

  for( int i0 = 0 ; i0 < n0 ; ++i0 ){
    mg0[i0]    = mg0vec.at(i0);
    exp0[i0]   = exp0vec.at(i0);
    obs0[i0]   = obs0vec.at(i0);
  }

  for( int i50 = 0 ; i50 < n50 ; ++i50 ){
    mg50[i50]  = mg50vec.at(i50);
    exp50[i50] = exp50vec.at(i50);
    obs50[i50] = obs50vec.at(i50);
  }

  // TGraph grobs0(nbins,mgval,obs0);    
  // TGraph grexp0(nbins,mgval,exp0);
  // TGraph grobs50(nbins,mgval,obs50);
  // TGraph grexp50(nbins,mgval,exp50);

  TGraph grobs0(n0,mg0,obs0);    
  TGraph grexp0(n0,mg0,exp0);
  TGraph grobs50(n50,mg50,obs50);    
  TGraph grexp50(n50,mg50,exp50);

  grobs0.SetName("grobs0");
  grexp0.SetName("grexp0");
  grobs50.SetName("grobs50");
  grexp50.SetName("grexp50");

  grobs0.SetTitle("grobs0");
  grexp0.SetTitle("grexp0");
  grobs50.SetTitle("grobs50");
  grexp50.SetTitle("grexp50");

  TFile* outfile = TFile::Open(Form("cards/%s/observed_limit.root",version),"RECREATE");
  grobs0.Write();
  grexp0.Write();
  grobs50.Write();
  grexp50.Write();
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
