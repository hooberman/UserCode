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
  
  float mgval[nbins];
  float exp0[nbins];
  float obs0[nbins];
  float exp50[nbins];
  float obs50[nbins];

  int mlbin = 1;

  for( int ml = 0 ; ml < 2 ; ml++ ){
    
    if( ml == 0 ) mlbin = 1;
    if( ml == 1 ) mlbin = 3;

    for( int mgbin = 6 ; mgbin < 18 ; mgbin++ ){

      float mg  = ( mgbin - 1 ) * 25;
      mgval[mgbin-6] = mg;

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
	  obs0[mgbin-6]  = mylimit.obs;
	  exp0[mgbin-6]  = mylimit.exp;
	}
	else if( mlbin == 3 ){
	  obs50[mgbin-6] = mylimit.obs;
	  exp50[mgbin-6] = mylimit.exp;
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

  TGraph grobs0(nbins,mgval,obs0);
  TGraph grexp0(nbins,mgval,exp0);
  TGraph grobs50(nbins,mgval,obs50);
  TGraph grexp50(nbins,mgval,exp50);

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
