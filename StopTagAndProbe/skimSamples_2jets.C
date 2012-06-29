#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

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

using namespace std;

void skim(char* version, char* sample);

void skimSamples_2jets(){

  skim("V00-00-03","data_DoubleElectron_May10");
  skim("V00-00-03","data_DoubleElectron_PRv4");
  skim("V00-00-03","data_DoubleElectron_Aug05");
  skim("V00-00-03","data_DoubleElectron_PRv6");
  skim("V00-00-03","data_DoubleElectron_B30");
  skim("V00-00-03","data_DoubleElectron_B34");

  skim("V00-00-03","data_SingleMu_May10");
  skim("V00-00-03","data_SingleMu_PRv4");
  skim("V00-00-03","data_SingleMu_Aug05");
  skim("V00-00-03","data_SingleMu_PRv6");
  skim("V00-00-03","data_SingleMu_B30");
  skim("V00-00-03","data_SingleMu_B34");

  skim("V00-00-03","dymm_testskim");

}



void skim(char* version, char* sample){

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  char* path                    = Form("smurf/%s",version);
  char* infilename		= Form("%s/%s.root",path,sample);

  //--------------------------------------------------
  // list of output files
  //--------------------------------------------------

  char* outfilename		= Form("%s/%s_2jets.root",path,sample);    // >=2 jets

  //--------------------------------------------------
  // list of cuts definining output files
  //--------------------------------------------------
  
  char* njets2                  = "njets>=2";

  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------

  cout << endl << endl;
  cout << "Reading in : " << infilename     << endl;
  cout << "Writing to : " << outfilename    << " " << njets2   << endl;

  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("leptons");
  chain->Add(infilename);

  //--------------------------------------------------------------
  // njets2
  //--------------------------------------------------------------

  TFile *file_njets2 = TFile::Open(outfilename, "RECREATE");
  assert( file_njets2 != 0 );
  TTree* tree_njets2 = chain->CopyTree( njets2 );
  tree_njets2->Write();
  file_njets2->Close();

  //--------------------------------------------------------------
  // dummy check
  //--------------------------------------------------------------

  TChain *chin  = new TChain("leptons");
  TChain *chout = new TChain("leptons");

  chin->Add(infilename);
  chout->Add(outfilename);

  cout << "Infile  total  entries " << chin->GetEntries()              << endl;
  cout << "Infile  njets2 entries " << chin->GetEntries(TCut(njets2))  << endl;
  cout << "Outfile total  entries " << chout->GetEntries()             << endl;
  cout << "Outfile njets2 entries " << chout->GetEntries(TCut(njets2)) << endl;


}
