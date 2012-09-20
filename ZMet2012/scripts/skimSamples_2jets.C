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

  skim("V00-01-06","data_53X");
  // skim("V00-01-04","data_2012C_53X");
  // skim("V00-01-04","ttbar_53X");
  // skim("V00-01-04","zjets_full_53X");
  // skim("V00-01-03","zjets_53X");
  // skim("V00-01-04","wz_53X");
  // skim("V00-01-04","wz2l2q_53X");
  // skim("V00-01-04","zz_53X");
  // skim("V00-01-04","zz4l_53X");
  // skim("V00-01-04","zz2l2q_53X");
  // skim("V00-01-04","ww_53X");
  // skim("V00-01-04","t_53X");
  // skim("V00-01-04","wjets_53X");
  // skim("V00-01-04","ttW_53X");
  // skim("V00-01-04","ttZ_53X");
  // skim("V00-01-04","VVV_53X");

  // skim("V00-01-00","DoubleElectron");
  // skim("V00-01-00","DoubleElectron_2012Cv2");

}



void skim(char* version, char* sample){

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  //char* path                    = Form("../output/%s",version);
  char* path                    = Form("../output/%s",version);
  char* infilename		= Form("%s/%s_baby.root",path,sample);

  //--------------------------------------------------
  // list of output files
  //--------------------------------------------------

  char* outfilename		= Form("%s/%s_baby_2jets.root",path,sample);    // >=2 jets

  //--------------------------------------------------
  // list of cuts definining output files
  //--------------------------------------------------
  
  char* njets2                  = "njets>=2 || njets40>=2";

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

  TChain *chain = new TChain("T1");
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

  TChain *chin  = new TChain("T1");
  TChain *chout = new TChain("T1");

  chin->Add(infilename);
  chout->Add(outfilename);

  cout << "Infile  total  entries " << chin->GetEntries()              << endl;
  cout << "Infile  njets2 entries " << chin->GetEntries(TCut(njets2))  << endl;
  cout << "Outfile total  entries " << chout->GetEntries()             << endl;
  cout << "Outfile njets2 entries " << chout->GetEntries(TCut(njets2)) << endl;


}
