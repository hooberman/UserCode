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

void skim(char* path, char* cut, char* label, char* sample);

void skimSamples(){

  char* path  = "../output/V00-01-07";
  char* cut   = "( (njets>=2 || njets40>=2) && pfmet>100.0 )";
  char* label = "2jets_met100";

  skim(path,cut,label,"data_53X_baby");
  skim(path,cut,label,"data_2012C_53X_baby");
  // skim(path,cut,label,"ttbar_53X_baby");
  // skim(path,cut,label,"zjets_full_53X_baby");
  // skim(path,cut,label,"ww_53X_baby");
  // skim(path,cut,label,"wz_53X_baby");
  // skim(path,cut,label,"wz2l2q_53X_baby");
  // skim(path,cut,label,"zz_53X_baby");
  // skim(path,cut,label,"zz4l_53X_baby");
  // skim(path,cut,label,"zz2l2q_53X_baby");
  // skim(path,cut,label,"t_53X_baby");
  // skim(path,cut,label,"ttW_53X_baby");
  // skim(path,cut,label,"ttZ_53X_baby");
  // skim(path,cut,label,"VVV_53X_baby");

}



void skim(char* path, char* cut, char* label, char* sample){

  //--------------------------------------------------
  // path of input and output files
  //--------------------------------------------------

  char* infilename		= Form("%s/%s.root"   ,path,sample);
  char* outfilename		= Form("%s/%s_%s.root",path,sample,label);

  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------

  cout << endl << endl;
  cout << "Reading in : " << infilename     << endl;
  cout << "Writing to : " << outfilename    << " " << cut   << endl;

  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("T1");
  chain->Add(infilename);

  //--------------------------------------------------
  // output file and tree
  //--------------------------------------------------

  TFile *file_cut = TFile::Open(outfilename, "RECREATE");
  assert( file_cut != 0 );
  TTree* tree_cut = chain->CopyTree( cut );
  tree_cut->Write();
  file_cut->Close();

  //--------------------------------------------------
  // dummy check
  //--------------------------------------------------

  TChain *chin  = new TChain("T1");
  TChain *chout = new TChain("T1");

  chin->Add(infilename);
  chout->Add(outfilename);

  cout << "Infile  total  entries " << chin->GetEntries()           << endl;
  cout << "Infile  cut    entries " << chin->GetEntries(TCut(cut))  << endl;
  cout << "Outfile total  entries " << chout->GetEntries()          << endl;
  cout << "Outfile cut    entries " << chout->GetEntries(TCut(cut)) << endl;


}
