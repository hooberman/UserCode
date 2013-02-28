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

void skim(char* path, char* cut, char* label, char* sample, char* treename);

void skimSamples(){

  char* path     = "../output/V00-02-13";
  char* cut      = "mg==200 && ml==0"; 
  char* label    = "mg200_ml0";
  char* treename = "T1";

  skim(path,cut,label,"wzsms_baby_oldIso"  ,treename);

  /*
  //char* path     = "../output/V00-02-13";
  //char* cut      = "( njets>=2 || njetsup>=2 || njetsdn>=2 )"; 

  char* path     = "../photon_output/V00-02-03";
  char* cut      = "( njets>=2 )"; 
  char* label    = "2jets";
  char* treename = "T1";

  // skim(path,cut,label,"data_53X_2012A_baby"  ,treename);
  // skim(path,cut,label,"data_53X_2012B_baby"  ,treename);
  // skim(path,cut,label,"data_53X_2012C_baby"  ,treename);
  // skim(path,cut,label,"data_53X_2012D_baby"  ,treename);

  // skim(path,cut,label,"zjets_53X_slim_baby"        , treename);
  // skim(path,cut,label,"zjets_small_53X_slim_baby"  , treename);
  // skim(path,cut,label,"ttbar_53X_slim_baby"        , treename);
  // skim(path,cut,label,"zz2l2q_53X_slim_baby"       , treename);
  // skim(path,cut,label,"zz2l2nu_53X_slim_baby"      , treename);
  // skim(path,cut,label,"zz4l_53X_slim_baby"         , treename);
  // skim(path,cut,label,"t_53X_slim_baby"            , treename);
  // skim(path,cut,label,"ttw_53X_slim_baby"          , treename);
  // skim(path,cut,label,"ttz_53X_slim_baby"          , treename);
  // skim(path,cut,label,"tbz_53X_slim_baby"          , treename);
  // skim(path,cut,label,"vvv_53X_slim_baby"          , treename);
  // skim(path,cut,label,"ww_53X_slim_baby"           , treename);
  // skim(path,cut,label,"wz3lnu_53X_slim_baby"       , treename);
  // skim(path,cut,label,"wz2l2q_53X_slim_baby"       , treename);

  // skim(path,cut,label,"data_53X_2012ALL_baby"      , treename);
  skim(path,cut,label,"data_53X_2012A_baby"        , treename);
  skim(path,cut,label,"data_53X_2012B_baby"        , treename);
  skim(path,cut,label,"data_53X_2012C_baby"        , treename);
  skim(path,cut,label,"data_53X_2012D_baby"        , treename);
  */

  // char* path     = "../photon_output/V00-02-00";
  // char* cut      = "( njets>=2 )"; 
  // char* label    = "2jets";
  // char* treename = "T1";

  // skim(path,cut,label,"data_53X_2012A_baby"  ,treename);
  // skim(path,cut,label,"data_53X_2012B_baby"  ,treename);
  // skim(path,cut,label,"data_53X_2012C_baby"  ,treename);
  // skim(path,cut,label,"data_53X_2012D_baby"  ,treename);

  /*
  char* path     = "../photon_output/V00-02-03";
  char* cut      = "(sqrt(  pow(jet1->mass(),2) + pow(jet2->mass(),2) + 2 * jet1->E() * jet2->E() - 2 * jet1->Px() * jet2->Px() - 2 * jet1->Py() * jet2->Py() - 2 * jet1->Pz() * jet2->Pz() ) > 70.0 && sqrt(  pow(jet1->mass(),2) + pow(jet2->mass(),2) + 2 * jet1->E() * jet2->E() - 2 * jet1->Px() * jet2->Px() - 2 * jet1->Py() * jet2->Py() - 2 * jet1->Pz() * jet2->Pz() ) < 110.0 )";
  char* label    = "mjj70to110";
  char* treename = "T1";

  skim(path,cut,label,"data_53X_2012ALL_baby_2jets"  ,treename);
  */

  // char* path     = "../output/V00-01-07";
  // char* cut      = "( (njets>=2 || njets40>=2) && pfmet>100.0 )";
  // char* label    = "2jets_met100";
  // char* treename = "T1";

  // char* path     = "/tas/benhoob/home/StopTagAndProbe2012/smurf";
  // //char* cut      = "(njets>=2)";
  // //char* label    = "2jets";
  // char* cut      = "(probe->pt()>100.0)";
  // char* label    = "probept100";
  // char* treename = "leptons";

  // char* path     = "/tas/benhoob/testFiles/T2tt_8TeV";
  // char* cut      = "(t1metphicorrmtdn>120.0||t1metphicorrmt>120.0||t1metphicorrmtup>120.0) && (t1metphicorrdn>100.0||t1metphicorr>100.0||t1metphicorrup>100.0) && (njetsDown>=4||njetsUp>=4)";
  // char* label    = "met100_mt120_njets4";
  // char* treename = "t";

  // char* path     = "/tas/benhoob/testFiles/T2tt_8TeV";
  // char* cut      = "mg==250 && ml==0";
  // char* label    = "250_0";
  // char* treename = "t";

  // skim(path,cut,label,"merged"  ,treename);
  // skim(path,cut,label,"merged_1",treename);
  // skim(path,cut,label,"merged_2",treename);
  // skim(path,cut,label,"merged_3",treename);
  // skim(path,cut,label,"merged_4",treename);
  // skim(path,cut,label,"merged_5",treename);

  // skim(path,cut,label,"SingleMu2012A_V00-00-03/merged_json",treename);
  // skim(path,cut,label,"SingleMu2012B_V00-00-03/merged_json",treename);
  // skim(path,cut,label,"SingleMu2012C_V00-00-03/merged_json",treename);
  // skim(path,cut,label,"SingleEl2012A_V00-00-03/merged_json",treename);
  // skim(path,cut,label,"SingleEl2012B_V00-00-03/merged_json",treename);
  // skim(path,cut,label,"SingleEl2012C_V00-00-03/merged_json",treename);
  // skim(path,cut,label,"ZJets_V00-00-03/merged"        ,treename);

  // skim(path,cut,label,"data_2012C_53X_baby");
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



void skim(char* path, char* cut, char* label, char* sample, char* treename){

  //--------------------------------------------------
  // path of input and output files
  //--------------------------------------------------

  char* infilename		= Form("%s/%s.root"   ,path,sample);
  char* outfilename		= Form("%s/%s_%s.root",path,sample,label);

  // char* infilename		= Form("%s/%s*njets4.root"   ,path,sample);
  // char* outfilename		= Form("%s/%s_%s.root",path,sample,label);

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

  TChain *chain = new TChain(treename);
  chain->Add(infilename);
  // cout << "Adding all unskimmed files" << endl;
  // chain->Add( Form("%s/merged.root",path) );
  // chain->Add( Form("%s/merged_1.root",path) );
  // chain->Add( Form("%s/merged_2.root",path) );
  // chain->Add( Form("%s/merged_3.root",path) );
  // chain->Add( Form("%s/merged_4.root",path) );
  // chain->Add( Form("%s/merged_5.root",path) );

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

  TChain *chin  = new TChain(treename);
  TChain *chout = new TChain(treename);

  chin->Add(infilename);
  chout->Add(outfilename);

  cout << "Infile  total  entries " << chin->GetEntries()           << endl;
  cout << "Infile  cut    entries " << chin->GetEntries(TCut(cut))  << endl;
  cout << "Outfile total  entries " << chout->GetEntries()          << endl;
  cout << "Outfile cut    entries " << chout->GetEntries(TCut(cut)) << endl;


}
