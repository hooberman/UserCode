#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TString.h"
#include "TROOT.h"
#include "../Tools/goodrun.cc"


void skimTree(char* runlist, char* inputFile, char* outputFile)
{
  
  cout << endl;
  cout << "Apply goodrun list : " << runlist    << endl;
  cout << "Input file         : " << inputFile  << endl;
  cout << "Output file        : " << outputFile << endl;

  set_goodrun_file(runlist);

  TFile* fin = new TFile(inputFile);
  TTree* ch  = (TTree*)fin->Get("T1"); 

  if (ch==0x0) return; 

  TFile *newfile= new TFile(outputFile,"recreate");

  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

  //unsigned int run_ = 0;
  //unsigned int lumi_ = 0;

  Int_t  run_;
  Int_t  lumi_;

  ch->SetBranchAddress( "run"           , &run_     );     
  ch->SetBranchAddress( "lumi"          , &lumi_     );     

  for(int ievt = 0; ievt < ch->GetEntries() ;ievt++) {
    ch->GetEntry(ievt); 
            
    if (!goodrun_json(run_, lumi_)) continue;

    evt_tree->Fill();
  }

  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
}  



void applyJSON(){

  char* runlist    = "../jsons/final_19p47fb_cms2.txt";

  char* inputFile  = "../output/V00-02-02/data_53X_2012A_baby_nojson.root";
  char* outputFile = "../output/V00-02-02/data_53X_2012A_baby.root";
  skimTree( runlist , inputFile , outputFile );

  inputFile  = "../output/V00-02-02/data_53X_2012B_baby_nojson.root";
  outputFile = "../output/V00-02-02/data_53X_2012B_baby.root";
  skimTree( runlist , inputFile , outputFile );

  inputFile  = "../output/V00-02-02/data_53X_2012C_baby_nojson.root";
  outputFile = "../output/V00-02-02/data_53X_2012C_baby.root";
  skimTree( runlist , inputFile , outputFile );

  inputFile  = "../output/V00-02-02/data_53X_2012D_baby_nojson.root";
  outputFile = "../output/V00-02-02/data_53X_2012D_baby.root";
  skimTree( runlist , inputFile , outputFile );

}
