#include "TChain.h"
#include "makeVictorTemplates.C"

void runVictorTemplates(){

  TChain* ch = new TChain("jetTree");

  string file = Form("qcd-mini-ntuple/qcd-templates*root");
  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  makeVictorTemplates* myLooper = new makeVictorTemplates();
  
  myLooper->ScanChain(ch);
  
}

