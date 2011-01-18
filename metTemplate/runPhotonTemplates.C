#include "TChain.h"
#include "makePhotonTemplates.C"

void runPhotonTemplates(){

  TChain* ch = new TChain("T1");

  string file = Form("photonTemplates/nov5th_v3/*_baby.root");
  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  makePhotonTemplates* myLooper = new makePhotonTemplates();
  
  myLooper->ScanChain(ch);
  
}

