#include "TChain.h"
#include "makePhotonTemplates.C"

void runPhotonTemplates( char* iter , char* sample ){

  TChain* ch = new TChain("T1");

  string file = Form("../photon_output/%s/%s_baby_2jets.root",iter,sample);

  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  makePhotonTemplates* myLooper = new makePhotonTemplates();
  
  myLooper->ScanChain( ch , iter , sample );
  
}

