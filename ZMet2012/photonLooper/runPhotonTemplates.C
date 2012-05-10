#include "TChain.h"
#include "makePhotonTemplates.C"

void runPhotonTemplates( char* iter ){

  TChain* ch = new TChain("T1");

  string file = Form("../templates/%s/Photon_baby.root",iter);

  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  makePhotonTemplates* myLooper = new makePhotonTemplates();
  
  myLooper->ScanChain( ch , iter );
  
}

