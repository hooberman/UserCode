#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");

  char* datapath = "/tas07/disk00/cms2/SpecializedSkims";

  ch->Add( Form( "%s/%s.root" , datapath , prefix ) );
  
  bool isData         = true;
  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  looper* myLooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET);
  
}





