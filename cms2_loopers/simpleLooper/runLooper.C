#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if( strcmp( prefix , "zmm" ) == 0 ){
    ch->Add("/tas/cms2/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root");
  }
  
  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  looper* mylooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  mylooper->ScanChain(ch, prefix, isData);
  
}





