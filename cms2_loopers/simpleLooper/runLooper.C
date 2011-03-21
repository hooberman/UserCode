#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if( strcmp( prefix , "zmm" ) == 0 ){
    ch->Add("/tas/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged_ntuple_1.root");
  }
  
  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  looper* mylooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  mylooper->ScanChain(ch, prefix, isData);
  
}





