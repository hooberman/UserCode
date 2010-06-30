#include "TChain.h"
#include "looper.C"
//#include "looper.h"

void runLooper(char* prefix , bool isData = false){

  TChain* ch = new TChain("Events");

  if( strcmp( prefix , "zee" ) == 0 ){
    ch->Add("/tas05/disk00/benhoob/tcmetTestFiles/output/PFstudies_zmm.root");
  }

  
  looper* myLooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData);
  
}

