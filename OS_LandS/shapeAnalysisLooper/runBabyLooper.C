#include "TChain.h"
#include "babylooper.C"

void runBabyLooper(char* version, char* prefix , bool isData = false ){

  TChain* ch = new TChain("t");

  string file = Form("output/%s/highpt/%s_smallTree.root",version,prefix);
  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  babylooper* myLooper = new babylooper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, version, prefix, isData);
  
}

