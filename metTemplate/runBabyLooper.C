#include "TChain.h"
#include "babylooper.C"

void runBabyLooper(char* iter, char* prefix , bool isData = true, 
                   babylooper::selectionType sel = babylooper::e_QCDSelection, 
                   bool makeTemplate = false){

  TChain* ch = new TChain("T1");

  string file = Form("output/%s/%s_baby.root",iter,prefix);
  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  babylooper* myLooper = new babylooper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, sel, makeTemplate);
  
}

