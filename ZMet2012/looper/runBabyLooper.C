#include "TChain.h"
#include "babylooper.C"

void runBabyLooper(char* Z_version, char* template_version, char* prefix , bool isData = true, 
                   babylooper::selectionType sel = babylooper::e_QCDSelection, 
                   bool makeTemplate = false){

  TChain* ch = new TChain("T1");

  string file = Form("../output/%s/%s_baby_2jets.root",Z_version,prefix);
  //string file = Form("../output/%s/%s_baby_tenPercent.root",Z_version,prefix);
  cout << "Adding " << file << endl;

  ch->Add( file.c_str() );

  babylooper* myLooper = new babylooper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, Z_version, template_version, prefix, isData, sel, makeTemplate);
  
}

