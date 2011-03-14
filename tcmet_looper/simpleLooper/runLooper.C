#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if( strcmp( prefix , "dyee" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-07-03_PFCandidates/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/ntuple*root");
  }

  else if( strcmp( prefix , "dymm" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-07-03_PFCandidates/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/ntuple*root");
  }
  
  else if( strcmp( prefix , "h130" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-07-03_PFCandidates/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/ntuple*root");
  }

  else if( strcmp( prefix , "dyee_nopu" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-06-18_PFCandidates/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-START38_V12-v1/ntuple*root");
  }

  else if( strcmp( prefix , "dymm_nopu" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-06-18_PFCandidates/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-START38_V12-v1/ntuple*root");
  }
  
  else if( strcmp( prefix , "h130_nopu" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-06-18_PFCandidates/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-START38_V12-v1/ntuple*root");
  }
  
  else if( strcmp( prefix , "test384" ) == 0 ){
    ch->Add("ntuple_V03_08_04.root");
  }
  
  else if( strcmp( prefix , "test412" ) == 0 ){
    ch->Add("test_3_11_1_patch2.root");
    isData = true;
  }
  
  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  looper* mylooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  mylooper->ScanChain(ch, prefix, isData);
  
}





