#include "TChain.h"
#include "makePhotonBabies.C"

void runMakePhotonBabies(char* prefix , bool isData = true, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //-----------------------------------------------------------------------------------

  if( strcmp( prefix , "Photon" ) == 0 ){
    //ch->Add("/tas/cms2/Photon_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-01-02_merged/merged_ntuple_161312_0.root");
    ch->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_1_2_patch1_V04-01-05/Photon_Run2011A-Apr22ReReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-05_merged/V04-01-05/merged*root");
  }

  //-----------------------------------------------------------------------------------

  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  //-----------------------------------------------------------------------------------
    
  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  makePhotonBabies* myLooper = new makePhotonBabies();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, -1 ,kFactor);
  
}

