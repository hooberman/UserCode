#include "TChain.h"
#include "makePhotonBabies.C"


void pickSkimIfExists( TChain *ch, const std::string& base )
{

  TChain *dummy = new TChain("Events");
  int nFiles = 0;
  if (dummy->Add(base.c_str())) {
    nFiles = ch->Add(base.c_str());
    std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
              << nFiles << " files" << std::endl;
  } else{
    std::cout << "Didn't find sample " << base << " , quitting" << std::endl;
    exit(0);
  }

  // be paranoid
  if (nFiles == 0) {
    std::cout << "ERROR: expected to read files " 
              << base.c_str() << "  but found none" << std::endl;
    assert(0);
  }

  return;
}


void runMakePhotonBabies(char* prefix , bool isData = true, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //-----------------------------------------------------------------------------------

  if( strcmp( prefix , "Photon" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-4/userdata/cms2/Photon_Run2011A-Apr22ReReco-v2_AOD/V04-01-05/merged*root");
    pickSkimIfExists(ch,"/nfs-4/userdata/cms2/Photon_Run2011A-PromptReco-v1_AOD/V04-01-02/merged_160329_161312.root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/imacneill/CMSSW_4_1_2_patch1_V04-01-03/Photon_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/merged*root");
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

