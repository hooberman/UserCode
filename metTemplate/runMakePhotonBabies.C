#include "TChain.h"
#include "makePhotonBabies.C"

void runMakePhotonBabies(char* prefix , bool isData = true, float kFactor = 1.){

  TChain* ch = new TChain("Events");
  
  //------------------------------------------------------------------------------------------------------------
  
  if( strcmp( prefix , "PhotonJet" ) == 0 ){
    //ch->Add("/tas/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged_ntuple.root");
    ch->Add("/tas/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*root");
    ch->Add("/tas/cms2/PhotonJet_Pt30_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*root");
    ch->Add("/tas/cms2/PhotonJet_Pt80_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*root");
    ch->Add("/tas/cms2/PhotonJet_Pt170_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "EG" ) == 0 ){
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/pfJetPt30Skim/skimmed_ntuple_144112_19.root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "Photon" ) == 0 ){
    ch->Add("/tas/cms2/Photon_Run2010B-PromptReco-v2_RECO/V03-06-09/merged_ntuple_147926_0.root");
    //ch->Add("/tas/cms2/Photon_Run2010B-PromptReco-v2_RECO/V03-06-14/merged*root");
    //ch->Add("/tas/cms2/Photon_Run2010B-PromptReco-v2_RECO/V03-06-09/merged*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  //------------------------------------------------------------------------------------------------------------
    
  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  makePhotonBabies* myLooper = new makePhotonBabies();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, -1 ,kFactor);
  
}

