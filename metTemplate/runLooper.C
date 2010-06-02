#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix , bool isData = true, bool makeMetTemplate = false){

  TChain* ch = new TChain("Events");

  if( strcmp( prefix , "PhotonJet_Pt15" ) == 0 ){
    ch->Add("/tas07/disk00/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root");
  }
  else if( strcmp( prefix , "JetMETTau" ) == 0 ){
    ch->Add("/tas07/disk00/cms2/SpecializedSkims/Commissioning10-SD_JetMETTau-v9_goodrunPfJetPt30.root");
  }
  else if( strcmp( prefix , "EG" ) == 0 ){
    ch->Add("/tas07/disk00/cms2/SpecializedSkims/Commissioning10-SD_EG-v9_goodrunPhotonGT10GeV.root");
  }

  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  looper* myLooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, makeMetTemplate);
  
}

