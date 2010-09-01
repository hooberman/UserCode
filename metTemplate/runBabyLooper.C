#include "TChain.h"
#include "babylooper.C"

void runBabyLooper(char* prefix , bool isData = true, babylooper::metAlgo algo = babylooper::e_makeTemplate){

  TChain* ch = new TChain("T1");

  if     ( strcmp( prefix , "PhotonJet_Pt15" ) == 0 )      ch->Add("root/PhotonJet_Pt15_baby.root");
  else if( strcmp( prefix , "JetMETTau_250nb" ) == 0 )     ch->Add("root/JetMETTau_250nb_baby.root");
  else if( strcmp( prefix , "EG" ) == 0 )                  ch->Add("root/EG_baby.root");
  else if( strcmp( prefix , "EG_1p9pb" ) == 0 )            ch->Add("root/EG_1p9pb_baby.root");
  else if( strcmp( prefix , "ZJets" ) == 0 )               ch->Add("root/ZJets_baby.root");
  else if( strcmp( prefix , "ZJets_1p9pb" ) == 0 )         ch->Add("root/ZJets_1p9pb_baby.root");
  else if( strcmp( prefix , "QCD_Pt15" ) == 0 )            ch->Add("root/QCD_Pt15_baby.root");
  else if( strcmp( prefix , "dilep" ) == 0 )               ch->Add("root/dilep_baby.root");
  else if( strcmp( prefix , "lepdata_1p9pb" ) == 0 )       ch->Add("root/lepdata_1p9pb_baby.root");
  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  babylooper* myLooper = new babylooper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, algo);
  
}

