#include "TChain.h"
#include "Z_looper.C"
//#include "Z_looper.h"

void runZLooper(char* prefix , bool isData = true, Z_looper::metAlgo algo = Z_looper::e_makeTemplate, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //------------------------------------------------------------------------------------------------------------

  if( strcmp( prefix , "testdata" ) == 0 ){
    ch->Add("/tas/benhoob/skims/cms2/ExpressPhysics_Run2010A-Express-v4_FEVT_dilepmet.root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ZJets" ) == 0 ){
    ch->Add("/tas/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "TTbar" ) == 0 ){
    ch->Add("/tas/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }


  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata" ) == 0 ){

    //electron
    ch->Add("/tas/cms2/EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/diLepPt1010skim/dilepton_skimmed_ntuple.root");
    ch->Add("/tas/cms2/EG_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-07/singleLepPt5Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-07/singleLepPt10Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed*root");
    
    //muon
    ch->Add("/tas/cms2/Mu_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/diLepPt1010skim/dilepton*root");
    ch->Add("/tas/cms2/Mu_Run2010A-Jul16thReReco-v1_RECO/V03-04-26-07/singleLepPt5Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-07/singleLepPt10Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed*root");
    
    
    //ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed_ntuple_143657_12.root");
  }


  //------------------------------------------------------------------------------------------------------------

  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  //------------------------------------------------------------------------------------------------------------
    
  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  Z_looper* myLooper = new Z_looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, algo, -1 ,kFactor);
  
}

