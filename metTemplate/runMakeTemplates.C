#include "TChain.h"
#include "makeTemplates.C"

void runMakeTemplates(char* prefix , bool isData = true,  
                      makeTemplates::selectionType sel = makeTemplates::e_QCDSelection, 
                      float kFactor = 1.){

  TChain* ch = new TChain("Events");


  //------------------------------------------------------------------------------------------------------------

  if( strcmp( prefix , "QCD_Pt15" ) == 0 ){
    ch->Add("/tas/cms2/QCD_Pt15_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "PhotonJet" ) == 0 ){
    ch->Add("/tas/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root");
    ch->Add("/tas/cms2/PhotonJet_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "QCD_Pt30" ) == 0 ){
    ch->Add("/tas/cms2/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }


  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "JetMETTau" ) == 0 ){
    ch->Add("/tas/cms2/MinimumBias_Commissioning10-SD_JetMETTau-Jun14thSkim_v1_RECO/V03-04-26-02/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/JetMETTau_Run2010A-Jun14thReReco_v2_RECO/V03-04-26-01/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/JetMETTau_Run2010A-Jul16thReReco-v1_RECO/V03-04-26-07/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/JetMETTau_Run2010A-PromptReco-v4_RECO/V03-04-26-07/pfJetPt30Skim/skimmed*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "EG" ) == 0 ){
 
    ch->Add("/tas/cms2/EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-07/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-25/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-07/pfJetPt30Skim/skimmed*root");
    ch->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/pfJetPt30Skim/skimmed*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  //------------------------------------------------------------------------------------------------------------
    
  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  makeTemplates* myLooper = new makeTemplates();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, sel, -1 ,kFactor);
  
}

