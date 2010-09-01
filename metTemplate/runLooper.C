#include "TChain.h"
#include "looper.C"
//#include "looper.h"

void runLooper(char* prefix , bool isData = true, looper::metAlgo algo = looper::e_makeTemplate, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //------------------------------------------------------------------------------------------------------------

  if( strcmp( prefix , "PhotonJet" ) == 0 ){
    ch->Add("/tas/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root");
    ch->Add("/tas/cms2/PhotonJet_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "JetMETTau_250nb" ) == 0 ){
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

  else if( strcmp( prefix , "ZJets" ) == 0 ){
    ch->Add("/tas/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "TTbar" ) == 0 ){
    ch->Add("/tas/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "QCD_Pt15" ) == 0 ){
    ch->Add("/tas/cms2/QCD_Pt15_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "QCD_Pt30" ) == 0 ){
    ch->Add("/tas/cms2/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "Wenu" ) == 0 ){
    ch->Add("/tas/cms2/Wenu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root");
  }

  //------------------------------------------------------------------------------------------------------------

  //top data skims
  else if ( strcmp( prefix , "dilep" ) == 0 ){
    ch->Add("/home/users/jmuelmen/CMSSW_3_6_1_patch4/src/CMS2/NtupleMacros/NtupleTools/dilep_skim_2.root");
    ch->Add("/nfs-3/userdata/fgolf/SSskims/data/skimmed_ntuple*.root");
  
    //Muon Prompt reco files
    ch->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed_ntuple_1426*.root");
    
    //EG prompt reco files      
    ch->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed_ntuple_1426*.root");
    
    //from hcal bad runs
    ch->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139783_0.root");
    ch->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139783_0.root");
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
  
  looper* myLooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, algo, -1 ,kFactor);
  
}

