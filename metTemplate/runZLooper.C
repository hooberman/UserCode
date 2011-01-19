#include "TChain.h"
#include "Z_looper.C"
//#include "Z_looper.h"

void runZLooper(char* prefix , bool isData = true, Z_looper::metAlgo algo = Z_looper::e_makeTemplate, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //------------------------------------------------------------------------------------------------------------
  if( strcmp( prefix , "testdata" ) == 0 ){
    ch->Add("/tas/benhoob/skims/cms2/zjetsmet/zjetsmet_allevents_oct15th.root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ZJets" ) == 0 ){
    //ch->Add("/tas/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
    ch->Add("/tas/cms2/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/merged*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "TTbar" ) == 0 ){
    //ch->Add("/tas/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root");
    ch->Add("/tas/cms2/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/diLepPt2020Skim/skim*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4" ) == 0 ){
    //ch->Add("/tas/cms2/LM4_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*root");
    ch->Add("/tas/cms2/LM4_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata" ) == 0 ){
       
    //    ch->Add("/tas/cms2/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple_143827_7.root");
    ch->Add("/tas/cms2/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root");
    ch->Add("/tas/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed*root");
    ch->Add("/tas/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root");
    
    ch->Add("/tas/cms2/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed*root");
    ch->Add("/tas/cms2/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root");
    
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata_skim" ) == 0 ){
    ch->Add("/tas/cms2/dilepSkim35pb/els.root");
    ch->Add("/tas/cms2/dilepSkim35pb/mus.root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "WJets" ) == 0 ){
    ch->Add("/tas/cms2/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/skimmed*root");
  }
  
  else if ( strcmp( prefix , "WW" ) == 0 ){
    //ch->Add("/tas/cms2/WW_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root");
    ch->Add("/tas/cms2/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/diLepPt2020/merged*root");
  }        
  
  else if ( strcmp( prefix , "WZ" ) == 0 ){
    //ch->Add("/tas/cms2/WZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root");
    ch->Add("/tas/cms2/WZTo3LNu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root");
  } 
  
  else if ( strcmp( prefix , "ZZ" ) == 0 ){
    //ch->Add("/tas/cms2/ZZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root");
    ch->Add("/tas/cms2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/diLepPt2020/merged*root");
  }
  
  else if ( strcmp( prefix , "tW" ) == 0 ){
    //ch->Add("/tas/cms2/SingleTop_sChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root");
    //ch->Add("/tas/cms2/SingleTop_tChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root");
    //ch->Add("/tas/cms2/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root");
    ch->Add("/tas/cms2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Fall10-START38_V12-v2/V03-06-14/diLepPt2020/merged*root");
    ch->Add("/tas/cms2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Fall10-START38_V12-v2/V03-06-17/merged*root");
    ch->Add("/tas/cms2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Fall10-START38_V12-v1/V03-06-17/merged*root");
  }
  
  //------------------------------------------------------------------------------------------------------------

  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  //------------------------------------------------------------------------------------------------------------
    
  bool calculateTCMET = true;  //recalculate tcmet on-the-fly?
  
  Z_looper* myLooper = new Z_looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, algo, -1 ,kFactor);
  
}




