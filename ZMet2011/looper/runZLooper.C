#include "TChain.h"
#include "Z_looper.C"
//#include "Z_looper.h"

void runZLooper(char* prefix , bool isData = true, Z_looper::metAlgo algo = Z_looper::e_makeTemplate, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //char* datapath = "/nfs-3/userdata/cms2";
  char* datapath = "/tas/cms2";

  //------------------------------------------------------------------------------------------------------------
  if( strcmp( prefix , "testdata" ) == 0 ){
    cout << "Sample not present on UAF!" << endl;
    exit(0);
    //ch->Add("/tas/benhoob/skims/cms2/zjetsmet/zjetsmet_allevents_oct15th.root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "dyee_spring11" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "dymm_spring11" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "tt_spring11" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ZJets" ) == 0 ){
    //ch->Add("Form("%s/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root",datapath));
    ch->Add(Form("%s/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/merged*root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "DY" ) == 0 ){
    ch->Add(Form("%s/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/merged*root",datapath));
    ch->Add(Form("%s/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "TTbar" ) == 0 ){
    //ch->Add("Form("%s/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root",datapath));
    ch->Add(Form("%s/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/diLepPt2020Skim/skim*root",datapath));
  }
  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM0" ) == 0 ){
    //ch->Add(Form("%s/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM1" ) == 0 ){
    //ch->Add(Form("%s/LM1_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM1_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM2" ) == 0 ){
    //ch->Add(Form("%s/LM2_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM2_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM3" ) == 0 ){
    //ch->Add(Form("%s/LM3_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM3_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4" ) == 0 ){
    //ch->Add(Form("%s/LM4_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM4_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM8" ) == 0 ){
    //ch->Add(Form("%s/LM8_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM8_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM9" ) == 0 ){
    //ch->Add(Form("%s/LM9_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM9_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata" ) == 0 ){
       
    //    ch->Add("Form("%s/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple_143827_7.root",datapath));
    ch->Add(Form("%s/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    
    ch->Add(Form("%s/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata_skim" ) == 0 ){
    //ch->Add(Form("%s/dilepSkim35pb/Sep17_promptv2/els.root",datapath));
    //ch->Add(Form("%s/dilepSkim35pb/Sep17_promptv2/mus.root",datapath));
    ch->Add(Form("%s/dilepSkim35pb/els.root",datapath));
    ch->Add(Form("%s/dilepSkim35pb/mus.root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata_skim_nov4" ) == 0 ){
    ch->Add(Form("%s/dilepSkim35pb_Nov4/els.root",datapath));
    ch->Add(Form("%s/dilepSkim35pb_Nov4/mus.root",datapath));
  }

  //------------------------------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "WJets" ) == 0 ){
    ch->Add(Form("%s/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/skimmed*root",datapath));
  }
  
  else if ( strcmp( prefix , "WW" ) == 0 ){
    //ch->Add(Form("%s/WW_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root",datapath));
    ch->Add(Form("%s/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  }        
  
  else if ( strcmp( prefix , "WZ" ) == 0 ){
    //ch->Add(Form("%s/WZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root",datapath));
    ch->Add(Form("%s/WZTo3LNu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  } 
  
  else if ( strcmp( prefix , "ZZ" ) == 0 ){
    //ch->Add(Form("%s/ZZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root",datapath));
    ch->Add(Form("%s/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  }
  
  else if ( strcmp( prefix , "tW" ) == 0 ){
    //ch->Add(Form("%s/SingleTop_sChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",datapath));
    //ch->Add(Form("%s/SingleTop_tChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",datapath));
    //ch->Add(Form("%s/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",datapath));
    ch->Add(Form("%s/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Fall10-START38_V12-v2/V03-06-14/diLepPt2020/merged*root",datapath));
    ch->Add(Form("%s/TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Fall10-START38_V12-v2/V03-06-17/merged*root",datapath));
    ch->Add(Form("%s/TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Fall10-START38_V12-v1/V03-06-17/merged*root",datapath));
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




