#include "TChain.h"
#include "Z_looper.C"
//#include "Z_looper.h"


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

void runZLooper(char* prefix , bool isData = true, float kFactor = 1.){

  TChain* ch = new TChain("Events");

  //----------------------------------------------------------------------------------------

  if( strcmp( prefix , "data" ) == 0 ){

    pickSkimIfExists(ch,"cms2_data/DoubleElectron_Run2011A-May10ReReco-v1_AOD/V04-02-15/DoubleElectronTriggerSkim/skim*root");
    pickSkimIfExists(ch,"cms2_data/DoubleMu_Run2011A-May10ReReco-v1_AOD/V04-02-15/DoubleMuTriggerSkim/skim*root");
    pickSkimIfExists(ch,"cms2_data/MuEG_Run2011A-May10ReReco-v1_AOD/V04-02-15/merged*root");

  }
  
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ttbar" ) == 0 ){
    pickSkimIfExists(ch,"cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets" ) == 0 ){
    pickSkimIfExists(ch,"cms2/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4" ) == 0 ){
    pickSkimIfExists(ch,"cms2/LM4_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*root");
  }
  
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM8" ) == 0 ){
    pickSkimIfExists(ch,"cms2/LM8_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*root");
  }
  
  //----------------------------------------------------------------------------------------

  else{
    cout << "ERROR: cannot find sample " << prefix << endl;
    exit(0);
  }

  //----------------------------------------------------------------------------------------
    
  bool calculateTCMET = false;  //recalculate tcmet on-the-fly?
  
  Z_looper* myLooper = new Z_looper();
  
  cout << "Running on sample " << prefix << endl;
  myLooper->ScanChain(ch, prefix, isData, calculateTCMET, -1 ,kFactor);
  
}






































  /*
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "dyee_spring11" ) == 0 ){
    ch->Add("/tas/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "dymm_spring11" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
  }
  else if( strcmp( prefix , "DY" ) == 0 ){
    ch->Add(Form("%s/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/merged*root",datapath));
    ch->Add(Form("%s/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "TTbar" ) == 0 ){
    //ch->Add("Form("%s/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*root",datapath));
    ch->Add(Form("%s/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/diLepPt2020Skim/skim*root",datapath));
  }
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM0" ) == 0 ){
    //ch->Add(Form("%s/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM1" ) == 0 ){
    //ch->Add(Form("%s/LM1_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM1_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM2" ) == 0 ){
    //ch->Add(Form("%s/LM2_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM2_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM3" ) == 0 ){
    //ch->Add(Form("%s/LM3_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM3_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4" ) == 0 ){
    //ch->Add(Form("%s/LM4_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM4_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM8" ) == 0 ){
    //ch->Add(Form("%s/LM8_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM8_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM9" ) == 0 ){
    //ch->Add(Form("%s/LM9_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
    ch->Add(Form("%s/LM9_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*root",datapath));
  }

  //----------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata" ) == 0 ){
       
    //    ch->Add("Form("%s/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple_143827_7.root",datapath));
    ch->Add(Form("%s/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    
    ch->Add(Form("%s/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed*root",datapath));
    ch->Add(Form("%s/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed*root",datapath));
    
  }

  //----------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata_skim" ) == 0 ){
    //ch->Add(Form("%s/dilepSkim35pb/Sep17_promptv2/els.root",datapath));
    //ch->Add(Form("%s/dilepSkim35pb/Sep17_promptv2/mus.root",datapath));
    ch->Add(Form("%s/dilepSkim35pb/els.root",datapath));
    ch->Add(Form("%s/dilepSkim35pb/mus.root",datapath));
  }

  //----------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "lepdata_skim_nov4" ) == 0 ){
    ch->Add(Form("%s/dilepSkim35pb_Nov4/els.root",datapath));
    ch->Add(Form("%s/dilepSkim35pb_Nov4/mus.root",datapath));
  }

  //----------------------------------------------------------------------------------------

  else if ( strcmp( prefix , "WJets" ) == 0 ){
    ch->Add(Form("%s/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/skimmed*root",datapath));
  }

  //----------------------------------------------------------------------------------------
  
  else if ( strcmp( prefix , "WW" ) == 0 ){
    //ch->Add(Form("%s/WW_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root",datapath));
    ch->Add(Form("%s/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  }        

  //----------------------------------------------------------------------------------------  

  else if ( strcmp( prefix , "WZ" ) == 0 ){
    //ch->Add(Form("%s/WZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root",datapath));
    ch->Add(Form("%s/WZTo3LNu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  } 

  //----------------------------------------------------------------------------------------
  
  else if ( strcmp( prefix , "ZZ" ) == 0 ){
    //ch->Add(Form("%s/ZZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root",datapath));
    ch->Add(Form("%s/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1/V03-06-14/merged*root",datapath));
  }

  //----------------------------------------------------------------------------------------
  
  else if ( strcmp( prefix , "tW" ) == 0 ){
    //ch->Add(Form("%s/SingleTop_sChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",datapath));
    //ch->Add(Form("%s/SingleTop_tChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",datapath));
    //ch->Add(Form("%s/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",datapath));
    ch->Add(Form("%s/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Fall10-START38_V12-v2/V03-06-14/diLepPt2020/merged*root",datapath));
    ch->Add(Form("%s/TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Fall10-START38_V12-v2/V03-06-17/merged*root",datapath));
    ch->Add(Form("%s/TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Fall10-START38_V12-v1/V03-06-17/merged*root",datapath));
  }
  
  //----------------------------------------------------------------------------------------
  */
