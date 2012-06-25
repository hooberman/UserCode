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

    // 2012A
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/store*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012A-PromptReco-v1_AOD/unmerged/store*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/MuEG_Run2012A-PromptReco-v1_AOD/unmerged/store*root");

    // 2012 B
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012B-PromptReco-v1_AOD/unmerged/store*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012B-PromptReco-v1_AOD/unmerged/store*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/MuEG_Run2012B-PromptReco-v1_AOD/unmerged/store*root");

  }
  
  //----------------------------------------------------------------------------------------
  else if( strcmp( prefix , "dataskim" ) == 0 ){    

    // 2012A
    pickSkimIfExists(ch,"/nfs-6/userdata/benhoob/ZMet2012/DoubleElectron_Run2012A-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/benhoob/ZMet2012/DoubleMu_Run2012A-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/benhoob/ZMet2012/MuEG_Run2012A-PromptReco-v1_AOD/V05-02-27/merged*root");

    // 2012B
    pickSkimIfExists(ch,"/nfs-6/userdata/benhoob/ZMet2012/DoubleElectron_Run2012B-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/benhoob/ZMet2012/DoubleMu_Run2012B-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/benhoob/ZMet2012/MuEG_Run2012B-PromptReco-v1_AOD/V05-02-27/merged*root");
  }
  
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged*1.root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged_ntuple.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "testfilter_newJEC" ) == 0 ){
    //pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012A-PromptReco-v1_AOD/merged/merged_ntuple_193334_0.root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012A-PromptReco-v1_AOD/merged/merged_ntuple_193334_0.root");

    pickSkimIfExists(ch,"/home/users/benhoob/filters/output/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012A-PromptReco-v1_AOD/merged/ZMet2012/merged_ntuple.root");
    pickSkimIfExists(ch,"/home/users/benhoob/filters/output/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012A-PromptReco-v1_AOD/merged/ZMet2012/merged_ntuple.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ttbar" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012//Summer12MC/TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012//Summer12MC/TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged_ntuple.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v3/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v3/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v3/V05-02-27/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ww" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "RelValZEE" ) == 0 ){
    pickSkimIfExists(ch,"RelValZEE.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "RelValZMM" ) == 0 ){
    pickSkimIfExists(ch,"RelValZMM.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "t" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
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
