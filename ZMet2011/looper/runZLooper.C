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

  if( strcmp( prefix , "singlemudata" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");

  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "data" ) == 0 ){

    //pickSkimIfExists(ch,"/nfs-7a/userdata/cms2/DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim/skimmed_ntuple_173438_0.root");
    //pickSkimIfExists(ch,"/nfs-7a/userdata/cms2/DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim/skimmed*root");
    //pickSkimIfExists(ch,"cms2_data/DoubleElectron_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skimmed_ntuple_999999_9_1.root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleMuTriggerSkim/skimmed_ntuple_999999_1.root");

    /*
    //---------------------------
    // May10 rereco
    //---------------------------

    pickSkimIfExists(ch,"cms2_data/DoubleElectron_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skim*root");
    pickSkimIfExists(ch,"cms2_data/DoubleMu_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skim*root");
    pickSkimIfExists(ch,"cms2_data/MuEG_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skim*root");
    
    //---------------------------
    // prompt reco v4
    //---------------------------
    
    pickSkimIfExists(ch,"cms2_data/DoubleElectron_Run2011A-PromptReco-v4_AOD/V04-02-20/DoubleElectronTriggerSkim/skim*root");
    pickSkimIfExists(ch,"cms2_data/DoubleMu_Run2011A-PromptReco-v4_AOD/V04-02-20/DoubleMuTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
    
    //---------------------------
    // august rereco
    //---------------------------
    
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleElectronTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleMuTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/MuEG_Run2011A-05Aug2011-v1_AOD/V04-02-30/SSignSkim/skim*root");
    
    //---------------------------
    // prompt reco v6
    //---------------------------
    
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleElectronTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/MuEG_Run2011A-PromptReco-v6_AOD/V04-02-30/SSignSkim/skim*root");
    
    //---------------------------
    // Run2011B prompt reco v1
    //---------------------------
    
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleElectron_Run2011B-PromptReco-v1_AOD/V04-02-30/DoubleElectronTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleMu_Run2011B-PromptReco-v1_AOD/V04-02-30/DoubleMuTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/MuEG_Run2011B-PromptReco-v1_AOD/V04-02-30/SSignSkim/skim*root");

    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleElectron_Run2011B-PromptReco-v1_AOD/V04-02-34/DoubleElectronTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/DoubleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/DoubleMuTriggerSkim/skim*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/MuEG_Run2011B-PromptReco-v1_AOD/V04-02-34/SSignSkim/skim*root");
    */
  }
  
  //----------------------------------------------------------------------------------------
  
  else if( strcmp( prefix , "ttbar" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------
  
  else if( strcmp( prefix , "ggmsb" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/PhysicsProcess_higgsino_mu_SLHA_new_v1_macneill-PhysicsProcess_higgsino_mu_SLHA_new_v1-73833d9410398befa066c3c6e1fda77f/VB04-02-29_FastSim_Rutgers/preprocessing/*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wzsms" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/macneill/CMS2_VB04-02-29_FastSim_TChizz/SMS-TChiwz_mNeutralino-100to500_mLSP-0to400_7TeV-Pythia6Z_Summer11-START42_V11-v1/*.root"); //'official'
    //pickSkimIfExists(ch,"/hadoop/cms/store/user/fgolf/CMS2_VB04-02-29_Fastsim/TChiwz/*.root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/user/fgolf/CMS2_VB04-02-29_Fastsim/TChiwz/TChiwz_400_50To425_200.lhe_50000.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zzsms" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/macneill/CMS2_VB04-02-29_FastSim_TChizz/SMS-TChizz_mNeutralino-100to500_mLSP-0to400_7TeV-Pythia6Z_Summer11-START42_V11-v1/*.root"); //'official'
    //pickSkimIfExists(ch,"/hadoop/cms/store/user/fgolf/CMS2_VB04-02-29_Fastsim/TChizz/*.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjetsS6_incomplete" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMS2_V04-02-29/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/ntuple_8*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "dyee" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "dymm" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-3/userdata/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "T5zz" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/warren/SMS-T5zz_x-05_Mgluino-150to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ZZZ" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/ZZZNoGstar_TuneZ2_7TeV-madgraphCMSSW42xPUv1_spadhi-ZZZNoGstar_TuneZ2_7TeV-madgraphCMSSW42xPUv1-9ab11d163a88ab8f3641ab081403ebc5/VB04-02-29_FastSim/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "T5zzl" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7a/userdata/cms2/SMS-T5zzl_Mgluino-150to1200_mLSP-50to1100_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "T5zzh" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/benhoob/CMS2_CMS2_VB04-02-29_Fastsim/SMS-T5zzh_mGluino-150to1200_mLSP-50to1100_7TeV-Pythia6Z_Summer11-START42_V11-v5/ntuple*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "T5zzgmsb" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/SMS-T5zzgmsb_mGluino-100to1200_mNLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "T5zzgmsb_hadoop" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/benhoob/CMS2_VB04-02-29_Fastsim/SMS-T5zzgmsb_mGluino-100to1200_mNLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/ntuple*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz_summer11_madgraph" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz_summer11_pythia" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz_spring11_pythia" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-3/userdata/cms2/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz_summer11_madgraph" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root"); 
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz_summer11_pythia" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz_spring11_pythia" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-3/userdata/cms2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "singletop" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/LM4_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4v2" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/LM4_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM8" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/LM8_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM4v2" ) == 0 ){
    pickSkimIfExists(ch,"/tas/cms2/LM4_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM8v2" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/LM8_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "LM9" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-7/userdata/cms2/LM9_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
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
