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
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012A-PromptReco-v1_AOD/merged/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012A-PromptReco-v1_AOD/merged/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/MuEG_Run2012A-PromptReco-v1_AOD/merged/merged*root");

    // 2012 B
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012B-PromptReco-v1_AOD/merged/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012B-PromptReco-v1_AOD/merged/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/MuEG_Run2012B-PromptReco-v1_AOD/merged/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "data_53X" ) == 0 ){    

    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012A-13Jul2012-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012A-13Jul2012-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012A-13Jul2012-v1_AOD/V05-03-13/merged*root");

    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012A-recover-06Aug2012-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012A-recover-06Aug2012-v1_AOD/V05-03-13/merged*root");

    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012B-13Jul2012-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012B-13Jul2012-v4_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012B-13Jul2012-v1_AOD/V05-03-13/merged*root");

  }
  
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "data_2012C_53X" ) == 0 ){    

    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012C-PromptReco-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012C-PromptReco-v1_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012C-PromptReco-v1_AOD/V05-03-13/merged*root");

    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012C-PromptReco-v2_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012C-PromptReco-v2_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012C-PromptReco-v2_AOD/V05-03-13/merged*root");

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

  else if( strcmp( prefix , "dataskim2010" ) == 0 ){    

    // 2012A
    //pickSkimIfExists(ch,"ZMet2012_pt2010/DoubleElectron_Run2012A-PromptReco-v1_AOD/V05-02-27/merged_ntuple_193621_0_skim.root");

    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012A-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012A-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012A-PromptReco-v1_AOD/V05-02-27/merged*root");
    
    // 2012B
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012B-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012B-PromptReco-v1_AOD/V05-02-27/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012B-PromptReco-v1_AOD/V05-02-27/merged*root");
  }
  
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "data2012c" ) == 0 ){    

    // pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-09/DoubleElectron_Run2012C-PromptReco-v2_AOD/merged/merged*root");
    // pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-09/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/merged*root");
    // pickSkimIfExists(ch,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-09/MuEG_Run2012C-PromptReco-v2_AOD/merged/merged*root");

    // 2012B
    pickSkimIfExists(ch,"ZMet2012_pt2010/DoubleElectron_Run2012C-PromptReco-v2_AOD/V05-03-09/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010/DoubleMu_Run2012C-PromptReco-v2_AOD/V05-03-09/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010/MuEG_Run2012C-PromptReco-v2_AOD/V05-03-09/merged*root");
  }
  
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets" ) == 0 ){
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged*1.root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/SingleOrDiLepton/merged_ntuple_1*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "t_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz2l2q_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz2l2q_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wjets_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/benhoob/ZMet2010_pt2010/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets_full_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/SingleOrDiLepton/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets_MET50_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/user/benhoob/MET50/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "data2012cv2" ) == 0 ){    

    // 2012C v2
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleElectron_Run2012C-PromptReco-v2_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/DoubleMu_Run2012C-PromptReco-v2_AOD/V05-03-13/merged*root");
    pickSkimIfExists(ch,"ZMet2012_pt2010_nfs/MuEG_Run2012C-PromptReco-v2_AOD/V05-03-13/merged*root");
  }
   
  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zjets_10to50" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-10To50filter_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
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
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged_ntuple.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ttbar_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ttbar_massiveb" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S6_START52_V9-v1/V05-02-28/SingleOrDiLepton/merged*root");
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
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged_ntuple.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/nfs-6/userdata/cms2/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wz_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "zz4l_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ttW_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "ttZ_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "VVV_53X" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZZNoGstarJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWZNoGstarJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "wzsms" ) == 0 ){
    pickSkimIfExists(ch,"/nfs-3/userdata/cms2/SMS-TChiwz_Mchargino-100to500_mLSP-0to400_8TeV-Pythia6Z_Summer12-START52_V9_FSIM-v1/V05-03-13/merged*root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/SMS-TChiwz_Mchargino-100to500_mLSP-0to400_8TeV-Pythia6Z_StoreResults-V2-PU_START52_V9_FastSim-v1/V05-02-28/merged*root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/SMS-TChiwz_Mchargino-100to500_mLSP-0to400_8TeV-Pythia6Z_StoreResults-V2-PU_START52_V9_FastSim-v1/V05-02-28/merged_ntuple_15.root");
  }

  //----------------------------------------------------------------------------------------
  
  else if( strcmp( prefix , "gmsb" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Higgsino_StoreResults-Higgsino-ebe1b6443ab75fb2cd06c2581e9a7621/V05-02-28/merged*root");
    //pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Higgsino_StoreResults-Higgsino-ebe1b6443ab75fb2cd06c2581e9a7621/V05-02-28/merged_ntuple.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "RelValZEE" ) == 0 ){
    pickSkimIfExists(ch,"/tas/benhoob/home/ntupling/CMSSW_5_3_2_patch4/src/CMS2/NtupleMaker/test/RelValZEE_53X.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "RelValZMM" ) == 0 ){
    pickSkimIfExists(ch,"/tas/benhoob/home/ntupling/CMSSW_5_3_2_patch4/src/CMS2/NtupleMaker/test/RelValZMM_53X.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "t" ) == 0 ){
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
    pickSkimIfExists(ch,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "DoubleElectron_199752" ) == 0 ){
    pickSkimIfExists(ch,"/tas/benhoob/testFiles/DoubleElectron_Run2012C-PromptReco-v2/V05-03-13/run199752/DoubleElectron_53X_199752.root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "DoubleMu_199752" ) == 0 ){
    pickSkimIfExists(ch,"/tas/benhoob/testFiles/DoubleMu_Run2012C-PromptReco-v2/V05-03-13/run199752/ntuple_*root");
  }

  //----------------------------------------------------------------------------------------

  else if( strcmp( prefix , "MuEG_199752" ) == 0 ){
    pickSkimIfExists(ch,"/tas/benhoob/testFiles/MuEG_Run2012C-PromptReco-v2/V05-03-13/run199752/ntuple_*root");
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
