#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h" 
#include "wwtypes.h"
#include "doAnalysis.h"
#include "LeptonTreeMaker.h"
#include "SmurfDataTypes.h"
#include "processLeptonTree.h"
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{

  // ------------------------------------------------------------------------
  // 		For BatchSubmission
  // ------------------------------------------------------------------------

  if (argc == 3) {
    SmurfTree::DataType dataType = (SmurfTree::DataType)atoi(argv[1]);
    std::string dataFile = argv[2];
    bool realData = false;
    if (dataType == 0) realData = true;
    const std::string cms2_json_file = "./files/Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2.jmu";
    processLeptonTree("test", dataType, dataFile, realData, "");
    return 0;
  }
	
  // ------------------------------------------------------------------------
  //  	For local submission
  // ------------------------------------------------------------------------

  TString goodrunlist = "Cert_160404-180252_7TeV_mergePromptMay10Aug5_JSON_goodruns.txt";

  int ijob     = atoi(argv[1]);
  int prescale = 1;

  
  // processLeptonTree("DoubleElectron_May10", SmurfTree::data, 
  // 		    "StopTNPSkim/DoubleElectron_Run2011A-May10ReReco-v1_AOD/V04-02-20/merged_ntuple_999999_10_skim.root", true, goodrunlist,prescale);

  if( ijob == 1 ){
    processLeptonTree("DoubleElectron_May10" , SmurfTree::data, "StopTNPSkim/DoubleElectron_Run2011A-May10ReReco-v1_AOD/V04-02-20/merged*root", true, goodrunlist,prescale);
    processLeptonTree("DoubleElectron_PRv4"  , SmurfTree::data, "StopTNPSkim/DoubleElectron_Run2011A-PromptReco-v4_AOD/V04-02-20/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("DoubleElectron_Aug05" , SmurfTree::data, "StopTNPSkim/DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/merged*root"  , true, goodrunlist,prescale);
    processLeptonTree("DoubleElectron_PRv6"  , SmurfTree::data, "StopTNPSkim/DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-20/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("DoubleElectron_B30"   , SmurfTree::data, "StopTNPSkim/DoubleElectron_Run2011B-PromptReco-v1_AOD/V04-02-30/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("DoubleElectron_B34"   , SmurfTree::data, "StopTNPSkim/DoubleElectron_Run2011B-PromptReco-v1_AOD/V04-02-34/merged*root" , true, goodrunlist,prescale);
  }

  else if( ijob == 1 ){
    processLeptonTree("SingleMu_May10" , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-20/merged*root", true, goodrunlist,prescale);
    processLeptonTree("SingleMu_PRv4"  , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-PromptReco-v4_AOD/V04-02-20/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_Aug05" , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/merged*root"  , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_PRv6"  , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-PromptReco-v6_AOD/V04-02-20/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_B30"   , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-30/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_B34"   , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/merged*root" , true, goodrunlist,prescale);
  }




  /*
  if( ijob == 1 ){
    cout << "Processing May10 data" << endl;
    processLeptonTree("May10", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root", true, goodrunlist,prescale);
  }
  
  else if( ijob == 2 ){
    cout << "Processing PRv4 data" << endl;
    processLeptonTree("PRv4", SmurfTree::data, "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 3 ){
    cout << "Processing Aug05 data" << endl;
    processLeptonTree("Aug05", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 4 ){
    cout << "Processing PRv6 data" << endl;
    processLeptonTree("PRv6", SmurfTree::data, "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 5 ){
    cout << "Processing 2011B-V33 data" << endl;
    processLeptonTree("2011B-V33", SmurfTree::data, "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 6 ){
    cout << "Processing 2011B-V34 data" << endl;
    processLeptonTree("2011B-V34", SmurfTree::data, "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 7 ){
    cout << "Processing DYJets MC" << endl;
    processLeptonTree("test", SmurfTree::dymm, "/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*root", false, "",prescale);
  }

  else if( ijob == 8 ){
    cout << "Processing May10 data skim" << endl;
    processLeptonTree("May10skim", SmurfTree::data, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root", true, goodrunlist,prescale);
  }
  
  else if( ijob == 9 ){
    cout << "Processing PRv4 data skim" << endl;
    processLeptonTree("PRv4skim", SmurfTree::data, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 10 ){
    cout << "Processing Aug05 data skim" << endl;
    processLeptonTree("Aug05skim", SmurfTree::data, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 11 ){
    cout << "Processing PRv6 data skim" << endl;
    processLeptonTree("PRv6skim", SmurfTree::data, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 12 ){
    cout << "Processing 2011B-V33 data skim" << endl;
    processLeptonTree("2011B-V33skim", SmurfTree::data, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 13 ){
    cout << "Processing 2011B-V34 data skim" << endl;
    processLeptonTree("2011B-V34skim", SmurfTree::data, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/SingleLeptonAndTwoJets/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 14 ){
    cout << "Processing DYJets MC skim" << endl;
    processLeptonTree("testskim", SmurfTree::dymm, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root", false, "",prescale);
  }

  else if( ijob == 15 ){
    cout << "Processing May10 DoubleElectron" << endl;
    processLeptonTree("DoubleElectron_May10", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleElectron_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 16 ){
    cout << "Processing PRv4 DoubleElectron" << endl;
    processLeptonTree("DoubleElectron_PRv4", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleElectron_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 17 ){
    cout << "Processing Aug05 DoubleElectron" << endl;
    processLeptonTree("DoubleElectron_Aug05", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 18 ){
    cout << "Processing PRv6 DoubleElectron" << endl;
    processLeptonTree("DoubleElectron_PRv6", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 19 ){
    cout << "Processing 2011B V30 DoubleElectron" << endl;
    processLeptonTree("DoubleElectron_2011B30", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root", true, goodrunlist,prescale);
  }

  else if( ijob == 20 ){
    cout << "Processing 2011B V34 DoubleElectron" << endl;
    processLeptonTree("DoubleElectron_2011B34", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/DoubleElectron_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged_ntuple_178365_0.root", true, goodrunlist,prescale);
  }
  */

  return 0; 
}
