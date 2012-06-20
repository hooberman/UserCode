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
  int prescale = 10;

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
    processLeptonTree("test", SmurfTree::dymm, "/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root", false, "",prescale);
  }


  return 0; 
}
