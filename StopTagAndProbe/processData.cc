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
  bool dodata = true;
  bool domc   = false;

  // 
  // Data
  //
  if(dodata){
    // TString goodrunlist = "Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON_cms2.txt";
    // processLeptonTree("test", SmurfTree::data, "/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/SingleElectron_Run2012B-PromptReco-v1_AOD/merged/merged_ntuple_195552_0.root", true, goodrunlist);
    
// /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/  
// /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/  
// /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/  
// /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/  
// /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ 
// /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/   

    TString goodrunlist = "Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON_cms2.txt";
    processLeptonTree("test", SmurfTree::data, "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged_ntuple_180072_0.root", true, goodrunlist);
  }

  // 
  // MC
  // 
  if(domc){
    processLeptonTree("test", SmurfTree::dyee, "/nfs-7/userdata/cms2/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/merged_ntuple_29*.root", false, "");
  }

  return 0; 
}
