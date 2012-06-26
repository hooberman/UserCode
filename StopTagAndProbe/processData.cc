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

  else if( ijob == 2 ){
    processLeptonTree("SingleMu_May10" , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/merged*root", true, goodrunlist,prescale);
    processLeptonTree("SingleMu_PRv4"  , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-PromptReco-v4_AOD/V04-02-33/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_Aug05" , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-05Aug2011-v1_AOD/V04-02-33/merged*root"  , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_PRv6"  , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011A-PromptReco-v6_AOD/V04-02-33/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_B30"   , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/merged*root" , true, goodrunlist,prescale);
    processLeptonTree("SingleMu_B34"   , SmurfTree::data, "StopTNPSkim/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/merged*root" , true, goodrunlist,prescale);
  }

  else if( ijob == 3 ){
    cout << "Processing DYJets MC" << endl;
    processLeptonTree("test", SmurfTree::dymm, "/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*root", false, "",prescale);
  }

  else if( ijob == 4 ){
    cout << "Processing DYJets MC skim" << endl;
    processLeptonTree("testskim", SmurfTree::dymm, "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root", false, "",prescale);
  }

  return 0; 
}
