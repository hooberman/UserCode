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

  TString goodrunlist = "Cert_190456-201678_8TeV_PromptReco_Collisions12_JSON_goodruns.txt";

  cout << "Using json: " << goodrunlist << endl;

  int ijob     = atoi(argv[1]);
  int prescale = 1;

  if( ijob == 1 ){
    cout << "Processing single muon data" << endl;

    processLeptonTree("SingleMu_2012A",SmurfTree::data,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/merged_ntuple_999999_11_9.root",true, goodrunlist,prescale);

    // processLeptonTree("SingleMu_2012A",SmurfTree::data,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/merged*root",true, goodrunlist,prescale);
    // processLeptonTree("SingleMu_2012B",SmurfTree::data,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012B-13Jul2012-v1_AOD/merged/merged*root",true, goodrunlist,prescale);
    // processLeptonTree("SingleMu_2012C",SmurfTree::data,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v2_AOD/merged/merged*root",true, goodrunlist,prescale);
  }

  if( ijob == 2 ){
    cout << "Processing single electron data" << endl;
    processLeptonTree("SingleEl_2012A",SmurfTree::data,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/merged_ntuple_999999_4.root",true, goodrunlist,prescale);

    // processLeptonTree("SingleEl_2012A",SmurfTree::data,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/merged*root",true, goodrunlist,prescale);
    // processLeptonTree("SingleEl_2012B",SmurfTree::data,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/merged*root",true, goodrunlist,prescale);
    // processLeptonTree("SingleEl_2012C",SmurfTree::data,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v2_AOD/merged/merged*root",true, goodrunlist,prescale);
  }

  else if( ijob == 3 ){
    cout << "Processing DYJets MC" << endl;
    //processLeptonTree("test",SmurfTree::dymm,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/SingleOrDiLepton/merged*root",false,"",prescale);
    processLeptonTree("test",SmurfTree::dymm,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/SingleOrDiLepton/merged_ntuple_1.root",false,"",prescale);
  }

  return 0; 
}
