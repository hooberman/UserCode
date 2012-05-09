
#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "looper.h"

#include <iostream>
#endif

void pickSkimIfExists( TChain *ch, const std::string& base, const std::string& skimPrefix = "" )
{
  TChain *dummy = new TChain("Events");

  int nFiles = 0;
  if (dummy->Add(base.c_str())) {
    nFiles = ch->Add(base.c_str());
    std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
              << nFiles << " files" << std::endl;
  } else
    std::cout << "Couldn't find " << base << std::endl;

  // be paranoid
  if (nFiles == 0) {
    std::cout << "ERROR: expected to read files " 
              << base.c_str() << "  but found none" << std::endl;
    //assert(0);
    exit(0);
  }

  return;
}

void doAll(bool skipFWLite = true)
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-00-07";
  const char* jsonfile   = "jsons/json_DCSONLY_190389_191859_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  // Load various tools  
  gROOT->ProcessLine(Form(".x setup.C(%d)", skipFWLite));

  // Load FWLite
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  // Load and compile the looping code
  gSystem->CompileMacro("looper.C","++k", "liblooper");
  
  looper* looper = new looper();

  //set looper parameters
  looper->set_susybaseline(0);
  //make baby ntuple
  looper->set_createTree(1);
  //use bitmask selection
  looper->set_useBitMask(0);
  //set version
  looper->set_version(version);
  //set json
  looper->set_json( jsonfile );
 
  //----------------------------------------
  // flags for files to run over
  //----------------------------------------

  bool runData          = 0;
  bool runElData        = 0;
  bool runMuData        = 0;
  bool runSingleElData  = 0;
  bool runSingleMuData  = 0;
  bool runElHadData     = 1;
  bool runDoubleElData  = 0;
  bool runDoubleMuData  = 0;
  bool runMuHadData     = 0;
  bool runPhotonData    = 0;

  //----------------------------------------------------------------------------------------------------------
  // all data
  //----------------------------------------------------------------------------------------------------------

  TChain* chData = new TChain("Events");

  if( runData ){
    cout << "Adding all data" << endl;

    pickSkimIfExists(chData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");  	  	 
    pickSkimIfExists(chData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root"); 	  	 
    pickSkimIfExists(chData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/MuEG_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/SingleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	 
    pickSkimIfExists(chData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/Photon_Run2012A-PromptReco-v1_AOD/unmerged/*root");	 
    pickSkimIfExists(chData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/MuHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	  	 
    pickSkimIfExists(chData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/SingleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  //----------------------------------------------------------------------------------------------------------
  // all electron data
  //----------------------------------------------------------------------------------------------------------

  TChain* chElData = new TChain("Events");

  if( runElData ){
    cout << "Adding all electron data" << endl;
    pickSkimIfExists(chElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root"); 	  	 
    pickSkimIfExists(chElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/SingleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chElData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	     
  }

  //----------------------------------------------------------------------------------------------------------
  // all muon data
  //----------------------------------------------------------------------------------------------------------

  TChain* chMuData = new TChain("Events");

  if( runMuData ){
    cout << "Adding all muon data" << endl;
    //pickSkimIfExists(chMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/SingleMu_Run2012A-PromptReco-v1_AOD/unmerged/store_data_Run2012A_SingleMu_AOD_PromptReco-v1_000_190_663_125F6FAA-C382-E111-8D4A-003048F1110E.root");
    //pickSkimIfExists(chMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/SingleMu_Run2012A-PromptReco-v1_AOD/unmerged/store_data_Run2012A_SingleMu_AOD_PromptReco-v1_000_190_663_BE881788-D182-E111-86AB-001D09F2462D.root");
    // pickSkimIfExists(chMuData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/MuEG_Run2012A-PromptReco-v1_AOD/unmerged/*root");

    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");  	  	 
    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/MuHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	  	 
    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/SingleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  //----------------------------------------------------------------------------------------------------------
  // single electron data
  //----------------------------------------------------------------------------------------------------------

  TChain* chSingleElData = new TChain("Events");

  if( runSingleElData ){
    cout << "Adding single electron data" << endl;
    pickSkimIfExists(chSingleElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/SingleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  //----------------------------------------------------------------------------------------------------------
  // single muon data
  //----------------------------------------------------------------------------------------------------------

  TChain* chSingleMuData = new TChain("Events");

  if( runSingleMuData ){
    cout << "Adding single muon data" << endl;
    pickSkimIfExists(chSingleMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/SingleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  //----------------------------------------------------------------------------------------------------------
  // ElHad data
  //----------------------------------------------------------------------------------------------------------

  TChain* chElHadData = new TChain("Events");

  if( runElHadData ){
    cout << "Adding electronhad data" << endl;
    //pickSkimIfExists(chElHadData,"/hadoop/cms/store/user/cwelke/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/store_data_Run2012A_ElectronHad_AOD_PromptReco-v1_000_191_718_C448B8E1-028C-E111-9E99-003048F110BE.root");
    //pickSkimIfExists(chElHadData,"/hadoop/cms/store/user/cwelke/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/*191_718*root");
    pickSkimIfExists(chElHadData,"/hadoop/cms/store/user/cwelke/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/*191_830*root");
    //pickSkimIfExists(chElHadData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	     
  }

  //----------------------------------------------------------------------------------------------------------
  // doubleElectron data
  //----------------------------------------------------------------------------------------------------------

  TChain* chDoubleElData = new TChain("Events");

  if( runDoubleElData ){
    cout << "Adding double electron data" << endl;
    pickSkimIfExists(chDoubleElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  //----------------------------------------------------------------------------------------------------------
  // doubleMu data
  //----------------------------------------------------------------------------------------------------------

  TChain* chDoubleMuData = new TChain("Events");

  if( runDoubleMuData ){
    cout << "Adding double muon data" << endl;
    pickSkimIfExists(chDoubleMuData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    //pickSkimIfExists(chDoubleMuData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleMu_Run2012A-PromptReco-v1_AOD/unmerged/store_data_Run2012A_DoubleMu_AOD_PromptReco-v1_000_191_090_7C2A77D9-2087-E111-877A-003048D2BC4C.root");
  }

  //----------------------------------------------------------------------------------------------------------
  // MuHad data
  //----------------------------------------------------------------------------------------------------------

  TChain* chMuHadData = new TChain("Events");

  if( runMuHadData ){
    cout << "Adding muhad data" << endl;
    pickSkimIfExists(chMuHadData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/MuHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }


  //----------------------------------------------------------------------------------------------------------
  // photon data
  //----------------------------------------------------------------------------------------------------------

  TChain* chPhotonData = new TChain("Events");

  if( runPhotonData ){
    cout << "Adding all photon data" << endl;

    pickSkimIfExists(chPhotonData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chPhotonData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/Photon_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  //----------------------------------------------------------------------------------------------------------
  // run on samples
  //----------------------------------------------------------------------------------------------------------

  if (runData) {
    cout << "Processing all data" << endl;
    looper->ScanChain(chData,"allData");
    cout << "Done processing all data" << endl;
  }

  if (runDoubleMuData) {
    cout << "Processing double muon data" << endl;
    looper->ScanChain(chDoubleMuData,"doubleMuData");
    cout << "Done processing double muon data" << endl;
  }

  if (runSingleMuData) {
    cout << "Processing single muon data" << endl;
    looper->ScanChain(chSingleMuData,"singleMuData");
    cout << "Done processing single muon data" << endl;
  }

  if (runSingleElData) {
    cout << "Processing single electron data" << endl;
    looper->ScanChain(chSingleElData,"singleElData");
    cout << "Done processing single electron data" << endl;
  }

  if (runElHadData) {
    cout << "Processing elhad data" << endl;
    looper->ScanChain(chElHadData,"elHadData");
    cout << "Done processing elhad data" << endl;
  }

  if (runDoubleElData) {
    cout << "Processing double electron data" << endl;
    looper->ScanChain(chDoubleElData,"doubleElData");
    cout << "Done processing double electron data" << endl;
  }

  if (runMuHadData) {
    cout << "Processing muhad data" << endl;
    looper->ScanChain(chMuHadData,"muHadData");
    cout << "Done processing muhad data" << endl;
  }

  if (runElData) {
    cout << "Processing all electron data" << endl;
    looper->ScanChain(chElData,"elData");
    cout << "Done processing all electron data" << endl;
  }

  if (runMuData) {
    cout << "Processing all muon data" << endl;
    looper->ScanChain(chMuData,"muData");
    cout << "Done processing all muon data" << endl;
  }

  if (runPhotonData) {
    cout << "Processing all photon data" << endl;
    looper->ScanChain(chPhotonData,"photonData");
    cout << "Done processing all photon data" << endl;
  }

  //----------------------------------------
  // save all the histograms
  //----------------------------------------

  const char* outFile  = Form("../output/%s/singleLepton.root", version );
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(outFile);
  deleteHistos();
  
  gSystem->Exit(0);
}
