
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
  
  const char* version    = "V00-00-00";
  const char* jsonfile   = "jsons/json_DCSONLY_Apr10_goodruns.txt";
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
  bool runElData        = 1;
  bool runMuData        = 1;
  bool runPhotonData    = 0;

  //----------------------------------------
  // add samples to TChains
  //----------------------------------------

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

  TChain* chElData = new TChain("Events");

  if( runElData ){
    cout << "Adding all electron data" << endl;

    pickSkimIfExists(chElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root"); 	  	 
    pickSkimIfExists(chElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/MuEG_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chElData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/SingleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chElData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/ElectronHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	 
  }

  TChain* chMuData = new TChain("Events");

  if( runMuData ){
    cout << "Adding all muon data" << endl;

    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");  	  	 
    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/MuEG_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/MuHad_Run2012A-PromptReco-v1_AOD/unmerged/*root");	  	 
    pickSkimIfExists(chMuData,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch3_V05-02-07/SingleMu_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }

  TChain* chPhotonData = new TChain("Events");

  if( runPhotonData ){
    cout << "Adding all photon data" << endl;

    pickSkimIfExists(chPhotonData,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/unmerged/*root");
    pickSkimIfExists(chPhotonData,"/hadoop/cms/store/user/macneill/CMSSW_5_2_3_patch3_V05-02-07/Photon_Run2012A-PromptReco-v1_AOD/unmerged/*root");
  }


  //----------------------------------------
  // run on samples
  //----------------------------------------

  if (runData) {
    cout << "Processing all data" << endl;
    looper->ScanChain(chData,"allData");
    cout << "Done processing all data" << endl;
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
