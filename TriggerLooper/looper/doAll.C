
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
  const char* jsonfile   = "jsons/Cert_160404-180252_7TeV_mergePromptMay10Aug5_JSON_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  //Load CORE stuff
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/utilities.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/SimpleFakeRate.cc+");
  gROOT->ProcessLine(".L ../CORE/mcSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MT2/MT2.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");  
  gROOT->ProcessLine(".L ../CORE/susySelections.cc+");
  gROOT->ProcessLine(".L ../CORE/mcSUSYkfactor.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerSuperModel.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
  //gROOT->ProcessLine(".L ../CORE/jetSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");

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

  bool runElHad     = 1;
  bool runMuHad     = 1;

  //----------------------------------------
  // add samples to TChains
  //----------------------------------------

  TChain* chElHad     = new  TChain("Events");

  if(runElHad){
    cout << "adding ElectronHad data" << endl;

    //pickSkimIfExists(chElHad,"/hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-34/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged_ntuple_180252_0.root");

    //pickSkimIfExists(chElHad,"/hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-34/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged_ntuple_179959_0.root");


    pickSkimIfExists(chElHad,"/hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged_ntuple_178*root");
    pickSkimIfExists(chElHad,"/hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-34/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged_ntuple*root");

  }

  TChain* chMuHad     = new  TChain("Events");

  if(runMuHad){
    cout << "adding MuHad data" << endl;

    //pickSkimIfExists(chMuHad,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-35/MuHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-35_merged/V04-02-35/merged_ntuple_180252_0.root");

    pickSkimIfExists(chMuHad,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged_ntuple_178*root");
    pickSkimIfExists(chMuHad,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-35/MuHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-35_merged/V04-02-35/merged*root");
  }

  //----------------------------------------
  // run on samples
  //----------------------------------------

  if (runElHad) {
    cout << "Processing ElHad data" << endl;
    looper->ScanChain(chElHad,"elhad");
    cout << "Done processing ElHad data" << endl;
  }

  if (runMuHad) {
    cout << "Processing MuHad data" << endl;
    looper->ScanChain(chMuHad,"muhad");
    cout << "Done processing MuHadHad data" << endl;
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
