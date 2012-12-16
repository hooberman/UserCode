
#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "singleLeptonLooper.h"

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
  
  const char* version    = "V00-01-06";
  const char* jsonfile   = "jsons/Cert_EPSFINAL_May10ReReco_v2_PromptReco_160404_167913_JSON_goodruns.txt";
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
  //gROOT->ProcessLine(".L ../CORE/topmass/ttdilepsolve.cpp+"); 

  // Load FWLite
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  // Load and compile the looping code
  gSystem->CompileMacro("singleLeptonLooper.C","++k", "libsingleLeptonLooper");
  
  singleLeptonLooper* looper = new singleLeptonLooper();

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

  // k-factors
  float kttall    = 1.;
  float kqcd      = 1.;  
  float kWjets    = 1.;  
  float kDYtot    = 1.;  
  float ktW       = 1.;

  // prescales
  int preqcd      = 1;
  int prettall    = 1;
  int preWjets    = 1;
  int preDYtot    = 1;
  int pretW       = 1;
 
  // flags for files to run over
  bool rundata     = 0;
  bool runttall    = 1;
  bool runWjets    = 0;
  bool runQCD      = 0;
  bool runtW       = 0;
  bool runDYtot    = 0;
  bool runT2tt     = 0; 
  bool runT2tt_few = 0;
  bool runT2bw     = 0;

  if( useMCSkims )  cout << "Using MC skims" << endl;
  else              cout << "Using full MC samples" << endl;

  //----------------------------------------
  // QCD
  //----------------------------------------

  TChain* chQCD = new  TChain("Events");

  if(runQCD){

    string skimdir = "/home/users/benhoob/filters/output/";

    //SingleLeptonSkim ntuples
    if( useMCSkims ){
      pickSkimIfExists(chQCD, skimdir+"QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
      pickSkimIfExists(chQCD, skimdir+"QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
      pickSkimIfExists(chQCD, skimdir+"QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
      pickSkimIfExists(chQCD, skimdir+"QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
      pickSkimIfExists(chQCD, skimdir+"QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
      pickSkimIfExists(chQCD, skimdir+"QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    }

    //full cms2 ntuples
    else{
      pickSkimIfExists(chQCD, "/nfs-7/userdata/cms2/QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/merged*root");
      pickSkimIfExists(chQCD, "/nfs-7/userdata/cms2/QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/merged*root");
      pickSkimIfExists(chQCD, "/nfs-7/userdata/cms2/QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/merged*root");
      pickSkimIfExists(chQCD, "/nfs-7/userdata/cms2/QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/merged*root");
      pickSkimIfExists(chQCD, "/nfs-7/userdata/cms2/QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/merged*root");
      pickSkimIfExists(chQCD, "/nfs-7/userdata/cms2/QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/merged*root");
    }
  }

  //----------------------------------------
  // ttbar
  //----------------------------------------

  TChain* chtopall = new TChain("Events");
  if (runttall) {
    pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged*root");
    //pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged_ntuple.root");
  }

  //----------------------------------------
  // W+jets
  //----------------------------------------

  TChain* chWjets = new  TChain("Events");
  if(runWjets){

    //SingleLeptonSkim ntuples
    if( useMCSkims ){
      pickSkimIfExists(chWjets,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged_ntuple_*.root");
      //pickSkimIfExists(chWjets,"/nfs-7/userdata/cms2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndJets/merged*root");//l+3j filtered
      //pickSkimIfExists(chWjets,"/hadoop/cms/store/user/imacneill/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndJets/merged*root");
    }

    //full cms2 ntuples
    else{
      pickSkimIfExists(chWjets,"/nfs-7/userdata/cms2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged*root");
    }

  }

  //----------------------------------------
  // DY
  //----------------------------------------

  TChain* chDYtot = new  TChain("Events");
  if(runDYtot){
    cout << "ERROR! NEED TO ADD SINGLE LEPTON SKIM OF DY SAMPLES!!!" << endl;
    exit(0);
  }

  //----------------------------------------
  // single top
  //----------------------------------------
  
  TChain* chtW = new  TChain("Events");
  if (runtW) {

    pickSkimIfExists(chtW,"/nfs-7/userdata/benhoob/cms2/T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/benhoob/cms2/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    
    pickSkimIfExists(chtW,"/nfs-7/userdata/benhoob/cms2/T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/benhoob/cms2/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

    pickSkimIfExists(chtW,"/nfs-7/userdata/benhoob/cms2/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/benhoob/cms2/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

  }

  //----------------------------------------
  // T2tt (a few sample points)
  //----------------------------------------

  TChain *chT2tt_few = new TChain("Events");
  if (runT2tt_few) {
    
    // these are the T2tt files with m(stop),m(LSP) = 400,100 OR 500,100 OR 500,300
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_37.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_46.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_66.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_71.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_84.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_85.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_88.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_147.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_149.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_151.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_191.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_193.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_212.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_213.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_222.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_225.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_232.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_268.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_328.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_390.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_395.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_396.root");
    
    // these are the T2tt files with m(stop),m(LSP) = 250,50 OR 300,50 OR 400,50 OR 400,200
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_1.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_5.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_6.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_16.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_62.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_63.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_70.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_91.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_105.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_149.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_150.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_160.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_171.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_243.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_246.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_264.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_274.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_312.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_315.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_325.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_326.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_333.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_335.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_349.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_351.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_369.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_388.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_394.root");

    //350/100 (from ATLAS paper)
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_122.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_274.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_286.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_326.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_347.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_364.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_391.root");
    pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_392.root");
  }

  //----------------------------------------
  // T2tt (a few sample points)
  //----------------------------------------

  TChain *chT2tt = new TChain("Events");
  if (runT2tt) {
    pickSkimIfExists(chT2tt,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged*root");
  }

  //----------------------------------------
  // T2bw
  //----------------------------------------

  TChain *chT2bw = new TChain("Events");
  if (runT2bw) {
    // pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged_ntuple_98.root");
    //    pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged_ntuple_99.root");
    //    pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged_ntuple_100.root");
    //    pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged_ntuple_101.root");
    pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root");
  }

  //----------------------------------------
  // DATA
  //----------------------------------------

  TChain* chdata     = new  TChain("Events");

  if(rundata){
    
    cout << "adding ElectronHad and MuHad or SingleMu data" << endl;
    
    //---------------------------
    // May10 rereco
    //---------------------------
   
    //pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root"); 
    //pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");    

    // pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");//l+3j filtered
    // pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");//l+3j filtered

    pickSkimIfExists(chdata,"/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");//l+2j filtered
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/ElectronHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");//l+2j filtered

    //---------------------------
    // prompt reco v4
    //---------------------------
    
    //pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
    //pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");

    // pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");//l+3j filtered
    // pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");//l+3j filtered

    pickSkimIfExists(chdata,"/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");//l+2j filtered
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/ElectronHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");//l+2j filtered
    
  }



  
  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  float lumi              = 1.0; 
  bool  calculateTCMET    = false; //redo tcmet calculation on the fly

  char* jetTypeStrings[3] = {"JPT", "calo","pfjet"};
  char* metTypeStrings[4] = {"tcmet", "muon", "muonjes","pfmet"};
  char* zvetoStrings[4]   = {"", "_allzveto", "_nozveto","_selectz"};
  char* frmodeStrings[2]  = {"QCDType","WjetsType"}; //e_qcd = 0, e_wjets
  bool doFakeApp          = false;

  singleLeptonLooper::TrigEnum trig;

  for (int jetTypeIdx = 2; jetTypeIdx < 3; ++jetTypeIdx)
    {
      for (int metTypeIdx = 3; metTypeIdx < 4; ++metTypeIdx)
	{
	  for (int zvetoIdx = 0; zvetoIdx < 1; ++zvetoIdx)
	    {
	      for (int frmodeIdx = 0; frmodeIdx < (2-(1*!doFakeApp)); ++frmodeIdx)
		//for (int frmodeIdx = 0; frmodeIdx < 1; ++frmodeIdx)
		{
                  
		  singleLeptonLooper::JetTypeEnum  jetType(jetTypeIdx);
		  singleLeptonLooper::MetTypeEnum  metType(metTypeIdx);
		  singleLeptonLooper::ZVetoEnum    zveto(zvetoIdx);
		  singleLeptonLooper::FREnum       frmode(frmodeIdx);

		  if( doFakeApp ){
		    if( frmodeIdx == 0 ) cout << "Doing double fake estimate" << endl;
		    if( frmodeIdx == 1 ) cout << "Doing single fake estimate" << endl;
		  }

		  //--------------------------------------------------------------------
		  if (runttall) {
		    cout << "Processing ttbar all.. " << endl;
		    looper->ScanChain(chtopall,"ttall", kttall, prettall, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing ttbar all.. " << endl;
		  }
		  //--------------------------------------------------------------------
		  if (rundata) {
		    cout << "Processing data" << endl;
		    looper->ScanChain(chdata,"data", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing data" << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runDYtot) {
		    cout << "Processing DY->all" << endl;
		    looper->ScanChain(chDYtot,"DYtot", kDYtot, preDYtot, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done rocessing DY->ee" << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runQCD) {
		    cout << "Processing QCD.. " << endl;
		    looper->ScanChain(chQCD,"qcd", kqcd, preqcd, lumi, jetType, metType, zveto,frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing  QCD.. " << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runWjets) {
		    cout << "Processing Wjets.." << endl;
		    looper->ScanChain(chWjets,"wjets", kWjets, preWjets, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing Wjets.." << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runtW) {
		    cout << "Processing tW" << endl;
		    looper->ScanChain(chtW,"tW", ktW, pretW, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing tW" << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runT2tt) {
		    cout << "Processing T2tt" << endl;
		    looper->ScanChain(chT2tt, "T2tt", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing T2tt" << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runT2tt_few) {
		    cout << "Processing T2tt_few" << endl;
		    looper->ScanChain(chT2tt_few, "T2tt_few", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing T2tt_few" << endl;
		  }
		  //--------------------------------------------------------------------
		  if (runT2bw) {
		    cout << "Processing T2bw all.. " << endl;
		    looper->ScanChain(chT2bw,"T2bw", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		    cout << "Done processing T2bw all.. " << endl;
		  }

		  // save all the histograms
		  const char* outFile;
		  if(doFakeApp) {
		    outFile = Form("../output/%s/singleLepton_%s_%s%s_%s_FakeApp.root", version,
				   jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx],zvetoStrings[zvetoIdx],frmodeStrings[frmode]);
		  }
		  else {
		    outFile = Form("../output/%s/singleLepton_%s_%s%s.root", version,
				   jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx],zvetoStrings[zvetoIdx]);
		  }
                  
		  //const char* outFile = Form("victory_baseline_genmetgt50_nosumjetptcut_%s_%s_pleasework_varbins.root", 
		  //jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx]);
		  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
		  rootdir->cd();
		  saveHist(outFile);
		  deleteHistos();
                  
		} // frmodeIdx
	    }//zvetoIdx
	} // metTypeIdx
    } // jetTypeIdx
  //}

  gSystem->Exit(0);
  
}
