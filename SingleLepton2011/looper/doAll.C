
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
  
  const char* version    = "V00-04-00";
  const char* jsonfile   = "jsons/Cert_160404-180252_7TeV_mergePromptMay10Aug5_JSON_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  // Load Everything                                                                                                                                                                                                          
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gSystem->Load("libsingleLeptonCORE.so");
  gSystem->Load("libsingleLeptonLooper.so");

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
  float kVV       = 1.;
  float kDYtot    = 1.;  
  float ktW       = 1.;

  // prescales
  int preqcd      = 1;
  int prettall    = 1;
  int preWjets    = 1;
  int preVV       = 1;
  int preDYtot    = 1;
  int pretW       = 1;
 
  // flags for files to run over
  bool rundata     = 0;
  bool runttall    = 0;
  bool runWjets    = 0;
  bool runVV       = 1;
  bool runQCD      = 0;
  bool runMuQCD    = 0;
  bool runtW       = 0;
  bool runDYtot    = 0;
  bool runT2tt     = 0; 
  bool runT2tt_few = 0;
  bool runT2bw     = 0;

  bool rundatamay10   = 0;
  bool rundataprv4    = 0;
  bool rundataaug05   = 0;
  bool rundataprv6    = 0;
  bool rundata2011b33 = 0;
  bool rundata2011b34 = 0;

  if( useMCSkims )  cout << "Using MC skims" << endl;
  else              cout << "Using full MC samples" << endl;

  //----------------------------------------
  // QCD
  //----------------------------------------

  TChain* chQCD = new  TChain("Events");

  if(runQCD){
    string skimdir = "/home/users/benhoob/filters/output/";
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
  }

  //----------------------------------------
  // Muon QCD
  //----------------------------------------

  TChain* chMuQCD = new  TChain("Events");

  if(runMuQCD){

    string skimdir = "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/";
    string cms2dir = "/nfs-7/userdata/cms2/";

    //SingleLeptonSkim ntuples
    if( useMCSkims ){
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    }
    //full cms2 ntuples  
    else{
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    }
  }

  //----------------------------------------
  // ttbar
  //----------------------------------------

  TChain* chtopall = new TChain("Events");
  if (runttall) {
    //    pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged_ntuple_35.root");
    pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged*root");
  }

  //----------------------------------------
  // W+jets
  //----------------------------------------

  TChain* chWjets = new  TChain("Events");
  if(runWjets){
    pickSkimIfExists(chWjets,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
  }

  //----------------------------------------                                                                                                                                                                
  // Diboson VV                                                                                                                                                                                             
  //----------------------------------------                                                                                                                                                                

  TChain* chVV = new  TChain("Events");
  if(runVV){
    pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    // samples in multi-lepton decay modes
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v2/V04-02-29/SingleLepton/merged*root");
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
    // pickSkimIfExists(chVV, "/nfs-6/userdata/cms2/ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/merged*root");
  }

  //----------------------------------------
  // DY
  //----------------------------------------

  TChain* chDYtot = new  TChain("Events");
  if(runDYtot){
    string dypath = "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/";
    pickSkimIfExists(chDYtot,dypath+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    // samples in specific decay modes
    // pickSkimIfExists(chDYtot,dypath+"DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    // pickSkimIfExists(chDYtot,dypath+"DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    // pickSkimIfExists(chDYtot,dypath+"DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    // pickSkimIfExists(chDYtot,dypath+"DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    // pickSkimIfExists(chDYtot,dypath+"DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    // pickSkimIfExists(chDYtot,dypath+"DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
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
    pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root");
  }

  //----------------------------------------
  // DATA: May10
  //----------------------------------------

  TChain* chdatamay10 = new  TChain("Events");

  if(rundatamay10){

    cout << "adding ElectronHad, MuHad, SingleMu May10 data" << endl;
    
    // May10
    pickSkimIfExists(chdatamay10,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdatamay10,"/nfs-7/userdata/cms2/MuHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdatamay10,"/nfs-3/userdata/cms2/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

  }

  //----------------------------------------
  // DATA: PRv4
  //----------------------------------------

  TChain* chdataprv4 = new  TChain("Events");

  if(rundataprv4){
    
    cout << "adding ElectronHad, MuHad, SingleMu PRv4 data" << endl;

    // PRv4
    pickSkimIfExists(chdataprv4,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdataprv4,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdataprv4,"/nfs-3/userdata/cms2/SingleMu_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");
  }

  //----------------------------------------
  // DATA: Aug05
  //----------------------------------------

  TChain* chdataaug05 = new  TChain("Events");

  if(rundataaug05){
    
    cout << "adding ElectronHad, MuHad, SingleMu Aug05 data" << endl;
    
    // Aug05
    pickSkimIfExists(chdataaug05,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdataaug05,"/nfs-7/userdata/cms2/MuHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdataaug05,"/nfs-3/userdata/cms2/SingleMu_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

  }

  //----------------------------------------
  // DATA: PRv6
  //----------------------------------------

  TChain* chdataprv6     = new  TChain("Events");

  if(rundataprv6){
    
    cout << "adding ElectronHad, MuHad, SingleMu PRv6 data" << endl;

    // PRv6
    pickSkimIfExists(chdataprv6,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdataprv6,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdataprv6,"/nfs-3/userdata/cms2/SingleMu_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");
  }

  //----------------------------------------
  // DATA: 2011B V33
  //----------------------------------------

  TChain* chdata2011b33 = new  TChain("Events");

  if(rundata2011b33){
    
    cout << "adding ElectronHad, MuHad, SingleMu 2011B V33 data" << endl;

    // 2011B
    pickSkimIfExists(chdata2011b33,"/nfs-7/userdata/cms2/ElectronHad_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata2011b33,"/nfs-7/userdata/cms2/MuHad_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata2011b33,"/nfs-3/userdata/cms2/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

  }

  //----------------------------------------
  // DATA: 2011B V34
  //----------------------------------------

  TChain* chdata2011b34 = new  TChain("Events");

  if(rundata2011b34){
    
    cout << "adding ElectronHad, MuHad, SingleMu 2011B V34 data" << endl;
    // 2011B
    pickSkimIfExists(chdata2011b34,"/nfs-7/userdata/cms2/ElectronHad_Run2011B-PromptReco-v1_AOD/V04-02-34/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata2011b34,"/nfs-7/userdata/cms2/MuHad_Run2011B-PromptReco-v1_AOD/V04-02-35/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata2011b34,"/nfs-3/userdata/cms2/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/SingleLeptonAndTwoJets/merged*root");
    
  }

  //----------------------------------------
  // DATA
  //----------------------------------------

  TChain* chdata     = new  TChain("Events");

  if(rundata){
    
    cout << "adding ElectronHad and MuHad or SingleMu data" << endl;

    //pickSkimIfExists(chdata,"/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_999999_0_skim.root");//l+2j filtered
    
    // May10
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-3/userdata/cms2/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");   
    // PRv4
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-3/userdata/cms2/SingleMu_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

    // Aug05
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-3/userdata/cms2/SingleMu_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

    // PRv6
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-3/userdata/cms2/SingleMu_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

    // 2011B
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-3/userdata/cms2/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged*root");

    // 2011B
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011B-PromptReco-v1_AOD/V04-02-34/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011B-PromptReco-v1_AOD/V04-02-35/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-3/userdata/cms2/SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/SingleLeptonAndTwoJets/merged*root");
  }

  
  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  float lumi              = 1.0; 
  
  //--------------------------------------------------------------------
  if (rundatamay10) {
    cout << "Processing datamay10" << endl;
    looper->ScanChain(chdatamay10,"datamay10", 1, 1, lumi);
    cout << "Done processing datamay10" << endl;
  }
  //--------------------------------------------------------------------
  if (rundataprv4) {
    cout << "Processing dataprv4" << endl;
    looper->ScanChain(chdataprv4,"dataprv4", 1, 1, lumi);
    cout << "Done processing dataprv4" << endl;
  }
  //--------------------------------------------------------------------
  if (rundataaug05) {
    cout << "Processing dataaug05" << endl;
    looper->ScanChain(chdataaug05,"dataaug05", 1, 1, lumi);
    cout << "Done processing dataaug05" << endl;
  }
		  //--------------------------------------------------------------------
  if (rundataprv6) {
    cout << "Processing dataprv6" << endl;
    looper->ScanChain(chdataprv6,"dataprv6", 1, 1, lumi);
    cout << "Done processing dataprv6" << endl;
  }
  //--------------------------------------------------------------------
  if (rundata2011b33) {
    cout << "Processing data2011b33" << endl;
    looper->ScanChain(chdata2011b33,"data2011b33", 1, 1, lumi);
    cout << "Done processing data2011b33" << endl;
  }
  //--------------------------------------------------------------------
  if (rundata2011b34) {
    cout << "Processing data2011b34" << endl;
    looper->ScanChain(chdata2011b34,"data2011b34", 1, 1, lumi);
    cout << "Done processing data2011b34" << endl;
  }
  //--------------------------------------------------------------------
  if (runttall) {
    cout << "Processing ttbar all.. " << endl;
    looper->ScanChain(chtopall,"ttall", kttall, prettall, lumi);
    cout << "Done processing ttbar all.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runDYtot) {
    cout << "Processing DY->all" << endl;
    looper->ScanChain(chDYtot,"DYtot", kDYtot, preDYtot, lumi);
    cout << "Done rocessing DY->ee" << endl;
  }
  //--------------------------------------------------------------------
  if (runQCD) {
    cout << "Processing QCD.. " << endl;
    looper->ScanChain(chQCD,"qcd", kqcd, preqcd, lumi);
    cout << "Done processing  QCD.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runMuQCD) {
    cout << "Processing Mu QCD.. " << endl;
    looper->ScanChain(chMuQCD,"muqcd", kqcd, preqcd, lumi);
    cout << "Done processing  QCD.. " << endl;
                  }
  //--------------------------------------------------------------------
  if (runWjets) {
    cout << "Processing Wjets.." << endl;
    looper->ScanChain(chWjets,"wjets", kWjets, preWjets, lumi);
    cout << "Done processing Wjets.." << endl;
  }
  //-------------------------------------------------------------------
  if (runVV) {
    cout << "Processing Diboson.." << endl;
    looper->ScanChain(chVV,"diboson", kVV, preVV, lumi);
    cout << "Done processing Diboson.." << endl;
  }
  //--------------------------------------------------------------------
  if (runtW) {
    cout << "Processing tW" << endl;
    looper->ScanChain(chtW,"tW", ktW, pretW, lumi);
    cout << "Done processing tW" << endl;
  }
  //--------------------------------------------------------------------
  if (runT2tt) {
    cout << "Processing T2tt" << endl;
    looper->ScanChain(chT2tt, "T2tt", 1, 1, lumi);
    cout << "Done processing T2tt" << endl;
  }
  //--------------------------------------------------------------------
  if (runT2tt_few) {
    cout << "Processing T2tt_few" << endl;
    looper->ScanChain(chT2tt_few, "T2tt_few", 1, 1, lumi);
    cout << "Done processing T2tt_few" << endl;
  }
  //--------------------------------------------------------------------
  if (runT2bw) {
    cout << "Processing T2bw all.. " << endl;
    looper->ScanChain(chT2bw,"T2bw", 1, 1, lumi);
    cout << "Done processing T2bw all.. " << endl;
  }
  //--------------------------------------------------------------------
  if (rundata) {
    cout << "Processing data" << endl;
    looper->ScanChain(chdata,"data", 1, 1, lumi);
    cout << "Done processing data" << endl;
  }
  //--------------------------------------------------------------------
  
  gSystem->Exit(0);
  
}
