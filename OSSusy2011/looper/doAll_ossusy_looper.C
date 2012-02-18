
#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "ossusy_looper.h"

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

void doAll_ossusy_looper(bool skipFWLite = true)
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version   = "V00-02-20";
  const char* jsonfile  = "jsons/Cert_160404-180252_7TeV_mergePromptMay10Aug5_JSON_goodruns.txt";

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
  gSystem->CompileMacro("ossusy_looper.C","++k", "libossusy_looper");
  
  ossusy_looper* looper = new ossusy_looper();

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

  //qcdpt15 correction
  // 	    * (8.158/8.762) // approximate ratio of cross section for pt-hat > 15 to 15 < pt-hat < 30
  // 	    * (490779./381623.); // ratio of number events with pt-hat < 30 to total with pt-hat > 15 ;; for the early data sample qcdPt15 and qcdPt30 MC
  float kqcdpt15  = (8.158/8.762)*(490779./381623.);  
  float kqcdpt30  = 1.;  
  float kqcd      = 1.;  
  float kttall    = 1.;  //157.5/165.0;  
  float ktt42     = 1.;  //157.5/165.0;  
  float kttpowheg = 1.;  //157.5/165.0;  
  float kttdil    = 1.;  //157.5/165.0;  
  float kttrelval = 1;  
  float kttem     = 1.;  //157.5/165.0;  
  float kttotr    = 1.;  //157.5/165.0;  
  float kVV       = 1.3;
  float kWW       = 1.;
  float kWZ       = 1.;
  float kZZ       = 1.;
  float kWjetsMG  = 1.;  // 11850 pb, 980000 events processed
  float kWjets    = 1.;  // 11850 pb, 980000 events processed
  float kWcharm   = 1.1;
  float kZjets    = 1.;  
  float kDYtot    = 1.;  
  float kDYee     = 1.;  // 1230 pb,  970360 events processed
  float kDYmm     = 1.;  // 1230 pb,  970360 events processed
  float kDYtautau = 1.;  // 1230 pb,  970360 events processed
  float kppMuX    = 1.;  // xsec/nevents
  float kEM       = 1.;
  float ktW       = 1.;  // the evtScale is all negative for some reason
  float kVQQ      = 1.;
  float kLM0      = 1.;
  float kLM1      = 1.;
  float kLM2      = 1.;
  float kLM3      = 1.;
  float kLM4      = 1.;
  float kLM5      = 1.;
  float kLM6      = 1.;
  float kLM7      = 1.;
  float kLM8      = 1.;
  float kLM9      = 1.;
  float kLM10     = 1.;
  float kLM11     = 1.;
  float kLM12     = 1.;
  float kLM13     = 1.;
  float kLMscan   = 1.;
  float kML1      = 1.;
  float kML2      = 1.;
  float kML3      = 1.;
  float kML4      = 1.;
  float kML5      = 1.;
  float kML6      = 1.;
  float kML7      = 1.;
  float kML8      = 1.;

  // Prescales
  int preqcdpt15  = 1;
  int preqcdpt30  = 1;
  int preqcd      = 1;
  int prettall    = 1;
  int prett42     = 1;
  int prettpowheg = 1;
  int prettdil    = 1;
  int prettem     = 1;
  int prettrelval = 1;
  int prettotr    = 1;
  int preVV       = 1;
  int preWW       = 1;
  int preWZ       = 1;
  int preZZ       = 1;
  int preWjets    = 1;
  int preWjetsMG  = 1;
  int preWcharm   = 1;
  int preZjets    = 1;
  int preDYee     = 1;
  int preDYtot    = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;
  int preppMuX    = 1;
  int preEM       = 1;
  int pretW       = 1;
  int preVQQ      = 1;
  int preLM0      = 1;
  int preLM1      = 1;
  int preLM2      = 1;
  int preLM3      = 1;
  int preLM4      = 1;
  int preLM5      = 1;
  int preLM6      = 1;
  int preLM7      = 1;
  int preLM8      = 1;
  int preLM9      = 1;
  int preLM10     = 1;
  int preLM11     = 1;
  int preLM12     = 1;
  int preLM13     = 1;
  int preML1      = 1;
  int preML2      = 1;
  int preML3      = 1;
  int preML4      = 1;
  int preML5      = 1;
  int preML6      = 1;
  int preML7      = 1;
  int preML8      = 1;
  int preLMscan   = 1;

  /*
  //Flags for files to run over
  bool rundata     = 0;
  bool rundata41   = 0;
  bool rundataskim = 0;
  bool runQCDpt15  = 0;
  bool runQCDpt30  = 0;
  bool runQCD      = 0;
  bool runttall    = 0;
  bool runtt42     = 0;
  bool runttpowheg = 1;
  bool runttdil    = 0;
  bool runttem     = 0;
  bool runttrelval = 0;
  bool runttotr    = 0;
  bool runVV       = 0;
  bool runWW       = 1;
  bool runWZ       = 1;
  bool runZZ       = 1;
  bool runWjets    = 1;
  bool runWjetsMG  = 0;
  bool runWcharm   = 0;
  bool runZjets    = 0;
  bool runDYtot    = 1;
  bool runDYee     = 0;
  bool runDYmm     = 0;
  bool runDYtautau = 0;
  bool runppMuX    = 0;
  bool runEM       = 0;
  bool runtW       = 1;
  bool runVQQ      = 0;
  bool runLM0      = 0;
  bool runLM1      = 0;
  bool runLM2      = 0;
  bool runLM3      = 0;
  bool runLM4      = 0;
  bool runLM5      = 0;
  bool runLM6      = 0;
  bool runLM7      = 0;
  bool runLM8      = 0;
  bool runLM9      = 0;
  bool runLM10     = 0;
  bool runLM11     = 0;
  bool runLM12     = 0;
  bool runLM13     = 0;
  bool runML1      = 0;
  bool runML2      = 0;
  bool runML3      = 0;
  bool runML4      = 0;
  bool runML5      = 0;
  bool runML6      = 0;
  bool runML7      = 0;
  bool runML8      = 0;
  bool runLMscan   = 0; 
  bool runT1lh     = 0;
  bool runT2tt     = 0;
  */
    
  //Flags for files to run over
  bool rundata     = 0;
  bool rundata41   = 0;
  bool rundataskim = 0;
  bool runQCDpt15  = 0;
  bool runQCDpt30  = 0;
  bool runQCD      = 0;
  bool runphotons  = 0;
  bool runttall    = 0;
  bool runttallPUS6= 0;
  bool runttpowheg = 0;
  bool runtt42     = 0;
  bool runttdil    = 0;
  bool runttrelval = 0;
  bool runttem     = 0;
  bool runttotr    = 0;
  bool runVV       = 0;
  bool runWW       = 0;
  bool runWZ       = 0;
  bool runZZ       = 0;
  bool runWjets    = 0;
  bool runWjetsMG  = 0;
  bool runWcharm   = 0;
  bool runZjets    = 0;
  bool runDYtot    = 0;
  bool runDYee     = 0;
  bool runDYmm     = 0;
  bool runDYtautau = 0;
  bool runppMuX    = 0;
  bool runEM       = 0;
  bool runtW       = 0;
  bool runVQQ      = 0;
  bool runLM0      = 0;
  bool runLM1      = 0;
  bool runLM1v2    = 0;
  bool runLM2      = 0;
  bool runLM3      = 0;
  bool runLM3v2    = 0;
  bool runLM4      = 0;
  bool runLM5      = 0;
  bool runLM6      = 0;
  bool runLM6v2    = 0;
  bool runLM7      = 0;
  bool runLM8      = 0;
  bool runLM9      = 0;
  bool runLM10     = 0;
  bool runLM11     = 0;
  bool runLM12     = 0;
  bool runLM13     = 0;
  bool runML1      = 0;
  bool runML2      = 0;
  bool runML3      = 0;
  bool runML4      = 0;
  bool runML5      = 0;
  bool runML6      = 0;
  bool runML7      = 0;
  bool runML8      = 0;
  bool runLMscan   = 0; 
  bool runLMscanFall11      = 0; 
  bool runLMscanFall11dil   = 0; 
  bool runLMscanFall11dil1  = 1; 
  bool runLMscanFall11dil2  = 0; 
  bool runLMscanFall11dil3  = 0; 
  bool runLMscanFall11dil4  = 0; 
  bool runLMscanFall11dil5  = 0; 
  bool runLMscanFall11dil6  = 0; 
  bool runLMscanFall11dil7  = 0; 
  bool runLMscanFall11dil8  = 0; 
  bool runLMscanFall11dil9  = 0; 
  bool runLMscanFall11dil10 = 0; 
  bool runT2tt     = 0;
  bool runT1lh     = 0;
  bool runZZZ      = 0;
  
  char* dir = "";

  bool useMCSkims = true;
  if( useMCSkims ){
    cout << "Using MC skims" << endl;
    dir = "met50/";
  }
  else{
    cout << "Using full MC samples" << endl;
  }

  TChain* chdataskim = new  TChain("Events");
  if(rundataskim){
    
    pickSkimIfExists(chdataskim,
                     "/tas/cms2/dilepSkim35pb/els.root",
                     "dataskim");

    pickSkimIfExists(chdataskim,
                     "/tas/cms2/dilepSkim35pb/mus.root",
                     "dataskim");
  }
  
  // TChain* chQCDpt15 = new  TChain("Events");
  // if(runQCDpt15){
  //   pickSkimIfExists(chQCDpt15, 
  //                    "/tas/cms2/QCD_Pt15_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/*root",
  //                    "QCDpt15");

  // }
  
  // TChain* chQCDpt30 = new  TChain("Events");
  // if(runQCDpt30){
  //   pickSkimIfExists(chQCDpt30, 
  //                    "/tas/cms2/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/*root",
  //                    "QCDpt30");
  // }

  TChain* chQCD = new  TChain("Events");
  if(runQCD){
    pickSkimIfExists(chQCD, 
		     "cms2/QCD_Pt_15to30_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-03/merged*root",
		     "QCD");

    pickSkimIfExists(chQCD, 
		     "cms2/QCD_Pt_30to50_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-03/merged*root",
		     "QCD");

    pickSkimIfExists(chQCD, 
		     "cms2/QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-03/merged*root",
		     "QCD");


  }

  TChain* chphotons = new  TChain("Events");
  if(runphotons){
    pickSkimIfExists(chphotons,"/nfs-7/userdata/cms2/G_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/merged*root"); 
  }

  TChain* chZZZ = new  TChain("Events");
  if(runZZZ){
    pickSkimIfExists(chZZZ,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/ZZZNoGstar_TuneZ2_7TeV-madgraphCMSSW42xPUv1_spadhi-ZZZNoGstar_TuneZ2_7TeV-madgraphCMSSW42xPUv1-9ab11d163a88ab8f3641ab081403ebc5/VB04-02-29_FastSim/merged*root");
  }

  TChain* chZjets = new  TChain("Events");
  if(runZjets){

    //pickSkimIfExists(chDYtot,"/nfs-7/userdata/cms2/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chZjets,"/nfs-7/userdata/cms2/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*root");
  }

  TChain* chtopall = new TChain("Events");
  if (runttall) {
    //pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root");
    pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  TChain* chtopallPUS6 = new TChain("Events");
  if (runttallPUS6) {
    //pickSkimIfExists(chtopallPUS6,"/nfs-6/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-29/merged_ntuple.*root");
    pickSkimIfExists(chtopallPUS6,"/nfs-6/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-29/merged*root");
  }

  TChain* chtop42 = new TChain("Events");
  if (runtt42) {
    pickSkimIfExists(chtop42,
		     "cms2/TT_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-12/merged*root",
		     "TTJets");
  }

  TChain* chtoppowheg = new TChain("Events");
  if (runttpowheg) {
    pickSkimIfExists(chtoppowheg, 
		     "/nfs-7/userdata/cms2/TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root",
		     "TTPowheg");
  }

  TChain* chtopdil = new TChain("Events");
  if (runttdil) {
    pickSkimIfExists(chtopdil,
		     "cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "TTJets");

  }

  TChain* chttrelval = new TChain("Events");
  if (runttrelval) {
    pickSkimIfExists(chttrelval,
		     "ttbar_relval_3_11.root",
                     "TTJets");
  }

  TChain* chtopem = new TChain("Events");
  if (runttem) {
    pickSkimIfExists(chtopem, 
		     "cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "TTJets");
  }

  TChain* chtopotr = new TChain("Events");
  if(runttotr){
    pickSkimIfExists(chtopotr, 
		     "cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "TTJets");

  }

  TChain* chvv = new TChain("Events");
  // if(runVV){
  //   pickSkimIfExists(chvv, 
  //                    Form("/tas/cms2/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Fall10-START38_V12-v1/V03-06-14/%smerged*root",dir),
  //                    "VV");
  // }

  TChain* chww = new TChain("Events");
  if(runWW){
    pickSkimIfExists(chww,"/nfs-7/userdata/cms2/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  TChain* chWZ = new TChain("Events");
  if(runWZ){
    pickSkimIfExists(chWZ,"/nfs-7/userdata/cms2/WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chWZ,"/nfs-7/userdata/cms2/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  TChain* chZZ = new TChain("Events");
  if(runZZ){
    pickSkimIfExists(chZZ,"/nfs-7/userdata/cms2/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root"); 
    pickSkimIfExists(chZZ,"/nfs-7/userdata/cms2/ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chZZ,"/nfs-7/userdata/cms2/ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  TChain* chWjets = new  TChain("Events");
  if(runWjets){
    pickSkimIfExists(chWjets,"/hadoop/cms/store/user/imacneill/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/DileptonHyp/merged*root");
  }

  TChain* chWjetsMG = new  TChain("Events");
  if(runWjetsMG){
    pickSkimIfExists(chWjetsMG, 
                     "cms2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
		     "WJets");
  }

  TChain* chWcharm = new TChain("Events");
  // if(runWcharm){
  //   pickSkimIfExists(chWcharm, 
  //                    //is this the correct sample????
  //  		     "/tas/cms2/WCJets_7TeV-madgraph_Spring10-START3X_V26-v1/V03-04-13-01/merged*root",
  //                    "Wc");
  // }


  TChain* chDYtot = new  TChain("Events");
  if(runDYtot){
    pickSkimIfExists(chDYtot,"/nfs-7/userdata/cms2/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/met50skim/merged*root");
    pickSkimIfExists(chDYtot,"/nfs-3/userdata/cms2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/met50skim/merged*root");
    pickSkimIfExists(chDYtot,"/nfs-7/userdata/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/met50skim/merged*root");
    pickSkimIfExists(chDYtot,"/nfs-3/userdata/cms2/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/met50skim/merged*root");
    pickSkimIfExists(chDYtot,"/nfs-3/userdata/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/met50skim/merged*root");
    pickSkimIfExists(chDYtot,"/nfs-3/userdata/cms2/DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/met50skim/merged*root");
    pickSkimIfExists(chDYtot,"/nfs-7/userdata/cms2/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/met50skim/merged*root");
  }
  
  TChain* chDYtautau = new  TChain("Events");
  
  if(runDYtautau){
  
    char* dypath = "/tas03/home/benhoob/OSSusy2011/filter/output";
    
    pickSkimIfExists(chDYtautau, 
		     Form("%s/DYToTauTau_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v2/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYtautau");
    
    pickSkimIfExists(chDYtautau, 
		     Form("%s/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYtautau");
    
    pickSkimIfExists(chDYtautau, 
		     Form("%s/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYtautau");
  }
  
  TChain* chDYee = new  TChain("Events");
 
  if(runDYee){
    
    char* dypath = "/tas03/home/benhoob/OSSusy2011/filter/output";
    
    pickSkimIfExists(chDYee, 
		     Form("%s/DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYee");

    pickSkimIfExists(chDYee, 
		     Form("%s/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYee");
    
    pickSkimIfExists(chDYee, 
		     Form("%s/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYee");
    
  }

  TChain* chDYmm = new  TChain("Events");

  if(runDYmm){
    
    char* dypath = "/tas03/home/benhoob/OSSusy2011/filter/output";
      
    pickSkimIfExists(chDYmm, 
		     Form("%s/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYmm");
      
    pickSkimIfExists(chDYmm, 
		     Form("%s/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYmm");
      
    pickSkimIfExists(chDYmm, 
		     Form("%s/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/met50skim/merged_ntuple.root",dypath),
		     "DYmm");
  }
  
  
  
  // ppMuX
  TChain* chppMuX = new  TChain("Events");
  if (runppMuX) {
    pickSkimIfExists(chppMuX, 
                     "/tas/cms2/InclusiveMu15_Spring10-START3X_V26_S09-v1/V03-04-13-07/",
                     "InclusiveMuPt15"); 
    // can try InclusiveMu5Pt50 .. figure out how to merge later
  }
  
  // ppEM
  //only pt80to170 sample is currently available!!!
  TChain* chEM =  new  TChain("Events");
  if (runEM) {
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
    pickSkimIfExists(chEM, 
                     "/tas/cms2/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                     "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
  }

  // tW
  TChain* chtW = new  TChain("Events");
  if (runtW) {
    pickSkimIfExists(chtW,"/nfs-7/userdata/cms2/T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/cms2/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/cms2/T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/cms2/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/cms2/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtW,"/nfs-7/userdata/cms2/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  // VQQ
  TChain* chVQQ = new TChain("Events");
  if (runVQQ) {
    pickSkimIfExists(chVQQ, 
                     //is this the correct sample???
                     "/tas/cms2/VqqJets-madgraph_Spring10-START3X_V26_S09-v1/merged*.root",
  		     //"data/VQQ-madgraph_Summer08_IDEAL_V11_redigi_v2/merged*.root", 
  		     "VQQ");
  }
  
  // LM0
  TChain *chLM0 = new TChain("Events");
  if (runLM0) {
    pickSkimIfExists(chLM0, 
		     "/nfs-7/userdata/cms2/LM0_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
                     "SUSY_LM0");

  }
  
  // LM1
  TChain *chLM1 = new TChain("Events");
  if (runLM1) {
    pickSkimIfExists(chLM1, 
		     "/nfs-7/userdata/cms2/LM1_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
		     "SUSY_LM1");
  }
  
  // LM1v2
  TChain *chLM1v2 = new TChain("Events");
  if (runLM1v2) {
    pickSkimIfExists(chLM1v2, 
		     "/nfs-7/userdata/cms2/LM1_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged*root",
		     "SUSY_LM1");
  }
  
  // LM2
  TChain *chLM2 = new TChain("Events");
  if (runLM2) {
    pickSkimIfExists(chLM2, 
		     "/nfs-7/userdata/cms2/LM2_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
                     "SUSY_LM2");
  }

  // LM3
  TChain *chLM3 = new TChain("Events");
  if (runLM3) {
    pickSkimIfExists(chLM3, 
		     "/nfs-7/userdata/cms2/LM3_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
                     "SUSY_LM3");
  }

  // LM3v2
  TChain *chLM3v2 = new TChain("Events");
  if (runLM3v2) {
    pickSkimIfExists(chLM3v2, 
		     "/nfs-7/userdata/cms2/LM3_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged*root",
                     "SUSY_LM3");
  }

  // LM4
  TChain *chLM4 = new TChain("Events");
  if (runLM4) {
    pickSkimIfExists(chLM4, 
		     "/nfs-7/userdata/cms2/LM4_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
                     "SUSY_LM4");
  }

  // LM5
  TChain *chLM5 = new TChain("Events");
  if (runLM5) {
    pickSkimIfExists(chLM5, 
		     "/nfs-7/userdata/cms2/LM5_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
                     "SUSY_LM5");
  }

  // LM6
  TChain *chLM6 = new TChain("Events");
  if (runLM6) {
    pickSkimIfExists(chLM6, 
		     "/nfs-7/userdata/cms2/LM6_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root",
                     "SUSY_LM6");
  }

  // LM6v2
  TChain *chLM6v2 = new TChain("Events");
  if (runLM6v2) {
    //pickSkimIfExists(chLM6v2,"/nfs-7/userdata/cms2/LM6_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged_ntuple.root");
    pickSkimIfExists(chLM6v2,"/nfs-7/userdata/cms2/LM6_SUSY_sftsht_7TeV-pythia6_Summer11-PU_S4_START42_V11-v2/V04-02-29/merged*root");
  }

  // LM7
  TChain *chLM7 = new TChain("Events");
  if (runLM7) {
    pickSkimIfExists(chLM7, 
                     "cms2/LM7_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM7");
  }

  // LM8
  TChain *chLM8 = new TChain("Events");
  if (runLM8) {
    pickSkimIfExists(chLM8, 
                     "cms2/LM8_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM8");
  }

  // LM9
  TChain *chLM9 = new TChain("Events");
  if (runLM9) {
    pickSkimIfExists(chLM9, 
                     "cms2/LM9_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM9");
  }

  // LM10
  TChain *chLM10 = new TChain("Events");
  if (runLM10) {
    pickSkimIfExists(chLM10, 
                     "cms2/LM10_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM10");
  }

  // LM11 
  TChain *chLM11 = new TChain("Events");
  if (runLM11) {
    pickSkimIfExists(chLM11, 
                     "cms2/LM11_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM11");
  }

  // LM12
  TChain *chLM12 = new TChain("Events");
  if (runLM12) {
    pickSkimIfExists(chLM12, 
                     "cms2/LM12_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM12");
  }

  // LM13
  TChain *chLM13 = new TChain("Events");
  if (runLM13) {
    pickSkimIfExists(chLM13, 
                     "cms2/LM13_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root",
                     "SUSY_LM13");
  }
  
  // ML1
  TChain *chML1 = new TChain("Events");
  if (runML1) {
    pickSkimIfExists(chML1, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML01_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML1");
  }

  // ML2
  TChain *chML2 = new TChain("Events");
  if (runML2) {
    pickSkimIfExists(chML2, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML02_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML2");
  }

  // ML3
  TChain *chML3 = new TChain("Events");
  if (runML3) {
    pickSkimIfExists(chML3, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML03_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML3");
  }

  // ML4
  TChain *chML4 = new TChain("Events");
  if (runML4) {
    pickSkimIfExists(chML4, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML04_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML4");
  }

  // ML5
  TChain *chML5 = new TChain("Events");
  if (runML5) {
    pickSkimIfExists(chML5, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML05_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML5");
  }

  // ML6
  TChain *chML6 = new TChain("Events");
  if (runML6) {
    pickSkimIfExists(chML6, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML06_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML6");
  }

  // ML7
  TChain *chML7 = new TChain("Events");
  if (runML7) {
    pickSkimIfExists(chML7, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML07_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML7");
  }
  
  // ML8
  TChain *chML8 = new TChain("Events");
  if (runML8) {
    pickSkimIfExists(chML8, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML08_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML8");
  }

  TChain *chLMscanFall11dil = new TChain("Events");
  if (runLMscanFall11dil) {

    //pickSkimIfExists(chLMscanFall11dil,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple*root");

    pickSkimIfExists(chLMscanFall11dil,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_147_1_sap.root");

  }

  TChain *chLMscanFall11dil1 = new TChain("Events");
  if (runLMscanFall11dil1) {
    pickSkimIfExists(chLMscanFall11dil1,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_11*root"); //111
    pickSkimIfExists(chLMscanFall11dil1,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_12*root"); //111
    pickSkimIfExists(chLMscanFall11dil1,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_13*root"); //111
  }

  TChain *chLMscanFall11dil2 = new TChain("Events");
  if (runLMscanFall11dil2) {
    pickSkimIfExists(chLMscanFall11dil2,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_14*root"); //111
    pickSkimIfExists(chLMscanFall11dil2,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_15*root"); //111
    pickSkimIfExists(chLMscanFall11dil2,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_16*root"); //111
  }

  TChain *chLMscanFall11dil3 = new TChain("Events");
  if (runLMscanFall11dil3) {
    pickSkimIfExists(chLMscanFall11dil3,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_17*root"); //111
    pickSkimIfExists(chLMscanFall11dil3,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_18*root"); //111
    pickSkimIfExists(chLMscanFall11dil3,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_19*root"); //111
  }

  TChain *chLMscanFall11dil4 = new TChain("Events");
  if (runLMscanFall11dil4) {
    pickSkimIfExists(chLMscanFall11dil4,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_20*root"); //111
    pickSkimIfExists(chLMscanFall11dil4,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_21*root"); //111
    pickSkimIfExists(chLMscanFall11dil4,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_22*root"); //111
  }

  TChain *chLMscanFall11dil5 = new TChain("Events");
  if (runLMscanFall11dil5) {
    pickSkimIfExists(chLMscanFall11dil5,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_23*root"); //111
    pickSkimIfExists(chLMscanFall11dil5,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_24*root"); //111
    pickSkimIfExists(chLMscanFall11dil5,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_25*root"); //111
  }

  TChain *chLMscanFall11dil6 = new TChain("Events");
  if (runLMscanFall11dil6) {
    pickSkimIfExists(chLMscanFall11dil6,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_26*root"); //111
    pickSkimIfExists(chLMscanFall11dil6,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_27*root"); //111
    pickSkimIfExists(chLMscanFall11dil6,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_28*root"); //111
  }

  TChain *chLMscanFall11dil7 = new TChain("Events");
  if (runLMscanFall11dil7) {
    pickSkimIfExists(chLMscanFall11dil7,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_3*root"); //295
  }

  TChain *chLMscanFall11dil8 = new TChain("Events");
  if (runLMscanFall11dil8) {
    pickSkimIfExists(chLMscanFall11dil8,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_29*root"); //111
    pickSkimIfExists(chLMscanFall11dil8,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_4*root"); //111
    pickSkimIfExists(chLMscanFall11dil8,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_5*root"); //111
  }

  TChain *chLMscanFall11dil9 = new TChain("Events");
  if (runLMscanFall11dil9) {
    pickSkimIfExists(chLMscanFall11dil9,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_6*root"); //111
    pickSkimIfExists(chLMscanFall11dil9,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_7*root"); //111
    pickSkimIfExists(chLMscanFall11dil9,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_8*root"); //111
  }

  TChain *chLMscanFall11dil10 = new TChain("Events");
  if (runLMscanFall11dil10) {
    pickSkimIfExists(chLMscanFall11dil10,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_9*root"); //111
    pickSkimIfExists(chLMscanFall11dil10,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_10*root"); //111
    pickSkimIfExists(chLMscanFall11dil10,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_1_*root"); //1
    pickSkimIfExists(chLMscanFall11dil10,"/nfs-7/userdata/cms2/mSUGRA_dilepton_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v6/VB04-02-29_Fastsim_mSUGRA_Dilep/preprocessing/ntuple_2_*root"); //1
  }

  TChain *chLMscanFall11 = new TChain("Events");
  if (runLMscanFall11) {

    pickSkimIfExists(chLMscanFall11,"/nfs-7/userdata/cms2/mSUGRA_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v5/VB04-02-29_Fastsim_mSUGRA/preprocessing/ntuple*root");


    // LM6
    //pickSkimIfExists(chLMscanFall11,"/nfs-7/userdata/cms2/mSUGRA_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v5/VB04-02-29_Fastsim_mSUGRA/preprocessing/ntuple_2596_*.root");
    //pickSkimIfExists(chLMscanFall11,"/nfs-7/userdata/cms2/mSUGRA_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v5/VB04-02-29_Fastsim_mSUGRA/preprocessing/ntuple_2597_*.root");

    //pickSkimIfExists(chLMscanFall11,"/nfs-7/userdata/cms2/mSUGRA_m0-220to3000_m12-100to1000_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v5/VB04-02-29_Fastsim_mSUGRA/preprocessing/ntuple_999_1_FO8.root");

  }

  // LMscan
  TChain *chLMscan = new TChain("Events");
  if (runLMscan) {

    pickSkimIfExists(chLMscan,
     		     "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/merged_ntuple*root",
		     "LMscan");

    /*   
    pickSkimIfExists(chLMscan,
		     "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/merged_ntuple_1_15.root",
                     "LMscan");
    
    pickSkimIfExists(chLMscan,
		     "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/merged_ntuple_1_19.root",
                     "LMscan");

    pickSkimIfExists(chLMscan,
		     "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/merged_ntuple_1_39.root",
                     "LMscan");

    pickSkimIfExists(chLMscan,
		     "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/merged_ntuple_1_42.root",
                     "LMscan");
    
    pickSkimIfExists(chLMscan,
		     "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/merged_ntuple_1_18.root",
                     "LMscan");

    */
  }

  TChain *chT2tt = new TChain("Events");
  if (runT2tt) {
    
    // pickSkimIfExists(chT2tt,
    // 		     "/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple.root",
    //                  "T2tt");  

    pickSkimIfExists(chT2tt,
    		     "/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged*root",
                     "T2tt");  
    
  }

  // LMscan
  TChain *chT1lh = new TChain("Events");
  if (runT1lh) {
    
    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged*root",
                     "T1lh");  
    /*
    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged_ntuple_3.root",
                     "T1lh");  

    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged_ntuple_4.root",
                     "T1lh");  

    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged_ntuple_5.root",
                     "T1lh");  

    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged_ntuple_7.root",
                     "T1lh");  

    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged_ntuple_9.root",
                     "T1lh");  

    pickSkimIfExists(chT1lh,
		     "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/merged_ntuple_10.root",
                     "T1lh");  
    */

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

  ossusy_looper::TrigEnum trig;

  TChain* chdata     = new  TChain("Events");
  TChain* chdata41   = new  TChain("Events");

  for( int pt = 0 ; pt < 1 ; ++pt ){

    //set trigger type
    if( pt == 0 ) trig = ossusy_looper::e_highpt;
    if( pt == 1 ) trig = ossusy_looper::e_lowpt;
  
    looper->set_trigger( trig );

    chdata->Reset();
    chdata41->Reset();
  
    if(rundata41){

      cout << "41X data obsolete! quitting" << endl;
      exit(0);

    }

    if(rundata){
      
      if( trig == ossusy_looper::e_highpt ){

	cout << "Doing high-pT dilepton trigger data" << endl;
	
	//---------------------------
	// May10 rereco
	//---------------------------

	pickSkimIfExists(chdata,"cms2_data/DoubleElectron_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skim*root");
	pickSkimIfExists(chdata,"cms2_data/DoubleMu_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skim*root");
	pickSkimIfExists(chdata,"cms2_data/MuEG_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skim*root");

	//---------------------------
	// prompt reco v4
	//---------------------------

	pickSkimIfExists(chdata,"cms2_data/DoubleElectron_Run2011A-PromptReco-v4_AOD/V04-02-20/DoubleElectronTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"cms2_data/DoubleMu_Run2011A-PromptReco-v4_AOD/V04-02-20/DoubleMuTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");

	//---------------------------
	// august rereco
	//---------------------------

	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleElectron_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleElectronTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleMu_Run2011A-05Aug2011-v1_AOD/V04-02-30/DoubleMuTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/MuEG_Run2011A-05Aug2011-v1_AOD/V04-02-30/SSignSkim/skim*root");

	//---------------------------
	// prompt reco v6
	//---------------------------
	
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleElectronTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleMu_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleMuTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/MuEG_Run2011A-PromptReco-v6_AOD/V04-02-30/SSignSkim/skim*root");
	
	//---------------------------
	// Run2011B prompt reco v1
	//---------------------------

	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleElectron_Run2011B-PromptReco-v1_AOD/V04-02-30/DoubleElectronTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleMu_Run2011B-PromptReco-v1_AOD/V04-02-30/DoubleMuTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/MuEG_Run2011B-PromptReco-v1_AOD/V04-02-30/SSignSkim/skim*root");

	//---------------------------
	// Run2011B prompt reco v1
	//---------------------------

	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleElectron_Run2011B-PromptReco-v1_AOD/V04-02-34/DoubleElectronTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/DoubleMu_Run2011B-PromptReco-v1_AOD/V04-02-34/DoubleMuTriggerSkim/skim*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/MuEG_Run2011B-PromptReco-v1_AOD/V04-02-34/SSignSkim/skim*root");


      }
      
      else if( trig == ossusy_looper::e_lowpt ){
	
	cout << "Doing dilepton-HT trigger data" << endl;

	//---------------------------
	// may10 rereco
	//---------------------------

	pickSkimIfExists(chdata,"cms2_data/ElectronHad_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skimmed*root");
	pickSkimIfExists(chdata,"cms2_data/MuHad_Run2011A-May10ReReco-v1_AOD/V04-02-20/SSignSkim/skimmed*root");

	//---------------------------
	// prompt reco v4
	//---------------------------

	pickSkimIfExists(chdata,"cms2_data/ElectronHad_Run2011A-PromptReco-v4_AOD/V04-02-20/SSignSkim/skimmed*root");
	pickSkimIfExists(chdata,"cms2_data/MuHad_Run2011A-PromptReco-v4_AOD/V04-02-20/SSignSkim/skimmed*root");

	//---------------------------
	// aug05 rereco
	//---------------------------

	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/ElectronHad_Run2011A-05Aug2011-v1_AOD/V04-02-29/SSignSkim/skimmed*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/MuHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SSignSkim/skimmed*root");

	//---------------------------
	// prompt reco v6
	//---------------------------

	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/ElectronHad_Run2011A-PromptReco-v6_AOD/V04-02-29/SSignSkim/skimmed*root");
	pickSkimIfExists(chdata,"/nfs-6/userdata/cms2/MuHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SSignSkim/skim*root");

	//---------------------------
	// 2011B
	//---------------------------

	pickSkimIfExists(chdata,"/hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");
	pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");



      }
	
    }

    for (int jetTypeIdx = 2; jetTypeIdx < 3; ++jetTypeIdx)
      {
	for (int metTypeIdx = 3; metTypeIdx < 4; ++metTypeIdx)
	  {
	    for (int zvetoIdx = 0; zvetoIdx < 1; ++zvetoIdx)
	      {
		for (int frmodeIdx = 0; frmodeIdx < (2-(1*!doFakeApp)); ++frmodeIdx)
		//for (int frmodeIdx = 0; frmodeIdx < 1; ++frmodeIdx)
		  {
                  
		    ossusy_looper::JetTypeEnum  jetType(jetTypeIdx);
		    ossusy_looper::MetTypeEnum  metType(metTypeIdx);
		    ossusy_looper::ZVetoEnum    zveto(zvetoIdx);
		    ossusy_looper::FREnum       frmode(frmodeIdx);

		    if( doFakeApp ){
		      if( frmodeIdx == 0 ) cout << "Doing double fake estimate" << endl;
		      if( frmodeIdx == 1 ) cout << "Doing single fake estimate" << endl;
		    }

		    if (rundataskim) {
		      cout << "Processing data skim" << endl;
		      looper->ScanChain(chdataskim,"dataskim", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing data skim" << endl;
		      hist::color("dataskim", kBlack);
		    }            
		    if (rundata) {
		      cout << "Processing data" << endl;
		      looper->ScanChain(chdata,"data", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing data" << endl;
		      hist::color("data", kBlack);
		    }
		    if (rundata41) {
		      cout << "Processing data 41X" << endl;
		      looper->ScanChain(chdata41,"data41", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing data 41X" << endl;
		      hist::color("data41", kBlack);
		    }
		    if (runttall) {
		      cout << "Processing ttbar all.. " << endl;
		      looper->ScanChain(chtopall,"ttall", kttall, prettall, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar all.. " << endl;
		      hist::color("ttall", kYellow);
		    }
		    if (runttallPUS6) {
		      cout << "Processing ttbar all PUS6.. " << endl;
		      looper->ScanChain(chtopallPUS6,"ttallPUS6", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar all PUS6.. " << endl;
		    }
		    if (runDYtot) {
		      cout << "Processing DY->all" << endl;
		      looper->ScanChain(chDYtot,"DYtot", kDYtot, preDYtot, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done rocessing DY->ee" << endl;
		      hist::color("DYtot", kMagenta);
		    }
		    if (runDYee) {
		      cout << "Processing DY->ee" << endl;
		      looper->ScanChain(chDYee,"DYee", kDYee, preDYee, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done rocessing DY->ee" << endl;
		      hist::color("DYee", kMagenta);
		    }
		    if (runDYmm) {
		      cout << "Processing DY->mm" << endl;
		      looper->ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing DY->mm" << endl;
		      hist::color("DYmm", kCyan);
		    }
		    if (runDYtautau) {
		      cout << "Processing DY->tautau" << endl;
		      looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing DY->tautau" << endl;
		      hist::color("DYtautau", kBlack);
		    }
		    if (runZjets) {
		      cout << "Processing Zjets" << endl;
		      looper->ScanChain(chZjets,"Zjets", kZjets, preZjets, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing Zjets" << endl;
		      hist::color("Zjets", kBlack);
		    }
		    if (runQCD) {
		      cout << "Processing QCD.. " << endl;
		      looper->ScanChain(chQCD,"qcd", kqcd, preqcd, lumi, jetType, metType, zveto,frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing  QCD.. " << endl;
		      hist::color("qcd", kOrange);
		    }
		    if (runphotons) {
		      cout << "Processing photons.. " << endl;
		      looper->ScanChain(chphotons,"photons", 1, 1, lumi, jetType, metType, zveto,frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing  photons.. " << endl;
		    }
		    if (runZZZ) {
		      cout << "Processing ZZZ.. " << endl;
		      looper->ScanChain(chZZZ,"ZZZ", 1, 1, lumi, jetType, metType, zveto,frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing  ZZZ.. " << endl;
		    }
		    if (runQCDpt15) {
		      cout << "Processing QCDpt15.. " << endl;
		      looper->ScanChain(chQCDpt15,"qcdpt15", kqcdpt15, preqcdpt15, lumi, jetType, metType, zveto,frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing  QCDpt15.. " << endl;
		      hist::color("qcdpt15", kOrange);
		    }
		    if (runQCDpt30) {
		      cout << "Processing QCDpt30.. " << endl;
		      looper->ScanChain(chQCDpt30,"qcdpt30", kqcdpt30, preqcdpt30, lumi, jetType, metType, zveto,frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing  QCDpt30.. " << endl;
		      hist::color("qcdpt30", kOrange);
		    }
		    if (runtt42) {
		      cout << "Processing ttbar 42.. " << endl;
		      looper->ScanChain(chtop42,"tt42", ktt42, prett42, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar 42.. " << endl;
		      hist::color("tt42", kYellow);
		    }
		    if (runttpowheg) {
		      cout << "Processing ttbar powheg.. " << endl;
		      looper->ScanChain(chtoppowheg,"ttpowheg", kttpowheg, prettpowheg, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar powheg.. " << endl;
		      hist::color("ttpowheg", kYellow);
		    }
		    if (runttdil) {
		      cout << "Processing ttbar dileptonic.. " << endl;
		      looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar dileptonic.. " << endl;
		      hist::color("ttdil", kYellow);
		    }
		    if (runttrelval) {
		      cout << "Processing ttbar relval.. " << endl;
		      looper->ScanChain(chttrelval,"ttrelval", kttrelval, prettrelval, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar dileptonic.. " << endl;
		      hist::color("ttrelval", kYellow);
		    }
		    if (runttem) {
		      cout << "Processing ttbar em.. " << endl;
		      looper->ScanChain(chtopem,"ttem", kttem, prettem, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar em.. " << endl;
		    }
		    if (runttotr) {
		      cout << "Processing ttbar no-dileptons.. " << endl;
		      looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ttbar no-dileptons.. " << endl;
		      hist::color("ttotr", 30);
		    }
		    if (runVV) {
		      cout << "Processing VV.." << endl;
		      looper->ScanChain(chvv,"vv", kVV, preVV, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing VV.." << endl;
		      hist::color("vv", kRed);
		    }
		    if (runWW) {
		      cout << "Processing WW.." << endl;
		      looper->ScanChain(chww,"ww", kWW, preWW, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing WW.." << endl;
		      hist::color("ww", kRed);
		    }
		    if (runWZ) {
		      cout << "Processing WZ.." << endl;
		      looper->ScanChain(chWZ,"wz", kWZ, preWZ, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing WZ.." << endl;
		      hist::color("wz", kBlue);
		    }
		    if (runZZ) {
		      cout << "Processing ZZ.." << endl;
		      looper->ScanChain(chZZ,"zz", kZZ, preZZ, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ZZ.." << endl;
		      hist::color("zz", kGreen);
		    }
		    if (runWjets) {
		      cout << "Processing Wjets.." << endl;
		      looper->ScanChain(chWjets,"wjets", kWjets, preWjets, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing Wjets.." << endl;
		      hist::color("wjets", 40);
		    }
		    if (runWjetsMG) {
		      cout << "Processing Wjets MG.." << endl;
		      looper->ScanChain(chWjetsMG,"wjetsMG", kWjetsMG, preWjetsMG, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing WjetsMG.." << endl;
		      hist::color("wjetsMG", 40);
		    }
		    if (runWcharm) {
		      cout << "Processing Wcharm.." << endl;
		      looper->ScanChain(chWcharm, "wcharm", kWcharm, preWcharm, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing Wcharm.." << endl;
		      hist::color("wcharm", 50);
		    }

		    if (runppMuX) {
		      cout << "Processing ppMuX" << endl;
		      looper->ScanChain(chppMuX,"ppMuX", kppMuX, preppMuX, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ppMuX" << endl;
		      hist::color("ppMuX", 51);
		    }
		    if (runEM) {
		      cout << "Processing EM" << endl;
		      looper->ScanChain(chEM,"EM", kEM, preEM, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing EM" << endl;
		      hist::color("EM", 49);
		    }
		    if (runtW) {
		      cout << "Processing tW" << endl;
		      looper->ScanChain(chtW,"tW", ktW, pretW, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing tW" << endl;
		      hist::color("tW", 63);
		    }
		    if (runVQQ) { 
		      cout << "Processing VQQ" << endl;
		      looper->ScanChain(chVQQ,"VQQ", kVQQ, preVQQ, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing VQQ" << endl;
		      hist::color("VQQ", 45);
		    }
		    if (runLM0) {
		      cout << "Processing LM0" << endl;
		      looper->ScanChain(chLM0, "LM0", kLM0, preLM0, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM0" << endl;
		      hist::color("LM0", kOrange);
		    }
		    if (runLM1) {
		      cout << "Processing LM1" << endl;
		      looper->ScanChain(chLM1, "LM1", kLM1, preLM1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM1" << endl;
		      hist::color("LM1", kOrange+1);
		    }
		    if (runLM1v2) {
		      cout << "Processing LM1v2" << endl;
		      looper->ScanChain(chLM1v2, "LM1v2", kLM1, preLM1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM1v2" << endl;
		    }
		    if (runLM2) {
		      cout << "Processing LM2" << endl;
		      looper->ScanChain(chLM2, "LM2", kLM2, preLM2, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM2" << endl;
		      hist::color("LM2", kOrange+2);
		    }
		    if (runLM3) {
		      cout << "Processing LM3" << endl;
		      looper->ScanChain(chLM3, "LM3", kLM3, preLM3, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM3" << endl;
		      hist::color("LM3", kOrange+3);
		    }
		    if (runLM3v2) {
		      cout << "Processing LM3v2" << endl;
		      looper->ScanChain(chLM3v2, "LM3v2", kLM3, preLM3, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM3v2" << endl;
		    }
		    if (runLM4) {
		      cout << "Processing LM4" << endl;
		      looper->ScanChain(chLM4, "LM4", kLM4, preLM4, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM4" << endl;
		      hist::color("LM4", kOrange+4);
		    }
		    if (runLM5) {
		      cout << "Processing LM5" << endl;
		      looper->ScanChain(chLM5, "LM5", kLM5, preLM5, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM5" << endl;
		      hist::color("LM5", kOrange+5);
		    }
		    if (runLM6) {
		      cout << "Processing LM6" << endl;
		      looper->ScanChain(chLM6, "LM6", kLM6, preLM6, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM6" << endl;
		      hist::color("LM6", kOrange+6);
		    }
		    if (runLM6v2) {
		      cout << "Processing LM6v2" << endl;
		      looper->ScanChain(chLM6v2, "LM6v2", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM6v2" << endl;
		    }
		    if (runLM7) {
		      cout << "Processing LM7" << endl;
		      looper->ScanChain(chLM7, "LM7", kLM7, preLM7, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM7" << endl;
		      hist::color("LM7", kOrange+7);
		    }
		    if (runLM8) {
		      cout << "Processing LM8" << endl;
		      looper->ScanChain(chLM8, "LM8", kLM8, preLM8, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM8" << endl;
		      hist::color("LM8", kOrange+8);
		    }
		    if (runLM9) {
		      cout << "Processing LM9" << endl;
		      looper->ScanChain(chLM9, "LM9", kLM9, preLM9, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM9" << endl;
		      hist::color("LM9", kOrange+9);
		    }
		    if (runLM10) {
		      cout << "Processing LM10" << endl;
		      looper->ScanChain(chLM10, "LM10", kLM10, preLM10, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM10" << endl;
		      hist::color("LM10", kOrange+10);
		    }
		    if (runLM11) { 
		      cout << "Processing LM11" << endl;
		      looper->ScanChain(chLM11, "LM11", kLM11, preLM11, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM11" << endl;
		      hist::color("LM11", kOrange-7);
		    }
		    if (runLM12) {
		      cout << "Processing LM12" << endl;
		      looper->ScanChain(chLM12, "LM12", kLM12, preLM12, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM12" << endl;
		      hist::color("LM12", kOrange-7);
		    }
		    if (runLM13) {
		      cout << "Processing LM13" << endl;
		      looper->ScanChain(chLM13, "LM13", kLM13, preLM13, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LM13" << endl;
		      hist::color("LM13", kOrange-7);
		    }
		    if (runML1) {
		      cout << "Processing ML1" << endl;
		      looper->ScanChain(chML1, "ML1", kML1, preML1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML1" << endl;
		    }
		    if (runML2) {
		      cout << "Processing ML2" << endl;
		      looper->ScanChain(chML2, "ML2", kML2, preML2, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML2" << endl;
		    }
		    if (runML3) {
		      cout << "Processing ML3" << endl;
		      looper->ScanChain(chML3, "ML3", kML3, preML3, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML3" << endl;
		    }
		    if (runML4) {
		      cout << "Processing ML4" << endl;
		      looper->ScanChain(chML4, "ML4", kML4, preML4, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML4" << endl;
		    }
		    if (runML5) {
		      cout << "Processing ML5" << endl;
		      looper->ScanChain(chML5, "ML5", kML5, preML5, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML5" << endl;
		    }
		    if (runML6) {
		      cout << "Processing ML6" << endl;
		      looper->ScanChain(chML6, "ML6", kML6, preML6, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML6" << endl;
		    }
		    if (runML7) {
		      cout << "Processing ML7" << endl;
		      looper->ScanChain(chML7, "ML7", kML7, preML7, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML7" << endl;
		    }
		    if (runML8) {
		      cout << "Processing ML8" << endl;
		      looper->ScanChain(chML8, "ML8", kML8, preML8, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing ML8" << endl;
		    }
		    if (runLMscan) {
		      cout << "Processing LMscan" << endl;
		      looper->ScanChain(chLMscan, "LMscan", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscan" << endl;
		      hist::color("LMscan", kOrange-7);
		    }
		    if (runLMscanFall11) {
		      cout << "Processing LMscanFall11" << endl;
		      looper->ScanChain(chLMscanFall11, "LMscanFall11", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11" << endl;
		    }
		    if (runLMscanFall11dil) {
		      cout << "Processing LMscanFall11 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil, "LMscanFall11dil", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil" << endl;
		    }
		    if (runLMscanFall11dil1) {
		      cout << "Processing LMscanFall111 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil1, "LMscanFall11dil1", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil1" << endl;
		    }
		    if (runLMscanFall11dil2) {
		      cout << "Processing LMscanFall112 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil2, "LMscanFall11dil2", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil2" << endl;
		    }
		    if (runLMscanFall11dil3) {
		      cout << "Processing LMscanFall113 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil3, "LMscanFall11dil3", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil3" << endl;
		    }
		    if (runLMscanFall11dil4) {
		      cout << "Processing LMscanFall114 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil4, "LMscanFall11dil4", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil4" << endl;
		    }
		    if (runLMscanFall11dil5) {
		      cout << "Processing LMscanFall115 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil5, "LMscanFall11dil5", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil5" << endl;
		    }
		    if (runLMscanFall11dil6) {
		      cout << "Processing LMscanFall116 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil6, "LMscanFall11dil6", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil6" << endl;
		    }
		    if (runLMscanFall11dil7) {
		      cout << "Processing LMscanFall117 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil7, "LMscanFall11dil7", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil7" << endl;
		    }
		    if (runLMscanFall11dil8) {
		      cout << "Processing LMscanFall118 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil8, "LMscanFall11dil8", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil8" << endl;
		    }
		    if (runLMscanFall11dil9) {
		      cout << "Processing LMscanFall119 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil9, "LMscanFall11dil9", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil9" << endl;
		    }
		    if (runLMscanFall11dil10) {
		      cout << "Processing LMscanFall1110 dilepton filter" << endl;
		      looper->ScanChain(chLMscanFall11dil10, "LMscanFall11dil10", kLMscan, preLMscan, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing LMscanFall11dil10" << endl;
		    }
		    if (runT1lh) {
		      cout << "Processing T1lh" << endl;
		      looper->ScanChain(chT1lh, "T1lh", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing T1lh" << endl;
		      hist::color("LMscan", kOrange-7);
		    }
		    if (runT2tt) {
		      cout << "Processing T2tt" << endl;
		      looper->ScanChain(chT2tt, "T2tt", 1, 1, lumi, jetType, metType, zveto, frmode, doFakeApp, calculateTCMET);
		      cout << "Done processing T2tt" << endl;
		    }

		    char* dir = "";
		    if( trig == ossusy_looper::e_lowpt  ) dir = "lowpt"  ;
		    if( trig == ossusy_looper::e_highpt ) dir = "highpt" ;
                  
		    // save all the histograms
		    if(doFakeApp) {
		      const char* outFile = Form("../output/%s/%s/ossusy_%s_%s%s_%s_FakeApp.root", version,dir,
						 jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx],zvetoStrings[zvetoIdx],frmodeStrings[frmode]);
		    }
		    else {
		      const char* outFile = Form("../output/%s/%s/ossusy_%s_%s%s.root", version,dir,
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
  }

  gSystem->Exit(0);
  
}
