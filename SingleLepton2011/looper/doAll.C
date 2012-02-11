
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
  
  const char* version    = "V00-03-00";
  const char* jsonfile   = "jsons/Cert_160404-180252_7TeV_mergePromptMay10Aug5_JSON_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable(); 

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
  bool rundata     = 1;
  bool runttall    = 0;
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
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
  }

  //----------------------------------------
  // ttbar
  //----------------------------------------

  TChain* chtopall = new TChain("Events");
  if (runttall) {
    pickSkimIfExists(chtopall,"/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged*root");
  }

  //----------------------------------------
  // W+jets
  //----------------------------------------

  TChain* chWjets = new  TChain("Events");
  if(runWjets){
    pickSkimIfExists(chWjets,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndJets/merged*root");
  }

  //----------------------------------------
  // DY
  //----------------------------------------

  TChain* chDYtot = new  TChain("Events");
  if(runDYtot){
    string dypath = "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/"
    pickSkimIfExists(chDYtot,dypath+"DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    pickSkimIfExists(chDYtot,dypath+"DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    pickSkimIfExists(chDYtot,dypath+"DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    pickSkimIfExists(chDYtot,dypath+"DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    pickSkimIfExists(chDYtot,dypath+"DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    pickSkimIfExists(chDYtot,dypath+"DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    pickSkimIfExists(chDYtot,dypath+"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
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
  // DATA
  //----------------------------------------

  TChain* chdata     = new  TChain("Events");

  if(rundata){
    
    cout << "adding ElectronHad and MuHad or SingleMu data" << endl;

    //pickSkimIfExists(chdata,"/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/SingleMu_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_999999_0_skim.root");//l+2j filtered
    
    // May10
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-May10ReReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");   
    // PRv4
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v4_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");

    // Aug05
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");

    // PRv6
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");

    // 2011B
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/merged*root");

    // 2011B
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/ElectronHad_Run2011B-PromptReco-v1_AOD/V04-02-34/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/nfs-7/userdata/cms2/MuHad_Run2011B-PromptReco-v1_AOD/V04-02-35/SingleLeptonAndJets/merged*root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
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
