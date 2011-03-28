#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;

  char* path = "/tas/cms2";
  //char* path = "/nfs-3/userdata/cms2";
  
  if( strcmp( prefix , "dyee" ) == 0 ){
    ch->Add(Form("%s/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged_ntuple*root",path));
  }

  else if( strcmp( prefix , "dymm" ) == 0 ){
    ch->Add(Form("%s/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged*root",path));
  }
  
  else if( strcmp( prefix , "h130" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged*root",path));
  }

  else if( strcmp( prefix , "tt" ) == 0 ){
    ch->Add(Form("%s/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged*root",path));
  }

  else if( strcmp( prefix , "dyee_spring11" ) == 0 ){
    //ch->Add("/tas/benhoob/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
    ch->Add("/tas/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/merged*root");
  }

  else if( strcmp( prefix , "dymm_spring11" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
  }

  else if( strcmp( prefix , "tt_spring11" ) == 0 ){
    //ch->Add("/tas/benhoob/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
    ch->Add("/tas/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-00/merged*root");
  }

  else if( strcmp( prefix , "h130_spring11" ) == 0 ){
    //ch->Add("/tas/benhoob/cms2/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-00-10/n*root");
    ch->Add("/tas/cms2/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-10/merged*root");
  }

  else if( strcmp( prefix , "dyee_winter10" ) == 0 ){
    ch->Add("/tas/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged*root");
  }

  else if( strcmp( prefix , "dymm_winter10" ) == 0 ){
    ch->Add("/tas/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V04-00-05-03/merged*root");
  }

  else if( strcmp( prefix , "data" ) == 0 ){
    //ch->Add(Form("%s/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/ntuple.root",path));
    //ch->Add(Form("%s/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/ntuple_1.root",path));
    ch->Add(Form("%s/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/ntuple*root",path));
    isData = true;
  }

  else if( strcmp( prefix , "dyee_373" ) == 0 ){
    ch->Add("/tas/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V03-07-03_PFCandidates/merged_ntuple*root");
  }

  else if( strcmp( prefix , "dymm_373" ) == 0 ){
    ch->Add("/tas/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V03-07-03_PFCandidates/merged_ntuple*root");
  }
  
  else if( strcmp( prefix , "h130_373" ) == 0 ){
    ch->Add("/tas/cms2/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Winter10-E7TeV_ProbDist_2011Flat_BX156_START39_V8-v1/V03-07-03_PFCandidates/merged*root");
  }

  else if( strcmp( prefix , "dyee_nopu" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-06-18_PFCandidates/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-START38_V12-v1/ntuple*root");
  }

  else if( strcmp( prefix , "dymm_nopu" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-06-18_PFCandidates/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall10-START38_V12-v1/ntuple*root");
  }
  
  else if( strcmp( prefix , "h130_nopu" ) == 0 ){
    ch->Add("/tas/benhoob/cms2/CMS2_V03-06-18_PFCandidates/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-START38_V12-v1/ntuple*root");
  }
  
  else if( strcmp( prefix , "test384" ) == 0 ){
    ch->Add("ntuple_V03_08_04.root");
  }
  
  else if( strcmp( prefix , "test412" ) == 0 ){
    ch->Add("test_3_11_1_patch2.root");
    isData = true;
  }
  
  else if( strcmp( prefix , "test" ) == 0 ){
    //ch->Add("ntuple_indicies.root");
    ch->Add("ntuple412_40004.root");
  }
  
  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  looper* mylooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  mylooper->ScanChain(ch, prefix, isData);
  
}





