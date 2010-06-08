#include "TChain.h"
#include "tcmetLooperTemplate.C"

void runTCMETLooperTemplate(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if     ( strcmp( prefix , "zee" )   == 0 ){
    //ch->Add( "/tas07/disk00/cms2/Zee_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root" );
    ch->Add( "/tas/cms2/Zee_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root" );
  }

  else if( strcmp( prefix , "zmm" )   == 0 ){
    //ch->Add( "/tas07/disk00/cms2/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root" );
    ch->Add( "/tas/cms2/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root" );
  }

  else if( strcmp( prefix , "ttbar" ) == 0 ){
    //ch->Add( "/tas07/disk00/cms2/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*root" );
    ch->Add( "/tas/cms2/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*root" );
  }

  else if( strcmp( prefix , "qcd" )   == 0 ){
    //ch->Add( "/tas07/disk00/cms2/QCD_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*root" );
    ch->Add( "/tas/cms2/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*root" );
  }

  else if( strcmp( prefix , "data" )  == 0 ){  
    ch->Add( "/tas01/disk01/cms2/MinimumBias_Commissioning10-GOODCOLL-v8/merged_ntuple*.root" );
    isData = true;
  }

  else if( strcmp( prefix , "datahighmet" )  == 0 ){  
    //ch->Add( "/tas01/disk02/benhoob/ntuples/ntuple_picked_event_only_tcMET_gt_45GeV_7TeV_HighStats.root" );
    ch->Add( "/tas03/home/benhoob/CMSSW_3_6_0_V03-04-09-01/src/CMS2/NtupleMaker/test/ntuple_picked_events_tcMET_gt_45GeV_Artur.root");
    isData = true;
  }

  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  tcmetLooperTemplate* looper = new tcmetLooperTemplate();
  
  cout << "Running on sample " << prefix << endl;
  looper->ScanChain(ch, prefix, isData);
  
}





