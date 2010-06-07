#include "TChain.h"
#include "tcmetLooperTemplate.C"

void runTCMETLooperTemplate(char* prefix){

  string mc_v26   = "/tas01/disk01/cms2/MinBias_Spring10-START3X_V26A_356ReReco-v1/V03-03-07/merged_ntuple*.root";
  string mc_we    = "/tas01/disk02/cms2/Wenu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root";
  string data_highmet = "/tas01/disk02/benhoob/ntuples/ntuple_picked_event_only_tcMET_gt_45GeV_7TeV_HighStats.root";
 
  string data_v8  = "/tas01/disk01/cms2/MinimumBias_Commissioning10-GOODCOLL-v8/merged_ntuple*.root";
  
  //Spring10 samples
  string mc_zee   = "/tas07/disk00/cms2/Zee_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root";
  string mc_zmm   = "/tas07/disk00/cms2/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*root";
  string mc_ttbar = "/tas07/disk00/cms2/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*root";
  string mc_qcd   = "/tas07/disk00/cms2/QCD_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*root";

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if     ( strcmp( prefix , "zee" )   == 0 )  ch->Add( mc_zee.c_str() );
  else if( strcmp( prefix , "zmm" )   == 0 )  ch->Add( mc_zmm.c_str() );
  else if( strcmp( prefix , "ttbar" ) == 0 )  ch->Add( mc_ttbar.c_str() );
  else if( strcmp( prefix , "qcd" )   == 0 )  ch->Add( mc_qcd.c_str() );
  else if( strcmp( prefix , "data" )  == 0 ){  
    ch->Add( data_v8.c_str() );
    isData = true;
  }
  else if( strcmp( prefix , "datahighmet" )  == 0 ){  
    ch->Add( data_highmet.c_str() );
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





