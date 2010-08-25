#include "TChain.h"
#include "looper.C"

void runLooper(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if( strcmp( prefix , "data" )  == 0 ){  
    ch->Add("/tas/benhoob/skims/cms2/ExpressPhysics_Run2010A-Express-v4_dilepskim.root");
    isData = true;
  } 

  else if( strcmp( prefix , "ZJets" )  == 0 ){  
    ch->Add("/tas/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root");
  }
   
  else if( strcmp( prefix , "TTBar" )  == 0 ){  
    ch->Add("/tas/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*root");
  }
   else if( strcmp( prefix , "Vqq" )  == 0 ){  
    ch->Add("/tas/cms2/VqqJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/merged*root");
  } 
 
  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  looper* mylooper = new looper();
  
  cout << "Running on sample " << prefix << endl;
  mylooper->ScanChain(ch, prefix, isData);
  
}





