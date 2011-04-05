{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  char* Z_version        = "V00-00-02";
  char* template_version = "V00-00-01";
  
  runBabyLooper(Z_version,template_version, "pr_data"   , true   , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "ttbar"     , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zjets"     , false  , babylooper::e_ZSelection  );

  
}
