{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  char* Z_version        = "V00-00-10";
  char* template_version = "V00-00-07";
  
  runBabyLooper(Z_version,template_version, "data"      , true   , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "ttbar"     , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zjets"     , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "LM4"       , false  , babylooper::e_ZSelection  );

  
}
