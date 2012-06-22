{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  char* Z_version        = "V00-00-16";
  char* template_version = "V00-00-11";
  
  runBabyLooper(Z_version,template_version, "dataskim"  , true   , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "wz"      , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zz"        , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "wz"        , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zjets"     , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "LM4"       , false  , babylooper::e_ZSelection  );

  
}

