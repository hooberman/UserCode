{

  //gSystem->Load("../../MiniFWLite/libMiniFWLite.so");

  #include "babylooper_edge.h"

  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runBabyLooper_edge.C+");

  char* Z_version        = "V00-01-04";
  char* template_version = "V00-01-00";
  
  runBabyLooper_edge(Z_version,template_version, "data_ALL_53X" , true , babylooper_edge::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "wz"          , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zz"          , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "wz"          , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zjets"       , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "LM4"         , false  , babylooper::e_ZSelection  );

  
}

