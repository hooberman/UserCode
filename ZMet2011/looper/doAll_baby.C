{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  char* Z_version        = "V00-00-02";
  char* template_version = "V00-00-01";
  
  //runBabyLooper(Z_version,template_version"EG"             , true  , babylooper::e_photonSelection , true );
  //runBabyLooper(Z_version,template_version"PhotonJet"      , false , babylooper::e_photonSelection , true  );
  

  runBabyLooper(Z_version,template_version, "pr_data"   , true   , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "ttbar"     , false  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "zjets"     , false  , babylooper::e_ZSelection  );

  //runBabyLooper(Z_version,template_version"lepdata_skim"   , true  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"lepdata"        , true  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "express_data"   , true  , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "dyee_spring11"  , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "dymm_spring11"  , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version, "tt_spring11"    , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"TTbar"          , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"ZJets"          , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"WJets"          , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"WW"             , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"WZ"             , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"ZZ"             , false , babylooper::e_ZSelection  );
  //runBabyLooper(Z_version,template_version"tW"             , false , babylooper::e_ZSelection  );
  
}
