{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  //char* iter = "V01-03_copy";
  //char* iter = "nov5th_v4";
  char* iter = "v7";
  
  //runBabyLooper(iter,"EG"             , true  , babylooper::e_photonSelection , true );
  //runBabyLooper(iter,"PhotonJet"      , false , babylooper::e_photonSelection , true  );
  
  runBabyLooper(iter,"lepdata_skim"        , true  , babylooper::e_ZSelection  );
  //   //runBabyLooper(iter,"lepdata"        , true  , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"TTbar"          , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"ZJets"          , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"WJets"          , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"WW"             , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"WZ"             , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"ZZ"             , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"tW"             , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"LM4"            , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"LM8"            , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"LM0"            , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"LM1"            , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"LM2"            , false , babylooper::e_ZSelection  );
//   runBabyLooper(iter,"LM3"            , false , babylooper::e_ZSelection  );
}
