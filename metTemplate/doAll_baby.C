{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  //runBabyLooper("V01-02","EG"             , true  , babylooper::e_photonSelection , true );
  //runBabyLooper("V01-02","PhotonJet"      , false , babylooper::e_photonSelection , true  );
  
  //runBabyLooper("V01-02","lepdata"        , true  , babylooper::e_ZSelection  );
  //runBabyLooper("V01-02","TTbar"          , false , babylooper::e_ZSelection  );
  runBabyLooper("V01-02","ZJets"          , false , babylooper::e_ZSelection  );
  //runBabyLooper("V01-02","LM4"            , false , babylooper::e_ZSelection  );

}
