{

#include "Z_looper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runZLooper.C+");
  
  //run on samples
  runZLooper("lepdata_skim"   , true ,  Z_looper::e_ZSelection          );
  runZLooper("lepdata"        , true ,  Z_looper::e_ZSelection          );
  runZLooper("ZJets"          , false , Z_looper::e_ZSelection , 1.27   );
  runZLooper("TTbar"          , false , Z_looper::e_ZSelection , 0.9545 );
  runZLooper("WJets"          , false , Z_looper::e_ZSelection , 1.116  );
  runZLooper("WW"             , false , Z_looper::e_ZSelection          );
  runZLooper("WZ"             , false , Z_looper::e_ZSelection          );
  runZLooper("ZZ"             , false , Z_looper::e_ZSelection          );
  runZLooper("tW"             , false , Z_looper::e_ZSelection          );

  //runZLooper("testdata"                          , true ,  Z_looper::e_ZSelection         );
  //runZLooper("LM4"                             , false , Z_looper::e_ZSelection , 1);

}
