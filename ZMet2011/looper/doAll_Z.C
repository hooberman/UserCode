{

#include "Z_looper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runZLooper.C+");
  
  //run on samples
  runZLooper("express_data"   , true ,  Z_looper::e_ZSelection );
  runZLooper("dyee_spring11"  , false , Z_looper::e_ZSelection );
  runZLooper("dymm_spring11"  , false , Z_looper::e_ZSelection );
  runZLooper("tt_spring11"    , false , Z_looper::e_ZSelection );
  runZLooper("zjets_spring11" , false , Z_looper::e_ZSelection );




//   runZLooper("DY"             , false , Z_looper::e_ZSelection );
//   runZLooper("WJets"          , false , Z_looper::e_ZSelection );
//   runZLooper("WW"             , false , Z_looper::e_ZSelection );
//   runZLooper("WZ"             , false , Z_looper::e_ZSelection );
//   runZLooper("ZZ"             , false , Z_looper::e_ZSelection );
//   runZLooper("tW"             , false , Z_looper::e_ZSelection );
//   runZLooper("LM0"            , false , Z_looper::e_ZSelection );
//   runZLooper("LM1"            , false , Z_looper::e_ZSelection );
//   runZLooper("LM2"            , false , Z_looper::e_ZSelection );
//   runZLooper("LM3"            , false , Z_looper::e_ZSelection );
//   runZLooper("LM4"            , false , Z_looper::e_ZSelection );
//   runZLooper("LM8"            , false , Z_looper::e_ZSelection );
//   runZLooper("LM9"            , false , Z_looper::e_ZSelection );
//   runZLooper("lepdata_skim"        , true ,  Z_looper::e_ZSelection );
//   runZLooper("lepdata_skim_nov4"   , true ,  Z_looper::e_ZSelection );
//   runZLooper("lepdata"        , true ,  Z_looper::e_ZSelection );
//   runZLooper("TTbar"          , false , Z_looper::e_ZSelection ); 
//   runZLooper("ZJets"          , false , Z_looper::e_ZSelection );
//   runZLooper("testdata"                          , true ,  Z_looper::e_ZSelection         );
  

}
