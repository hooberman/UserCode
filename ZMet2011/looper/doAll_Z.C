{
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/utilities.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
  gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/mcSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/susySelections.cc+");
  //gROOT->ProcessLine(".L ../CORE/jetSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/mcSUSYkfactor.cc+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");
  
  runZLooper("wzsms"                 , false );
  runZLooper("zzsms"                 , false );
  runZLooper("ggmsb"                 , false );
  //runZLooper("data"                  , true  );
  // runZLooper("ttbar"                 , false );
  // runZLooper("zjets"                 , false );
  // runZLooper("zjetsS6_incomplete"    , false );
  // runZLooper("dyee"                  , false );
  // runZLooper("dymm"                  , false );
  // runZLooper("wz_summer11_madgraph"  , false );
  // runZLooper("zz_summer11_madgraph"  , false );
  // runZLooper("LM4"                   , false );
  // runZLooper("LM4v2"                 , false );
  // runZLooper("LM8"                   , false );
  // runZLooper("LM4v2"                 , false );
  // runZLooper("LM8v2"                 , false );
  // runZLooper("LM9"                   , false );
  // runZLooper("singletop"             , false );
  // runZLooper("T5zzgmsb"              , false );
  // runZLooper("T5zzgmsb_hadoop"       , false );
  //runZLooper("T5zz"                  , false );
  // runZLooper("ZZZ"                   , false );
  // runZLooper("T5zzh"                 , false );
  // runZLooper("T5zzl"                 , false );
  // runZLooper("wz_summer11_pythia"    , false );
  // runZLooper("zz_summer11_pythia"    , false );

}




