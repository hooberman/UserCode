{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  //runZLooper("dataskim"                , true  );
  runZLooper("zjets"                   , false );  
  runZLooper("ttbar"                   , false );  
  //runZLooper("wz"                      , false );  
  //runZLooper("zz"                      , false );  
  //runZLooper("ww"                      , false );  
  //runZLooper("t"                       , false );  

  //runZLooper("testfilter_newJEC"     , true  );  
  //runZLooper("RelValZEE"             , false );  
  //runZLooper("RelValZMM"             , false );  
  //runZLooper("data"                  , true  );
}




