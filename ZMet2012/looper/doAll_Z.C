{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  // runZLooper("dataskim"                , true  );
  // runZLooper("wz"                      , false );  
  // runZLooper("zz"                      , false );  
  // runZLooper("ttbar"                   , false );  
  // runZLooper("ww"                      , false );  
  // runZLooper("t"                       , false );  
  // runZLooper("zjets"                   , false );  
  runZLooper("wzsms"     , false  );  
  runZLooper("gmsb"      , false  );  

  //runZLooper("testfilter_newJEC"     , true  );  
  //runZLooper("RelValZEE"             , false );  
  //runZLooper("RelValZMM"             , false );  
  //runZLooper("data"                  , true  );
}




