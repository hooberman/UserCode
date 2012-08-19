{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  runZLooper("ttbar_massiveb"          , false );  
  runZLooper("ww"                      , false );  
  runZLooper("t"                       , false );  
  runZLooper("wz"                      , false );  
  runZLooper("zz"                      , false );  
  runZLooper("ttbar"                   , false );  
  runZLooper("dataskim2010"            , true  );
  runZLooper("zjets"                   , false );  

  // runZLooper("data"                  , true  );
  // runZLooper("dataskim"              , true  );
  // runZLooper("data2012c"             , true  );
  // runZLooper("wzsms"                 , false );  
  // runZLooper("gmsb"                  , false );  
  // runZLooper("testfilter_newJEC"     , true  );  
  // runZLooper("RelValZEE"             , false );  
  // runZLooper("RelValZMM"             , false );  

}




