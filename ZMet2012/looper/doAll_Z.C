{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  // runZLooper("dataskim2010"            , true  );
  // runZLooper("data2012cv2"             , true  );
  // runZLooper("ttbar_massiveb"          , false );  
  // runZLooper("ww"                      , false );  
  // runZLooper("t"                       , false );  
  // runZLooper("wz"                      , false );  
  // runZLooper("zz"                      , false );  
  // runZLooper("ttbar"                   , false );  
  // runZLooper("zjets"                   , false );  
  // runZLooper("zjets_10to50"            , false );  

  // runZLooper("data"                    , true  );
  // runZLooper("dataskim"                , true  );
  // runZLooper("data2012c"               , true  );
  runZLooper("wzsms"                   , false );  
  // runZLooper("gmsb"                    , false );  
  // runZLooper("testfilter_newJEC"       , true  );  

  //--------------------------------------------------
  // samples for sync exercise
  //--------------------------------------------------

  // runZLooper("RelValZEE"                  , false );  
  // runZLooper("RelValZMM"                  , false );
  // runZLooper("MuEG_199752"                , true  );
  // runZLooper("DoubleElectron_199752"      , true  );
  // runZLooper("DoubleMu_199752"            , true  );


}




