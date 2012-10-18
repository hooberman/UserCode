{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  // runZLooper("data_53X"                , true  );  
  // runZLooper("data_2012C_53X"          , true  );  
  // runZLooper("data_53X_edgeSync"       , true  );  
  // runZLooper("wz_53X"                  , false );  
  // runZLooper("zz_53X"                  , false );  
  // runZLooper("gmsb"                    , false );  
  runZLooper("gmsb_526"                   , false );  
  // runZLooper("wzsms"                   , false );  
  // runZLooper("wjets_53X"               , false );  
  // runZLooper("ww_53X"                  , false );  
  // runZLooper("t_53X"                   , false );  
  // runZLooper("wz2l2q_53X"              , false );  
  // runZLooper("zz2l2q_53X"              , false );  

  // runZLooper("zz4l_53X"                , false );  
  // runZLooper("ttW_53X"                 , false );  
  // runZLooper("ttZ_53X"                 , false );  
  // runZLooper("VVV_53X"                 , false );
  // runZLooper("zjets_full_53X"          , false );  
  // runZLooper("ttbar_53X"               , false );  
  // runZLooper("zjets_53X"               , false );  
  // runZLooper("zjets_MET50_53X"         , false );  


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
  // runZLooper("wzsms"                   , false );  
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




