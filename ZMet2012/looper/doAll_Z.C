{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  //runZLooper("data"                  , true  );
  //runZLooper("zjets"                 , false );  
  //runZLooper("ttbar"                 , false );  
  //runZLooper("wzmg"                    , false );  
  //runZLooper("zz"                      , false );  
  //runZLooper("t"                       , false );  
  //runZLooper("ww"                      , false );  
  //runZLooper("testfilter_newJEC"         , true );  
  runZLooper("RelValZEE"         , false );  
  runZLooper("RelValZMM"         , false );  
}




