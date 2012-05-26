{
  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  //runZLooper("data"                  , true  );
  runZLooper("zjets"                 , false );  
  //runZLooper("ttbar"                 , false );  

}




