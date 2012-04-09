{
  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  runZLooper("zjets"                  , false );  
  // runZLooper("data"                  , true  );

}




