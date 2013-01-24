int doSample_Z(std::string sample){

  gROOT->ProcessLine(".L histtools.C+");
  //gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L runZLooper.C+");

  bool isData = TString(sample).Contains("data");
  runZLooper(sample.c_str(),isData);

  return 0;
}




