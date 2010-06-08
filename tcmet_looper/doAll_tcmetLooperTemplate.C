{

  gROOT->ProcessLine(".L histtools.C++");
  gROOT->ProcessLine(".L runTCMETLooperTemplate.C++");
  gSystem->Load("/tas03/home/benhoob/CMSSW_3_6_0_V03-04-09-01/src/CMS2/NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");

  //runTCMETLooperTemplate("data");
  //runTCMETLooperTemplate("zmm");
  runTCMETLooperTemplate("zee");
  //runTCMETLooperTemplate("datahighmet");
  //runTCMETLooperTemplate("qcd");
  //runTCMETLooperTemplate("ttbar");
  
}
