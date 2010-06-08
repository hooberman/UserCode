{

  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  //gROOT->ProcessLine(".L Tools/goodrun.cc+");
  //gROOT->ProcessLine(".L CORE/utilities.cc+");
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runTCMETLooperTemplate.C+");
  gSystem->Load("libMiniFWLite.so");





  //choose samples to run over----------------------

  //runTCMETLooperTemplate("data");
  //runTCMETLooperTemplate("zmm");
  //runTCMETLooperTemplate("zee");
  runTCMETLooperTemplate("datahighmet");
  //runTCMETLooperTemplate("qcd");
  //runTCMETLooperTemplate("ttbar");
  
  //------------------------------------------------
  
}
