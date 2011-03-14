{

  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L CORE/susySelections.cc+");
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runLooper.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  
  //choose samples to run over----------------------
  runLooper("dyee");
  runLooper("dymm");
  runLooper("h130");
  //runLooper("dyee_nopu");
  //runLooper("dymm_nopu");
  //runLooper("h130_nopu");
  //runLooper("test384");
  //runLooper("test412");
  
}
