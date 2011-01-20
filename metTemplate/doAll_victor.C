{

  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runVictorTemplates.C+");
  
  runVictorTemplates();
  
}
