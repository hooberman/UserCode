{

  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runPhotonTemplates.C+");
  
  runPhotonTemplates();
  
}
