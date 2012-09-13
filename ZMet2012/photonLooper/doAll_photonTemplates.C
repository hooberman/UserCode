{

  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runPhotonTemplates.C+");
  gSystem->Load("../../MiniFWLite/libMiniFWLite.so");
  
  //runPhotonTemplates("V00-00-11","Photon");
  runPhotonTemplates("V00-01-00","DoubleElectron");
  
}
