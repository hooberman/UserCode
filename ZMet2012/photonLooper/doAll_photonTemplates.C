{

  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runPhotonTemplates.C+");
  
  //runPhotonTemplates("V00-00-11","Photon");
  runPhotonTemplates("V00-00-11","DoubleElectron");
  
}
