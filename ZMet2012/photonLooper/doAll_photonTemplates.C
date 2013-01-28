{

  gROOT->ProcessLine(".L histtools.C+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runPhotonTemplates.C+");
  
  //runPhotonTemplates("V00-00-07","Photon");
  runPhotonTemplates("V00-00-07","DoubleElectron");
  
}
