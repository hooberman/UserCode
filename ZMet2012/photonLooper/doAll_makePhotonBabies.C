{
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runMakePhotonBabies.C+");
  
  //--photon+jets---
  runMakePhotonBabies("Photon"                , true  );   //Photon data
}
