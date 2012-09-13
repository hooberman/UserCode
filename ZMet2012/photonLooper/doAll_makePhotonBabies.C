{
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runMakePhotonBabies.C+");
  
  //--photon+jets---
  //runMakePhotonBabies("DoubleElectron"        , true  );   // DoubleElectron data
  runMakePhotonBabies("DoubleElectron_2012Cv2"  , true  );   // DoubleElectron data
  //runMakePhotonBabies("Photon"                , true  );   // Photon data
}
