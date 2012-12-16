{
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runMakePhotonBabies.C+");

  runMakePhotonBabies("data_53X_2012A"             , true  );  
  // runMakePhotonBabies("data_53X_2012B"             , true  );  
  // runMakePhotonBabies("data_53X_2012C"             , true  );  
  // runMakePhotonBabies("data_53X_2012D"             , true  );  
  
  //--photon+jets---
  //runMakePhotonBabies("DoubleElectron"          , true  );   // DoubleElectron data
  //runMakePhotonBabies("DoubleElectron_2012Cv2"  , true  );   // DoubleElectron data
  //runMakePhotonBabies("Photon"                  , true  );   // Photon data
}
