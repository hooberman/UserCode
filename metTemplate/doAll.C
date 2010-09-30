{

#include "looper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runLooper.C+");
  


  //---make templates---
  //runLooper("JetMETTau_250nb"                 , true  , looper::e_makeTemplate     );
  //runLooper("QCD_Pt15"                        , false , looper::e_makeTemplate  );
  
  //--photon+jets---
  //runLooper("EG"                              , true  , looper::e_photonSelection  );
  //runLooper("PhotonJet"                       , false , looper::e_photonSelection  );
  //runLooper("Wenu"                            , false , looper::e_photonSelection  );
  //runLooper("QCD_Pt15"                        , false , looper::e_photonSelection  );
  //runLooper("QCD_Pt30"                        , false , looper::e_photonSelection  );
  
  //---Z+jets---
  //runLooper("lepdata"                         , true ,  looper::e_ZSelection );
  //runLooper("ZJets"                           , false , looper::e_ZSelection  , 1.27);
  //runLooper("TTbar"                           , false , looper::e_ZSelection );

  runLooper("testdata"                        , true , looper::e_ZSelection );
}
