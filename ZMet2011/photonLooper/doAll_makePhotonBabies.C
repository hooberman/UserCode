{

#include "makePhotonBabies.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runMakePhotonBabies.C+");
  
  //--photon+jets---
  //runMakePhotonBabies("PhotonJet"             , false );   //PhotonJet MC
  //runMakePhotonBabies("EG"                    , true  );   //EG data
  runMakePhotonBabies("Photon"                , true  );   //Photon data
}
