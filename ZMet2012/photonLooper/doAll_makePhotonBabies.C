{

#include "makePhotonBabies.h"


  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/utilities.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
  gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/susySelections.cc+");
  gROOT->ProcessLine(".L ../CORE/photonSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/jetSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/mcSUSYkfactor.cc+");
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runMakePhotonBabies.C+");
  
  //--photon+jets---
  runMakePhotonBabies("Photon"                , true  );   //Photon data
}
