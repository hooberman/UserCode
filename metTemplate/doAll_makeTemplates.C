{

#include "makeTemplates.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L runMakeTemplates.C+");
  
  //---make templates---
  //runMakeTemplates("JetMETTau"                       , true  , makeTemplates::e_QCDSelection );
  //runMakeTemplates("QCD_Pt15"                        , false , makeTemplates::e_QCDSelection );
  //runMakeTemplates("QCD_Pt30"                        , false , makeTemplates::e_QCDSelection );
  
  //--photon+jets---
  runMakeTemplates("PhotonJet"                       , false , makeTemplates::e_photonSelection );
  runMakeTemplates("EG"                              , true  , makeTemplates::e_photonSelection );  
}
