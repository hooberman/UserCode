{
  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/utilities.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/SimpleFakeRate.cc+");
  gROOT->ProcessLine(".L CORE/mcSelections.cc+");
  gROOT->ProcessLine(".L CORE/MT2/MT2.cc+");
  gROOT->ProcessLine(".L CORE/triggerUtils.cc+");  
  gROOT->ProcessLine(".L CORE/susySelections.cc+");
  gROOT->ProcessLine(".L CORE/mcSUSYkfactor.cc+");
  gROOT->ProcessLine(".L CORE/triggerSuperModel.cc+");
  gROOT->ProcessLine(".L CORE/triggerUtils.cc+");
  gROOT->ProcessLine(".L CORE/ttbarSelections.cc+");
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runLooper.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
  
  //choose samples to run over----------------------
  runLooper("zmm");
  
}
