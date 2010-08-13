{

#include "looper.h"

  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L CMS2.cc+");
  gROOT->ProcessLine(".L runLooper.C+");
  
  gSystem->Load("/tas03/home/benhoob/CMSSW_3_6_1_patch4/src/CMS2/NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");

  runLooper("data",true);
  runLooper("zmm");
  runLooper("zee");
  runLooper("qcd"); 
  runLooper("ttbar");

}
