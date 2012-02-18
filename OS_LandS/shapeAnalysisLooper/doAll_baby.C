{

#include "babylooper.h"
  
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runBabyLooper.C+");
  
  char* version        = "V00-02-12";
  
  runBabyLooper(version, "LMscanFall11dil" , false );
  
}
