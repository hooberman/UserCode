

{

#include "looper.h"

gROOT->ProcessLine(".L histtools.C+");
gROOT->ProcessLine(".L runLooper.C+");
gROOT->ProcessLine(".L CORE/CMS2.cc+");

gSystem->Load("/tas03/home/benhoob/CMSSW_3_3_6/src/metTemplate/Tools/MiniFWLite/libMiniFWLite.so");

//runLooper("may5PDSkim2_SD_EG_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTau_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTauMonitor_goodrunPfJetPt15");

//void runLooper(char* prefix , bool isData = true, metAlgo algo){

//runLooper("JetMETTau"                       , true  , looper::e_makeTemplate     );
//runLooper("QCD_Pt15"                        , false , looper::e_makeTemplate  );
runLooper("EG"                              , true  , looper::e_photonSelection  );
//runLooper("PhotonJet_Pt15"                  , false , looper::e_photonSelection  );
//runLooper("ZJets"                           , false , looper::e_ZSelection  );


}
