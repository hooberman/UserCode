

{

#include "looper.h"

gROOT->ProcessLine(".L histtools.C+");
gROOT->ProcessLine(".L CORE/CMS2.cc+");
gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
gROOT->ProcessLine(".L runLooper.C+");

//runLooper("may5PDSkim2_SD_EG_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTau_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTauMonitor_goodrunPfJetPt15");

//void runLooper(char* prefix , bool isData = true, metAlgo algo){

runLooper("JetMETTau_250nb"                 , true  , looper::e_makeTemplate     );
runLooper("QCD_Pt15"                        , false , looper::e_makeTemplate  );
runLooper("EG"                              , true  , looper::e_photonSelection  );
runLooper("PhotonJet_Pt15"                  , false , looper::e_photonSelection  );
runLooper("ZJets"                           , false , looper::e_ZSelection  );


}
