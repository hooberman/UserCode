{

#include "looper.h"

gROOT->ProcessLine(".L histtools.C+");
gROOT->ProcessLine(".L CORE/CMS2.cc+");
gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
gROOT->ProcessLine(".L runLooper.C+");

//void runLooper(char* prefix , bool isData = true, metAlgo algo){

//runLooper("JetMETTau_250nb"                 , true  , looper::e_makeTemplate     );
//runLooper("QCD_Pt15"                        , false , looper::e_makeTemplate  );
runLooper("EG"                              , true  , looper::e_photonSelection  );
//runLooper("PhotonJet_Pt15"                  , false , looper::e_photonSelection  );
//runLooper("ZJets"                           , false , looper::e_ZSelection  );

}
