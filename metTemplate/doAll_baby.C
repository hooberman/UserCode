{

#include "babylooper.h"

gROOT->ProcessLine(".L histtools.C+");
//gROOT->ProcessLine(".L CORE/CMS2.cc+");
gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");
gROOT->ProcessLine(".L runBabyLooper.C+");

//void runBabyLooper(char* prefix , bool isData = true, metAlgo algo){

//runBabyLooper("JetMETTau_250nb"                 , true  , babylooper::e_makeTemplate     );
//runBabyLooper("QCD_Pt15"                        , false , babylooper::e_makeTemplate  );
//runBabyLooper("EG"                              , true  , babylooper::e_photonSelection  );
//runBabyLooper("EG_1p9pb"                        , true  , babylooper::e_photonSelection  );
//runBabyLooper("PhotonJet_Pt15"                  , false , babylooper::e_photonSelection  );
//runBabyLooper("PhotonJet_1p9pb"                 , false , babylooper::e_photonSelection  );
//runBabyLooper("ZJets"                           , false , babylooper::e_ZSelection  );
//runBabyLooper("ZJets_1p9pb"                     , false , babylooper::e_ZSelection  );
//runBabyLooper("lepdata_1p9pb"                   , true , babylooper::e_ZSelection  );

runBabyLooper("lepdata_2p8pb_ttdilsync_ttdilskims"        , true ,  babylooper::e_ZSelection  );
runBabyLooper("ZJets_2p8pb_ttdilsync_ttdilskims"          , false , babylooper::e_ZSelection  );
runBabyLooper("TTbar_2p8pb_ttdilsync_ttdilskims"          , false , babylooper::e_ZSelection  );

//runBabyLooper("makeTemplates_EG_ttdilsync"                 , true   , babylooper::e_photonSelection  , true);
//runBabyLooper("makeTemplates_PhotonJet_ttdilsync"           , false  , babylooper::e_photonSelection  , true);
}
