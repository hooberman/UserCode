{

gROOT->ProcessLine(".L histtools.C++");
gROOT->ProcessLine(".L runLooper.C++");

gSystem->Load("/tas03/home/benhoob/CMSSW_3_3_6/src/metTemplate/Tools/MiniFWLite/libMiniFWLite.so");

//runLooper("may5PDSkim2_SD_EG_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTau_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTauMonitor_goodrunPfJetPt15");

//void runLooper(char* prefix , bool isData = true, bool makeMetTemplate = false){
//runLooper("JetMETTau"                       , true  , true  );
runLooper("EG"                              , true  , false );
runLooper("PhotonJet_Pt15"                  , false , false );


}
