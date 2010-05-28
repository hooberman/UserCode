{

gROOT->ProcessLine(".L histtools.C++");
gROOT->ProcessLine(".L runLooper.C++");

runLooper("may5PDSkim2_SD_EG_goodrun");
//runLooper("may5PDSkim2_SD_JetMETTau_goodrun");
runLooper("may5PDSkim2_SD_JetMETTauMonitor_goodrunPfJetPt15");
//runLooper("Commissioning10-SD_JetMETTauMonitor-v9_goodrunPfJetPt15");
//runLooper("Commissioning10-SD_EG-v9_goodrunPhotonGT10GeV");
//runLooper("Commissioning10-SD_JetMETTau-v9_goodrunPfJetPt30");

}
