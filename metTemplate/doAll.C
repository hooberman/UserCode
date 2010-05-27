{

gROOT->ProcessLine(".L histtools.C++");
gROOT->ProcessLine(".L runLooper.C++");

runLooper("MinBias");
//runLooper("JetMETTau");
//runLooper("JetMETTauMonitor");
//runLooper("EG");

}
