{

gROOT->ProcessLine(".L histtools.C++");
gROOT->ProcessLine(".L runTCMETLooperTemplate.C++");


//runTCMETLooperTemplate("data");
//runTCMETLooperTemplate("zmm");
//runTCMETLooperTemplate("zee");
//runTCMETLooperTemplate("datahighmet");
//runTCMETLooperTemplate("qcd");
runTCMETLooperTemplate("ttbar");

}
