{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L Core/StopTree.h+");
    gROOT->ProcessLine(".L libStopTreeLooper.so");

    StopTreeLooper *looper = new StopTreeLooper();

    // 
    // data
    //

    TChain *ch_data = new TChain("t");
    ch_data->Add("/tas/vimartin/SingleLepton2011/output/V00-04-11_superskim/data_smallTree.root");
    looper->setOutFileName("data_histos.root"); 
    looper->loop(ch_data, "data");

    //
    // mc
    //

    // TChain *ch_ttsl = new TChain("t");
    // ch_ttsl->Add("/tas/vimartin/SingleLepton2011/output/V00-04-11_superskim/ttsl_smallTree.root");
    // looper->setOutFileName("ttsl_histos.root");
    // looper->loop(ch_ttsl, "ttsl");

    // TChain *ch_ttdl = new TChain("t");
    // ch_ttdl->Add("/tas/vimartin/SingleLepton2011/output/V00-04-11_superskim/ttdl_smallTree.root");
    // looper->setOutFileName("ttdl_histos.root");
    // looper->loop(ch_ttdl, "ttdl");

    delete looper;
    delete ch_data;
    // delete ch_ttsl;
    // delete ch_ttdl;

}


