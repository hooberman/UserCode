{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

    MyScanChain *looper = new MyScanChain();

    TString ntuple_location = "/tas/cms2/";
    //TString ntuple_location = "/nfs-3/userdata/cms2/";

    // DYEE
    TChain *chain_dyee   = new TChain("Events");
    chain_dyee   ->Add(ntuple_location + "DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root");

    // DYMM
    TChain *chain_dymumu = new TChain("Events");
    chain_dymumu->Add(ntuple_location + "DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root");
    
    cout << "Doing DY->mm" << endl;
    looper->ScanChain(false, "dymm", chain_dymumu, 1.0);
    cout << "Doing DY->ee" << endl;
    looper->ScanChain(false, "dyee", chain_dyee, 1.0);
   
    // tidy up
    delete looper;
    delete chain_dymumu;
    delete chain_dyee;

}

