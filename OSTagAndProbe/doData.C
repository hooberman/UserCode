{
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

    MyScanChain *looper = new MyScanChain();

    //
    // 2011 data
    //

    TChain *chain_data2011_skim = new TChain("Events");

    chain_data2011_skim->Add("cms2/DoubleElectron_Run2011A-Apr22ReReco-v2_AOD/V04-01-05/tagAndProbeSkim/skimmed*root");
    chain_data2011_skim->Add("cms2/DoubleElectron_Run2011A-PromptReco-v1_AOD/V04-00-13/tagAndProbeSkim_merged/merged_160329_161312.root");
    chain_data2011_skim->Add("cms2/DoubleMu_Run2011A-PromptReco-v1_AOD/V04-00-13/tagAndProbeSkim_merged/merged_160329_161312.root");
    chain_data2011_skim->Add("cms2/DoubleElectron_Run2011A-PromptReco-v2_AOD/V04-01-03/tagAndProbeSkim/skimmed*root");
    chain_data2011_skim->Add("cms2/DoubleMu_Run2011A-PromptReco-v2_AOD/V04-01-03/tagAndProbeSkim/skimmed*root");
    
    //TChain *chain_data2011_express = new TChain("Events");
    //chain_data2011_express->Add("cms2/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/ntuple*.root");

    //
    // set the good run list
    //
    //looper->setGoodRunList("runlists/Cert_TopOct15_Merged_135821-147454_allPVT.txt");
    // for 2011 unofficial so far
    //looper->setGoodRunList("runlists/json_DCSONLY_ManualCert.txt");
    looper->setGoodRunList("Cert_160404-163757_7TeV_PromptReco_Collisions11_JSON_goodruns.txt");

    // do the looping
    looper->ScanChain(true, "data_skim", chain_data2011_skim);
    
    //
    // tidy up
    //

    delete looper;
    //delete chain_data2011_express;
    delete chain_data2011_skim;

}

