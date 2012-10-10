#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TFile.h" 
#include "../../../Smurf/Core/SmurfTree.h"
#include "LeptonTreeMaker.h"
#include "SmurfDataTypes.h"
#include "processLeptonTree.h"
 

int processLeptonTree(TString outfileid, enum SmurfTree::DataType sample, TString file, bool realData, TString goodrunlist, int prescale)
{

//    gSystem->Load("libTree.so");
//    gSystem->Load("libPhysics.so");
//    gSystem->Load("libEG.so");
//    gSystem->Load("libMathCore.so");

//    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
//    gSystem->Load("libCMS2NtupleMacrosLooper.so");
//    gSystem->Load("libCMS2NtupleMacrosCORE.so");

    //
    // const config parameters
    //

    const bool lockToCoreSelectors 	= false;
    const bool useLHeleId 		= false;
    const bool useMVAeleId 		= false;
    const bool doDYNNLOw 		= false;
    //const unsigned int prescale 	= 1;
    const double integratedLumi 	= 1000.0; // pb^1  

    //
    // create looper
    //

    std::cout << "going to loop" << std::endl;
    LeptonTreeMaker *looper = new LeptonTreeMaker(lockToCoreSelectors, useLHeleId, useMVAeleId, doDYNNLOw, prescale, realData);
    looper->SetBatchMode(true);

    //
    // set up chain
    //

    /*
    int count = 0;
    TFile *f;
    while ((f = TFile::Open(file)) == 0){
      std::cout << "Couldn't open file, waiting 60s and then retrying..." << std::endl;
      if (++count > 3){
	std::cout << "Failed to open file after 3 attempts, quitting" << std::endl;
	return 999;
      }
      gSystem->Sleep(60 * 1000); // In miliseconds
    }

    TChain *chain = (TChain*) f->Get("Events");
    */

    TChain *chain = new TChain("Events");
    chain->Add(file);

    std::cout << "Entries in chain: " << chain->GetEntries() << std::endl;

    if( chain->GetEntries() == 0 ){
      std::cout << "0 entries in chain --> QUITTING!!!" << std::endl;
      return 666;
    }

    //
    // loop
    //
    looper->ScanChain(outfileid, chain, sample, integratedLumi, -1, -1, false, realData, goodrunlist);

    //
    // tidy up
    //

    delete looper;
    delete chain;

    return 0;
}
