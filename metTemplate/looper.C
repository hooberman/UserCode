#include "looper.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>


#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>

#include "CORE/CMS2.h"
#include "CORE/CMS2.cc"
#include "CORE/metSelections.cc"
#include "CORE/trackSelections.cc"
#include "CORE/eventSelections.cc"
//#include "CORE/electronSelections.cc"
#include "CORE/muonSelections.cc"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"


void looper::ScanChain (TChain* chain, const char* prefix, bool isData, bool calculateTCMET, int nEvents){

  set_goodrun_file("goodruns_official_0526.txt");

  bookHistos();

  // make a baby ntuple
  stringstream babyfilename;
  babyfilename << prefix << "_baby.root";
  MakeBabyNtuple(babyfilename.str().c_str());

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  //pass fail counters
  int nPassGoodRun    = 0;
  int nPassBPTX       = 0;
  int nPassBSC        = 0;
  int nPassBeamHalo   = 0;
  int nPassGoodTracks = 0;
  int nPassGoodVertex = 0;

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())){
    
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();

    for (unsigned int event = 0; event < nEvents; ++event){
        
      cms2.GetEntry(event);
      ++nEventsTotal;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
            
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //
      // APPLY BASIC EVENT SELECTIONS
      // TRIGGER, TRACKING AND VERTEX
      //
      if (!isData || goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))
        nPassGoodRun++;
      else continue;
          
      // determine if current event passes BPTX triggers
      if (cleaning_BPTX( isData ))
        nPassBPTX++;
      else continue;

      // determine if current event passes BSC triggers
      if (cleaning_BSC())
        nPassBSC++;
      else continue;

      // determine if current event passes beam halo triggers
      if (cleaning_beamHalo())
        nPassBeamHalo++;
      else continue;

      // determine if current event is a beam scraping event
      if (cleaning_goodTracks())
        nPassGoodTracks++;
      else continue;

      // determine if current event has a good vertex
      if (cleaning_goodVertex())
        nPassGoodVertex++;
      else continue;
          
      //
      // N.B. BABY NTUPLE IS FILLED
      // FOR EACH EVENT
      //
      InitBabyNtuple();

      // event stuff
      run_    = cms2.evt_run();
      lumi_   = cms2.evt_lumiBlock();
      event_  = cms2.evt_event();

      //met quantities--------------------------------------------------------------------------
      //this is probably more info than necessary, can probably reduce number of met branches

      //calomet
      met_       = cms2.evt_met();
      metphi_    = cms2.evt_metPhi();
      sumet_     = cms2.evt_sumet();

      // pf met stuff
      pfmet_    = cms2.evt_pfmet();
      pfmetphi_ = cms2.evt_pfmetPhi();
      pfsumet_  = cms2.evt_pfsumet();

      //muon-corrected met stuff
      mumet_    = cms2.evt_metMuonCorr();
      mumetphi_ = cms2.evt_metMuonCorrPhi();
      musumet_  = cms2.evt_sumetMuonCorr();

      //muon-corrected JES met stuff
      mujesmet_    = cms2.evt_metMuonJESCorr();
      mujesmetphi_ = cms2.evt_metMuonJESCorrPhi();
      mujessumet_  = -9999.; //cms2.evt_sumetMuonJESCorr(); //branch doesn't exist!!!

      // genmet stuff
      if (!isData){
        genmet_     = cms2.gen_met();
        genmetphi_  = cms2.gen_metPhi();
        gensumet_   = cms2.gen_sumEt();
      }
      
      if(calculateTCMET){
        
        // calculate tcmet on-the-fly
        bool usePV            = false;
        bool useHFcleaning    = false;
        bool useHCALcleaning  = false;
        bool useECALcleaning  = false;
        
        if ( isData ){
          useHFcleaning   = true;
          useHCALcleaning = true;
          useECALcleaning = true;
        }
        
        metStruct structMET = correctedTCMET(usePV, useHFcleaning, useHCALcleaning, useECALcleaning);
        
        tcmet_     = structMET.met;
        tcmetphi_  = structMET.metphi;
        tcsumet_   = structMET.sumet;
      }else{
        
        //just take tcmet from event
        tcmet_     = evt_tcmet();
        tcmetphi_  = evt_tcmetPhi();
        tcsumet_   = evt_tcsumet();

      }


      //photon quantities-----------------------------------------------------------------------
      nPhotons_    = 0;
      maxPhotonPt_ = 0.;

      //count photons pt > 10 GeV
      for (unsigned int iphoton = 0 ; iphoton < photons_p4().size() ; iphoton++) {

        LorentzVector vphoton = photons_p4().at(iphoton);

        if(vphoton.pt() > 10){
          nPhotons_++;
          if(vphoton.pt() > maxPhotonPt_) maxPhotonPt_ = vphoton.pt();
        }

      }

      //jet quantities--------------------------------------------------------------------------
      
      nJets_      = 0;
      sumJetPt_   = 0.;

      //count pfjets with pt > 20 and eta < 3
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
    
        LorentzVector vjet = pfjets_p4().at(ijet);
  
        //check if jet has overlapping photon with pt > 10 GeV
        bool overlapPhoton = false;

        for (unsigned int iphoton = 0 ; iphoton < photons_p4().size() ; iphoton++ ){
          LorentzVector vphoton  = photons_p4().at(iphoton);
          if (dRbetweenVectors(vjet, vphoton) < 0.5 && vphoton.pt() > 10) overlapPhoton = true;
        }

        if(overlapPhoton) continue;

        if ( vjet.pt() > 20. && fabs( vjet.eta() ) < 3.){
          nJets_++;
          sumJetPt_ += vjet.pt();
        }

        
      }

      //----------------------------------------------------------------------------------------
      
      FillBabyNtuple();
      
    } // end loop over events
  } // end loop over files
  
  cout << "\n\n********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events that pass good run/lumi: " << nPassGoodRun 
       << " (" << 100*(double)nPassGoodRun/nEventsTotal << "% of total)" << endl;
  cout << "Total number of events that pass BPTX trigger: " << nPassBPTX
       << " (" << 100*(double)nPassBPTX/nPassGoodRun << "%)" << endl;
  cout << "Total number of events that pass BSC trigger: " << nPassBSC
       << " (" << 100*(double)nPassBSC/nPassBPTX << "%)" << endl;
  cout << "Total number of events that pass BeamHalo trigger: " << nPassBeamHalo
       << " (" << 100*(double)nPassBeamHalo/nPassBSC << "%)" << endl;
  cout << "Total number of events that pass tracking cuts: " << nPassGoodTracks
       << " (" << 100*(double)nPassGoodTracks/nPassBeamHalo << "%)" << endl;
  cout << "Total number of events that pass vertex cuts: " << nPassGoodVertex
       << " (" << 100*(double)nPassGoodVertex/nPassGoodTracks << "%)" << endl;
  cout << endl << endl;
  
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  CloseBabyNtuple();

  // make histos rootfile
  stringstream rootfilename;
  rootfilename << prefix << "_histos.root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(rootfilename.str().c_str());
  deleteHistos();
  
} // end ScanChain

void looper::InitBabyNtuple ()
{
  // event stuff
  run_    = -999999;
  lumi_   = -999999;
  event_  = -999999;
 
  // genmet stuff
  genmet_     = -999999.;
  genmetphi_  = -999999.;
  gensumet_   = -999999.;

  // pfmet stuff
  pfmet_     = -999999.;
  pfmetphi_  = -999999.;
  pfsumet_   = -999999.;

  // calomet stuff
  met_          = -999999.;
  metphi_       = -999999.;
  sumet_        = -999999.;

  // muon-corrected calomet stuff
  mumet_        = -999999.;
  mumetphi_     = -999999.;
  musumet_      = -999999.;

  // calomet stuff
  mujesmet_     = -999999.;
  mujesmetphi_  = -999999.;
  mujessumet_   = -999999.;

  // tcmet stuff
  tcmet_        = -999999.;
  tcmetphi_     = -999999.;
  tcsumet_      = -999999.;

}

void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  //book histos here
 
}


void looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("tree", "A Baby Ntuple");

  babyTree_->Branch("run",          &run_,          "run/I"  );
  babyTree_->Branch("lumi",         &lumi_,         "lumi/I" );
  babyTree_->Branch("event",        &event_,        "event/I");
  babyTree_->Branch("pfmet",        &pfmet_,        "pfmet/F"   );
  babyTree_->Branch("pfmetphi",     &pfmetphi_,     "pfmetphi/F");
  babyTree_->Branch("pfsumet",      &pfsumet_,      "pfsumet/F" );
  babyTree_->Branch("met",          &met_,          "met/F"      );
  babyTree_->Branch("metphi",       &metphi_,       "metphi/F"   );
  babyTree_->Branch("sumet",        &sumet_,        "sumet/F"    );
  babyTree_->Branch("mumet",        &mumet_,        "mumet/F"      );
  babyTree_->Branch("mumetphi",     &mumetphi_,     "mumetphi/F"   );
  babyTree_->Branch("musumet",      &musumet_,      "musumet/F"    );
  babyTree_->Branch("mujesmet",     &mujesmet_,     "mujesmet/F"      );
  babyTree_->Branch("mujesmetphi",  &mujesmetphi_,  "mujesmetphi/F"   );
  babyTree_->Branch("mujessumet",   &mujessumet_,   "mujessumet/F"    );
  babyTree_->Branch("genmet",       &genmet_,       "genmet/F"   );
  babyTree_->Branch("genmetphi",    &genmetphi_,    "genmetphi/F");
  babyTree_->Branch("gensumet",     &gensumet_,     "gensumet/F" );
  babyTree_->Branch("tcmet",        &tcmet_,        "tcmet/F"      );
  babyTree_->Branch("tcmetphi",     &tcmetphi_,     "tcmetphi/F"   );
  babyTree_->Branch("tcsumet",      &tcsumet_,      "tcsumet/F"    );
  babyTree_->Branch("nphotons",     &nPhotons_,     "nphotons/I"    );
  babyTree_->Branch("maxphotonpt",  &maxPhotonPt_,  "maxphotonpt/I"    );
  babyTree_->Branch("njets",        &nJets_,        "njets/I"    );
  babyTree_->Branch("sumjetpt",     &sumJetPt_,     "sumjetpt/I"    );
  
}



void looper::FillBabyNtuple ()
{
  babyTree_->Fill();
}

void looper::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}




