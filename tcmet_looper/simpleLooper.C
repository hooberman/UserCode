#include "simpleLooper.h"
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
#include "CORE/trackSelections.h"
#include "CORE/metSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

char* iter = "default";

bool usePV    = false;
bool makebaby = false;
bool debug    = false;

// inline double fround(double n, double d){
//   return floor(n * pow(10., d) + .5) / pow(10., d);
// }

void simpleLooper::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  bookHistos();

  //set_goodrun_file("Cert_132440-136119_7TeV_May27thReReco_Collisions10_JSON.txt");
  ofile.open( Form( "output/%s_%s_events.txt" , prefix , iter) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  MakeBabyNtuple( Form( "output/%s_%s_baby.root" , prefix , iter) );

  //pass fail counters
  int nPassGoodRun    = 0;
  int nPassBPTX       = 0;
  int nPassBSC        = 0;
  int nPassBeamHalo   = 0;
  int nPassGoodTracks = 0;
  int nPassGoodVertex = 0;

  if( debug ) cout << "Begin looping over files" << endl;

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next()))
    {
      TFile f(currentFile->GetTitle());
      TTree *tree = (TTree*)f.Get("Events");
      cms2.Init(tree);

      // event loop
      unsigned int nEvents = tree->GetEntries();

      for (unsigned int event = 0; event < nEvents; ++event)
        {
          if( debug ) cout << "Event " << event << endl;

          cms2.GetEntry(event);
          ++nEventsTotal;

          // progress feedback to user
          if (nEventsTotal % 1000 == 0)
            {
              // xterm magic from L. Vacavant and A. Cerri
              if (isatty(1))
                {
                  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                         "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
                  fflush(stdout);
                }
            }
	  
	  cout << __LINE__ << endl;
          //
          // APPLY BASIC EVENT SELECTIONS
          // TRIGGER, TRACKING AND VERTEX
          //
          
	  //if (!isData || goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))
	  //  nPassGoodRun++;
	  //else continue;
          
	  // determine if current event passes BPTX triggers
	  // if (cleaning_BPTX( isData ))
// 	    nPassBPTX++;
// 	  else continue;
          //	  cout << __LINE__ << endl;
	  // determine if current event passes BSC triggers
	  //if (cleaning_BSC())
	  //  nPassBSC++;
	  //else continue;
          	  cout << __LINE__ << endl;
	  // determine if current event passes beam halo triggers
	  // if (cleaning_beamHalo())
// 	    nPassBeamHalo++;
// 	  else continue;
          	  cout << __LINE__ << endl;
	  // determine if current event is a beam scraping event
	  // if (cleaning_goodTracks())
// 	    nPassGoodTracks++;
// 	  else continue;
          	  cout << __LINE__ << endl;
		  cout << "rho " << cms2.vtxs_position()[0].Rho() << endl;
  // determine if current event has a good vertex
	   if (cleaning_goodVertex())
 	    nPassGoodVertex++;
 	  else continue;

          	  cout << __LINE__ << endl;
	  
	  if( debug ) cout << "Pass event selection" << endl;
          
          
          //
          // N.B. BABY NTUPLE IS FILLED
          // FOR EACH EVENT
          //
          InitBabyNtuple();

          // event stuff
          run_    = cms2.evt_run();
          lumi_   = cms2.evt_lumiBlock();
          event_  = cms2.evt_event();

          // pf met stuff
          pfmet_    = cms2.evt_pfmet();
          pfmetphi_ = cms2.evt_pfmetPhi();
          pfsumet_  = cms2.evt_pfsumet();

          // raw  tcmet stuff
          rawtcmet_    = cms2.evt_tcmet();
          rawtcmetphi_ = cms2.evt_tcmetPhi();
          rawtcsumet_  = cms2.evt_tcsumet();
          
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
            //genmet_     = cms2.gen_met();
            //genmetphi_  = cms2.gen_metPhi();
            //gensumet_   = cms2.gen_sumEt();
            genmet_     = 0;
            genmetphi_  = 0;
            gensumet_   = 0;
          }

          //calomet
          met_       = cms2.evt_met();
          metphi_    = cms2.evt_metPhi();
          sumet_     = cms2.evt_sumet();

          fillUnderOverFlow( hmumet , mumet_ );
          fillUnderOverFlow( hmujesmet , mujesmet_ );
          fillUnderOverFlow( htcmet , tcmet_ );
          fillUnderOverFlow( hrawtcmet , rawtcmet_ );
  
          fillUnderOverFlow( hdmumet , mumet_ - genmet_ );
          fillUnderOverFlow( hdmujesmet , mujesmet_ - genmet_ );
          fillUnderOverFlow( hdtcmet , tcmet_ - genmet_ );
          fillUnderOverFlow( hdrawtcmet , rawtcmet_ - genmet_ );

	  fillUnderOverFlow( hdtcmet_mumet    , tcmet_ - mumet_ );
 	  fillUnderOverFlow( hdtcmet_mujesmet , tcmet_ - mujesmet_ );

 	  fillUnderOverFlow( htcmetNew ,    cms2.evtNew_tcmet() );
 	  fillUnderOverFlow( hpfmet ,       cms2.evt_pfmet() );
 	  fillUnderOverFlow( htcmetNewPFC , cms2.evtNewPFC_tcmet() );
   
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
  //stringstream histfile;
  //histfile << prefix << "_histos.root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  //saveHist( histfile.str().c_str() );
  saveHist( Form( "output/%s_%s_histos.root" , prefix , iter ) );
  deleteHistos();
  
} // end ScanChain

void simpleLooper::InitBabyNtuple ()
{
  // event stuff
  run_    = -999999;
  lumi_   = -999999;
  event_  = -999999;
  hase_   = false;

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

  // calomet stuff
  rawtcmet_     = -999999.;
  rawtcmetphi_  = -999999.;
  rawtcsumet_   = -999999.;

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
  ectcmet_      = -999999.;
  ectcmetphi_   = -999999.;
  ectcsumet_    = -999999.;

  // electron stuff
  eldr_   = -999999.;
  elpt_   = -999999.;
  eleta_  = -999999.;


}

void simpleLooper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hmumet    = new TH1F("hmumet",   "",100,0,100);
  hmujesmet = new TH1F("hmujesmet","",100,0,100);
  htcmet    = new TH1F("htcmet",   "",100,0,100);
  hrawtcmet = new TH1F("hrawtcmet","",100,0,100);

  hdmumet    = new TH1F("hdmumet",   "",400,-200,200);
  hdmujesmet = new TH1F("hdmujesmet","",400,-200,200);
  hdtcmet    = new TH1F("hdtcmet",   "",400,-200,200);
  hdrawtcmet = new TH1F("hdrawtcmet","",400,-200,200);

  hdtcmet_mumet    = new TH1F("hdtcmet_mumet","",   400,-200,200);
  hdtcmet_mujesmet = new TH1F("hdtcmet_mujesmet","",400,-200,200);

  htcmetNew    = new TH1F("htcmetNew",     "calo-tcmet",100,0,100);
  hpfmet       = new TH1F("hpfmet",        "pfmet",100,0,100);
  htcmetNewPFC = new TH1F("htcmetNewPFC",  "pfc-tcmet",100,0,100);
}


void simpleLooper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  eventTree_ = new TTree("Events", "Events Tree");

  eventTree_->Branch("run"              , &run_              , "run/I"  );
  eventTree_->Branch("lumi"             , &lumi_             , "lumi/I" );
  eventTree_->Branch("event"            , &event_            , "event/I");
  eventTree_->Branch("hase"             , &hase_             , "hase/O" );

  eventTree_->Branch("pfmet"            , &pfmet_            , "pfmet/F"   );
  eventTree_->Branch("pfmetphi"         , &pfmetphi_         , "pfmetphi/F");
  eventTree_->Branch("pfsumet"          , &pfsumet_          , "pfsumet/F" );

  eventTree_->Branch("met"              , &met_              , "met/F"      );
  eventTree_->Branch("metphi"           , &metphi_           , "metphi/F"   );
  eventTree_->Branch("sumet"            , &sumet_            , "sumet/F"    );

  eventTree_->Branch("mumet"            , &mumet_            , "mumet/F"      );
  eventTree_->Branch("mumetphi"         , &mumetphi_         , "mumetphi/F"   );
  eventTree_->Branch("musumet"          , &musumet_          , "musumet/F"    );

  eventTree_->Branch("mujesmet"         , &mujesmet_         , "mujesmet/F"      );
  eventTree_->Branch("mujesmetphi"      , &mujesmetphi_      , "mujesmetphi/F"   );
  eventTree_->Branch("mujessumet"       , &mujessumet_       , "mujessumet/F"    );

  eventTree_->Branch("genmet"           , &genmet_           , "genmet/F"   );
  eventTree_->Branch("genmetphi"        , &genmetphi_        , "genmetphi/F");
  eventTree_->Branch("gensumet"         , &gensumet_         , "gensumet/F" );
 
  eventTree_->Branch("tcmet"            , &tcmet_            , "tcmet/F"      );
  eventTree_->Branch("tcmetphi"         , &tcmetphi_         , "tcmetphi/F"   );
  eventTree_->Branch("tcsumet"          , &tcsumet_          , "tcsumet/F"    );

  eventTree_->Branch("rawtcmet"         , &rawtcmet_         , "rawtcmet/F"      );
  eventTree_->Branch("rawtcmetphi"      , &rawtcmetphi_      , "rawtcmetphi/F"   );
  eventTree_->Branch("rawtcsumet"       , &rawtcsumet_       , "rawtcsumet/F"    );

  eventTree_->Branch("raw35Xtcmet"      , &raw35Xtcmet_      , "raw35Xtcmet/F"      );
  eventTree_->Branch("raw35Xtcmetphi"   , &raw35Xtcmetphi_   , "raw35Xtcmetphi/F"   );
  eventTree_->Branch("raw35Xtcsumet"    , &raw35Xtcsumet_    , "raw35Xtcsumet/F"    );

  //rootdir->cd();
  babyFile_->cd();
  trackTree_ = new TTree("Tracks", "Tracks Tree");

  trackTree_->Branch("pt"           , &trk_pt_           , "pt/F"           );
  trackTree_->Branch("d0vtx"        , &trk_d0vtx_        , "d0vtx/F"        );
  trackTree_->Branch("d0corr"       , &trk_d0corr_       , "d0corr/F"       );
  trackTree_->Branch("nhits"        , &trk_nhits_        , "nhits/I"        );
  trackTree_->Branch("chi2"         , &trk_chi2_         , "chi2/F"         );
  trackTree_->Branch("ndf"          , &trk_ndf_          , "ndf/I"          );
  trackTree_->Branch("pterr"        , &trk_pterr_        , "pterr/F"        );
  trackTree_->Branch("phi"          , &trk_phi_          , "phi/F"          );
  trackTree_->Branch("eta"          , &trk_eta_          , "eta/F"          );
  trackTree_->Branch("qual"         , &trk_qual_         , "qual/I"         );
  trackTree_->Branch("pass"         , &trk_pass_         , "pass/I"         );
  trackTree_->Branch("algo"         , &trk_algo_         , "algo/I"         );
  trackTree_->Branch("chg"          , &trk_chg_          , "chg/I"          );
  trackTree_->Branch("nlayers"      , &trk_nlayers_      , "nlayers/I"      );
  trackTree_->Branch("isdup"        , &trk_isdup_        , "isdup/I"        );
  trackTree_->Branch("passnlayers"  , &trk_passnlayers_  , "passnlayers/I"  );

  trackTree_->Branch("vtxz"      , &trk_vtxz_      , "vtxz/F"       );
  trackTree_->Branch("z0"        , &trk_z0_        , "z0/F"       );
  trackTree_->Branch("z0corr"    , &trk_z0corr_    , "z0corr/F"       );

  trackTree_->Branch("nlayersp"  , &trkp_nlayers_  , "nlayersp/I"      );
  trackTree_->Branch("nhitsp"    , &trkp_nhits_    , "nhitsp/I"      );
  trackTree_->Branch("etap"      , &trkp_eta_      , "etap/F"      );
  trackTree_->Branch("chgp"      , &trkp_chg_      , "chgp/I"      );
  trackTree_->Branch("ptp"       , &trkp_pt_       , "ptp/F"      );
  trackTree_->Branch("algop"     , &trkp_algo_     , "algop/I"      );
  trackTree_->Branch("dr"        , &trkp_dr_       , "dr/F"      );
  trackTree_->Branch("dcotth"    , &trkp_dcotth_   , "dcotth/F"      );
  trackTree_->Branch("dphi"      , &trkp_dphi_     , "dphi/F"      );

}

//--------------------------------------------------------------------

void simpleLooper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

//--------------------------------------------------------------------

void simpleLooper::CloseBabyNtuple (){

  babyFile_->cd();
  eventTree_->Write();
  trackTree_->Write();
  babyFile_->Close();
  
}
