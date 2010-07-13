#include "tcmetLooperTemplate.h"
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

// #include "CORE/CMS2.h"
// #include "CORE/trackSelections.h"
// #include "CORE/metSelections.h"
// #include "CORE/eventSelections.h"
// #include "CORE/electronSelectionsParameters.h"
// #include "CORE/electronSelections.h"
// #include "CORE/muonSelections.h"
// #include "Tools/goodrun.cc"
// #include "CORE/utilities.cc"
// #include "histtools.h"

#include "CORE/CMS2.cc"
#include "CORE/trackSelections.cc"
#include "CORE/metSelections.cc"
#include "CORE/eventSelections.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/electronSelections.cc"
#include "CORE/muonSelections.cc"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.C"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

//char* iter = "default";
//char* iter = "oldTQ";
char* iter = "dupveto";
//char* iter = "nlayerscut_5_7";
//char* iter = "trkpt_lt10eta25";
//char* iter = "eta25";

bool usePV    = false;
bool makebaby = false;
bool debug    = false;

// inline double fround(double n, double d){
//   return floor(n * pow(10., d) + .5) / pow(10., d);
// }

void tcmetLooperTemplate::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  bookHistos();

  set_goodrun_file("jsonlist_132440_132697.txt");
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

          //
          // APPLY BASIC EVENT SELECTIONS
          // TRIGGER, TRACKING AND VERTEX
          //
          if( strcmp( prefix , "dipion" ) !=0 ){
            
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
          }
          
          if( strcmp( prefix , "zee")   == 0 || strcmp( prefix , "zmm")   == 0 || strcmp( prefix , "ttbar") == 0 ){
            if (!isGoodDilepton())    continue;
          } 

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

          // raw  tcmet35X stuff
          raw35Xtcmet_    = -9999;
          raw35Xtcmetphi_ = -9999;
          raw35Xtcsumet_  = -9999;
          //raw35Xtcmet_    = cms2.evt35X_tcmet();
          //raw35Xtcmetphi_ = cms2.evt35X_tcmetPhi();
          //raw35Xtcsumet_  = cms2.evt35X_tcsumet();
          
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

          //calomet
          met_       = cms2.evt_met();
          metphi_    = cms2.evt_metPhi();
          sumet_     = cms2.evt_sumet();
          
          // calculate tcmet on-the-fly
          bool electronVetoCone = true;
          bool usePV            = false;
          bool useHFcleaning    = false;
          bool useHCALcleaning  = false;
          bool useECALcleaning  = false;
          
          if ( isData ){
            useHFcleaning   = true;
            useHCALcleaning = true;
            useECALcleaning = true;
          }

          metStruct structMET = correctedTCMET( usePV, electronVetoCone, useHFcleaning, useHCALcleaning, useECALcleaning );
            
          //tcmet
          tcmet_     = structMET.met;
          tcmetphi_  = structMET.metphi;
          tcsumet_   = structMET.sumet;

          if( makebaby || tcmet_ > 45 ){
            
//             cout << "|" << setw(12)  << cms2.evt_run()         << setw(4)
//                  << "|" << setw(12)  << cms2.evt_lumiBlock()   << setw(4)
//                  << "|" << setw(12)  << cms2.evt_event()       << setw(4)
//                  << "|" << setw(12)  << fround( tcmet_ , 1 )   << setw(4) << "|" << endl;
            
            structMET = correctedTCMET( usePV, electronVetoCone, useHFcleaning, useHCALcleaning, useECALcleaning , true, ofile);          
            
            tcmet_     = structMET.met;
            tcmetphi_  = structMET.metphi;
            tcsumet_   = structMET.sumet;
            
	    eventTree_->Fill();
          }

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

          
          bool badvtx = false;
          
          if (cms2.vtxs_isFake()[0])                     badvtx = true;
          if (cms2.vtxs_ndof()[0] < 4.)                  badvtx = true;
          if (cms2.vtxs_position()[0].Rho() > 2.0)       badvtx = true;
          if (fabs(cms2.vtxs_position()[0].Z()) > 15.0)  badvtx = true;
      
          if( badvtx )  trk_vtxz_ = -9999;
          else          trk_vtxz_ = cms2.vtxs_position()[0].Z();
        

	  if( makebaby || tcmet_ > 45 ){

            duplicateTracks.clear();
            findDuplicateTracks();
      
	    for( unsigned int i = 0; i < cms2.trks_trk_p4().size(); i++ ) {

	      if( isMuon( i ) )      continue;
	      if( isElectron( i ) )  continue;
	  	      
	      //if (useElectronVetoCone && closeToElectron(i))  continue;
	      
	      trk_pass_     = isGoodTrack( i, usePV ) ? 1 : 0;            
	      trk_pt_       = cms2.trks_trk_p4().at(i).pt();
	      trk_d0vtx_    = cms2.trks_d0vtx().at(i);
	      trk_d0corr_   = cms2.trks_d0corr().at(i);
	      trk_nhits_    = cms2.trks_validHits().at(i);
	      trk_chi2_     = cms2.trks_chi2().at(i);
	      trk_ndf_      = (int) cms2.trks_ndof().at(i);
	      trk_pterr_    = cms2.trks_ptErr().at(i);
	      trk_phi_      = cms2.trks_trk_p4().at(i).phi();
	      trk_eta_      = cms2.trks_trk_p4().at(i).eta();
	      trk_qual_     = isTrackQuality( i , (1 << highPurity) ) ? 1 : 0;
              trk_algo_     = cms2.trks_algo().at(i);
              trk_chg_      = cms2.trks_charge().at(i);
              trk_nlayers_  = cms2.trks_nlayers().at(i);
              trk_z0_       = cms2.trks_z0().at(i);
              trk_z0corr_   = cms2.trks_z0corr().at(i);
              trk_isdup_    = 0;

              for( unsigned int idup = 0 ; idup < duplicateTracks.size() ; ++idup ){
                if( i == duplicateTracks.at(idup) ) trk_isdup_ = 1;
              }

              trk_passnlayers_ = 1;
              if( trk_algo_ < 8 && trks_nlayers().at(i) < 5 ) trk_passnlayers_ = 0;
              if( trk_algo_ > 7 && trks_nlayers().at(i) < 7 ) trk_passnlayers_ = 0;

              //find closest track
              float drmin = 100;
              int   jmin  = -1;
              for( unsigned int j = 0; j < cms2.trks_trk_p4().size(); j++ ) {

                if( isMuon( i ) )      continue;
                if( isElectron( i ) )  continue;
                if( i == j )           continue;

                float dr = dRbetweenVectors(cms2.trks_trk_p4().at(i),cms2.trks_trk_p4().at(j));
                if( dr < drmin ){
                  drmin = dr;
                  jmin  = j;
                }
              }
              
              if( jmin > -1 ){
                trkp_nlayers_  =  cms2.trks_nlayers().at(jmin);
                trkp_nhits_    =  cms2.trks_validHits().at(jmin);
                trkp_eta_      =  cms2.trks_trk_p4().at(jmin).eta();
                trkp_chg_      =  cms2.trks_charge().at(jmin);
                trkp_dr_       =  drmin;
                trkp_pt_       =  cms2.trks_trk_p4().at(jmin).pt();
                trkp_algo_     =  cms2.trks_algo().at(jmin);
                trkp_dcotth_   =  1./tan(cms2.trks_trk_p4().at(i).theta()) - 1./tan(cms2.trks_trk_p4().at(jmin).theta());
                trkp_dphi_     =  fabs( cms2.trks_trk_p4().at(i).phi() - cms2.trks_trk_p4().at(jmin).phi() );
                if( trkp_dphi_ > TMath::Pi() ) trkp_dphi_ = TMath::TwoPi() - trkp_dphi_;
              }
              else{
                trkp_nlayers_  = -9999;
                trkp_nhits_    = -9999;
                trkp_eta_      = -9999.;
                trkp_chg_      = -9999;
                trkp_dr_       = -9999.;
                trkp_pt_       = -9999.;
                trkp_algo_     = -9999;
                trkp_dcotth_   = -9999.;
                trkp_dphi_     = -9999.;
              }

	      trackTree_->Fill();
            }
          }
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

void tcmetLooperTemplate::InitBabyNtuple ()
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

void tcmetLooperTemplate::bookHistos(){

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

  
}


void tcmetLooperTemplate::MakeBabyNtuple (const char* babyFileName)
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

void tcmetLooperTemplate::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

//--------------------------------------------------------------------

void tcmetLooperTemplate::CloseBabyNtuple (){

  babyFile_->cd();
  eventTree_->Write();
  trackTree_->Write();
  babyFile_->Close();
  
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isTruthZee () {
  
  if ( cms2.hyp_p4().size() != 1)              return false;
  if ( cms2.els_p4().size() != 2)              return false;
  if ( abs( cms2.els_mc3_id().at(0) ) != 11 )  return false; 
  if ( abs( cms2.els_mc3_id().at(1) ) != 11 )  return false; 

  return true;

}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isTruthZmm () {
  
  if ( cms2.hyp_p4().size() != 1)              return false;
  if ( cms2.mus_p4().size() != 2)              return false;
  if ( abs( cms2.mus_mc3_id().at(0) ) != 13 )  return false; 
  if ( abs( cms2.mus_mc3_id().at(1) ) != 13 )  return false; 

  return true;

}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::jetVeto () {

  for (unsigned int ijet = 0; ijet < cms2.jets_p4().size(); ijet++) {
    
    LorentzVector vjet = cms2.jets_p4().at(ijet);
    LorentzVector vlt  = cms2.hyp_lt_p4()[0];
    LorentzVector vll  = cms2.hyp_ll_p4()[0];

    if (dRbetweenVectors(vjet, vll) < 0.5) continue;
    if (dRbetweenVectors(vjet, vlt) < 0.5) continue;
    
    LorentzVector bah = cms2.jets_p4().at(ijet);
    
    if ( bah.pt() * cms2.jets_cor().at(ijet) > 20. && fabs( bah.eta() ) < 3.)
      return true;
  
  }
  
  return false;
  
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isGoodDilepton ()
{
  //cout <<"tcmetLooperTemplate::isGoodDilepton"<<endl;
  bool isgood = false;

  for (unsigned int i = 0; i < cms2.hyp_p4().size(); ++i){
    
    if (cms2.hyp_lt_charge()[i] * cms2.hyp_ll_charge()[i] > 0)           continue;
    if (min(cms2.hyp_lt_p4()[i].pt(), cms2.hyp_ll_p4()[i].pt()) < 20.)   continue;
    
    //muon ID
    if (abs(cms2.hyp_ll_id()[i]) == 13  && (! (fabs(cms2.hyp_ll_p4()[i].eta()) < 2.4 && muonId(cms2.hyp_ll_index()[i]))))   continue;
    if (abs(cms2.hyp_lt_id()[i]) == 13  && (! (fabs(cms2.hyp_lt_p4()[i].eta()) < 2.4 && muonId(cms2.hyp_lt_index()[i]))))   continue;
    
    //cand01
    //cout << "ELECTRON SELECTION TURNED OFF!!!" << endl;
    //exit(0);

    if (abs(cms2.hyp_ll_id()[i]) == 11  && (! pass_electronSelection( cms2.hyp_ll_index()[i] , electronSelection_cand01 ))) continue;
    if (abs(cms2.hyp_lt_id()[i]) == 11  && (! pass_electronSelection( cms2.hyp_lt_index()[i] , electronSelection_cand01 ))) continue;
    
    isgood = true;
  }
  
  return isgood;
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isGoodZee ()
{
  bool isgood = false;

  for (unsigned int i = 0; i < cms2.hyp_p4().size(); ++i)
    {
      if (cms2.hyp_type()[i] != 3)
        continue;

      if (cms2.hyp_lt_charge()[i] * cms2.hyp_ll_charge()[i] > 0)
        continue;

      if (min(cms2.hyp_lt_p4()[i].pt(), cms2.hyp_ll_p4()[i].pt()) < 20.)
        continue;

      if (! pass_electronSelection( cms2.hyp_ll_index()[i] , electronSelection_cand01 ))
        continue;

      if (! pass_electronSelection( cms2.hyp_lt_index()[i] , electronSelection_cand01 ))
        continue;
      
      //cout << "ELECTRON SELECTION TURNED OFF!!!" << endl;
      //exit(0);

      float dphi = fabs( cms2.hyp_lt_p4()[i].phi() - cms2.hyp_ll_p4()[i].phi() );
      if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;

      if( fabs( dphi - TMath::Pi() ) > 0.3 ) continue;

      isgood = true;
    }

  return isgood;
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isGoodZmm ()
{
  bool isgood = false;

  for (unsigned int i = 0; i < cms2.hyp_p4().size(); ++i)
    {
      if (cms2.hyp_type()[i] != 0)
        continue;

      if (cms2.hyp_lt_charge()[i] * cms2.hyp_ll_charge()[i] > 0)
        continue;

      if (min(cms2.hyp_lt_p4()[i].pt(), cms2.hyp_ll_p4()[i].pt()) < 20.)
        continue;

      if (!muonId(cms2.hyp_lt_index()[i]))
        continue;

      if (!muonId(cms2.hyp_ll_index()[i]))
        continue;

      isgood = true;
    }

  return isgood;
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isGoodTrack( int index, bool usePV ) {
  
     float corrected_d0 = cms2.trks_d0corr().at(index);

     if(usePV)
	  corrected_d0 = cms2.trks_d0vtx().at(index);

     if( cms2.trks_algo().at(index) < 8 ) {

	  float d0cut = min( sqrt( pow(0.015,2) + pow(0.5/cms2.trks_trk_p4().at(index).pt(),2) ), 0.3);
          if( d0cut > 0.3 ) d0cut = 0.3;
	  if( fabs( corrected_d0 ) > d0cut ) return false;
     }
     else {
	  if( cms2.trks_validHits().at(index) < 11 )                              return false;
	  if( cms2.trks_chi2().at(index) / cms2.trks_ndof().at(index) > 3 )            return false;
	  if( cms2.trks_ptErr().at(index) / cms2.trks_trk_p4().at(index).pt() > 0.10 ) return false;
     }

     if( cms2.trks_trk_p4().at(index).pt() > 100 )                           return false;
     if( fabs( cms2.trks_trk_p4().at(index).eta() ) > 2.65 )                 return false;
     if( cms2.trks_validHits().at(index) < 6 )                               return false;
     if( cms2.trks_chi2().at(index) / cms2.trks_ndof().at(index) > 5 )            return false;
     if( cms2.trks_ptErr().at(index) / cms2.trks_trk_p4().at(index).pt() > 0.20 ) return false;
     if( !isTrackQuality( index, (1 << highPurity) ) )                  return false;

     return true;
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isMuon( int index ) {

  for( unsigned int i = 0; i < cms2.mus_p4().size(); i++ ) {

    if( cms2.mus_trkidx().at(i) == index ) return true;
  }

  return false;
}

//--------------------------------------------------------------------

bool tcmetLooperTemplate::isElectron( int index ) {

  for( unsigned int i = 0; i < cms2.els_p4().size(); i++ ) {

    if( cms2.els_trkidx().at(i) == index && cms2.els_hOverE().at(i) < 0.1 ) return true;
  }

  return false;
}

//--------------------------------------------------------------------

void tcmetLooperTemplate::findDuplicateTracks(){

  for( unsigned int i = 0; i < cms2.trks_trk_p4().size(); i++ ) {

    //if( trks_trk_p4().at(i).pt() < 5. ) continue;

    for( unsigned int j = i + 1 ; j < cms2.trks_trk_p4().size(); j++ ) {

      if( cms2.trks_charge().at(i) * cms2.trks_charge().at(j) < 0 ) continue;
      //if( trks_trk_p4().at(j).pt() < 5. )                 continue;
      
      float dphi = fabs( cms2.trks_trk_p4().at(i).phi() - cms2.trks_trk_p4().at(j).phi() );
      if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;

      float dcotth = fabs( 1./tan( cms2.trks_trk_p4().at(i).theta() ) - 1./tan( cms2.trks_trk_p4().at(j).theta() ) );

      if( dphi   > 0.03 )      continue;
      if( dcotth > 0.0006 )    continue;

      int iVeto = vetoTrack( i , j );

      duplicateTracks.push_back(iVeto);
      
    }
  }
}

//--------------------------------------------------------------------

int tcmetLooperTemplate::vetoTrack( int i , int j ){

  //given 2 tracks, decide which one to veto

  if     ( trks_validHits().at(i) < trks_validHits().at(j) )                               return i;
  else if( trks_validHits().at(i) > trks_validHits().at(j) )                               return j;
  else if( trks_chi2().at(i) / trks_ndof().at(i) > trks_chi2().at(j) / trks_ndof().at(j) ) return i;  
  else if( trks_chi2().at(i) / trks_ndof().at(i) < trks_chi2().at(j) / trks_ndof().at(j) ) return j;  
  else if( trks_ptErr().at(i) > trks_ptErr().at(j) )                                       return i; 
  else if( trks_ptErr().at(i) < trks_ptErr().at(j) )                                       return j; 
  return i;

}

//--------------------------------------------------------------------
