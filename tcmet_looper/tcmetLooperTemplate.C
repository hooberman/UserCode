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

#include "../CORE/CMS2.h"
#include "../CORE/CMS2.cc"
#include "../CORE/metSelections.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/muonSelections.cc"
#include "../Tools/goodrun.cc"
#include "../CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

void tcmetLooperTemplate::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  bookHistos();

  //save met printouts to txt file
  stringstream s1;
  s1 << prefix << "_events.txt";
  ofile.open(s1.str().c_str());

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  // make a baby ntuple
  stringstream s2;
  s2 << prefix << "_baby.root";
  MakeBabyNtuple(s2.str().c_str());

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
  while ((currentFile = (TFile*)fileIter.Next()))
    {
      TFile f(currentFile->GetTitle());
      TTree *tree = (TTree*)f.Get("Events");
      cms2.Init(tree);

      // event loop
      unsigned int nEvents = tree->GetEntries();

      for (unsigned int event = 0; event < nEvents; ++event)
        {
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
          
          if( strcmp( prefix , "zee")   == 0 || 
              strcmp( prefix , "zmm")   == 0 || 
              strcmp( prefix , "ttbar") == 0 ){
            if (!isGoodDilepton())    continue;
          } 

          /*
          //for Z->ee event, require exactly 2 truth-matched electrons
          if (!isTruthZee())   continue;
          
          // for Z->mm event, require exactly 2 truth-matched muons
          if (!isTruthZmm())   continue;
          
          // for Z->ee/Z->mm event, veto if >=1 jet with corrected pt > 20 GeV, eta < 3
          if (jetVeto())       continue;
          */


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

          // pf met stuff
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
            genmet_     = cms2.gen_met();
            genmetphi_  = cms2.gen_metPhi();
            gensumet_   = cms2.gen_sumEt();
          }
          
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
          
          metStruct structMET = correctedTCMET(electronVetoCone, usePV, useHFcleaning, useHCALcleaning, useECALcleaning);
          
          //calomet
          met_       = cms2.evt_met();
          metphi_    = cms2.evt_metPhi();
          sumet_     = cms2.evt_sumet();

          //tcmet
          tcmet_     = structMET.met;
          tcmetphi_  = structMET.metphi;
          tcsumet_   = structMET.sumet;

          FillBabyNtuple();

          hmumet->Fill(mumet_);
          hmujesmet->Fill(mujesmet_);
          htcmet->Fill(tcmet_);
          hrawtcmet->Fill(rawtcmet_);

          hdmumet->Fill(mumet_-genmet_);
          hdmujesmet->Fill(mujesmet_-genmet_);
          hdtcmet->Fill(tcmet_-genmet_);
          hdrawtcmet->Fill(rawtcmet_-genmet_);

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
  stringstream s3;
  s3 << prefix << "_histos.root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(s3.str().c_str());
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
  
}


void tcmetLooperTemplate::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("tree", "A Baby Ntuple");

  babyTree_->Branch("run"  , &run_  , "run/I"  );
  babyTree_->Branch("lumi" , &lumi_ , "lumi/I" );
  babyTree_->Branch("event", &event_, "event/I");
  babyTree_->Branch("hase" , &hase_ , "hase/O" );

  babyTree_->Branch("pfmet"   , &pfmet_   , "pfmet/F"   );
  babyTree_->Branch("pfmetphi", &pfmetphi_, "pfmetphi/F");
  babyTree_->Branch("pfsumet" , &pfsumet_ , "pfsumet/F" );

  babyTree_->Branch("met"      , &met_      , "met/F"      );
  babyTree_->Branch("metphi"   , &metphi_   , "metphi/F"   );
  babyTree_->Branch("sumet"    , &sumet_    , "sumet/F"    );

  babyTree_->Branch("mumet"    , &mumet_    , "mumet/F"      );
  babyTree_->Branch("mumetphi" , &mumetphi_ , "mumetphi/F"   );
  babyTree_->Branch("musumet"  , &musumet_  , "musumet/F"    );

  babyTree_->Branch("mujesmet"      , &mujesmet_      , "mujesmet/F"      );
  babyTree_->Branch("mujesmetphi"   , &mujesmetphi_   , "mujesmetphi/F"   );
  babyTree_->Branch("mujessumet"    , &mujessumet_    , "mujessumet/F"    );

  babyTree_->Branch("genmet"      , &genmet_      , "genmet/F"   );
  babyTree_->Branch("genmetphi"   , &genmetphi_   , "genmetphi/F");
  babyTree_->Branch("gensumet"    , &gensumet_    , "gensumet/F" );

  babyTree_->Branch("tcmet"      , &tcmet_      , "tcmet/F"      );
  babyTree_->Branch("tcmetphi"   , &tcmetphi_   , "tcmetphi/F"   );
  babyTree_->Branch("tcsumet"    , &tcsumet_    , "tcsumet/F"    );

  babyTree_->Branch("rawtcmet"      , &rawtcmet_      , "rawtcmet/F"      );
  babyTree_->Branch("rawtcmetphi"   , &rawtcmetphi_   , "rawtcmetphi/F"   );
  babyTree_->Branch("rawtcsumet"    , &rawtcsumet_    , "rawtcsumet/F"    );

}



void tcmetLooperTemplate::FillBabyNtuple ()
{
  babyTree_->Fill();
}

void tcmetLooperTemplate::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}

bool tcmetLooperTemplate::isTruthZee () {
  
  if ( cms2.hyp_p4().size() != 1)              return false;
  if ( cms2.els_p4().size() != 2)              return false;
  if ( abs( cms2.els_mc3_id().at(0) ) != 11 )  return false; 
  if ( abs( cms2.els_mc3_id().at(1) ) != 11 )  return false; 

  return true;

}

bool tcmetLooperTemplate::isTruthZmm () {
  
  if ( cms2.hyp_p4().size() != 1)              return false;
  if ( cms2.mus_p4().size() != 2)              return false;
  if ( abs( cms2.mus_mc3_id().at(0) ) != 13 )  return false; 
  if ( abs( cms2.mus_mc3_id().at(1) ) != 13 )  return false; 

  return true;

}

bool tcmetLooperTemplate::jetVeto () {

  for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
    
    LorentzVector vjet = jets_p4().at(ijet);
    LorentzVector vlt  = hyp_lt_p4()[0];
    LorentzVector vll  = hyp_ll_p4()[0];

    if (dRbetweenVectors(vjet, vll) < 0.5) continue;
    if (dRbetweenVectors(vjet, vlt) < 0.5) continue;
    
    LorentzVector bah = jets_p4().at(ijet);
    
    if ( bah.pt() * jets_cor().at(ijet) > 20. && fabs( bah.eta() ) < 3.)
      return true;
  
  }
  
  return false;
  
}


bool tcmetLooperTemplate::isGoodDilepton ()
{
  //cout <<"tcmetLooperTemplate::isGoodDilepton"<<endl;
  bool isgood = false;

  for (unsigned int i = 0; i < cms2.hyp_p4().size(); ++i){
    
    if (cms2.hyp_lt_charge()[i] * cms2.hyp_ll_charge()[i] > 0)           continue;
    if (min(cms2.hyp_lt_p4()[i].pt(), cms2.hyp_ll_p4()[i].pt()) < 20.)   continue;
    
    //muon ID
    if (abs(hyp_ll_id()[i]) == 13  && (! (fabs(hyp_ll_p4()[i].eta()) < 2.4 && muonId(hyp_ll_index()[i]))))   continue;
    if (abs(hyp_lt_id()[i]) == 13  && (! (fabs(hyp_lt_p4()[i].eta()) < 2.4 && muonId(hyp_lt_index()[i]))))   continue;
    
    //cand01
    if (abs(hyp_ll_id()[i]) == 11  && (! electronSelection_cand01(hyp_ll_index()[i])))   continue;
    if (abs(hyp_lt_id()[i]) == 11  && (! electronSelection_cand01(hyp_lt_index()[i])))   continue;
    
    isgood = true;
  }
  
  return isgood;
}

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

      if (!electronSelection_cand01(cms2.hyp_lt_index()[i]))
        continue;

      if (!electronSelection_cand01(cms2.hyp_ll_index()[i]))
        continue;
      
      float dphi = fabs( cms2.hyp_lt_p4()[i].phi() - cms2.hyp_ll_p4()[i].phi() );
      if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;

      if( fabs( dphi - TMath::Pi() ) > 0.3 ) continue;

      isgood = true;
    }

  return isgood;
}


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

bool tcmetLooperTemplate::isGoodTrack( int index, bool usePV ) {
  
  float corrected_d0 = cms2.trks_d0corr().at(index);

  if(usePV)
    corrected_d0 = cms2.trks_d0vtx().at(index);

  if( cms2.trks_algo().at(index) < 8 ) {

    float d0cut = min( sqrt( pow(0.015,2) + pow(0.5/cms2.trks_trk_p4().at(index).pt(),2) ), 0.3);

    if( corrected_d0 > d0cut ) return false;
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



















/*


bool badeveto = false;

if( ectcmet_ - tcmet_ > 10 ){

ofile << endl << "---------------------------------------------------------------------------------" << endl;
            
int width = 10;

ofile << "|" << setw(width) << "type" << setw(width) 
<< "|" << setw(width) << "met"  << setw(width) 
<< "|" << setw(width) << "metx" << setw(width) 
<< "|" << setw(width) << "mety" << setw(width) << "|" << endl;
                 
         
ofile << "|" << setw(width) << "standard"        << setw(width) 
<< "|" << setw(width) << tcmet_            << setw(width) 
<< "|" << setw(width) << structMET.metx    << setw(width) 
<< "|" << setw(width) << structMET.mety    << setw(width) << "|" << endl; 
         
ofile << "|" << setw(width) << "eVeto"           << setw(width) 
<< "|" << setw(width) << ectcmet_          << setw(width) 
<< "|" << setw(width) << tcmetStruct.metx  << setw(width) 
<< "|" << setw(width) << tcmetStruct.mety  << setw(width) << "|" << endl; 
              
          
badeveto = true;
            
metStruct tempMet;
if (isData)
tempMet = correctedTCMET(false, false, true,  true,  true,  false, ofile);
else
tempMet = correctedTCMET(false, false, false, false, false, false, ofile);

}



if(badeveto){

ofile << " electron pt " << elspt  
<< " duplicate pt " << elpt_ 
<< " dR " << dR << endl;

}






*/
