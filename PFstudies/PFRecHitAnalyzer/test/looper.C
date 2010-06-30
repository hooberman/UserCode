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

#include "CMS2.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

bool debug = false;

using namespace tas;

//thresholds-----------------------------------------
float eb_threshold  = 0.07;
float ee_threshold  = 0.3;
float hb_threshold  = 0.7;
float he_threshold  = 0.8;
float hfh_threshold = 1.7;
float hfe_threshold = 0.8;
//--------------------------------------------------

void looper::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  bookHistos();

  // make a baby ntuple
  stringstream babyfilename;
  babyfilename << prefix << "_baby.root";
  MakeBabyNtuple( Form("%s_baby.root", prefix ) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  if(debug) cout << "Begin file loop" << endl;

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())){
    
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();
    ///nEvents = 100000;

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

    
      InitBabyNtuple();

      fillUnderOverFlow( hgenmet             ,    genmet()  );
      fillUnderOverFlow( hcalomet            ,    calomet()  );
      fillUnderOverFlow( hpfrmet_nothresh    ,    pfmet()  );

      float met_x    = 0;
      float met_y    = 0;
      float sumet    = 0;
      float ebmet_x  = 0;
      float ebmet_y  = 0;
      float ebsumet  = 0;
      float eemet_x  = 0;
      float eemet_y  = 0;
      float eesumet  = 0;
      float hbmet_x  = 0;
      float hbmet_y  = 0;
      float hbsumet  = 0;
      float hemet_x  = 0;
      float hemet_y  = 0;
      float hesumet  = 0;
      float hfhmet_x = 0;
      float hfhmet_y = 0;
      float hfhsumet = 0;
      float hfemet_x = 0;
      float hfemet_y = 0;
      float hfesumet = 0;
      
      //ECAL Barrel

      for( unsigned int i = 0 ; i < pf_ebrechit_e().size() ; ++i ){

        if( pf_ebrechit_e().at(i) < eb_threshold ) continue;

        sumet     += pf_ebrechit_et().at(i);
        ebsumet   += pf_ebrechit_et().at(i);
        met_x     -= pf_ebrechit_et().at(i) * cos( pf_ebrechit_phi().at(i) );
        met_y     -= pf_ebrechit_et().at(i) * sin( pf_ebrechit_phi().at(i) );
        ebmet_x   -= pf_ebrechit_et().at(i) * cos( pf_ebrechit_phi().at(i) );
        ebmet_y   -= pf_ebrechit_et().at(i) * sin( pf_ebrechit_phi().at(i) );
      }

      //ECAL Endcap

      for( unsigned int i = 0 ; i < pf_eerechit_e().size() ; ++i ){

        if( pf_eerechit_e().at(i) < ee_threshold ) continue;

        sumet     += pf_eerechit_et().at(i);
        eesumet   += pf_eerechit_et().at(i);
        met_x     -= pf_eerechit_et().at(i) * cos( pf_eerechit_phi().at(i) );
        met_y     -= pf_eerechit_et().at(i) * sin( pf_eerechit_phi().at(i) );
        eemet_x   -= pf_eerechit_et().at(i) * cos( pf_eerechit_phi().at(i) );
        eemet_y   -= pf_eerechit_et().at(i) * sin( pf_eerechit_phi().at(i) );
      }

      //HCAL Barrel

      for( unsigned int i = 0 ; i < pf_hbrechit_e().size() ; ++i ){

        if( pf_hbrechit_e().at(i) < hb_threshold ) continue;

        sumet     += pf_hbrechit_et().at(i);
        hbsumet   += pf_hbrechit_et().at(i);
        met_x     -= pf_hbrechit_et().at(i) * cos( pf_hbrechit_phi().at(i) );
        met_y     -= pf_hbrechit_et().at(i) * sin( pf_hbrechit_phi().at(i) );
        hbmet_x   -= pf_hbrechit_et().at(i) * cos( pf_hbrechit_phi().at(i) );
        hbmet_y   -= pf_hbrechit_et().at(i) * sin( pf_hbrechit_phi().at(i) );
      }

      //HCAL Endcap

      for( unsigned int i = 0 ; i < pf_herechit_e().size() ; ++i ){

        if( pf_herechit_e().at(i) < he_threshold ) continue;

        sumet     += pf_herechit_et().at(i);
        hesumet   += pf_herechit_et().at(i);
        met_x     -= pf_herechit_et().at(i) * cos( pf_herechit_phi().at(i) );
        met_y     -= pf_herechit_et().at(i) * sin( pf_herechit_phi().at(i) );
        hemet_x   -= pf_herechit_et().at(i) * cos( pf_herechit_phi().at(i) );
        hemet_y   -= pf_herechit_et().at(i) * sin( pf_herechit_phi().at(i) );
      }

      //HF Hadronic

      for( unsigned int i = 0 ; i < pf_hfhrechit_e().size() ; ++i ){

        if( pf_hfhrechit_e().at(i) < hfh_threshold ) continue;

        sumet      += pf_hfhrechit_et().at(i);
        hfhsumet   += pf_hfhrechit_et().at(i);
        met_x      -= pf_hfhrechit_et().at(i) * cos( pf_hfhrechit_phi().at(i) );
        met_y      -= pf_hfhrechit_et().at(i) * sin( pf_hfhrechit_phi().at(i) );
        hfhmet_x   -= pf_hfhrechit_et().at(i) * cos( pf_hfhrechit_phi().at(i) );
        hfhmet_y   -= pf_hfhrechit_et().at(i) * sin( pf_hfhrechit_phi().at(i) );
      }

      //HF EM

      for( unsigned int i = 0 ; i < pf_hferechit_e().size() ; ++i ){

        if( pf_hferechit_e().at(i) < hfe_threshold ) continue;

        sumet      += pf_hferechit_et().at(i);
        hfesumet   += pf_hferechit_et().at(i);
        met_x      -= pf_hferechit_et().at(i) * cos( pf_hferechit_phi().at(i) );
        met_y      -= pf_hferechit_et().at(i) * sin( pf_hferechit_phi().at(i) );
        hfemet_x   -= pf_hferechit_et().at(i) * cos( pf_hferechit_phi().at(i) );
        hfemet_y   -= pf_hferechit_et().at(i) * sin( pf_hferechit_phi().at(i) );
      }

      float pfrmet   = sqrt( pow( met_x , 2 ) + pow( met_y , 2 ) );
      
      fillUnderOverFlow( hpfrmet     ,    pfrmet  );
      
      FillBabyNtuple();
      
    } // end loop over events
  } // end loop over files
  

  
  CloseBabyNtuple();

  // make histos rootfile
  stringstream rootfilename;
  rootfilename << prefix << "_histos.root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(rootfilename.str().c_str());
  deleteHistos();
  
} // end ScanChain


void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  //declare histos
  hgenmet  = new TH1F("hgenmet","Gen-Level MET",100,0,100);
  hcalomet = new TH1F("hcalomet","caloMET",100,0,100);
  hpfrmet  = new TH1F("hpfrmet","PF RecHit MET",100,0,100);
  hpfrmet_nothresh  = new TH1F("hpfrmet_nothresh","PF RecHit MET (no thresholds)",100,0,100);


}


void looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("tree", "A Baby Ntuple");

  //declare branches
  babyTree_->Branch("run",          &run_,          "run/I"  );

  
}



void looper::FillBabyNtuple ()
{
  babyTree_->Fill();
}


void looper::InitBabyNtuple (){
}

void looper::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}


void looper::fillUnderOverFlow(TH1F *h1, float value, float weight)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}


// cout << "EE:     " << sqrt( pow( eemet_x , 2 )  + pow( eemet_y , 2 ) )  << "  " << pf_ee_met()  << endl;
// cout << "EB:     " << sqrt( pow( ebmet_x , 2 )  + pow( ebmet_y , 2 ) )  << "  " << pf_eb_met()  << endl;
// cout << "HE:     " << sqrt( pow( hemet_x , 2 )  + pow( hemet_y , 2 ) )  << "  " << pf_he_met()  << endl;
// cout << "HB:     " << sqrt( pow( hbmet_x , 2 )  + pow( hbmet_y , 2 ) )  << "  " << pf_hb_met()  << endl;
// cout << "HFH:    " << sqrt( pow( hfhmet_x , 2 ) + pow( hfhmet_y , 2 ) ) << "  " << pf_hfh_met() << endl;
// cout << "HFE:    " << sqrt( pow( hfemet_x , 2 ) + pow( hfemet_y , 2 ) ) << "  " << pf_hfe_met() << endl;
// cout << "TOT:    " << sqrt( pow( met_x , 2 )    + pow( met_y , 2 ) )    << "  " << pfmet()      << endl;
