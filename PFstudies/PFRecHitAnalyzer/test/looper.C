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

//#include "CMS2.cc"
#include "CMS2.h"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
#include "TProfile.h"

bool debug = false;

using namespace tas;


inline double fround(double n, double d)
{
  return floor(n * pow(10., d) + .5) / pow(10., d);
}


looper::metStruct looper::ScanChain (TChain* chain, const char* prefix, const char* suffix, bool isData, int nEvents,
                                     float eb_threshold, float ee_threshold,  float hb_threshold,
                                     float he_threshold, float hfh_threshold, float hfe_threshold,
                                     float hfshort_threshold, float hflong_threshold){
  
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

      //met
      fillUnderOverFlow( hgenmet             ,    genmet()     );
      fillUnderOverFlow( hcalo_met            ,    calomet()    );
      fillUnderOverFlow( hpfr_met_nothresh    ,    pfmet()      );
      fillUnderOverFlow( hpfc_met             ,    pfclusmet()  );

      //sumet
      fillUnderOverFlow( hgensumet             ,    gensumet()     );
      fillUnderOverFlow( hcalo_sumet            ,    calosumet()    );
      fillUnderOverFlow( hpfr_sumet_nothresh    ,    pfsumet()      );
      fillUnderOverFlow( hpfc_sumet             ,    pfclussumet()  );

      //sumet in subdetectors (calomet)
      fillUnderOverFlow( hcalo_ebsumet      ,    calo_eb_sumet()    );
      fillUnderOverFlow( hcalo_eesumet      ,    calo_ee_sumet()    );
      fillUnderOverFlow( hcalo_hbsumet      ,    calo_hb_sumet()    );
      fillUnderOverFlow( hcalo_hesumet      ,    calo_he_sumet()    );
      fillUnderOverFlow( hcalo_hfhsumet     ,    calo_hfh_sumet()   );
      fillUnderOverFlow( hcalo_hfesumet     ,    calo_hfe_sumet()   );
      
      //sumet in subdetectors (PFClusters)
      fillUnderOverFlow( hpfc_ebsumet      ,    pfclus_eb_sumet()   );
      fillUnderOverFlow( hpfc_eesumet      ,    pfclus_ee_sumet()   );
      fillUnderOverFlow( hpfc_hbsumet      ,    pfclus_hb_sumet()   );
      fillUnderOverFlow( hpfc_hesumet      ,    pfclus_he_sumet()   );
      fillUnderOverFlow( hpfc_hfhsumet     ,    pfclus_hfh_sumet()  );
      fillUnderOverFlow( hpfc_hfesumet     ,    pfclus_hfe_sumet()  );

      //met in subdetectors (calomet)
//       fillUnderOverFlow( hcalo_ebmet      ,    calo_eb_met()    );
//       fillUnderOverFlow( hcalo_eemet      ,    calo_ee_met()    );
//       fillUnderOverFlow( hcalo_hbmet      ,    calo_hb_met()    );
//       fillUnderOverFlow( hcalo_hemet      ,    calo_he_met()    );
//       fillUnderOverFlow( hcalo_hfhmet     ,    calo_hfh_met()   );
//       fillUnderOverFlow( hcalo_hfemet     ,    calo_hfe_met()   );

      //met in subdetectors (PFClusters)
      fillUnderOverFlow( hpfc_ebmet      ,    pfclus_eb_met()   );
      fillUnderOverFlow( hpfc_eemet      ,    pfclus_ee_met()   );
      fillUnderOverFlow( hpfc_hbmet      ,    pfclus_hb_met()   );
      fillUnderOverFlow( hpfc_hemet      ,    pfclus_he_met()   );
      fillUnderOverFlow( hpfc_hfhmet     ,    pfclus_hfh_met()  );
      fillUnderOverFlow( hpfc_hfemet     ,    pfclus_hfe_met()  );

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

        fillUnderOverFlow( hebe ,  pf_ebrechit_e().at(i) );
        if( pf_ebrechit_e().at(i) < eb_threshold ) continue;

        sumet     += pf_ebrechit_et().at(i);
        ebsumet   += pf_ebrechit_et().at(i);
        met_x     -= pf_ebrechit_et().at(i) * cos( pf_ebrechit_phi().at(i) );
        met_y     -= pf_ebrechit_et().at(i) * sin( pf_ebrechit_phi().at(i) );
        ebmet_x   -= pf_ebrechit_et().at(i) * cos( pf_ebrechit_phi().at(i) );
        ebmet_y   -= pf_ebrechit_et().at(i) * sin( pf_ebrechit_phi().at(i) );
      }

      float ebmet = sqrt( ebmet_x * ebmet_x + ebmet_y * ebmet_y );

      //ECAL Endcap

      for( unsigned int i = 0 ; i < pf_eerechit_e().size() ; ++i ){

        fillUnderOverFlow( heee ,  pf_eerechit_e().at(i) );
        if( pf_eerechit_e().at(i) < ee_threshold ) continue;

        sumet     += pf_eerechit_et().at(i);
        eesumet   += pf_eerechit_et().at(i);
        met_x     -= pf_eerechit_et().at(i) * cos( pf_eerechit_phi().at(i) );
        met_y     -= pf_eerechit_et().at(i) * sin( pf_eerechit_phi().at(i) );
        eemet_x   -= pf_eerechit_et().at(i) * cos( pf_eerechit_phi().at(i) );
        eemet_y   -= pf_eerechit_et().at(i) * sin( pf_eerechit_phi().at(i) );
      }

      float eemet = sqrt( eemet_x * eemet_x + eemet_y * eemet_y );

      //HCAL Barrel

      for( unsigned int i = 0 ; i < pf_hbrechit_e().size() ; ++i ){

        fillUnderOverFlow( hhbe ,  pf_hbrechit_e().at(i) );
        if( pf_hbrechit_e().at(i) < hb_threshold ) continue;

        sumet     += pf_hbrechit_et().at(i);
        hbsumet   += pf_hbrechit_et().at(i);
        met_x     -= pf_hbrechit_et().at(i) * cos( pf_hbrechit_phi().at(i) );
        met_y     -= pf_hbrechit_et().at(i) * sin( pf_hbrechit_phi().at(i) );
        hbmet_x   -= pf_hbrechit_et().at(i) * cos( pf_hbrechit_phi().at(i) );
        hbmet_y   -= pf_hbrechit_et().at(i) * sin( pf_hbrechit_phi().at(i) );
      }

      float hbmet = sqrt( hbmet_x * hbmet_x + hbmet_y * hbmet_y );

      //HCAL Endcap

      for( unsigned int i = 0 ; i < pf_herechit_e().size() ; ++i ){

        fillUnderOverFlow( hhee ,  pf_herechit_e().at(i) );
        if( pf_herechit_e().at(i) < he_threshold ) continue;

        sumet     += pf_herechit_et().at(i);
        hesumet   += pf_herechit_et().at(i);
        met_x     -= pf_herechit_et().at(i) * cos( pf_herechit_phi().at(i) );
        met_y     -= pf_herechit_et().at(i) * sin( pf_herechit_phi().at(i) );
        hemet_x   -= pf_herechit_et().at(i) * cos( pf_herechit_phi().at(i) );
        hemet_y   -= pf_herechit_et().at(i) * sin( pf_herechit_phi().at(i) );
      }

      float hemet = sqrt( hemet_x * hemet_x + hemet_y * hemet_y );


      assert( pf_hfhrechit_e().size() ==  pf_hferechit_e().size() );
      
      //HF Hadronic

      for( unsigned int i = 0 ; i < pf_hfhrechit_e().size() ; ++i ){

        fillUnderOverFlow( hhfhe ,  pf_hfhrechit_e().at(i) );
        if( pf_hfhrechit_e().at(i) < hfh_threshold ) continue;

        float eshort = pf_hfhrechit_e().at(i) / 2.;

        if( eshort < hfshort_threshold ) continue;

        sumet      += pf_hfhrechit_et().at(i);
        hfhsumet   += pf_hfhrechit_et().at(i);
        met_x      -= pf_hfhrechit_et().at(i) * cos( pf_hfhrechit_phi().at(i) );
        met_y      -= pf_hfhrechit_et().at(i) * sin( pf_hfhrechit_phi().at(i) );
        hfhmet_x   -= pf_hfhrechit_et().at(i) * cos( pf_hfhrechit_phi().at(i) );
        hfhmet_y   -= pf_hfhrechit_et().at(i) * sin( pf_hfhrechit_phi().at(i) );
      }

      float hfhmet = sqrt( hfhmet_x * hfhmet_x + hfhmet_y * hfhmet_y );

      //HF EM

      for( unsigned int i = 0 ; i < pf_hferechit_e().size() ; ++i ){

        fillUnderOverFlow( hhfee ,  pf_hferechit_e().at(i) );
        if( pf_hferechit_e().at(i) < hfe_threshold ) continue;

        float eshort = pf_hfhrechit_e().at(i) / 2.;
        float elong = eshort + pf_hferechit_e().at(i);

        if( eshort < hfshort_threshold ) continue;
        if( elong  < hflong_threshold  ) continue;
        
        sumet      += pf_hferechit_et().at(i);
        hfesumet   += pf_hferechit_et().at(i);
        met_x      -= pf_hferechit_et().at(i) * cos( pf_hferechit_phi().at(i) );
        met_y      -= pf_hferechit_et().at(i) * sin( pf_hferechit_phi().at(i) );
        hfemet_x   -= pf_hferechit_et().at(i) * cos( pf_hferechit_phi().at(i) );
        hfemet_y   -= pf_hferechit_et().at(i) * sin( pf_hferechit_phi().at(i) );
      }
      
      float hfemet = sqrt( hfemet_x * hfemet_x + hfemet_y * hfemet_y );

      float pfrmet   = sqrt( pow( met_x , 2 ) + pow( met_y , 2 ) );
      
      //met/sumet (PFRecHits)
      fillUnderOverFlow( hpfr_met       ,    pfrmet  );
      fillUnderOverFlow( hpfr_sumet     ,    sumet   );

      //sumet in subdetectors (PFRecHits)
      fillUnderOverFlow( hpfr_ebsumet      ,    ebsumet   );
      fillUnderOverFlow( hpfr_eesumet      ,    eesumet   );
      fillUnderOverFlow( hpfr_hbsumet      ,    hbsumet   );
      fillUnderOverFlow( hpfr_hesumet      ,    hesumet   );
      fillUnderOverFlow( hpfr_hfhsumet     ,    hfhsumet  );
      fillUnderOverFlow( hpfr_hfesumet     ,    hfesumet  );

      //met in subdetectors (PFRecHits)
      fillUnderOverFlow( hpfr_ebmet      ,    ebmet   );
      fillUnderOverFlow( hpfr_eemet      ,    eemet   );
      fillUnderOverFlow( hpfr_hbmet      ,    hbmet   );
      fillUnderOverFlow( hpfr_hemet      ,    hemet   );
      fillUnderOverFlow( hpfr_hfhmet     ,    hfhmet  );
      fillUnderOverFlow( hpfr_hfemet     ,    hfemet  );

      //met x/y components (PFRecHits)
      fillUnderOverFlow( hdpfr_met         ,    pfrmet - genmet()   );
      fillUnderOverFlow( hpfr_metxy        ,    met_x     );
      fillUnderOverFlow( hpfr_metxy        ,    met_y     );
      fillUnderOverFlow( hpfr_ebmetxy      ,    ebmet_x   );
      fillUnderOverFlow( hpfr_ebmetxy      ,    ebmet_y   );
      fillUnderOverFlow( hpfr_eemetxy      ,    eemet_x   );
      fillUnderOverFlow( hpfr_eemetxy      ,    eemet_y   );
      fillUnderOverFlow( hpfr_hbmetxy      ,    hbmet_x   );
      fillUnderOverFlow( hpfr_hbmetxy      ,    hbmet_y   );
      fillUnderOverFlow( hpfr_hemetxy      ,    hemet_x   );
      fillUnderOverFlow( hpfr_hemetxy      ,    hemet_y   );
      fillUnderOverFlow( hpfr_hfhmetxy     ,    hfhmet_x  );
      fillUnderOverFlow( hpfr_hfhmetxy     ,    hfhmet_y  );
      fillUnderOverFlow( hpfr_hfemetxy     ,    hfemet_x  );
      fillUnderOverFlow( hpfr_hfemetxy     ,    hfemet_y  );
      
      tmet->Fill( gensumet() , sumet / gensumet() );
      

      FillBabyNtuple();
      
    } // end loop over events
  } // end loop over files
  

  CloseBabyNtuple();

  metStruct mystruct;

  mystruct.dmet_mean  = hdpfr_met->GetMean(1);
  mystruct.dmet_rms   = hdpfr_met->GetRMS(1);
  mystruct.met_rms    = hpfr_metxy->GetRMS(1);

  mystruct.ebmet_rms     = hpfr_ebmetxy->GetRMS(1);
  mystruct.eemet_rms     = hpfr_eemetxy->GetRMS(1);
  mystruct.hbmet_rms     = hpfr_hbmetxy->GetRMS(1);
  mystruct.hemet_rms     = hpfr_hemetxy->GetRMS(1);
  mystruct.hfhmet_rms    = hpfr_hfhmetxy->GetRMS(1);
  mystruct.hfemet_rms    = hpfr_hfemetxy->GetRMS(1);

  // make histos rootfile
  stringstream rootfilename;
  rootfilename << "root/" << prefix << suffix << "_histos.root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(rootfilename.str().c_str());
  deleteHistos();
 
  
  return mystruct;
  
} // end ScanChain


void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  
  float max1 = 500;

  //declare histos
  hgenmet  = new TH1F("hgenmet","Gen-Level MET",100,0,max1);
  hcalo_met = new TH1F("hcalo_met","caloMET",100,0,max1);
  hpfr_met  = new TH1F("hpfr_met","PF RecHit MET",100,0,max1);
  hpfc_met  = new TH1F("hpfc_met","PF Cluster MET",100,0,max1);
  hpfr_met_nothresh  = new TH1F("hpfr_met_nothresh","PF RecHit MET (no thresholds)",100,0,max1);

  hgensumet  = new TH1F("hgensumet","Gen-Level sumET",100,0,1000);
  hcalo_sumet = new TH1F("hcalo_sumet","calo sumET",100,0,1000);
  hpfr_sumet  = new TH1F("hpfr_sumet","PF RecHit sumET",100,0,1000);
  hpfc_sumet  = new TH1F("hpfc_sumet","PF Cluster sumET",100,0,1000);
  hpfr_sumet_nothresh  = new TH1F("hpfr_sumet_nothresh","PF RecHit sumET (no thresholds)",100,0,1000);

  hcalo_ebsumet      = new TH1F("hcalo_ebsumet","EB caloMET", 100,0,max1);
  hcalo_eesumet      = new TH1F("hcalo_eesumet","EE caloMET", 100,0,max1);
  hcalo_hbsumet      = new TH1F("hcalo_hbsumet","HB caloMET", 100,0,max1);
  hcalo_hesumet      = new TH1F("hcalo_hesumet","HE caloMET", 100,0,max1);
  hcalo_hfhsumet     = new TH1F("hcalo_hfhsumet","HFH caloMET",100,0,max1);
  hcalo_hfesumet     = new TH1F("hcalo_hfesumet","HFE caloMET",100,0,max1);

  hcalo_ebmet      = new TH1F("hcalo_ebmet","EB caloMET", 100,0,max1);
  hcalo_eemet      = new TH1F("hcalo_eemet","EE caloMET", 100,0,max1);
  hcalo_hbmet      = new TH1F("hcalo_hbmet","HB caloMET", 100,0,max1);
  hcalo_hemet      = new TH1F("hcalo_hemet","HE caloMET", 100,0,max1);
  hcalo_hfhmet     = new TH1F("hcalo_hfhmet","HFH caloMET",100,0,max1);
  hcalo_hfemet     = new TH1F("hcalo_hfemet","HFE caloMET",100,0,max1);

  hpfr_ebsumet      = new TH1F("hpfr_ebsumet","EB pfrsumET", 100,0,max1);
  hpfr_eesumet      = new TH1F("hpfr_eesumet","EE pfrsumET", 100,0,max1);
  hpfr_hbsumet      = new TH1F("hpfr_hbsumet","HB pfrsumET", 100,0,max1);
  hpfr_hesumet      = new TH1F("hpfr_hesumet","HE pfrsumET", 100,0,max1);
  hpfr_hfhsumet     = new TH1F("hpfr_hfhsumet","HFH pfrsumET",100,0,max1);
  hpfr_hfesumet     = new TH1F("hpfr_hfesumet","HFE pfrsumET",100,0,max1);

  hpfr_ebmet      = new TH1F("hpfr_ebmet","EB pfrmet", 100,0,max1);
  hpfr_eemet      = new TH1F("hpfr_eemet","EE pfrmet", 100,0,max1);
  hpfr_hbmet      = new TH1F("hpfr_hbmet","HB pfrmet", 100,0,max1);
  hpfr_hemet      = new TH1F("hpfr_hemet","HE pfrmet", 100,0,max1);
  hpfr_hfhmet     = new TH1F("hpfr_hfhmet","HFH pfrmet",100,0,max1);
  hpfr_hfemet     = new TH1F("hpfr_hfemet","HFE pfrmet",100,0,max1);

  hpfc_ebsumet      = new TH1F("hpfc_ebsumet","EB pfcsumET", 100,0,max1);
  hpfc_eesumet      = new TH1F("hpfc_eesumet","EE pfcsumET", 100,0,max1);
  hpfc_hbsumet      = new TH1F("hpfc_hbsumet","HB pfcsumET", 100,0,max1);
  hpfc_hesumet      = new TH1F("hpfc_hesumet","HE pfcsumET", 100,0,max1);
  hpfc_hfhsumet     = new TH1F("hpfc_hfhsumet","HFH pfcsumET",100,0,max1);
  hpfc_hfesumet     = new TH1F("hpfc_hfesumet","HFE pfcsumET",100,0,max1);

  hpfc_ebmet      = new TH1F("hpfc_ebmet","EB pfcmet", 100,0,max1);
  hpfc_eemet      = new TH1F("hpfc_eemet","EE pfcmet", 100,0,max1);
  hpfc_hbmet      = new TH1F("hpfc_hbmet","HB pfcmet", 100,0,max1);
  hpfc_hemet      = new TH1F("hpfc_hemet","HE pfcmet", 100,0,max1);
  hpfc_hfhmet     = new TH1F("hpfc_hfhmet","HFH pfcmet",100,0,max1);
  hpfc_hfemet     = new TH1F("hpfc_hfemet","HFE pfcmet",100,0,max1);

  hpfr_metxy        = new TH1F("hpfr_metxy",   "pfrMET(x/y)",    200,-max1,max1);
  hpfr_ebmetxy      = new TH1F("hpfr_ebmetxy", "EB pfrMET(x/y)", 200,-max1,max1);
  hpfr_eemetxy      = new TH1F("hpfr_eemetxy", "EE pfrMET(x/y)", 200,-max1,max1);
  hpfr_hbmetxy      = new TH1F("hpfr_hbmetxy", "HB pfrMET(x/y)", 200,-max1,max1);
  hpfr_hemetxy      = new TH1F("hpfr_hemetxy", "HE pfrMET(x/y)", 200,-max1,max1);
  hpfr_hfhmetxy     = new TH1F("hpfr_hfhmetxy","HFH pfrMET(x/y)",200,-max1,max1);
  hpfr_hfemetxy     = new TH1F("hpfr_hfemetxy","HFE pfrMET(x/y)",200,-max1,max1);
  hdpfr_met         = new TH1F("hdpfr_met",    "pfrMET - genMET",200,-max1,max1);

  hebe             = new TH1F("hebe", "EB rechit energy", 1000,0,10);
  heee             = new TH1F("heee", "EE rechit energy", 1000,0,10);
  hhbe             = new TH1F("hhbe", "HB rechit energy", 1000,0,10);
  hhee             = new TH1F("hhee", "HE rechit energy", 1000,0,10);
  hhfhe            = new TH1F("hhfhe","HFH rechit energy",1000,0,10);
  hhfee            = new TH1F("hhfee","HFE rechit energy",1000,0,10);

//   hclusebe             = new TH1F("hclusebe", "EB cluster energy", 1000,0,10);
//   hcluseee             = new TH1F("hcluseee", "EE cluster energy", 1000,0,10);
//   hclushbe             = new TH1F("hclushbe", "HB cluster energy", 1000,0,10);
//   hclushee             = new TH1F("hclushee", "HE cluster energy", 1000,0,10);
//   hclushfhe            = new TH1F("hclushfhe","HFH cluster energy",1000,0,10);
//   hclushfee            = new TH1F("hclushfee","HFE cluster energy",1000,0,10);
  
  tmet             = new TProfile("tmet","",200,0,2000,0,2);
  tmet->GetXaxis()->SetTitle("gen sumET (GeV)");
  tmet->GetYaxis()->SetTitle("pfrechit sumET / gen sumET");
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








//thresholds-----------------------------------------

//calotower-based
// float eb_threshold  = 0.18;
// float ee_threshold  = 0.6;
// float hb_threshold  = 0.8;
// float he_threshold  = 1.5;
// float hfh_threshold = 1.7;
// float hfe_threshold = 8.;

//pfrechit-based
// float eb_threshold  = 0.2;
// float ee_threshold  = 0.45;
// float hb_threshold  = 0.4;
// float he_threshold  = 0.4;
// float hfh_threshold = 2.4;
// float hfe_threshold = 0.8;

//--------------------------------------------------

/*
  cout << endl << endl 
       << "----------------------------------------------------------" << endl;

  cout << "|" << setw(15) << ""           << setw(4)
       << "|" << setw(15) << "Mean"       << setw(4)
       << "|" << setw(15) << "RMS"        << setw(4) << "|"  << endl; 

  cout << "|" << setw(15) << "d(met)" << setw(4)
       << "|" << setw(15) << fround(hdpfrmet->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hdpfrmet->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_metxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_metxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "EB met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_ebmetxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_ebmetxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "EE met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_eemetxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_eemetxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "HB met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_hbmetxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_hbmetxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "HE met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_hemetxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_hemetxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "HFH met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_hfhmetxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_hfhmetxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "|" << setw(15) << "HFE met(x/y)" << setw(4)
       << "|" << setw(15) << fround(hpfr_hfemetxy->GetMean(1),2) << setw(4)
       << "|" << setw(15) << fround(hpfr_hfemetxy->GetRMS(1),1)  << setw(4) << "|" << endl; 

  cout << "----------------------------------------------------------" << endl << endl << endl;
*/
