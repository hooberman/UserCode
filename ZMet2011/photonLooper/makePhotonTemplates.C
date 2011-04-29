#include "makePhotonTemplates.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "../CORE/metTemplatesSelections.cc"
#include "TChain.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

bool debug = true;

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void makePhotonTemplates::ScanChain ( TChain* chain , char* iter ){

  
  int npass = 0;
  bookHistos();
  
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;

  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  if(debug) cout << "Begin file loop" << endl;

  // file loop
  TIter fileIter(listOfFiles);

  TFile* currentFile = 0;

  while ((currentFile = (TFile*)fileIter.Next())){

    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("T1");
    //cms2.Init(tree);
    setBranches(tree);

    // event loop
    //unsigned int nEvents = tree->GetEntries();
    nEvents = tree->GetEntries();

    for (unsigned int event = 0 ; event < nEvents; ++event){
   
      tree->GetEntry(event);
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

      //kinematic selection
      if( nJets_ < 2 )        continue;

      //cleaning
      if( jet_pt_ - etg_ < -5 )                                      continue; // cleaning
      if( fabs( etag_ ) > 2 )                                        continue; // photon eta < 2
      if( pfjetid_ != 1 )                                            continue; // pass PFJetID
      if( hlt20_ < 1 && hlt30_ < 1 && hlt50_ < 1 && hlt75_ < 1 )     continue; // require trig

      int iJetBin          = getJetBin       ( nJets_    );
      int iSumJetPtBin     = getSumJetPtBin  ( sumJetPt_ );
      int iBosonPtBin      = getBosonPtBin   ( etg_      );
      int iVtxBin          = getVtxBin       ( nvtx_     );
      float templateWeight = 1;

      //fill templates binned by njets, sumjetpt, boson pt        
      fillUnderOverFlow( tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  tcmet_    , templateWeight );
      fillUnderOverFlow( pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  pfmet_    , templateWeight );
      
      //fill templates binned by njets, sumjetpt, nVtx
      fillUnderOverFlow( tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ]    ,  tcmet_    , templateWeight );
      fillUnderOverFlow( pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ]    ,  pfmet_    , templateWeight );
    
      //fill templates binned by njets, sumjetpt
      fillUnderOverFlow( tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
      fillUnderOverFlow( pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );
      
      ++npass;

      int iTrigBin = -1;


      if( hlt75_ > 0 ){
        
        templateWeight = hlt75_;
        iTrigBin = 3;

        fillUnderOverFlow( tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );

        fillUnderOverFlow( hphotonPt70 , etg_ , templateWeight );
      }

      else if( hlt50_ > 0 ){
	
        templateWeight = hlt50_;
        iTrigBin = 2;

        fillUnderOverFlow( tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );

        fillUnderOverFlow( hphotonPt50 , etg_ , templateWeight );
      }

      else if( hlt30_ > 0 ){
        
        templateWeight = hlt30_;
        iTrigBin = 1;

        fillUnderOverFlow( tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );

        fillUnderOverFlow( hphotonPt30 , etg_ , templateWeight );
      }

      else if( hlt20_ > 0 ){
        
        templateWeight = hlt20_;
        iTrigBin = 0;

        fillUnderOverFlow( tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );

        fillUnderOverFlow( hphotonPt20 , etg_ , templateWeight );
      }




    } // end loop over events
  } // end loop over files
      
  cout << npass << " events passing selection" << endl;
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  //-------------------------------------------
  // normalize templates
  //-------------------------------------------
   
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        for( int iBosonPtBin = 0 ; iBosonPtBin < nBosonPtBins ; iBosonPtBin++ ){
          
          float scale = tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
          if( scale > 0 )  tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
          
          scale = pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
          if( scale > 0 )  pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
          
        }
      }
    }

    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        for( int iVtxBin = 0 ; iVtxBin < nVtxBins ; iVtxBin++ ){
          
          float scale = tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Integral();
          if( scale > 0 )  tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Scale ( 1. / scale );
          
          scale = pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Integral();
          if( scale > 0 )  pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Scale ( 1. / scale );
          
        }
      }
    }

    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
        float scale = tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      
        scale = pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
        
        
      }
    }

    for( int iTrigBin = 0 ; iTrigBin < 4 ; ++iTrigBin ){
      for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
        for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
          
          float scale = tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Integral();
          if( scale > 0 )  tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
          
          scale = pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Integral();
          if( scale > 0 )  pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
          
        }
      }
    }


    // make histos rootfile
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    cout << "Writing templates to " << Form("../templates/%s/photon_templates.root",iter) << endl;
    saveHist(Form("../templates/%s/photon_templates.root",iter));
    deleteHistos();
  
} // end ScanChain


void makePhotonTemplates::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

void makePhotonTemplates::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hphotonPt20  = new TH1F("hphotonPt20", "",500,0,500);
  hphotonPt30  = new TH1F("hphotonPt30", "",500,0,500);
  hphotonPt50  = new TH1F("hphotonPt50", "",500,0,500);
  hphotonPt70  = new TH1F("hphotonPt70", "",500,0,500);

  int maxmet = 200;

  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      for( int iBosonPtBin = 0 ; iBosonPtBin < nBosonPtBins ; iBosonPtBin++ ){
        
        tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("tcmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                     Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                          bosonPtString(iBosonPtBin).c_str()),maxmet,0,maxmet);
        
        pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("pfmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                     Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                          bosonPtString(iBosonPtBin).c_str()),maxmet,0,maxmet);
        
        tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();
        pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();
        
        tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
        pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
        
      }
    }
  }

  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      for( int iVtxBin = 0 ; iVtxBin < nVtxBins ; iVtxBin++ ){
        
        tcmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin] = new TH1F(Form("tcmetTemplate_njets_ht_nvtx_%i_%i_%i",iJetBin,iSumJetPtBin,iVtxBin),
                                                                               Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                                    nVTXString(iVtxBin).c_str()),maxmet,0,maxmet);
                
        pfmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin] = new TH1F(Form("pfmetTemplate_njets_ht_nvtx_%i_%i_%i",iJetBin,iSumJetPtBin,iVtxBin),
                                                                               Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                                    nVTXString(iVtxBin).c_str()),maxmet,0,maxmet);
        
        tcmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->Sumw2();
        pfmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->Sumw2();
        
        tcmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->GetXaxis()->SetTitle("tcmet (GeV)");
        pfmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
        
      }
    }
  }
  
  
  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      
      
      tcmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                               Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
      
      pfmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("pfmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                               Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
      
      tcmetTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
      pfmetTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
      
      tcmetTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
      pfmetTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
      
    }
  }
  
  char* trigName[4]={"HLT20","HLT30","HLT50","HLT75"};

  for( int iTrigBin = 0 ; iTrigBin < 4 ; iTrigBin++ ){
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
        
        tcmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetTemplate_photon_%i_%i_%i",iTrigBin,iJetBin,iSumJetPtBin),
									 Form("%s, %s, %s",trigName[iTrigBin],
									      jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
	
        pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin] = new TH1F(Form("pfmetTemplate_photon_%i_%i_%i",iTrigBin,iJetBin,iSumJetPtBin),
									 Form("%s, %s, %s",trigName[iTrigBin],
									      jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
        
        tcmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->Sumw2();
        pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->Sumw2();
        
        tcmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
        pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
        
      }
    }
  }

}


void makePhotonTemplates::setBranches (TTree* tree){
  
  tree->SetBranchAddress("tcmet"	,       &tcmet_         );
  tree->SetBranchAddress("pfmet"	,       &pfmet_         );
  tree->SetBranchAddress("njets"	,       &nJets_         );
  tree->SetBranchAddress("sumjetpt"	,       &sumJetPt_      );
  tree->SetBranchAddress("nvtx"		,       &nvtx_          );
  tree->SetBranchAddress("jetpt"	,       &jet_pt_        );  
  tree->SetBranchAddress("etag"		,       &etag_          ); 
  tree->SetBranchAddress("etg"		,       &etg_           ); 
  tree->SetBranchAddress("pfjetid"	,       &pfjetid_       );     
  tree->SetBranchAddress("hlt20"	,	&hlt20_		);
  tree->SetBranchAddress("hlt30"	,	&hlt30_		);
  tree->SetBranchAddress("hlt50"	,	&hlt50_		);
  tree->SetBranchAddress("hlt75"	,	&hlt75_		);
  tree->SetBranchAddress("hlt125"	,	&hlt125_	);

}
