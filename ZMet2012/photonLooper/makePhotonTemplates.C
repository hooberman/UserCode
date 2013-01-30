#include "makePhotonTemplates.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#ifndef __CINT__
#include "../CORE/metTemplatesSelections.cc"
#endif

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

const bool debug          =  true;
const bool vtxreweight    =  true;
const bool bveto          =  true;
const bool mjjcut         =  true;
const bool pt40           = false;

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void makePhotonTemplates::ScanChain ( TChain* chain , char* iter , char* sample ){

  bool useHGGTriggers = true;
  //if( TString(sample).Contains("DoubleElectron") ) useHGGTriggers = true;

  cout << "Sample : " << sample << endl;
  if( useHGGTriggers ) cout << "Using H->gg triggers" << endl;
  else                 cout << "Using standard triggers" << endl;

  if( vtxreweight ) cout << "Doing vtx reweighting" << endl;
  else              cout << "NO vtx reweighting"    << endl;

  if( bveto )       cout << "Doing b-veto" << endl;
  else              cout << "NO b-veto"    << endl;

  if( mjjcut )      cout << "Doing mjj cut" << endl;
  else              cout << "NO mjj cut"    << endl;

  if( pt40 )        cout << "Require >=2 40 GeV jets with HT > 100 GeV"        << endl;
  else              cout << "DO NOT require >=2 40 GeV jets with HT > 100 GeV" << endl;

  int npass = 0;
  bookHistos();
  
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;

  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  if(debug) cout << "Begin file loop" << endl;

  TH1F* reweightHist[5];

  if( vtxreweight ){ 

    char* vtxfile = (char*) "vtxreweight_DoubleElectron_19fb.root";
    //if( TString(sample).Contains("DoubleElectron") ) vtxfile = (char*) "vtxreweight_DoubleElectron_9p2fb.root";
    //if( TString(sample).Contains("DoubleElectron") ) vtxfile = (char*) "vtxreweight_DoubleElectron_9p2fb_2012AB.root";
    //if( TString(sample).Contains("DoubleElectron") ) vtxfile = (char*) "vtxreweight_DoubleElectron_9p2fb_2012C.root";
    //if( TString(sample).Contains("DoubleElectron") ) vtxfile = (char*) "vtxreweight_DoubleElectron_19fb.root";

    cout << "Using vtx reweighting file " << vtxfile << endl;

    TFile *reweightFile = reweightFile = TFile::Open(vtxfile);
    reweightHist[0] = (TH1F*) reweightFile->Get("hratio20");
    reweightHist[1] = (TH1F*) reweightFile->Get("hratio30");
    reweightHist[2] = (TH1F*) reweightFile->Get("hratio50");
    reweightHist[3] = (TH1F*) reweightFile->Get("hratio70");
    reweightHist[4] = (TH1F*) reweightFile->Get("hratio90");
  }

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

      int vtxbin = -1;
      if( nvtx_ < 5 )                 vtxbin = 1;
      if( nvtx_ >=  5 && nvtx_ < 10 ) vtxbin = 2;
      if( nvtx_ >= 10 && nvtx_ < 15 ) vtxbin = 3;
      if( nvtx_ >= 15 && nvtx_ < 20 ) vtxbin = 4;
      if( nvtx_ >= 20 && nvtx_ < 25 ) vtxbin = 5;
      if( nvtx_ >= 25 && nvtx_ < 30 ) vtxbin = 6;
      if( nvtx_ >= 30 )               vtxbin = 7;

      int h20 = hlt20_;
      int h30 = hlt30_;
      int h50 = hlt50_;
      int h75 = hlt75_;
      int h90 = hlt90_;
      
      if( useHGGTriggers ){
	h20 = hgg22_;
	h30 = hgg36_;
	h50 = hgg50_;
	h75 = hgg75_;
	h90 = hgg90_;
      }

      //float mjj = ( *jet1_ + *jet2_ ).mass();
      //if( mjjcut && ( mjj < 70.0 || mjj > 110.0 ) )         continue; // apply dijet mass cut

      // event selection 
      if( nJets_ < 2 )                                      continue; // >=2 jets
      if( etg_ <  20 )                                      continue; // photon pt  > 20 GeV
      if( fabs( etag_ ) > 2 )                               continue; // photon eta < 2
      if( hoe_ > 0.1 )                                      continue; // H/E < 0.1
      if( photon_pixelseed_ == 1 )                          continue; // veto pixel match
      if( jetneutralemfrac_ < 0.7 )                         continue; // jet neutral EM fraction cut
      if( jet_pt_     - etg_ < -5 )                         continue; // pfjet cleaning
      //if( calojet_pt_ - etg_ < -5 )                         continue; // calojet cleaning
      if( elveto_ == 1 )                                    continue; // remove photons with nearby electrons
      if( maxleppt_ > 20.0 )                                continue; // veto leptons pt > 20 GeV
      if( acos( cos( phig_ - pfmetphi_ ) ) < 0.14 )         continue; // kill photons aligned with MET
      if( bveto && nbm_ > 0 )                               continue; // apply b-veto 
      //if( bveto && nbl_ > 0 )                               continue; // apply b-veto 
      if( pt40 && ( nJets40_ < 2 || ht40_ < 100.0 ) )       continue; // require 2 pT > 40 GeV jets with HT > 100 GeV
      if( !(csc_==0 && hbhe_==1 && hcallaser_==1 && ecaltp_==1 && trkfail_==1 && eebadsc_==1 && hbhenew_==1) ) continue; // MET filters
      //if( run_ > 196531 )                                   continue; // 2012 A+B
      //if( run_ <= 196531 )                                  continue; // 2012 C
      //if( run_ >= 197556 && run_ <= 198913 )                continue; // veto 2012C-PromptReco-v1
      // //if( pfjetid_ != 1 )                                                     continue; // pass PFJetID
      
      //if( h20 < 1 && h30 < 1 && h50 < 1 && h75 < 1 && h90 < 1 )                    continue; // require trig
      if( h20 < 1 && h50 < 1 && h75 < 1 && h90 < 1 )                    continue; // require trig

      int iJetBin;
      int iSumJetPtBin;

      if( pt40 ){
	iJetBin          = getJetBin       ( nJets40_    );
	iSumJetPtBin     = getSumJetPtBin  ( ht40_       );
      }
      else{
	iJetBin          = getJetBin       ( nJets_  );
	iSumJetPtBin     = getSumJetPtBin  ( ht_     );
      }

      int iBosonPtBin      = getBosonPtBin   ( etg_      );
      int iVtxBin          = getVtxBin       ( nvtx_     );
      float templateWeight = 1;

      /*
      //fill templates binned by njets, sumjetpt, boson pt        
      fillUnderOverFlow( tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  tcmet_    , templateWeight );
      fillUnderOverFlow( pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  pfmet_    , templateWeight );
      
      //fill templates binned by njets, sumjetpt, nVtx
      fillUnderOverFlow( tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ]    ,  tcmet_    , templateWeight );
      fillUnderOverFlow( pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ]    ,  pfmet_    , templateWeight );
    
      //fill templates binned by njets, sumjetpt
      fillUnderOverFlow( tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
      fillUnderOverFlow( pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );
      */

      ++npass;

      int iTrigBin = -1;

      if( h20 > 0 )  fillUnderOverFlow( hphotonPt20 , etg_ , h20 );
      if( h30 > 0 )  fillUnderOverFlow( hphotonPt30 , etg_ , h30 );
      if( h50 > 0 )  fillUnderOverFlow( hphotonPt50 , etg_ , h50 );
      if( h75 > 0 )  fillUnderOverFlow( hphotonPt70 , etg_ , h75 );
      if( h90 > 0 )  fillUnderOverFlow( hphotonPt90 , etg_ , h90 );

      if     ( h90 > 0 )  templateWeight = h90;
      else if( h75 > 0 )  templateWeight = h75;
      else if( h50 > 0 )  templateWeight = h50;
      else if( h20 > 0 )  templateWeight = h20;

      if( etg_ > 90.0 ){
	//if( h90 > 0 ){
        //templateWeight = h90;
        iTrigBin = 4;

	fillUnderOverFlow( hphotonPt90_exc , etg_ , templateWeight );
	fillUnderOverFlow( hphotonAll , etg_  , templateWeight );
	fillUnderOverFlow( hnvtxPt90  , nvtx_ , templateWeight );
	fillUnderOverFlow( hnvtxAll   , nvtx_ , templateWeight );

	if( vtxreweight ) templateWeight *= reweightHist[4]->GetBinContent(vtxbin);

	//cout << "nvtx vtxbin weight " << nvtx_ << " " << vtxbin << " " << reweightHist[4]->GetBinContent(vtxbin) << endl;

        fillUnderOverFlow(   tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_         , templateWeight );
        fillUnderOverFlow(   pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_         , templateWeight );
        fillUnderOverFlow( t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmett1_       , templateWeight );
        fillUnderOverFlow( t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] ,  pfmett1new_    , templateWeight );
      }
      
      else if( etg_ > 75.0 ){
	//else if( h75 > 0 ){
        //templateWeight = h75;
        iTrigBin = 3;

	fillUnderOverFlow( hphotonPt70_exc , etg_ , templateWeight );
	fillUnderOverFlow( hphotonAll , etg_  , templateWeight );
	fillUnderOverFlow( hnvtxPt70  , nvtx_ , templateWeight );
	fillUnderOverFlow( hnvtxAll   , nvtx_ , templateWeight );

	if( vtxreweight ) templateWeight *= reweightHist[3]->GetBinContent(vtxbin);

        fillUnderOverFlow(   tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_         , templateWeight );
        fillUnderOverFlow(   pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_         , templateWeight );
        fillUnderOverFlow( t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmett1_       , templateWeight );
        fillUnderOverFlow( t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] ,  pfmett1new_    , templateWeight );
      }

      else if( etg_ > 50.0 ){
	//else if( h50 > 0 ){
        //templateWeight = h50;
        iTrigBin = 2;

	fillUnderOverFlow( hphotonPt50_exc , etg_ , templateWeight );
	fillUnderOverFlow( hphotonAll , etg_  , templateWeight );
	fillUnderOverFlow( hnvtxPt50  , nvtx_ , templateWeight );
	fillUnderOverFlow( hnvtxAll   , nvtx_ , templateWeight );

	if( vtxreweight ) templateWeight *= reweightHist[2]->GetBinContent(vtxbin);

        fillUnderOverFlow(   tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_         , templateWeight );
        fillUnderOverFlow(   pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_         , templateWeight );
        fillUnderOverFlow( t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmett1_       , templateWeight );
        fillUnderOverFlow( t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] ,  pfmett1new_    , templateWeight );
      }

      else if( etg_ > 36.0 ){
	//else if( h30 > 0 ){
	//templateWeight = h30;
	iTrigBin = 1;

	fillUnderOverFlow( hphotonPt30_exc , etg_ , templateWeight );
	fillUnderOverFlow( hphotonAll , etg_  , templateWeight );
	fillUnderOverFlow( hnvtxPt30  , nvtx_ , templateWeight );
	fillUnderOverFlow( hnvtxAll   , nvtx_ , templateWeight );

	if( vtxreweight ) templateWeight *= reweightHist[1]->GetBinContent(vtxbin);

	fillUnderOverFlow(   tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_         , templateWeight );
	fillUnderOverFlow(   pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_         , templateWeight );
	fillUnderOverFlow( t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmett1_       , templateWeight );
	fillUnderOverFlow( t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] ,  pfmett1new_    , templateWeight );
      }

      else if( etg_ > 20.0 ){
	//else if( h20 > 0 ){
        //templateWeight = h20;
        iTrigBin = 0;

	fillUnderOverFlow( hphotonPt20_exc , etg_ , templateWeight );
	fillUnderOverFlow( hphotonAll , etg_  , templateWeight );
	fillUnderOverFlow( hnvtxPt20  , nvtx_ , templateWeight );
	fillUnderOverFlow( hnvtxAll   , nvtx_ , templateWeight );

	if( vtxreweight ) templateWeight *= reweightHist[0]->GetBinContent(vtxbin);

        fillUnderOverFlow(   tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  tcmet_         , templateWeight );
        fillUnderOverFlow(   pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmet_         , templateWeight );
        fillUnderOverFlow( t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]    ,  pfmett1_       , templateWeight );
        fillUnderOverFlow( t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] ,  pfmett1new_    , templateWeight );
      }

      else{
	cout << "NO TRIGGERS PASS!!!" << endl;
	exit(0);
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

    for( int iTrigBin = 0 ; iTrigBin < 5 ; ++iTrigBin ){
      for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
        for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
          
          float scale = tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Integral();
          if( scale > 0 )  tcmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
          
          scale = pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Integral();
          if( scale > 0 )  pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );

          scale = t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Integral();
          if( scale > 0 )  t1pfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );

          scale = t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Integral();
          if( scale > 0 )  t1newpfmetTemplate_photon[ iTrigBin ][ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
          
        }
      }
    }


    // make histos rootfile
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();

    char* vtxchar = "";
    if( vtxreweight ) vtxchar = "_vtxreweight";

    char* bvetochar = "";
    if( bveto ) bvetochar = "_bveto";

    char* mjjchar = "";
    if( mjjcut ) mjjchar = "_mjjcut";

    char* pt40char = "";
    //if( pt40 ) pt40char = "_pt40";
    //if( pt40 ) pt40char = "_pt40_2012AB";
    if( pt40 ) pt40char = "_pt40_2012C";

    cout << "Writing templates to " << Form("../photon_output/%s/%s_templates%s%s%s%s.root",iter,sample,vtxchar,bvetochar,mjjchar,pt40char) << endl;
    saveHist(Form("../photon_output/%s/%s_templates%s%s%s%s.root",iter,sample,vtxchar,bvetochar,mjjchar,pt40char));

    //deleteHistos();

    // TFile* fout = TFile::Open(Form("../photon_output/%s/%s_templateHistos.root",iter,sample),"RECREATE");
    // fout->cd();
    // hphotonPt20->Write();
    // hphotonPt30->Write();
    // hphotonPt50->Write();
    // hphotonPt70->Write();
    // hphotonPt90->Write();
    // hphotonAll->Write();
    // fout->Close();
  
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
  hphotonPt90  = new TH1F("hphotonPt90", "",500,0,500);
  hphotonAll   = new TH1F("hphotonAll" , "",500,0,500);

  hphotonPt20_exc  = new TH1F("hphotonPt20_exc", "",500,0,500);
  hphotonPt30_exc  = new TH1F("hphotonPt30_exc", "",500,0,500);
  hphotonPt50_exc  = new TH1F("hphotonPt50_exc", "",500,0,500);
  hphotonPt70_exc  = new TH1F("hphotonPt70_exc", "",500,0,500);
  hphotonPt90_exc  = new TH1F("hphotonPt90_exc", "",500,0,500);

  hnvtxPt20  = new TH1F("hnvtxPt20", "",50,0,50);
  hnvtxPt30  = new TH1F("hnvtxPt30", "",50,0,50);
  hnvtxPt50  = new TH1F("hnvtxPt50", "",50,0,50);
  hnvtxPt70  = new TH1F("hnvtxPt70", "",50,0,50);
  hnvtxPt90  = new TH1F("hnvtxPt90", "",50,0,50);
  hnvtxAll   = new TH1F("hnvtxAll" , "",50,0,50);

  int maxmet = 350;
  if( bveto ) maxmet = 250;

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
  
  char* trigName[5]={(char*)"HLT20",(char*)"HLT30",(char*)"HLT50",(char*)"HLT75",(char*)"HLT90"};

  for( int iTrigBin = 0 ; iTrigBin < 5 ; iTrigBin++ ){
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
        
        tcmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetTemplate_photon_%i_%i_%i",iTrigBin,iJetBin,iSumJetPtBin),
									 Form("%s, %s, %s",trigName[iTrigBin],
									      jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
	
        pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin] = new TH1F(Form("pfmetTemplate_photon_%i_%i_%i",iTrigBin,iJetBin,iSumJetPtBin),
									 Form("%s, %s, %s",trigName[iTrigBin],
									      jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);

        t1pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin] = new TH1F(Form("t1pfmetTemplate_photon_%i_%i_%i",iTrigBin,iJetBin,iSumJetPtBin),
									   Form("%s, %s, %s",trigName[iTrigBin],
										jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
        

        t1newpfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin] = new TH1F(Form("t1newpfmetTemplate_photon_%i_%i_%i",iTrigBin,iJetBin,iSumJetPtBin),
									      Form("%s, %s, %s",trigName[iTrigBin],
										   jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
        
        tcmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->Sumw2();
        pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->Sumw2();
        t1pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->Sumw2();
        t1newpfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->Sumw2();
        
        tcmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
        pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
        t1pfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("type1 pfmet 42X (GeV)");          
        t1newpfmetTemplate_photon[iTrigBin][iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("type1 pfmet 52X(GeV)");          
        
      }
    }
  }

}


void makePhotonTemplates::setBranches (TTree* tree){
  
  tree->SetBranchAddress("run"	               ,        &run_                   );
  tree->SetBranchAddress("tcmet"	       ,        &tcmet_                 );
  tree->SetBranchAddress("pfmet"	       ,        &pfmet_                 );
  tree->SetBranchAddress("pfmett1new"	       ,        &pfmett1new_            );
  tree->SetBranchAddress("pfmett1"	       ,        &pfmett1_               );
  tree->SetBranchAddress("pfmetphi"	       ,        &pfmetphi_              );
  tree->SetBranchAddress("njets"	       ,        &nJets_                 );
  tree->SetBranchAddress("njets40"	       ,        &nJets40_               );
  tree->SetBranchAddress("nbl"	               ,        &nbl_                   );
  tree->SetBranchAddress("nbm"	               ,        &nbm_                   );
  tree->SetBranchAddress("nbt"	               ,        &nbt_                   );
  // tree->SetBranchAddress("nbcsvl"              ,        &nbcsvl_                );
  // tree->SetBranchAddress("nbcsvm"              ,        &nbcsvm_                );
  // tree->SetBranchAddress("nbcsvt"              ,        &nbcsvt_                );
  tree->SetBranchAddress("ht"	               ,        &ht_                    );
  tree->SetBranchAddress("ht40"	               ,        &ht40_                  );
  tree->SetBranchAddress("nvtx"		       ,        &nvtx_                  );
  tree->SetBranchAddress("jetpt"	       ,        &jet_pt_                );
  tree->SetBranchAddress("calojetpt"	       ,	&calojet_pt_		);  
  tree->SetBranchAddress("etag"		       ,        &etag_                  ); 
  tree->SetBranchAddress("phig"		       ,        &phig_                  ); 
  tree->SetBranchAddress("etg"		       ,        &etg_                   ); 
  tree->SetBranchAddress("pfjetid"	       ,        &pfjetid_               );     
  tree->SetBranchAddress("hlt20" 	       ,	&hlt20_		        );
  tree->SetBranchAddress("hlt30"	       ,	&hlt30_		        );
  tree->SetBranchAddress("hlt50"	       ,	&hlt50_		        );
  tree->SetBranchAddress("hlt75"	       ,	&hlt75_		        );
  tree->SetBranchAddress("hlt90"   	       ,	&hlt90_		        );
  tree->SetBranchAddress("hgg22"	       ,	&hgg22_		        );
  tree->SetBranchAddress("hgg36"	       ,	&hgg36_	        	);
  tree->SetBranchAddress("hgg50"       	       ,	&hgg50_		        );
  tree->SetBranchAddress("hgg75"	       ,	&hgg75_		        );
  tree->SetBranchAddress("hgg90"	       ,	&hgg90_		        );
  tree->SetBranchAddress("hoe"	               ,	&hoe_	        	);
  tree->SetBranchAddress("photon_pixelseed"    ,	&photon_pixelseed_	);
  tree->SetBranchAddress("maxleppt"	       ,	&maxleppt_		);
  tree->SetBranchAddress("elveto"	       ,	&elveto_		);
  tree->SetBranchAddress("jetneutralemfrac"    ,        &jetneutralemfrac_      );
  tree->SetBranchAddress("csc"                 ,        &csc_                   );
  tree->SetBranchAddress("hbhe"                ,        &hbhe_                  );
  tree->SetBranchAddress("hcallaser"           ,        &hcallaser_             );
  tree->SetBranchAddress("ecaltp"              ,        &ecaltp_                );
  tree->SetBranchAddress("trkfail"             ,        &trkfail_               );
  tree->SetBranchAddress("eebadsc"             ,        &eebadsc_               );
  tree->SetBranchAddress("hbhenew"             ,        &hbhenew_               );
  //tree->SetBranchAddress("jet1"                ,        &jet1Ptr_               );
  //tree->SetBranchAddress("jet2"                ,        &jet2Ptr_               );
  //tree->SetBranchAddress("jet1"                ,        &jet1_                  );
  //tree->SetBranchAddress("jet2"                ,        &jet2_                  );

}
