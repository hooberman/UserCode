#include "babylooper.h"
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
// #include "CORE/metSelections.cc"
// #include "CORE/trackSelections.cc"
// #include "CORE/eventSelections.cc"

// #include "CORE/electronSelections.cc"
// #include "CORE/electronSelectionsParameters.cc"
// #include "CORE/muonSelections.cc"
// #include "Tools/goodrun.cc"
// #include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

enum metType   { e_tcmet = 0, e_pfmet = 1 };

//------------USER PARAMS----------
bool debug              = false;
metType myMetType       = e_tcmet;
const int nJetBins      = 11;
const int nSumJetPtBins = 23;
//---------------------------------

using namespace std;

int getJetBin( int njets ){
  
  int bin = njets;
  if( bin >= nJetBins ) bin = nJetBins - 1;

  return bin;
}

string jetString( int bin ){

  stringstream s;
  s << "nJets = " << bin ;

  return s.str();

}

int getSumJetPtBin( float x ){

  //bins array defines the sumJetPt binning
  float bins[nSumJetPtBins+1]={0,25,50,75,100,125,150,175,200,250,300,350,400,450,
                               500,600,700,800,900,1000,2000,3000,4000,5000};
  
  if( x < bins[0] )              return 0;
  if( x >= bins[nSumJetPtBins] ) return nSumJetPtBins - 1;

  int ptbin = -1;

  for( int ibin = 0 ; ibin < nSumJetPtBins+1 ; ibin++){
    if( x >= bins[ibin] && x< bins[ibin+1] ){
      ptbin = ibin;
      break;
    }
  }

  if( ptbin == -1 ) 
    cout << "ERROR CANNOT FIND BIN FOR SUMJETPT " << x << endl;

  return ptbin;
}

string sumJetPtString( int bin ){
  
  float bins[nSumJetPtBins+1]={0,25,50,75,100,125,150,175,200,250,300,350,400,450,
                               500,600,700,800,900,1000,2000,3000,4000,5000};

  stringstream s;
  s << bins[bin] << " < sumJetPt < " << bins[bin+1] << " GeV";

  return s.str();
}

void babylooper::ScanChain (TChain* chain, const char* prefix, bool isData, bool calculateTCMET, metAlgo algo, int nEvents){

  int npass = 0;
  algo_ = algo;

  if( algo_ == e_makeTemplate )    cout << "metAlgo makeTemplate" << endl;
  if( algo_ == e_photonSelection ) cout << "metAlgo photonSelection" << endl;
  if( algo_ == e_ZSelection )      cout << "metAlgo ZSelection" << endl;

  TFile *metTemplateFile;
  string metTemplateString = "";
  
  if( algo_ != e_makeTemplate){
    metTemplateString = "_dataTemplate";
    if( myMetType == e_tcmet )
      metTemplateFile = TFile::Open("root/JetMETTau_250nb_babyhistos_tcmet.root");
    else if( myMetType == e_pfmet )
      metTemplateFile = TFile::Open("root/JetMETTau_250nb_babyhistos_pfmet.root");
    else{
      cout << "UNDEFINED METTYPE!!!" << endl;
      exit(0);
    }
  }

  bookHistos();

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
    TTree *tree = (TTree*)f.Get("T1");
    //cms2.Init(tree);
    setBranches(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();
 
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

      weight_ = 1;
      //require at least 2 jets
      if ( nJets_ < 2 )     continue;

      //require leading jet pt > 50
      if( jetmax_pt_ < 50 ) continue;

      float theMet = -1;
      if     ( myMetType == e_tcmet ) theMet = tcmet_;
      else if( myMetType == e_pfmet ) theMet = pfmet_;

      //fill met template-----------------------------------------------------------------------
      
      if( algo_ == e_makeTemplate ) {

        if( HLT_Jet15U_ == 0 ) continue;
        
        int iJetBin      = getJetBin( nJets_ );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        
        if( debug ) {
          cout << "nJets " << nJets_ << " sumJetPt " << sumJetPt_ << endl;
          cout << "iJetBin " << iJetBin << " iSumJetPtBin " << iSumJetPtBin << " met " << theMet << endl;
        }

        fillUnderOverFlow( metTemplate[ iJetBin ][ iSumJetPtBin ] ,  theMet , weight_ );
        
      }

      
      //apply good photon selection-------------------------------------------------------------
  
      if( algo_ == e_photonSelection ) {

        if( HLT_Photon10_L1R_ == 0         &&         
            HLT_Photon15_L1R_ == 0         && 
            HLT_Photon10_Cleaned_L1R_ == 0 && 
            HLT_Photon15_Cleaned_L1R_ == 0  ) continue;

        if ( etag_ > 1 )                                 continue;
        //if ( jet_eta_ > 1 )                              continue;
        if ( etg_ < 20 )                                 continue;
        if ( (1.-r4_) < 0.05 )                           continue;
        if ( hoe_ > 0.1 )                                continue;
        if ( jet_dr_ > 0.5 )                             continue;
        //if ( jet_neu_emfrac_ + jet_chg_emfrac_< 0.95 )   continue; 
        //if ( drel_ < 0.5 )                               continue;
        if ( jet_neu_emfrac_ < 0.95 )                    continue; 
        if ( ( 1 - jet_neu_emfrac_ ) * etg_ > 1 )        continue;
        //if ( sumJetPt_ < 200. )                          continue;
        
        //dphixmet_  = deltaPhi( tcmetphi_ , phig_ );
        //metPar_    = theMet * cos( dphixmet_ );
        //metPerp_   = theMet * sin( dphixmet_ );
    
      }
      
      //cout << "passz " << passz_ << endl;
      //cout << "pdgid " << pdgid_ << endl;
      //cout << "ptll "  << ptll_  << " ptlt " << ptlt_ << endl;
      //cout << "passe " << passe_ttbarV1_ << " passm " << passm_nomttbar_ << endl;

      //apply Z selection-----------------------------------------------------------------------

      if( algo_ == e_ZSelection ) {

        //if( pdgid_ == 13 && ( flagll_ != 5 || flaglt_ != 5 ) )    continue;
        if( passz_ != 1 )                     continue;
        if( pdgid_ != 11 && pdgid_ != 13 )    continue;
        if( ptll_ < 10 || ptlt_ < 10 )        continue;

        //if( pdgid_ == 11 && !( passe_cand01_ ==1 ) )         continue;
        //if( pdgid_ == 13 && !( passm_nom_    ==1 ) )  continue;
        //if( pdgid_ == 13 && !( passm_nomttbarV2_    ==1 ) )  continue;
        if( pdgid_ == 11 && !( passe_ttbarV1_  == 1 ) )      continue;
        if( pdgid_ == 13 && !( passm_nomttbar_ == 1 ) )  continue;

      }
   
      npass++;

      //fill predicted and observed met histos--------------------------------------------------

      if( algo_ != e_makeTemplate){
        
        int iJetBin      = getJetBin( nJets_ );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        TH1F* hmet     = (TH1F*) metTemplateFile->Get(Form("metTemplate_%i_%i",iJetBin,iSumJetPtBin));
        TH1F* hmetPar  = (TH1F*) metTemplateFile->Get(Form("metParTemplate_%i_%i",iJetBin,iSumJetPtBin));
        TH1F* hmetPerp = (TH1F*) metTemplateFile->Get(Form("metPerpTemplate_%i_%i",iJetBin,iSumJetPtBin));
        
        fillUnderOverFlow( metObserved_njets[iJetBin] , theMet , weight_ );
        metPredicted_njets[iJetBin]->Add( hmet );
        
        fillUnderOverFlow( metObserved , theMet , weight_ );
        metPredicted->Add( hmet );

        fillUnderOverFlow( metParObserved , metPar_ , weight_ );
        metParPredicted->Add( hmetPar );

        fillUnderOverFlow( metPerpObserved , metPerp_ , weight_ );
        metPerpPredicted->Add( hmetPerp );

        if( algo_ == e_ZSelection ) {
          if( pdgid_ == 11 ){
            fillUnderOverFlow( metObserved_ee , theMet , weight_ );
            metPredicted_ee->Add( hmet );
          }
          if( pdgid_ == 13 ){
            fillUnderOverFlow( metObserved_mm , theMet , weight_ );
            metPredicted_mm->Add( hmet );
          }
        }

        if( algo_ == e_photonSelection ){
          if( etg_ < 40 ){
            fillUnderOverFlow( metObserved_ptlt40 , theMet , weight_ );
            metPredicted_ptlt40->Add( hmet );
          }
          if( etg_ > 40 && etg_ < 60 ){
            fillUnderOverFlow( metObserved_pt40_60 , theMet , weight_ );
            metPredicted_pt40_60->Add( hmet );
          }
          if( etg_ > 60 ){
            fillUnderOverFlow( metObserved_ptgt60 , theMet , weight_ );
            metPredicted_ptgt60->Add( hmet );
          }
          if( etg_ < 50 ){
            fillUnderOverFlow( metObserved_ptlt50 , theMet , weight_ );
            metPredicted_ptlt50->Add( hmet );
          }else{
            fillUnderOverFlow( metObserved_ptgt50 , theMet , weight_ );
            metPredicted_ptgt50->Add( hmet );
          }
        }
        
        delete hmet;

      }
    } // end loop over events
  } // end loop over files
  
  cout << npass << " events passing selection" << endl;
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  
  //normalize met templates
  if( algo_ == e_makeTemplate ) {
    
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
       
        float scale = metTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        
        if( scale > 0 )
          metTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );

        scale = metParTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        
        if( scale > 0 )
          metParTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );

        scale = metPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        
        if( scale > 0 )
          metPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
    
      }
    }
  }

  // make histos rootfile
  stringstream rootfilename;
  rootfilename << "root/" << prefix << metTemplateString << "_babyhistos.root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(rootfilename.str().c_str());
  deleteHistos();
  
} // end ScanChain

float babylooper::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}


void babylooper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

void babylooper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  metObserved  = new TH1F("metObserved", "Observed MET",500,0,500);
  metPredicted = new TH1F("metPredicted","Predicted MET",500,0,500);

  metObserved_ptlt40     = new TH1F("metObserved_ptlt40", "Observed MET (p_{T}<40 GeV)",500,0,500);
  metPredicted_ptlt40    = new TH1F("metPredicted_ptlt40","Predicted MET  (p_{T}<40 GeV)",500,0,500);

  metObserved_ptlt50     = new TH1F("metObserved_ptlt50", "Observed MET (p_{T}<50 GeV)",500,0,500);
  metPredicted_ptlt50    = new TH1F("metPredicted_ptlt50","Predicted MET  (p_{T}<50 GeV)",500,0,500);

  metObserved_ptgt50     = new TH1F("metObserved_ptgt50", "Observed MET (p_{T}>50 GeV)",500,0,500);
  metPredicted_ptgt50    = new TH1F("metPredicted_ptgt50","Predicted MET  (p_{T}>50 GeV)",500,0,500);

  metObserved_pt40_60    = new TH1F("metObserved_pt40_60", "Observed MET (40<p_{T}<60 GeV)",500,0,500);
  metPredicted_pt40_60   = new TH1F("metPredicted_pt40_60","Predicted MET  (40<p_{T}<60 GeV)",500,0,500);

  metObserved_ptgt60     = new TH1F("metObserved_ptgt60", "Observed MET (p_{T}>60 GeV)",500,0,500);
  metPredicted_ptgt60    = new TH1F("metPredicted_ptgt60","Predicted MET  (p_{T}>60 GeV)",500,0,500);

  metObserved_ee  = new TH1F("metObserved_ee", "Observed MET (ee)",500,0,500);
  metPredicted_ee = new TH1F("metPredicted_ee","Predicted MET (ee)",500,0,500);

  metObserved_mm  = new TH1F("metObserved_mm", "Observed MET (#mu#mu)",500,0,500);
  metPredicted_mm = new TH1F("metPredicted_mm","Predicted MET (#mu#mu)",500,0,500);

  metParObserved  = new TH1F("metParObserved", "Observed MET (Parallel)",1000,-500,500);
  metParPredicted = new TH1F("metParPredicted","Predicted MET (Parallel)",1000,-500,500);

  metPerpObserved  = new TH1F("metPerpObserved", "Observed MET (Perpendicular)",500,0,500);
  metPerpPredicted = new TH1F("metPerpPredicted","Predicted MET (Perpendicular)",500,0,500);
 
  metObserved->Sumw2();
  metPredicted->Sumw2();

  metObserved_ptlt40->Sumw2();
  metPredicted_ptlt40->Sumw2();

  metObserved_pt40_60->Sumw2();
  metPredicted_pt40_60->Sumw2();

  metObserved_ptgt60->Sumw2();
  metPredicted_ptgt60->Sumw2();

  metObserved_ee->Sumw2();
  metPredicted_ee->Sumw2();

  metObserved_mm->Sumw2();
  metPredicted_mm->Sumw2();

  metParObserved->Sumw2();
  metParPredicted->Sumw2();

  metPerpObserved->Sumw2();
  metPerpPredicted->Sumw2();

  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){

    metObserved_njets[iJetBin]  = new TH1F(Form("metObserved_njets%i",iJetBin), Form("Observed MET NJets %i", iJetBin),500,0,500);
    metPredicted_njets[iJetBin] = new TH1F(Form("metPredicted_njets%i",iJetBin),Form("Predicted MET NJets %i",iJetBin),500,0,500);
    
    metObserved_njets[iJetBin] ->Sumw2();
    metPredicted_njets[iJetBin]->Sumw2();
  }
  

  if( algo_ == e_makeTemplate ) {
    
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
        metTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("metTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                          Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);


        metParTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("metParTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                             Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),1000,-500,500);

        metPerpTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("metPerpTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                              Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
       
        metTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
        metParTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
        metPerpTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 

        metTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("met");
        metParTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("met_{#parallel}");
        metPerpTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("met_{#perp}");
      }
    }
  }
  
}


void babylooper::setBranches (TTree* tree){
  tree->SetBranchAddress("run",          &run_          );
  tree->SetBranchAddress("weight",       &weight_       );
  tree->SetBranchAddress("lumi",         &lumi_         );
  tree->SetBranchAddress("event",        &event_        );
  tree->SetBranchAddress("pfmet",        &pfmet_        );
  tree->SetBranchAddress("pfmetphi",     &pfmetphi_     );
  tree->SetBranchAddress("pfsumet",      &pfsumet_      );
  tree->SetBranchAddress("met",          &met_          );
  tree->SetBranchAddress("metphi",       &metphi_       );
  tree->SetBranchAddress("sumet",        &sumet_        );
  tree->SetBranchAddress("mumet",        &mumet_        );
  tree->SetBranchAddress("mumetphi",     &mumetphi_     );
  tree->SetBranchAddress("musumet",      &musumet_      );
  tree->SetBranchAddress("mujesmet",     &mujesmet_     );
  tree->SetBranchAddress("mujesmetphi",  &mujesmetphi_  );
  tree->SetBranchAddress("mujessumet",   &mujessumet_   );
  tree->SetBranchAddress("genmet",       &genmet_       );
  tree->SetBranchAddress("genmetphi",    &genmetphi_    );
  tree->SetBranchAddress("gensumet",     &gensumet_     );
  tree->SetBranchAddress("tcmet",        &tcmet_        );

  tree->SetBranchAddress("dphixmet",     &dphixmet_     );
  tree->SetBranchAddress("metpar",       &metPar_       );
  tree->SetBranchAddress("metperp",      &metPerp_      );


  tree->SetBranchAddress("tcmetphi",     &tcmetphi_     );     
  tree->SetBranchAddress("tcsumet",      &tcsumet_      );     
  tree->SetBranchAddress("njets",        &nJets_        );     
  tree->SetBranchAddress("sumjetpt",     &sumJetPt_     );     
  tree->SetBranchAddress("vecjetpt",     &vecJetPt_     );     

  //photon stuff
  tree->SetBranchAddress("ng",      &nPhotons_          ); 
  tree->SetBranchAddress("etg",     &etg_               ); 
  tree->SetBranchAddress("etag",    &etag_              ); 
  tree->SetBranchAddress("phig",    &phig_              ); 
  tree->SetBranchAddress("hoe",     &hoe_               ); 
  tree->SetBranchAddress("eciso",   &eciso_             ); 
  tree->SetBranchAddress("hciso",   &hciso_             ); 
  tree->SetBranchAddress("tkiso",   &tkiso_             );   
  tree->SetBranchAddress("swiss",   &swiss_             ); 
  tree->SetBranchAddress("seed",    &seed_              ); 
  tree->SetBranchAddress("s4",      &s4_                ); 
  tree->SetBranchAddress("r4",      &r4_                ); 

  //photon-matched jet stuff
  tree->SetBranchAddress("jetdr",                 &jet_dr_                );           
  tree->SetBranchAddress("jetpt",                 &jet_pt_                );   
  tree->SetBranchAddress("jeteta",                &jet_eta_               ); 
  tree->SetBranchAddress("jetenergy",             &jet_energy_            ); 
  tree->SetBranchAddress("jetchargedemfrac",      &jet_chg_emfrac_        ); 
  tree->SetBranchAddress("jetchargedhadfrac",     &jet_chg_hadfrac_       ); 
  tree->SetBranchAddress("jetneutralemfrac",      &jet_neu_emfrac_        ); 
  tree->SetBranchAddress("jetneutralhadfrac",     &jet_neu_hadfrac_       ); 
  //tree->SetBranchAddress("jetncharged",           &jet_nchg_              ); 
  //tree->SetBranchAddress("jetnmuon",              &jet_nmuon_             ); 
  //tree->SetBranchAddress("jetnneutral",           &jet_nneu_              ); 
  tree->SetBranchAddress("jetdphimet",            &jet_dphimet_           ); 
  tree->SetBranchAddress("jetdpt",                &jet_dpt_               ); 
  tree->SetBranchAddress("jetdrgen",              &jet_drgen_             ); 
  tree->SetBranchAddress("drel",                  &drel_                  ); 

  tree->SetBranchAddress("maxjetpt",              &jetmax_pt_             ); 
  tree->SetBranchAddress("maxjetdphimet",         &jetmax_dphimet_        ); 

  //trigger
  tree->SetBranchAddress("HLT_Jet15U",            &HLT_Jet15U_            ); 
  tree->SetBranchAddress("HLT_Jet30U",            &HLT_Jet30U_            ); 

  tree->SetBranchAddress("HLT_Photon10_L1R",      &HLT_Photon10_L1R_      ); 
  tree->SetBranchAddress("HLT_Photon15_L1R",      &HLT_Photon15_L1R_      ); 

  tree->SetBranchAddress("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_  ); 
  tree->SetBranchAddress("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_  ); 
  tree->SetBranchAddress("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_  ); 

  //Z stuff
  tree->SetBranchAddress("passz",             &passz_             );
  tree->SetBranchAddress("pdgid",             &pdgid_             );
  tree->SetBranchAddress("passm_nom",         &passm_nom_         );
  tree->SetBranchAddress("passm_nomttbar",    &passm_nomttbar_    );  
  tree->SetBranchAddress("passm_nomttbarV2",  &passm_nomttbarV2_  );  
  tree->SetBranchAddress("passe_ttbar",       &passe_ttbar_       );
  tree->SetBranchAddress("passe_ttbarV1",     &passe_ttbarV1_     );
  tree->SetBranchAddress("passe_cand01",      &passe_cand01_      );  
  tree->SetBranchAddress("ptll",              &ptll_              );  
  tree->SetBranchAddress("ptlt",              &ptlt_              ); 
  tree->SetBranchAddress("dilmass",           &dilmass_           ); 
  tree->SetBranchAddress("flagll",            &flagll_            );  
  tree->SetBranchAddress("flaglt",            &flaglt_            ); 

}




