#include "makeTemplates.h"
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
#include "CORE/metSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/muonSelections.h"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"
#include "CORE/ttbarSelections.cc"
#include "CORE/jetSelections.cc"
#include "CORE/triggerUtils.cc"


#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

using namespace tas;

//--------------------------------------------------------------------

bool debug                = false;
const int nJetBins        = 11;
const int nSumJetPtBins   = 23;
const int nBosonPtBins    = 4;

float lumi                = 0.0027945;
char* iter                = "V01-01";
char* jsonfilename        = "Cert_TopAug30_Merged_135059-144114_recover_noESDCS_goodruns.txt";

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

int getBosonPtBin( float x ){

  if     ( x < 40.           ) return 0;
  else if( x > 40 && x < 60  ) return 1;
  else if( x > 60 && x < 100 ) return 2;
  else if( x > 100           ) return 3;
  else { 
    cout << "Error could not find boson pt bin" << endl;
    exit(0);
  }
}

//--------------------------------------------------------------------

string bosonPtString( int bin ){

  if     ( bin == 0 ) return "boson p_{T} < 40 GeV";
  else if( bin == 1 ) return "40 < boson p_{T} < 60 GeV";
  else if( bin == 2 ) return "60 < boson p_{T} < 100 GeV";
  else if( bin == 3 ) return "boson p_{T} > 100 GeV";
  else{
    cout << "Error unrecognized boson pt bin " << bin << endl;
    exit(0);
  }

}

//--------------------------------------------------------------------
  
int getJetBin( int njets ){
  
  int bin = njets;
  if( bin >= nJetBins ) bin = nJetBins - 1;
  
  return bin;
}

//--------------------------------------------------------------------

string jetString( int bin ){

  stringstream s;
  s << "nJets = " << bin ;

  return s.str();

}

//--------------------------------------------------------------------

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

//--------------------------------------------------------------------

string sumJetPtString( int bin ){
  
  float bins[nSumJetPtBins+1]={0,25,50,75,100,125,150,175,200,250,300,350,400,450,
                               500,600,700,800,900,1000,2000,3000,4000,5000};

  stringstream s;
  s << bins[bin] << " < sumJetPt < " << bins[bin+1] << " GeV";

  return s.str();
}

//--------------------------------------------------------------------

void makeTemplates::ScanChain (TChain* chain, const char* prefix, bool isData, 
                               bool calculateTCMET, selectionType mySelectionType, int nEvents, float kFactor){

  selection_ = mySelectionType;

  if     ( selection_ == e_QCDSelection    ) cout << "using QCD selection" << endl;
  else if( selection_ == e_photonSelection ) cout << "using photon selection" << endl;
  else if( selection_ == e_ZSelection      ) cout << "using Z selection" << endl;
  else{
    cout << "Unrecognized selectionType" << endl;
    exit(0);
  }

  set_goodrun_file( jsonfilename );
  
  if( isData ){
    ofile_tcmet.open(  Form( "output/%s/%s_tcmetprintout.txt" , iter , prefix  ) );
    ofile_events.open( Form( "output/%s/%s_highmetevents.txt" , iter , prefix  ) );


    ofile_events << "|" << setw(8)  << "run"          << setw(4) 
                 << "|" << setw(6)  << "lumi"         << setw(4) 
                 << "|" << setw(12) << "event"        << setw(4) 
                 << "|" << setw(6)  << "njets"        << setw(4) 
                 << "|" << setw(6)  << "nbtags"       << setw(4) 
                 << "|" << setw(8)  << "tcmet"        << setw(4) 
                 << "|" << setw(8)  << "pfmet"        << setw(4) 
                 << "|" << setw(8)  << "dphi"         << setw(4) << "|" << endl; 
  }


  bookHistos();

  // make a baby ntuple
  //stringstream babyfilename;
  //babyfilename << prefix << "_baby.root";
  MakeBabyNtuple( Form("output/%s/%s_baby.root", iter , prefix ) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  //pass fail counters
  int nSkip_els_conv_dist = 0;
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
 
    for (unsigned int event = 0 ; event < nEvents; ++event){
      
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
      
      if( isData ) {
	DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
	if (is_duplicate(id) )
	  continue;
      }
      
      //skip events with bad els_conv_dist 
      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
      }
      
      if( skipEvent ){
        cout << "SKIPPING EVENT WITH BAD ELS_CONV_DIST" << endl;       
        nSkip_els_conv_dist++;
        continue;
      }
      
      //good run+event selection-----------------------------------------------------------

      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      if( !cleaning_standardAugust2010( isData) )                    continue;

      
      if(debug) cout << "Pass event selection" << endl;

      InitBabyNtuple();

      // event stuff
      run_    = cms2.evt_run();
      lumi_   = cms2.evt_lumiBlock();
      event_  = cms2.evt_event();

      weight_ = 1.;
      pthat_  = -1;
      if( !isData ){
        weight_ = cms2.evt_scale1fb() * kFactor * lumi;
        pthat_  = cms2.genps_pthat();
      }

      float maxphotonpt = -1;
     
      //stitching together PhotonJet_Pt15 and PhotonJet_Pt30
      if( strcmp( prefix , "PhotonJet") == 0 ){
        if( TString(currentFile->GetTitle()).Contains("PhotonJet_Pt15") ){
          if( cms2.genps_pthat() > 30. ) continue;
        }   
        for( unsigned int ig = 0 ; ig < photons_p4().size() ; ++ig ){
          if( photons_p4().at(ig).pt() > maxphotonpt ) 
            maxphotonpt = photons_p4().at(ig).pt();
        }
 
        fillUnderOverFlow( hgenps_pthat  , genps_pthat() , weight_ );
        fillUnderOverFlow( hphotonpt     , maxphotonpt   , weight_ );
      }
      
      

      //access HLT triggers------------------------------------------------------------------

      for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_L1Jet6U" ) == 0 ){
          if( passHLTTrigger("HLT_L1Jet6U") )                  HLT_L1Jet6U_ = 1;
          else                                                 HLT_L1Jet6U_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_L1Jet10U" ) == 0 ){
          if( passHLTTrigger("HLT_L1Jet10U") )                 HLT_L1Jet10U_ = 1;
          else                                                 HLT_L1Jet10U_ = 0;
        }
        
        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Jet15U" ) == 0 ){
          if( passHLTTrigger("HLT_Jet15U") )                   HLT_Jet15U_ = 1;
          else                                                 HLT_Jet15U_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Jet30U" ) == 0 ){
          if( passHLTTrigger("HLT_Jet30U") )                   HLT_Jet30U_ = 1;
          else                                                 HLT_Jet30U_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "L1_SingleEG5" ) == 0 ){
          if( passHLTTrigger("L1_SingleEG5") )                 L1_SingleEG5_ = 1;
          else                                                 L1_SingleEG5_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon10_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon10_L1R") )             HLT_Photon10_L1R_ = 1;
          else                                                 HLT_Photon10_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon15_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon15_L1R") )             HLT_Photon15_L1R_ = 1;
          else                                                 HLT_Photon15_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon10_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon10_Cleaned_L1R") )     HLT_Photon10_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon10_Cleaned_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon15_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon15_Cleaned_L1R") )     HLT_Photon15_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon15_Cleaned_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon20_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon20_Cleaned_L1R") )     HLT_Photon20_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon20_Cleaned_L1R_ = 0;
        }
      }      

            
      // calomet, pfmet, genmet
      met_       = cms2.evt_met();
      metphi_    = cms2.evt_metPhi();
      sumet_     = cms2.evt_sumet();

      pfmet_    = cms2.evt_pfmet();
      pfmetphi_ = cms2.evt_pfmetPhi();
      pfsumet_  = cms2.evt_pfsumet();

      if (!isData){
        genmet_     = cms2.gen_met();
        genmetphi_  = cms2.gen_metPhi();
        gensumet_   = cms2.gen_sumEt();
      }
      
      //tcmet stuff------------------------------------------------------------------- 

      // out-of-the-box  tcmet stuff
      tcmet_    = evt_tcmet();
      tcmetphi_ = evt_tcmetPhi();
      tcsumet_  = evt_tcsumet();

      if( calculateTCMET ){
        
        metStruct tcmetNewStruct = correctedTCMET();
        tcmetNew_     = tcmetNewStruct.met;
        tcmetphiNew_  = tcmetNewStruct.metphi;
        tcsumetNew_   = tcmetNewStruct.sumet;
        
      }else{
        
        tcmetNew_    = -9999;
        tcmetphiNew_ = -9999;
        tcsumetNew_  = -9999;
       
      }

      nGoodVertex_ = 0;
      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
        if(isGoodVertex(v)) nGoodVertex_++;
      }

                  
      //photon stuff------------------------------------------------------------------ 

      int ijetg    = -1; //index of jet matched to photon

      if( selection_ == e_photonSelection ){
        
        nPhotons_           =  0;
        float maxPhotonPt   = -1;
        int   igmax         = -1;
          
        if( photons_p4().size() == 0 ) continue;
          
        //count photons pt > 10 GeV
        for (unsigned int iphoton = 0 ; iphoton < photons_p4().size() ; iphoton++) {
            
          LorentzVector vphoton = photons_p4().at(iphoton);
            
          if(vphoton.pt() > 10){
            nPhotons_++;
            if(vphoton.pt() > maxPhotonPt){
              maxPhotonPt  = vphoton.pt();
              igmax        = iphoton;
            }
          }
        }
                
        if( igmax < 0 ) continue;

        etg_        = photons_p4()[igmax].pt();
        etag_       = photons_p4()[igmax].eta();
        phig_       = photons_p4()[igmax].phi();
        hoe_        = photons_hOverE()[igmax];
        swiss_      = photons_swissSeed()[igmax];
        int scind   = photons_scindex()[igmax] ;
        seed_       = scs_eSeed()[scind] ;
        s4_         = swiss_ - seed_ ;
        r4_         = 1 - s4_ / seed_ ;

        //eciso_ = photons_ecalIso()[igmax];
        //hciso_ = photons_hcalIso()[igmax];
        //tkiso_ = photons_tkIsoSolid()[igmax];
                
        if( isData ){
          photon_pixelseed_        = photons_haspixelSeed()[igmax] ? 1 : 0;
          photon_e15_              = photons_e1x5()[igmax];      
          photon_e25max_           = photons_e2x5Max()[igmax];      
          photon_e33_              = photons_e3x3()[igmax];           
          photon_e55_              = photons_e5x5()[igmax];           
          photon_ecalIso03_        = photons_ecalIso03()[igmax];      
          photon_ecalIso04_        = photons_ecalIso04()[igmax];      
          photon_hcalIso03_        = photons_hcalIso03()[igmax];      
          photon_hcalIso04_        = photons_hcalIso04()[igmax];      
          photon_ntkIsoHollow03_   = photons_ntkIsoHollow03()[igmax];
          photon_ntkIsoHollow04_   = photons_ntkIsoHollow04()[igmax];
          photon_ntkIsoSolid03_    = photons_ntkIsoSolid03()[igmax]; 
          photon_ntkIsoSolid04_    = photons_ntkIsoSolid04()[igmax]; 
          photon_sigmaEtaEta_      = photons_sigmaEtaEta()[igmax];    
          photon_sigmaIEtaIEta_    = photons_sigmaIEtaIEta()[igmax];
          photon_tkisoHollow03_    = photons_tkIsoHollow03()[igmax];
          photon_tkisoHollow04_    = photons_tkIsoHollow04()[igmax]; 
          photon_tkisoSolid03_     = photons_tkIsoSolid03()[igmax];   
          photon_tkisoSolid04_     = photons_tkIsoSolid04()[igmax];  
        }
      
        jet_dr_      = 10000;
        ijetg    = -1;

        //find jet matched to to photon
        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vg   = photons_p4().at(igmax);
          
          if( vjet.pt()  < 10  )       continue;
          if( vjet.eta() > 2.5 )       continue;
          if( !passesPFJetID(ijet) )   continue;

          float dr = dRbetweenVectors(vjet, vg);
         
          if( dr < jet_dr_ ){
            jet_dr_ = dr;
            ijetg   = ijet;
          }
        }
        
        if( ijetg < 0 ) continue;
        
        jet_pt_             = pfjets_p4().at(ijetg).pt();
        jet_energy_         = pfjets_p4().at(ijetg).energy();
        jet_eta_            = pfjets_p4().at(ijetg).eta();
        jet_chg_emfrac_     = pfjets_chargedEmE().at(ijetg)     / jet_energy_;
        jet_chg_hadfrac_    = pfjets_chargedHadronE().at(ijetg) / jet_energy_;
        jet_neu_emfrac_     = pfjets_neutralEmE().at(ijetg)     / jet_energy_;
        jet_neu_hadfrac_    = pfjets_neutralHadronE().at(ijetg) / jet_energy_;
        jet_nchg_           = pfjets_chargedMultiplicity().at(ijetg);
        jet_nmuon_          = pfjets_muonMultiplicity().at(ijetg);
        jet_nneu_           = pfjets_neutralMultiplicity().at(ijetg);
        jet_dphimet_        = deltaPhi( pfjets_p4().at(ijetg).phi() , tcmetphi_);

      }
              
      //jet stuff--------------------------------------------------------------------- 

                    
      nJets_        = 0;
      sumJetPt_     = 0.;
      nJets40_      = 0;
      sumJetPt10_   = 0.;
      nbtags_       = 0;

      LorentzVector jetSystem(0.,0.,0.,0.);        
      float maxcosdphi  = -99;
      //int   imaxcosdphi = -1;
      int   imaxjet     = -1;
      float maxpt       = -1;
                      
      //loop over pfjets pt > 30 GeV |eta| < 2.5
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
        
        if( !passesPFJetID(ijet) )               continue;

        //skip jet matched to photon
        if( selection_ == e_photonSelection && (int)ijet == ijetg ) continue;
 
        LorentzVector vjet      = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
        if( fabs( vjet.eta() ) > 2.5 )           continue;

        if ( vjet.pt() > 10. ){
          sumJetPt10_ += vjet.pt();
        }
        if ( vjet.pt() > 15. ){
          sumJetPt_ += vjet.pt();
          jetSystem += vjet;
        }
              
        if( vjet.pt() < 30. )                    continue;

        //find max jet pt
        if( pfjets_p4().at(ijet).pt() > maxpt ){
          maxpt   = pfjets_p4().at(ijet).pt();
          imaxjet = ijet;
        }

        //find jet (anti-)aligned with tcmet
        if( fabs( cos( tcmetphi_ - vjet.phi() ) ) > maxcosdphi ){
          maxcosdphi  = fabs( cos( tcmetphi_ - vjet.phi() ) );
          //imaxcosphi  = ijet;
          dphijetmet_ = fabs( tcmetphi_ - vjet.phi() );
          if( dphijetmet_ > TMath::Pi() ) dphijetmet_ = TMath::TwoPi() - dphijetmet_;
        }
        
        //find closest calojet to use btagging info
        float dRmin    = 100;
        int   iCaloJet = -1;
          
        for( unsigned int iC = 0 ; iC < jets_p4().size() ; iC++ ){
            
          LorentzVector vcalojet = jets_p4().at(iC);
          if( vcalojet.pt() * jets_cor().at(iC) < 10 ) continue;
            
          float dR = dRbetweenVectors(vjet, vcalojet);
          if( dR < dRmin ){
            dRmin = dR;
            iCaloJet = iC;
          }
        }
                  
        if( iCaloJet > -1 ){
          if( jets_simpleSecondaryVertexHighEffBJetTag().at(iCaloJet) > 1.74 ) ++nbtags_;
          //if( jets_trackCountingHighEffBJetTag().at(iCaloJet) > 1.7 ) ++nbtags_;
        }

        if ( vjet.pt() > 30. ) nJets_++;
        if ( vjet.pt() > 40. ) nJets40_++;
          
      }
                            
      jetmax_pt_ = -1;

      if( imaxjet > -1 ){
        jetmax_pt_       = pfjets_p4().at(imaxjet).pt();
        jetmax_dphimet_  = deltaPhi( pfjets_p4().at(imaxjet).phi() , tcmetphi_);
      }

      vecJetPt_ = jetSystem.pt();
      
      //fill histos and ntuple----------------------------------------------------------- 
              
      FillBabyNtuple();
      
      //fillHistos( htcmet            , tcmet_           , weight_ , leptype_ , nJets_ );
      //fillHistos( htcmetNew         , tcmetNew_        , weight_ , leptype_ , nJets_ );
      //fillHistos( hpfmet            , pfmet_           , weight_ , leptype_ , nJets_  );
              
      if( isData && ( tcmet_ > 30 || pfmet_ > 30 ) ){

        metStruct dummyStruct = correctedTCMET( true, ofile_tcmet );

        ofile_events << "|" << setw(8)  << evt_run()                   << setw(4) 
                     << "|" << setw(6)  << evt_lumiBlock()             << setw(4) 
                     << "|" << setw(12) << evt_event()                 << setw(4) 
                     << "|" << setw(6)  << nJets_                      << setw(4) 
                     << "|" << setw(6)  << nbtags_                     << setw(4) 
                     << "|" << setw(8)  << fround(tcmet_,1)            << setw(4) 
                     << "|" << setw(8)  << fround(pfmet_,1)            << setw(4) 
                     << "|" << setw(8)  << fround(dphijetmet_,2)       << setw(4) << "|" << endl; 
       
      }
      
      //selection for met templates----------------------------------------------------- 
      
      //QCD selection
      if( selection_ == e_QCDSelection ) {
        
        if( !(HLT_Jet15U_ == 1) )  continue;
        if( jetmax_pt_ < 50  )     continue;
        if( nJets_ < 2 )           continue;
        
      }

      //photon+jets selection
      else if( selection_ == e_photonSelection ){
        
        if( !(HLT_Photon10_L1R_ == 1         || HLT_Photon15_L1R_ == 1 || 
              HLT_Photon10_Cleaned_L1R_ == 1 || HLT_Photon15_Cleaned_L1R_ == 1)  ) continue;
        if( jetmax_pt_ < 30  )  continue;
        if( nJets_ < 1 )        continue;

        //Ben's photon-object selection
        if ( etag_ > 1 )                                 continue;
        if ( etg_ < 20 )                                 continue;
        if ( (1.-r4_) < 0.05 )                           continue;
        if ( hoe_ > 0.1 )                                continue;
        if ( jet_neu_emfrac_ < 0.95 )                    continue; 
          
        /*
        //Warren's photon selection
        if ( etag_ > 1 )                                 continue;
        if ( etg_ < 20 )                                 continue;
        if ( (1.-r4_) < 0.05 )                           continue;
        if ( hoe_ > 0.05 )                               continue;
        if ( photon_hcalIso03_ > 2.4 && photon_hcalIso03_ / etg_ > 0.05)        continue;              
        if ( photon_ecalIso03_ > 1.7 && photon_ecalIso03_ / etg_ > 0.05)        continue;       
        if ( photon_ntkIsoSolid03_ > 2 || photon_tkisoSolid03_ / etg_ > 0.1 )   continue;
        */
        
        int iJ = nJets_;
        if( iJ > 4 ) iJ = 4;
        fillUnderOverFlow( hetg[iJ]  , etg_  , weight_ );
        fillUnderOverFlow( hetg[0]   , etg_  , weight_ );
      }

      else{
        cout << "UNRECOGNIZED SELECTION ENUM" << selection_ << endl;
        exit(0);
      }

      int iJetBin      = getJetBin( nJets_ );
      int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
      int iBosonPtBin  = getBosonPtBin( etg_ );
  
      //fill templates binned by njets, sumjetpt, boson pt        
      fillUnderOverFlow( tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  tcmet_    , weight_ );
      fillUnderOverFlow( pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  pfmet_    , weight_ );
      fillUnderOverFlow( tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] ,  tcmetNew_ , weight_ );
      
      //fill templates binned by njets, sumjetpt
      fillUnderOverFlow( tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , weight_ );
      fillUnderOverFlow( pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , weight_ );
      fillUnderOverFlow( tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] ,  tcmetNew_ , weight_ );
      
      /*
        float dphi_tcmet    = deltaPhi( tcmetphi_ , phig_ );
        float tcmetPar      = tcmet_ * cos( dphi_tcmet );
        float tcmetPerp     = tcmet_ * sin( dphi_tcmet );

        float dphi_tcmetNew    = deltaPhi( tcmetphiNew_, phig_ );
        float tcmetNewPar      = tcmetNew_ * cos( dphi_tcmetNew );
        float tcmetNewPerp     = tcmetNew_ * sin( dphi_tcmetNew );

        float dphi_pfmet    = deltaPhi( pfmetphi_ , phig_ );
        float pfmetPar      = pfmet_ * cos( dphi_pfmet );
        float pfmetPerp     = pfmet_ * sin( dphi_pfmet );

        fillUnderOverFlow( tcmetParTemplate[ iJetBin ][ iSumJetPtBin ]  , tcmetPar   , weight_ );
        fillUnderOverFlow( tcmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] , tcmetPerp  , weight_ );
     
        fillUnderOverFlow( tcmetNewParTemplate[ iJetBin ][ iSumJetPtBin ]  , tcmetNewPar   , weight_ );
        fillUnderOverFlow( tcmetNewPerpTemplate[ iJetBin ][ iSumJetPtBin ] , tcmetNewPerp  , weight_ );

        fillUnderOverFlow( pfmetParTemplate[ iJetBin ][ iSumJetPtBin ]  , pfmetPar   , weight_ );
        fillUnderOverFlow( pfmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] , pfmetPerp  , weight_ );
      */

    } // end loop over events
  } // end loop over files

  if( nSkip_els_conv_dist > 0 ){
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist branch" << endl;
  }

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  CloseBabyNtuple();

  //normalize met templates
  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      for( int iBosonPtBin = 0 ; iBosonPtBin < nBosonPtBins ; iBosonPtBin++ ){
        
        float scale = tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
        if( scale > 0 )  tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
        
        scale = tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
        if( scale > 0 )  tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
        
        scale = pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
        if( scale > 0 )  pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
        
      }
    }
  }
  
  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      
      float scale = tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
      if( scale > 0 )  tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );

      scale = tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
      if( scale > 0 )  tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      
      scale = pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
      if( scale > 0 )  pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
     

      /*
        scale = tcmetParTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetParTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      
        scale = tcmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      
        scale = tcmetNewParTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetNewParTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      
        scale = tcmetNewPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetNewPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
     
        scale = pfmetParTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  pfmetParTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      
        scale = pfmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  pfmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
      */

    }
  }
  
  // make histos rootfile
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form("output/%s/%s_templates.root", iter , prefix ) );
  deleteHistos();
  
} // end ScanChain

//--------------------------------------------------------------------

float makeTemplates::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//--------------------------------------------------------------------

void makeTemplates::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

//--------------------------------------------------------------------

void makeTemplates::InitBabyNtuple (){

  //triggers
  HLT_L1Jet6U_          = -1;
  HLT_L1Jet10U_         = -1;
  HLT_Jet15U_           = -1;
  HLT_Jet30U_           = -1;
  L1_SingleEG5_         = -1;
  HLT_Photon10_L1R_     = -1;
  HLT_Photon15_L1R_     = -1;
  HLT_Photon10_Cleaned_L1R_     = -1;
  HLT_Photon15_Cleaned_L1R_     = -1;
  HLT_Photon20_Cleaned_L1R_     = -1;

  // event stuff
  run_          = -999999;
  lumi_         = -999999;
  event_        = -999999;
  weight_       = -999999.;
  pthat_        = -999999.;
  nGoodVertex_  = -999999;
  leptype_      = -999999;

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
  dphixmet_     = -999999.;
  metPar_       = -999999.;
  metPerp_      = -999999.;

  tcmet_        = -999999.;
  tcmetphi_     = -999999.;
  tcsumet_      = -999999.;

  tcmetNew_     = -999999.;
  tcsumetNew_   = -999999.;
  tcmetphiNew_  = -999999.;

  nJets_        = -999999;
  sumJetPt_     = -999999;
  vecJetPt_     = -999999;
  nJets40_      = -999999;
  sumJetPt10_   = -999999;

  nbtags_       = -999999;
  dphijetmet_   = -999999;

  //leading jet stuff
  jetmax_pt_        = -999999;
  jetmax_dphimet_   = -999999;


  //photon stuff
  nPhotons_ = -999999;
  etg_      = -999999.;
  etag_     = -999999.;
  phig_     = -999999.;
  hoe_      = -999999.;
  eciso_    = -999999.;
  hciso_    = -999999.;
  tkiso_    = -999999.;
  swiss_    = -999999.;
  seed_     = -999999.;
  s4_       = -999999.;
  r4_       = -999999.;

  //more photon stuff
  photon_pixelseed_        = -999999;
  photon_e15_              = -999999.;
  photon_e25max_           = -999999.;
  photon_e33_              = -999999.;
  photon_e55_              = -999999.;
  photon_ecalIso03_        = -999999.;
  photon_ecalIso04_        = -999999.;
  photon_hcalIso03_        = -999999.;
  photon_hcalIso04_        = -999999.;
  photon_ntkIsoHollow03_   = -999999.;
  photon_ntkIsoHollow04_   = -999999.;
  photon_ntkIsoSolid03_    = -999999.;
  photon_ntkIsoSolid04_    = -999999.;
  photon_sigmaEtaEta_      = -999999.;
  photon_sigmaIEtaIEta_    = -999999.;
  photon_tkisoHollow03_    = -999999.;
  photon_tkisoHollow04_    = -999999.;
  photon_tkisoSolid03_     = -999999.;
  photon_tkisoSolid04_     = -999999.;
                                  
  //photon-matched jet stuff
  jet_eta_          = -999999.;  
  jet_energy_       = -999999.;  
  jet_pt_           = -999999.;  
  jet_chg_emfrac_   = -999999.;  
  jet_chg_hadfrac_  = -999999.;  
  jet_neu_emfrac_   = -999999.;  
  jet_neu_hadfrac_  = -999999.;  
  jet_nchg_         = -999999;      
  jet_nmuon_        = -999999;  
  jet_nneu_         = -999999;  
  jet_dphimet_      = -999999.;  
  jet_dpt_          = -999999.;  

}

//--------------------------------------------------------------------

void makeTemplates::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hgenps_pthat = new TH1F("hgenps_pthat","",100,0,100);
  hphotonpt    = new TH1F("hphotonpt","",100,0,100);

  hgenps_pthat->GetXaxis()->SetTitle("gen p_{T}(hat) (GeV)");
  hphotonpt->GetXaxis()->SetTitle("max photon p_{T} (GeV)");

  char* pttitle[5]={"all jets","1 jet","2 jet","3 jet","#geq 4 jet"};

  for( int iJ = 0 ; iJ < 5 ; iJ++ ){
    hetg[iJ] = new TH1F(Form("hetg_%i",iJ),pttitle[iJ],200,0,200);
    hetg[iJ]->GetXaxis()->SetTitle("photon p_{T} (GeV)");
  }


  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      for( int iBosonPtBin = 0 ; iBosonPtBin < nBosonPtBins ; iBosonPtBin++ ){

        tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("tcmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                     Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                          bosonPtString(iBosonPtBin).c_str()),500,0,500);

        tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("tcmetNewTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                        Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                             bosonPtString(iBosonPtBin).c_str()),500,0,500);

        pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("pfmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                     Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                          bosonPtString(iBosonPtBin).c_str()),500,0,500);
        
        tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();
        tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();
        pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();

        tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
        tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("tcmetNew (GeV)");
        pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          

      }
    }
  }


  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      

      tcmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                               Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
        
      tcmetNewTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetNewTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                  Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
        
      pfmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("pfmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                               Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
        
      tcmetTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
      tcmetNewTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
      pfmetTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
        
      tcmetTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
      tcmetNewTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmetNew (GeV)");
      pfmetTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
        
    }
  }


  /*
  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
      //tcmet
      tcmetTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("tcmetTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                          Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
  
      tcmetParTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("tcmetParTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                             Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),1000,-500,500);
  
      tcmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("tcmetPerpTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                              Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
  
      tcmetTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
      tcmetParTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
      tcmetPerpTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
  
      tcmetTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("tcmet");
      tcmetParTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("tcmet_{#parallel}");
      tcmetPerpTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("tcmet_{#perp}");

      //tcmetNew
      tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("tcmetNewTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                             Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
  
      tcmetNewParTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("tcmetNewParTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                                Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),1000,-500,500);
  
      tcmetNewPerpTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("tcmetNewPerpTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                                 Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
  
      tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
      tcmetNewParTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
      tcmetNewPerpTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
  
      tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("tcmetNew");
      tcmetNewParTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("tcmetNew_{#parallel}");
      tcmetNewPerpTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("tcmetNew_{#perp}");
  
      //pfmet
      pfmetTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("pfmetTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                          Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
  
      pfmetParTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("pfmetParTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                             Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),1000,-500,500);
  
      pfmetPerpTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("pfmetPerpTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                              Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
  
      pfmetTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
      pfmetParTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
      pfmetPerpTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
  
      pfmetTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("pfmet");
      pfmetParTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("pfmet_{#parallel}");
      pfmetPerpTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("pfmet_{#perp}");
    }
  }
  */
  
  char* leptype[4]   = {"ee", "mm", "em", "all"};
  char* jetbin[4]    = {"0j", "1j", "geq2j", "allj"};

  char* leptype_title[4]   = {"ee", "#mu#mu", "e#mu", "all leptons"};
  char* jetbin_title[4]    = {"0 jets", "1 jet", "#geq 2 jets", "all jets"};

  for (int i = 0; i < 4; i++) {
   
    hdilMass[i] = new TH1F(Form("hdilMass_%s",leptype[i]),  leptype_title[i],   150,0,300);
    hdilMass[i]->GetXaxis()->SetTitle("M(ll) (GeV)");
 
    for (int j = 0; j < 4; j++) {

      char* suffix       = Form("%s_%s",leptype[i],jetbin[j]);
      char* suffix_title = Form("%s %s",leptype_title[i],jetbin_title[j]);
    
      htcmet[i][j]    = new TH1F(Form("htcmet_%s",suffix),    suffix_title, 100,0,100);
      htcmetNew[i][j] = new TH1F(Form("htcmetNew_%s",suffix), suffix_title, 100,0,100);
      hpfmet[i][j]    = new TH1F(Form("hpfmet_%s",suffix),    suffix_title, 100,0,100);
      htcmet[i][j]->GetXaxis()->SetTitle("tcmet (GeV)");
      htcmetNew[i][j]->GetXaxis()->SetTitle("tcmetNew (GeV)");
      hpfmet[i][j]->GetXaxis()->SetTitle("pfmet (GeV)");
    }
  }
}
 
//--------------------------------------------------------------------

void makeTemplates::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("T1", "A Baby Ntuple");

  //event stuff
  babyTree_->Branch("run",          &run_,          "run/I"  );
  babyTree_->Branch("lumi",         &lumi_,         "lumi/I" );
  babyTree_->Branch("event",        &event_,        "event/I");
  babyTree_->Branch("nvtx",         &nGoodVertex_,  "nvtx/I");
  babyTree_->Branch("weight",       &weight_,       "weight/F");
  babyTree_->Branch("pthat",        &pthat_,        "pthat/F");

  //met stuff
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
  babyTree_->Branch("dphixmet",     &dphixmet_,      "dphixmet/F"    );
  babyTree_->Branch("metpar",       &metPar_,        "metpar/F"      );
  babyTree_->Branch("metperp",      &metPerp_,       "metperp/F"     );
  babyTree_->Branch("tcmet",        &tcmet_,        "tcmet/F"      );
  babyTree_->Branch("tcmetphi",     &tcmetphi_,     "tcmetphi/F"   );
  babyTree_->Branch("tcsumet",      &tcsumet_,      "tcsumet/F"    );
  babyTree_->Branch("tcmetNew",     &tcmetNew_,     "tcmetNew/F"      );
  babyTree_->Branch("tcmetphiNew",  &tcmetphiNew_,  "tcmetphiNew/F"   );
  babyTree_->Branch("tcsumetNew",   &tcsumetNew_,   "tcsumetNew/F"    );

  //jet stuff
  babyTree_->Branch("njets",          &nJets_,          "njets/I"       );
  babyTree_->Branch("njets40",        &nJets40_,        "njets40/I"     );
  babyTree_->Branch("sumjetpt",       &sumJetPt_,       "sumjetpt/F"    );
  babyTree_->Branch("sumjetpt10",     &sumJetPt10_,     "sumjetpt10/F"    );
  babyTree_->Branch("vecjetpt",       &vecJetPt_,       "vecjetpt/F"    );
  babyTree_->Branch("nbtags",         &nbtags_,         "nbtags/I");
  babyTree_->Branch("ndphijetmet",    &dphijetmet_,     "dphijetmet/F");
  babyTree_->Branch("maxjetpt",       &jetmax_pt_,      "maxjetpt/F");
  babyTree_->Branch("maxjetdphimet",  &jetmax_dphimet_, "maxjetdphimet/F");
                                   
  //trigger stuff
  babyTree_->Branch("HLT_Jet15U",                    &HLT_Jet15U_,                   "HLT_Jet15U/I");
  babyTree_->Branch("HLT_Jet30U",                    &HLT_Jet30U_,                   "HLT_Jet30U/I");
  babyTree_->Branch("HLT_Photon10_L1R",              &HLT_Photon10_L1R_,             "HLT_Photon10_L1R/I");
  babyTree_->Branch("HLT_Photon15_L1R",              &HLT_Photon15_L1R_,             "HLT_Photon15_L1R/I");
  babyTree_->Branch("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_,     "HLT_Photon10_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_,     "HLT_Photon15_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_,     "HLT_Photon20_Cleaned_L1R/I");  


  //photon stuff
  babyTree_->Branch("ng",      &nPhotons_, "ng/I"); 
  babyTree_->Branch("etg",     &etg_,      "etg/F");	   
  babyTree_->Branch("etag",    &etag_,     "etag/F");	   
  babyTree_->Branch("phig",    &phig_,     "phig/F");
  babyTree_->Branch("hoe",     &hoe_,      "hoe/F");	   
  babyTree_->Branch("eciso",   &eciso_,    "eciso/F");	   
  babyTree_->Branch("hciso",   &hciso_,    "hciso/F");	   
  babyTree_->Branch("tkiso",   &tkiso_,    "tkiso/F");
  babyTree_->Branch("swiss",   &swiss_,    "swiss/F");
  babyTree_->Branch("seed",    &seed_,     "seed/F");
  babyTree_->Branch("s4",      &s4_,       "s4/F");
  babyTree_->Branch("r4",      &r4_,       "r4/F");

  //more photon stuff
  babyTree_->Branch("photon_pixelseed",         &photon_pixelseed_,         "photon_pixelseed/I");         
  babyTree_->Branch("photon_e15",               &photon_e15_,               "photon_e15/F");                
  babyTree_->Branch("photon_e25max",            &photon_e25max_,            "photon_e25max/F");             
  babyTree_->Branch("photon_e33",               &photon_e33_,               "photon_e33/F");                
  babyTree_->Branch("photon_e55",               &photon_e55_,               "photon_e55/F");                
  babyTree_->Branch("photon_ecalIso03",         &photon_ecalIso03_,         "photon_ecalIso03/F");          
  babyTree_->Branch("photon_ecalIso04",         &photon_ecalIso04_,         "photon_ecalIso04/F");          
  babyTree_->Branch("photon_hcalIso03",         &photon_hcalIso03_,         "photon_hcalIso03/F");          
  babyTree_->Branch("photon_hcalIso04",         &photon_hcalIso04_,         "photon_hcalIso04/F");          
  babyTree_->Branch("photon_ntkIsoHollow03",    &photon_ntkIsoHollow03_,    "photon_ntkIsoHollow03/F");     
  babyTree_->Branch("photon_ntkIsoHollow04",    &photon_ntkIsoHollow04_,    "photon_ntkIsoHollow04/F");     
  babyTree_->Branch("photon_ntkIsoSolid03",     &photon_ntkIsoSolid03_,     "photon_ntkIsoSolid03/F");      
  babyTree_->Branch("photon_ntkIsoSolid04",     &photon_ntkIsoSolid04_,     "photon_ntkIsoSolid04/F");      
  babyTree_->Branch("photon_sigmaEtaEta",       &photon_sigmaEtaEta_,       "photon_sigmaEtaEta/F");        
  babyTree_->Branch("photon_sigmaIEtaIEta",     &photon_sigmaIEtaIEta_,     "photon_sigmaIEtaIEta/F");      
  babyTree_->Branch("photon_tkisoHollow03",     &photon_tkisoHollow03_,     "photon_tkisoHollow03/F");      
  babyTree_->Branch("photon_tkisoHollow04",     &photon_tkisoHollow04_,     "photon_tkisoHollow04/F");      
  babyTree_->Branch("photon_tkisoSolid03",      &photon_tkisoSolid03_,      "photon_tkisoSolid03/F");      
  babyTree_->Branch("photon_tkisoSolid04",      &photon_tkisoSolid04_,      "photon_tkisoSolid04/F");           
                                                                            
  //photon-matched jet stuff
  babyTree_->Branch("jetdr",                 &jet_dr_,               "jetdr/F");
  babyTree_->Branch("jetpt",                 &jet_pt_,               "jetpt/F");
  babyTree_->Branch("jeteta",                &jet_eta_,              "jeteta/F");
  babyTree_->Branch("jetenergy",             &jet_energy_,           "jetenergy/F");
  babyTree_->Branch("jetchargedemfrac",      &jet_chg_emfrac_,       "jetchargedemfrac/F");
  babyTree_->Branch("jetchargedhadfrac",     &jet_chg_hadfrac_,      "jetchargedhadfrac/F");
  babyTree_->Branch("jetneutralemfrac",      &jet_neu_emfrac_,       "jetneutralemfrac/F");
  babyTree_->Branch("jetneutralhadfrac",     &jet_neu_hadfrac_,      "jetneutralhadfrac/F");
  babyTree_->Branch("jetncharged",           &jet_nchg_,             "jetncharged/I");
  babyTree_->Branch("jetnmuon",              &jet_nmuon_,            "jetnmuon/I");
  babyTree_->Branch("jetnneutral",           &jet_nneu_,             "jetnneutral/I");
  babyTree_->Branch("jetdphimet",            &jet_dphimet_,          "jetdphimet/F");
  babyTree_->Branch("jetdpt",                &jet_dpt_,              "jetdpt/F");

}

//--------------------------------------------------------------------

void makeTemplates::FillBabyNtuple ()
{
  babyTree_->Fill();
}

//--------------------------------------------------------------------

void makeTemplates::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}

//--------------------------------------------------------------------

void makeTemplates::fillHistos(TH1F *h1[4],float value, float weight, int myType)
{

  fillUnderOverFlow(h1[myType], value, weight);      
  fillUnderOverFlow(h1[3],      value, weight);      
}

//--------------------------------------------------------------------

void makeTemplates::fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{

  if( nJetsIdx > 2 ) nJetsIdx = 2;
  
  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      

}

//--------------------------------------------------------------------
