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
#include "CORE/metSelections.cc"
#include "CORE/trackSelections.cc"
#include "CORE/eventSelections.cc"

#include "CORE/electronSelections.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/muonSelections.cc"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

bool debug = false;

const int nJetBins      = 11;
const int nSumJetPtBins = 23;

using namespace tas;

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

void looper::ScanChain (TChain* chain, const char* prefix, bool isData, bool calculateTCMET, metAlgo algo, int nEvents){

  algo_ = algo;

  if( algo_ == e_makeTemplate )    cout << "metAlgo makeTemplate" << endl;
  if( algo_ == e_photonSelection ) cout << "metAlgo photonSelection" << endl;
  if( algo_ == e_ZSelection )      cout << "metAlgo ZSelection" << endl;

  TFile *metTemplateFile;
  string metTemplateString = "";
  
  DorkyEventIdentifier dei;

  if( algo_ != e_makeTemplate){
    
    metTemplateString = "_dataTemplate";
    metTemplateFile = TFile::Open("root/JetMETTau_250nb.root");
    //metTemplateFile = TFile::Open("root/JetMETTau_metTemplate_maxjetpt30.root");
  
    //metTemplateString = "_mcTemplate";
    //metTemplateFile = TFile::Open("root/QCD_Pt15_metTemplate.root");
  }


  //set_goodrun_file("goodruns_official_0526.txt");
  //set_goodrun_file("/home/users/kalavase/UserCode/Kalavase/TTbarThesis/840_nb/Cert_TopAug13_Merged_135059-142664_goodruns.txt");
  set_goodrun_file("Cert_TopAug13_Merged_135059-142664_goodruns.txt");

  bookHistos();

  // make a baby ntuple
  stringstream babyfilename;
  babyfilename << prefix << "_baby.root";
  MakeBabyNtuple( Form("root/%s_baby.root", prefix ) );

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

      //duplicate event veto CHECK THIS
      bool myduplicate = dei.is_duplicate(DorkyEvent());
      if(myduplicate) continue;
      
      //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Kalavase/TTbarThesis/ScanChain.C?revision=1.33&view=markup
     
      //basic event selection----------------------------------------------------------------
      
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
      
      if ( !cleaning_standard( isData ) ) continue;

      if(debug) cout << "Pass event selection" << endl;
   
      InitBabyNtuple();

      // event stuff
      run_    = cms2.evt_run();
      lumi_   = cms2.evt_lumiBlock();
      event_  = cms2.evt_event();

     
      //access HLT triggers------------------------------------------------------------------

      for( unsigned int i = 0 ; i < hlt_trigNames().size() ; i++ ){

        if( strcmp(hlt_trigNames().at(i) , "HLT_L1Jet6U" ) == 0 ){
          if( passHLTTrigger("HLT_L1Jet6U") )          HLT_L1Jet6U_ = 1;
          else                                         HLT_L1Jet6U_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_L1Jet10U" ) == 0 ){
          if( passHLTTrigger("HLT_L1Jet10U") )         HLT_L1Jet10U_ = 1;
          else                                         HLT_L1Jet10U_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Jet15U" ) == 0 ){
          if( passHLTTrigger("HLT_Jet15U") )           HLT_Jet15U_ = 1;
          else                                         HLT_Jet15U_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Jet30U" ) == 0 ){
          if( passHLTTrigger("HLT_Jet30U") )           HLT_Jet30U_ = 1;
          else                                         HLT_Jet30U_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "L1_SingleEG5" ) == 0 ){
          if( passHLTTrigger("L1_SingleEG5") )         L1_SingleEG5_ = 1;
          else                                         L1_SingleEG5_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Photon10_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon10_L1R") )     HLT_Photon10_L1R_ = 1;
          else                                         HLT_Photon10_L1R_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Photon15_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon15_L1R") )     HLT_Photon15_L1R_ = 1;
          else                                         HLT_Photon15_L1R_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Photon10_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon10_Cleaned_L1R") )     HLT_Photon10_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon10_Cleaned_L1R_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Photon15_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon15_Cleaned_L1R") )     HLT_Photon15_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon15_Cleaned_L1R_ = 0;
        }

        if( strcmp(hlt_trigNames().at(i) , "HLT_Photon20_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon20_Cleaned_L1R") )     HLT_Photon20_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon20_Cleaned_L1R_ = 0;
        }
      }

      //met quantities--------------------------------------------------------------------------
      
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
        //genmet_     = cms2.gen_met();
        //genmetphi_  = cms2.gen_metPhi();
        //gensumet_   = cms2.gen_sumEt();
      }
      
      // if(calculateTCMET){
        
//         // calculate tcmet on-the-fly
//         bool usePV            = false;
//         bool useHFcleaning    = false;
//         bool useHCALcleaning  = false;
//         bool useECALcleaning  = false;
        
//         if ( isData ){
//           useHFcleaning   = true;
//           useHCALcleaning = true;
//           useECALcleaning = true;
//         }
        
//         metStruct structMET = correctedTCMET(usePV, useHFcleaning, useHCALcleaning, useECALcleaning);
        
//         tcmet_     = structMET.met;
//         tcmetphi_  = structMET.metphi;
//         tcsumet_   = structMET.sumet;

//       }else{
        
      //just take tcmet from event
      tcmet_     = evt_tcmet();
      tcmetphi_  = evt_tcmetPhi();
      tcsumet_   = evt_tcsumet();

      //}
      

      //photon quantities-----------------------------------------------------------------------

      if(debug) cout << "Get photon quantities" << endl;

      if(debug) cout << "run lumi event " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
      
      nPhotons_         =  0;
      float maxPhotonPt = -1;
      int igmax         = -1;

      if( photons_p4().size() == 0 && algo_ == e_photonSelection ) continue;
      
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

      //if( igmax < 0 && algo_ == e_photonSelection ) continue;

      if( igmax > -1 ){
        etg_   = photons_p4()[igmax].pt();
        etag_  = photons_p4()[igmax].eta();
        phig_  = photons_p4()[igmax].phi();
        hoe_   = photons_hOverE()[igmax];
        //eciso_ = photons_ecalIso()[igmax];
        //hciso_ = photons_hcalIso()[igmax];
        //tkiso_ = photons_tkIsoSolid()[igmax];
        swiss_ = photons_swissSeed()[igmax];
        
        int scind = photons_scindex()[igmax] ;
        seed_  = scs_eSeed()[scind] ;
        s4_ = swiss_ - seed_ ;
        r4_ = 1 - s4_ / seed_ ;

        drel_ = 1000;

        for( unsigned int iel = 0 ; iel < els_p4().size() ; iel++ ){
        
          if( els_p4().at(iel).pt() < 5 ) continue;

          float deta = etag_ - trks_trk_p4().at(iel).eta();
          float dphi = fabs( phig_ - trks_trk_p4().at(iel).phi() );
          if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
          float dr = sqrt( deta * deta + dphi * dphi );
          if( dr < drel_ ) drel_ = dr;

        }
      }
      
      //photon-matched jet quantities---------------------------------------------------------------

      if(debug) cout << "Get jet quantities" << endl;

      jet_dr_      = 10000;
      int   ijetg  = -1;

      //find jet corresponding to photon
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

        float etajet = pfjets_p4().at(ijet).eta();
        
        if( etajet > 3 ) continue;
        
        float phijet = pfjets_p4().at(ijet).phi();
        float deta   = etajet - etag_;
        float dphi   = phijet - phig_;
        if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
        
        float deltaR = sqrt( deta*deta + dphi*dphi );

        if( deltaR < jet_dr_ ){
          jet_dr_ = deltaR;
          ijetg   = ijet;
        }
      }

    
      if( ijetg < 0 && algo_ == e_photonSelection ) continue;

      if( ijetg > -1 ){
      
        jet_pt_             = pfjets_p4().at(ijetg).pt();
        jet_energy_         = pfjets_p4().at(ijetg).energy();
        jet_eta_            = pfjets_p4().at(ijetg).eta();

        //energy component fractions (add protection for index out of range)
        jet_chg_emfrac_     = pfjets_chargedEmE().at(ijetg)     / jet_energy_;
        jet_chg_hadfrac_    = pfjets_chargedHadronE().at(ijetg) / jet_energy_;
        jet_neu_emfrac_     = pfjets_neutralEmE().at(ijetg)     / jet_energy_;
        jet_neu_hadfrac_    = pfjets_neutralHadronE().at(ijetg) / jet_energy_;
        
        //multiplicities
        jet_nchg_           = pfjets_chargedMultiplicity().at(ijetg);
        jet_nmuon_          = pfjets_muonMultiplicity().at(ijetg);
        jet_nneu_           = pfjets_neutralMultiplicity().at(ijetg);
        
        //deltaPhi( jet - met )
        jet_dphimet_          = deltaPhi( pfjets_p4().at(ijetg).phi() , tcmetphi_);
      
//         if(!isData){
          
//           //deltaR match to genjet
//           int iMin    = -1;
//           float dRmin = -1;
          
//           for (unsigned int igenjet = 0 ; igenjet < genjets_p4().size() ; igenjet++ ){
            
//             LorentzVector vgenjet = genjets_p4().at(igenjet);
            
//             float dR = dRbetweenVectors(pfjets_p4().at(ijetg), vgenjet);
            
//             if(dR < dRmin){
//               iMin = igenjet;
//               dRmin = dR;
//             }
//           }
          
//           if(iMin > -1){
//             jet_dpt_   = jet_pt_ - genjets_p4().at(iMin).pt();
//             jet_drgen_ = dRmin;
//           }
//         }
      }

      //find leading jet------------------------------------------------------------------------

      int imaxjet = -1;
      float maxpt = -1;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

        if( fabs( pfjets_p4().at(ijet).eta() ) > 5.) continue;
        
        if( pfjets_p4().at(ijet).pt() > maxpt ){
          maxpt = pfjets_p4().at(ijet).pt();
          imaxjet = ijet;
        }

      }
      
      if( imaxjet > -1 ){
        jetmax_pt_       = pfjets_p4().at(imaxjet).pt();
        jetmax_dphimet_  = deltaPhi( pfjets_p4().at(imaxjet).phi() , tcmetphi_);
      }
      
      //loop over pfjets, find nJets and sumJetPt-----------------------------------------------

      if(debug) cout << "Get nJets and sumJetPt" << endl;

      nJets_      = 0;
      sumJetPt_   = 0.;
      LorentzVector jetSystem(0.,0.,0.,0.);

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

        LorentzVector vjet = pfjets_p4().at(ijet);
        
        bool skipJet = false;

        //skip jet matched to photon
        if( (int)ijet == ijetg && algo_ == e_photonSelection ){
          skipJet = true;
        }
       
        //skip all hyp leptons
        if( algo_ == e_ZSelection ) { 

          for( unsigned int hypIdx = 0 ; hypIdx < hyp_p4().size() ; hypIdx++ ) {
         
            LorentzVector vlt  = hyp_lt_p4()[hypIdx];
            LorentzVector vll  = hyp_ll_p4()[hypIdx];
            
            if (dRbetweenVectors(vjet, vll) < 0.4) skipJet = true;
            if (dRbetweenVectors(vjet, vlt) < 0.4) skipJet = true;
         
          }
        }
  
        if( skipJet ) continue;

        if( fabs( vjet.eta() ) < 5.){
       
          if ( vjet.pt() > 30. ){
            nJets_++;
          }        
          if ( vjet.pt() > 15. ){
            sumJetPt_ += vjet.pt();
            jetSystem += vjet;
          }
        }
      }

      vecJetPt_ = jetSystem.pt();

      if( algo_ == e_ZSelection ){
      
        if( hyp_p4().size() != 1 ) continue;       
      
        passz_ = passZSelection() ? 1 : 0;
      
        pdgid_ = 0;
      
        if( hyp_type()[0] == 0 ){
          pdgid_ = 13; 
          passm_nom_        = passMuon_Nominal()      ? 1 : 0;
          passm_nomttbar_   = passMuon_NominalTTbar() ? 1 : 0;
        }
        if( hyp_type()[0] == 3 ){
          pdgid_ = 11;
          passe_ttbar_      = passElectron_ttbar()    ? 1 : 0;
          passe_cand01_     = passElectron_cand01()   ? 1 : 0;
        }

        ptll_             = hyp_ll_p4()[0].pt();
        ptlt_             = hyp_lt_p4()[0].pt();
        
      }
      

      //fill baby ntuple
      if ( nJets_ < 1 )     continue;
      if( jetmax_pt_ < 30 ) continue;
      FillBabyNtuple();

      //require at least 2 jets
      if ( nJets_ < 2 )     continue;

      //require leading jet pt > 50
      if( jetmax_pt_ < 50 ) continue;

      dphixmet_  = deltaPhi( tcmetphi_ , jetSystem.phi() );
      metPar_    = tcmet_ * cos( dphixmet_ );
      metPerp_   = tcmet_ * sin( dphixmet_ );
      
      //fill met template-----------------------------------------------------------------------
      
      if( algo_ == e_makeTemplate ) {

        if( HLT_Jet15U_ == 0 ) continue;
        
        //dphixmet_  = deltaPhi( tcmetphi_ , jetSystem.phi() );
        //metPar_    = -1 * tcmet_ * cos( dphixmet_ );
        //metPerp_   = tcmet_ * sin( dphixmet_ );
        
        int iJetBin      = getJetBin( nJets_ );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        
        if( debug ) {
          cout << "nJets " << nJets_ << " sumJetPt " << sumJetPt_ << endl;
          cout << "iJetBin " << iJetBin << " iSumJetPtBin " << iSumJetPtBin << " tcmet " << tcmet_ << endl;
        }

        metTemplate[ iJetBin ][ iSumJetPtBin ]->Fill( tcmet_ );
        metParTemplate[ iJetBin ][ iSumJetPtBin ]->Fill( metPar_ );
        metPerpTemplate[ iJetBin ][ iSumJetPtBin ]->Fill( metPerp_ );
        
      }

      
      //apply good photon selection-------------------------------------------------------------
  
      if( algo_ == e_photonSelection ) {

        if( HLT_Photon10_L1R_ == 0         &&         
            HLT_Photon15_L1R_ == 0         && 
            HLT_Photon10_Cleaned_L1R_ == 0 && 
            HLT_Photon15_Cleaned_L1R_ == 0  ) continue;
          
        if ( etg_ < 20 )                                 continue;
        if ( (1.-r4_) < 0.05 )                           continue;
        if ( hoe_ > 0.1 )                                continue;
        if ( jet_dr_ > 0.5 )                             continue;
        if ( jet_neu_emfrac_ + jet_chg_emfrac_< 0.95 )   continue; 
        if ( drel_ < 0.3 )                               continue;
        //if ( jet_neu_emfrac_ < 0.95 )                    continue; 
        //if ( sumJetPt_ < 200. )                          continue;
        
        //dphixmet_  = deltaPhi( tcmetphi_ , phig_ );
        //metPar_    = tcmet_ * cos( dphixmet_ );
        //metPerp_   = tcmet_ * sin( dphixmet_ );
        
      }

      //apply Z selection-----------------------------------------------------------------------

      if( algo_ == e_ZSelection ) {

        //if( vecJetPt_  < 100)                                       continue;
        if( hyp_p4().size() != 1 )                                  continue;
        if( hyp_lt_id()[0] * hyp_ll_id()[0] > 0 )                   continue;
        if( hyp_type()[0]==1 || hyp_type()[0]==2)                   continue;
        if( hyp_p4()[0].mass() < 76. || hyp_p4()[0].mass() > 106.)  continue;
        if( hyp_ll_p4()[0].pt() < 10 )                              continue;
        if( hyp_lt_p4()[0].pt() < 10 )                              continue;
        //if ( sumJetPt_ < 200. )                                     continue;

        //muon ID
        if (abs(hyp_ll_id()[0]) == 13  && (! (fabs(hyp_ll_p4()[0].eta()) < 2.4 && muonId(hyp_ll_index()[0]))))   continue;
        if (abs(hyp_lt_id()[0]) == 13  && (! (fabs(hyp_lt_p4()[0].eta()) < 2.4 && muonId(hyp_lt_index()[0]))))   continue;
        
        //cand01
        if (abs(hyp_ll_id()[0]) == 11  && (! pass_electronSelection( hyp_ll_index()[0] , electronSelection_cand01 ))) continue;
        if (abs(hyp_lt_id()[0]) == 11  && (! pass_electronSelection( hyp_lt_index()[0] , electronSelection_cand01 ))) continue;

        //ttbar muon ID
        //if (abs(hyp_ll_id()[0]) == 13  && (! (fabs(hyp_ll_p4()[0].eta()) < 2.4 && muonId(hyp_ll_index()[0],NominalTTbar))))   continue;
        //if (abs(hyp_lt_id()[0]) == 13  && (! (fabs(hyp_lt_p4()[0].eta()) < 2.4 && muonId(hyp_lt_index()[0],NominalTTbar))))   continue;
        
        //ttbar electron ID
        //if (abs(hyp_ll_id()[0]) == 11  && (! pass_electronSelection( hyp_ll_index()[0] , electronSelection_ttbar ))) continue;
        //if (abs(hyp_lt_id()[0]) == 11  && (! pass_electronSelection( hyp_lt_index()[0] , electronSelection_ttbar ))) continue;

        //dphixmet_  = deltaPhi( tcmetphi_ , hyp_p4()[0].phi() );
        //metPar_    = tcmet_ * cos( dphixmet_ );
        //metPerp_   = tcmet_ * sin( dphixmet_ );
        
      }

     

      //fill predicted and observed met histos--------------------------------------------------

      if( algo_ != e_makeTemplate){
        
        int iJetBin      = getJetBin( nJets_ );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        TH1F* hmet     = (TH1F*) metTemplateFile->Get(Form("metTemplate_%i_%i",iJetBin,iSumJetPtBin));
        TH1F* hmetPar  = (TH1F*) metTemplateFile->Get(Form("metParTemplate_%i_%i",iJetBin,iSumJetPtBin));
        TH1F* hmetPerp = (TH1F*) metTemplateFile->Get(Form("metPerpTemplate_%i_%i",iJetBin,iSumJetPtBin));
        
        metObserved_njets[iJetBin]->Fill( tcmet_ );
        metPredicted_njets[iJetBin]->Add( hmet );
        
        metObserved->Fill( tcmet_ );
        metPredicted->Add( hmet );

        metParObserved->Fill( metPar_ );
        metParPredicted->Add( hmetPar );

        metPerpObserved->Fill( metPerp_ );
        metPerpPredicted->Add( hmetPerp );
        
        delete hmet;

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
  rootfilename << "root/" << prefix << metTemplateString << ".root";

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist(rootfilename.str().c_str());
  deleteHistos();
  
} // end ScanChain

float looper::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

bool looper::passMuon_NominalTTbar(){

  //ttbar muon ID
  if (abs(hyp_ll_id()[0]) == 13  && (! (fabs(hyp_ll_p4()[0].eta()) < 2.4 && muonId(hyp_ll_index()[0],NominalTTbar))))   return false;
  if (abs(hyp_lt_id()[0]) == 13  && (! (fabs(hyp_lt_p4()[0].eta()) < 2.4 && muonId(hyp_lt_index()[0],NominalTTbar))))   return false;
  return true;

}

bool looper::passMuon_Nominal(){

  //nominal muon ID
  if (abs(hyp_ll_id()[0]) == 13  && (! (fabs(hyp_ll_p4()[0].eta()) < 2.4 && muonId(hyp_ll_index()[0]))))   return false;
  if (abs(hyp_lt_id()[0]) == 13  && (! (fabs(hyp_lt_p4()[0].eta()) < 2.4 && muonId(hyp_lt_index()[0]))))   return false;
  return true;

}

bool looper::passElectron_ttbar(){

  //ttbar electron ID
  if (abs(hyp_ll_id()[0]) == 11  && (! pass_electronSelection( hyp_ll_index()[0] , electronSelection_ttbar ))) return false;
  if (abs(hyp_lt_id()[0]) == 11  && (! pass_electronSelection( hyp_lt_index()[0] , electronSelection_ttbar ))) return false;
  return true;

}

bool looper::passElectron_cand01(){

  //cand01 electron ID
  if (abs(hyp_ll_id()[0]) == 11  && (! pass_electronSelection( hyp_ll_index()[0] , electronSelection_cand01 ))) return false;
  if (abs(hyp_lt_id()[0]) == 11  && (! pass_electronSelection( hyp_lt_index()[0] , electronSelection_cand01 ))) return false;
  return true;

}

bool looper::passZSelection(){
  
  if( hyp_p4().size() != 1 )                                  return false;
  if( hyp_lt_id()[0] * hyp_ll_id()[0] > 0 )                   return false;
  if( hyp_type()[0]==1 || hyp_type()[0]==2)                   return false;
  if( hyp_p4()[0].mass() < 76. || hyp_p4()[0].mass() > 106.)  return false;

  return true;

}

void looper::InitBabyNtuple (){

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
  dphixmet_     = -999999.;
  metPar_       = -999999.;
  metPerp_      = -999999.;
  tcmetphi_     = -999999.;
  tcsumet_      = -999999.;

  nJets_        = -999999;
  sumJetPt_     = -999999;
  vecJetPt_     = -999999;

  //photon stuff
  nPhotons_ = 0;
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
  drel_     = -999999.;
                                  
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
  jet_drgen_        = -999999.;  

  //leading jet stuff
  jetmax_pt_        = -999999;
  jetmax_dphimet_   = -999999;

  //Z stuff
  passz_           = -999999;
  passe_ttbar_     = -999999;
  passe_cand01_    = -999999;
  passm_nomttbar_  = -999999;
  passm_nom_       = -999999;
  pdgid_           = -999999;
  ptll_            = -999999;
  ptlt_            = -999999;
}

void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  metObserved  = new TH1F("metObserved", "Observed MET",500,0,500);
  metPredicted = new TH1F("metPredicted","Predicted MET",500,0,500);

  metParObserved  = new TH1F("metParObserved", "Observed MET (Parallel)",1000,-500,500);
  metParPredicted = new TH1F("metParPredicted","Predicted MET (Parallel)",1000,-500,500);

  metPerpObserved  = new TH1F("metPerpObserved", "Observed MET (Perpendicular)",500,0,500);
  metPerpPredicted = new TH1F("metPerpPredicted","Predicted MET (Perpendicular)",500,0,500);
 
  metObserved->Sumw2();
  metPredicted->Sumw2();

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


void looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("T1", "A Baby Ntuple");

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

  babyTree_->Branch("dphixmet",     &dphixmet_,      "dphixmet/F"    );
  babyTree_->Branch("metpar",       &metPar_,        "metpar/F"      );
  babyTree_->Branch("metperp",      &metPerp_,       "metperp/F"     );


  babyTree_->Branch("tcmetphi",     &tcmetphi_,     "tcmetphi/F"   );
  babyTree_->Branch("tcsumet",      &tcsumet_,      "tcsumet/F"    );
  babyTree_->Branch("njets",        &nJets_,        "njets/I"    );
  babyTree_->Branch("sumjetpt",     &sumJetPt_,     "sumjetpt/F"    );
  babyTree_->Branch("vecjetpt",     &vecJetPt_,     "vecjetpt/F"    );

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
  babyTree_->Branch("jetdrgen",              &jet_drgen_,            "jetdrgen/F");
  babyTree_->Branch("drel",                  &drel_,                 "drel/F");

  babyTree_->Branch("maxjetpt",              &jetmax_pt_,            "maxjetpt/F");
  babyTree_->Branch("maxjetdphimet",         &jetmax_dphimet_,       "maxjetdphimet/F");

  //trigger
  //babyTree_->Branch("HLT_Jet6U",             &HLT_Jet6U_,            "HLT_Jet6U/I");
  //babyTree_->Branch("HLT_Jet10U",            &HLT_Jet10U_,           "HLT_Jet10U/I");
  babyTree_->Branch("HLT_Jet15U",            &HLT_Jet15U_,           "HLT_Jet15U/I");
  babyTree_->Branch("HLT_Jet30U",            &HLT_Jet30U_,           "HLT_Jet30U/I");

  babyTree_->Branch("HLT_Photon10_L1R",      &HLT_Photon10_L1R_,     "HLT_Photon10_L1R/I");
  babyTree_->Branch("HLT_Photon15_L1R",      &HLT_Photon15_L1R_,     "HLT_Photon15_L1R/I");

  babyTree_->Branch("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_,     "HLT_Photon10_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_,     "HLT_Photon15_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_,     "HLT_Photon20_Cleaned_L1R/I");  

  //Z stuff
  babyTree_->Branch("passz",           &passz_,           "passz/I");  
  babyTree_->Branch("pdgid",           &pdgid_,           "pdgid/I");  
  babyTree_->Branch("passm_nom",       &passm_nom_,       "passm_nom/I");  
  babyTree_->Branch("passm_nomttbar",  &passm_nomttbar_,  "passm_nomttbar/I");  
  babyTree_->Branch("passe_ttbar",     &passe_ttbar_,     "passe_ttbar/I");  
  babyTree_->Branch("passe_cand01",    &passe_cand01_,    "passe_cand01/I");  
  babyTree_->Branch("ptll",            &ptll_,            "ptll/F");  
  babyTree_->Branch("ptlt",            &ptlt_,            "ptlt/F");  

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







        
//         if( pfjets_p4().size() != pfjets_chargedEmE().size() && debug ){
//           cout << evt_dataset() << endl;
//           cout << "run lumi event : " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
//           cout << "pfjets_p4().size() " << pfjets_p4().size() << " pfjets_chargedEmE().size() " << pfjets_chargedEmE().size() << endl;
//         }


        
//         //add protection for index out of range
//         if( ijetg < pfjets_chargedEmE().size() ){
