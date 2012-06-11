#include "makePhotonBabies.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <set>

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDatabasePDG.h"
#include <sstream>

#include "../CORE/CMS2.cc"
#ifndef __CINT__
#include "../CORE/utilities.cc"
//#include "../CORE/ssSelections.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/MITConversionUtilities.cc"
#include "../CORE/muonSelections.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/ttbarSelections.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/metSelections.cc"
#include "../CORE/jetSelections.cc"
#include "../CORE/photonSelections.cc"
#include "../CORE/triggerUtils.cc"
#include "../CORE/triggerSuperModel.cc"
#include "../CORE/mcSelections.cc"
#include "../CORE/susySelections.cc"
#include "../CORE/mcSUSYkfactor.cc"
#include "../CORE/SimpleFakeRate.cc"
#include "../Tools/goodrun.cc"
#include "../Tools/vtxreweight.cc"
#include "../Tools/msugraCrossSection.cc"
#include "../Tools/bTagEff_BTV.cc"
#endif

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

using namespace tas;

//--------------------------------------------------------------------

const bool debug                = false;
const float lumi                = 1.0;
const char* iter                = "V00-00-09";
const char* jsonfilename        = "../jsons/Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_goodruns.txt";

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

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

//--------------------------------------------------------------------

int makePhotonBabies::passThisHLTTrigger( char* hltname ){

  if( debug) cout << "Checking for pattern " << hltname << endl;

  //-------------------------------------------------------
  // First check if trigger is present. If not, return -1.
  //-------------------------------------------------------

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( hltname ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger ){
    if( debug ) cout << "Did not find trigger" << endl;
    return -1;
  }

  //----------------------------------------------------
  // Now check if trigger passed. If so return prescale
  //----------------------------------------------------

  if( debug) cout << "Found trigger " << exact_hltname << endl;

  if( !passHLTTrigger( exact_hltname ) ){
    if( debug ) cout << "Trigger did not pass" << endl;
    return 0;
  }
  else{
    int PS = HLT_prescale( exact_hltname );
    if( debug ) cout << "Trigger passed, prescale " << PS << endl;
    return PS;
  }

  //-------------------
  // shouldn't get here
  //-------------------

  return -2;
}

//--------------------------------------------------------------------

void makePhotonBabies::ScanChain (TChain* chain, const char* prefix, bool isData, 
                                  bool calculateTCMET, int nEvents, float kFactor){

  set_goodrun_file( jsonfilename );
  cout << "Using json " << jsonfilename << endl;

  //------------------------------------------------------------------------------------------------------
  // load here the on-the-fly corrections/uncertainties L1FastL2L3 (MC) and L1FastL2L3Residual (DATA)
  // corrections are stored in jet_corrected_pfL1FastJetL2L3
  // uncertainties are stored in pfUncertainty
  //------------------------------------------------------------------------------------------------------

  std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

  jetcorr_filenames_pfL1FastJetL2L3.clear();

  if ( TString(prefix).Contains("data") ) {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L1FastJet.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2Relative.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L3Absolute.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2L3Residual.txt");
  } 
  else {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L1FastJet.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L2Relative.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L3Absolute.txt");
  }

  jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  bookHistos();

  int npass = 0;
  
  //---------------------
  // make a baby ntuple
  //---------------------
  
  MakeBabyNtuple( Form("../photon_output/%s/%s_baby.root", iter , prefix ) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  if(debug) cout << "Begin file loop" << endl;

  //------------------
  // begin file loop
  //------------------

  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())){

    cout << currentFile->GetTitle() << endl;

    TFile* f = new TFile(currentFile->GetTitle());

    if( !f || f->IsZombie() ) {
      cout << "Skipping bad input file: " << currentFile->GetTitle() << endl;
      continue; //exit(1);                                                                                             
    }

    TTree *tree = (TTree*)f->Get("Events");

    if( tree == 0 ){
      cout << "Can't get tree from file " << currentFile->GetTitle() << endl;
      continue;
    }

    cms2.Init(tree);

    unsigned int nEvents = tree->GetEntries();
 
    for (unsigned int event = 0 ; event < nEvents; ++event){
        
      cms2.GetEntry(event);
      ++nEventsTotal;

      //-------------------------------
      // progress feedback to user
      //-------------------------------

      if (nEventsTotal % 1000 == 0){

        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }
      
      //------------------------------
      // duplicate event cleaning
      //------------------------------
      
      if( isData ) {
	DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
	if (is_duplicate(id) )
	  continue;
      }
      
      //--------------------------
      // good run+event selection
      //--------------------------

      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      if( !cleaning_goodVertexApril2011() )                          continue;

      if(debug) cout << "Pass event selection" << endl;

      InitBabyNtuple();

      //----------------
      // event stuff
      //----------------

      strcpy(dataset_, cms2.evt_dataset().at(0).Data());
      run_     = cms2.evt_run();
      lumi_    = cms2.evt_lumiBlock();
      event_   = cms2.evt_event();

      weight_ = 1.;
      pthat_  = -1;
      
      if( !isData ){
        weight_ = cms2.evt_scale1fb() * kFactor * lumi;
        pthat_  = cms2.genps_pthat();
      }

      //-------------------------------------------------------
      // stitching together Photon MC
      //-------------------------------------------------------
      
      float maxphotonpt = -1;

      if( strcmp( prefix , "PhotonJet") == 0 ){

        if( TString(currentFile->GetTitle()).Contains("PhotonJet_Pt15") ){
          if( cms2.genps_pthat() > 30. ) continue;
        }
        if( TString(currentFile->GetTitle()).Contains("PhotonJet_Pt30") ){
          if( cms2.genps_pthat() > 80. ) continue;
        }
        if( TString(currentFile->GetTitle()).Contains("PhotonJet_Pt80") ){
          if( cms2.genps_pthat() > 170. ) continue;
        }

        for( unsigned int ig = 0 ; ig < photons_p4().size() ; ++ig ){
          if( photons_p4().at(ig).pt() > maxphotonpt ) 
            maxphotonpt = photons_p4().at(ig).pt();
        }
 
        fillUnderOverFlow( hgenps_pthat  , genps_pthat() , weight_ );
        fillUnderOverFlow( hphotonpt     , maxphotonpt   , weight_ );
      }
            
      //-------------------------
      // access HLT triggers
      //-------------------------
     
      hlt20_  = passThisHLTTrigger( "HLT_Photon20_CaloIdVL_IsoL_v"  );
      hlt30_  = passThisHLTTrigger( "HLT_Photon30_CaloIdVL_IsoL_v"  );
      hlt50_  = passThisHLTTrigger( "HLT_Photon50_CaloIdVL_IsoL_v"  );
      hlt75_  = passThisHLTTrigger( "HLT_Photon75_CaloIdVL_IsoL_v"  );      
      hlt90_  = passThisHLTTrigger( "HLT_Photon90_CaloIdVL_IsoL_v"  );
      hlt135_ = passThisHLTTrigger( "HLT_Photon135_v"               );
      hlt150_ = passThisHLTTrigger( "HLT_Photon150_v"               );
      hlt160_ = passThisHLTTrigger( "HLT_Photon160_v"               );

      hgg22_  = passThisHLTTrigger( "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v" );
      hgg36_  = passThisHLTTrigger( "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v" );
      hgg50_  = passThisHLTTrigger( "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v" );
      hgg75_  = passThisHLTTrigger( "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v" );
      hgg90_  = passThisHLTTrigger( "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v" );

      rho_ = evt_ww_rho_vor(); // TO BE REPLACED
      //rho_ = evt_kt6pf_foregiso_rho();

      maxleppt_ = 0;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                                 continue;
	if( !pass_electronSelection( iel , electronSelection_el_OSV2 , false , false ) ) continue;
	if( els_p4().at(iel).pt() > maxleppt_ ) maxleppt_ = els_p4().at(iel).pt();
      }
              
      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 10 )           continue;
	if( !muonId( imu , OSZ_v4 ))               continue;
	if( mus_p4().at(imu).pt() > maxleppt_ ) maxleppt_ = mus_p4().at(imu).pt();
      }
              
      //-------------------------
      // MET variables
      //-------------------------

      met_        = cms2.evt_met();
      metphi_     = cms2.evt_metPhi();
      sumet_      = cms2.evt_sumet();

      pfmet_      = cms2.evt_pfmet();
      pfmetphi_   = cms2.evt_pfmetPhi();
      pfsumet_    = cms2.evt_pfsumet();

      pfmett1_    = -99; //cms2.evt_pfmet_type1cor();
      pfmetphit1_ = -99; //cms2.evt_pfmetPhi_type1cor();

      tcmet_     = evt_tcmet();
      tcmetphi_  = evt_tcmetPhi();
      tcsumet_   = evt_tcsumet();

      if (!isData){
        genmet_     = cms2.gen_met();
        genmetphi_  = cms2.gen_metPhi();
        gensumet_   = cms2.gen_sumEt();
      }
                                    
      //------------------------
      // vertex stuff
      //------------------------

      nGoodVertex_ = 0;

      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
        if(isGoodVertex(v)) nGoodVertex_++;
      }

      //---------------------------------
      // photon stuff
      //---------------------------------

      int   ijetg         = -1; //index of jet matched to photon
      float maxPhotonPt   = -1;
      int   igmax         = -1;

      nPhotons_           =  0;
      
      if( photons_p4().size() == 0 ) continue;
      
      //---------------------------
      //count photons pt > 10 GeV
      //---------------------------
      
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
      
      ijetg = isGoodEMObject( igmax );
      
      if( ijetg < 0 ) continue;
      
      etg_        = photons_p4().at(igmax).pt();
      etag_       = photons_p4().at(igmax).eta();
      phig_       = photons_p4().at(igmax).phi();
      hoe_        = photons_hOverE().at(igmax);
      swiss_      = photons_swissSeed().at(igmax);
      int scind   = photons_scindex().at(igmax) ;
      
      if( scind > - 1 ){
        seed_       = scs_eSeed().at(scind) ;
        s4_         = swiss_ - seed_ ;
        r4_         = 1 - s4_ / seed_ ;
	photon_sigmaIPhiIPhi_ = scs_sigmaIPhiIPhi()[scind]; 
     }else{
        seed_                 = -9999.;
        s4_                   = -9999.;
        r4_                   = -9999.;
	photon_sigmaIPhiIPhi_ = -9999.;
      }
      
      photon_scidx_            = photons_scindex().at(igmax);
      photon_pixelseed_        = photons_haspixelSeed().at(igmax) ? 1 : 0;
      photon_e15_              = photons_e1x5().at(igmax);      
      photon_e25max_           = photons_e2x5Max().at(igmax);      
      photon_e33_              = photons_e3x3().at(igmax);           
      photon_e55_              = photons_e5x5().at(igmax);           
      photon_ecalIso03_        = photons_ecalIso03().at(igmax);      
      photon_ecalIso04_        = photons_ecalIso04().at(igmax);      
      photon_hcalIso03_        = photons_hcalIso03().at(igmax);      
      photon_hcalIso04_        = photons_hcalIso04().at(igmax);      
      photon_ntkIsoHollow03_   = photons_ntkIsoHollow03().at(igmax);
      photon_ntkIsoHollow04_   = photons_ntkIsoHollow04().at(igmax);
      photon_ntkIsoSolid03_    = photons_ntkIsoSolid03().at(igmax); 
      photon_ntkIsoSolid04_    = photons_ntkIsoSolid04().at(igmax); 
      photon_sigmaEtaEta_      = photons_sigmaEtaEta().at(igmax);    
      photon_sigmaIEtaIEta_    = photons_sigmaIEtaIEta().at(igmax);
      photon_tkisoHollow03_    = photons_tkIsoHollow03().at(igmax);
      photon_tkisoHollow04_    = photons_tkIsoHollow04().at(igmax); 
      photon_tkisoSolid03_     = photons_tkIsoSolid03().at(igmax);   
      photon_tkisoSolid04_     = photons_tkIsoSolid04().at(igmax);  

      elveto_ = 0;
      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                                 continue;
	if( dRbetweenVectors( photons_p4().at(igmax) , els_p4().at(iel) ) > 0.2 )        continue;
	elveto_ = 1;
      }
      
      LorentzVector myvjet = pfjets_p4().at(ijetg);
      LorentzVector myvg   = photons_p4().at(igmax);
      
      jet_dr_             = dRbetweenVectors(myvjet, myvg);
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
      jet_pfjetid_        = passesPFJetID( ijetg ) ? 1 : 0;
      

      // get pt of closest calojet, within dr < 0.3
      calojet_pt_          = -1;

      float mindr = 100;

      for( int ic = 0 ; ic < jets_p4().size() ; ic++ ){
	float dr = dRbetweenVectors( myvg, jets_p4().at(ic) );

	if( dr > 0.3 ) continue;

	if( dr < mindr ){
	  mindr = dr;
	  calojet_pt_ = jets_p4().at(ic).pt();
	}

      }

      //--------------------
      // jet stuff
      //--------------------

      nJets10_      = 0;
      nJets15_      = 0;
      nJets20_      = 0;
      nJets_        = 0;
      nJets40_      = 0;
      ht_           = 0.;
      ht10_         = 0.;
      nbl_          = 0;
      nbm_          = 0;
      nbt_          = 0;

      LorentzVector jetSystem(0.,0.,0.,0.);        
      float maxcosdphi  = -99;

      VofP4 goodJets;
      goodJets.clear();

      failjetid_ =  0;
      maxemf_    = -1;
      
      //-----------------------------------------
      // loop over pfjets pt > 30 GeV |eta| < 3.0
      //-----------------------------------------

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

	//---------------------------------------------------------------------------
        // skip jet matched to photon
	//---------------------------------------------------------------------------

        if( (int)ijet == ijetg ) continue;
 
	if( fabs( pfjets_p4().at(ijet).eta() ) > 5.0 ) continue;

	//---------------------------------------------------------------------------
	// get total correction: L1FastL2L3 for MC, L1FastL2L3Residual for data
	// and require |eta| < 3.0
	//---------------------------------------------------------------------------

	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();

	LorentzVector vjet = corr * pfjets_p4().at(ijet);

	if( fabs(vjet.eta()) > 3.0 ) continue;

	//---------------------------------------------------------------------------
        // PFJetID
	//---------------------------------------------------------------------------

        if( !passesPFJetID(ijet) ){
          failjetid_ = 1;
          continue;
        }

	//---------------------------------------------------------------------------
        // HT variables
	//---------------------------------------------------------------------------

        if ( vjet.pt() > 10. ){
	  nJets10_++;
          ht10_ += vjet.pt();
	}

        if ( vjet.pt() > 15. ){
	  nJets15_++;
          ht_ += vjet.pt();
          jetSystem += vjet;

          float emfrac = pfjets_neutralEmE().at(ijet) / pfjets_p4().at(ijet).energy();
          if( emfrac > maxemf_ ) maxemf_ = emfrac;
	}

        if ( vjet.pt() > 20. ) nJets20_++;
        if ( vjet.pt() > 30. ) nJets_++;
        if ( vjet.pt() > 40. ) nJets40_++;
              
        if( vjet.pt() < 30. )                    continue;

	goodJets.push_back(vjet);

        //find jet (anti-)aligned with tcmet
        if( fabs( cos( evt_pfmetPhi() - vjet.phi() ) ) > maxcosdphi ){
          maxcosdphi  = fabs( cos( evt_pfmetPhi() - vjet.phi() ) );
          dphijetmet_ = fabs( evt_pfmetPhi() - vjet.phi() );
          if( dphijetmet_ > TMath::Pi() ) dphijetmet_ = TMath::TwoPi() - dphijetmet_;
        }
        
	if( pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.244 ) nbl_ ++;
	if( pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.679 ) nbm_ ++;
	if( pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.898 ) nbt_ ++;

      }

      vecJetPt_ = jetSystem.pt();

      sort(goodJets.begin()    , goodJets.end()    , sortByPt);

      if( goodJets.size()  > 0 ) jet1_   = &(goodJets.at(0));
      if( goodJets.size()  > 1 ) jet2_   = &(goodJets.at(1));
      if( goodJets.size()  > 2 ) jet3_   = &(goodJets.at(2));
      if( goodJets.size()  > 3 ) jet4_   = &(goodJets.at(3));
      
      //-------------------------
      // fill histos and ntuple
      //-------------------------
              
      npass++;
      FillBabyNtuple();

    } // end loop over events

    delete f;

  } // end loop over files

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  CloseBabyNtuple();

  cout << "numEvents passing selection " << npass << endl;

  //deleteHistos();
  
} // end ScanChain

//--------------------------------------------------------------------

float makePhotonBabies::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//--------------------------------------------------------------------

void makePhotonBabies::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

//--------------------------------------------------------------------

void makePhotonBabies::InitBabyNtuple (){

  // pfjets
  jet1_                         = 0;
  jet2_                         = 0;
  jet3_                         = 0;
  jet4_                         = 0;

  // trigger stuff
  hlt20_			= -9999;
  hlt30_			= -9999;
  hlt50_			= -9999;
  hlt75_			= -9999;
  hlt90_			= -9999;
  hlt135_			= -9999;
  hlt150_			= -9999;
  hlt160_			= -9999;

  // event stuff
  run_				= -999999;
  memset(dataset_, '\0', 200);
  lumi_				= -999999;
  event_			= -999999;
  weight_			= -999999.;
  pthat_			= -999999.;
  nGoodVertex_			= -999999;
  leptype_			= -999999;

  // genmet stuff
  genmet_			= -999999.;
  genmetphi_			= -999999.;
  gensumet_			= -999999.;

  // pfmet stuff
  pfmet_			= -999999.;
  pfmet_type1_pt30_		= -999999.;
  pfmet_type1_pt15_		= -999999.;
  pfmetphi_			= -999999.;
  pfsumet_			= -999999.;
  pfmett1_			= -999999.;
  pfmetphit1_			= -999999.;

  // calomet stuff
  met_				= -999999.;
  metphi_			= -999999.;
  sumet_			= -999999.;

  // muon-corrected calomet stuff
  mumet_			= -999999.;
  mumetphi_			= -999999.;
  musumet_			= -999999.;

  // calomet stuff
  mujesmet_			= -999999.;
  mujesmetphi_			= -999999.;
  mujessumet_			= -999999.;

  // tcmet stuff
  dphixmet_			= -999999.;
  metPar_			= -999999.;
  metPerp_			= -999999.;

  tcmet_			= -999999.;
  tcmetphi_			= -999999.;
  tcsumet_			= -999999.;

  tcmetNew_			= -999999.;
  tcmetNew_type1_pt30_		= -999999.;
  tcmetNew_type1_pt15_		= -999999.;
  tcsumetNew_			= -999999.;
  tcmetphiNew_			= -999999.;

  nJets_			= -999999;
  ht_		        	= -999999;
  vecJetPt_			= -999999;
  nJets40_			= -999999;
  nJets10_			= -999999;
  nJets15_			= -999999;
  nJets20_			= -999999;
  ht10_	         		= -999999;

  nbl_  			= -999999;
  nbm_  			= -999999;
  nbt_  			= -999999;
  dphijetmet_			= -999999;

  //leading jet stuff
  jetmax_pt_			= -999999;
  jetmax_dphimet_		= -999999;

  failjetid_			= -999999;
  maxemf_			= -999999.;

  //photon stuff
  nPhotons_			= -999999;
  etg_				= -999999.;
  etag_				= -999999.;
  phig_				= -999999.;
  hoe_				= -999999.;
  eciso_			= -999999.;
  hciso_			= -999999.;
  tkiso_			= -999999.;
  swiss_			= -999999.;
  seed_				= -999999.;
  s4_				= -999999.;
  r4_				= -999999.;

  //more photon stuff
  photon_scidx_			= -999999;
  photon_pixelseed_		= -999999;
  photon_e15_			= -999999.;
  photon_e25max_		= -999999.;
  photon_e33_			= -999999.;
  photon_e55_			= -999999.;
  photon_ecalIso03_		= -999999.;
  photon_ecalIso04_		= -999999.;
  photon_hcalIso03_		= -999999.;
  photon_hcalIso04_		= -999999.;
  photon_ntkIsoHollow03_	= -999999.;
  photon_ntkIsoHollow04_	= -999999.;
  photon_ntkIsoSolid03_		= -999999.;
  photon_ntkIsoSolid04_		= -999999.;
  photon_sigmaEtaEta_		= -999999.;
  photon_sigmaIEtaIEta_		= -999999.;
  photon_sigmaIPhiIPhi_		= -999999.;
  photon_tkisoHollow03_		= -999999.;
  photon_tkisoHollow04_		= -999999.;
  photon_tkisoSolid03_		= -999999.;
  photon_tkisoSolid04_		= -999999.;
                                  
  //photon-matched jet stuff
  jet_eta_			= -999999.;  
  jet_energy_			= -999999.;  
  jet_pt_			= -999999.;  
  calojet_pt_			= -999999.;  
  jet_chg_emfrac_		= -999999.;  
  jet_chg_hadfrac_		= -999999.;  
  jet_neu_emfrac_		= -999999.;  
  jet_neu_hadfrac_		= -999999.;  
  jet_nchg_			= -999999;      
  jet_nmuon_			= -999999;  
  jet_nneu_			= -999999;  
  jet_dphimet_			= -999999.;  
  jet_pfjetid_			= -999999;  
  jet_dpt_			= -999999.;  

}

//--------------------------------------------------------------------

void makePhotonBabies::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hgenps_pthat = new TH1F("hgenps_pthat","",500,0,500);
  hphotonpt    = new TH1F("hphotonpt","",500,0,500);

  hgenps_pthat->GetXaxis()->SetTitle("gen p_{T}(hat) (GeV)");
  hphotonpt->GetXaxis()->SetTitle("max photon p_{T} (GeV)");

  char* pttitle[5]={"all jets","1 jet","2 jet","3 jet","#geq 4 jet"};

  for( int iJ = 0 ; iJ < 5 ; iJ++ ){
    hetg[iJ] = new TH1F(Form("hetg_%i",iJ),pttitle[iJ],200,0,200);
    hetg[iJ]->GetXaxis()->SetTitle("photon p_{T} (GeV)");
  }

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

void makePhotonBabies::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("T1", "A Baby Ntuple");


  //event stuff
  babyTree_->Branch("dataset"	,	&dataset_          ,	"dataset[200]/C"   );
  babyTree_->Branch("run"	,       &run_              ,	"run/I"            );
  babyTree_->Branch("lumi"	,       &lumi_             ,	"lumi/I"           );
  babyTree_->Branch("event"	,       &event_            ,	"event/I"          );
  babyTree_->Branch("nvtx"	,       &nGoodVertex_      ,	"nvtx/I"           );
  babyTree_->Branch("weight"	,       &weight_           ,	"weight/F"         );
  babyTree_->Branch("pthat"	,       &pthat_            ,	"pthat/F"          );
  babyTree_->Branch("failjetid"	,	&failjetid_        ,	"failjetid/I"      );
  babyTree_->Branch("maxemf"	,       &maxemf_           ,	"maxemf/F"         );
  babyTree_->Branch("maxleppt"	,	&maxleppt_         ,	"maxleppt/F"       );
  babyTree_->Branch("elveto"	,	&elveto_           ,	"elveto/I"         );

  //met stuff
  babyTree_->Branch("pfmet"			,       &pfmet_                ,	 "pfmet/F"			);
  babyTree_->Branch("pfmett1"			,       &pfmett1_              ,	 "pfmett1/F"			);
  babyTree_->Branch("pfmet_type1_pt30"		,       &pfmet_type1_pt30_     ,	 "pfmet_type1_pt30/F"		);
  babyTree_->Branch("pfmet_type1_pt15"		,       &pfmet_type1_pt15_     ,	 "pfmet_type1_pt15/F"		);
  babyTree_->Branch("pfmetphi"			,	&pfmetphi_             ,	 "pfmetphi/F"			);
  babyTree_->Branch("pfmetphit1"		,	&pfmetphit1_           ,	 "pfmetphit1/F"			);
  babyTree_->Branch("pfsumet"			,	&pfsumet_              ,	 "pfsumet/F"			);
  babyTree_->Branch("met"			,       &met_                  ,	 "met/F"			);
  babyTree_->Branch("metphi"			,       &metphi_               ,	 "metphi/F"			);
  babyTree_->Branch("sumet"			,       &sumet_                ,	 "sumet/F"			);
  babyTree_->Branch("mumet"			,       &mumet_                ,	 "mumet/F"			);
  babyTree_->Branch("mumetphi"			,	&mumetphi_             ,	 "mumetphi/F"			);
  babyTree_->Branch("musumet"			,	&musumet_              ,	 "musumet/F"			);
  babyTree_->Branch("mujesmet"			,	&mujesmet_             ,	 "mujesmet/F"			);
  babyTree_->Branch("mujesmetphi"		,	&mujesmetphi_          ,	 "mujesmetphi/F"		);
  babyTree_->Branch("mujessumet"		,	&mujessumet_           ,	 "mujessumet/F"			);
  babyTree_->Branch("genmet"			,       &genmet_               ,	 "genmet/F"			);
  babyTree_->Branch("genmetphi"			,	&genmetphi_            ,	 "genmetphi/F"			);
  babyTree_->Branch("gensumet"			,	&gensumet_             ,	 "gensumet/F"			);
  babyTree_->Branch("dphixmet"			,	&dphixmet_             ,	 "dphixmet/F"			);
  babyTree_->Branch("metpar"			,       &metPar_               ,	 "metpar/F"			);
  babyTree_->Branch("metperp"			,	&metPerp_              ,	 "metperp/F"			);
  babyTree_->Branch("tcmet"			,       &tcmet_                ,	 "tcmet/F"			);
  babyTree_->Branch("tcmetphi"			,	&tcmetphi_             ,	 "tcmetphi/F"			);
  babyTree_->Branch("tcsumet"			,	&tcsumet_              ,	 "tcsumet/F"			);
  babyTree_->Branch("tcmetNew"			,	&tcmetNew_             ,	 "tcmetNew/F"			);
  babyTree_->Branch("tcmetNew_type1_pt30"	,	&tcmetNew_type1_pt30_  ,	 "tcmetNew_type1_pt30/F"	);
  babyTree_->Branch("tcmetNew_type1_pt15"	,	&tcmetNew_type1_pt15_  ,	 "tcmetNew_type1_pt15/F"	);
  babyTree_->Branch("tcmetphiNew"		,	&tcmetphiNew_          ,	 "tcmetphiNew/F"		);
  babyTree_->Branch("tcsumetNew"		,	&tcsumetNew_           ,	 "tcsumetNew/F"			);

  //jet stuff
  babyTree_->Branch("njets"			,       &nJets_                ,         "njets/I"		);
  babyTree_->Branch("njets10"			,       &nJets10_              ,	 "njets10/I"		);
  babyTree_->Branch("njets15"			,       &nJets15_              ,	 "njets15/I"		);
  babyTree_->Branch("njets20"			,       &nJets20_              ,	 "njets20/I"		);
  babyTree_->Branch("njets40"			,       &nJets40_              ,	 "njets40/I"		);
  babyTree_->Branch("ht"			,       &ht_                   ,	 "ht/F"		        );
  babyTree_->Branch("ht10"	         	,	&ht10_                 ,	 "ht10/F"		);
  babyTree_->Branch("vecjetpt"			,       &vecJetPt_             ,	 "vecjetpt/F"		);
  babyTree_->Branch("nbl"			,       &nbl_                  ,         "nbl/I"		);
  babyTree_->Branch("nbm"			,       &nbm_                  ,         "nbm/I"		);
  babyTree_->Branch("nbt"			,       &nbt_                  ,         "nbt/I"		);
  babyTree_->Branch("ndphijetmet"		,	&dphijetmet_           ,	 "dphijetmet/F"		);
  babyTree_->Branch("maxjetpt"			,       &jetmax_pt_            ,	 "maxjetpt/F"		);
  babyTree_->Branch("maxjetdphimet"		,	&jetmax_dphimet_       ,	 "maxjetdphimet/F"	);
                                   
  //trigger stuff
  babyTree_->Branch("hlt20"			,	&hlt20_  ,  "hlt20/I"    );  
  babyTree_->Branch("hlt30"			,	&hlt30_  ,  "hlt30/I"    );  
  babyTree_->Branch("hlt50"			,	&hlt50_  ,  "hlt50/I"    );  
  babyTree_->Branch("hlt75"			,	&hlt75_  ,  "hlt60/I"    );  
  babyTree_->Branch("hlt90"			,	&hlt90_  ,  "hlt90/I"    );  
  babyTree_->Branch("hlt135"			,	&hlt135_ ,  "hlt135/I"   );  
  babyTree_->Branch("hlt150"			,	&hlt150_ ,  "hlt150/I"   );  
  babyTree_->Branch("hlt160"			,	&hlt160_ ,  "hlt160/I"   );  

  babyTree_->Branch("hgg22"			,	&hgg22_  ,  "hgg22/I"    );  
  babyTree_->Branch("hgg36"			,	&hgg36_  ,  "hgg36/I"    );  
  babyTree_->Branch("hgg50"			,	&hgg50_  ,  "hgg50/I"    );  
  babyTree_->Branch("hgg75"			,	&hgg75_  ,  "hgg75/I"    );  
  babyTree_->Branch("hgg90"			,	&hgg90_  ,  "hgg90/I"    );  

  babyTree_->Branch("rho"			,	&rho_    ,  "rho/F"      );  

  //photon stuff
  babyTree_->Branch("ng"			,	&nPhotons_, "ng/I"); 
  babyTree_->Branch("etg"			,	&etg_,      "etg/F");	   
  babyTree_->Branch("etag"			,	&etag_,     "etag/F");	   
  babyTree_->Branch("phig"			,	&phig_,     "phig/F");
  babyTree_->Branch("hoe"			,	&hoe_,      "hoe/F");	   
  babyTree_->Branch("eciso"			,	&eciso_,    "eciso/F");	   
  babyTree_->Branch("hciso"			,	&hciso_,    "hciso/F");	   
  babyTree_->Branch("tkiso"			,	&tkiso_,    "tkiso/F");
  babyTree_->Branch("swiss"			,	&swiss_,    "swiss/F");
  babyTree_->Branch("seed"			,	&seed_,     "seed/F");
  babyTree_->Branch("s4"			,	&s4_,       "s4/F");
  babyTree_->Branch("r4"			,	&r4_,       "r4/F");

  //more photon stuff
  babyTree_->Branch("photon_scidx"		,       &photon_scidx_,             "photon_scidx/I");         
  babyTree_->Branch("photon_pixelseed"		,       &photon_pixelseed_,         "photon_pixelseed/I");         
  babyTree_->Branch("photon_e15"		,       &photon_e15_,               "photon_e15/F");                
  babyTree_->Branch("photon_e25max"		,       &photon_e25max_,            "photon_e25max/F");             
  babyTree_->Branch("photon_e33"		,       &photon_e33_,               "photon_e33/F");                
  babyTree_->Branch("photon_e55"		,       &photon_e55_,               "photon_e55/F");                
  babyTree_->Branch("photon_ecalIso03"		,       &photon_ecalIso03_,         "photon_ecalIso03/F");          
  babyTree_->Branch("photon_ecalIso04"		,       &photon_ecalIso04_,         "photon_ecalIso04/F");          
  babyTree_->Branch("photon_hcalIso03"		,       &photon_hcalIso03_,         "photon_hcalIso03/F");          
  babyTree_->Branch("photon_hcalIso04"		,       &photon_hcalIso04_,         "photon_hcalIso04/F");          
  babyTree_->Branch("photon_ntkIsoHollow03"	,	&photon_ntkIsoHollow03_,    "photon_ntkIsoHollow03/F");     
  babyTree_->Branch("photon_ntkIsoHollow04"	,	&photon_ntkIsoHollow04_,    "photon_ntkIsoHollow04/F");     
  babyTree_->Branch("photon_ntkIsoSolid03"	,	&photon_ntkIsoSolid03_,     "photon_ntkIsoSolid03/F");      
  babyTree_->Branch("photon_ntkIsoSolid04"	,	&photon_ntkIsoSolid04_,     "photon_ntkIsoSolid04/F");      
  babyTree_->Branch("photon_sigmaEtaEta"	,       &photon_sigmaEtaEta_,       "photon_sigmaEtaEta/F");        
  babyTree_->Branch("photon_sigmaIEtaIEta"	,	&photon_sigmaIEtaIEta_,     "photon_sigmaIEtaIEta/F");      
  babyTree_->Branch("photon_sigmaIPhiIPhi"	,	&photon_sigmaIPhiIPhi_,     "photon_sigmaIPhiIPhi/F");      
  babyTree_->Branch("photon_tkisoHollow03"	,	&photon_tkisoHollow03_,     "photon_tkisoHollow03/F");      
  babyTree_->Branch("photon_tkisoHollow04"	,	&photon_tkisoHollow04_,     "photon_tkisoHollow04/F");      
  babyTree_->Branch("photon_tkisoSolid03"	,	&photon_tkisoSolid03_,      "photon_tkisoSolid03/F");      
  babyTree_->Branch("photon_tkisoSolid04"	,	&photon_tkisoSolid04_,      "photon_tkisoSolid04/F");           
                                                                            
  //photon-matched jet stuff
  babyTree_->Branch("jetdr"			,       &jet_dr_,               "jetdr/F");
  babyTree_->Branch("jetpt"			,       &jet_pt_,               "jetpt/F");
  babyTree_->Branch("calojetpt"			,       &calojet_pt_,           "calojetpt/F");
  babyTree_->Branch("pfjetid"			,       &jet_pfjetid_,          "pfjetid/I");
  babyTree_->Branch("jeteta"			,       &jet_eta_,              "jeteta/F");
  babyTree_->Branch("jetenergy"			,       &jet_energy_,           "jetenergy/F");
  babyTree_->Branch("jetchargedemfrac"		,	&jet_chg_emfrac_,       "jetchargedemfrac/F");
  babyTree_->Branch("jetchargedhadfrac"		,	&jet_chg_hadfrac_,      "jetchargedhadfrac/F");
  babyTree_->Branch("jetneutralemfrac"		,	&jet_neu_emfrac_,       "jetneutralemfrac/F");
  babyTree_->Branch("jetneutralhadfrac"		,	&jet_neu_hadfrac_,      "jetneutralhadfrac/F");
  babyTree_->Branch("jetncharged"		,       &jet_nchg_,             "jetncharged/I");
  babyTree_->Branch("jetnmuon"			,       &jet_nmuon_,            "jetnmuon/I");
  babyTree_->Branch("jetnneutral"		,       &jet_nneu_,             "jetnneutral/I");
  babyTree_->Branch("jetdphimet"		,       &jet_dphimet_,          "jetdphimet/F");
  babyTree_->Branch("jetdpt"			,       &jet_dpt_,              "jetdpt/F");

  babyTree_->Branch("jet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet1_	);
  babyTree_->Branch("jet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet2_	);
  babyTree_->Branch("jet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet3_	);
  babyTree_->Branch("jet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet4_	);

}

//--------------------------------------------------------------------

void makePhotonBabies::FillBabyNtuple ()
{
  babyTree_->Fill();
}

//--------------------------------------------------------------------

void makePhotonBabies::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}

//--------------------------------------------------------------------

void makePhotonBabies::fillHistos(TH1F *h1[4],float value, float weight, int myType)
{

  fillUnderOverFlow(h1[myType], value, weight);      
  fillUnderOverFlow(h1[3],      value, weight);      
}

//--------------------------------------------------------------------

void makePhotonBabies::fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{

  if( nJetsIdx > 2 ) nJetsIdx = 2;
  
  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      

}

//--------------------------------------------------------------------
