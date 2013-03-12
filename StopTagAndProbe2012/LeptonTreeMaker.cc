//now make the source file
#include "LeptonTreeMaker.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTreeCache.h"
#include "TPRegexp.h"

#include "analysisObjects.h"
#include "analysisTools.h"
#include "analysisSelections.h"

using namespace std;

#include "math.h"

//#include "../Tools/MuonIDMVA.h"
#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
//#include "../Tools/EGammaMvaEleEstimator.h"
//#include "../Tools/MuonMVAEstimator.h"
//#include "../Tools/MuonEffectiveArea.h"

//#ifndef __CINT__
#include "../CORE/CMS2.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/MITConversionUtilities.h"
#include "../CORE/muonSelections.h"
#include "../CORE/susySelections.h"
//#endif

MuonIDMVA* muonIdMVA_leptree = 0;
EGammaMvaEleEstimator* egammaMvaEleEstimator_leptree = 0;
MuonMVAEstimator* muonMVAEstimator_leptree = 0;

bool verbose = true;

struct myTrackIso {

  //defaultValue 
  float iso_dr03_dz005_pt00;

  // iso sum options
  float isoDir_dr03_dz005_pt00;

  // r04 cone option
  float iso_dr04_dz005_pt00;

  // veto cone
  float iso_dr01503_dz005_pt00;
  float iso_dr0503_dz005_pt00;

  // dz variation
  float iso_dr03_dz000_pt00;
  float iso_dr03_dz010_pt00;
  float iso_dr03_dz020_pt00;

  //pt Variation
  float iso_dr03_dz005_pt01;
  float iso_dr03_dz005_pt02;
  float iso_dr03_dz005_pt03;
  float iso_dr03_dz005_pt04;
  float iso_dr03_dz005_pt05;
  float iso_dr03_dz005_pt06;
  float iso_dr03_dz005_pt07;
  float iso_dr03_dz005_pt08;
  float iso_dr03_dz005_pt09;
  float iso_dr03_dz005_pt10;

};

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
//
// constructor and destructor
//
inline bool sortByPt(const LorentzVector &vec1, const LorentzVector &vec2 ) {
    return vec1.pt() > vec2.pt();
}


double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

struct myTrackIso trackIso_struct( int thisPf , float coneR = 0.3 , float dz_thresh = 0.05 , bool dovtxcut=false , float pt_thresh = 0.0 ){

  dovtxcut=false;

  struct myTrackIso iso;

  // no cuts  
  iso.iso_dr03_dz005_pt00=0.; // this one will be the default
  
  // dz variation                                                                                                                         
  iso.iso_dr03_dz020_pt00=0.;
  iso.iso_dr03_dz010_pt00=0.;
  iso.iso_dr03_dz000_pt00=0.;
  
  // isolation sum option
  iso.isoDir_dr03_dz005_pt00=0.;

  // cone04
  iso.iso_dr04_dz005_pt00=0.;

  // veto cone
  iso.iso_dr01503_dz005_pt00=0.; 
  iso.iso_dr0503_dz005_pt00=0.;

  //pt Variation                                                                                                                      
  iso.iso_dr03_dz005_pt01=0.;
  iso.iso_dr03_dz005_pt02=0.;
  iso.iso_dr03_dz005_pt03=0.;
  iso.iso_dr03_dz005_pt04=0.;
  iso.iso_dr03_dz005_pt05=0.;
  iso.iso_dr03_dz005_pt06=0.;
  iso.iso_dr03_dz005_pt07=0.;
  iso.iso_dr03_dz005_pt08=0.;
  iso.iso_dr03_dz005_pt09=0.;
  iso.iso_dr03_dz005_pt10=0.;
  


  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate
                                                                                                               
    if(cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals                                                                                                                          
    // when we want to find the second isolated muon and electron
    // we do not use the electron and muon in the isolation sum,                                                                                              
    // to avoid overlap with the other lepton in the event  
    if((abs(cms2.pfcands_particleId().at(thisPf))==13 || abs(cms2.pfcands_particleId().at(thisPf))==11) && abs(cms2.pfcands_particleId().at(ipf))==13) continue;                                                                                         
    if((abs(cms2.pfcands_particleId().at(thisPf))==13 || abs(cms2.pfcands_particleId().at(thisPf))==11) && abs(cms2.pfcands_particleId().at(ipf))==11) continue;                                                                                         

    //----------------------------------------                                                                                                                                                   
    // find closest PV and dz w.r.t. that PV                                                                                                                                                     
   //----------------------------------------                                                                                                                                                   

    float mindz = 999.;
    int vtxi    = -1;

    if (dovtxcut) {
      for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {


	int itrk = cms2.pfcands_trkidx().at(ipf);

	if( itrk >= (int)cms2.trks_trk_p4().size() || itrk < 0 ){
	  //note: this should only happen for electrons which do not have a matched track                                                                                                            
	  //currently we are just ignoring these guys                                                                                                                                                
	  continue;
	}

	////////

        if(!isGoodVertex(ivtx)) continue;

        float mydz = trks_dz_pv(itrk,ivtx).first;
	//        fillOverFlow( h_dz_vtx_trk , mydz );

        if (fabs(mydz) < fabs(mindz)) {
          mindz = mydz;
          vtxi = ivtx;
        }

      }

      //----------------------------------------------------------------------------                                                                                                               
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons                                                                                                                  
      //----------------------------------------------------------------------------                                                                                                               

      if ( vtxi != 0 )     continue;
    } else {
      //      mindz = trks_dz_pv(itrk,0).first;

      int itrk = -1;

      if (abs(cms2.pfcands_particleId().at(ipf))!=11) {
        itrk = cms2.pfcands_trkidx().at(ipf);
        if( itrk >= (int)cms2.trks_trk_p4().size() || itrk < 0 ) continue;
        mindz=trks_dz_pv(itrk,0).first;
      }

      if (abs(cms2.pfcands_particleId().at(ipf))==11 && cms2.pfcands_pfelsidx().at(ipf)>=0) {
        itrk = cms2.els_gsftrkidx().at(cms2.pfcands_pfelsidx().at(ipf));
        if( itrk >= (int)cms2.gsftrks_p4().size() || itrk < 0 ) continue;
        mindz=gsftrks_dz_pv(itrk,0).first;
      }


    }

    //---------------------------------------                                                                                                                                     
    // passes cuts, add up isolation value                                                                                                                                        
    //---------------------------------------                                                                                                                                                    

    coneR=0.3;    
    double dr=ROOT::Math::VectorUtil::DeltaR( cms2.pfcands_p4().at(ipf) , cms2.pfcands_p4().at(thisPf) );
    if( dr > coneR ) continue; // skip pfcands outside the cone                                     

    // this is the default
    if( cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt00+= cms2.pfcands_p4().at(ipf).pt();

    // cone 04
    //    if( pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt00+= pfcands_p4().at(ipf).pt();

    // veto Cone 0.05 // this hould be ok for the taus
    if(cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.05 && dr >= 0.05) iso.iso_dr0503_dz005_pt00 += cms2.pfcands_p4().at(ipf).pt();

    // veto Cone 0.015 // this should be ok for the electron in the endcap
    if(cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.05 && dr >= 0.015) iso.iso_dr01503_dz005_pt00 += cms2.pfcands_p4().at(ipf).pt();

    // this is the iso-sum option
    if( cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.05) iso.isoDir_dr03_dz005_pt00 +=cms2.pfcands_p4().at(ipf).pt()*(1-3*dr);

    // some dz variation
    if( cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.00) iso.iso_dr03_dz000_pt00+= cms2.pfcands_p4().at(ipf).pt(); 
    if( cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.10) iso.iso_dr03_dz010_pt00+= cms2.pfcands_p4().at(ipf).pt();     
    if( cms2.pfcands_p4().at(ipf).pt()>=0.0 && fabs(mindz) <= 0.20) iso.iso_dr03_dz020_pt00+= cms2.pfcands_p4().at(ipf).pt();
    
    // some pt variation
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.1 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt01+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.2 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt02+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.3 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt03+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.4 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt04+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.5 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt05+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.6 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt06+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.7 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt07+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.8 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt08+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 0.9 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt09+= cms2.pfcands_p4().at(ipf).pt();
    if(cms2.pfcands_p4().at(ipf).pt() >= 1.0 && fabs(mindz) <= 0.05) iso.iso_dr03_dz005_pt10+= cms2.pfcands_p4().at(ipf).pt();    
      
  } // end loop of the cand 

  return iso;
}


float trackIso( int thisPf , float coneR = 0.3 , float dz_thresh = 0.05 , bool dovtxcut = false , float minpt = 0.0 , float dRveto = 0.0 , float detaveto = 0.0 ){

  float iso = 0.0;

  // cout << "candidate: pt eta phi ID  " 
  //      << Form("%.1f",cms2.pfcands_p4().at(thisPf).pt())  << " " 
  //      << Form("%.2f",cms2.pfcands_p4().at(thisPf).eta()) << " " 
  //      << Form("%.2f",cms2.pfcands_p4().at(thisPf).phi()) << " "
  //      << cms2.pfcands_particleId().at(thisPf)            << " " 
  //      << endl;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                          ) continue; // skip this PFCandidate
    if( cms2.pfcands_charge().at(ipf) == 0     ) continue; // skip neutrals
    if( cms2.pfcands_p4().at(ipf).pt() < minpt ) continue; // minimum pt for tracks in trk iso

    if( dRbetweenVectors( cms2.pfcands_p4().at(ipf) , cms2.pfcands_p4().at(thisPf) ) > coneR    ) continue; // only consider tracks within dR < coneR
    if( dRbetweenVectors( cms2.pfcands_p4().at(ipf) , cms2.pfcands_p4().at(thisPf) ) < dRveto   ) continue; // veto tracks inside dRveto cone
    if( fabs( cms2.pfcands_p4().at(ipf).eta() - cms2.pfcands_p4().at(thisPf).eta() ) < detaveto ) continue; // veto tracks in phi strip

    int itrk = cms2.pfcands_trkidx().at(ipf);
    
    if( itrk >= (int)cms2.trks_trk_p4().size() || itrk < 0 ){
      //note: this should only happen for electrons which do not have a matched track
      //currently we are just ignoring these guys
      continue;
    }
    
    //----------------------------------------
    // find closest PV and dz w.r.t. that PV
    //----------------------------------------
    
    float mindz = 999.;
    int vtxi    = -1;
      
    if (dovtxcut) {
      for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
	
	if(!isGoodVertex(ivtx)) continue;
	
	float mydz = dz_trk_vtx(itrk,ivtx);
	
	if (fabs(mydz) < fabs(mindz)) {
	  mindz = mydz;
	  vtxi = ivtx;
	}
	
      }
      

    //----------------------------------------------------------------------------
    // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
    //----------------------------------------------------------------------------
    
      if ( vtxi != 0 )     continue;
    } else {
      mindz = dz_trk_vtx(itrk,0);
    }
    if ( fabs(mindz) > dz_thresh )     continue;

    //---------------------------------------
    // passes cuts, add up isolation value
    //---------------------------------------

    // cout << "track: pt eta phi ID dR    " 
    // 	 << Form("%.1f",cms2.pfcands_p4().at(ipf).pt())  << " " 
    // 	 << Form("%.2f",cms2.pfcands_p4().at(ipf).eta()) << " " 
    // 	 << Form("%.2f",cms2.pfcands_p4().at(ipf).phi()) << " "
    // 	 << cms2.pfcands_particleId().at(ipf)            << " " 
    // 	 << Form("%.3f",dRbetweenVectors( cms2.pfcands_p4().at(ipf) , cms2.pfcands_p4().at(thisPf))) << " "
    // 	 << endl;

    iso += cms2.pfcands_p4().at(ipf).pt();

  }

  //cout << "trackiso " << iso << endl;
  return iso;
}
/*
bool overlapMuon_ZMet2012_v1(int index , float ptcut = 10.0 ){

  for( unsigned int imu = 0 ; imu < cms2.mus_p4().size(); ++imu ){

    float dr = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4().at(index) , cms2.mus_p4().at(imu) );
    
    if( dr > 0.1                           ) continue;
    if( cms2.mus_p4().at(imu).pt() < ptcut ) continue;
    if( !muonId( imu , ZMet2012_v1 )       ) continue;
    
    return true;
  }

  return false;

}

bool passElectronSelection_Stop2012_v2(int index, bool vetoTransition, bool eta24, bool useOldIsolation ){

  if( vetoTransition && fabs(cms2.els_etaSC()[index]) > 1.4442 && fabs(cms2.els_etaSC()[index]) < 1.566 ) return false;
  if( eta24 && fabs(cms2.els_p4()[index].eta()) > 2.4 )                                                   return false;
  if( overlapMuon_ZMet2012_v1(index,10.0) )                                                               return false;

  electronIdComponent_t answer_loose_2012 = electronId_WP2012_v2(index, MEDIUM, useOldIsolation);
  if ((answer_loose_2012 & wp2012::PassAllWP2012Cuts) == PassAllWP2012Cuts)  return true;
  
  return false;
}
*/

int findTriggerIndex(TString trigName)
{
    vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
    vector<TString>::const_iterator end_it   = cms2.hlt_trigNames().end();
    vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
    if(found_it != end_it) return found_it - begin_it;
    return -1;
}

TString triggerName(TString triggerPattern){

  //-------------------------------------------------------
  // get exact trigger name corresponding to given pattern
  //-------------------------------------------------------

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < cms2.hlt_trigNames().size() ; ++itrig ){
    if( TString( cms2.hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = cms2.hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger) return "TRIGGER_NOT_FOUND";

  return exact_hltname;

}

int objectPassTrigger(const LorentzVector &obj, char* trigname, float ptmin = 0, float drmax = 0.1){

  TString exact_trigname = triggerName(trigname);

  if( exact_trigname.Contains("TRIGGER_NOT_FOUND") ){
    cout << "ERROR! no trigger matching pattern " << trigname << endl;
    return 0;
  }

  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(exact_trigname)];

  // cout << "Trigger pattern " << trigname       << endl;
  // cout << "Exact name      " << exact_trigname << endl;
  // cout << "Num objects     " << trigp4.size()  << endl;

  if( trigp4.size() == 0 ) return 0;

  for (unsigned int i = 0; i < trigp4.size(); ++i){
    if ( trigp4[i].Pt() < ptmin ) continue;
    float dr = dRbetweenVectors(trigp4[i], obj);
    if (dr < drmax) return cms2.hlt_prescales().at(findTriggerIndex(exact_trigname));
  }

  return 0;

}


LeptonTreeMaker::LeptonTreeMaker(bool lockToCoreSelectors, bool useLHeleId, bool useMVAeleId, bool doDYNNLOw, unsigned int prescale, bool realData) 
  : lockToCoreSelectors_(lockToCoreSelectors), useLHeleId_(useLHeleId), useMVAeleId_(useMVAeleId), 
    doDYNNLOw_(doDYNNLOw), prescale_(prescale), realData_(realData)
{
  std::cout << "LeptonTreeMaker::LeptonTreeMaker" << std::endl;
  std::cout << "\t lockToCoreSelectors: " << lockToCoreSelectors << std::endl;
  std::cout << "\t useLHeleId         : " << useLHeleId << std::endl;
  std::cout << "\t useMVAeleId        : " << useMVAeleId << std::endl;
  std::cout << "\t doDYNNLOw          : " << doDYNNLOw << std::endl;
  std::cout << "\t prescale           : " << prescale << std::endl;
  std::cout << "\t realData           : " << realData << std::endl;

  //
  // set up jet corrections
  //

  // new 53X
  char* dataJEC = (char*) "GR_P_V39_AN3";
  char* mcJEC   = (char*) "DESIGN53_V15";

  std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  try {
    jetcorr_filenames_pfL1FastJetL2L3.clear();
    if (realData) {
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L1FastJet_AK5PF.txt"    , dataJEC ));
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L2Relative_AK5PF.txt"   , dataJEC ));
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L3Absolute_AK5PF.txt"   , dataJEC ));
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L2L3Residual_AK5PF.txt" , dataJEC ));

    } else {
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L1FastJet_AK5PF.txt"  , mcJEC ));
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L2Relative_AK5PF.txt" , mcJEC ));
      jetcorr_filenames_pfL1FastJetL2L3.push_back  (Form("jetCorrections/%s_L3Absolute_AK5PF.txt" , mcJEC ));    

    }
    jet_corrector_pfL1FastJetL2L3_ = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);
  } catch (...){
    cout << "\nFailed to setup correctors needed to get Jet Enetry Scale. Abort\n" << endl;
    assert(0);
  }

  //
  // other parameters
  //

  batchMode_ = false;

}

LeptonTreeMaker::~LeptonTreeMaker() {
  std::cout << "LeptonTreeMaker::~LeptonTreeMaker" << std::endl;
  //if (electronIdMVA_ != 0) delete electronIdMVA_;
  //if (jet_corrector_pfL1FastJetL2L3_ != 0) delete jet_corrector_pfL1FastJetL2L3_;
}

void LeptonTreeMaker::ScanChain(TString outfileid,
				TChain* chain, 
				SmurfTree::DataType sample,
				double integratedLumi, // in unit of pb^-1, if negative the weight is 1.
				double xsec,           // in unit of pb, if negative take it from evt_xsec_excl*evt_kfactor
				int nProcessedEvents,  // if negative, take it from evt_nEvts
				bool identifyEvents, 
				bool realData,
				TString cms2_json_file)
{

  std::string prefix = SmurfTree::name(sample);
  std::cout << sample << ", " << prefix << std::endl;

  if ( chain->GetListOfFiles()->GetEntries()==0 ){
    printf("\nERROR: chain is empty for sample: %s\n\n",prefix.c_str());
    assert(0);
  }

  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  gErrorIgnoreLevel = -1;
  unsigned int nEventsTotal = 0;

  // make smurf ntuples
  gSystem->MakeDirectory("smurf");
  //gSystem->MakeDirectory("smurf/V00-00-00");
  TFile* fSmurf = TFile::Open(Form("smurf/%s_%s.root",prefix.c_str(), outfileid.Data()),"RECREATE");
  assert(fSmurf);
  LeptonTree leptonTree;
  leptonTree.CreateTree();
  leptonTree.tree_->SetDirectory(fSmurf);

  // triggers
  HLT_Ele27_WP80_tag_    = 0;
  HLT_Ele27_WP80_probe_  = 0;
  HLT_IsoMu24_tag_       = 0;
  HLT_IsoMu24_probe_     = 0;

  vtxweight_ = 0.0;

  nbl_ = 0;
  nbm_ = 0;

  tagiso_   = 0.0;
  probeiso_ = 0.0;

  // leptype_ = 0;
  // elid_    = 0;
  // eliso_   = 0;
  // muid_    = 0;
  // muiso_   = 0;

  muispf_    = 0;
  mustahits_ = 0;
  munsegs_   = 0;
  munpix_    = 0;
  munlayers_ = 0;
  mutrk_     = 0;
  muglb_     = 0;
  mud0_      = 0.0;
  mudz_      = 0.0;
  muchi2ndf_ = 0.0;

  probepassid_  = 0;
  probepassiso_ = 0;

  probepfpt_    = 0.0;

  leptonTree.tree_->Branch("muispf"	       , &muispf_		        ,"muispf/i"    );
  leptonTree.tree_->Branch("mustahits"	       , &mustahits_		        ,"mustahits/i" );
  leptonTree.tree_->Branch("munsegs"	       , &munsegs_		        ,"munsegs/i"   );
  leptonTree.tree_->Branch("munpix"	       , &munpix_		        ,"munpix/i"    );
  leptonTree.tree_->Branch("munlayers"	       , &munlayers_		        ,"munlayers/i" );
  leptonTree.tree_->Branch("mutrk"	       , &mutrk_		        ,"mutrk/i"     );
  leptonTree.tree_->Branch("muglb"	       , &muglb_		        ,"muglb/i"     );
  leptonTree.tree_->Branch("mud0"	       , &mud0_		        	,"mud0/F"      );
  leptonTree.tree_->Branch("mudz"	       , &mudz_		        	,"mudz/F"      );
  leptonTree.tree_->Branch("muchi2ndf"	       , &muchi2ndf_		        ,"muchi2ndf/F" );

  leptonTree.tree_->Branch("probepassid"       , &probepassid_		        ,"probepassid/i" );
  leptonTree.tree_->Branch("probepassiso"      , &probepassiso_		        ,"probepassiso/i");

  leptonTree.tree_->Branch("probepfpt"         , &probepfpt_		        ,"probepfpt/F");

  leptonTree.tree_->Branch("HLT_Ele27_WP80_tag"	       , &HLT_Ele27_WP80_tag_		        	,"HLT_Ele27_WP80_tag/i");
  leptonTree.tree_->Branch("HLT_Ele27_WP80_probe"      , &HLT_Ele27_WP80_probe_		                ,"HLT_Ele27_WP80_probe/i");

  leptonTree.tree_->Branch("HLT_IsoMu24_tag"	       , &HLT_IsoMu24_tag_			        ,"HLT_IsoMu24_tag/i");
  leptonTree.tree_->Branch("HLT_IsoMu24_probe"         , &HLT_IsoMu24_probe_		                ,"HLT_IsoMu24_probe/i");

  leptonTree.tree_->Branch("vtxweight"                 , &vtxweight_	                 	        ,"vtxweight/F");
  leptonTree.tree_->Branch("drprobe"                   , &drprobe_	                 	        ,"drprobe/F");
  leptonTree.tree_->Branch("tkisoold"                  , &tkiso_old_	                 	        ,"tkisoold/F");
  leptonTree.tree_->Branch("tkisonew"                  , &tkiso_new_	                 	        ,"tkisonew/F");
  leptonTree.tree_->Branch("tkisonewnoveto"            , &tkiso_new_noveto_	          	        ,"tkisonewnoveto/F");
  leptonTree.tree_->Branch("tkisonewpt1"               , &tkiso_new_pt1_	            	        ,"tkisonewpt1/F");
  leptonTree.tree_->Branch("tkisonewpt2"               , &tkiso_new_pt2_	            	        ,"tkisonewpt2/F");
  leptonTree.tree_->Branch("tkisonewpt3"               , &tkiso_new_pt3_	            	        ,"tkisonewpt3/F");
  leptonTree.tree_->Branch("tkisonewpt4"               , &tkiso_new_pt4_	            	        ,"tkisonewpt4/F");
  leptonTree.tree_->Branch("tkisonewpt5"               , &tkiso_new_pt5_	            	        ,"tkisonewpt5/F");
  leptonTree.tree_->Branch("isoch"                     , &isoch_	                	        ,"isoch/F");
  leptonTree.tree_->Branch("tkisov4"                   , &tkisov4_	                 	        ,"tkisov4/F");

  leptonTree.tree_->Branch("nbl"                       , &nbl_	            	                        ,"nbl/I");
  leptonTree.tree_->Branch("nbm"                       , &nbm_	            	                        ,"nbm/I");

  leptonTree.tree_->Branch("tagiso"                    , &tagiso_	            	                ,"tagiso/F");
  leptonTree.tree_->Branch("probeiso"                  , &probeiso_	            	                ,"probeiso/F");




  // leptonTree.tree_->Branch("leptype"                , &leptype_	                    	        ,"leptype/I");
  // leptonTree.tree_->Branch("elid"                   , &elid_	            	                ,"elid/I");
  // leptonTree.tree_->Branch("eliso"                  , &eliso_	            	                ,"eliso/I");
  // leptonTree.tree_->Branch("muid"                   , &muid_	            	                ,"muid/I");
  // leptonTree.tree_->Branch("muiso"                  , &muiso_	            	                ,"muiso/I");

  std::map<unsigned int, std::set<unsigned int> > runList;

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int nFailedIdentification = 0;
  int nFilteredOut = 0;

  int i_permille_old = 0;

  if(realData) {
    if(cms2_json_file=="") {  
      std::cout<<"\n WARNING: Running on the real Data, but not JSON file!"<<std::endl;
      // assert(0);
    } else {
      //set_goodrun_file_json(cms2_json_file);
      set_goodrun_file(cms2_json_file);

    }
  }

  set_vtxreweight_rootfile("vtxreweight_Summer12MC_PUS10_19fb_Zselection.root",true);

  // vector<TH2D*>     fDYNNLOKFactorHists;           //vector of hist for Drell-Yan NNLO Kfactor
  // if (doDYNNLOw_ && (sample==SmurfTree::dyee || sample==SmurfTree::dymm || sample==SmurfTree::dytt) ) {
  //     TFile *tmpFile = new TFile( "files/fewz_powheg_weights_stepwise_2011_fine7.root", "READ");
  //     const int nMassBins = 41;
  //     for(int i=0; i<nMassBins; i++){         
  //         TString hname = TString::Format("weight_%02d",i+1);
  //         TH2D *tmpHist = (TH2D*)tmpFile->Get(hname);
  //         tmpHist->SetDirectory(0);
  //         fDYNNLOKFactorHists.push_back(tmpHist);
  //     }
  //     tmpFile->Close();
  //     delete tmpFile;
  // }

  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
    //        printf("current file: %s (%s), %s\n", currentFile->GetName(), 
    // 	      currentFile->GetTitle(), currentFile->IsA()->GetName());

    // if( verbose ) cout << currentFile->GetTitle() << endl;
    // if( TString(currentFile->GetTitle()).Contains("SingleMu_Run2011A-05Aug2011-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_999999_8_skim.root")  ) continue;
    // if( TString(currentFile->GetTitle()).Contains("SingleMu_Run2011A-PromptReco-v6_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_172949_0_skim.root") ) continue;
    // if( TString(currentFile->GetTitle()).Contains("SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_176201_0_skim.root") ) continue;
    // if( TString(currentFile->GetTitle()).Contains("SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_176304_0_skim.root") ) continue;
    // if( TString(currentFile->GetTitle()).Contains("SingleMu_Run2011B-PromptReco-v1_AOD/V04-02-33/SingleLeptonAndTwoJets/merged_ntuple_176548_0_skim.root") ) continue;


    TFile *f = TFile::Open(currentFile->GetTitle()); 
    assert(f);
    TTree *tree = (TTree*)f->Get("Events");
    assert(tree);
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    cms2.Init(tree);  // set branch addresses for TTree tree
    TStopwatch t;

    //Event Loop
    unsigned int nEvents = tree->GetEntries();

    for( unsigned int event = 0; event < nEvents; ++event) {
      tree->LoadTree(event);
      cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
      if (cms2.evt_event() % prescale_ !=0) continue;

      // bool passevent = false;
      // if( cms2.evt_run() == 166438 && cms2.evt_lumiBlock() == 738 && cms2.evt_event() == 831532430 ) passevent = true;
      // if (!passevent) continue;

      // Select the good runs from the json file
      if(realData && cms2_json_file!="") {
	if( !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      }
      runList[cms2.evt_run()].insert(cms2.evt_lumiBlock());

      double weight = 1.0;
      if ( !realData && integratedLumi>0 ){
	double mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
	weight = integratedLumi * mcweight * (xsec>0?xsec:cms2.evt_xsec_excl()*cms2.evt_kfactor()*cms2.evt_filt_eff()) /
	  (nProcessedEvents>0?nProcessedEvents:cms2.evt_nEvts());
      }       
      ++nEventsTotal;
      //
      // basic event preselection
      //

      //if (cms2.trks_d0().size() == 0) continue;  // needed to get rid of back Monte Carlo events in CMSSW_2_X analysis
      if (cms2.hyp_type().size() == 0)       continue; // skip events without hypothesis
      if ( !cleaning_goodVertexApril2011() ) continue; // require >=1 good vertex
      //
      // duplicate removal
      //

      EventIdentifier id(cms2,realData); // <<<<--------CRASH HERE
      if (is_duplicate(id)) {
	duplicates_total_n++;
	if(!realData) duplicates_total_weight += cms2.evt_scale1fb();
	//cout << "Duplicate event found. Run: " << cms2.evt_run() << ", Event:" << cms2.evt_event() << ", Lumi: " << cms2.evt_lumiBlock() << endl;
	continue;
      }
      //
      // counter
      //

      if (!batchMode_) {
	int i_permille = (int)floor(10 * nEventsTotal / float(nEventsChain));
	if (i_permille != i_permille_old) {
	  // xterm magic from L. Vacavant and A. Cerri
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	  fflush(stdout);
	  i_permille_old = i_permille;
	}
      }

      //
      // filter out invalid contributions
      // to the sample being analyzed
      //

      if ( identifyEvents ){
	// check if we know what we are looking at
	if ( ! isIdentified(sample) ) nFailedIdentification++;

	// filter by process
	if ( ! filterByProcess(sample) ) {
	  nFilteredOut++;
	  continue;
	}
      }

      //
      // fill lepton tree
      //

      /*
      unsigned int eventSelection = 0;
      //
      // electron FR
      //
      TPMERegexp regexp_ele8	    (	"HLT_Ele8_CaloIdL_TrkIdVL_v(\\d+)", "o"	     );
      TPMERegexp regexp_ele8trkidiso  (	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)", "o"			);
      TPMERegexp regexp_ele8trkidiso_jet30	(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)", "o"	);
      TPMERegexp regexp_ele17trkidiso			(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)", "o"			);
      TPMERegexp regexp_ele17trkidiso_jet30	(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)", "o"	);

      unsigned int eleFRMask = LeptonTree::QCDFakeEle;
      if (passedTriggerRegExp(regexp_ele8) 				||            		
      passedTriggerRegExp(regexp_ele8trkidiso) 		||            
      passedTriggerRegExp(regexp_ele8trkidiso_jet30) 	||      
      passedTriggerRegExp(regexp_ele17trkidiso) 		||           
      passedTriggerRegExp(regexp_ele17trkidiso_jet30) )     	eventSelection |= LeptonTree::QCDFakeEle;
      if (eventSelection & eleFRMask)                 			MakeElectronFakeRateTree(leptonTree, weight, sample, eventSelection);

      //
      // muon FR
      //
      TPMERegexp regexp_mu8	("HLT_Mu8_v(\\d+)", "o"	);
      TPMERegexp regexp_mu17	("HLT_Mu17_v(\\d+)", "o");

      unsigned int muonFRMask = LeptonTree::QCDFakeMu;
      if (	passedTriggerRegExp(regexp_mu8) || 
      passedTriggerRegExp(regexp_mu17) ) 		eventSelection |= LeptonTree::QCDFakeMu;
      if (eventSelection & muonFRMask)               	MakeMuonFakeRateTree(leptonTree, weight, sample, eventSelection);
      */

      //-------------------------------
      // electron tag and probe
      //-------------------------------

      if (cms2.els_p4().size() > 1)    MakeElectronTagAndProbeTree(leptonTree, weight, sample);

      //-------------------------------
      // muon tag and probe tree
      //-------------------------------

      if (cms2.mus_p4().size() > 1)    MakeMuonTagAndProbeTree(leptonTree, weight, sample);

    }

    //
    // timing information
    //

    t.Stop();

    if( verbose ){
      printf("Finished processing file: %s\n",currentFile->GetTitle());
      printf("Real time: %u events / %f s = %e event/s\n", nEvents, 
	     t.RealTime(), nEvents / t.RealTime());
      printf("CPU time: %u events / %f s = %e event/s\n", nEvents, 
	     t.CpuTime(), nEvents / t.CpuTime());
      printf("Total duplicate count: %d.  Total weight %f\n",   
	     duplicates_total_n, duplicates_total_weight);
    }

    //
    // tidy up
    //

    delete tree;
    delete f;
  }

  //
  // print monitoring information
  //

  if ( nEventsChain != nEventsTotal ) {
    printf("ERROR: number of events from files (%d) is not equal to total number"
	   " of events (%d)\n", nEventsChain, nEventsTotal);
  }
  printf("Total number of skipped events due to bad identification: %d (%0.0f %%)\n",   
	 nFailedIdentification, nFailedIdentification*100.0/(nEventsChain+1e-5));
  printf("Total number of filtered out events: %d (%0.0f %%)\n",   
	 nFilteredOut, nFilteredOut*100.0/(nEventsChain+1e-5));

  //
  // write smurf tree
  // and tidy up
  //

  fSmurf->cd(); 
  leptonTree.tree_->Write();
  fSmurf->Close();

}


void LeptonTreeMaker::MakeElectronRecoTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample)
{
  /*
  // loop on tags
  for (int el = 0; el < cms2.els_p4().size(); ++el) {

  // loop on probes
  for (int sc = 0; sc < cms2.scs_energy().size(); ++sc) {

  // check probe and tag do not overlap
  // this is based on comparing the edm::Ref 
  // in the supercluster and electron supercluster
  if (cms2.scs_elsidx()[sc] == el) continue;

  // basic kinematics
  if (fabs(cms2.scs_eta()[sc]) > 2.5) continue;
  if (cms2.scs_energy()[sc]/cosh(cms2.scs_eta()[sc]) < 10.0) continue;

  // fill the leptonTree
  SetCommonTreeVariables(leptonTree, weight, sample);
  leptonTree.pt_          = cms2.scs_energy()[sc]/cosh(cms2.scs_eta()[sc]);
  leptonTree.eta_         = cms2.scs_eta()[sc];
  leptonTree.phi_         = cms2.scs_phi()[sc];
  leptonTree.q_           = 0;

  // away jet
  //const std::vector<JetPair> &jets = getJets(jetType(), cms2.els_p4()[el], cms2.scs_p4()[sc], 0, 5.0, true, false);
  //leptonTree.leadingAwayJetPt_ =  -999.;
  //if (jets.size() > 0) leptonTree.leadingAwayJetPt_ = jets[0].first.Pt();

  // tag and probe mass
  LorentzVector tag = cms2.els_p4()[el];
  LorentzVector probe = cms2.scs_p4()[sc];
  LorentzVector tp = tag + probe;
  leptonTree.tagAndProbeMass_ = tp.M();

  // lepton selection
  leptonTree.leptonSelection_                                       |= LeptonTree::PassEleSC;
  if (cms2.scs_elsidx()[sc] != -1) leptonTree.leptonSelection_      |= LeptonTree::PassEleReco;

  // event selection
  leptonTree.eventSelection_                                        |= LeptonTree::ZeeTagAndProbe;

  // fill it
  leptonTree.tree_->Fill();
  }

  }
  */
}


void LeptonTreeMaker::MakeElectronTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample)
{

  bool isData = ( sample == SmurfTree::data );

  //----------------------------------------------------
  // store all good leptons for jet-lepton cleaning
  //----------------------------------------------------

  VofP4 goodLeptons;
            
  for( unsigned int iel = 0 ; iel < cms2.els_p4().size(); ++iel ){
    if( cms2.els_p4().at(iel).pt() < 10 )                                         continue;
    if( !passElectronSelection_Stop2012_v3(iel,true,true,false) )                 continue; // stop 2012 ID/iso
    goodLeptons.push_back( cms2.els_p4().at(iel) );
  }

  for( unsigned int imu = 0 ; imu < cms2.mus_p4().size(); ++imu ){
    if( cms2.mus_p4().at(imu).pt() < 10 )      continue;
    if( !muonId( imu , ZMet2012_v1 ))          continue;
    goodLeptons.push_back( cms2.mus_p4().at(imu) );
  }  

  VofP4 jets;


  sort(jets.begin(), jets.end(), sortByPt);

  // loop on tags
  for (unsigned int tag = 0; tag < cms2.els_p4().size(); ++tag) {

    // check ID criteria on tag
    //if (!goodElectronIsolated(tag, useLHeleId_, useMVAeleId_, electronIdMVA_, lockToCoreSelectors_)) continue;
    //if (! goodElectronIsolated(tag, useLHeleId_, useMVAeleId_, egammaMvaEleEstimator_leptree, lockToCoreSelectors_)) continue;
    //if( !pass_electronSelection( tag , electronSelection_ssV5 , false , false ) ) continue;

    if( !passElectronSelection_Stop2012_v3( tag , true , true , false ) )   continue;
    if (cms2.els_p4()[tag].Pt() < 20.)                                      continue;
    if (fabs(cms2.els_etaSC()[tag]) > 2.4)                                  continue;
    
    // loop on probes
    for (unsigned int probe = 0; probe < cms2.els_p4().size(); ++probe) {

      // check probe and tag do not overlap
      if (tag == probe) continue;

      // basic probe denominator
      if (cms2.els_p4()[probe].Pt() < 5.0)      continue;
      if (fabs(cms2.els_etaSC()[probe]) > 2.5)  continue;

      // fill the tree - event general variables
      SetCommonTreeVariables(leptonTree, weight, sample);
      leptonTree.eventSelection_ = LeptonTree::ZeeTagAndProbe;
      //leptype_ = 0;

      // fill the tree - probe specific variables
      leptonTree.probe_           = cms2.els_p4()[probe];
      leptonTree.qProbe_          = cms2.els_charge()[probe];
      leptonTree.tag_             = cms2.els_p4()[tag];
      leptonTree.qTag_            = cms2.els_charge()[tag];

      leptonTree.tagAndProbeMass_ = (cms2.els_p4()[probe] + cms2.els_p4()[tag]).M();
           
      HLT_Ele27_WP80_tag_	  = isData ? cms2.els_HLT_Ele27_WP80()[tag]   : 1;
      HLT_Ele27_WP80_probe_       = isData ? cms2.els_HLT_Ele27_WP80()[probe] : 1;

      probeiso_                   = electronIsoValuePF2012_FastJetEffArea_v3( probe , 0.3 , 0 , false );
      tagiso_                     = electronIsoValuePF2012_FastJetEffArea_v3( tag   , 0.3 , 0 , false );

      // fill the tree - criteria the probe passed 
      // const std::vector<JetPair> &jets = getJets(jetType(), cms2.els_p4()[tag], cms2.els_p4()[probe], 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
      // if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
      // if (jets.size()>1)	leptonTree.jet2_ = jets.at(1).first;
      // if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;

      //----------------------------------------------------
      // count jets
      //----------------------------------------------------

      jets.clear();
      if( jets.size() > 0 ) cout << "ERROR! jets.size() " << jets.size() << endl;

      nbl_ = 0;
      nbm_ = 0;

      for (unsigned int ijet = 0 ; ijet < cms2.pfjets_p4().size() ; ijet++) {
    
	// get L1FastL2L3Residual total correction
	jet_corrector_pfL1FastJetL2L3_->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3_->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3_->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3_->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3_->getCorrection();

	LorentzVector vjet = corr * cms2.pfjets_p4().at(ijet);

	if(       vjet.pt()    < 30. )           continue;
	if( fabs( vjet.eta() ) > 2.5 )           continue;

	if( dRbetweenVectors( vjet , cms2.els_p4()[probe] ) < 0.4 ) continue;
	if( dRbetweenVectors( vjet , cms2.els_p4()[tag]   ) < 0.4 ) continue;

	// lepton-jet overlap removal
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
	          
	// PFJetID
	if( !passesPFJetID(ijet) ) continue;
    
	if( cms2.pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.244 ) nbl_++;
	if( cms2.pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.679 ) nbm_++;

	jets.push_back(vjet);
      }

      // store jets quantities
      if( jets.size() > 0 ) leptonTree.jet1_ = jets.at(0);
      if( jets.size() > 1 ) leptonTree.jet2_ = jets.at(1);
      if( jets.size() > 2 ) leptonTree.jet3_ = jets.at(2);
      if( jets.size() > 3 ) leptonTree.jet4_ = jets.at(3);

      leptonTree.njets_ = jets.size();
			
      // FO
      // if (ElectronFOV4(probe))
      // 	leptonTree.leptonSelection_     |= LeptonTree::PassEleFO;
			
      // ID
      //if (goodElectronWithoutIsolation(probe, useLHeleId_, useMVAeleId_, egammaMvaEleEstimator_leptree))
      //if( pass_electronSelection( probe , electronSelection_ssV5_noIso , false , false ) ) 
      if( passElectronSelection_Stop2012_v3_NoIso( probe , true , true , false ) )
	leptonTree.leptonSelection_     |= LeptonTree::PassEleID;

      // ISO
      //if (ww_elIso(probe))
      //if( pass_electronSelection( probe , electronSelection_ssV5_iso , false , false ) ) 
      if( passElectronSelection_Stop2012_v3_Iso( probe , true , true , false ) )
	leptonTree.leptonSelection_     |= LeptonTree::PassEleIso;

      vtxweight_ = vtxweight(isData);

      float mindr   = 100;
      int   pfindex = -1;

      // find pfcandidate index matching probe muon
      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( cms2.pfcands_p4().at(ipf).pt() < 5.0         ) continue;
	if( cms2.pfcands_charge().at(ipf) == 0           ) continue;
	
	float dr = dRbetweenVectors( cms2.pfcands_p4().at(ipf) , cms2.els_p4()[probe] );

	if( dr < mindr ){
	  mindr   = dr;
	  pfindex = ipf;
	}
      }

      drprobe_         = mindr;

      //cout << endl << endl;
      //cout << "pT(probe)  dR " << cms2.els_p4()[probe].pt() << " " << mindr << endl;

      if( pfindex > -1 ){
	tkiso_old_       = trackIso(pfindex, 0.3, 0.20, true  , 0.0 , 0.07 , 0.025 );
	tkiso_new_       = trackIso(pfindex, 0.3, 0.05, false , 0.0 , 0.07 , 0.025 );
	tkiso_new_noveto_= trackIso(pfindex, 0.3, 0.05, false , 0.0 , 0.0  , 0.0   );
	tkiso_new_pt1_   = trackIso(pfindex, 0.3, 0.05, false , 0.1 , 0.07 , 0.025 );
	tkiso_new_pt2_   = trackIso(pfindex, 0.3, 0.05, false , 0.2 , 0.07 , 0.025 );
	tkiso_new_pt3_   = trackIso(pfindex, 0.3, 0.05, false , 0.3 , 0.07 , 0.025 );
	tkiso_new_pt4_   = trackIso(pfindex, 0.3, 0.05, false , 0.4 , 0.07 , 0.025 );
	tkiso_new_pt5_   = trackIso(pfindex, 0.3, 0.05, false , 0.5 , 0.07 , 0.025 );

	struct myTrackIso myTrackIso=trackIso_struct(pfindex);
	tkisov4_         = myTrackIso.iso_dr03_dz010_pt00;
      }
      else{
	tkiso_old_       =  -1;
	tkiso_new_       =  -1;
	tkiso_new_noveto_=  -1;
	tkiso_new_pt1_   =  -1;
	tkiso_new_pt2_   =  -1;
	tkiso_new_pt3_   =  -1;
	tkiso_new_pt4_   =  -1;
	tkiso_new_pt5_   =  -1;
	tkisov4_         =  -1;
      }

      isoch_           = electronIsoValuePF( probe , 0 , 0.3, 99999. , 0.1 , 0.07 , 0.025 , 0.025 , 0 );

      //float isoch2 = cms2.els_iso03_pf_ch().at(probe);
      //cout << "charged pfiso " << isoch_ << endl;

      leptonTree.mt_   = sqrt( 2 * cms2.els_p4().at(tag).pt() * cms2.evt_pfmet() * ( 1 - cos( cms2.els_p4().at(tag).phi() - cms2.evt_pfmetPhi() ) ) );
      
      // fill it
      leptonTree.tree_->Fill();

    }

  }
}


void LeptonTreeMaker::MakeMuonTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample) {

  bool isData = ( sample == SmurfTree::data );

  //----------------------------------------------------
  // store all good leptons for jet-lepton cleaning
  //----------------------------------------------------

  VofP4 goodLeptons;
            
  for( unsigned int iel = 0 ; iel < cms2.els_p4().size(); ++iel ){
    if( cms2.els_p4().at(iel).pt() < 10 )                                         continue;
    if( !passElectronSelection_Stop2012_v3(iel,true,true,false) )                 continue;
    goodLeptons.push_back( cms2.els_p4().at(iel) );
  }

  for( unsigned int imu = 0 ; imu < cms2.mus_p4().size(); ++imu ){
    if( cms2.mus_p4().at(imu).pt() < 10 )      continue;
    if( !muonId( imu , ZMet2012_v1 ))          continue;
    goodLeptons.push_back( cms2.mus_p4().at(imu) );
  }  

  //----------------------------------------------------
  // count jets
  //----------------------------------------------------

  VofP4 jets;
	
  //std::vector<Int_t> nullMu;  // null identified muons // FIXME
  //std::vector<Int_t> nullEle; // null identified electrons  // FIXME

  //----------------------------------------------------
  // loop over tags
  //----------------------------------------------------

  for (unsigned int tag = 0; tag < cms2.mus_p4().size(); ++tag) {

    // check ID criteria on tag
    // APPLY FULL TAG CRITERIA: USING STOP MUON SELECTION
    // if (!goodMuonIsolated(tag, lockToCoreSelectors_, false, muonIdMVA_leptree, muonMVAEstimator_leptree,  nullMu, nullEle) ) continue;
    
    if( cms2.mus_p4()[tag].Pt() < 20.0)          continue; // pt cut
    if( fabs(cms2.mus_p4()[tag].Eta()) > 2.4)    continue; // eta cut
    if( !muonId( tag , ZMet2012_v1 )   )         continue; // full ID and iso

    // loop on probes
    for (unsigned int probe = 0; probe < cms2.mus_p4().size(); ++probe) {

      // check probe and tag do not overlap
      if (tag == probe) continue;

      // basic probe denominator
      if (cms2.mus_p4()[probe].Pt() < 5.0)                 continue; // pt cut
      if (fabs(cms2.mus_p4()[probe].Eta()) > 2.4)          continue; // eta cut
      //if (((cms2.mus_type()[probe]) & (1<<1)) == 0)        continue; // global muon

      SetCommonTreeVariables(leptonTree, weight, sample);

      muispf_    = cms2.mus_pid_PFMuon().at(probe);                               // pf muon
      muchi2ndf_ = cms2.mus_gfit_chi2().at(probe)/cms2.mus_gfit_ndof().at(probe); // chi2/ndf
      mustahits_ = cms2.mus_gfit_validSTAHits().at(probe);                        // standalone hits
      munsegs_   = cms2.mus_numberOfMatchedStations().at(probe);                  // muon segments
      mutrk_     = (((cms2.mus_type()[probe]) & (1<<2)) == 0) ? 0 : 1;            // tracker muon
      muglb_     = (((cms2.mus_type()[probe]) & (1<<1)) == 0) ? 0 : 1;            // tracker muon

      // require muon has track index, and at least 1 good vertex
      int trkidx = cms2.mus_trkidx().at(probe);
      int vtxidx = firstGoodVertex();

      //if ( trkidx < 0 ) continue;
      //if ( vtxidx < 0 ) continue;

      if( trkidx >= 0 ){
	munpix_    = cms2.trks_valid_pixelhits().at(trkidx);                        // pixel hits
	munlayers_ = cms2.trks_nlayers().at(trkidx);                                // nlayers
      }
      else{
	munpix_    = -1;
	munlayers_ = -1;
      }

      if( trkidx >=0 && vtxidx >= 0 ){
	mud0_      = fabs(trks_d0_pv(trkidx, vtxidx).first);                        // d0
	mudz_      = fabs(trks_dz_pv(trkidx, vtxidx).first);                        // dz
      }
      else{
	mud0_      = -1.0;
	mudz_      = -1.0;
      }

      probepfpt_ = -1;
      int ipf1 = cms2.mus_pfmusidx().at(probe);
      if( ipf1 >= 0 ) probepfpt_ = cms2.pfmus_p4().at(ipf1).pt();

      // fill the tree - event general variables
      leptonTree.eventSelection_ = LeptonTree::ZmmTagAndProbe;
      //leptype_ = 1;

      // fill the tree - probe specific variables
      leptonTree.probe_       = cms2.mus_p4()[probe];
      leptonTree.qProbe_      = cms2.mus_charge()[probe];
      leptonTree.tag_         = cms2.mus_p4()[tag];
      leptonTree.qTag_        = cms2.mus_charge()[tag];

      leptonTree.tagAndProbeMass_ = (cms2.mus_p4()[probe] + cms2.mus_p4()[tag]).M();

      // 2011: REWRITE MATCHING OBJECT TO MUON
      if( cms2.evt_run() >= 193806 ){
	HLT_IsoMu24_tag_        =	isData ? objectPassTrigger( cms2.mus_p4()[tag]   , (char*) "HLT_IsoMu24_v" ) : 1;
	HLT_IsoMu24_probe_ 	=	isData ? objectPassTrigger( cms2.mus_p4()[probe] , (char*) "HLT_IsoMu24_v" ) : 1;
      }
      else{
	HLT_IsoMu24_tag_        =	isData ? objectPassTrigger( cms2.mus_p4()[tag]   , (char*) "HLT_IsoMu24_eta2p1_v" ) : 1;
	HLT_IsoMu24_probe_ 	=	isData ? objectPassTrigger( cms2.mus_p4()[probe] , (char*) "HLT_IsoMu24_eta2p1_v" ) : 1;
      }

      probeiso_                 = muonIsoValuePF2012_deltaBeta(probe);
      tagiso_                   = muonIsoValuePF2012_deltaBeta(tag);

      //----------------------------------------------------
      // count jets
      //----------------------------------------------------

      jets.clear();
      if( jets.size() > 0 ) cout << "ERROR! jets.size() " << jets.size() << endl;

      nbl_ = 0;
      nbm_ = 0;

      for (unsigned int ijet = 0 ; ijet < cms2.pfjets_p4().size() ; ijet++) {
    
	// get L1FastL2L3Residual total correction
	jet_corrector_pfL1FastJetL2L3_->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3_->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3_->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3_->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3_->getCorrection();

	LorentzVector vjet = corr * cms2.pfjets_p4().at(ijet);

	if(       vjet.pt()    < 30. )           continue;
	if( fabs( vjet.eta() ) > 2.5 )           continue;

	if( dRbetweenVectors( vjet , cms2.mus_p4()[probe] ) < 0.4 ) continue;
	if( dRbetweenVectors( vjet , cms2.mus_p4()[tag]   ) < 0.4 ) continue;

	// lepton-jet overlap removal
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
	          
	// PFJetID
	if( !passesPFJetID(ijet) ) continue;

	if( cms2.pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.244 ) nbl_++;
	if( cms2.pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.679 ) nbm_++;

	jets.push_back(vjet);
      }

      // store jets quantities
      if( jets.size() > 0 ) leptonTree.jet1_ = jets.at(0);
      if( jets.size() > 1 ) leptonTree.jet2_ = jets.at(1);
      if( jets.size() > 2 ) leptonTree.jet3_ = jets.at(2);
      if( jets.size() > 3 ) leptonTree.jet4_ = jets.at(3);

      leptonTree.njets_ = jets.size();

      // // FO
      // if ( fakableMuon(probe, MuFOV2, muonMVAEstimator_leptree,  nullMu, nullEle) ) 
      // 	leptonTree.leptonSelection_     |= LeptonTree::PassMuFO;

      // ID
      //if ( goodMuonWithoutIsolation(probe,false, muonIdMVA_leptree) )

      probepassid_  = 0;
      probepassiso_ = 0;

      if( muonIdNotIsolated( probe , ZMet2012_v1 )){
	leptonTree.leptonSelection_     |= LeptonTree::PassMuID;
	probepassid_ = 1;
      }

      // ISO
      //if ( muonIsoValue(probe,false) < 0.15 )
      if( muonIsoValuePF2012_deltaBeta(probe) < 0.15 ){
	leptonTree.leptonSelection_     |= LeptonTree::PassMuIso;
	probepassiso_ = 1;
      }

      vtxweight_ = vtxweight(isData);

      float mindr   = 100;
      int   pfindex = -1;

      // find pfcandidate index matching probe muon
      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {
	if( cms2.pfcands_p4().at(ipf).pt() < 5.0         ) continue;
	if( cms2.pfcands_charge().at(ipf) == 0           ) continue;

	float dr = dRbetweenVectors( cms2.pfcands_p4().at(ipf) , cms2.mus_p4()[probe] );

	if( dr < mindr ){
	  mindr   = dr;
	  pfindex = ipf;
	}
      }

      drprobe_         = mindr;

      if( pfindex > -1 ){
	tkiso_old_       = trackIso(pfindex, 0.3, 0.20, true       );
	tkiso_new_       = trackIso(pfindex, 0.3, 0.05, false      );
	tkiso_new_noveto_= trackIso(pfindex, 0.3, 0.05, false      );
	tkiso_new_pt1_   = trackIso(pfindex, 0.3, 0.05, false , 0.1);
	tkiso_new_pt2_   = trackIso(pfindex, 0.3, 0.05, false , 0.2);
	tkiso_new_pt3_   = trackIso(pfindex, 0.3, 0.05, false , 0.3);
	tkiso_new_pt4_   = trackIso(pfindex, 0.3, 0.05, false , 0.4);
	tkiso_new_pt5_   = trackIso(pfindex, 0.3, 0.05, false , 0.5);

	struct myTrackIso myTrackIso=trackIso_struct(pfindex);
	tkisov4_         = myTrackIso.iso_dr03_dz010_pt00;
      }
      else{
	tkiso_old_       =  -1;
	tkiso_new_       =  -1;
	tkiso_new_noveto_=  -1;
	tkiso_new_pt1_   =  -1;
	tkiso_new_pt2_   =  -1;
	tkiso_new_pt3_   =  -1;
	tkiso_new_pt4_   =  -1;
	tkiso_new_pt5_   =  -1;
	tkisov4_         =  -1;
      }

      isoch_           = muonIsoValuePF( probe , 0 , 0.3, 1.0, 0.1, 0);

      leptonTree.mt_   = sqrt( 2 * cms2.mus_p4().at(tag).pt() * cms2.evt_pfmet() * ( 1 - cos( cms2.mus_p4().at(tag).phi() - cms2.evt_pfmetPhi() ) ) );

      // fill it
      leptonTree.tree_->Fill();
    }

  }
}

//
// fake rates
//

// void LeptonTreeMaker::MakeMuonFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection)
// {
//   std::vector<Int_t> nullMu; // null identified muons // FIXME
//   std::vector<Int_t> nullEle; // null identified electrons  // FIXME

//   //
//   // count FO
//   //

//   unsigned int fo = 0;
//   unsigned int nfo = 0;
//   for (unsigned int i = 0; i < cms2.mus_p4().size(); ++i) {
//     if ( fakableMuon(i, MuFOV2, muonMVAEstimator_leptree,  nullMu, nullEle) ) { ++nfo; fo = i; }
//   }

//   //
//   // require exactly one FO
//   //

//   if (nfo == 1) {

//     //
//     // check event properties
//     //

//     if (metValue() > 20.0) return;

//     // 
//     // fill the tree
//     //

//     SetCommonTreeVariables(leptonTree, weight, sample);

//     // fill the tree - probe specific variables
//     leptonTree.probe_       = cms2.mus_p4()[fo];
//     leptonTree.qProbe_      = cms2.mus_charge()[fo];

//     LorentzVector nullvector(0.0, 0.0, 0.0, 0.0);

//     leptonTree.eventSelection_ = eventSelection;

//     // FR triggers
//     //		if( passedTrigger("HLT_Mu8_v16") ) leptonTree.fr_Mu8_ = true;
//     //		if( passedTrigger("HLT_Mu17_v3") ) leptonTree.fr_Mu17_ = true;


//     leptonTree.leptonSelection_     |= LeptonTree::PassMuFO;
//     if (goodMuonWithoutIsolation(fo,false, muonIdMVA_leptree))
//       leptonTree.leptonSelection_     |= LeptonTree::PassMuID;
//     if (ww_muIso(fo, muonMVAEstimator_leptree, nullMu,  nullEle))
//       leptonTree.leptonSelection_     |= LeptonTree::PassMuIso;

//     const std::vector<JetPair> &jets = getJets(jetType(), cms2.mus_p4()[fo], nullvector, 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
//     if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
//     if (jets.size()>1)	leptonTree.jet2_ = jets.at(1).first;
//     if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;

//     leptonTree.tree_->Fill(); 

//   }

// }

// void LeptonTreeMaker::MakeElectronFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection)
// {

//   //
//   // count FO
//   //

//   // 1 and only 1 very loose
//   unsigned int fo = 0;
//   unsigned int nfo = 0;
//   for (unsigned int i = 0; i < cms2.els_p4().size(); ++i) {
//     if (ElectronFOV4(i)) { 	++nfo; fo = i; }
//   }

//   //
//   // require exactly one FO
//   //

//   if (nfo == 1) {

//     //
//     // check event properties
//     //

//     if (metValue() > 20.0) return;

//     // 
//     // fill the tree
//     //

//     SetCommonTreeVariables(leptonTree, weight, sample);

//     // fill the tree - probe specific variables
//     leptonTree.probe_       = cms2.els_p4()[fo];
//     leptonTree.qProbe_      = cms2.els_charge()[fo];
//     LorentzVector nullvector(0.0, 0.0, 0.0, 0.0);
//     leptonTree.eventSelection_ = eventSelection;
//     /*
//     // FR triggers
//     if( passedTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v2")	|| 
//     passedTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v3")	|| 
//     passedTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v4") )
//     leptonTree.fr_Ele8_CaloIdL_TrkIdVL_ = true;
//     if( passedTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12")	||
//     passedTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13")	||
//     passedTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14") )
//     leptonTree.fr_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_ = true;
//     if( passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3")	|| 
//     passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4") 	||
//     passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") )
//     leptonTree.fr_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_ = true;
//     if( passedTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v3")	||
//     passedTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4") 	||
//     passedTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5") )	
//     leptonTree.fr_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30 = true;
//     if( passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v3") ||
//     passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4") ||
//     passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5") ) 
//     leptonTree.fr_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30 = true;
//     */
//     const std::vector<JetPair> &jets = getJets(jetType(), cms2.els_p4()[fo], nullvector, 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
//     if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
//     if (jets.size()>1)  leptonTree.jet2_ = jets.at(1).first;
//     if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;
//     leptonTree.leptonSelection_     	|= LeptonTree::PassEleFO;
//     if (goodElectronWithoutIsolation(fo, useLHeleId_, useMVAeleId_, egammaMvaEleEstimator_leptree))
//       leptonTree.leptonSelection_     |= LeptonTree::PassEleID;
//     if (ww_elIso(fo))
//       leptonTree.leptonSelection_     |= LeptonTree::PassEleIso;
//     leptonTree.tree_->Fill();

//   }

// }


void LeptonTreeMaker::SetCommonTreeVariables(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample)
{

  leptonTree.InitVariables();
  leptonTree.run_         = cms2.evt_run();
  leptonTree.event_       = cms2.evt_event();
  leptonTree.lumi_        = cms2.evt_lumiBlock();
  leptonTree.nvtx_        = nGoodVertex();
  leptonTree.scale1fb_    = weight;
  leptonTree.met_         = cms2.evt_pfmet();

  //leptonTree.rho_         = cms2.evt_ww_rho_vor();

  muispf_       = -1;
  muchi2ndf_    = -1.0;
  mustahits_    = -1;
  munsegs_      = -1;
  munpix_       = -1;
  munlayers_    = -1;
  mutrk_        = -1;
  muglb_        = -1;
  mud0_         = -1.0;
  mudz_         = -1.0;
  probepfpt_    = -1;
  probepassid_  = -1;
  probepassiso_ = -1;

}

float LeptonTreeMaker::GetAwayJetPt(LorentzVector lep1, LorentzVector lep2)
{
  const std::vector<JetPair> &jets = getJets(jetType(), lep1, lep2, 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
  if (jets.size() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lep1,jets[0].first))>1.0) return jets[0].first.Pt();
  else return -999.9;
}


