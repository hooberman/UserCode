#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/LorentzVector.h"
#include "histtools.h"
#include "looper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"

#include "../CORE/CMS2.cc"
#ifndef __CINT__
#include "../CORE/utilities.cc"
#include "../CORE/ssSelections.cc"
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

bool verbose        = false;
bool doTenPercent   = false;
bool debug          = false;

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

//--------------------------------------------------------------------

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

//--------------------------------------------------------------------

int findTriggerIndex(TString trigName)
{
    vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
    vector<TString>::const_iterator end_it = hlt_trigNames().end();
    vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
    if(found_it != end_it) return found_it - begin_it;
    return -1;
}

//--------------------------------------------------------------------

float getMinDeltaRBetweenObjects( TString trigname , int id1 , int id2 , bool verbose = false ){ 
  
  //-------------------------------------------------------------------------------------------
  // this function returns the minimum deltaR between 2 triggers objects with ID's id1 and id2
  //-------------------------------------------------------------------------------------------

  //cout << "Checking trigger " << trigname << " " << passHLTTrigger(trigname) << endl;

  // first, get p4 and ID vectors
  int trigindex = findTriggerIndex(trigname);

  if( trigindex < 0 ) return -1.0; //ERROR! didn't find this trigger

  std::vector<int>           trigId = cms2.hlt_trigObjs_id()[findTriggerIndex(trigname)];
  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(trigname)];

  //cout << "number of objects " << trigId.size() << endl;

  assert( trigId.size() == trigp4.size() );
  if( trigId.size() == 0 ) return -2.0;

  if( verbose ){
    cout << endl;
    cout << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
    cout << trigname << " pass? " << passHLTTrigger(trigname) << endl;
    
    cout << "|" << setw(12) << "index" << setw(4) 
	 << "|" << setw(12) << "ID"    << setw(4) 
	 << "|" << setw(12) << "pt"    << setw(4) 
	 << "|" << setw(12) << "eta"   << setw(4) 
	 << "|" << setw(12) << "phi"   << setw(4) 
	 << "|" << endl;

    for(unsigned int i = 0 ; i < trigId.size() ; ++i ){
      
      cout << "|" << setw(12) << i << setw(4) 
	   << "|" << setw(12) << trigId.at(i) << setw(4) 
	   << "|" << setw(12) << Form("%.2f",trigp4.at(i).pt())  << setw(4) 
	   << "|" << setw(12) << Form("%.2f",trigp4.at(i).eta()) << setw(4) 
	   << "|" << setw(12) << Form("%.2f",trigp4.at(i).phi()) << setw(4) 
	   << "|" << endl;
    }
  }


  // store p4's of objects with ID's id1 and id2  
  VofP4 obs1;
  VofP4 obs2;
  
  for (int i = 0; i < trigp4.size(); ++i){
    if( trigId.at(i) == id1 ) obs1.push_back( trigp4.at(i) );
    if( trigId.at(i) == id2 ) obs2.push_back( trigp4.at(i) );
  }

  // compute min deltaR between id1 and id2 objects
  double drmin = 100.0;
  int i1min = -1;
  int i2min = -1;

  for( unsigned int i1 = 0 ; i1 < obs1.size() ; i1++ ){
    for( unsigned int i2 = 0 ; i2 < obs2.size() ; i2++ ){
      float dr12 = dRbetweenVectors( obs1.at(i1) , obs2.at(i2) );
      if( dr12 < drmin ){
	drmin = dr12;
	i1min = i1;
	i2min = i2;
      }
    }
  }

  if( verbose ) cout << "i1min i2min dR " << i1min << " " << i2min << " " << drmin << endl;

  return drmin;
}



bool objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt) 
{

  float drMin = 999.99;
  for (size_t i = 0; i < trigObjs.size(); ++i)
    {
      if (trigObjs[i].Pt() < pt) continue;
      float dr = dRbetweenVectors(trigObjs[i], obj);
      if (dr < drMin) drMin = dr;
    }

  if (drMin < 0.1) return true;
  return false;

}

bool objectPassTrigger(const LorentzVector &obj, char* trigname, float ptmin, int id, float drmax){

  TString exactTriggerName = triggerName( trigname );

  if( !passHLTTrigger( exactTriggerName ) ) return false;

  if( exactTriggerName.Contains("TRIGGER_NOT_FOUND") ){
    cout << "Error in objectPassTrigger: couldn't find trigger matched to : " << trigname << endl;
    return false;
  }

  std::vector<int>           trigId = cms2.hlt_trigObjs_id()[findTriggerIndex(exactTriggerName)];
  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(exactTriggerName)];

  assert( trigId.size() == trigp4.size() );
  if( trigId.size() == 0 ) return false;

  float drMin = 999.99;

  for (int i = 0; i < trigp4.size(); ++i){
    if ( trigp4[i].Pt() < ptmin ) continue;
    if ( trigId[i] != id        ) continue;
    float dr = dRbetweenVectors(trigp4[i], obj);
    if (dr < drMin) drMin = dr;
  }

  if (drMin < drmax) return true;
  return false;

}


float getTriggerObjectPt(char* trigname, int id){

  int index = findTriggerIndex(trigname);

  if( index < 0 ){
    cout << "ERROR! can't find trigger " << trigname << " in run " << evt_run() << ", quitting" << endl;
    exit(0);
  }

  std::vector<int>           trigId = cms2.hlt_trigObjs_id()[findTriggerIndex(trigname)];
  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(trigname)];

  assert( trigId.size() == trigp4.size() );
  if( trigId.size() == 0 ) return false;

  float ptmax = -2;

  for (int i = 0; i < trigp4.size(); ++i){
    if ( trigId[i] != id        ) continue;
    if( trigp4.at(i).pt() > ptmax ) ptmax = trigp4.at(i).pt();
  }

  return ptmax;
}

//--------------------------------------------------------------------

looper::looper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
}

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

void looper::InitBaby(){

  memset(dataset_, '\0', 500);
  eledijet_hltele_		= 0; 
  eletrijet_hltele_		= 0; 
  mudijet_hltmu_		= 0; 
  mutrijet_hltmu_		= 0; 
  lep1_				= 0;
  lep2_				= 0;
  lep3_				= 0;
  lep4_				= 0;
  pjet1_			= 0;
  pjet2_			= 0;
  pjet3_			= 0;
  pjet4_			= 0;
  pjetnoiso1_			= 0;
  pjetnoiso2_			= 0;
  pjetnoiso3_			= 0;
  pjetnoiso4_			= 0;
  cjet1_			= 0;
  cjet2_			= 0;
  cjet3_			= 0;
  cjet4_			= 0;
  elnoiso1_                     = 0;
  elnoiso2_                     = 0;
  elnoiso3_                     = 0;
  elnoiso4_                     = 0;
  munoiso1_                     = 0;
  munoiso2_                     = 0;
  munoiso3_                     = 0;
  munoiso4_                     = 0;
  pjet1_res_			= -999.;
  pjet2_res_			= -999.;
  pjet3_res_			= -999.;
  pjet4_res_			= -999.;
  pjet1_L1Fast_			= -999.;
  pjet2_L1Fast_			= -999.;
  pjet3_L1Fast_			= -999.;
  pjet4_L1Fast_			= -999.;
  pjet1_L2L3_			= -999.;
  pjet2_L2L3_			= -999.;
  pjet3_L2L3_			= -999.;
  pjet4_L2L3_			= -999.;
  eledijet_n82_			= -9;
  eletrijet_n82_		= -9;
  eledijet_n85_			= -9;
  eletrijet_n85_		= -9;
  mudijet_n83_			= -9;
  mutrijet_n83_			= -9;
  mudijet_n85_			= -9;
  mutrijet_n85_			= -9;
  eledijet_trigmindr_ejet_	= 999.;
  eletrijet_trigmindr_ejet_	= 999.;
  mudijet_trigmindr_mujet_	= 999.;
  mutrijet_trigmindr_mujet_	= 999.;
  eledijet_trigdr_pjet1_	= 999.;
  eledijet_trigdr_pjet2_	= 999.;
  eledijet_trigdr_pjet3_	= 999.;
  eledijet_trigdr_pjet4_	= 999.;
  eletrijet_trigdr_pjet1_	= 999.;
  eletrijet_trigdr_pjet2_	= 999.;
  eletrijet_trigdr_pjet3_	= 999.;
  eletrijet_trigdr_pjet4_	= 999.;
  mudijet_trigdr_pjet1_		= 999.;
  mudijet_trigdr_pjet2_		= 999.;
  mudijet_trigdr_pjet3_		= 999.;
  mudijet_trigdr_pjet4_		= 999.;
  mutrijet_trigdr_pjet1_	= 999.;
  mutrijet_trigdr_pjet2_	= 999.;
  mutrijet_trigdr_pjet3_	= 999.;
  mutrijet_trigdr_pjet4_	= 999.;
  run_				= -999;
  event_			= -999;
  lumi_				= -999;
}

//--------------------------------------------------------------------

void looper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

int looper::ScanChain(TChain* chain, char *prefix){

  if( debug )  cout << __LINE__ << ": start ScanChain" << endl;

  //------------------------------------------------------------------------------------------------------
  // load here the on-the-fly corrections/uncertainties L1FastL2L3 (MC) and L1FastL2L3Residual (DATA)
  // corrections are stored in jet_corrected_pfL1FastJetL2L3
  // uncertainties are stored in pfUncertainty
  //------------------------------------------------------------------------------------------------------

  std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

  jetcorr_filenames_pfL1FastJetL2L3.clear();

  jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L1FastJet.txt");
  jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2Relative.txt");
  jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L3Absolute.txt");
  jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2L3Residual.txt");
  
  jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);


  set_goodrun_file( g_json );

  bool isData = true;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  BookHistos(prefix);
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  int nSkip_els_conv_dist = 0;

  float netot  = 0.;
  float nmtot  = 0.;
  float nepass = 0.;
  float nmpass = 0.;

  if(g_createTree) makeTree(prefix);

  bool hasJptBtagBranch = true;

  char* thisFile = "blah";

  if( debug )  cout << __LINE__ << ": begin file loop" << endl;

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    if( !f || f->IsZombie() ) {
      cout << "Skipping bad input file: " << currentFile->GetTitle() << endl;
      continue; //exit(1);                                                                                             
    }

    if( strcmp(thisFile,currentFile->GetTitle()) != 0 ){
      thisFile = (char*) currentFile->GetTitle();
      cout << thisFile << endl;
    }

    TTree *tree = (TTree*)f->Get("Events");

    //Matevz
    //TTreeCache::SetLearnEntries(100);
    //tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      if( doTenPercent ){
	if( !(nEventsTotal%10==0) ) continue;
      }

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //Matevz
      //tree->LoadTree(z);

      cms2.GetEntry(z);

      if( debug )  cout << __LINE__ << ": got entry" << endl;

      InitBaby();

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      if( !cleaning_goodDAVertexApril2011() )                        continue;
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------
 
      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
          skipEvent = true;
        }
      }
             
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }
   
      //---------------------------------------------
      // find leptons passing analysis selection
      //---------------------------------------------

      if( debug )  cout << __LINE__ << ": do lepton selection" << endl;

      VofP4 goodLeptons;
      vector<int> lepId;
      vector<int> lepIndex;

      ngoodlep_ = 0;
      ngoodel_  = 0;
      ngoodmu_  = 0;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                              continue;
	if( !pass_electronSelection( iel , electronSelection_ssV5 , false , false ) ) continue;
	goodLeptons.push_back( els_p4().at(iel) );
	lepId.push_back( els_charge().at(iel) * 11 );
	lepIndex.push_back(iel);
	ngoodel_++;
	ngoodlep_++;
      }

      nosel_ = 0;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 20 )                                                 continue;
	if( !pass_electronSelection( iel , electronSelection_el_OSV3 , false , false ) ) continue;
	nosel_++;
      }

          
      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 5 )            continue;
	if( !muonId( imu , OSGeneric_v3 ))         continue;
	goodLeptons.push_back( mus_p4().at(imu) );
	lepId.push_back( mus_charge().at(imu) * 13 );
	lepIndex.push_back(imu);
	ngoodmu_++;
	ngoodlep_++;
      }  

      sort( goodLeptons.begin(), goodLeptons.end(), sortByPt);

      if( ngoodlep_ > 0 ) 	lep1_ = &( goodLeptons.at(0) );
      if( ngoodlep_ > 1 ) 	lep2_ = &( goodLeptons.at(1) );
      if( ngoodlep_ > 2 ) 	lep3_ = &( goodLeptons.at(2) );
      if( ngoodlep_ > 3 ) 	lep4_ = &( goodLeptons.at(3) );

      dilmass_ = -1;

      if( ngoodlep_ > 1 ){
	dilmass_ = (*lep1_ + *lep2_).mass();	
      }

      VofP4 goodElectronsNoIso;
      vector<int> elnoisoIndex;

      nelnoiso_ = 0;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 80 )                                                       continue;
	if( !pass_electronSelection( iel , electronSelection_ssV5_noIso , false , false ) ) continue;
	goodElectronsNoIso.push_back( els_p4().at(iel) );
	elnoisoIndex.push_back(iel);
	nelnoiso_ ++;
      }

      VofP4 goodMuonsNoIso;
      vector<int> munoisoIndex;

      elnoiso1_mt_ = -1;

      if( nelnoiso_ > 0 ){
 	elnoiso1_           = &( goodElectronsNoIso.at(0) );
	elnoiso1_wp80_      = objectPassTrigger( *elnoiso1_ , "HLT_Ele27_WP80_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso1_top_       = objectPassTrigger( *elnoiso1_ , "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso1_iso_       = electronIsolation_rel_v1      ( elnoisoIndex.at(0) , true );
	elnoiso1_isofj_     = electronIsolation_rel_FastJet ( elnoisoIndex.at(0) , true );
	elnoiso1_isovtx_    = electronIsolation_cor_rel_v1  ( elnoisoIndex.at(0) , true );
	elnoiso1_isopf_     = electronIsoValuePF            ( elnoisoIndex.at(0) , 0    );
	elnoiso1_isopffj03_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(0) , 0.3 , 0);
	elnoiso1_isopffj04_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(0) , 0.4 , 0);
	elnoiso1_mt_        = sqrt( 2 * evt_pfmet() * (*elnoiso1_).pt() * ( 1 - cos( evt_pfmetPhi() - (*elnoiso1_).eta() ) ) );
	elnoiso1_d0pv_      = electron_d0PV_smurfV3( elnoisoIndex.at(0) );
	elnoiso1_d0bs_      = els_d0corr().at(elnoisoIndex.at(0));
      }
      if( nelnoiso_ > 1 ){
 	elnoiso2_           = &( goodElectronsNoIso.at(1) );
	elnoiso2_wp80_      = objectPassTrigger( *elnoiso2_ , "HLT_Ele27_WP80_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso2_top_       = objectPassTrigger( *elnoiso2_ , "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso2_iso_       = electronIsolation_rel_v1      ( elnoisoIndex.at(1) , true );
	elnoiso2_isofj_     = electronIsolation_rel_FastJet ( elnoisoIndex.at(1) , true );
	elnoiso2_isovtx_    = electronIsolation_cor_rel_v1  ( elnoisoIndex.at(1) , true );
	elnoiso2_isopf_     = electronIsoValuePF            ( elnoisoIndex.at(1) , 0    );
	elnoiso2_isopffj03_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(1) , 0.3 , 0);
	elnoiso2_isopffj04_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(1) , 0.4 , 0);
      }
      if( nelnoiso_ > 2 ){
 	elnoiso3_           = &( goodElectronsNoIso.at(2) );
	elnoiso3_wp80_      = objectPassTrigger( *elnoiso3_ , "HLT_Ele27_WP80_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso3_top_       = objectPassTrigger( *elnoiso3_ , "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso3_iso_       = electronIsolation_rel_v1      ( elnoisoIndex.at(2) , true );
	elnoiso3_isofj_     = electronIsolation_rel_FastJet ( elnoisoIndex.at(2) , true );
	elnoiso3_isovtx_    = electronIsolation_cor_rel_v1  ( elnoisoIndex.at(2) , true );
	elnoiso3_isopf_     = electronIsoValuePF            ( elnoisoIndex.at(2) , 0    );
	elnoiso3_isopffj03_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(2) , 0.3 , 0);
	elnoiso3_isopffj04_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(2) , 0.4 , 0);
      }
      if( nelnoiso_ > 3 ){
 	elnoiso4_           = &( goodElectronsNoIso.at(3) );
	elnoiso4_wp80_      = objectPassTrigger( *elnoiso4_ , "HLT_Ele27_WP80_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso4_top_       = objectPassTrigger( *elnoiso4_ , "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v" , 20.0 , 82 , 0.2 ) ? 1 : 0;
	elnoiso4_iso_       = electronIsolation_rel_v1      ( elnoisoIndex.at(3) , true );
	elnoiso4_isofj_     = electronIsolation_rel_FastJet ( elnoisoIndex.at(3) , true );
	elnoiso4_isovtx_    = electronIsolation_cor_rel_v1  ( elnoisoIndex.at(3) , true );
	elnoiso4_isopf_     = electronIsoValuePF            ( elnoisoIndex.at(3) , 0    );
	elnoiso4_isopffj03_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(3) , 0.3 , 0);
	elnoiso4_isopffj04_ = electronIsoValuePF2012_FastJetEffArea( elnoisoIndex.at(3) , 0.4 , 0);
      }
      
      nmunoiso_ = 0;
      
      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 8 )                                    continue;
	if( !muonIdNotIsolated( imu , OSGeneric_v3 ))                      continue;
	goodMuonsNoIso.push_back( mus_p4().at(imu) );
	munoisoIndex.push_back(imu);
	nmunoiso_ ++;
      }

      munoiso1_mt_ = -1;

      if( nmunoiso_ > 0 ){
 	munoiso1_        = &( goodMuonsNoIso.at(0) );
	munoiso1_mu24_   = objectPassTrigger( *munoiso1_ , "HLT_IsoMu24_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso1_mu30_   = objectPassTrigger( *munoiso1_ , "HLT_IsoMu30_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso1_iso_    = muonIsoValue         ( munoisoIndex.at(0) , false );
	munoiso1_isofj_  = muonIsoValue_FastJet ( munoisoIndex.at(0) , false );
	munoiso1_isovtx_ = muonCorIsoValue      ( munoisoIndex.at(0) , false );
	munoiso1_isopf_  = muonIsoValuePF       ( munoisoIndex.at(0) , 0     );
	munoiso1_mt_     = sqrt( 2 * evt_pfmet() * (*munoiso1_).pt() * ( 1 - cos( evt_pfmetPhi() - (*munoiso1_).eta() ) ) );
	munoiso1_d0pv_   = mud0PV_smurfV3(munoisoIndex.at(0));
	munoiso1_d0bs_   = mus_d0corr().at(munoisoIndex.at(0));
      }
      if( nmunoiso_ > 1 ){
 	munoiso2_        = &( goodMuonsNoIso.at(1) );
	munoiso2_mu24_   = objectPassTrigger( *munoiso2_ , "HLT_IsoMu24_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso2_mu30_   = objectPassTrigger( *munoiso2_ , "HLT_IsoMu30_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso2_iso_    = muonIsoValue         ( munoisoIndex.at(1) , false );
	munoiso2_isofj_  = muonIsoValue_FastJet ( munoisoIndex.at(1) , false );
	munoiso2_isovtx_ = muonCorIsoValue      ( munoisoIndex.at(1) , false );
	munoiso2_isopf_  = muonIsoValuePF       ( munoisoIndex.at(1) , 0     );
      }
      if( nmunoiso_ > 2 ){
 	munoiso3_        = &( goodMuonsNoIso.at(2) );
	munoiso3_mu24_   = objectPassTrigger( *munoiso3_ , "HLT_IsoMu24_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso3_mu30_   = objectPassTrigger( *munoiso3_ , "HLT_IsoMu30_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso3_iso_    = muonIsoValue         ( munoisoIndex.at(2) , false );
	munoiso3_isofj_  = muonIsoValue_FastJet ( munoisoIndex.at(2) , false );
	munoiso3_isovtx_ = muonCorIsoValue      ( munoisoIndex.at(2) , false );
	munoiso3_isopf_  = muonIsoValuePF       ( munoisoIndex.at(2) , 0     );
      }
      if( nmunoiso_ > 3 ){
 	munoiso4_        = &( goodMuonsNoIso.at(3) );
	munoiso4_mu24_   = objectPassTrigger( *munoiso4_ , "HLT_IsoMu24_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso4_mu30_   = objectPassTrigger( *munoiso4_ , "HLT_IsoMu30_eta2p1_v" , 20.0 , 83 , 0.2 ) ? 1 : 0;
	munoiso4_iso_    = muonIsoValue         ( munoisoIndex.at(3) , false );
	munoiso4_isofj_  = muonIsoValue_FastJet ( munoisoIndex.at(3) , false );
	munoiso4_isovtx_ = muonCorIsoValue      ( munoisoIndex.at(3) , false );
	munoiso4_isopf_  = muonIsoValuePF       ( munoisoIndex.at(3) , 0     );
      }


      if( debug )  cout << __LINE__ << ": do jet selection" << endl;

      //-------------------------------------
      // jet counting
      //-------------------------------------
      //Add 2 jet corrections per jet
      //Min dR between e at HLT and jet object at HLT
      VofP4 vpfjets_p4;
      vpfjets_p4.clear();

      VofP4 vpfjetsnoiso_p4;
      vpfjetsnoiso_p4.clear();

      njets_      = 0.;
      ht_         = 0.;
      njetsnoiso_ = 0.;
      htnoiso_    = 0.;

      int   imaxjet   = -1;
      float maxjetpt  = -1.;

      vector<int> goodjets;
      goodjets.clear();
      
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

	if( fabs( pfjets_p4().at(ijet).eta() ) > 5.0 ) continue;

	//---------------------------------------------------------------------------
	// get total correction: L1FastL2L3 for MC, L1FastL2L3Residual for data
	//---------------------------------------------------------------------------

	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();

	LorentzVector vjet   = corr * pfjets_p4().at(ijet);
      
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( !passesPFJetID(ijet) )     continue;
	if( fabs( vjet.eta() ) > 2.5 ) continue;
	if( vjet.pt() < 30 )           continue;
        if( !passesPFJetID(ijet) )     continue;

	njets_++;
	ht_ += vjet.pt();

	vpfjets_p4.push_back( vjet );

	//--------------------------------------------------------------------------
	// overlap removal with non-isolated leptons
	//--------------------------------------------------------------------------

	rejectJet = false;
	for( int ilep = 0 ; ilep < goodElectronsNoIso.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodElectronsNoIso.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	for( int ilep = 0 ; ilep < goodMuonsNoIso.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodMuonsNoIso.at(ilep) )     < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;

	njetsnoiso_++;
	htnoiso_ += vjet.pt();

	vpfjetsnoiso_p4.push_back( vjet );


      }
      
      sort( vpfjets_p4.begin(), vpfjets_p4.end(), sortByPt);

      if( njets_ > 0 ) 	pjet1_ = &( vpfjets_p4.at(0) );
      if( njets_ > 1 ) 	pjet2_ = &( vpfjets_p4.at(1) );
      if( njets_ > 2 ) 	pjet3_ = &( vpfjets_p4.at(2) );
      if( njets_ > 3 ) 	pjet4_ = &( vpfjets_p4.at(3) );

      sort( vpfjetsnoiso_p4.begin(), vpfjetsnoiso_p4.end(), sortByPt);

      if( njetsnoiso_ > 0 ) 	pjetnoiso1_ = &( vpfjetsnoiso_p4.at(0) );
      if( njetsnoiso_ > 1 ) 	pjetnoiso2_ = &( vpfjetsnoiso_p4.at(1) );
      if( njetsnoiso_ > 2 ) 	pjetnoiso3_ = &( vpfjetsnoiso_p4.at(2) );
      if( njetsnoiso_ > 3 ) 	pjetnoiso4_ = &( vpfjetsnoiso_p4.at(3) );
      
      if( njets_ > 0 ) {
	int i_j1 = getJetIndex(vpfjets_p4.at(0),jet_corrector_pfL1FastJetL2L3);
      
	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(i_j1)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(i_j1).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(i_j1).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	if( factors.size() == 4 ){
	  pjet1_L1Fast_ = factors.at(0);
	  pjet1_L2L3_   = factors.at(2) / factors.at(0);
	  pjet1_res_    = factors.at(3) / factors.at(2);
	}

      }
      
      if( njets_ > 1 ) {
	int i_j2 = getJetIndex(vpfjets_p4.at(1),jet_corrector_pfL1FastJetL2L3);

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(i_j2)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(i_j2).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(i_j2).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	if( factors.size() == 4 ){
	  pjet2_L1Fast_ = factors.at(0);
	  pjet2_L2L3_   = factors.at(2) / factors.at(0);
	  pjet2_res_    = factors.at(3) / factors.at(2);
	}

      }
      
      if( njets_ > 2 ) {
	int i_j3 = getJetIndex(vpfjets_p4.at(2),jet_corrector_pfL1FastJetL2L3);

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(i_j3)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(i_j3).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(i_j3).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	if( factors.size() == 4 ){
	  pjet3_L1Fast_ = factors.at(0);
	  pjet3_L2L3_   = factors.at(2) / factors.at(0);
	  pjet3_res_    = factors.at(3) / factors.at(2);
	}

      }
      
      if( njets_ > 3 ) {
	int i_j4 = getJetIndex(vpfjets_p4.at(3),jet_corrector_pfL1FastJetL2L3);

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(i_j4)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(i_j4).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(i_j4).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	if( factors.size() == 4 ){
	  pjet4_L1Fast_ = factors.at(0);
	  pjet4_L2L3_   = factors.at(2) / factors.at(0);
	  pjet4_res_    = factors.at(3) / factors.at(2);
	}

      }
      
      //------------------------------------------
      // count calojets
      //------------------------------------------

      VofP4 vcalojets_p4;
      ncjets_ = 0;
      htc_    = 0.0;

      for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
	
	LorentzVector vjet = jets_p4().at(ijet) * jets_corL1FastL2L3().at(ijet);

	if( !passesCaloJetID( vjet ) )         continue;	
	if( fabs( vjet.eta() ) > 2.5 )         continue;
	if( vjet.pt() < 20           )         continue;

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;

	htc_ += vjet.pt();
	ncjets_ ++;
	
	vcalojets_p4.push_back(vjet);
      }
      
      sort( vcalojets_p4.begin(), vcalojets_p4.end(), sortByPt);

      if( ncjets_ > 0 ) 	cjet1_ = &( vcalojets_p4.at(0) );
      if( ncjets_ > 1 ) 	cjet2_ = &( vcalojets_p4.at(1) );
      if( ncjets_ > 2 ) 	cjet3_ = &( vcalojets_p4.at(2) );
      if( ncjets_ > 3 ) 	cjet4_ = &( vcalojets_p4.at(3) );

      //----------------------------
      // MET flavors
      //----------------------------

      pfmet_    = evt_pfmet();
      pfmetphi_ = evt_pfmetPhi();
      pfsumet_  = evt_pfsumet();

      //----------------------------------------
      // nvertex variables
      //----------------------------------------

      nvtx_ = 0;
    
      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
	if(isGoodVertex(v)) ++nvtx_;
      }

      ndavtx_ = 0;
    
      for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
	if(isGoodDAVertex(v)) ++ndavtx_;
      }
      
      strcpy(dataset_, cms2.evt_dataset().at(0).Data());  //dataset name
      run_          = evt_run();                    //run
      lumi_         = evt_lumiBlock();              //lumi
      event_        = evt_event();                  //event

      //----------------------------------------
      // triggers
      //----------------------------------------

      //-------------------------------------------
      //  passTriggerPrescale functions returns:
      // -1: no matching trigger found
      //  0: trigger didn't pass
      //  1: trigger passed, un-prescaled
      //  N: trigger passed, prescale N
      //-------------------------------------------

      // top electron+jets triggers
      eltrijet_             = passTriggerPrescale("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v");
      eltrijetbackup_       = passTriggerPrescale("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet50_Jet40_Jet30_v");
      eldijet_              = passTriggerPrescale("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralPFJet30_v");
      eljet_                = passTriggerPrescale("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v");
      elnoisotrijet_        = passTriggerPrescale("HLT_Ele25_CaloIdVT_TrkIdT_TriCentralPFJet30_v");
      elnoisotrijetbackup_  = passTriggerPrescale("HLT_Ele25_CaloIdVT_TrkIdT_TriCentralPFJet50_Jet40_Jet30_v");

      // top muon+jets triggers
      mutrijet_             = passTriggerPrescale("HLT_Iso10Mu20_eta2p1_TriCentralPFJet30_v");
      mutrijetbackup_       = passTriggerPrescale("HLT_Iso10Mu20_eta2p1_CentralPFJet50_Jet40_Jet30_v");
      mudijet_              = passTriggerPrescale("HLT_Iso10Mu20_eta2p1_DiCentralPFJet30_v");
      mujet_                = passTriggerPrescale("HLT_Iso10Mu20_eta2p1_CentralPFJet30_v");
      munoisotrijet_        = passTriggerPrescale("HLT_Mu20_eta2p1_TriCentralPFJet30_v");
      munoisotrijetbackup_  = passTriggerPrescale("HLT_Mu20_eta2p1_CentralPFJet50_Jet40_Jet30_v");

      // non-isolated dilepton-HT triggers
      eeht175_              = passTriggerPrescale("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT175_v");
      eeht225_              = passTriggerPrescale("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT225_v");
      mmht175_              = passTriggerPrescale("HLT_DoubleMu8_Mass8_PFHT175_v");
      mmht225_              = passTriggerPrescale("HLT_DoubleMu8_Mass8_PFHT225_v");
      emht175_              = passTriggerPrescale("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT175_v");
      emht225_              = passTriggerPrescale("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT225_v");

      // isolated dilepton-HT triggers
      mmisoht175_           = passTriggerPrescale("HLT_DoubleRelIso1p0Mu5_Mass8_PFHT175_v");
      mmisoht225_           = passTriggerPrescale("HLT_DoubleRelIso1p0Mu5_Mass8_PFHT225_v");
      emisoht175_           = passTriggerPrescale("HLT_RelIso1p0Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT175_v");
      emisoht225_           = passTriggerPrescale("HLT_RelIso1p0Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT225_v");

      // isolated single muon triggers
      isomu20_              = passTriggerPrescale("HLT_IsoMu20_eta2p1_v");
      isomu24_              = passTriggerPrescale("HLT_IsoMu24_eta2p1_v");
      isomu30_              = passTriggerPrescale("HLT_IsoMu30_eta2p1_v");
      isomu34_              = passTriggerPrescale("HLT_IsoMu34_eta2p1_v");
      isomu40_              = passTriggerPrescale("HLT_IsoMu40_eta2p1_v");

      // non-isolated single muon triggers
      mu24_                 = passTriggerPrescale("HLT_Mu24_eta2p1_v");
      mu30_                 = passTriggerPrescale("HLT_Mu30_eta2p1_v");
      mu40_                 = passTriggerPrescale("HLT_Mu40_eta2p1_v");
      mu50_                 = passTriggerPrescale("HLT_Mu50_eta2p1_v");

      mu8_                  = passTriggerPrescale("HLT_Mu8_v");
      mu17_                 = passTriggerPrescale("HLT_Mu17_v");

      // single-electron triggers
      el27wp80_             = passTriggerPrescale("HLT_Ele27_WP80_v");
      el27wp70_             = passTriggerPrescale("HLT_Ele27_WP70_v");
      el27_                 = passTriggerPrescale("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
      el32_                 = passTriggerPrescale("HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

      // multi-jet triggers
      quadjet70_            = passTriggerPrescale("HLT_QuadJet70_v");
      quadjet80_            = passTriggerPrescale("HLT_QuadJet80_v");
      quadjet90_            = passTriggerPrescale("HLT_QuadJet90_v");

      // single photon triggers
      photon20_             = passTriggerPrescale("HLT_Photon20_CaloIdVL_IsoL_v");
      photon30_             = passTriggerPrescale("HLT_Photon30_CaloIdVL_IsoL_v");
      photon50_             = passTriggerPrescale("HLT_Photon50_CaloIdVL_IsoL_v");
      photon75_             = passTriggerPrescale("HLT_Photon75_CaloIdVL_IsoL_v");
      photon90_             = passTriggerPrescale("HLT_Photon90_CaloIdVL_IsoL_v");
      photon135_            = passTriggerPrescale("HLT_Photon135_v");
      photon150_            = passTriggerPrescale("HLT_Photon150_v");
      photon160_            = passTriggerPrescale("HLT_Photon160_v");

      // higgs single photon triggers
      hphoton22_            = passTriggerPrescale("HLT_Photon22_R9Id90_HE10_Iso40_v");
      hphoton36_            = passTriggerPrescale("HLT_Photon36_R9Id90_HE10_Iso40_v");
      hphoton50_            = passTriggerPrescale("HLT_Photon50_R9Id90_HE10_Iso40_v");
      hphoton75_            = passTriggerPrescale("HLT_Photon75_R9Id90_HE10_Iso40_v");
      hphoton90_            = passTriggerPrescale("HLT_Photon90_R9Id90_HE10_Iso40_v");

      // single electron utility triggers
      el8_                  = passTriggerPrescale("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
      el8jet30_             = passTriggerPrescale("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
      el17_                 = passTriggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
      el17jet30_            = passTriggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
      el8vl_                = passTriggerPrescale("HLT_Ele8_CaloIdL_CaloIsoVL_v");
      el17vl_               = passTriggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_v");
      el8noiso_             = passTriggerPrescale("HLT_Ele8_CaloIdT_TrkIdVL_v");
      el8noisojet30_        = passTriggerPrescale("HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v");
      el30noiso_            = passTriggerPrescale("HLT_Ele30_CaloIdVT_TrkIdT_v");

      // Higgs dilepton triggers
      ee_                   = passTriggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
      mmtrk_                = passTriggerPrescale("HLT_Mu17_TkMu8_v");
      mm_                   = passTriggerPrescale("HLT_Mu17_Mu8_v");
      em_                   = passTriggerPrescale("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
      me_                   = passTriggerPrescale("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

     
      mindrej_ = getMinDeltaRBetweenObjects( triggerName("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30") , 82 , 85 ); // 82 = electron  85 = jet 
      mindrmj_ = getMinDeltaRBetweenObjects( triggerName("HLT_IsoMu17_eta2p1_TriCentralPFJet30")                         , 83 , 85 ); // 83 = muon      85 = jet 
      
      // store pt of electron matched to electron-dijet trigger
      elptmatch_ = -1;

      if( evt_run() >= 178420 && evt_run() <= 179889 ){
	elptmatch_ = getTriggerObjectPt( "HLT_Ele27_WP80_DiCentralPFJet25_v4" , 82);
      }
      else if( evt_run() >= 179959 && evt_run() <= 180291 ){
	elptmatch_ = getTriggerObjectPt( "HLT_Ele27_WP80_DiCentralPFJet25_v5" , 82);
      }
      
      // now get trigger objects and store minimum jet dR

      outTree->Fill();
      
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  if(g_createTree) closeTree();
  
  //already_seen.clear();

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  return 0;

}

//--------------------------------------------------------------------
 
void looper::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()

//--------------------------------------------------------------------

void looper::makeTree(char *prefix ){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  outFile   = new TFile(Form("../output/%s/%s.root",g_version,prefix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in looper.h
  outTree->Branch("ele8dijet30",     &ele8dijet30_,      "ele8dijet30/I");
  outTree->Branch("mindrej",         &mindrej_,          "mindrej/F");
  outTree->Branch("nosel",           &nosel_,            "nosel/I");
  outTree->Branch("mindrmj",         &mindrmj_,          "mindrmj/F");
  outTree->Branch("njets",           &njets_,            "njets/I");
  outTree->Branch("njetsnoiso",      &njetsnoiso_,       "njetsnoiso/I");
  outTree->Branch("ncjets",          &ncjets_,           "ncjets/I");
  outTree->Branch("ht",              &ht_,               "ht/F");
  outTree->Branch("htnoiso",         &htnoiso_,          "htnoiso/F");
  outTree->Branch("htc",             &htc_,              "htc/F");
  outTree->Branch("pfmet",           &pfmet_,            "pfmet/F");
  outTree->Branch("mmht150",         &mmht150_,          "mmht150/I");
  outTree->Branch("eeht150",         &eeht150_,          "ht150/I");
  outTree->Branch("emht150",         &emht150_,          "emht150/I");
  outTree->Branch("pfmetphi",        &pfmetphi_,         "pfmetphi/F");
  outTree->Branch("elptmatch",       &elptmatch_,        "elptmatch/F");
  outTree->Branch("pfsumet",         &pfsumet_,          "pfsumet/F");
  outTree->Branch("dataset",         &dataset_,          "dataset[200]/C");
  outTree->Branch("run",             &run_,              "run/I");
  outTree->Branch("lumi",            &lumi_,             "lumi/I");
  outTree->Branch("event",           &event_,            "event/I");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("ngoodmu",         &ngoodmu_,          "ngoodmu/I");
  outTree->Branch("nvtx",            &nvtx_,             "nvtx/I");
  outTree->Branch("ndavtx",          &ndavtx_,           "ndavtx/I");
  outTree->Branch("dilmass",         &dilmass_,          "dilmass/F");

  // outTree->Branch("eledijet_hg",     &eledijet_hg_,      "eledijet_hg/I");
  // outTree->Branch("eledijetmht15",   &eledijetmht15_,    "eledijetmht15/I");
  // outTree->Branch("eledijetmht25",   &eledijetmht25_,    "eledijetmht25/I");
  // outTree->Branch("eledijet",        &eledijet_,         "eledijet/I");
  // outTree->Branch("ele27dijet25",    &ele27dijet25_,     "ele27dijet25/I");
  // outTree->Branch("eletrijet",       &eletrijet_,        "eletrijet/I");
  // outTree->Branch("elequadjet",      &elequadjet_,       "elequadjet/I");
  // outTree->Branch("mudijet",         &mudijet_,          "mudijet/I");
  // outTree->Branch("mudijetmht15",    &mudijetmht15_,     "mudijetmht15/I");
  // outTree->Branch("mudijetmht25",    &mudijetmht25_,     "mudijetmht25/I");
  // outTree->Branch("mutrijet",        &mutrijet_,         "mutrijet/I");
  // outTree->Branch("muquadjet",       &muquadjet_,        "muquadjet/I");

  outTree->Branch("pjet1_res",       &pjet1_res_,        "pjet1_res/F");
  outTree->Branch("pjet2_res",       &pjet2_res_,        "pjet2_res/F");
  outTree->Branch("pjet3_res",       &pjet3_res_,        "pjet3_res/F");
  outTree->Branch("pjet4_res",       &pjet4_res_,        "pjet4_res/F");
  outTree->Branch("pjet1_L1Fast",    &pjet1_L1Fast_,     "pjet1_L1Fast/F");
  outTree->Branch("pjet2_L1Fast",    &pjet2_L1Fast_,     "pjet2_L1Fast/F");
  outTree->Branch("pjet3_L1Fast",    &pjet3_L1Fast_,     "pjet3_L1Fast/F");
  outTree->Branch("pjet4_L1Fast",    &pjet4_L1Fast_,     "pjet4_L1Fast/F");
  outTree->Branch("pjet1_L2L3",      &pjet1_L2L3_,       "pjet1_L2L3/F");
  outTree->Branch("pjet2_L2L3",      &pjet2_L2L3_,       "pjet2_L2L3/F");
  outTree->Branch("pjet3_L2L3",      &pjet3_L2L3_,       "pjet3_L2L3/F");
  outTree->Branch("pjet4_L2L3",      &pjet4_L2L3_,       "pjet4_L2L3/F");

  /*
  outTree->Branch("eledijet_n82",    &eledijet_n82_,     "eledijet_n82/I");
  outTree->Branch("eletrijet_n82",   &eletrijet_n82_,    "eletrijet_n82/I");
  outTree->Branch("mudijet_n83",     &mudijet_n83_,      "mudijet_n83/I");
  outTree->Branch("mutrijet_n83",    &mutrijet_n83_,     "mutrijet_n83/I");
  outTree->Branch("eledijet_n85",    &eledijet_n85_,     "eledijet_n85/I");
  outTree->Branch("eletrijet_n85",   &eletrijet_n85_,    "eletrijet_n85/I");
  outTree->Branch("mudijet_n85",     &mudijet_n85_,      "mudijet_n85/I");
  outTree->Branch("mutrijet_n85",    &mutrijet_n85_,     "mutrijet_n85/I");
  outTree->Branch("eledijet_trigmindr_ejet",    &eledijet_trigmindr_ejet_,     "eledijet_trigmindr_ejet/F");
  outTree->Branch("eletrijet_trigmindr_ejet",   &eletrijet_trigmindr_ejet_,    "eletrijet_trigmindr_ejet/F");
  outTree->Branch("mudijet_trigmindr_mujet",     &mudijet_trigmindr_mujet_,      "mudijet_trigmindr_mujet/F");
  outTree->Branch("mutrijet_trigmindr_mujet",    &mutrijet_trigmindr_mujet_,     "mutrijet_trigmindr_mujet/F");
  outTree->Branch("eledijet_trigdr_pjet1",        &eledijet_trigdr_pjet1_,         "eledijet_trigdr_pjet1/F");
  outTree->Branch("eledijet_trigdr_pjet2",        &eledijet_trigdr_pjet2_,         "eledijet_trigdr_pjet2/F");
  outTree->Branch("eledijet_trigdr_pjet3",        &eledijet_trigdr_pjet3_,         "eledijet_trigdr_pjet3/F");
  outTree->Branch("eledijet_trigdr_pjet4",        &eledijet_trigdr_pjet4_,         "eledijet_trigdr_pjet4/F");
  outTree->Branch("eletrijet_trigdr_pjet1",       &eletrijet_trigdr_pjet1_,        "eletrijet_trigdr_pjet1/F");
  outTree->Branch("eletrijet_trigdr_pjet2",       &eletrijet_trigdr_pjet2_,        "eletrijet_trigdr_pjet2/F");
  outTree->Branch("eletrijet_trigdr_pjet3",       &eletrijet_trigdr_pjet3_,        "eletrijet_trigdr_pjet3/F");
  outTree->Branch("eletrijet_trigdr_pjet4",       &eletrijet_trigdr_pjet4_,        "eletrijet_trigdr_pjet4/F");
  outTree->Branch("mudijet_trigdr_pjet1",         &mudijet_trigdr_pjet1_,          "mudijet_trigdr_pjet1/F");
  outTree->Branch("mudijet_trigdr_pjet2",         &mudijet_trigdr_pjet2_,          "mudijet_trigdr_pjet2/F");
  outTree->Branch("mudijet_trigdr_pjet3",         &mudijet_trigdr_pjet3_,          "mudijet_trigdr_pjet3/F");
  outTree->Branch("mudijet_trigdr_pjet4",         &mudijet_trigdr_pjet4_,          "mudijet_trigdr_pjet4/F");
  outTree->Branch("mutrijet_trigdr_pjet1",        &mutrijet_trigdr_pjet1_,         "mutrijet_trigdr_pjet1/F");
  outTree->Branch("mutrijet_trigdr_pjet2",        &mutrijet_trigdr_pjet2_,         "mutrijet_trigdr_pjet2/F");
  outTree->Branch("mutrijet_trigdr_pjet3",        &mutrijet_trigdr_pjet3_,         "mutrijet_trigdr_pjet3/F");
  outTree->Branch("mutrijet_trigdr_pjet4",        &mutrijet_trigdr_pjet4_,         "mutrijet_trigdr_pjet4/F");
  outTree->Branch("eledijet_hltele"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &eledijet_hltele_	);
  outTree->Branch("eletrijet_hltele"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &eletrijet_hltele_);
  outTree->Branch("mudijet_hltmu"       , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mudijet_hltmu_	);
  outTree->Branch("mutrijet_hltmu"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mutrijet_hltmu_  );
  */
  
  outTree->Branch("cjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet1_	);
  outTree->Branch("cjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet2_	);
  outTree->Branch("cjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet3_	);
  outTree->Branch("cjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet4_	);

  outTree->Branch("pjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet1_	);
  outTree->Branch("pjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet2_	);
  outTree->Branch("pjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet3_	);
  outTree->Branch("pjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet4_	);

  outTree->Branch("pjetnoiso1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjetnoiso1_	);
  outTree->Branch("pjetnoiso2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjetnoiso2_	);
  outTree->Branch("pjetnoiso3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjetnoiso3_	);
  outTree->Branch("pjetnoiso4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjetnoiso4_	);

  outTree->Branch("lep1"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_	);
  outTree->Branch("lep2"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_	);
  outTree->Branch("lep3"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep3_	);
  outTree->Branch("lep4"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep4_	);

  outTree->Branch("elnoiso1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &elnoiso1_	);
  outTree->Branch("elnoiso2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &elnoiso2_	);
  outTree->Branch("elnoiso3" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &elnoiso3_	);
  outTree->Branch("elnoiso4" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &elnoiso4_	);

  outTree->Branch("munoiso1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &munoiso1_	);
  outTree->Branch("munoiso2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &munoiso2_	);
  outTree->Branch("munoiso3" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &munoiso3_	);
  outTree->Branch("munoiso4" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &munoiso4_	);

  outTree->Branch("elnoiso1mt"               , &elnoiso1_mt_             ,  "elnoiso1mt/F"            );             
  outTree->Branch("munoiso1mt"               , &munoiso1_mt_             ,  "munoiso1mt/F"            );             

  outTree->Branch("munoiso1d0pv"             , &munoiso1_d0pv_           ,  "munoiso1d0pv/F"          );             
  outTree->Branch("munoiso1d0bs"             , &munoiso1_d0bs_           ,  "munoiso1d0ns/F"          );             
  outTree->Branch("elnoiso1d0pv"             , &elnoiso1_d0pv_           ,  "elnoiso1d0pv/F"          );             
  outTree->Branch("elnoiso1d0bs"             , &elnoiso1_d0bs_           ,  "elnoiso1d0ns/F"          );             

  // top electron+jets triggers
  outTree->Branch("eltrijet"                 , &eltrijet_                ,  "eltrijet/I"              );             
  outTree->Branch("eltrijetbackup"           , &eltrijetbackup_          ,  "eltrijetbackup/I"        );             
  outTree->Branch("eldijet"                  , &eldijet_                 ,  "eldijet/I"               );             
  outTree->Branch("eljet"                    , &eljet_                   ,  "eljet/I"                 );             
  outTree->Branch("eltrijet"                 , &eltrijet_                ,  "eltrijet/I"              );             
  outTree->Branch("elnoisotrijet"            , &elnoisotrijet_           ,  "elnoisotrijet/I"         );             
  outTree->Branch("elnoisotrijetbackup"      , &elnoisotrijetbackup_     ,  "enoisoltrijetbackup/I"   );             
					       
  // top muon+jets triggers		       
  outTree->Branch("mutrijet"                 , &mutrijet_                ,  "mutrijet/I"              );             
  outTree->Branch("mutrijetbackup"           , &mutrijetbackup_          ,  "mutrijetbackup/I"        );             
  outTree->Branch("mudijet"                  , &mudijet_                 ,  "mudijet/I"               );             
  outTree->Branch("mujet"                    , &mujet_                   ,  "mujet/I"                 );             
  outTree->Branch("mutrijet"                 , &mutrijet_                ,  "mutrijet/I"              );             
  outTree->Branch("munoisotrijet"            , &munoisotrijet_           ,  "munoisotrijet/I"         );             
  outTree->Branch("munoisotrijetbackup"      , &munoisotrijetbackup_     ,  "enoisoltrijetbackup/I"   );             
					       
  // non-isolated dilepton-HT triggers	       
  outTree->Branch("eeht175"                  , &eeht175_                 ,  "eeht175/I"               );             
  outTree->Branch("eeht225"                  , &eeht225_                 ,  "eeht225/I"               );             
  outTree->Branch("mmht175"                  , &mmht175_                 ,  "mmht175/I"               );             
  outTree->Branch("mmht225"                  , &mmht225_                 ,  "mmht225/I"               );             
  outTree->Branch("emht175"                  , &emht175_                 ,  "emht175/I"               );             
  outTree->Branch("emht225"                  , &emht225_                 ,  "emht225/I"               );             

  // isolated dilepton-HT triggers
  outTree->Branch("mmisoht175"               , &mmisoht175_              ,  "mmisoht175/I"            );             
  outTree->Branch("mmisoht225"               , &mmisoht225_              ,  "mmisoht225/I"            );             
  outTree->Branch("emisoht175"               , &emisoht175_              ,  "emisoht175/I"            );             
  outTree->Branch("emisoht225"               , &emisoht225_              ,  "emisoht225/I"            );             
					       
  // isolated single muon triggers	       
  outTree->Branch("isomu20"                  , &isomu20_                 ,  "isomu20/I"               );             
  outTree->Branch("isomu24"                  , &isomu24_                 ,  "isomu24/I"               );             
  outTree->Branch("isomu30"                  , &isomu30_                 ,  "isomu30/I"               );             
  outTree->Branch("isomu34"                  , &isomu34_                 ,  "isomu34/I"               );             
  outTree->Branch("isomu40"                  , &isomu40_                 ,  "isomu40/I"               );             
					       
  // non-isolated single muon triggers
  outTree->Branch("mu24"                     , &mu24_                    ,  "mu24/I"                  );             
  outTree->Branch("mu30"                     , &mu30_                    ,  "mu30/I"                  );             
  outTree->Branch("mu40"                     , &mu40_                    ,  "mu40/I"                  );             
  outTree->Branch("mu50"                     , &mu50_                    ,  "mu50/I"                  );             

  outTree->Branch("mu8"                      , &mu8_                     ,  "mu8/I"                   );             
  outTree->Branch("mu17"                     , &mu17_                    ,  "mu17/I"                  );             
					       
  // single-electron triggers		       
  outTree->Branch("el27wp80"                 , &el27wp80_                ,  "el27wp80/I"              );             
  outTree->Branch("el27wp70"                 , &el27wp70_                ,  "el27wp"                  );
  outTree->Branch("el27"                     , &el27_                    ,  "el27/I"                  );             
  outTree->Branch("el32"                     , &el32_                    ,  "el32/I"                  );             
					       
  // multi-jet triggers
  outTree->Branch("quadjet70"                , &quadjet70_               ,  "quadjet70/I"             );             
  outTree->Branch("quadjet80"                , &quadjet80_               ,  "quadjet80/I"             );             
  outTree->Branch("quadjet90"                , &quadjet90_               ,  "quadjet90/I"             );             
					       
  // single photon triggers		       
  outTree->Branch("photon20"                 , &photon20_                ,  "photon20/I"              );             
  outTree->Branch("photon30"                 , &photon30_                ,  "photon30/I"              );             
  outTree->Branch("photon50"                 , &photon50_                ,  "photon50/I"              );             
  outTree->Branch("photon75"                 , &photon75_                ,  "photon75/I"              );             
  outTree->Branch("photon90"                 , &photon90_                ,  "photon90/I"              );             
  outTree->Branch("photon135"                , &photon135_               ,  "photon135/I"             );             
  outTree->Branch("photon150"                , &photon150_               ,  "photon150/I"             );             
  outTree->Branch("photon160"                , &photon160_               ,  "photon160/I"             );             
					       
  // higgs single photon triggers	       
  outTree->Branch("hphoton22"                , &hphoton22_               ,  "hphoton22/I"             );             
  outTree->Branch("hphoton36"                , &hphoton36_               ,  "hphoton36/I"             );             
  outTree->Branch("hphoton50"                , &hphoton50_               ,  "hphoton50/I"             );             
  outTree->Branch("hphoton75"                , &hphoton75_               ,  "hphoton75/I"             );             
  outTree->Branch("hphoton90"                , &hphoton90_               ,  "hphoton90/I"             );             
					       
  // single electron utility triggers	       
  outTree->Branch("el8"                      , &el8_                     ,  "el8/I"                   );             
  outTree->Branch("el8jet30"                 , &el8jet30_                ,  "el8jet30/I"              );             
  outTree->Branch("el17"                     , &el17_                    ,  "el17/I"                  );             
  outTree->Branch("el17jet30"                , &el17jet30_               ,  "el17jet30/I"             );             
  outTree->Branch("el8vl"                    , &el8vl_                   ,  "el8vl/I"                 );             
  outTree->Branch("el17vl"                   , &el17vl_                  ,  "el17vl/I"                );             
  outTree->Branch("el8noiso"                 , &el8noiso_                ,  "el8noiso/I"              );             
  outTree->Branch("el8noisojet30"            , &el8noisojet30_           ,  "el8noisojet30/I"         );             
  outTree->Branch("el30noiso"                , &el30noiso_               ,  "el30noiso/I"             );             

  // Higgs dilepton triggers		       
  outTree->Branch("ee"                       , &ee_                      ,  "ee/I"                    );             
  outTree->Branch("mm"                       , &mm_                      ,  "mm/I"                    );             
  outTree->Branch("em"                       , &em_                      ,  "em/I"                    );             
  outTree->Branch("me"                       , &me_                      ,  "me/I"                    );             
  outTree->Branch("mmtrk"                    , &mmtrk_                   ,  "mmtrk/I"                 );             
					       
  outTree->Branch("nelnoiso"                 , &nelnoiso_                ,  "nelnoiso/I"              );             
  outTree->Branch("nmunoiso"                 , &nmunoiso_                ,  "nmunoiso/I"              );             

  outTree->Branch("elnoiso1_wp80"            , &elnoiso1_wp80_           ,  "elnoiso1_wp80/I"         );             
  outTree->Branch("elnoiso1_top"             , &elnoiso1_top_            ,  "elnoiso1_top/I"          );             
  outTree->Branch("elnoiso1_iso"             , &elnoiso1_iso_            ,  "elnoiso1_iso/F"          );             
  outTree->Branch("elnoiso1_isofj"           , &elnoiso1_isofj_          ,  "elnoiso1_isofj/F"        );             
  outTree->Branch("elnoiso1_isovtx"          , &elnoiso1_isovtx_         ,  "elnoiso1_isovtx/F"       );             
  outTree->Branch("elnoiso1_isopf"           , &elnoiso1_isopf_          ,  "elnoiso1_isopf/F"        );             
  outTree->Branch("elnoiso1_isopffj03"       , &elnoiso1_isopffj03_      ,  "elnoiso1_isopffj03/F"    );             
  outTree->Branch("elnoiso1_isopffj04"       , &elnoiso1_isopffj04_      ,  "elnoiso1_isopffj04/F"    );             
					       
  outTree->Branch("elnoiso2_wp80"            , &elnoiso2_wp80_           ,  "elnoiso2_wp80/I"         );             
  outTree->Branch("elnoiso2_top"             , &elnoiso2_top_            ,  "elnoiso2_top/I"          );             
  outTree->Branch("elnoiso2_iso"             , &elnoiso2_iso_            ,  "elnoiso2_iso/F"          );             
  outTree->Branch("elnoiso2_isofj"           , &elnoiso2_isofj_          ,  "elnoiso2_isofj/F"        );             
  outTree->Branch("elnoiso2_isovtx"          , &elnoiso2_isovtx_         ,  "elnoiso2_isovtx/F"       );             
  outTree->Branch("elnoiso2_isopf"           , &elnoiso2_isopf_          ,  "elnoiso2_isopf/F"        );             
  outTree->Branch("elnoiso2_isopffj03"       , &elnoiso2_isopffj03_      ,  "elnoiso2_isopffj03/F"    );             
  outTree->Branch("elnoiso2_isopffj04"       , &elnoiso2_isopffj04_      ,  "elnoiso2_isopffj04/F"    );             
					       
  outTree->Branch("munoiso1_mu24"            , &munoiso1_mu24_           ,  "munoiso1_mu24/I"         );             
  outTree->Branch("munoiso1_mu30"            , &munoiso1_mu30_           ,  "munoiso1_mu30/I"         );             
  outTree->Branch("munoiso1_iso"             , &munoiso1_iso_            ,  "munoiso1_iso/F"          );             
  outTree->Branch("munoiso1_isofj"           , &munoiso1_isofj_          ,  "munoiso1_isofj/F"        );             
  outTree->Branch("munoiso1_isovtx"          , &munoiso1_isovtx_         ,  "munoiso1_isovtx/F"       );             
  outTree->Branch("munoiso1_isopf"           , &munoiso1_isopf_          ,  "munoiso1_isopf/F"        );             
					       
  outTree->Branch("munoiso2_mu24"            , &munoiso2_mu24_           ,  "munoiso2_mu24/I"         );             
  outTree->Branch("munoiso2_mu30"            , &munoiso2_mu30_           ,  "munoiso2_mu30/I"         );             
  outTree->Branch("munoiso2_iso"             , &munoiso2_iso_            ,  "munoiso2_iso/F"          );             
  outTree->Branch("munoiso2_isofj"           , &munoiso2_isofj_          ,  "munoiso2_isofj/F"        );             
  outTree->Branch("munoiso2_isovtx"          , &munoiso2_isovtx_         ,  "munoiso2_isovtx/F"       );             
  outTree->Branch("munoiso2_isopf"           , &munoiso2_isopf_          ,  "munoiso2_isopf/F"        );             

  
}

//--------------------------------------------------------------------

vector<int> goodDAVertices(){

  vector<int> myVertices;
  myVertices.clear();
  
  for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
    if( !isGoodDAVertex(v) ) continue;
    myVertices.push_back(v);
  }
  
  return myVertices;
}

//--------------------------------------------------------------------

float looper::dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}


//--------------------------------------------------------------------

float looper::getMinDR(int type1, int type2, std::vector<int> hltid, std::vector<LorentzVector> hltobj) {
  
  float mindr = 9999.;
  for( unsigned int ihlt1 = 0 ; ihlt1 < hltobj.size() ; ++ihlt1 ){
    if ( hltid.at(ihlt1)!=type1 ) continue;
    for( unsigned int ihlt2 = 0 ; ihlt2 < hltobj.size() ; ++ihlt2 ){
      if ( hltid.at(ihlt2)!=type2 ) continue;
      float dr = dRbetweenVectors( hltobj.at(ihlt1) , hltobj.at(ihlt2) );
      if ( mindr > dr ) mindr = dr;
      }
  }
  return mindr;
}

//--------------------------------------------------------------------

int looper::getJetIndex(LorentzVector jet, FactorizedJetCorrector* jet_corrector) {

  int matchindex = -1;
  for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

    jet_corrector->setRho   ( cms2.evt_ww_rho_vor()           );
    jet_corrector->setJetA  ( cms2.pfjets_area().at(ijet)     );
    jet_corrector->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
    jet_corrector->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
    double corr = jet_corrector->getCorrection();
    
    LorentzVector vjet   = corr * pfjets_p4().at(ijet);

    //LorentzVector vjet      = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet);
    if ( dRbetweenVectors( vjet , jet ) > 0.001 ) continue;
    if ( abs(jet.pt()-vjet.pt())        > 0.001)  continue;

    matchindex = ijet;
    break;
  }

  if( matchindex < 0 ){
    cout << __FILE__ << " " << __LINE__ << " : ERROR! can't find matched jet, quitting" << endl;
    exit(0);
  }
  
  return matchindex;
}
