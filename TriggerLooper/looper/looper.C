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
#include "../CORE/CMS2.h"
#include "../CORE/metSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/mcSelections.h"
#include "../CORE/muonSelections.h"
#include "../Tools/goodrun.cc"
#include "../CORE/utilities.cc"
#include "../CORE/ttbarSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/mcSUSYkfactor.h"
#include "../CORE/triggerSuperModel.h"
#include "../CORE/triggerUtils.h"
#include "../Tools/vtxreweight.cc"
#include "../Tools/msugraCrossSection.cc"

bool verbose        = false;
bool doTenPercent   = false;

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

  std::vector<int>           trigId = cms2.hlt_trigObjs_id()[findTriggerIndex(trigname)];
  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(trigname)];

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

TString triggerName(TString triggerPattern){

  //-------------------------------------------------------
  // get exact trigger name corresponding to given pattern
  //-------------------------------------------------------

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger) return "TRIGGER_NOT_FOUND";

  return exact_hltname;

}

//--------------------------------------------------------------------

bool passUnprescaledHLTTriggerPattern(const char* arg){

  //---------------------------------------------
  // Check if trigger is unprescaled and passes
  //---------------------------------------------

  TString HLTTriggerPattern(arg);
  TString HLTTrigger = triggerName( HLTTriggerPattern );

  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND")){
    return false;
  }
  return passUnprescaledHLTTrigger( HLTTrigger );
}

//--------------------------------------------------------------------

bool passHLTTriggerPattern(const char* arg){

  //---------------------------------------------
  // Check if trigger is unprescaled and passes
  //---------------------------------------------

  TString HLTTriggerPattern(arg);
  TString HLTTrigger = triggerName( HLTTriggerPattern );

  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND")){
    return false;
  }
  return passHLTTrigger( HLTTrigger );
}

//--------------------------------------------------------------------

bool passElDijetMHT(){



  return false;
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

bool passesCaloJetID (const LorentzVector &jetp4)
{
  int jet_idx = -1;
  double minDR = 999;

  for (unsigned int i = 0; i < cms2.jets_p4().size(); i++)
    {
      double deltaR = ROOT::Math::VectorUtil::DeltaR(jetp4, cms2.jets_p4()[i]);

      if (deltaR < minDR)
	{
	  minDR = deltaR;
	  jet_idx = i;
	}
    }

  if (jet_idx < 0)
    return false;

  if (cms2.jets_emFrac()[jet_idx] < 0.01 || cms2.jets_fHPD()[jet_idx] > 0.98 || cms2.jets_n90Hits()[jet_idx] < 2)
    return false;

  return true;
}

//--------------------------------------------------------------------

bool passesPFJetID(unsigned int pfJetIdx) {

  float pfjet_chf_  = cms2.pfjets_chargedHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  float pfjet_cef_  = cms2.pfjets_chargedEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  int   pfjet_cm_   = cms2.pfjets_chargedMultiplicity()[pfJetIdx];
  int   pfjet_mult_ = pfjet_cm_ + cms2.pfjets_neutralMultiplicity()[pfJetIdx] + cms2.pfjets_muonMultiplicity()[pfJetIdx];

  if (pfjet_nef_ >= 0.99)
    return false;
  if (pfjet_nhf_ >= 0.99)
    return false;
  if (pfjet_mult_ < 2)
    return false;

  if (fabs(cms2.pfjets_p4()[pfJetIdx].eta()) < 2.4)
    {
      if (pfjet_chf_ < 1e-6)
	return false;
      if (pfjet_cm_ < 1)
	return false;
      if (pfjet_cef_ >= 0.99)
	return false;
    }

  return true;
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

  eledijet_hltele_ = 0; 
  eletrijet_hltele_ = 0; 
  mudijet_hltmu_ = 0; 
  mutrijet_hltmu_ = 0; 
  pjet1_ = 0;
  pjet2_ = 0;
  pjet3_ = 0;
  pjet4_ = 0;
  cjet1_ = 0;
  cjet2_ = 0;
  cjet3_ = 0;
  cjet4_ = 0;
  pjet1_L1Fast_ = -999.;
  pjet2_L1Fast_ = -999.;
  pjet3_L1Fast_ = -999.;
  pjet4_L1Fast_ = -999.;
  pjet1_L2L3_ = -999.;
  pjet2_L2L3_ = -999.;
  pjet3_L2L3_ = -999.;
  pjet4_L2L3_ = -999.;
  eledijet_n82_ = -9;
  eletrijet_n82_ = -9;
  eledijet_n85_ = -9;
  eletrijet_n85_ = -9;
  mudijet_n83_ = -9;
  mutrijet_n83_ = -9;
  mudijet_n85_ = -9;
  mutrijet_n85_ = -9;
  eledijet_trigmindr_ejet_ = 999.;
  eletrijet_trigmindr_ejet_ = 999.;
  mudijet_trigmindr_mujet_ = 999.;
  mutrijet_trigmindr_mujet_ = 999.;
  eledijet_trigdr_pjet1_ = 999.;
  eledijet_trigdr_pjet2_ = 999.;
  eledijet_trigdr_pjet3_ = 999.;
  eledijet_trigdr_pjet4_ = 999.;
  eletrijet_trigdr_pjet1_ = 999.;
  eletrijet_trigdr_pjet2_ = 999.;
  eletrijet_trigdr_pjet3_ = 999.;
  eletrijet_trigdr_pjet4_ = 999.;
  mudijet_trigdr_pjet1_ = 999.;
  mudijet_trigdr_pjet2_ = 999.;
  mudijet_trigdr_pjet3_ = 999.;
  mudijet_trigdr_pjet4_ = 999.;
  mutrijet_trigdr_pjet1_ = 999.;
  mutrijet_trigdr_pjet2_ = 999.;
  mutrijet_trigdr_pjet3_ = 999.;
  mutrijet_trigdr_pjet4_ = 999.;
  run_     = -999;
  event_   = -999;
  lumi_    = -999;
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
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

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
      tree->LoadTree(z);

      cms2.GetEntry(z);

      InitBaby();

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset() << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      //if( evt_run() < 178420 || evt_run() > 180291 ) continue;
      //if( evt_run() < 179959 || evt_run() > 180291 ) continue;

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
	if( mus_p4().at(imu).pt() < 10 )           continue;
	if( !muonId( imu , OSGeneric_v3 ))         continue;
	goodLeptons.push_back( mus_p4().at(imu) );
	lepId.push_back( mus_charge().at(imu) * 13 );
	lepIndex.push_back(imu);
	ngoodmu_++;
	ngoodlep_++;
      }  

      //std::vector<int> mutrigId = cms2.hlt_trigObjs_id()[findTriggerIndex("HLT_IsoMu17_eta2p1_DiCentralPFJet25_v5")];
      //std::vector<int> eltrigId = cms2.hlt_trigObjs_id()[findTriggerIndex("HLT_Ele27_WP80_DiCentralPFJet25_v5")];
      // std::vector<int> eltrigId = cms2.hlt_trigObjs_id()[findTriggerIndex("HLT_Ele100_CaloIdVT_TrkIdT_v3")];

      // if( eltrigId.size()>0 ){
      // 	for( unsigned int i = 0 ; i < eltrigId.size() ; i++ ){
      // 	  cout << i << " " << eltrigId.at(i) << endl;
      // 	}
      // }

      //------------------------------------------------
      // trigger study: turn-on curve for jet triggers
      //------------------------------------------------

      /*
      trgjet_  = 0;
      passtrg_ = -1;

      if( doTriggerStudy ){

	// only consider muon events
	if( abs(id1_) == 11 ) continue;

	// run range for which HLT_Mu17_CentralJet30 is un-prescaled
	if( evt_run() < 160329 || evt_run() > 163261 ) continue;

	// require event passes HLT_Mu20_v1
	if( !passUnprescaledHLTTrigger("HLT_Mu20_v1") ) continue;
 
	// require found trigger object
	std::vector<LorentzVector> muonObj;
	muonObj = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu20_v1")];

	if( muonObj.size() == 0 ) continue;

	// require offline muon matched to muon trigger object (dR < 0.1)
	if( dRbetweenVectors( *lep1_ , muonObj.at(0) ) > 0.1 ) continue;

	// now get trigger objects for HLT_Mu20_CentralJet30
	std::vector<LorentzVector> muonJetObj;
	if( evt_run() >= 160329 && evt_run() <= 161176 )
	  muonJetObj = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu17_CentralJet30_v1")];
	else if( evt_run() >= 161210 && evt_run() <= 163261 )
	  muonJetObj = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu17_CentralJet30_v2")];

	// did event pass HLT_Mu17_CentralJet30?
	passtrg_ = 0;
	if( evt_run() >= 160329 && evt_run() <= 161176 )
	  passtrg_ = passUnprescaledHLTTrigger("HLT_Mu17_CentralJet30_v1") ? 1 : 0;
	else if( evt_run() >= 161210 && evt_run() <= 163261 )
	  passtrg_ = passUnprescaledHLTTrigger("HLT_Mu17_CentralJet30_v2") ? 1 : 0;

	// now find highest pt offline jet matched to HLT jet
	int offlinejet = -1;
	int njets      =  0;
	float maxpt    = -1;

        for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
          
          LorentzVector vjet = jets_p4().at(ijet) * jets_corL1FastL2L3().at(ijet);

	  if( dRbetweenVectors( vjet , *lep1_ ) < 0.4 ) continue;
	  //if( vjet.pt() < 20. )                                          continue;

	  // if event passed trigger, require jet is matched to HLT object
	  if( passtrg_ == 1 ){
	    
	    bool HLTmatch = false;

	    for( unsigned int ihlt = 0 ; ihlt < muonJetObj.size() ; ++ihlt ){
	      
	      // exclude HLT muon object
	      if( dRbetweenVectors( vjet , muonObj.at(0) ) < 0.1 )       continue;
	    
	      // dr match to HLT jet object
	      if( dRbetweenVectors( vjet , muonJetObj.at(ihlt) ) > 0.4 ) continue;

	      HLTmatch = true;
	    }

	    if( !HLTmatch ) continue;
	  }

	  njets++;

	  if( vjet.pt() > maxpt ){
	    maxpt      = vjet.pt();
	    offlinejet = ijet;
	  }

        }

	//if( passtrg_ == 0 ) cout << "FAIL TRIGGER" << endl;
	if( offlinejet < 0 ) continue;
	//if( njets > 1      ) continue;

	// now store offline jet p4 and whether muon-jet trigger passed
	trgjet_  = &(jets_p4().at(offlinejet) * jets_corL1FastL2L3().at(offlinejet));

	outTree->Fill();
	continue;

	// cout << "muonJetObj.size() " << muonJetObj.size() << endl;
	// cout << "muon trigger       : pt, eta, phi    " << muonObj.at(0).pt() << ", " << muonObj.at(0).eta() << ", " << muonObj.at(0).phi() << endl;

	// if( muonJetObj.size() > 1 ){
	//   cout << "muon-jet trigger 1 : pt, eta, phi    " << muonJetObj.at(0).pt() << ", " << muonJetObj.at(0).eta() << ", " << muonJetObj.at(0).phi() << endl;
	//   cout << "muon-jet trigger 2 : pt, eta, phi    " << muonJetObj.at(1).pt() << ", " << muonJetObj.at(1).eta() << ", " << muonJetObj.at(1).phi() << endl;

	//
	//   float dr2 = dRbetweenVectors( muonObj.at(0) , muonJetObj.at(1) );

	//   cout << "dr1 " << dr1 << " dr2 " << dr2 << endl;
	// }



	//cout << endl << endl;
	//cout << "trigger objects " << muonObj.size() << endl;
	//if( muonObj.size() > 0 ) cout << "pt, eta, phi    " << muonObj.at(0).pt() << ", " << muonObj.at(0).eta() << ", " << muonObj.at(0).phi() << endl;
	//cout << "muon " << lep1_->pt() << ", " << lep1_->eta() << ", " << lep1_->phi() << endl;

      }
      */

      //-------------------------------------
      // jet counting
      //-------------------------------------
      //Add 2 jet corrections per jet
      //Min dR between e at HLT and jet object at HLT
      VofP4 vpfjets_p4;

      njets_ = 0.;
      ht_    = 0.;

      int   imaxjet   = -1;
      float maxjetpt  = -1.;

      vector<int> goodjets;
      goodjets.clear();

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

	LorentzVector vjet      = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet);

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( !passesPFJetID(ijet) )     continue;
	if( fabs( vjet.eta() ) > 2.5 ) continue;
	if( vjet.pt() < 30 )           continue;

	njets_++;
	ht_ += vjet.pt();

	vpfjets_p4.push_back( vjet );
      }

      sort( vpfjets_p4.begin(), vpfjets_p4.end(), sortByPt);

      if( njets_ > 0 ) 	pjet1_ = &( vpfjets_p4.at(0) );
      if( njets_ > 1 ) 	pjet2_ = &( vpfjets_p4.at(1) );
      if( njets_ > 2 ) 	pjet3_ = &( vpfjets_p4.at(2) );
      if( njets_ > 3 ) 	pjet4_ = &( vpfjets_p4.at(3) );

      if( njets_ > 0 ) {
	int i_j1 = getJetIndex(vpfjets_p4.at(0));
	pjet1_L1Fast_ = pfjets_corL1FastL2L3().at(i_j1)/pfjets_cor().at(i_j1);
	pjet1_L2L3_ = pfjets_cor().at(i_j1);
      }
      if( njets_ > 1 ) {
	int i_j2 = getJetIndex(vpfjets_p4.at(1));
	pjet2_L1Fast_ = pfjets_corL1FastL2L3().at(i_j2)/pfjets_cor().at(i_j2);
	pjet2_L2L3_ = pfjets_cor().at(i_j2);
      }
      if( njets_ > 2 ) {
	int i_j3 = getJetIndex(vpfjets_p4.at(2));
	pjet3_L1Fast_ = pfjets_corL1FastL2L3().at(i_j3)/pfjets_cor().at(i_j3);
	pjet3_L2L3_ = pfjets_cor().at(i_j3);
      }
      if( njets_ > 3 ) {
	int i_j4 = getJetIndex(vpfjets_p4.at(3));
	pjet4_L1Fast_ = pfjets_corL1FastL2L3().at(i_j4)/pfjets_cor().at(i_j4);
	pjet4_L2L3_ = pfjets_cor().at(i_j4);
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

      strcpy(dataset_, cms2.evt_dataset().Data());  //dataset name
      run_          = evt_run();                    //run
      lumi_         = evt_lumiBlock();              //lumi
      event_        = evt_event();                  //event

      //----------------------------------------
      // triggers
      //----------------------------------------

      eledijetmht15_ = passUnprescaledHLTTriggerPattern("HLT_Ele27_WP80_DiCentralPFJet25_PFMHT15_v")            		? 1 : 0; // 178420-180291
      eledijetmht25_ = passUnprescaledHLTTriggerPattern("HLT_Ele32_WP80_DiCentralPFJet25_PFMHT25_v")            		? 1 : 0; // 178420-180291
      eledijet_hg_   = passHLTTriggerPattern("HLT_Ele27_WP80_DiCentralPFJet25_v")                               		? 1 : 0; // 178420-180291
      eledijet_      = passHLTTriggerPattern("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralPFJet30_v")   		? 1 : 0; // 178420-180291 
      eletrijet_     = passUnprescaledHLTTriggerPattern("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v")  	? 1 : 0; // 178420-180291 
      elequadjet_    = passUnprescaledHLTTriggerPattern("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_QuadCentralPFJet30_v") 	? 1 : 0; // 178420-180291 

      ele27dijet25_  = passHLTTriggerPattern("HLT_Ele27_CaloIdVT_TrkIdT_DiCentralPFJet25_v") ? 1 : 0;

      mudijetmht15_ = passUnprescaledHLTTriggerPattern("HLT_IsoMu17_eta2p1_DiCentralPFJet25_PFMHT15_v")         		? 1 : 0; // 178420-180291
      mudijetmht25_ = passUnprescaledHLTTriggerPattern("HLT_IsoMu17_eta2p1_DiCentralPFJet25_PFMHT25_v")         		? 1 : 0; // 178420-180291
      mudijet_      = passHLTTriggerPattern("HLT_IsoMu17_eta2p1_DiCentralPFJet25_v")                            		? 1 : 0; // 178420-180291
      mutrijet_     = passUnprescaledHLTTriggerPattern("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v")                           	? 1 : 0; // 178420-180291
      muquadjet_    = passUnprescaledHLTTriggerPattern("HLT_IsoMu17_eta2p1_QuadCentralPFJet30_v")                          	? 1 : 0; // 178420-180291

      ele8dijet30_  = passHLTTriggerPattern("HLT_Ele8_CaloIdT_TrkIdT_DiJet30_v") ? 1 : 0;

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
      bool passtrig = (eledijet_==1 || eletrijet_==1 || mudijet_==1 || mutrijet_==1) ? true: false;
      if ( passtrig ) {

	std::vector<LorentzVector> muDiJetObj;
	std::vector<LorentzVector> muTriJetObj;
	std::vector<LorentzVector> eleDiJetObj;
	std::vector<LorentzVector> eleTriJetObj;
	std::vector<int>  muDiJetId;
	std::vector<int>  muTriJetId;
	std::vector<int>  eleDiJetId;
	std::vector<int>  eleTriJetId;

	//Need to find the exact triggers used for different runs
	TString HLTTrigger_EleDiJet  = triggerName( "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralPFJet30_v" );
	TString HLTTrigger_EleTriJet = triggerName( "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v" );
	TString HLTTrigger_MuDiJet   = triggerName( "HLT_IsoMu17_eta2p1_DiCentralPFJet25_v" );
	TString HLTTrigger_MuTriJet  = triggerName( "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v" );
	//	TString triggerName(TString triggerPattern)
	eleDiJetObj  = cms2.hlt_trigObjs_p4()[findTriggerIndex(HLTTrigger_EleDiJet)];
	eleTriJetObj = cms2.hlt_trigObjs_p4()[findTriggerIndex(HLTTrigger_EleTriJet)];
	eleDiJetId   = cms2.hlt_trigObjs_id()[findTriggerIndex(HLTTrigger_EleDiJet)];
	eleTriJetId  = cms2.hlt_trigObjs_id()[findTriggerIndex(HLTTrigger_EleTriJet)];
	muDiJetObj  = cms2.hlt_trigObjs_p4()[findTriggerIndex(HLTTrigger_MuDiJet)];
	muTriJetObj = cms2.hlt_trigObjs_p4()[findTriggerIndex(HLTTrigger_MuTriJet)];
	muDiJetId   = cms2.hlt_trigObjs_id()[findTriggerIndex(HLTTrigger_MuDiJet)];
	muTriJetId  = cms2.hlt_trigObjs_id()[findTriggerIndex(HLTTrigger_MuTriJet)];

	// now find HLT jet matched to offline jet
	//loop over selected jets 
	for (unsigned int ijet = 0; ijet < vpfjets_p4.size(); ijet++) {
	  
	  LorentzVector vjet = vpfjets_p4.at(ijet);
	  // dijet triggers
	  if ( eledijet_==1 ) {
	    eledijet_trigmindr_ejet_ = getMinDR(82, 85, eleDiJetId, eleDiJetObj);
	    float minhltdr = 999.;
	    eledijet_n82_ = 0;
	    eledijet_n85_ = 0;
	    for( unsigned int ihlt = 0 ; ihlt < eleDiJetObj.size() ; ++ihlt ){
	      if ( eleDiJetId.at(ihlt)==82 ) {
		eledijet_n82_++;
		if (eledijet_n82_==1) eledijet_hltele_ = &( eleDiJetObj.at(ihlt) );
		if ( (eleDiJetObj.at(ihlt)).pt() > eledijet_hltele_->pt() ) 
		  eledijet_hltele_ = &( eleDiJetObj.at(ihlt) );
	      }
	      if ( eleDiJetId.at(ihlt)==85 ) eledijet_n85_++;
	      //Only look for jet objects
	      if ( eleDiJetId.at(ihlt)!=85 ) continue;
	      // dr match to HLT jet object
	      float hltdr = dRbetweenVectors( vjet , eleDiJetObj.at(ihlt) );
	      if ( minhltdr > hltdr ) minhltdr = hltdr;
	    } // end loop over trig objects

	    if (ijet==0) eledijet_trigdr_pjet1_ = minhltdr;
	    else if (ijet==1) eledijet_trigdr_pjet2_ = minhltdr;
	    else if (ijet==2) eledijet_trigdr_pjet3_ = minhltdr;
	    else if (ijet==3) eledijet_trigdr_pjet4_ = minhltdr;
	  } 
	  if ( mudijet_==1 ) {
	    mudijet_trigmindr_mujet_ = getMinDR(83, 85, muDiJetId, muDiJetObj);
	    float minhltdr = 999.;
	    mudijet_n83_ = 0;
	    mudijet_n85_ = 0;
	    for( unsigned int ihlt = 0 ; ihlt < muDiJetObj.size() ; ++ihlt ){
	      if ( muDiJetId.at(ihlt)==83 ) {
		mudijet_n83_++;
		if (mudijet_n83_==1) mudijet_hltmu_ = &( muDiJetObj.at(ihlt) );
		if ( muDiJetObj.at(ihlt).pt() > mudijet_hltmu_->pt() )
		  mudijet_hltmu_ = &( muDiJetObj.at(ihlt) );
	      }
	      if ( muDiJetId.at(ihlt)==85 ) mudijet_n85_++;
	      //Only look for jet objects
	      if ( muDiJetId.at(ihlt)!=85 ) continue;
	      // dr match to HLT jet object
	      float hltdr = dRbetweenVectors( vjet , muDiJetObj.at(ihlt) );
	      if ( minhltdr > hltdr ) minhltdr = hltdr;
	    } // end loop over trig objects

	    if (ijet==0) mudijet_trigdr_pjet1_ = minhltdr;
	    else if (ijet==1) mudijet_trigdr_pjet2_ = minhltdr;
	    else if (ijet==2) mudijet_trigdr_pjet3_ = minhltdr;
	    else if (ijet==3) mudijet_trigdr_pjet4_ = minhltdr;
	  }// end dijet triggers
	  
	  // trijet triggers
	  if ( eletrijet_==1 ) {
	    eletrijet_trigmindr_ejet_ = getMinDR(82, 85, eleTriJetId, eleTriJetObj);
	    float minhltdr = 999.;
	    eletrijet_n82_ = 0;
	    eletrijet_n85_ = 0;
	    for( unsigned int ihlt = 0 ; ihlt < eleTriJetObj.size() ; ++ihlt ){
	      if ( eleTriJetId.at(ihlt)==82 ) {
		eletrijet_n82_++;
		if (eletrijet_n82_==1) eletrijet_hltele_ = &( eleTriJetObj.at(ihlt) );
		if ( eleTriJetObj.at(ihlt).pt() > eletrijet_hltele_->pt() )
		  eletrijet_hltele_ = &( eleTriJetObj.at(ihlt) );
	      }
	      if ( eleTriJetId.at(ihlt)==85 ) eletrijet_n85_++;
	      //Only look for jet objects
	      if ( eleTriJetId.at(ihlt)!=85 ) continue;
	      // dr match to HLT jet object
	      float hltdr = dRbetweenVectors( vjet , eleTriJetObj.at(ihlt) );
	      if ( minhltdr > hltdr ) minhltdr = hltdr;
	    } // end loop over trig objects

	    if (ijet==0) eletrijet_trigdr_pjet1_ = minhltdr;
	    else if (ijet==1) eletrijet_trigdr_pjet2_ = minhltdr;
	    else if (ijet==2) eletrijet_trigdr_pjet3_ = minhltdr;
	    else if (ijet==3) eletrijet_trigdr_pjet4_ = minhltdr;
	  } 
	  if ( mutrijet_==1 ) {
	    mutrijet_trigmindr_mujet_ = getMinDR(83, 85, muTriJetId, muTriJetObj);
	    float minhltdr = 999.;
	    mutrijet_n83_ = 0;
	    mutrijet_n85_ = 0;
	    for( unsigned int ihlt = 0 ; ihlt < muTriJetObj.size() ; ++ihlt ){
	      if ( muTriJetId.at(ihlt)==83 ) {
		mutrijet_n83_++;
		if (mutrijet_n83_==1) mutrijet_hltmu_ = &( muTriJetObj.at(ihlt) );
		if ( muTriJetObj.at(ihlt).pt() > mutrijet_hltmu_->pt() )
		  mutrijet_hltmu_ = &( muTriJetObj.at(ihlt) );
	      }
	      if ( muTriJetId.at(ihlt)==85 ) mutrijet_n85_++;
	      //Only look for jet objects
	      if ( muTriJetId.at(ihlt)!=85 ) continue;
	      // dr match to HLT jet object
	      float hltdr = dRbetweenVectors( vjet , muTriJetObj.at(ihlt) );
	      if ( minhltdr > hltdr ) minhltdr = hltdr;
	    } // end loop over trig objects

	    if (ijet==0) mutrijet_trigdr_pjet1_ = minhltdr;
	    else if (ijet==1) mutrijet_trigdr_pjet2_ = minhltdr;
	    else if (ijet==2) mutrijet_trigdr_pjet3_ = minhltdr;
	    else if (ijet==3) mutrijet_trigdr_pjet4_ = minhltdr;
	  }// end trijet triggers
	} // end loop over jets
      }// end check for pass trigger

      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  if(g_createTree) closeTree();
  
  already_seen.clear();

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
  outTree->Branch("ncjets",          &ncjets_,           "ncjets/I");
  outTree->Branch("ht",              &ht_,               "ht/F");
  outTree->Branch("htc",             &htc_,              "htc/F");
  outTree->Branch("pfmet",           &pfmet_,            "pfmet/F");
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
  outTree->Branch("eledijet_hg",     &eledijet_hg_,      "eledijet_hg/I");
  outTree->Branch("eledijetmht15",   &eledijetmht15_,    "eledijetmht15/I");
  outTree->Branch("eledijetmht25",   &eledijetmht25_,    "eledijetmht25/I");
  outTree->Branch("eledijet",        &eledijet_,         "eledijet/I");
  outTree->Branch("ele27dijet25",    &ele27dijet25_,     "ele27dijet25/I");
  outTree->Branch("eletrijet",       &eletrijet_,        "eletrijet/I");
  outTree->Branch("elequadjet",      &elequadjet_,       "elequadjet/I");
  outTree->Branch("mudijet",         &mudijet_,          "mudijet/I");
  outTree->Branch("mudijetmht15",    &mudijetmht15_,     "mudijetmht15/I");
  outTree->Branch("mudijetmht25",    &mudijetmht25_,     "mudijetmht25/I");
  outTree->Branch("mutrijet",        &mutrijet_,         "mutrijet/I");
  outTree->Branch("muquadjet",       &muquadjet_,        "muquadjet/I");
  outTree->Branch("pjet1_L1Fast",    &pjet1_L1Fast_,     "pjet1_L1Fast/F");
  outTree->Branch("pjet2_L1Fast",    &pjet2_L1Fast_,     "pjet2_L1Fast/F");
  outTree->Branch("pjet3_L1Fast",    &pjet3_L1Fast_,     "pjet3_L1Fast/F");
  outTree->Branch("pjet4_L1Fast",    &pjet4_L1Fast_,     "pjet4_L1Fast/F");
  outTree->Branch("pjet1_L2L3",      &pjet1_L2L3_,       "pjet1_L2L3/F");
  outTree->Branch("pjet2_L2L3",      &pjet2_L2L3_,       "pjet2_L2L3/F");
  outTree->Branch("pjet3_L2L3",      &pjet3_L2L3_,       "pjet3_L2L3/F");
  outTree->Branch("pjet4_L2L3",      &pjet4_L2L3_,       "pjet4_L2L3/F");
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
  outTree->Branch("cjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet1_	);
  outTree->Branch("cjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet2_	);
  outTree->Branch("cjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet3_	);
  outTree->Branch("cjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet4_	);
  outTree->Branch("pjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet1_	);
  outTree->Branch("pjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet2_	);
  outTree->Branch("pjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet3_	);
  outTree->Branch("pjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet4_	);

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

int looper::getJetIndex(LorentzVector jet) {

  int matchindex = -1;
  for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

    LorentzVector vjet      = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet);
    if (dRbetweenVectors( vjet , jet ) > 0.4 ) continue;
    if ( abs(jet.pt()-vjet.pt()) > 5) continue;
    matchindex = ijet;
    break;
  }
  
  return matchindex;
}
