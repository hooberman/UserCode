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

#include "../CORE/CMS2.cc"
#ifndef __CINT__
#include "../CORE/trackSelections.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/muonSelections.cc"
#include "../Tools/goodrun.cc"
#include "../CORE/triggerUtils.cc"
#include "../CORE/jetSelections.cc"
#include "../Tools/vtxreweight.cc"
#endif

bool verbose      = false;
bool doTenPercent = false;

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

//mSUGRA scan parameters-----------------------------

const bool  generalLeptonVeto = true;
const int   nm0points    = 100;
const float m0min        = 20.;
const float m0max        = 2020.;
const int   nm12points   = 38;
const float m12min       = 20.;
const float m12max       = 780.;

//---------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1.);
void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
//void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
void fillOverFlow(TH1F *h1, float value, float weight = 1.);
void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx);
void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx);
void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue,  int myType, int nJetsIdx);

//--------------------------------------------------------------------

looper::looper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
  initialized = false;
}

//--------------------------------------------------------------------

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

//--------------------------------------------------------------------


//-------------------------------------------------------
// get exact trigger name corresponding to given pattern
//-------------------------------------------------------
TString triggerName(TString triggerPattern){

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



//---------------------------------------------
// Check if trigger is unprescaled and passes
//---------------------------------------------
bool passUnprescaledHLTTriggerPattern(const char* arg){

  //cout << "Checking for pattern " << arg << endl;

  TString HLTTriggerPattern(arg);

  TString HLTTrigger = triggerName( HLTTriggerPattern );

  //cout << "Found trigger " << HLTTrigger << endl;

  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND")){
    //cout << "Didn't find trigger!" << endl;
    return false;
  }
  return passUnprescaledHLTTrigger( HLTTrigger );

}

//---------------------------------------------
// single muon triggers for lljj bump search
//---------------------------------------------

bool passMuMuJJTrigger_v1( bool isData ) {

  if( isData ){
    
    //-----------------------------------------------------------------------------
    if (evt_run() >= 160329 && evt_run() <= 163261){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu15_v5") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 163269 && evt_run() <= 164236){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v2") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 165088 && evt_run() <= 165887){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v4") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() == 166346 ){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v6") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 165922 && evt_run() <= 167043){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v5") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 167078 && evt_run() <= 170053){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v7") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 170071 && evt_run() <= 173198){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v8") )          return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 173212 && evt_run() <= 178380){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v3") )   return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 178420 && evt_run() <= 179889){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v6") )   return true;
    }
    //-----------------------------------------------------------------------------
    else if (evt_run() >= 179959 && evt_run() <= 180291){
      if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v7") )   return true;
    }
    //-----------------------------------------------------------------------------
  }

  else{
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu24_v") )  return true;
  }

  return false;
}

/*****************************************************************************************/
//passes the OS SUSY trigger selection 2011
/*****************************************************************************************/

bool passSUSYTrigger2011_v1( bool isData , int hypType , bool highpt ) {
  
  //----------------------------------------
  // no trigger requirements applied to MC
  //----------------------------------------
  
  if( !isData ) return true; 
  
  //---------------------------------
  // triggers for lepton-HT datasets
  //---------------------------------
  
  if( !highpt ) {
  
    //mm
    if( hypType == 0 ){
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT150_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu3_HT160_v") )   return true;
    }
    
    //em
    else if( hypType == 1 || hypType == 2 ){
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v") )     return true; 
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v") )     return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v") )     return true; 
      if( passUnprescaledHLTTriggerPattern("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v") )     return true;
    }
    
    //ee
    else if( hypType == 3 ){
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v") )   return true;
    }
  }
  
  //---------------------------------
  // triggers for dilepton datasets
  //---------------------------------
  
  else{
  
    //mm
    if( hypType == 0 ){
      if( passUnprescaledHLTTriggerPattern("HLT_DoubleMu7_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu13_Mu7_v" ) )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu13_Mu8_v" ) )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v" ) )   return true;
    }
    
    //em
    else if( hypType == 1 || hypType == 2 ){
      if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdL_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdL_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v") )   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v") )   return true;

    }
    
    //ee
    else if( hypType == 3 ){
      if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v") )                                   return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v") ) return true;
      if( passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) return true;
    }                                     
  }        
  
  return false;
    
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

int getIndexFromM0(float m0){
  
  float binsize = (m0max - m0min) / (float) nm0points;
  int index     = (int)(m0 - m0min) / binsize;
  return index;
    
}

//--------------------------------------------------------------------

int getIndexFromM12(float m12){
  
  float binsize = (m12max - m12min) / (float) nm12points;
  int index     = (int)(m12 - m12min) / binsize;
  return index;
    
}

//--------------------------------------------------------------------

void looper::InitBaby(){
  

  mllgen_	= -1;
  pthat_	= -1;
  qscale_	= -1;
  genmet_       = -9999;
  gensumet_     = -9999;
  genmetphi_    = -9999;   
  m0_		= -9999;
  m12_		= -9999;
  w1_		= -999;
  w2_		= -999;
  acc_2010_	= -999;
  acc_highmet_  = -999;
  acc_highht_	= -999;
  nels_		= -1;
  nmus_		= -1;
  ntaus_	= -1;
  nleps_	= -1;
  ptjetraw_	= -9999.;
  ptjet23_	= -9999.;
  ptjetF23_	= -9999.;
  ptjetO23_	= -9999.;
  cosphijz_	= -9999.;
  dilep_	= 0;

  lep1_		= 0;
  lep2_		= 0;

  jet1_		= 0;
  jet2_		= 0;
  jet3_		= 0;
  jet4_		= 0;
}

//--------------------------------------------------------------------

int getProcessType(char *prefix)
{
  int proc = -1;

  if(strcmp(prefix,"data")   == 0) proc = 0;
  if(strcmp(prefix,"Zjets")  == 0) proc = 1;
  if(strcmp(prefix,"ttdil")  == 0) proc = 2;
  if(strcmp(prefix,"ttotr")  == 0) proc = 3;
  if(strcmp(prefix,"ww")     == 0) proc = 4;
  if(strcmp(prefix,"wz")     == 0) proc = 5;
  if(strcmp(prefix,"zz")     == 0) proc = 6;
  if(strcmp(prefix,"wjets")  == 0) proc = 7;
  if(strcmp(prefix,"tW")     == 0) proc = 8;
  if(strcmp(prefix,"LM0")    == 0) proc = 10;
  if(strcmp(prefix,"LM1")    == 0) proc = 11;
  if(strcmp(prefix,"LM2")    == 0) proc = 12;
  if(strcmp(prefix,"LM3")    == 0) proc = 13;
  if(strcmp(prefix,"LM4")    == 0) proc = 14;
  if(strcmp(prefix,"LM5")    == 0) proc = 15;
  if(strcmp(prefix,"LM6")    == 0) proc = 16;
  if(strcmp(prefix,"LM7")    == 0) proc = 17;
  if(strcmp(prefix,"LM8")    == 0) proc = 18;
  if(strcmp(prefix,"LM9")    == 0) proc = 19;
  if(strcmp(prefix,"LM10")   == 0) proc = 20;
  if(strcmp(prefix,"LM11")   == 0) proc = 21;
  if(strcmp(prefix,"LM12")   == 0) proc = 22;

  return proc;
}

//--------------------------------------------------------------------

void looper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
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

bool looper::objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt) 
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

//--------------------------------------------------------------------

float getPingMass( int thisrun , int thislumi , unsigned int thisevent ){

  ifstream ifile("ping.txt");

  int run;
  int lumi;
  long event;
  float mass;
  string dummy;

  float pingmass = -1.;

  while( ifile.good() ){

    ifile >> run >> lumi >> event >> mass;

    if( run==thisrun && lumi==thislumi && event==thisevent ){
      //cout << "Found! " << run << " " << lumi << " " << event << " mass " << mass << endl;
      pingmass = mass;
    }

  }

  //if( pingmass < 0 ) cout << "ERROR! didn't find Ping mass! " << thisrun << " " << thislumi << " " << thisevent << " <<<-----------------------------------" << endl;
  ifile.close();
  return pingmass;

}

//--------------------------------------------------------------------

int looper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale, float lumi,
		      JetTypeEnum jetType, MetTypeEnum metType, ZVetoEnum zveto, FREnum frmode, bool doFakeApp, bool calculateTCMET)

{


  // Jet Corrections
  std::vector<std::string> jetcorr_pf_L2L3_filenames;
  jetcorr_pf_L2L3_filenames.clear();

  FactorizedJetCorrector *jet_pf_L2L3corrector;

  if( !initialized ){

    jetcorr_pf_L2L3_filenames.push_back("../CORE/jetcorr/data/START42_V13_AK5PF_L2L3Residual.txt");
    jet_pf_L2L3corrector = makeJetCorrector(jetcorr_pf_L2L3_filenames);

    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    //set vtx reweighting hist
    set_vtxreweight_rootfile("vtxreweight_Summer11MC_PUS4_4p7fb_Zselection.root");

    initialized = true;
  }




  bool isData = false;
  if( TString(prefix).Contains("data")  ){
    cout << "DATA!!!" << endl;
    doTenPercent = false;
    isData = true;
  }

  cout << "Doing ten percent only? " << doTenPercent << endl;

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

  int   nSkip_els_conv_dist = 0;
  float nevents             = 0.0;

  int ntot  = 0;
  int npass = 0;

  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

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

      InitBaby();

      pingmass_ = -1;
      if( TString(prefix).Contains("ping") ) pingmass_ = getPingMass(evt_run(),evt_lumiBlock(),evt_event());

      //------------------------------------------------------------------------
      // idiot check: make sure that HLT_IsoMu24 is present for all MC events
      //------------------------------------------------------------------------

      if( !isData ){
	TString trigname = triggerName("HLT_IsoMu24");

	if( trigname.Contains("TRIGGER_NOT_FOUND") ){
	  cout << "Didn't find trigger HLT_IsoMu24!!!" << endl;
	  cout << "Event " << z                                               << endl;
	  cout << "File  " << currentFile->GetTitle()                         << endl;
	  cout << evt_dataset() << " " 
	       << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	}
      }

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset() << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }
      
      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      if( !cleaning_goodDAVertexApril2011() )                        continue;
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;

      //--------------------------------
      // require trigger
      //--------------------------------

      if( !passMuMuJJTrigger_v1( isData ) ) continue;

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

      /*
	for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                                 continue;
	if( !pass_electronSelection( iel , electronSelection_el_OSV3 , false , false ) ) continue;
	goodLeptons.push_back( els_p4().at(iel) );
	lepId.push_back( els_charge().at(iel) * 11 );
	lepIndex.push_back(iel);
	ngoodel_++;
	ngoodlep_++;
	
	//cout << "Found electron " << ngoodlep_ << " pt " << els_p4().at(iel).pt() << endl;
	}
      */

      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 20 )                         continue;
	if( !muonIdNotIsolated( imu , muonSelection_mumujj ))    continue;

	goodLeptons.push_back( mus_p4().at(imu) );
	lepId.push_back( mus_charge().at(imu) * 13 );
	lepIndex.push_back(imu);
	ngoodmu_++;
	ngoodlep_++;

	//cout << "Found muon " << ngoodlep_ << " pt " << mus_p4().at(imu).pt() << endl;
      }  

      ntot++;

      // REQUIRE AT LEAST 1 GOOD LEPTON!!!
      if( goodLeptons.size() < 1 ){
	//cout << "Fail >=1 good lepton requirement" << endl;
	continue;
      }

      npass++;

      std::vector<LorentzVector> trigObjs;

      diltrig_ = passSUSYTrigger2011_v1( isData , 0 , true ) ? 1 : 0;

      if( isData ){

	//-----------------------------------------------------------------------------
	if (evt_run() >= 160329 && evt_run() <= 163261){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu15_v5")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 163269 && evt_run() <= 164236){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu24_v2")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 165088 && evt_run() <= 165887){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu24_v4")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() == 166346 ){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu24_v6")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 165922 && evt_run() <= 167043){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu24_v5")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 167078 && evt_run() <= 170053){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu24_v7")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 170071 && evt_run() <= 173198){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu24_v8")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 173212 && evt_run() <= 178380){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu30_eta2p1_v3")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 178420 && evt_run() <= 179889){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu30_eta2p1_v6")];
	}
	//-----------------------------------------------------------------------------
	else if (evt_run() >= 179959 && evt_run() <= 180291){
	  trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_IsoMu30_eta2p1_v7")];
	}
	//-----------------------------------------------------------------------------
	if( trigObjs.size() == 0 ){
	  cout << "ERROR! didn't find matching trigger objects" << endl;
	  continue;
	}
	//-----------------------------------------------------------------------------
      }

      else{

	TString trigname = triggerName("HLT_IsoMu24");
	trigObjs = cms2.hlt_trigObjs_p4()[findTriggerIndex(trigname)];

	if( trigObjs.size() == 0 ){
	  cout << "ERROR! didn't find matching trigger objects" << endl;
	  continue;
	}
      }

      //-----------------------------------------------------------------------
      // find leading lepton, require pt > 25 GeV and match to trigger object
      //-----------------------------------------------------------------------

      float maxpt   = -1;
      int   imaxpt  = -1;
      int   imaxpt2 = -1;
      
      for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){

	if( goodLeptons.at(ilep).pt() < 25 )                             continue;
	if( !objectPassTrigger( goodLeptons.at(ilep) , trigObjs ,25 ) )  continue;

	if( goodLeptons.at(ilep).pt() > maxpt ){
	  maxpt  = goodLeptons.at(ilep).pt();
	  imaxpt = ilep;
	}
      }

      //---------------------------------------------------------
      // require pt > 25 GeV lepton matched to trigger object
      //---------------------------------------------------------

      if( imaxpt < 0 ){
	hping->Fill(pingmass_);
	//cout << "didn't find muon matched to trigger object" << endl;
	continue;
      }

      id1_       = lepId.at(imaxpt);
      lep1_      = &goodLeptons.at(imaxpt);
      int index1 = lepIndex.at(imaxpt);
      iso1_      = muonIsoValue(index1,false);
      passid1_   = muonIdNotIsolated( index1 , OSGeneric_v3 ) ? 1 : 0;

      lep1trk_   = &(mus_trk_p4().at(imaxpt));
      lep1glb_   = &(mus_gfit_p4().at(imaxpt));
      lep1sta_   = &(mus_sta_p4().at(imaxpt));

      //cout << "Leading lepton: pt " << lep1_->pt() << " id " << id1_ << endl;

      //---------------------------------------------
      // find 2nd leading lepton (if >=2 leptons)
      //---------------------------------------------

      id2_       = -999;
      lep2_      = 0;
      int index2 = -1;
      dilmass_   = -999;
      iso2_      = -1.;
      passid2_   = -1;
      
      lep2trk_   = 0;
      lep2glb_   = 0;
      lep2sta_   = 0;

      if( ngoodlep_ > 1 ){

	maxpt = -1;
	
	for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  
	  if( ilep == imaxpt ) continue;
	  
	  if( goodLeptons.at(ilep).pt() > maxpt ){
	    maxpt   = goodLeptons.at(ilep).pt();
	    imaxpt2 = ilep;
	  }
	}
	
	if( imaxpt2 < 0 ){
	  cout << "ERROR! QUITTING imaxpt2 " << imaxpt2 << endl;
	  exit(3);
	}
	
	id2_       = lepId.at(imaxpt2);
	lep2_      = &goodLeptons.at(imaxpt2);
	index2     = lepIndex.at(imaxpt2);
	dilmass_   = (goodLeptons.at(imaxpt) + goodLeptons.at(imaxpt2)).mass();
	iso2_      = muonIsoValue(index2,false);
	passid2_   = muonIdNotIsolated( index2, OSGeneric_v3 ) ? 1 : 0;

	lep2trk_   = &(mus_trk_p4().at(imaxpt2));
	lep2glb_   = &(mus_gfit_p4().at(imaxpt2));
	lep2sta_   = &(mus_sta_p4().at(imaxpt2));

      }

      //cout << "2nd deading lepton: pt " << lep2_->pt() << " id " << id2_ << " mass " << dilmass_ << endl;
	  
      //--------------------------------
      // get MC quantities
      //--------------------------------
      
      int nels       =  0;
      int nmus       =  0;
      int ntaus      =  0;
      int nleps      =  0;
      float dilptgen = -1;
      /*
      if( !isData ){

	w1_ = -999;
	//w1_     = leptonOrTauIsFromW( index1 , id1_ , true );
	pthat_  = genps_pthat();
	qscale_ = genps_qScale();
	
	//splitting ttbar into ttdil/ttotr
	//nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
	nleps = leptonGenpCount(nels, nmus, ntaus);
	
	nels_  = nels;
	nmus_  = nmus;
	ntaus_ = ntaus;
	nleps_ = nleps;

	if( strcmp(prefix,"ttem")  == 0 && ( nels + nmus ) != 2 ) continue;
	if( strcmp(prefix,"ttdil") == 0 && nleps != 2           ) continue;
	if( strcmp(prefix,"ttotr") == 0 && nleps == 2           ) continue;
	
	LorentzVector vdilepton(0,0,0,0);
        
	for ( int igen = 0 ; igen < genps_id().size() ; igen++ ) { 
	  if ( abs( cms2.genps_id().at(igen) ) == 11) vdilepton += genps_p4().at(igen); 
	  if ( abs( cms2.genps_id().at(igen) ) == 13) vdilepton += genps_p4().at(igen); 
	}
        
	if( nels + nmus == 2) dilptgen = vdilepton.pt();
        
	if ( strcmp(prefix , "DYee"     ) == 0 &&  nels  != 2  ) continue;
	if ( strcmp(prefix , "DYmm"     ) == 0 &&  nmus  != 2  ) continue;
	if ( strcmp(prefix , "DYtautau" ) == 0 &&  ntaus != 2  ) continue;
	
	//splice together the DY samples - if its madgraph, then we do nothing
	if(TString(prefix).Contains("DY") && TString(evt_dataset()).Contains("madgraph") == false) {	
	  bool doNotContinue = false;
	  for(unsigned int i = 0; i < genps_p4().size(); i++){
	    if(abs(genps_id()[i]) == 23 && genps_p4()[i].M() > 50.)
	      doNotContinue = true;
	  }
	  if(doNotContinue)
	    continue;	
	}
	
	//extract pthat
	if(TString(prefix).Contains("DY")){
	  int nz = 0;
	  for(unsigned int i = 0; i < genps_p4().size(); i++){
	    if(abs(genps_id()[i]) == 23){
	      mllgen_ = genps_p4()[i].M();
	      nz++;
	    }
	  }
	  if(nz != 1 ) cout << "ERROR NZ " << nz << endl;
	}
      }
      */

      //pfjets
      VofP4 vpfjets_p4;

      njets_     = 0;
      ht_        = 0.;
      nbtags17_  = 0;
      nbtags20_  = 0;
      nbtags33_  = 0;
      nbtags20_24_  = 0;

      //---------------------
      // apply residual JEC
      //---------------------

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
	float         jet_cor = jetCorrection(cms2.pfjets_p4().at(ijet), jet_pf_L2L3corrector);
	LorentzVector vjet    = pfjets_corL1FastL2L3().at(ijet) * jet_cor * pfjets_p4().at(ijet);

	if( generalLeptonVeto ){
	  bool rejectJet = false;
	  for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	    if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.5 ) rejectJet = true;  
	  }
	  if( rejectJet ) continue;
	}

	if( vjet.pt() < 25         ) continue;
	if( fabs(vjet.eta()) > 3.0 ) continue;

	njets_++;
	ht_ += vjet.pt();

	vpfjets_p4.push_back( vjet );

	if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 1.7 ){
	  nbtags17_++;
	}

	if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 2.0 ){
	  nbtags20_++;
	}

	if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 2.0 && abs( vjet.eta() ) < 2.4 ){
	  nbtags20_24_++;
	}

	if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 3.3 ){
	  nbtags33_++;
	}
      }
	  
      //-----------------------------
      // DON'T apply residual JEC
      //-----------------------------

      njetsuncor_ = 0;
      VofP4 vpfjets_uncor_p4;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
	LorentzVector vjet    = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet);

	if( generalLeptonVeto ){
	  bool rejectJet = false;
	  for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	    if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.5 ) rejectJet = true;  
	  }
	  if( rejectJet ) continue;
	}

	if( vjet.pt() < 25         ) continue;
	if( fabs(vjet.eta()) > 3.0 ) continue;

	njetsuncor_++;
	vpfjets_uncor_p4.push_back( vjet );

      }
	  
      //--------------------------------------
      // DON'T apply any jet corrections
      //--------------------------------------

      njetsplain_ = 0;
      VofP4 vpfjets_plain_p4;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
	LorentzVector vjet    = pfjets_p4().at(ijet);

	if( generalLeptonVeto ){
	  bool rejectJet = false;
	  for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	    if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.5 ) rejectJet = true;  
	  }
	  if( rejectJet ) continue;
	}

	if( vjet.pt() < 25         ) continue;
	if( fabs(vjet.eta()) > 3.0 ) continue;

	njetsplain_++;
	vpfjets_plain_p4.push_back( vjet );

      }
	  
      //------------------------------------
      // require >=2 jets!!!!!!
      //------------------------------------

      // if( njets_ < 2 ){
      // 	cout << "less than 2 jets " << njets_ << endl;
      // 	hping->Fill(pingmass_);
      // 	continue;
      // }

      jet1_  = 0;
      jet2_  = 0;
      jet3_  = 0;
      jet4_  = 0;

      mmjj_      = -1;
      mmjjtrk_   = -1;
      mmjjglb_   = -1;
      mmjjsta_   = -1;
      mjjj_      = -1;
      mmjjuncor_ = -1;

      //--------------------------------------
      // corrected jets
      //--------------------------------------

      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vpfjets_p4_sorted(vpfjets_p4);
      sort(vpfjets_p4_sorted.begin(), vpfjets_p4_sorted.end(), sortByPt);

      if( njets_ > 0 ) jet1_  = &vpfjets_p4_sorted.at(0);
      if( njets_ > 1 ) jet2_  = &vpfjets_p4_sorted.at(1);
      if( njets_ > 2 ) jet3_  = &vpfjets_p4_sorted.at(2);
      if( njets_ > 3 ) jet4_  = &vpfjets_p4_sorted.at(3);

      if( njets_ >= 2 ){

	jetcor1_ = jetCorrection( vpfjets_p4_sorted.at(0), jet_pf_L2L3corrector);
	jetcor2_ = jetCorrection( vpfjets_p4_sorted.at(1), jet_pf_L2L3corrector);
	
	if( ngoodlep_ > 1 ){
	  mmjj_    = (*lep1_    + *lep2_    + *jet1_ + *jet2_).mass();
	  mmjjtrk_ = (*lep1trk_ + *lep2trk_ + *jet1_ + *jet2_).mass();
	  mmjjglb_ = (*lep1glb_ + *lep2glb_ + *jet1_ + *jet2_).mass();
	  mmjjsta_ = (*lep1sta_ + *lep2sta_ + *jet1_ + *jet2_).mass();
	}
	
	if( njets_ > 2    ) mjjj_ = (*lep1_+*jet1_+*jet2_+*jet3_).mass();
      }

      //--------------------------------------
      // UN-corrected jets
      //--------------------------------------

      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vpfjets_uncor_p4_sorted(vpfjets_uncor_p4);
      sort(vpfjets_uncor_p4_sorted.begin(), vpfjets_uncor_p4_sorted.end(), sortByPt);

      if( njetsuncor_ >= 2 && ngoodlep_ > 1 ){
	mmjjuncor_ = (*lep1_    + *lep2_    + vpfjets_uncor_p4_sorted.at(0) + vpfjets_uncor_p4_sorted.at(1)).mass();
      }

      if( isData ) mmjjdef_  = mmjj_;
      else         mmjjdef_  = mmjjuncor_;
     
      
      //--------------------------------------
      // plain jets
      //--------------------------------------

      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vpfjets_plain_p4_sorted(vpfjets_plain_p4);
      sort(vpfjets_plain_p4_sorted.begin(), vpfjets_plain_p4_sorted.end(), sortByPt);

      mmjjc_  = -1;
      cor1_   = -1;
      cor2_   = -1;
      metnew_ = -1;

      if( njetsplain_ >= 2 && ngoodlep_ > 1 ){

	float a = vpfjets_plain_p4_sorted.at(0).x();
	float b = vpfjets_plain_p4_sorted.at(1).x();
	float c = vpfjets_plain_p4_sorted.at(0).y();
	float d = vpfjets_plain_p4_sorted.at(1).y();
	
	float det = a * d - b * c;

	if( fabs(det) > 1e-10 ){
	  
	  float ap = d / det;
	  float bp = -1 * b / det;
	  float cp = -1 * c / det;
	  float dp = a / det;

	  float metx = evt_pfmet() * cos( evt_pfmetPhi() );
	  float mety = evt_pfmet() * sin( evt_pfmetPhi() );

	  cor1_ = 1 + ap * metx + bp * mety;
	  cor2_ = 1 + cp * metx + dp * mety;

	  if( cor1_ > 0 && cor2_ > 0 ){
	   
	    mmjjc_ = (*lep1_    + *lep2_    + cor1_ * vpfjets_plain_p4_sorted.at(0) + cor2_ * vpfjets_plain_p4_sorted.at(1)).mass();

	    float metxnew = metx - (cor1_-1) * vpfjets_plain_p4_sorted.at(0).x() - (cor2_-1) * vpfjets_plain_p4_sorted.at(1).x();
	    float metynew = mety - (cor1_-1) * vpfjets_plain_p4_sorted.at(0).y() - (cor2_-1) * vpfjets_plain_p4_sorted.at(1).y();
	    metnew_  = sqrt( metxnew * metxnew + metynew * metynew );
	    //cout << "new met " << metnew << endl;

	  }
	  else{
	    mmjjc_ = -3;
	  }
	}

	else{
	  mmjjc_ = -2;
	}


      }




      pfmet_    = evt_pfmet();
      pfmetphi_ = evt_pfmetPhi();
      pfsumet_  = evt_pfsumet();

      if( isData ){
	weight_ = 1;
      }

      else{

	weight_ = kFactor * evt_scale1fb() * lumi;

	if( doTenPercent )	  weight_ *= 10;
      }

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
      ndavtxweight_ = vtxweight(isData,true);
      hbhe_         = evt_hbheFilter();
      
      nevents += weight_;

      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << "nevents " << nevents << endl;
  cout << endl;
  cout << ntot << " total events, " << npass << " pass events" << endl;

  if(g_createTree) closeTree();
  
  //already_seen.clear();
  
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  //delete d_llsol; //REPLACETOPMASS

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
  hping = new TH1F("hping","hping",2000,0,2000);

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()

//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight)
{
  float maxx = h2->GetXaxis()->GetXmax();
  float minx = h2->GetXaxis()->GetXmin();
  float maxy = h2->GetYaxis()->GetXmax();
  float miny = h2->GetYaxis()->GetXmin();

  if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (xvalue < minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
  if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
  if (yvalue < miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

  h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

// void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue)
// {
//   float maxx = h2->GetXaxis()->GetXmax();
//   float minx = h2->GetXaxis()->GetXmin();
//   float maxy = h2->GetYaxis()->GetXmax();
//   float miny = h2->GetYaxis()->GetXmin();

//   if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
//   if (xvalue < minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
//   if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
//   if (yvalue < miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

//   h2->Fill(xvalue, yvalue);
// }

//--------------------------------------------------------------------

void fillOverFlow(TH1F *h1, float value, float weight)
{
  float max = h1->GetXaxis()->GetXmax();
  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight)
{
  float maxx = h2->GetXaxis()->GetXmax();
  float maxy = h2->GetYaxis()->GetXmax();

  if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());

  h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{
  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      
}

//--------------------------------------------------------------------

void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx)
{
  fillUnderOverFlow(h2[myType][nJetsIdx], xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[myType][3],        xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[3][nJetsIdx],      xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[3][3],             xvalue, yvalue, weight);      
}

//--------------------------------------------------------------------

void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue, int myType, int nJetsIdx)
{
  h2[myType][nJetsIdx] -> Fill(xvalue, yvalue);      
  h2[myType][3]        -> Fill(xvalue, yvalue);      
  h2[3][nJetsIdx]      -> Fill(xvalue, yvalue);      
  h2[3][3]             -> Fill(xvalue, yvalue);      
}

//--------------------------------------------------------------------

//--------------------------------------------------------------------

void looper::makeTree(char *prefix, bool doFakeApp, FREnum frmode ){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();


  //char* dir = "";
  //if     ( g_trig == e_lowpt  ) dir = "lowpt";
  //else if( g_trig == e_highpt ) dir = "highpt";

  //Super compressed ntuple here
  char* frsuffix = "";
  if( doFakeApp ){
    if ( frmode == e_qcd   ) frsuffix = "_doubleFake";
    if ( frmode == e_wjets ) frsuffix = "_singleFake";
  }

  char* tpsuffix = "";
  if( doTenPercent ) tpsuffix = "_tenPercent";

  outFile   = new TFile(Form("../output/%s/%s_smallTree%s%s.root",g_version,prefix,frsuffix,tpsuffix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in looper.h
  outTree->Branch("diltrig",         &diltrig_,          "diltrig/I");
  outTree->Branch("pingmass",        &pingmass_,         "pingmass/F");
  outTree->Branch("metnew",          &metnew_,           "metnew/F");
  outTree->Branch("mmjjdef",         &mmjjdef_,          "mmjjdef/F");
  outTree->Branch("mmjjc",           &mmjjc_,            "mmjjc/F");
  outTree->Branch("cor1",            &cor1_,             "cor1/F");
  outTree->Branch("cor2",            &cor2_,             "cor2/F");
  outTree->Branch("njetsuncor",      &njetsuncor_,       "njetsuncor/I");
  outTree->Branch("njetsplain",      &njetsplain_,       "njetsplain/I");
  outTree->Branch("jetcor1",         &jetcor1_,          "jetcor1/F");
  outTree->Branch("jetcor2",         &jetcor2_,          "jetcor2/F");
  outTree->Branch("hbhe",            &hbhe_,             "hbhe/I");
  outTree->Branch("json",            &json_,             "json/I");
  outTree->Branch("weight",          &weight_,           "weight/F");
  outTree->Branch("pthat",           &pthat_,            "pthat/F");
  outTree->Branch("qscale",          &qscale_,           "qscale/F");
  outTree->Branch("mllgen",          &mllgen_,           "mllgen/F");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("ngoodmu",         &ngoodmu_,          "ngoodmu/I");
  outTree->Branch("dilmass",         &dilmass_,          "dilmass/F");
  outTree->Branch("genmet",          &genmet_,           "genmet/F");
  outTree->Branch("genmet",          &genmet_,           "genmet/F");
  outTree->Branch("genmet",          &genmet_,           "genmet/F");
  outTree->Branch("pfmet",           &pfmet_,            "pfmet/F");
  outTree->Branch("pfmetveto",       &pfmetveto_,        "pfmetveto/F");
  outTree->Branch("pfmetsig",        &pfmetsig_,         "pfmetsig/F");
  outTree->Branch("pfmetphi",        &pfmetphi_,         "pfmetphi/F");
  outTree->Branch("pfsumet",         &pfsumet_,          "pfsumet/F");
  outTree->Branch("mucormet",        &mucormet_,         "mucormet/F");
  outTree->Branch("mucorjesmet",     &mucorjesmet_,      "mucorjesmet/F");
  outTree->Branch("tcmet35X",        &tcmet_35X_,        "tcmet35X/F");
  outTree->Branch("tcmetevent",      &tcmet_event_,      "tcmetevent/F");
  outTree->Branch("tcmetlooper",     &tcmet_looper_,     "tcmetlooper/F");
  outTree->Branch("tcmetphi",        &tcmetphi_,         "tcmetphi/F");
  outTree->Branch("tcsumet",         &tcsumet_,          "tcsumet/F");
  outTree->Branch("tcmetUp",         &tcmetUp_,          "tcmetUp/F");
  outTree->Branch("tcmetDown",       &tcmetDown_,        "tcmetDown/F");
  outTree->Branch("tcmetTest",       &tcmetTest_,        "tcmetTest/F");
  outTree->Branch("pfmetUp",         &pfmetUp_,          "pfmetUp/F");
  outTree->Branch("pfmetDown",       &pfmetDown_,        "pfmetDown/F");
  outTree->Branch("pfmetTest",       &pfmetTest_,        "pfmetTest/F");
  outTree->Branch("sumjetpt",        &sumjetpt_,         "sumjetpt/F");
  outTree->Branch("dileta",          &dileta_,           "dileta/F");
  outTree->Branch("dilpt",           &dilpt_,            "dilpt/F");
  outTree->Branch("dildphi",         &dildphi_,          "dildphi/F");
  outTree->Branch("njets",           &njets_,            "njets/I");
  outTree->Branch("ngenjets",        &ngenjets_,         "ngenjets/I");
  outTree->Branch("njpt",            &njpt_,             "njpt/I");
  outTree->Branch("npfjets25",       &npfjets25_,        "npfjets25/I");
  outTree->Branch("npfjets40",       &npfjets40_,        "npfjets40/I");
  outTree->Branch("npfjetspv",       &npfjetspv_,        "npfjetspv/I");
  outTree->Branch("njetsUp",         &njetsUp_,          "njetsUp/I");
  outTree->Branch("njetsDown",       &njetsDown_,        "njetsDown/I");
  outTree->Branch("nvtx",            &nvtx_,             "nvtx/I");
  outTree->Branch("ndavtx",          &ndavtx_,           "ndavtx/I");
  outTree->Branch("ndavtxweight",    &ndavtxweight_,     "ndavtxweight/F");
  outTree->Branch("nbtags17",        &nbtags17_,         "nbtags17/I");
  outTree->Branch("nbtags20",        &nbtags20_,         "nbtags20/I");
  outTree->Branch("nbtags2024",      &nbtags20_24_,      "nbtags2024/I");
  outTree->Branch("nbtags33",        &nbtags33_,         "nbtags33/I");
  outTree->Branch("m0",              &m0_,               "m0/F");
  outTree->Branch("m12",             &m12_,              "m12/F");
  outTree->Branch("id1",             &id1_,              "id1/I");
  outTree->Branch("id2",             &id2_,              "id2/I");
  outTree->Branch("w1",              &w1_,               "w1/I");
  outTree->Branch("w2",              &w2_,               "w2/I");
  outTree->Branch("iso1",            &iso1_,             "iso1/F");
  outTree->Branch("iso2",            &iso2_,             "iso2/F");
  outTree->Branch("passid1",         &passid1_,          "passid1/I");
  outTree->Branch("passid2",         &passid2_,          "passid2/I");
  outTree->Branch("mt",              &mt_,               "mt/F");
  outTree->Branch("dataset",         &dataset_,          "dataset[200]/C");
  outTree->Branch("run",             &run_,              "run/I");
  outTree->Branch("lumi",            &lumi_,             "lumi/I");
  outTree->Branch("event",           &event_,            "event/I");
  outTree->Branch("ht",              &ht_,               "ht/F");  
  outTree->Branch("htjpt",           &htjpt_,            "htjpt/F");  
  outTree->Branch("htpf25",          &htpf25_,           "htpf25/F");  
  outTree->Branch("htpf40",          &htpf40_,           "htpf40/F");  
  outTree->Branch("htpfpv",          &htpfpv_,           "htpfpv/F");  
  outTree->Branch("nels",            &nels_,             "nels/I");  
  outTree->Branch("nmus",            &nmus_,             "nmus/I");  
  outTree->Branch("ntaus",           &ntaus_,            "ntaus/I");  
  outTree->Branch("nleps",           &nleps_,            "nleps/I");  
  outTree->Branch("dphijm",          &dphijm_,           "dphijm/F");  
  outTree->Branch("ptjetraw",        &ptjetraw_,         "ptjetraw/F");  
  outTree->Branch("mmjj",            &mmjj_,             "mmjj/F");  
  outTree->Branch("mmjjuncor",       &mmjjuncor_,        "mmjjuncor/F");  
  outTree->Branch("mmjjtrk",         &mmjjtrk_,          "mmjjtrk/F");  
  outTree->Branch("mmjjglb",         &mmjjglb_,          "mmjjglb/F");  
  outTree->Branch("mmjjsta",         &mmjjsta_,          "mmjjsta/F");  
  outTree->Branch("mjjj",            &mjjj_,             "mjjj/F");  

  outTree->Branch("lep1trk" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1trk_	);
  outTree->Branch("lep1glb" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1glb_	);
  outTree->Branch("lep1sta" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1sta_	);
  outTree->Branch("lep2trk" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2trk_	);
  outTree->Branch("lep2glb" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2glb_	);
  outTree->Branch("lep2sta" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2sta_	);

  outTree->Branch("lep1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_	);
  outTree->Branch("lep2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_	);
  outTree->Branch("jet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet1_	);
  outTree->Branch("jet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet2_	);
  outTree->Branch("jet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet3_	);
  outTree->Branch("jet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet4_	);


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
