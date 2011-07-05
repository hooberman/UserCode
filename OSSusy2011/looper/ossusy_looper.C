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
#include "ossusy_looper.h"
#include "getMt2.C"
#include "TTreeCache.h"
#include "../CORE/CMS2.h"
#include "../CORE/metSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/mcSelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/MT2/MT2.h"
#include "../Tools/goodrun.cc"
#include "../CORE/utilities.cc"
#include "../CORE/ttbarSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/mcSUSYkfactor.h"
#include "../CORE/triggerSuperModel.h"
//#include "../CORE/jetSelections.h"
#include "../Tools/vtxreweight.cc"


bool verbose = false;

//#include "../CORE/topmass/getTopMassEstimate.icc" // REPLACETOPMASS
//#include "../CORE/triggerUtils.cc"

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

//mSUGRA scan parameters-----------------------------

const bool  generalLeptonVeto = true;
const int   nm0points    = 81;
const float m0min        = 0.;
const float m0max        = 4050.;
const int   nm12points   = 26;
const float m12min       = 100.;
const float m12max       = 620.;

TH1F* hsusydilPt[nm0points][nm12points]; 
TH1F* hsusytcmet[nm0points][nm12points]; 
TH2F* hsusy_met_sumjetpt[nm0points][nm12points]; 

//TH1F* msugra_

//---------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1.);
void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
//void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
void fillOverFlow(TH1F *h1, float value, float weight = 1.);
void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx);
void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx);
void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue,  int myType, int nJetsIdx);
float returnSigma(float sumJetPt, ossusy_looper::MetTypeEnum metType);
float returnBias(float sumJetPt, ossusy_looper::MetTypeEnum metType);

//--------------------------------------------------------------------

void checkElectron( int elidx ){

  cout << "Check electron" << endl;
  cout << "Pass all    " << pass_electronSelection( elidx , electronSelection_el_OSV3			) << endl;
  cout << "Pass ID     " << pass_electronSelection( elidx , electronSelection_el_OSV3_noiso		) << endl;
  cout << "Pass iso    " << pass_electronSelection( elidx , electronSelection_el_OSV3_iso		) << endl;
  cout << "VBTF90      " << pass_electronSelection( elidx , 1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL	) << endl;
  cout << "PV          " << pass_electronSelection( elidx , 1ll<<ELEIP_PV_OSV2				) << endl;
  cout << "nomuon      " << pass_electronSelection( elidx , 1ll<<ELENOMUON_010				) << endl;
  cout << "hitpattern  " << pass_electronSelection( elidx , 1ll<<ELENOTCONV_HITPATTERN			) << endl;
  cout << "convrej     " << pass_electronSelection( elidx , 1ll<<ELENOTCONV_DISTDCOT002			) << endl;
  cout << "pt10        " << pass_electronSelection( elidx , 1ll<<ELEPT_010				) << endl;
  cout << "eta25       " << pass_electronSelection( elidx , 1ll<<ELEETA_250				) << endl;
  cout << "transition  " << pass_electronSelection( elidx , 1ll<<ELE_NOT_TRANSITION			) << endl;
  cout << "HLT iso     " << pass_electronSelection( elidx , 1ll<<ELEISO_ECAL_RELNT020_NPS		) << endl;
  cout << "offline iso " << pass_electronSelection( elidx , 1ll<<ELEISO_RELNT015			) << endl;

}

void checkMuon( int muidx ){

  cout << "Check muon" << endl;
  cout << "Pass all  " <<  muonId(muidx , OSGeneric_v3)                                            << endl;
  cout << "Pass ID   " <<  muonIdNotIsolated(muidx , OSGeneric_v3 )                                << endl;
  cout << "Pass iso  " <<  ( muonIsoValue(muidx,false) < 0.15 )                                    << endl;
  cout << "eta24     " <<  ( TMath::Abs(cms2.mus_p4()[muidx].eta()) < 2.4)                         << endl;
  cout << "chi2/ndf  " <<  ( cms2.mus_gfit_chi2().at(muidx)/cms2.mus_gfit_ndof().at(muidx) < 10)   << endl;
  cout << "global    " <<  ( ((cms2.mus_type().at(muidx)) & (1<<1)) != 0)                          << endl;
  cout << "tracker   " <<  ( ((cms2.mus_type().at(muidx)) & (1<<2)) != 0)                          << endl;
  cout << "nhits     " <<  ( cms2.mus_validHits().at(muidx) > 10)                                  << endl;
  cout << "stahits   " <<  ( cms2.mus_gfit_validSTAHits().at(muidx) != 0)                          << endl;
  cout << "d0PV      " <<  ( TMath::Abs(mud0PV_smurfV3(muidx)) < 0.02)                             << endl;
  cout << "dzPV      " <<  ( TMath::Abs(mudzPV_smurfV3(muidx)) < 1  )                              << endl;
  cout << "dpt/pt    " <<  ( cms2.mus_ptErr().at(muidx)/cms2.mus_p4().at(muidx).pt()<0.1)          << endl;

}

//--------------------------------------------------------------------

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
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

void ossusy_looper::InitBaby(){

  w1_		= -999;
  w2_		= -999;

  acc_2010_	= -999;
  acc_highmet_  = -999;
  acc_highht_	= -999;

  dilep_	= 0;
  jet_		= 0;
  lep1_		= 0;
  lep2_		= 0;

  nels_		= -1;
  nmus_		= -1;
  ntaus_	= -1;

  ptjetraw_	= -9999.;
  ptjet23_	= -9999.;
  ptjetF23_	= -9999.;
  ptjetO23_	= -9999.;
  cosphijz_	= -9999.;

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

pair<float, float> ScaleMET( pair<float, float> p_met, LorentzVector p4_dilep, double rescale = 1.0 ){
  float met = p_met.first;
  float metPhi = p_met.second;
  float metx = met*cos(metPhi);
  float mety = met*sin(metPhi);

  float lepx = p4_dilep.Px();
  float lepy = p4_dilep.Py();
      
  //hadronic component of MET (well, mostly), scaled
  float metHx = (metx + lepx)*rescale;
  float metHy = (mety + lepy)*rescale;
  float metNewx = metHx - lepx;
  float metNewy = metHy - lepy;
  float metNewPhi = atan2(metNewy, metNewx);
      
  pair<float, float> p_met2 = make_pair(sqrt(metNewx*metNewx + metNewy*metNewy), metNewPhi);
  return p_met2;
}

//--------------------------------------------------------------------

void ossusy_looper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

//--------------------------------------------------------------------

bool nkcut(const unsigned value,const int n,
	   const int x0 =-1, const int x1 =-1, const int x2 =-1, const int x3 =-1,
	   const int x4 =-1, const int x5 =-1, const int x6 =-1, const int x7 =-1,
	   const int x8 =-1, const int x9 =-1, const int x10=-1, const int x11=-1,
	   const int x12=-1, const int x13=-1, const int x14=-1, const int x15=-1,
	   const int x16=-1, const int x17=-1, const int x18=-1, const int x19=-1,
	   const int x20=-1, const int x21=-1, const int x22=-1, const int x23=-1,
	   const int x24=-1, const int x25=-1, const int x26=-1, const int x27=-1,
	   const int x28=-1, const int x29=-1, const int x30=-1)
{
  //if (value<0) return false;
  for (int i=0;i<n;++i) {
    if (i==x0 ||i==x1 ||i==x2 ||i==x3 ||i==x4 ||
        i==x5 ||i==x6 ||i==x7 ||i==x8 ||i==x9 ||
        i==x10||i==x11||i==x12||i==x13||i==x14||
        i==x15||i==x16||i==x17||i==x18||i==x19||
        i==x20||i==x21||i==x22||i==x23||i==x24||
        i==x25||i==x26||i==x27||i==x28||i==x29||i==x30) continue;
    if (value&(1<<i)) {} else { return false;}
  }
  return true;
}

//--------------------------------------------------------------------

ossusy_looper::ossusy_looper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
  initialized = false;

}

//--------------------------------------------------------------------

int ossusy_looper::PassGenSelectionOS( bool isData , float metcut , float htcut , float ycut ){
  
  if( isData ) return -999.;
  
  //---------------------------------------------
  // get gen leptons
  //---------------------------------------------

  VofP4 mcLeptons;
  vector<int> mcId;

  //first, find electrons and muons
  for (size_t i = 0; i < cms2.genps_id().size(); ++i){
    
    //electron or muon
    if (!( abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13))       continue;

    //pt > 10 GeV, |eta| < 2.5
    if ( cms2.genps_p4()[i].Pt() < 10.0 || fabs(cms2.genps_p4()[i].Eta()) > 2.5)  continue;

    mcId.push_back(      cms2.genps_id()[i] );
    mcLeptons.push_back( cms2.genps_p4()[i] );
  }

  //second, look for tau->electron/muon
  for (size_t i = 0; i < cms2.genps_id().size(); ++i){
   
    //tau
    if (!( abs(cms2.genps_id()[i]) == 15 ) ) continue;

    //did this tau decay leptonically?
    bool lepTauDecay = false;

    for(unsigned int k = 0; k < cms2.genps_lepdaughter_id()[i].size(); k++) {
      int daughter = abs(cms2.genps_lepdaughter_id()[i][k]);

      if( daughter == 12 || daughter == 16 ) lepTauDecay = true;
    }

    //if tau decayed leptonically, find daughter electron/muon
    if( !lepTauDecay ) continue;

    for(unsigned int k = 0; k < cms2.genps_lepdaughter_id()[i].size(); k++) {
      int daughter = abs(cms2.genps_lepdaughter_id()[i][k]);

      if( ! ( daughter == 11 || daughter == 13) ) continue;

      //pt > 10 GeV, |eta| < 2.5
      if ( cms2.genps_lepdaughter_p4()[i][k].Pt() < 10.0 || fabs(cms2.genps_lepdaughter_p4()[i][k].Eta()) > 2.5)  continue;

      mcId.push_back( genps_lepdaughter_id()[i][k] );
      mcLeptons.push_back( genps_lepdaughter_p4()[i][k] );
    }

  }
  
  if( mcLeptons.size() < 2 ) return -1;

  //---------------------------------------------
  // look for OS pt > (20,10) GeV pair, Z-veto
  //---------------------------------------------

  bool foundPair = false;

  for( unsigned int i = 0 ; i < mcLeptons.size() ; ++i ){

    for( unsigned int j = i + 1 ; j < mcLeptons.size() ; ++j ){

      //20,10
      if( max( mcLeptons[i].pt() , mcLeptons[j].pt() ) < 20 ) continue;
      if( min( mcLeptons[i].pt() , mcLeptons[j].pt() ) < 10 ) continue;

      //OS
      if ( mcId[i] * mcId[j] > 0 )                            continue;

      //SF?
      bool SF = ( abs( mcId[i] ) == abs( mcId[j] ) );

      //Z mass veto SF pairs
      float dilmass = ( mcLeptons[i] + mcLeptons[j] ).mass();
      if( SF && dilmass > 76.0 && dilmass < 106. ) continue;
	
      //found OS pair!
      foundPair = true;
	     
    }
  }

  if( !foundPair ) return -2;

  //---------------------------------------------
  // get MC HT, njets
  //---------------------------------------------
    
  int   njets   = 0;
  float ht      = 0.;

  for (size_t j = 0; j < cms2.evt_ngenjets(); ++j) {

    if (cms2.genjets_p4()[j].Pt() < 30.0)       continue;
    if (fabs(cms2.genjets_p4()[j].Eta()) > 3.0) continue;
    bool clean = true;
    for ( size_t i = 0; i < mcLeptons.size(); ++i) 
      {
	if (ROOT::Math::VectorUtil::DeltaR(cms2.genjets_p4()[j], mcLeptons[i]) < 0.4) {
	  clean = false;
	  break;
	}
      }
    if (clean){
      njets ++;
      ht += genjets_p4()[j].Pt();
    }
  }
  
  //-----------------------------------
  // calculate gen-level MET, y
  //-----------------------------------
    
  float met = gen_met();
  float y   = ht > 0 ? met/sqrt(ht) : 0;
      
  if( njets < 2      ) return -3;
  if(    ht < htcut  ) return -4;
  if(   met < metcut ) return -5;
  if(     y < ycut   ) return -6;

  //pass!!!
  return 1;
    
}


int ossusy_looper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale, float lumi,
                             JetTypeEnum jetType, MetTypeEnum metType, ZVetoEnum zveto, FREnum frmode, bool doFakeApp, bool calculateTCMET)
{


  if( !initialized ){
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    set_vtxreweight_rootfile("vtxreweight_Spring11MC_336pb_Zselection.root",true);

    initialized = true;
  }

  float minpt  = -1;
  float maxpt  = -1;
  char* dir    = "";
  float htcut  = -1;
  bool  highpt = false;

  if( g_trig == e_lowpt ){
    cout << "Doing 10,5 selection" << endl;
    minpt = 5.;
    maxpt = 10.;
    htcut = 200.;
    dir   = "lowpt";
  }

  else if( g_trig == e_highpt ){
    cout << "Doing 20,10 selection" << endl;
    minpt = 10.;
    maxpt = 20.;
    htcut = 100.;
    dir   = "highpt";
    highpt = true;
  }
  
  bool isLM = TString(prefix).Contains("LM");
  int nSS = 0;
  int nOS = 0;

  bool isData = false;
  if( TString(prefix).Contains("data")  ){
    cout << "DATA!!!" << endl;
    isData = true;
  }
  // instanciate topmass solver REPLACETOPMASS
  //ttdilepsolve * d_llsol = new ttdilepsolve;


  //instantiate SimpleFakeRate class for electrons and muons
  //this is the default, can change it below if needed
  SimpleFakeRate* mufr = 0;
  SimpleFakeRate* elfr = 0;

  if(doFakeApp) {

    std::cout<<"**************************"<<std::endl;
    std::cout<<"Running FR application job"<<std::endl;
    std::cout<<"**************************"<<std::endl;

    if(isData) {
      std::cout<<"Using data derived FR files"<<std::endl;
      mufr = new SimpleFakeRate("fr_os7June2011.root", "fr_mu_OSGV3" );
      elfr = new SimpleFakeRate("fr_os7June2011.root", "fr_el_OSGV3" );
    }
    else {
      std::cout<<"Using data derived FR files"<<std::endl;
      std::cout<<"CURRENTLY USING DATA FR FOR MC FIXME!!!!!" <<std::endl;
      mufr = new SimpleFakeRate("fr_os7June2011.root", "fr_mu_OSGV3" );
      elfr = new SimpleFakeRate("fr_os7June2011.root", "fr_el_OSGV3" );
    }
  }

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

  int nGoodEl = 0;
  int nGoodMu = 0;
  int nGoodEM = 0;
  int nSkip_els_conv_dist = 0;

  int nevent = 0;
  float nee = 0.;
  float nmm = 0.;
  float nem = 0.;

  float neetot = 0.;
  float nmmtot = 0.;
  float nemtot = 0.;

  int ngen          = 1;

  if     ( TString(prefix).Contains("LM1") ) ngen = 219190;
  else if( TString(prefix).Contains("LM3") ) ngen = 220000;
  else if( TString(prefix).Contains("LM6") ) ngen = 219190;
  else if( TString(prefix).Contains("LM")  ){
    cout << "Setting LM ngen = 220000" << endl;
    ngen = 220000;
  }

  int nacc_2010     = 0;
  int nacc_highmet  = 0;
  int nacc_highht   = 0;

  int nreco_2010    = 0;
  int nreco_highmet = 0;
  int nreco_highht  = 0;

  int nreco_noacc_2010    = 0;
  int nreco_noacc_highmet = 0;
  int nreco_noacc_highht  = 0;

  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  bool hasJptBtagBranch = true;

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    if( !f || f->IsZombie() ) {
      cout << "Skipping bad input file: " << currentFile->GetTitle() << endl;
      continue; //exit(1);                                                                                             
    }

    TTree *tree = (TTree*)f->Get("Events");

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
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

      //Matevz
      tree->LoadTree(z);

      cms2.GetEntry(z);

      InitBaby();

      if( !cleaning_goodDAVertexApril2011() )                        continue;
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;

      float pthat_cutoff = 30.;
      if (strcmp( prefix , "qcdpt15" ) == 0 && genps_pthat() > pthat_cutoff) {
        continue;
      }

      // skip duplicates
      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }
   
      //skip events with bad els_conv_dist 
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

      //-------------------------------
      // get acceptance for LM points
      //-------------------------------

      if( TString(prefix).Contains("LM") ){

	acc_2010_    = 0;
	acc_highmet_ = 0;
	acc_highht_  = 0;

	acc_2010_    = PassGenSelectionOS( isData ,  -1 , 300 , 8.5 );
	acc_highmet_ = PassGenSelectionOS( isData , 275 , 300 ,  -1 );
	acc_highht_  = PassGenSelectionOS( isData , 200 , 600 ,  -1 );

	//2010
	if( PassGenSelectionOS( isData ,   -1 , 300 , 8.5 ) == 1 ){
	  nacc_2010++;
	}

	//high MET
	if( PassGenSelectionOS( isData ,  275 , 300 ,  -1 ) == 1 ){
	  nacc_highmet++;
	}

	//high HT
	if( PassGenSelectionOS( isData ,  200 , 600 ,  -1 ) == 1 ){
	  nacc_highht++;
	}

      }

      //goodrun list + event cleaning
      json_ = 1;



      //find good hyps, store in v_goodHyps
      vector<unsigned int> v_goodHyps;
      v_goodHyps.clear();
      vector<unsigned int> v_goodZHyps;
      vector<float> v_weights;
      v_weights.clear();
      v_goodZHyps.clear();

      bool foundMu_ll[20];
      bool foundMu_lt[20];
      bool foundEl_ll[20];
      bool foundEl_lt[20];
   
      VofP4 goodLeptons;

      ngoodlep_ = 0;
      ngoodel_  = 0;
      ngoodmu_  = 0;
  
      if( generalLeptonVeto ){
          
        for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
          if( els_p4().at(iel).pt() < 10 )                                                 continue;
          if( !pass_electronSelection( iel , electronSelection_el_OSV3 , false , false ) ) continue;
          goodLeptons.push_back( els_p4().at(iel) );
          ngoodel_++;
          ngoodlep_++;
        }
          
        for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
          if( mus_p4().at(imu).pt() < 10 )           continue;
          if( !muonId( imu , OSGeneric_v3 ))         continue;
          goodLeptons.push_back( mus_p4().at(imu) );
          ngoodmu_++;
          ngoodlep_++;
        }  
      }

      for(unsigned int i = 0; i < hyp_p4().size(); ++i) {

	if( verbose ){
	  cout << endl << "--------------------------------" << endl;
	  cout << "hyp       " << i << endl;
	  cout << "lep ll ID " << hyp_ll_id()[i] << " pt " << hyp_ll_p4()[i].pt() << endl;
	  cout << "lep ll ID " << hyp_lt_id()[i] << " pt " << hyp_lt_p4()[i].pt() << endl;
	  cout << "mass      " << hyp_p4()[i].mass() << endl;
	}

        if( !passSUSYTrigger2011_v1( isData , hyp_type()[i] , highpt ) ) continue;

        //OS, pt > (20,10) GeV, dilmass > 10 GeV
        if( hyp_lt_id()[i] * hyp_ll_id()[i] > 0 )                               continue;
        if( TMath::Max( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < maxpt )   continue;
        if( TMath::Min( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < minpt )   continue;
        if( hyp_p4()[i].mass() < 12 )                                           continue;

        float FRweight = 1;

	if( verbose ){
	  cout << "pass OS pt > (20,10) GeV, M > 12 GeV" << endl;
	  
	  if( abs(hyp_ll_id()[i]) == 11 ){
	    cout << "Check ll electron" << endl;
	    checkElectron( hyp_ll_index()[i] );
	  }else{
	    cout << "Check ll muon" << endl;
	    checkMuon( hyp_ll_index()[i] );
	  }

	  if( abs(hyp_lt_id()[i]) == 11 ){
	    cout << "Check lt electron" << endl;
	    checkElectron( hyp_lt_index()[i] );
	  }else{
	    cout << "Check lt muon" << endl;
	    checkMuon( hyp_lt_index()[i] );
	  }

	}
        if(doFakeApp) {
          FRweight = getFRWeight(i, mufr, elfr, frmode, isData); 
          
          // FRweight < -1 --> leptons don't satisfy FO selections
          if(FRweight < -1.) continue;

	  // std::cout << "hyp " << i << " fake rate weight = " << FRweight << std::endl;
          // v_goodHyps.push_back(hypIdx);
          // v_weights.push_back(FRweight); 
        }
        else{
          
          //muon ID
          if (abs(hyp_ll_id()[i]) == 13  && !( muonId(hyp_ll_index()[i] , OSGeneric_v3 ) ) )   continue;
          if (abs(hyp_lt_id()[i]) == 13  && !( muonId(hyp_lt_index()[i] , OSGeneric_v3 ) ) )   continue;
          
          //OSV3
          if (abs(hyp_ll_id()[i]) == 11  && !( pass_electronSelection( hyp_ll_index()[i] , electronSelection_el_OSV3  ))) continue;
          if (abs(hyp_lt_id()[i]) == 11  && !( pass_electronSelection( hyp_lt_index()[i] , electronSelection_el_OSV3  ))) continue;
          
        }

	v_goodHyps.push_back( i );
        v_weights.push_back( FRweight ); // this has to be multipiled to the orig weight later on! (FRweight is == 1 for std. run)

	if( verbose ){
	  cout << "Found good hyp!" << endl;
	  cout << "--------------------------------" << endl << endl;
	}

        if( hyp_p4()[i].mass() > 76. && hyp_p4()[i].mass() < 106. ){
          v_goodZHyps.push_back(i);
        }        
      }
           
      //loop over Z hypotheses
      if( v_goodZHyps.size() > 0 ){
        
        unsigned int zhyp = selectBestZHyp(v_goodZHyps);

        float weight = 1;
        if( isData ){
          weight = 1;
        }else{
          //weight = kFactor * evt_scale1fb() * lumi * triggerSuperModelEffic( zhyp );
          weight = kFactor * evt_scale1fb() * lumi;
        }
          
        //store dilepton type in myType
        int myType = 99;
        if (hyp_type()[zhyp] == 3)                              myType = 0; // ee
        if (hyp_type()[zhyp] == 0)                              myType = 1; // mm
        if (hyp_type()[zhyp] == 1 || hyp_type()[zhyp] == 2)     myType = 2; // em
        if (myType == 99) { 
          cout << "Skipping unknown dilepton type = " << hyp_type()[zhyp] << endl;
          continue;
        }

        if( myType == 0 ) nGoodEl+=weight_;
        if( myType == 1 ) nGoodMu+=weight_;
        if( myType == 2 ) nGoodEM+=weight_;
        
        int njets = 0;
        
        for (unsigned int ijet = 0; ijet < jpts_p4().size(); ijet++) {
          
          LorentzVector vjet = jpts_p4().at(ijet) * jpts_cor().at(ijet); 
          LorentzVector vlt  = hyp_lt_p4()[zhyp];
          LorentzVector vll  = hyp_ll_p4()[zhyp];
          
          if( dRbetweenVectors(vjet, vll) < 0.4) continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4) continue;
          if( vjet.pt() < 30.          )         continue;
          if( fabs( vjet.eta() ) > 3.0 )         continue;
          if( !passesCaloJetID( vjet ) )         continue;

          if( generalLeptonVeto ){
            bool rejectJet = false;
            for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
              if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
            }
            if( rejectJet ) continue;
          }
          
          njets++;
        }
          
        if( njets > 2 ) njets = 2;

	pair<float, float> p_met_Z; //met and met phi
        
        p_met_Z = getMet( "tcMET_looper" , zhyp);
        float tcmet_looper_Z = p_met_Z.first;
        
        p_met_Z = getMet( "tcMET" , zhyp);
        float tcmet_event_Z = p_met_Z.first;
        
        p_met_Z = getMet( "pfMET" , zhyp);
        float pfmet_Z = p_met_Z.first;
        
        p_met_Z = getMet( "muCorMET" , zhyp);
        float mucormet_Z = p_met_Z.first;
        
        //        p_met_Z = getMet( "muCorJESMET" , zhyp);
        float mucorjesmet_Z = 1234. ;//p_met_Z.first;
      
        fillHistos(hdilMass_Z      , hyp_p4()[zhyp].mass() , weight , myType , njets);
        fillHistos(htcmet_event_Z  , tcmet_event_Z         , weight , myType , njets);
        fillHistos(htcmet_looper_Z , tcmet_looper_Z        , weight , myType , njets);
        fillHistos(hpfmet_Z        , pfmet_Z               , weight , myType , njets);
        fillHistos(hmucormet_Z     , mucormet_Z            , weight , myType , njets);
        fillHistos(hmucorjesmet_Z  , mucorjesmet_Z         , weight , myType , njets);
        
      }
  


      //skip events with no good hyps
      if( v_goodHyps.size() == 0 ) continue;

      //ttbar hyp disambiguation
      //returns the index of the best hypothesis in the vector of hypotheses
      unsigned int goodHyp = selectHypByHighestSumPt(v_goodHyps);
      vector<unsigned int>::const_iterator goodHyp_it = find(v_goodHyps.begin(), v_goodHyps.end(), goodHyp);
      if(goodHyp_it == v_goodHyps.end()) {
        cout << "The weight index does not correspond to the index of the best hypothesis!!!! Something is wrong" 
             << "We will quit" << endl;
        return -999;
      }
      
      //clear this vector and put in the goodHyp vector in here so we can then save some space
      //and loop over this vector below. Useful when we're not using the hypDisambiguation
      //get the index of the goodHyp in the vector of goodHyps
      unsigned int goodHyp_idx = goodHyp_it - v_goodHyps.begin();
      //get the weight of the corresponding goodHyp
      float goodHyp_weight = v_weights[goodHyp_idx];
      v_goodHyps.clear();
      v_weights.clear();
      v_goodHyps.push_back(goodHyp);
      v_weights.push_back(goodHyp_weight);

      if( v_goodHyps.size() != 1 ){
        cout << "Error, nhyps = " << v_goodHyps.size() << ", this shouldn't happen!!!!" << endl;
        exit(0);
      }

      //loop over hyps (only 1 if hyp disambiguation is performed)
  
      for(unsigned int i = 0 ; i < v_goodHyps.size() ; ++i ){

        unsigned int hypIdx = v_goodHyps.at(i);

        //store dilepton type in myType
        int myType = 99;
        if (hyp_type()[hypIdx] == 3)                              myType = 0; // ee
        if (hyp_type()[hypIdx] == 0)                              myType = 1; // mm
        if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2)   myType = 2; // em
        if (myType == 99) {
          cout << "Skipping unknown dilepton type = " << hyp_type()[hypIdx] << endl;
          continue;
        }
      
        if( hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106. ){
         
        }
      
        if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 ) nSS++;
        else nOS++;

        int nels = 0;
        int nmus  = 0;
        int ntaus = 0;
        int nleps = 0;
        
        float dilptgen = -1;
        mllgen_ = -1;
	pthat_  = -1;
	qscale_ = -1;
        if( !isData ){

	  pthat_  = genps_pthat();
	  qscale_ = genps_qScale();

          //splitting ttbar into ttdil/ttotr
          //nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
          nleps = leptonGenpCount(nels, nmus, ntaus);

	  nels_  = nels;
	  nmus_  = nmus;
	  ntaus_ = ntaus;

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


        //for tt, check if 2 leptons are from W's
        //if(strcmp(prefix,"ttdil") == 0 && ttbarconstituents(hypIdx) != 1 ) continue;
        //if(strcmp(prefix,"ttotr") == 0 && ttbarconstituents(hypIdx) == 1 ) continue;

        //cout<<"pass lepton/trigger selection"<<endl;
        // check if it's a correct genp-event (deprecated)
        //std::string prefixStr(prefix);
        //if (prefixStr == "ttdil"    && genpCountPDGId(11,13,15) != 2) continue;
        //if (prefixStr == "ttotr"    && genpCountPDGId(11,13,15) == 2) continue;
        //if (prefixStr == "DYee"     && genpCountPDGId(11)       != 2) continue;
        //if (prefixStr == "DYmm"     && genpCountPDGId(13)       != 2) continue;
        //if (prefixStr == "DYtautau" && genpCountPDGId(15)       != 2) continue;

      
        int id_lt = hyp_lt_id()[hypIdx];
        int id_ll = hyp_ll_id()[hypIdx];

        // met and jet cuts down below...

        /*
          if (dilTruthMatch) {
          //this better be in the selections.cc
          bool isTrueLepton_ll = false;
          bool isTrueLepton_lt = false;
          isTrueLepton_ll = ( (abs(hyp_ll_id()[hypIdx]) == abs(hyp_ll_mc_id()[hypIdx]) &&
          abs(hyp_ll_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
          )
          || (hyp_ll_mc_id()[hypIdx]==22 && 
          TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hyp_ll_p4()[hypIdx],hyp_ll_mc_p4()[hypIdx])) <0.05
          && abs(hyp_ll_id()[hypIdx]) == abs(hyp_ll_mc_motherid()[hypIdx])
          )
          );
          isTrueLepton_lt = ( (abs(hyp_lt_id()[hypIdx]) == abs(hyp_lt_mc_id()[hypIdx]) &&
          abs(hyp_lt_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
          )
          || (hyp_lt_mc_id()[hypIdx]==22 && 
          TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hyp_lt_p4()[hypIdx],hyp_lt_mc_p4()[hypIdx])) <0.05
          && abs(hyp_lt_id()[hypIdx]) == abs(hyp_lt_mc_motherid()[hypIdx])
          )
          );
          if (!isTrueLepton_lt && !isTrueLepton_ll) continue;
          }
        */

        // jet counting
          	
        //calojets
        VofP4 vjets_noetacut_p4;
        VofP4 vjets_p4;

        for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
          
          LorentzVector vjet = jets_p4().at(ijet) * jets_cor().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          
          if( generalLeptonVeto ){
            bool rejectJet = false;
            for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
              if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
            }
            if( rejectJet ) continue;
          }

          if( dRbetweenVectors(vjet, vll) < 0.4) continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4) continue;
          if( vjet.pt() < 30.          )         continue;
          if( !passesCaloJetID( vjet ) )         continue;

          if( fabs( vjet.eta() ) > 2.5 )         continue;
          vjets_p4.push_back( vjet );
        }
   

        //jpts
        VofP4 vjpts_p4;
        VofP4 vjpts_btag_p4;
        LorentzVector  vjpts_p4_tot(0,0,0,0);
        VofP4 vjpts_noetacut_p4;

        //njetsUp_      = 0;
        //njetsDown_    = 0;
        //sumjetptUp_   = 0.;
        //sumjetptDown_ = 0.;

	njpt_  = 0;
	htjpt_ = 0.;

        for (unsigned int ijet = 0; ijet < jpts_p4().size(); ijet++) {
          
          LorentzVector vjet     = jpts_p4().at(ijet) * jpts_corL1FastL2L3().at(ijet); 
          //LorentzVector vjetUp   = jpts_p4().at(ijet) * jpts_corL1FastL2L3().at(ijet) * 1.05; 
          //LorentzVector vjetDown = jpts_p4().at(ijet) * jpts_corL1FastL2L3().at(ijet) * 0.95; 
          LorentzVector vlt      = hyp_lt_p4()[hypIdx];
          LorentzVector vll      = hyp_ll_p4()[hypIdx];

          if( generalLeptonVeto ){
            bool rejectJet = false;
            for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
              if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
            }
            if( rejectJet ) continue;
          }
          
          if( dRbetweenVectors(vjet, vll) < 0.4) continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4) continue;
          if( !passesCaloJetID( vjet ) )         continue;

          // if( vjetUp.pt() > 30. && fabs( vjetUp.eta() ) < 3.0 ){
          //   njetsUp_++;
          //   sumjetptUp_ += vjetUp.pt();
          // }

          // if( vjetDown.pt() > 30. && fabs( vjetDown.eta() ) < 3.0 ){
          //   njetsDown_++;
          //   sumjetptDown_ += vjetDown.pt();
          // }
          
          if( vjet.pt() < 30.          )         continue;

          vjpts_noetacut_p4.push_back( vjet );

          if( fabs( vjet.eta() ) > 3.0 )         continue;

	  njpt_++;
	  htjpt_ += vjet.pt();

          vjpts_p4.push_back( vjet );
          vjpts_p4_tot += vjet;
          
          float drmin_calojet = 1000;
          int   icalo = -1;

          for (unsigned int icalojet = 0; icalojet < jets_p4().size(); icalojet++) {
            LorentzVector vcalojet = jets_p4().at(icalojet);
            float dr = dRbetweenVectors( vjet , vcalojet );
            if( dr < drmin_calojet ){
              drmin_calojet = dr;
              icalo = icalojet;
            }
          }

          if( icalo > -1 && drmin_calojet < 0.3 ){
            if( jets_simpleSecondaryVertexHighEffBJetTag().at(icalo) > 1.74 ){
              vjpts_btag_p4.push_back( vjet );
            }
          }
        }
    
        //pfjets
        VofP4 vpfjets_p4;

	npfjets_ = 0.;
	htpf_    = 0.;
	nbtags_  = 0;

	nbtagstcl_  = 0;
	nbtagstcm_  = 0;

	npfjets40_  = 0;
	htpf40_     = 0.;

	npfjetspv_  = 0;
	htpfpv_     = 0.;

	htpf25_     = 0.;
	npfjets25_  = 0;
	njets15_    = 0;

	njetsUp_   = 0;
	njetsDown_ = 0;
	htUp_      = 0.;
	htDown_    = 0.;

	int   imaxjet   = -1;
	float maxjetpt  = -1.;

	vector<int> goodjets;
	goodjets.clear();

	jetid_   = 1;
	jetid30_ = 1;

        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet      = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vjetUp    = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet) * 1.075;
          LorentzVector vjetDown  = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet) * 0.925;
          LorentzVector vlt       = hyp_lt_p4()[hypIdx];
          LorentzVector vll       = hyp_ll_p4()[hypIdx];

          if( generalLeptonVeto ){
            bool rejectJet = false;
            for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
              if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
            }
            if( rejectJet ) continue;
          }
          
          if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
          if( !passesPFJetID(ijet) ){
	    jetid_ = 0;
	    if( vjet.pt() > 30 && fabs( vjet.eta() ) < 3.0 ) jetid30_ = 0;
	    continue;
	  }

	  if( vjet.pt() > 15 && fabs( vjet.eta() ) < 3.0 ){ 
	    njets15_++;
	  }

          if( vjetUp.pt() > 30. && fabs( vjetUp.eta() ) < 3.0 ){
            njetsUp_++;
            htUp_ += vjetUp.pt();
          }

          if( vjetDown.pt() > 30. && fabs( vjetDown.eta() ) < 3.0 ){
            njetsDown_++;
            htDown_ += vjetDown.pt();
          }

          if( vjet.pt() < 30. )                    continue;

          vjets_noetacut_p4.push_back( vjet );

          if( fabs( vjet.eta() ) > 3.0 )           continue;

	  if( fabs( vjet.eta() ) < 2.5 ){
	    npfjets25_ ++;
	    htpf25_  += vjet.pt();
	  }

	  if( vjet.pt() > 40. ){
	    npfjets40_ ++;
	    htpf40_    += vjet.pt();
	  }

	  if( jetFromSignalPV( ijet , 0 , 2 ) ){
	    npfjetspv_ ++;
	    htpfpv_    += vjet.pt();
	  }

	  npfjets_++;
	  htpf_ += vjet.pt();

          vpfjets_p4.push_back( vjet );
	  goodjets.push_back(ijet);

	  if( pfjets_simpleSecondaryVertexHighEffBJetTag().at(ijet) > 1.74 ){
	    nbtags_++;
	  }

	  if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 1.7 ){
	    nbtagstcl_++;
	  }

	  if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 3.3 ){
	    nbtagstcm_++;
	  }

	  if( vjet.pt() > maxjetpt ){
	    maxjetpt = vjet.pt();
	    imaxjet  = ijet;
	  }
        }

	if( imaxjet > -1 ){ 
	  jet_ = &(pfjets_corL1FastL2L3().at(imaxjet) * pfjets_p4().at(imaxjet));

	  LorentzVector vjetraw = pfjets_p4().at(imaxjet);

	  ptjetraw_     = vjetraw.pt();
	  ptjet23_      = pfjets_cor().at(imaxjet)           * vjetraw.pt();
	  ptjetF23_     = pfjets_corL1FastL2L3().at(imaxjet) * vjetraw.pt();
	  ptjetO23_     = pfjets_corL1L2L3().at(imaxjet)     * vjetraw.pt();
	  cosphijz_     = -1 * cos( vjetraw.phi() - hyp_p4()[hypIdx].phi() );
	  
          LorentzVector vjet = pfjets_corL1FastL2L3().at(imaxjet) * pfjets_p4().at(imaxjet);
	  dphijm_ = acos(cos(vjet.phi()-evt_pfmetPhi()));
	}

	//---------------------------------
	// L1offset jets
	//---------------------------------

	htoffset_    = 0.;
	njetsoffset_ = 0;

        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_corL1L2L3().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];

          if( generalLeptonVeto ){
            bool rejectJet = false;
            for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
              if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
            }
            if( rejectJet ) continue;
          }
          
          if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
          if( fabs( vjet.eta() ) > 3.0 )           continue;
          if( !passesPFJetID(ijet) )               continue;
          if( vjet.pt() < 30. )                    continue;

	  njetsoffset_++;
	  htoffset_ += vjet.pt();

	}

 
	//---------------------------------
	// uncorrected jets
	//---------------------------------

	htuncor_    = 0.;
	njetsuncor_ = 0;

        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];

          if( generalLeptonVeto ){
            bool rejectJet = false;
            for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
              if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
            }
            if( rejectJet ) continue;
          }
          
          if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
          if( fabs( vjet.eta() ) > 3.0 )           continue;
          if( !passesPFJetID(ijet) )               continue;
          if( vjet.pt() < 15. )                    continue;



	  if( vjet.pt() < 30. )                    continue;

	  njetsuncor_++;
	  htuncor_ += vjet.pt();

	}

	ngenjets_ = 0;
	htgen_    = 0;

	if( !isData ){
	  for (unsigned int igjet = 0 ; igjet < genjets_p4().size() ; igjet++) {
	    
	    LorentzVector vgjet = genjets_p4().at(igjet);
	    LorentzVector vlt   = hyp_lt_p4()[hypIdx];
	    LorentzVector vll   = hyp_ll_p4()[hypIdx];

	    if( generalLeptonVeto ){
	      bool rejectJet = false;
	      for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
		if( dRbetweenVectors( vgjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	      }
	      if( rejectJet ) continue;
	    }
	    
	    if( vgjet.pt() < 30.                   )  continue;
	    if( fabs( vgjet.eta() ) > 3.0          )  continue;
	    if( dRbetweenVectors(vgjet, vll) < 0.4 )  continue;
	    if( dRbetweenVectors(vgjet, vlt) < 0.4 )  continue;
	    
	    ngenjets_++;
	    htgen_ += vgjet.pt();
	  }
	}

        // sumjetpt, meff calculation
        
        //calojets
        float sumjetpt_jets_p4 = 0.;
        float meff_jets_p4     = 0.;

        for(unsigned int ijet = 0; ijet < vjets_p4.size(); ijet++) {
          sumjetpt_jets_p4 += vjets_p4.at(ijet).Pt();
          meff_jets_p4 +=     vjets_p4.at(ijet).Pt();
        }

        //jpts
        float sumjetpt_jpts_p4 = 0.;
        float meff_jpts_p4     = 0.;
 
        for(unsigned int ijet = 0; ijet < vjpts_p4.size(); ijet++) {
          sumjetpt_jpts_p4 += vjpts_p4.at(ijet).Pt();
          meff_jpts_p4     += vjpts_p4.at(ijet).Pt();
        }
        
        //pfjets
        float sumjetpt_pfjets_p4 = 0.;
        float meff_pfjets_p4     = 0.;

        for(unsigned int ijet = 0; ijet < vpfjets_p4.size(); ijet++) {
          sumjetpt_pfjets_p4 += vpfjets_p4.at(ijet).Pt();
          meff_pfjets_p4     += vpfjets_p4.at(ijet).Pt();
        }

        float pt_lt  = hyp_lt_p4()[hypIdx].pt();
        float pt_ll  = hyp_ll_p4()[hypIdx].pt();
     
        //genmet stuff
        float genmet    = -9999;
        float gensumet  = -9999;
        float genmetphi = -9999;
     
        if( !isData ){
          genmet     = gen_met();
          gensumet   = gen_sumEt();
          genmetphi  = gen_metPhi();

        }

        //ttdil tcmet definition
	pair<float, float> p_met; //met and met phi
        p_met = getMet( "tcMET"    , hypIdx);
        
        float tcmet    = p_met.first;
        float tcmetphi = p_met.second;
        float tcsumet  = 1;

        meff_jets_p4      += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_jets_p4      += evt_pfmet();
        
        meff_jpts_p4      += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_jpts_p4      += evt_pfmet();
        
        meff_pfjets_p4    += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_pfjets_p4    += evt_pfmet();

        // choose which met type to use
        float theMet = -999999.;
        if      (metType == e_tcmet)    theMet = tcmet;
        else if (metType == e_muon)     theMet = evt_metMuonCorr();
        else if (metType == e_muonjes)  theMet = evt_metMuonJESCorr();
        else if (metType == e_pfmet)    theMet = evt_pfmet();
        else {
          std::cout << "UNRECOGNIZED METTYPE!  Learn to program..." << std::endl;
          exit(1);
        }

        // choose which jet type to use
        float theSumJetPt = -999999.;
        int theNJets      = -999999;
        float vecjetpt    = -999999.;
        int theNBtags     = -999999;

        if (jetType == e_JPT) {
          theSumJetPt = sumjetpt_jpts_p4;
          theNJets = vjpts_p4.size();
          vecjetpt = vjpts_p4_tot.pt();
          //theNBtags = vjpts_btag_p4.size();
        } 

        else if (jetType == e_calo) {
          theSumJetPt = sumjetpt_jets_p4;
          theNJets = vjets_p4.size();
        } 

        else if (jetType == e_pfjet) {
          theSumJetPt = sumjetpt_pfjets_p4;
          theNJets = vpfjets_p4.size();
	  theNBtags     = nbtags_;
        }
        
        else {
          std::cout << "UNRECOGNIZED JETTYPE!  Learn to program..." << std::endl;
          exit(1);
        }
  
     
        //Probably should switch all of this stuff to pfjets!!!!!!!!!!!!!!!!!!!

        // Mt2 with all the jets, no cut on eta
        float thisMt2    = 99999.;
        float thisMet    = tcmet;
        float thisMetPhi = tcmetphi;

        // Using code in the CORE
        float mt2core = MT2(evt_pfmet(), evt_pfmetPhi(), hyp_ll_p4()[hypIdx], hyp_lt_p4()[hypIdx], 0., false);

        float mt2jcore = -1.;
        if (vjets_noetacut_p4.size() > 1)
          mt2jcore = MT2J(evt_pfmet(), evt_pfmetPhi(), hyp_ll_p4()[hypIdx], hyp_lt_p4()[hypIdx], vjets_noetacut_p4);

        // Custom Mt2
        for (unsigned int i=0; i<vjets_noetacut_p4.size(); i++) {
          LorentzVector vj1 = vjets_noetacut_p4.at(i);
          for (unsigned int j=i+1; j<vjets_noetacut_p4.size(); j++) {
            LorentzVector vj2 = vjets_noetacut_p4.at(j);

            LorentzVector v1 = vj1 + hyp_lt_p4()[hypIdx];
            LorentzVector v2 = vj2 + hyp_ll_p4()[hypIdx];
            float mt2 = getMt2(v1, v2, thisMet, thisMetPhi);
            if (mt2 < thisMt2) thisMt2 = mt2;

            v1 = vj1 + hyp_ll_p4()[hypIdx];
            v2 = vj2 + hyp_lt_p4()[hypIdx];
            mt2 = getMt2(v1, v2, thisMet, thisMetPhi);
            if (mt2 < thisMt2) thisMt2 = mt2;
          }
        }
        float mt2j = thisMt2;

        m_events.insert(pair<int,int>(evt_event(), 1));


        
        // The event weight including the kFactor (scaled to 100 pb-1)
        float weight = -1.;
        if(strcmp(prefix,"LMscan") == 0){
          //weight = kFactor * sparm_xsec() * 100. / 10000.; //xsec * lumi (100/pb) / nevents (10000)
          weight = 1;
          cout << "CURRENTLY NOT SET UP TO READ IN SUSY SCAN XSEC, SETTING WEIGHT = 1" << endl;
        }else if( isData ){
          weight = 1;
        }else{
          //weight = kFactor * evt_scale1fb() * lumi * triggerSuperModelEffic( hypIdx );
          weight = kFactor * evt_scale1fb() * lumi;

          if( TString(prefix).Contains("LM") ){
            if( strcmp( prefix , "LM0" )  == 0 ) weight *= kfactorSUSY( "lm0" );
            if( strcmp( prefix , "LM1" )  == 0 ) weight *= kfactorSUSY( "lm1" );
            if( strcmp( prefix , "LM2" )  == 0 ) weight *= kfactorSUSY( "lm2" );
            if( strcmp( prefix , "LM3" )  == 0 ) weight *= kfactorSUSY( "lm3" );
            if( strcmp( prefix , "LM4" )  == 0 ) weight *= kfactorSUSY( "lm4" );
            if( strcmp( prefix , "LM5" )  == 0 ) weight *= kfactorSUSY( "lm5" );
            if( strcmp( prefix , "LM6" )  == 0 ) weight *= kfactorSUSY( "lm6" );
            if( strcmp( prefix , "LM7" )  == 0 ) weight *= kfactorSUSY( "lm7" );
            if( strcmp( prefix , "LM8" )  == 0 ) weight *= kfactorSUSY( "lm8" );
            if( strcmp( prefix , "LM9" )  == 0 ) weight *= kfactorSUSY( "lm9" );
            if( strcmp( prefix , "LM10" ) == 0 ) weight *= kfactorSUSY( "lm10");
            if( strcmp( prefix , "LM11" ) == 0 ) weight *= kfactorSUSY( "lm11");
            if( strcmp( prefix , "LM12" ) == 0 ) weight *= kfactorSUSY( "lm12");
            if( strcmp( prefix , "LM13" ) == 0 ) weight *= kfactorSUSY( "lm13");
          }
        }

        if( doFakeApp ) {  // multiply orig weight with FR hyp weight (1 for std running, FRweight for FR run)
          weight *= v_weights.at(i);
	  //cout << "weight " << weight << endl;
        }

        // This isn't quite right, and works only if both em and ppmux are in play
        // ibl: is this still applicable? Check! 100302
        if ((! strcmp(prefix, "ppMuX") || ! strcmp(prefix,"EM")) && 
            (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2)) weight *= 0.5;

        VofP4 *new_jets_p4 =  &vjets_p4;
        VofP4 *new_jpts_p4 =  &vjpts_p4;
        int new_njets =       vjets_p4.size();
        int new_njpts =       vjpts_p4.size();
        int nJetsIdx =        min(new_njpts, 2);
        //int nJptsIdx =        min(new_njpts, 2);

        //extra variables for baby ntuple
        int pass   = ( theSumJetPt > htcut && theNJets >= 2 && theMet > 50. && id_lt * id_ll < 0 ) ? 1 : 0;
        int passz  = (passZSelection ( hypIdx ) || vetoZmumuGamma( hypIdx ) ) ? 1 : 0;
        float etaZ = hyp_p4()[hypIdx].eta();
        float m0   = -9999.;
        float m12  = -9999.;

        int nvtx = 0;
    
        for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
          if(isGoodVertex(v)) ++nvtx;
        }

        int ndavtx = 0;
    
        for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
          if(isGoodDAVertex(v)) ++ndavtx;
        }
             
        
        if(strcmp(prefix,"LMscan") == 0){
          //m0  = sparm_m0();
          //m12 = sparm_m12();
          m0  = -1;
          m12 = -1;
          cout << "CURRENTLY NOT SET UP TO READ IN SUSY PARMS, SETTING M0 = M1/2 = -1" << endl;
        }

        //hyp lepton pt
        float ptll     = hyp_ll_p4()[hypIdx].pt();
        float ptlt     = hyp_lt_p4()[hypIdx].pt();
        float etall    = hyp_ll_p4()[hypIdx].eta();
        float etalt    = hyp_lt_p4()[hypIdx].eta();
        float phill    = hyp_ll_p4()[hypIdx].phi();
        float philt    = hyp_lt_p4()[hypIdx].phi();

        float ptl1    = ( ptlt > ptll ) ? ptlt  : ptll; 
        float ptl2    = ( ptlt > ptll ) ? ptll  : ptlt; 
        float etal1   = ( ptlt > ptll ) ? etalt : etall; 
        float etal2   = ( ptlt > ptll ) ? etall : etalt; 
        float phil1   = ( ptlt > ptll ) ? philt : phill; 
        float phil2   = ( ptlt > ptll ) ? phill : philt; 
        float philep  = ( ptlt > ptll ) ? hyp_lt_trk_p4()[hypIdx].phi() : hyp_ll_trk_p4()[hypIdx].phi() ;

        //find leading jet pT
        int   imax  = -9999;
        float ptmax = -9999;
              
        for(unsigned int ijet = 0 ; ijet < vjpts_p4.size() ; ++ijet){
          if(vjpts_p4.at(ijet).pt() > ptmax){
            ptmax = vjpts_p4.at(ijet).pt();
            imax  = ijet;
          }
        } 
        
        //find 2nd leading jet pT
        int imax2    = -9999;
        float ptmax2 = -9999;
              
        for(unsigned int ijet = 0 ; ijet < vjpts_p4.size() ; ++ijet){
          if(vjpts_p4.at(ijet).pt() > ptmax2 && ijet!=imax){
            ptmax2 = vjpts_p4.at(ijet).pt();
            imax2  = ijet;
          }
        } 

        //tranverse mass leading lepton & met
        float deltaPhiLepMet = fabs( philep - tcmetphi );
        if( deltaPhiLepMet > TMath::Pi() ) deltaPhiLepMet = TMath::TwoPi() - deltaPhiLepMet;
        float mt = sqrt( 2 * ( ptl1 * tcmet * (1 - cos(deltaPhiLepMet) ) ) );

        //delta phi btw leptons
        double dphilep = fabs(hyp_lt_p4()[hypIdx].phi() - hyp_ll_p4()[hypIdx].phi());
        if (dphilep > TMath::Pi()) dphilep = TMath::TwoPi() - dphilep;


	/*
	  //REPLACETOPMASS
        // calculate the top mass
        float topMass = getTopMassEstimate(d_llsol, hypIdx, vpfjets_p4, evt_pfmet(), evt_pfmetPhi());
        if(topMass != -999 && 42 != 42) std::cout<<"And top mass from exteral: "<<topMass<<std::endl;

        vector<float> topMassAllComb;
        VofP4 vjpts_p4_Comb;
        int topMassCounter = 0;
        for(int jet1 =  0; jet1 < vjpts_p4.size(); ++jet1) {
          for(int jet2 =  jet1+1; jet2 < vjpts_p4.size(); ++jet2) {
            vjpts_p4_Comb.clear();
            vjpts_p4_Comb.push_back(vjpts_p4.at(jet1));
            vjpts_p4_Comb.push_back(vjpts_p4.at(jet2));
            topMassAllComb.push_back( getTopMassEstimate(d_llsol, hypIdx, vjpts_p4_Comb, tcmet, tcmetphi) );
            //             std::cout<<
            //               "We have "<<vjpts_p4.size()<<
            //               " jets. This is combination: "<<jet1<<
            //               " with "<<jet2<<
            //               "  and topmass "<<topMassAllComb.at(topMassCounter)<<
            //               std::endl;
            ++topMassCounter;
          }
        }
	*/

        //get various met types
        
        float tcmet_35X = -9999;

        p_met = getMet( "tcMET_looper"    , hypIdx);
        float tcmet_looper = p_met.first;

        p_met = getMet( "tcMET"    , hypIdx);
        float tcmet_event = p_met.first;

        p_met = getMet( "pfMET"    , hypIdx);
        float pfmet    = p_met.first;
	float pfmetphi = p_met.second;

        pair<float, float> pfmetUp   = ScaleMET( p_met , hyp_p4().at(hypIdx) , 1.075 );
        pair<float, float> pfmetDown = ScaleMET( p_met , hyp_p4().at(hypIdx) , 0.925 );
        pair<float, float> pfmetTest = ScaleMET( p_met , hyp_p4().at(hypIdx) , 1.000 );

        p_met = getMet( "muCorMET"    , hypIdx);
        float mucormet = p_met.first;

	pair<float, float> p_tcmet; //met and met phi
        p_tcmet = getMet( "tcMET"    , hypIdx);
          
        pair<float, float> metUp   = ScaleMET( p_tcmet , hyp_p4().at(hypIdx) , 1.075 );
        pair<float, float> metDown = ScaleMET( p_tcmet , hyp_p4().at(hypIdx) , 0.925 );
        pair<float, float> metTest = ScaleMET( p_tcmet , hyp_p4().at(hypIdx) , 1.000 );

        tcmet00_ = smearMet( theMet , theSumJetPt , 1.00 );
        tcmet10_ = smearMet( theMet , theSumJetPt , 1.10 );
        tcmet20_ = smearMet( theMet , theSumJetPt , 1.20 );
        tcmet30_ = smearMet( theMet , theSumJetPt , 1.30 );
        tcmet40_ = smearMet( theMet , theSumJetPt , 1.40 );
        tcmet50_ = smearMet( theMet , theSumJetPt , 1.50 );

        //        p_met = getMet( "muCorJESMET"    , hypIdx);
        float mucorjesmet = 1234. ;//p_met.first;

        //fill tree for baby ntuple 
        if(g_createTree){

          costhetaweight_ = -3;
          //if(strcmp(prefix,"ttdil") == 0 )
          //  costhetaweight_ = getCosThetaStarWeight();
        
          mullgen_      = foundMu_ll[hypIdx] ? 1 : 0;
          multgen_      = foundMu_lt[hypIdx] ? 1 : 0;
          //mull_         = (abs(hyp_ll_id()[hypIdx]) == 13  && (! muonId(hyp_ll_index()[hypIdx] , NominalTTbarV2 ) ) ) ? 0 : 1;
          //mult_         = (abs(hyp_lt_id()[hypIdx]) == 13  && (! muonId(hyp_lt_index()[hypIdx] , NominalTTbarV2 ) ) ) ? 0 : 1;
          nlep_         = nels + nmus;
          tcmet_looper_ = tcmet_looper;
          tcmet_35X_    = tcmet_35X;
          tcmet_event_  = tcmet_event;
          tcmetUp_      = metUp.first;
          tcmetDown_    = metDown.first;
          tcmetTest_    = metTest.first;
          pfmetUp_      = pfmetUp.first;
          pfmetDown_    = pfmetDown.first;
          pfmetTest_    = pfmetTest.first;
          pfmet_        = pfmet;
          pfmetsig_     = -999;//evt_pfmetSignificance();
          pfmetphi_     = pfmetphi;
          mucormet_     = mucormet;
          mucorjesmet_  = mucorjesmet;
          genmet_       = genmet;                       //generated met from neutrinos/LSP
          weight_       = weight;                       //event weight
          //smeff_        = isData ? 1 : triggerSuperModelEffic( hypIdx ); //trigger supermodel efficiency
          smeff_        = 1;
          proc_         = getProcessType(prefix);       //integer specifying sample
          topmass_      = -999;//topMass;                      //topepton mass //REPLACE TOPMASS
          dilmass_      = hyp_p4()[hypIdx].mass();      //dilepton mass
          dilpt_        = hyp_p4()[hypIdx].pt();        //dilepton pT
          dileta_       = hyp_p4()[hypIdx].eta();       //dilepton eta
          dildphi_      = dphilep;                      //dilepton delta phi
          tcmet_        = tcmet;                        //tcmet
          tcsumet_      = tcsumet;                      //tcsumet
          tcmetphi_     = tcmetphi;                     //tcmetphi
          sumjetpt_     = theSumJetPt;                  //scalar sum jet pt
          mt2_          = mt2core;                      //mt2 leptonic
          mt2j_         = mt2j;                         //mt2 with jets
          mt2jcore_     = mt2jcore;                     //mt2 with jets (core)
          njets_        = theNJets;                     //njets w pt>30 and |eta|<2.5
          nvtx_         = nvtx;                         //number of good vertices in this event
          ndavtx_       = ndavtx;                       //number of good DA vertices in this event
          vecjetpt_     = vecjetpt;                     //vector sum jet pt
          pass_         = pass;                         //pass kinematic cuts
          passz_        = passz;                        //pass Z selection
          m0_           = m0;                           //mSUGRA m0
          m12_          = m12;                          //mSUGRA m1/2
          ptl1_         = ptl1;                         //highest pT lepton
          ptl2_         = ptl2;                         //2nd highest pT lepton
          ptj1_         = ptmax;                        //leading jet
          ptj2_         = ptmax2;                       //2nd leading jet
          etal1_        = etal1;                        //highest pT lepton
          etal2_        = etal2;                        //2nd highest pT lepton
          phil1_        = phil1;                        //highest phi lepton
          phil2_        = phil2;                        //2nd highest phi lepton
          meff_         = meff_pfjets_p4;               //effective mass
          mt_           = mt;                           //transverse mass of leading lepton+met
	  y_		= theMet / sqrt( theSumJetPt ); //y=MET/sqrt(HT)
	  ht_		= theSumJetPt;                  //HT
          strcpy(dataset_, cms2.evt_dataset().Data());  //dataset name
          run_          = evt_run();                    //run
          lumi_         = evt_lumiBlock();              //lumi
          event_        = evt_event();                  //event
	  ndavtxweight_ = vtxweight(isData,true);
	  hbhe_         = evt_hbheFilter();

          k_				= 1;
          if( strcmp( prefix , "LM0"  )  == 0 ) k_ = kfactorSUSY( "lm0"  );
          if( strcmp( prefix , "LM1"  )  == 0 ) k_ = kfactorSUSY( "lm1"  );
          if( strcmp( prefix , "LM2"  )  == 0 ) k_ = kfactorSUSY( "lm2"  );
          if( strcmp( prefix , "LM3"  )  == 0 ) k_ = kfactorSUSY( "lm3"  );
          if( strcmp( prefix , "LM4"  )  == 0 ) k_ = kfactorSUSY( "lm4"  );
          if( strcmp( prefix , "LM5"  )  == 0 ) k_ = kfactorSUSY( "lm5"  );
          if( strcmp( prefix , "LM6"  )  == 0 ) k_ = kfactorSUSY( "lm6"  );
          if( strcmp( prefix , "LM7"  )  == 0 ) k_ = kfactorSUSY( "lm7"  );
          if( strcmp( prefix , "LM8"  )  == 0 ) k_ = kfactorSUSY( "lm8"  );
          if( strcmp( prefix , "LM9"  )  == 0 ) k_ = kfactorSUSY( "lm9"  );
          if( strcmp( prefix , "LM10" )  == 0 ) k_ = kfactorSUSY( "lm10" );
          if( strcmp( prefix , "LM11" )  == 0 ) k_ = kfactorSUSY( "lm11" );
          if( strcmp( prefix , "LM12" )  == 0 ) k_ = kfactorSUSY( "lm12" );

          float dzcut  = 0.1; // dz(trk,vtx) requirement
          float etacut = 3.0; // neutral PFCandidate eta requirement

          //met built from charged PFCandidates
          pair<float, float> p_trkmet    = PFCandidateMET( 0, hypIdx, goodjets, dzcut, 1.e10  , etacut ,  true , false );
          trkmet_      = p_trkmet.first;
          trkmetphi_   = p_trkmet.second;
          trkmetproj_  = projectedMET( trkmet_ , trkmetphi_ , hypIdx );

          //met built from charged PFCandidates and neutral PFCandidates pt > 4 GeV
          pair<float, float> p_trkmet4   = PFCandidateMET( 0, hypIdx, goodjets, dzcut,    4.  , etacut ,  true , false );
          trkmet4_      = p_trkmet.first;
          trkmet4phi_   = p_trkmet.second;
          trkmet4proj_  = projectedMET( trkmet4_ , trkmet4phi_ , hypIdx );

          //met built from charged PFCandidates and neutral PFCandidates pt > 8 GeV
          pair<float, float> p_trkmet8   = PFCandidateMET( 0, hypIdx, goodjets, dzcut,    8.  , etacut ,  true , false );
          trkmet8_      = p_trkmet8.first;
          trkmet8phi_   = p_trkmet8.second;
          trkmet8proj_  = projectedMET( trkmet8_ , trkmet8phi_ , hypIdx );

          //met built from jets and charged PFCandidates
          pair<float, float> p_trkjetmet = PFCandidateMET( 0, hypIdx, goodjets, dzcut,  1.e10 , etacut ,  true ,  true );
          trkjetmet_      = p_trkjetmet.first;
          trkjetmetphi_   = p_trkjetmet.second;
          trkjetmetproj_  = projectedMET( trkjetmet_ , trkjetmetphi_ , hypIdx );          

	  //--------------------------
	  // leading lepton = ll
	  //--------------------------

	  int index1 = -1;
	  int index2 = -1;

          if( hyp_ll_p4().at(hypIdx).pt() > hyp_lt_p4().at(hypIdx).pt() ){

	    index1 = hyp_ll_index()[hypIdx];
	    index2 = hyp_lt_index()[hypIdx];

            lep1_ = &hyp_ll_p4().at(hypIdx);
            lep2_ = &hyp_lt_p4().at(hypIdx);
	    id1_  = hyp_ll_id()[hypIdx];
	    id2_  = hyp_lt_id()[hypIdx];

	  }

	  //--------------------------
	  // leading lepton = lt
	  //--------------------------

	  else{

	    index1 = hyp_lt_index()[hypIdx];
	    index2 = hyp_ll_index()[hypIdx];

            lep1_ = &hyp_lt_p4().at(hypIdx);
            lep2_ = &hyp_ll_p4().at(hypIdx);
	    id1_  = hyp_lt_id()[hypIdx];
	    id2_  = hyp_ll_id()[hypIdx];

	  }

	  if( !isData ){
	    w1_          = leptonOrTauIsFromW( index1 , id1_ , isLM );
	    w2_          = leptonOrTauIsFromW( index2 , id2_ , isLM );
	  }

	  if( abs(id1_) == 11 ){
	    iso1_   = electronIsolation_rel   ( index1 , true ); //truncated
	    isont1_ = electronIsolation_rel_v1( index1 , true ); //non-truncated
	  }
	  else if( abs(id1_) == 13 ){
	    iso1_   = muonIsoValue( index1 , true  ); //truncated 
	    isont1_ = muonIsoValue( index1 , false ); //non-truncated
	  }
	  
	  if( abs(id2_) == 11 ){
	    iso2_   = electronIsolation_rel   ( index2 , true ); //truncated
	    isont2_ = electronIsolation_rel_v1( index2 , true ); //non-truncated
	  }
	  else if( abs(id2_) == 13 ){
	    iso2_   = muonIsoValue( index2 , true  ); //truncated 
	    isont2_ = muonIsoValue( index2 , false ); //non-truncated
	  }
	  
	  dilep_   = &hyp_p4().at(hypIdx);

          leptype_ = -1;
          if (hyp_type()[hypIdx] == 3) leptype_ = 0; // ee
          if (hyp_type()[hypIdx] == 0) leptype_ = 1; // mm
          if (hyp_type()[hypIdx] == 1) leptype_ = 2; // em
          if( hyp_type()[hypIdx] == 2) leptype_ = 2; // em
                
	  trgeff_ = 1;
	  if(!isData){
	    if( leptype_ == 0 ) trgeff_ = 1.00;
	    if( leptype_ == 1 ) trgeff_ = 0.90;
	    if( leptype_ == 2 ) trgeff_ = 0.95;
	  }

          outTree->Fill();
        }

	if( TString(prefix).Contains("LM") ){

	  //2010 signal region
	  if( PassGenSelectionOS( isData ,   -1 , 300 , 8.5 ) == 1 ){
	    if( passz == 0 && npfjets_ >= 2 && htpf_ > 300. && pfmet_ > 50. && y_ > 8.5 ) 
	      nreco_2010++;
	  }else{
	    if( passz == 0 && npfjets_ >= 2 && htpf_ > 300. && pfmet_ > 50. && y_ > 8.5 ) 
	      nreco_noacc_2010++;
	  }

	  //high MET signal region
	  if( PassGenSelectionOS( isData ,  275 , 300 ,  -1 ) == 1 ){
	    if( passz == 0 && npfjets_ >= 2 && htpf_ > 300. && pfmet_ > 275. ) 
	      nreco_highmet++;
	  }else{
	    if( passz == 0 && npfjets_ >= 2 && htpf_ > 300. && pfmet_ > 275. ) 
	      nreco_noacc_highmet++;

	    // if( PassGenSelectionOS( isData ,  275 , 300 ,  -1 ) == -1 ){
	    //   dumpDocLines(true);
	    // }

	  }

	  //high HT signal region
	  if( PassGenSelectionOS( isData ,  200 , 600 ,  -1 ) == 1){
	    if( passz == 0 && npfjets_ >= 2 && htpf_ > 600. && pfmet_ > 200. ) 
	      nreco_highht++;
	  }else{
	    if( passz == 0 && npfjets_ >= 2 && htpf_ > 600. && pfmet_ > 200. ) 
	      nreco_noacc_highht++;
	  }
	}

	if     ( leptype_ == 0 ) neetot += weight;
	else if( leptype_ == 1 ) nmmtot += weight;
	else if( leptype_ == 2 ) nemtot += weight;

        //selection (continue statements)-------------------------------
        if(!g_useBitMask){
        
          // all these cuts are ok for the FR application

          if (theSumJetPt < htcut)   continue;
          if (theNJets < 2)          continue;
          if (theMet < 50.)          continue;
      
          //if( myType == 2 && tcmet < 20 )                    continue; 
          //if( ( myType == 0 || myType == 1 ) && tcmet < 30 ) continue; 
         
          
          if (zveto == e_standard) {
          
            //veto same-flavor OS dileptons in Z mass window
            if ( ( hyp_type()[hypIdx] == 3 || hyp_type()[hypIdx] == 0 ) 
                 && hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106.) continue;

            if( vetoZmumuGamma( hypIdx ) ) continue;     
          }
        
          else if(zveto == e_allzveto){
            //veto dileptons in Z mass window regardless of flavor
            if (hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106.)   continue;

            if( vetoZmumuGamma( hypIdx ) ) continue;  
          }
        
          else if(zveto == e_nozveto){
            //no Z veto
          }
        
          else if(zveto == e_selectz){
            //require same-flavor OS dileptons in Z mass window
            if ( hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2 || 
                 hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.) continue;
          }
        
          else{
            cout<<"UNRECOGNIZED ZVETO"<<endl;
            exit(0);
          }
          

          fillHistos(hdilMass,  hyp_p4()[hypIdx].mass()  , weight, myType, nJetsIdx);
          fillHistos(htcmet  ,  tcmet            , weight, myType, nJetsIdx);
          fillHistos(hmetmuon,  evt_metMuonCorr(), weight, myType, nJetsIdx);
          fillHistos(hsumJetPt, sumjetpt_jets_p4 , weight, myType, nJetsIdx);
          fillHistos(hsumJptPt, sumjetpt_jpts_p4 , weight, myType, nJetsIdx);
          hnJet[myType]   ->Fill(new_njets,     weight);
          hnJet[3]        ->Fill(new_njets,     weight);
          hnJpt[myType]   ->Fill(new_njpts,     weight);
          hnJpt[3]        ->Fill(new_njpts,     weight);
          hnBtagJpt[myType]   ->Fill( theNBtags,     weight);
          hnBtagJpt[3]        ->Fill( theNBtags,     weight);

          
          //Fill susyscan histos
          if(strcmp(prefix,"LMscan") == 0){
                    
            //cout << " m0 "     << sparm_m0()  << " index " << getIndexFromM0(m0) 
            //     << " m1/2 "   << sparm_m12() << " index " << getIndexFromM12(m12)
            //     << " weight " << weight << endl;
                    
            fillUnderOverFlow( hsusydilPt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , hyp_p4()[hypIdx].pt(), weight);
            fillUnderOverFlow( hsusytcmet[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , tcmet, weight);
            hsusy_met_sumjetpt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)]->Fill(sumjetpt_jets_p4 , tcmet/sqrt(tcsumet), weight);
                    
          }
        }

        //selection (bitmask)-------------------------------------------

        if(g_useBitMask){

          const int ncut = 5;
          bool cut[ncut];
          for(int ic=0;ic<ncut;ic++)cut[ic]=false;
          unsigned cutbit=0;

          cut[0] = (id_lt * id_ll < 0); 
        
          if (zveto == e_standard) {
            //veto same-flavor OS dileptons in Z mass window
            cut[1] = (hyp_type()[hypIdx]==3 || hyp_type()[hypIdx]==0) ? 
              (hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.) : true; 

            if( vetoZmumuGamma( hypIdx ) ) cut[1] = false;
          }
          else if(zveto == e_allzveto){
            //veto OS dileptons in Z mass window regardless of flavor
            cut[1] = (hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.);

            if( vetoZmumuGamma( hypIdx ) ) cut[1] = false;
          }
          else if(zveto == e_nozveto){
            //no Z-veto
            cut[1] = true;
          }
          else if(zveto == e_selectz){
            //require same-flavor OS dileptons in Z mass window
            cut[1] = ( ( hyp_type()[hypIdx] == 3 || hyp_type()[hypIdx] == 0 ) && hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106. );
          }

          cut[2] = theMet > 50.;
          cut[3] = theSumJetPt > htcut;
          cut[4] = theNJets    > 1;

          for (int icut=0;icut<ncut;++icut) {
            if (cut[icut]) cutbit+=(1<<icut);
          }

          //Fill N-1 histos
          //NOTE: if(nkcut(cutbit,ncut,X)) means: if event passes all cuts *except* cut X,
          //where X refers to cut[X] = 'some requirement' above

          if(nkcut(cutbit,ncut,1)){//dilepton mass
            fillHistos(hdilMass,  hyp_p4()[hypIdx].mass()  , weight, myType, nJetsIdx);
          }
          if(nkcut(cutbit,ncut)){//dilepton mass
            fillHistos(hdilMass_cut,  hyp_p4()[hypIdx].mass()  , weight, myType, nJetsIdx);
          }
          if(nkcut(cutbit,ncut,2)){//met
            fillHistos(htcmet  ,  tcmet            , weight, myType, nJetsIdx);
            fillHistos(hmetmuon,  evt_metMuonCorr(), weight, myType, nJetsIdx);

            fillHistos(hetaZ_tcmet,          fabs(etaZ) , tcmet,               weight, myType, nJetsIdx);
            fillHistos(hetaZ_tcmetsqrtsumet, fabs(etaZ) , tcmet/sqrt(tcsumet), weight, myType, nJetsIdx);
            fillHistos(hetaZ_tcmetsumet,     fabs(etaZ) , tcmet/tcsumet,       weight, myType, nJetsIdx);

            float x = tcmet;
            float y = hyp_p4()[hypIdx].pt();

            //full region
            fillHistos(hmt2j_all,      mt2jcore, weight, myType, nJetsIdx);
            fillHistos(hmet_dilpt_all, x, y,     weight, myType, nJetsIdx);

            //signal region
            if( x > 100 && y > 50 && x > y){
              fillHistos(hmt2j_signal,      mt2jcore, weight, myType, nJetsIdx);
              fillHistos(hmet_dilpt_signal, x, y,     weight, myType, nJetsIdx);
            }
            //control region
            if( y > 100 && x > 50 && y > x){
              fillHistos(hmt2j_control,      mt2jcore, weight, myType, nJetsIdx);
              fillHistos(hmet_dilpt_control, x, y,     weight, myType, nJetsIdx);
            }

          }
          if(nkcut(cutbit,ncut,3)){//sumjetpt
            fillHistos(hsumJetPt, sumjetpt_jets_p4 , weight, myType, nJetsIdx);
            fillHistos(hsumJptPt, sumjetpt_jpts_p4 , weight, myType, nJetsIdx);
          }

          if(nkcut(cutbit,ncut)){//sumjetpt
            fillHistos(hsumJetPt_cut, sumjetpt_jets_p4 , weight, myType, nJetsIdx);
            fillHistos(hsumJptPt_cut, sumjetpt_jpts_p4 , weight, myType, nJetsIdx);
          }

          if(nkcut(cutbit,ncut,4)){//njets
            hnJet[myType]   ->Fill(new_njets,     weight);
            hnJet[3]        ->Fill(new_njets,     weight);
            hnJpt[myType]   ->Fill(new_njpts,     weight);
            hnJpt[3]        ->Fill(new_njpts,     weight);
            hnBtagJpt[myType]   ->Fill( theNBtags,     weight);
            hnBtagJpt[3]        ->Fill( theNBtags,     weight);
          }
                  
          if(nkcut(cutbit,ncut,2,3)){
            
            fillHistos( habcd_tprof_nopresel,   theSumJetPt , tcmet/sqrt(theSumJetPt),               myType, nJetsIdx);
            fillHistos( habcd_nopresel,         theSumJetPt , tcmet/sqrt(theSumJetPt),               weight, myType, nJetsIdx);

            fillHistos(hsumJetPt_tcmet,          theSumJetPt , tcmet,               weight, myType, nJetsIdx);
            fillHistos(hsumJetPt_tcmetsqrtsumet, theSumJetPt , tcmet/sqrt(tcsumet), weight, myType, nJetsIdx);
            fillHistos(hsumJetPt_tcmetsumet,     theSumJetPt , tcmet/tcsumet,       weight, myType, nJetsIdx);

            fillHistos(hsumJetPt_tcmetsqrtsumet_prof, theSumJetPt , tcmet/sqrt(tcsumet), myType, nJetsIdx);
            fillHistos(htcsumet_tcmet_prof,           tcsumet ,          tcmet,               myType, nJetsIdx);
            fillHistos(hgensumet_genmet_prof,         gensumet   ,       genmet,              myType, nJetsIdx);

          }
                  
          //Fill susyscan histos
          if(strcmp(prefix,"LMscan") == 0){
          
            //cout << " m0 "     << sparm_m0()  << " index " << getIndexFromM0(m0) 
            //     << " m1/2 "   << sparm_m12() << " index " << getIndexFromM12(m12)
            //     << " weight " << weight << endl;
          
            if(nkcut(cutbit,ncut)){
              fillUnderOverFlow( hsusydilPt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , hyp_p4()[hypIdx].pt(), weight);
              fillUnderOverFlow( hsusytcmet[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , tcmet, weight);
            }

            if(nkcut(cutbit,ncut,2,3)){
              hsusy_met_sumjetpt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)]->Fill(theSumJetPt , tcmet/sqrt(tcsumet), weight);
            }
          }

          if(nkcut(cutbit,ncut,2)){//met
            fillHistos( hdpfmet_genmet , genmet , ( pfmet - genmet ) / genmet , 1 , myType , nJetsIdx );
            fillHistos( tdpfmet_genmet , genmet , ( pfmet - genmet ) / genmet ,     myType , nJetsIdx );
            
            fillHistos( hdtcmetevent_genmet , genmet , ( tcmet_event - genmet ) / genmet , 1 , myType , nJetsIdx );
            fillHistos( tdtcmetevent_genmet , genmet , ( tcmet_event - genmet ) / genmet ,     myType , nJetsIdx );
            
            fillHistos( hdtcmetlooper_genmet , genmet , ( tcmet_looper - genmet ) / genmet , 1 , myType , nJetsIdx );
            fillHistos( tdtcmetlooper_genmet , genmet , ( tcmet_looper - genmet ) / genmet ,     myType , nJetsIdx );
            
            fillHistos( hdmucormet_genmet , genmet , ( mucormet - genmet ) / genmet , 1 , myType , nJetsIdx );
            fillHistos( tdmucormet_genmet , genmet , ( mucormet - genmet ) / genmet ,     myType , nJetsIdx );
            
            fillHistos( hdmucorjesmet_genmet , genmet , ( mucorjesmet - genmet ) / genmet , 1 , myType , nJetsIdx );
            fillHistos( tdmucorjesmet_genmet , genmet , ( mucorjesmet - genmet ) / genmet ,     myType , nJetsIdx );
          }
        
          if(!nkcut(cutbit,ncut)) continue; //continue if event doesn't pass ALL cuts
        }

        if( ZVetoGeneral() ){
          //cout << "veto Z event!!" << endl;
          //printEventInfo();
        }

        //-------------------------------------------------------------
        // Lots of histograms
        //-------------------------------------------------------------
              
        if     ( myType == 0 ) nee += weight;
        else if( myType == 1 ) nmm += weight;
        else if( myType == 2 ) nem += weight;
        else{ cout << "UNKNOWN TYPE " << myType << endl; exit(0); }

        hyield->Fill(0.5,        weight);
        hyield->Fill(1.5+myType, weight);

        //hyield_weight->Fill(0.5,        weight / triggerSuperModelEffic( hypIdx ));
        //hyield_weight->Fill(1.5+myType, weight / triggerSuperModelEffic( hypIdx ));
        hyield_weight->Fill(0.5,        weight );
        hyield_weight->Fill(1.5+myType, weight );

        if( theSumJetPt > 300. && theMet / sqrt( theSumJetPt ) > 8.5 ){
          hyieldsig->Fill(0.5,        weight);
          hyieldsig->Fill(1.5+myType, weight);
        }

        hyield_unweighted->Fill(0.5);
        hyield_unweighted->Fill(1.5+myType);

        if( theSumJetPt > 300. && theMet / sqrt( theSumJetPt ) > 8.5 ){
          hyieldsig_unweighted->Fill(0.5);
          hyieldsig_unweighted->Fill(1.5+myType);
        }

        fillHistos( htcmet_sqrtht,               tcmet/sqrt(theSumJetPt),               weight, myType, nJetsIdx);
        fillHistos( habcd,         theSumJetPt , tcmet/sqrt(theSumJetPt),               weight, myType, nJetsIdx);
        if( theSumJetPt > 150 && tcmet/sqrt(theSumJetPt) > 4.5 )
          fillHistos( habcd_tprof,   theSumJetPt , tcmet/sqrt(theSumJetPt),               myType, nJetsIdx);

        fillHistos(hetaz, fabs(etaZ) , weight, myType, nJetsIdx);
        fillHistos(hmt2jcore, mt2jcore, weight, myType, nJetsIdx);
        fillHistos(hmt2core,  mt2core,  weight, myType, nJetsIdx);
        fillHistos(hmt2j, mt2j, weight, myType, nJetsIdx);
        fillHistos(hmt,   mt,   weight, myType, nJetsIdx);
        //fillHistos(hsumJetPt, sumjetpt_jets_p4, weight, myType, nJetsIdx);
        fillHistos(hDtcmetgenmetVsumJetPt, theSumJetPt, tcmet-genmet, weight, myType, nJetsIdx);
        fillHistos(hDmetmuonjesgenmetVsumJetPt, theSumJetPt, evt_metMuonJESCorr()-genmet, weight, myType, nJetsIdx);
        fillHistos(hmeffJet, meff_jets_p4, weight, myType, nJetsIdx);
        //fillHistos(hsumJptPt, sumjetpt_jpts_p4, weight, myType, nJetsIdx);
        fillHistos(hmeffJPT, meff_jpts_p4, weight, myType, nJetsIdx);

        //fillHistos(hetaZ_tcmet, etaZ, tcmet, weight, myType, nJetsIdx);
        //fillHistos(hetaZ_tcmetsqrtsumet, etaZ, tcmet*sqrt(theSumJetPt), weight, myType, nJetsIdx);
        //fillHistos(hetaZ_tcmetsumet, etaZ, tcmet*theSumJetPt, weight, myType, nJetsIdx);
                

        // jet count
        //hnJet[myType]->Fill(new_njets, weight);
        //hnJet[3]->Fill(new_njets, weight);
        //hnJpt[myType]->Fill(new_njpts, weight);
        //hnJpt[3]->Fill(new_njpts, weight);
        //hnHypJet[myType]->Fill(new_hyp_njets, weight);
        //hnHypJet[3]->Fill(new_hyp_njets, weight);

        // lepton Pt
        if (abs(id_lt) == 11)  fillHistos(helePt, pt_lt, weight, myType, nJetsIdx);
        if (abs(id_ll) == 11)  fillHistos(helePt, pt_ll, weight, myType, nJetsIdx);
        if (abs(id_lt) == 13)  fillHistos(hmuPt,  pt_lt, weight, myType, nJetsIdx);
        if (abs(id_ll) == 13)  fillHistos(hmuPt,  pt_ll, weight, myType, nJetsIdx);

        if( pt_ll > ptlt ){
          fillHistos(hminLepPt,   pt_lt,                     weight, myType, nJetsIdx);
          fillHistos(hminLepEta,  hyp_lt_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
          fillHistos(hmaxLepPt,   pt_ll,                     weight, myType, nJetsIdx);
          fillHistos(hmaxLepEta,  hyp_ll_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
        }else{
          fillHistos(hminLepPt,   pt_ll,                     weight, myType, nJetsIdx);
          fillHistos(hminLepEta,  hyp_ll_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
          fillHistos(hmaxLepPt,   pt_lt,                     weight, myType, nJetsIdx);
          fillHistos(hmaxLepEta,  hyp_lt_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
        }

        // lepton Phi
        if (abs(id_lt) == 11) helePhi[myType][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[myType][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[myType][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[myType][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 11) helePhi[myType][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[myType][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[myType][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[myType][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 11) helePhi[3][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[3][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[3][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[3][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 11) helePhi[3][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[3][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[3][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[3][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);

        // dilepton mass
        //fillHistos(hdilMass, hyp_p4()[hypIdx].mass(), weight, myType, nJetsIdx);
        // top mass
	/* REPLACE TOPMASS
        if(topMass > 0.) fillHistos(htopMass, topMass, weight, myType, nJetsIdx);

        for(int imass = 0; imass < topMassAllComb.size(); ++imass) {
          if( topMassAllComb.at(imass) > 0. ) fillHistos(htopMassAllComb, topMassAllComb.at(imass), weight, myType, nJetsIdx);
        }
	*/

        // delta phi btw leptons
        double dphi = fabs(hyp_lt_p4()[hypIdx].phi() - hyp_ll_p4()[hypIdx].phi());
        if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
        fillHistos(hdphiLep, dphi, weight, myType, nJetsIdx);
        
        float dr_lep = dRbetweenVectors( hyp_lt_p4()[hypIdx] , hyp_ll_p4()[hypIdx] );
        fillHistos(hdrLep, dr_lep , weight, myType, nJetsIdx);

        // dphill vs mll, i.e. the 2d correlation between the previous two variables
        fillHistos(hdphillvsmll, hyp_p4()[hypIdx].mass(), dphi, weight, myType, nJetsIdx);

        // lepton Eta
        if (abs(id_lt) == 11)  fillHistos(heleEta, hyp_lt_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
        if (abs(id_ll) == 11)  fillHistos(heleEta, hyp_ll_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
        if (abs(id_lt) == 13)  fillHistos(hmuEta,  hyp_lt_p4()[hypIdx].eta(), weight, myType, nJetsIdx);
        if (abs(id_ll) == 13)  fillHistos(hmuEta,  hyp_ll_p4()[hypIdx].eta(), weight, myType, nJetsIdx);

        // dilepton pt
        fillHistos(hdilPt, hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        //dilepton pt with z-veto applied
        if(!passZSelection(hypIdx))
          fillHistos(hdilPt_zveto, hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);
    
        // smeared dilepton pt
        float sigma = returnSigma(theSumJetPt, metType);
        float bias  = returnBias(theSumJetPt, metType);
        float smear = random3_->Gaus(0, sigma);
        float dilPtSmeared = bias * hyp_p4()[hypIdx].pt() + smear;
        fillHistos(hdilPtSmeared, dilPtSmeared, weight, myType, nJetsIdx);

        // Gen Met and Met Phi
        fillHistos(hgenmet, genmet, weight, myType, nJetsIdx);
        fillHistos(hgenmetPhi, genmetphi , weight, myType, nJetsIdx);

        // Met and Met phi
        //fillHistos(hmetmuon, evt_metMuonCorr(), weight, myType, nJetsIdx);
        fillHistos(hmetmuonPhi, evt_metMuonCorrPhi(), weight, myType, nJetsIdx);

        // pat Met and Met phi
        fillHistos(hmetmuonjes, evt_metMuonJESCorr(), weight, myType, nJetsIdx);
        fillHistos(hmetmuonjesPhi, evt_metMuonJESCorrPhi(), weight, myType, nJetsIdx);

        // tc Met and Met phi
        //fillHistos(htcmet, tcmet, weight, myType, nJetsIdx);
        fillHistos(htcmetPhi, tcmetphi, weight, myType, nJetsIdx);

        // pf Met and Met phi
        fillHistos(hpfmet, evt_pfmet(), weight, myType, nJetsIdx);
        fillHistos(hpfmetPhi, evt_pfmetPhi(), weight, myType, nJetsIdx);

        // Met vs dilepton Pt
        fillHistos(hmetmuonVsDilepPt, evt_metMuonCorr(), hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        // pat  Met vs dilepton Pt
        fillHistos(hmetmuonjesVsDilepPt, evt_metMuonJESCorr(), hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        // tc  Met vs dilepton Pt
        fillHistos(htcmetVsDilepPt, tcmet, hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        // Met over dilepton Pt vs deltaphi btw the two
        double dphi2 = fabs(hyp_p4()[hypIdx].phi() - evt_metMuonCorrPhi() );
        if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
        dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
        fillHistos(hmetmuonOverPtVsDphi, evt_metMuonCorr()/hyp_p4()[hypIdx].pt(), dphi2, weight, myType, nJetsIdx);

        // pat Met over dilepton Pt vs deltaphi btw the two
        dphi2 = fabs(hyp_p4()[hypIdx].phi() - evt_metMuonJESCorrPhi());
        if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
        dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
        fillHistos(hmetmuonjesOverPtVsDphi, evt_metMuonJESCorr()/hyp_p4()[hypIdx].pt(), dphi2, weight, myType, nJetsIdx);

        // tc Met over dilepton Pt vs deltaphi btw the two
        dphi2 = fabs(hyp_p4()[hypIdx].phi() - tcmetphi);
        if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
        dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
        fillHistos(htcmetOverPtVsDphi, tcmet/hyp_p4()[hypIdx].pt(), dphi2, weight, myType, nJetsIdx);

        // Make a vector of jets sorted by pt and fill jet histograms
        if (new_njets > 0) {
          //vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > my_jets_p4(*new_jets_p4);
          vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > my_jets_p4(*new_jets_p4);
          sort(my_jets_p4.begin(), my_jets_p4.end(), sortByPt);

          fillHistos(hptJet1, my_jets_p4[0].Pt(), weight, myType, nJetsIdx);
          fillHistos(hetaJet1, my_jets_p4[0].Eta(), weight, myType, nJetsIdx);

          if (new_njets > 1) {
            fillHistos(hptJet2, my_jets_p4[1].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJet2, my_jets_p4[1].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njets > 2) {
            fillHistos(hptJet3, my_jets_p4[2].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJet3, my_jets_p4[2].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njets > 3) {
            fillHistos(hptJet4, my_jets_p4[3].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJet4, my_jets_p4[3].Eta(), weight, myType, nJetsIdx);
          }
        }

        // Make a vector of jets sorted by pt and fill jet histograms
        if (new_njpts > 0) {
          vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > my_jpts_p4(*new_jpts_p4);
          sort(my_jpts_p4.begin(), my_jpts_p4.end(), sortByPt);

          fillHistos(hptJpt1,  my_jpts_p4[0].Pt(),  weight, myType, nJetsIdx);
          fillHistos(hetaJpt1, my_jpts_p4[0].Eta(), weight, myType, nJetsIdx);

          if (new_njpts > 1) {
            fillHistos(hptJpt2,  my_jpts_p4[1].Pt(),  weight, myType, nJetsIdx);
            fillHistos(hetaJpt2, my_jpts_p4[1].Eta(), weight, myType, nJetsIdx);
            fillHistos(hdrJ1J2, dRbetweenVectors( my_jpts_p4[0] , my_jpts_p4[1] ) , weight, myType, nJetsIdx);
          }
          if (new_njpts > 2) {
            fillHistos(hptJpt3,  my_jpts_p4[2].Pt(),  weight, myType, nJetsIdx);
            fillHistos(hetaJpt3, my_jpts_p4[2].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njpts > 3) {
            fillHistos(hptJpt4,  my_jpts_p4[3].Pt(),  weight, myType, nJetsIdx);
            fillHistos(hetaJpt4, my_jpts_p4[3].Eta(), weight, myType, nJetsIdx);
          }
        }
        
        VofP4 *new_btag_jpts_p4 =  &vjpts_btag_p4;
        int    nb               =   vjpts_btag_p4.size();

        // Make a vector of jets sorted by pt and fill jet histograms
        if ( nb > 0) {
         
          vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > my_btag_jpts_p4(*new_btag_jpts_p4);
          sort(my_btag_jpts_p4.begin(), my_btag_jpts_p4.end(), sortByPt);
          
          fillHistos(hptBtagJpt1,  my_btag_jpts_p4[0].Pt(),  weight, myType, nJetsIdx);

          if( nb > 1 ){
            fillHistos(hptBtagJpt2,  my_btag_jpts_p4[1].Pt(),  weight, myType, nJetsIdx);
          }

          if( nb > 2 ){
            fillHistos(hptBtagJpt3,  my_btag_jpts_p4[2].Pt(),  weight, myType, nJetsIdx);
          }

          if( nb > 3 ){
            fillHistos(hptBtagJpt4,  my_btag_jpts_p4[3].Pt(),  weight, myType, nJetsIdx);
          }

        }
      }
    } // entries

    delete f;
  } // currentFile
  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << endl;
  cout << "-------------------" << endl;
  cout << "| Dilepton yields |" << endl;
  cout << "-------------------" << endl;
  cout << "nee " << neetot << endl;
  cout << "nmm " << nmmtot << endl;
  cout << "nem " << nemtot << endl;
  cout << "tot " << neetot+nmmtot+nemtot << endl;

  cout << endl;
  cout << "-----------------------" << endl;
  cout << "| Preselection yields |" << endl;
  cout << "-----------------------" << endl;
  cout << "nee " << nee << endl;
  cout << "nmm " << nmm << endl;
  cout << "nem " << nem << endl;
  cout << "tot " << nee+nmm+nem << endl;
  cout << endl;

  if( TString(prefix).Contains("LM") ){
    cout << endl;
    cout << "N(gen)              " << ngen << endl;
    cout << endl;
    cout << "N(acc)  2010        " << nacc_2010        << endl;
    cout << "N(reco) 2010        " << nreco_2010       << endl;
    cout << "N(reco) 2010 noacc  " << nreco_noacc_2010 << endl;
    cout << Form("acceptance %.3f efficiency %.2f",(float)nacc_2010/(float)ngen,(float)nreco_2010/(float)nacc_2010) << endl;
    cout << endl;
    cout << "N(acc)  high MET        " << nacc_highmet        << endl;
    cout << "N(reco) high MET        " << nreco_highmet       << endl;
    cout << "N(reco) high MET noacc  " << nreco_noacc_highmet << endl;
    cout << Form("acceptance %.3f efficiency %.2f",(float)nacc_highmet/(float)ngen,(float)nreco_highmet/(float)nacc_highmet) << endl;
    cout << endl;
    cout << "N(acc)  high HT       " << nacc_highht        << endl;
    cout << "N(reco) high HT       " << nreco_highht       << endl;
    cout << "N(reco) high HT noacc " << nreco_noacc_highht << endl;
    cout << Form("acceptance %.3f efficiency %.2f",(float)nacc_highht/(float)ngen,(float)nreco_highht/(float)nacc_highht) << endl;
  }


  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  //delete d_llsol; //REPLACETOPMASS

  return 0;

}

//--------------------------------------------------------------------

bool ossusy_looper::passZSelection(int hypIdx)
{
  //require OS leptons, same flavor, in Z mass window
  if(hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                   return false;
  if(hyp_type()[hypIdx]==1 || hyp_type()[hypIdx]==2)                   return false;
  if(hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  return false;

  return true;
}

//--------------------------------------------------------------------

bool ossusy_looper::passTrigger(int dilType)
{
  //  bool hlt_ele15_lw_l1r = passHLTTrigger("HLT_Ele15_SW_L1R");
  bool hlt_ele15_lw_l1r = passHLTTrigger("HLT_Ele15_LW_L1R");
  bool hltMu9           = passHLTTrigger("HLT_Mu9");

  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;

  return true;
}

//--------------------------------------------------------------------
 
void ossusy_looper::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hyield = new TH1F(Form("%s_yield",prefix),Form("%s Event Yields",prefix),4,0,4);
  hyield->GetXaxis()->SetTitle("dil type");
  hyield->GetXaxis()->SetBinLabel(1,"all");
  hyield->GetXaxis()->SetBinLabel(2,"ee");
  hyield->GetXaxis()->SetBinLabel(3,"mm");
  hyield->GetXaxis()->SetBinLabel(4,"em");
  hyield->Sumw2();

  hyieldsig = new TH1F(Form("%s_yieldsig",prefix),Form("%s Event Yields (Signal Region)",prefix),4,0,4);
  hyieldsig->GetXaxis()->SetTitle("dil type");
  hyieldsig->GetXaxis()->SetBinLabel(1,"all");
  hyieldsig->GetXaxis()->SetBinLabel(2,"ee");
  hyieldsig->GetXaxis()->SetBinLabel(3,"mm");
  hyieldsig->GetXaxis()->SetBinLabel(4,"em");
  hyieldsig->Sumw2();

  hyield_unweighted = new TH1F(Form("%s_yield_unweighted",prefix),Form("%s Event Yields, Un-weighted",prefix),4,0,4);
  hyield_unweighted->GetXaxis()->SetTitle("dil type");
  hyield_unweighted->GetXaxis()->SetBinLabel(1,"all");
  hyield_unweighted->GetXaxis()->SetBinLabel(2,"ee");
  hyield_unweighted->GetXaxis()->SetBinLabel(3,"mm");
  hyield_unweighted->GetXaxis()->SetBinLabel(4,"em");
  hyield_unweighted->Sumw2();

  hyieldsig_unweighted = new TH1F(Form("%s_yieldsig_unweighted",prefix),Form("%s Event Yields (Signal Region), Un-weighted",prefix),4,0,4);
  hyieldsig_unweighted->GetXaxis()->SetTitle("dil type");
  hyieldsig_unweighted->GetXaxis()->SetBinLabel(1,"all");
  hyieldsig_unweighted->GetXaxis()->SetBinLabel(2,"ee");
  hyieldsig_unweighted->GetXaxis()->SetBinLabel(3,"mm");
  hyieldsig_unweighted->GetXaxis()->SetBinLabel(4,"em");
  hyieldsig_unweighted->Sumw2();

  hyield_weight = new TH1F(Form("%s_yield_weight",prefix),Form("%s Event Yields (No trigger weight)",prefix),4,0,4);
  hyield_weight->GetXaxis()->SetTitle("dil type");
  hyield_weight->GetXaxis()->SetBinLabel(1,"all");
  hyield_weight->GetXaxis()->SetBinLabel(2,"ee");
  hyield_weight->GetXaxis()->SetBinLabel(3,"mm");
  hyield_weight->GetXaxis()->SetBinLabel(4,"em");
  hyield_weight->Sumw2();

  char jetbins[5][7]    = {"0", "1", "2", "3", "#geq 4"};
  char suffixall[4][4]  = {"ee", "mm", "em", "all"};
  char njetCh[4][5]     = {"0j", "1j", "2j", "allj"};
  //char btagpoints[3][7] = {"loose", "medium", "tight"};

  //double binedges1000[11] = {0.,  20.,  40.,  60.,  80., 100., 120., 150.,  200.,  300.,  500.};
  //double binedges1500[11] = {0., 100., 200., 300., 400., 500., 600., 800., 1000., 1200., 1500.};
  double binedges1500[6] = {0., 100., 200., 400., 800., 1500.};
  //double binedges2000[11] = {0., 100., 200., 300., 400., 500., 600., 800., 1000., 1500., 2000.};

  if(strcmp("LMscan",prefix)==0){
    for(int im0 = 0 ; im0 < nm0points ; im0++){
      for(int im12 = 0 ; im12 < nm12points ; im12++){
          
        hsusydilPt[im0][im12] = new TH1F(Form("susy_hdilPt_m0_%i_m12_%i",im0,im12),
                                         Form("susy_hdilPt_m0_%i_m12_%i",im0,im12),60,0.,300.);
        hsusydilPt[im0][im12] -> GetXaxis()->SetTitle("Pt (GeV)");
          
        hsusytcmet[im0][im12] = new TH1F(Form("susy_htcmet_m0_%i_m12_%i",im0,im12),
                                         Form("susy_htcmet_m0_%i_m12_%i",im0,im12),60,0.,300.);
        hsusytcmet[im0][im12] -> GetXaxis()->SetTitle("tcmet (GeV)");

        hsusy_met_sumjetpt[im0][im12] = new TH2F(Form("susy_hmet_sumjetpt_m0_%i_m12_%i",im0,im12),
                                                 Form("susy_hmet_sumjetpt_m0_%i_m12_%i",im0,im12),200,0,2000,500,0,50);

        hsusy_met_sumjetpt[im0][im12] -> GetXaxis()->SetTitle("sumJetPt (GeV)");
        hsusy_met_sumjetpt[im0][im12] -> GetYaxis()->SetTitle("tcmet / #sqrt{tcsumet} (GeV^{1/2})");
      }
    }
  }


  for (int i = 0; i < 4; i++) {
    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffixall[i]),Form("%s_nJet_%s",prefix,suffixall[i]),10,0,10);	
    hnJet[i]->GetXaxis()->SetTitle("nJets");

    hnJpt[i] = new TH1F(Form("%s_hnJpt_%s",prefix,suffixall[i]),Form("%s_nJpt_%s",prefix,suffixall[i]),10,0,10);
    hnJpt[i]->GetXaxis()->SetTitle("nJpts");

    hnBtagJpt[i] = new TH1F(Form("%s_hnBtagJpt_%s",prefix,suffixall[i]),Form("%s_nBtagJpt_%s",prefix,suffixall[i]),10,0,10);
    hnBtagJpt[i]->GetXaxis()->SetTitle("nBtagJpts");

    hnHypJet[i] = new TH1F(Form("%s_hnHypJet_%s",prefix,suffixall[i]),Form("%s_nHypJet_%s",prefix,suffixall[i]),10,0,10);
    hnHypJet[i]->GetXaxis()->SetTitle("nHypJets");

    for(int k = 0; k < 5; k++) {
      hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnJet[i]->GetXaxis()->SetLabelSize(0.07);

      hnJpt[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnJpt[i]->GetXaxis()->SetLabelSize(0.07);

      hnHypJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnHypJet[i]->GetXaxis()->SetLabelSize(0.07);
    }

    for (int j = 0; j < 4; j++) {
      char suffix[7];
      sprintf(suffix, "%s_%s", njetCh[j], suffixall[i]);

      hdilMass_Z[i][j]   = new TH1F(Form("%s_hdilMass_Z_%s",prefix,suffix),
                                    Form("%s_dilMass_Z_%s" ,prefix,suffix),300,0,300);

      htcmet_event_Z[i][j]   = new TH1F(Form("%s_htcmet_event_Z_%s",prefix,suffix),
                                        Form("%s_tcmet_event_Z_%s" ,prefix,suffix),60,0,300);
      
      htcmet_looper_Z[i][j]   = new TH1F(Form("%s_htcmet_looper_Z_%s",prefix,suffix),
                                         Form("%s_tcmet_looper_Z_%s" ,prefix,suffix),60,0,300);
      
      hpfmet_Z[i][j]   = new TH1F(Form("%s_hpfmet_Z_%s",prefix,suffix),
                                  Form("%s_pfmet_Z_%s" ,prefix,suffix),60,0,300);
      
      hmucormet_Z[i][j]   = new TH1F(Form("%s_hmucormet_Z_%s",prefix,suffix),
                                     Form("%s_mucormet_Z_%s" ,prefix,suffix),60,0,300);
      
      hmucorjesmet_Z[i][j]   = new TH1F(Form("%s_hmucorjesmet_Z_%s",prefix,suffix),
                                        Form("%s_mucorjesmet_Z_%s" ,prefix,suffix),60,0,300);
      
      hetaz[i][j]   = new TH1F(Form("%s_hetaz_%s",prefix,suffix),
                               Form("%s_etaz_%s" ,prefix,suffix),1000,0,5);
            
      hmt[i][j]   = new TH1F(Form("%s_hmt_%s",prefix,suffix),
                             Form("%s_mt_%s" ,prefix,suffix),1000,0,500);

      habcd[i][j]   = new TH2F(Form("%s_habcd_%s",prefix,suffix),
                               Form("%s_abcd_%s" ,prefix,suffix),1500,0,1500,300,0,30);

      habcd[i][j]->Sumw2();

      habcd_nopresel[i][j]   = new TH2F(Form("%s_habcd_nopresel_%s",prefix,suffix),
                                        Form("%s_abcd_nopresel_%s" ,prefix,suffix),1500,0,1500,300,0,30);

      //delta(met-genmet)/genmet TH2, TProfile

      //event-level tcmet
      hdtcmetevent_genmet[i][j]   = new TH2F(Form("%s_hdtcmetevent_genmet_%s",prefix,suffix),
                                             Form("%s_hdtcmetevent_genmet_%s",prefix,suffix),50,0,500,50,-1,1); 

      tdtcmetevent_genmet[i][j]   = new TProfile(Form("%s_tdtcmetevent_genmet_%s",prefix,suffix),
                                                 Form("%s_tdtcmetevent_genmet_%s",prefix,suffix),50,0,500,-1,1); 

      hdtcmetevent_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      hdtcmetevent_genmet[i][j]->GetYaxis()->SetTitle("(tcmetevent-genmet)/genmet");
      tdtcmetevent_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      tdtcmetevent_genmet[i][j]->GetYaxis()->SetTitle("(tcmetevent-genmet)/genmet");

      //looper-level tcmet
      hdtcmetlooper_genmet[i][j]   = new TH2F(Form("%s_hdtcmetlooper_genmet_%s",prefix,suffix),
                                              Form("%s_hdtcmetlooper_genmet_%s",prefix,suffix),50,0,500,50,-1,1); 
      
      tdtcmetlooper_genmet[i][j]   = new TProfile(Form("%s_tdtcmetlooper_genmet_%s",prefix,suffix),
                                                  Form("%s_tdtcmetlooper_genmet_%s",prefix,suffix),50,0,500,-1,1); 
      
      hdtcmetlooper_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      hdtcmetlooper_genmet[i][j]->GetYaxis()->SetTitle("(tcmetlooper-genmet)/genmet");
      tdtcmetlooper_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      tdtcmetlooper_genmet[i][j]->GetYaxis()->SetTitle("(tcmetlooper-genmet)/genmet");

      //muon-corrected calomet
      hdmucormet_genmet[i][j]   = new TH2F(Form("%s_hdmucormet_genmet_%s",prefix,suffix),
                                           Form("%s_hdmucormet_genmet_%s",prefix,suffix),50,0,500,50,-1,1); 
      
      tdmucormet_genmet[i][j]   = new TProfile(Form("%s_tdmucormet_genmet_%s",prefix,suffix),
                                               Form("%s_tdmucormet_genmet_%s",prefix,suffix),50,0,500,-1,1); 
      
      hdmucormet_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      hdmucormet_genmet[i][j]->GetYaxis()->SetTitle("(mucormet-genmet)/genmet");
      tdmucormet_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      tdmucormet_genmet[i][j]->GetYaxis()->SetTitle("(mucormet-genmet)/genmet");

      //muon-corrected Type1 calomet
      hdmucorjesmet_genmet[i][j]   = new TH2F(Form("%s_hdmucorjesmet_genmet_%s",prefix,suffix),
                                              Form("%s_hdmucorjesmet_genmet_%s",prefix,suffix),50,0,500,50,-1,1); 
      
      tdmucorjesmet_genmet[i][j]   = new TProfile(Form("%s_tdmucorjesmet_genmet_%s",prefix,suffix),
                                                  Form("%s_tdmucorjesmet_genmet_%s",prefix,suffix),50,0,500,-1,1); 
      
      hdmucorjesmet_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      hdmucorjesmet_genmet[i][j]->GetYaxis()->SetTitle("(mucorjesmet-genmet)/genmet");
      tdmucorjesmet_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      tdmucorjesmet_genmet[i][j]->GetYaxis()->SetTitle("(mucorjesmet-genmet)/genmet");

      //pfmet
      hdpfmet_genmet[i][j]   = new TH2F(Form("%s_hdpfmet_genmet_%s",prefix,suffix),
                                        Form("%s_hdpfmet_genmet_%s",prefix,suffix),50,0,500,50,-1,1); 
      
      tdpfmet_genmet[i][j]   = new TProfile(Form("%s_tdpfmet_genmet_%s",prefix,suffix),
                                            Form("%s_tdpfmet_genmet_%s",prefix,suffix),50,0,500,-1,1); 
      
      hdpfmet_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      hdpfmet_genmet[i][j]->GetYaxis()->SetTitle("(pfmet-genmet)/genmet");
      tdpfmet_genmet[i][j]->GetXaxis()->SetTitle("genmet (GeV)");
      tdpfmet_genmet[i][j]->GetYaxis()->SetTitle("(pfmet-genmet)/genmet");
            

      Double_t xbins[66];

      for( unsigned int ibin = 0 ; ibin < 51 ; ++ibin ) xbins[ibin]    = 10  * ibin;
      for( unsigned int ibin = 0 ; ibin < 10 ; ++ibin ) xbins[ibin+51] = 50  * ibin + 550;
      for( unsigned int ibin = 0 ; ibin <  5 ; ++ibin ) xbins[ibin+61] = 100 * ibin + 1100;

      //for( int ibin = 0 ; ibin < 66 ; ++ibin) cout << xbins[ibin] << endl;
      //exit(0);

      //habcd_tprof[i][j]   = new TProfile(Form("%s_habcd_tprof_%s",prefix,suffix),
      //                                   Form("%s_abcd_tprof_%s" ,prefix,suffix),150,0,1500,0,30);
      habcd_tprof[i][j]   = new TProfile(Form("%s_habcd_tprof_%s",prefix,suffix),
                                         Form("%s_abcd_tprof_%s" ,prefix,suffix),65,xbins,0,30);

      habcd_tprof_nopresel[i][j]   = new TProfile(Form("%s_habcd_tprof_nopresel_%s",prefix,suffix),
                                                  Form("%s_abcd_tprof_nopresel_%s" ,prefix,suffix),150,0,1500,0,30);
            
      habcd[i][j]->GetXaxis()->SetTitle("sumJetPt (GeV)");
      habcd[i][j]->GetYaxis()->SetTitle("tcmet/#sqrt{sumJetPt} (GeV^{1/2})");
            
      habcd_tprof[i][j]->GetXaxis()->SetTitle("sumJetPt (GeV)");
      habcd_tprof[i][j]->GetYaxis()->SetTitle("tcmet/#sqrt{sumJetPt} (GeV^{1/2})");

      habcd_nopresel[i][j]->GetXaxis()->SetTitle("sumJetPt (GeV)");
      habcd_nopresel[i][j]->GetYaxis()->SetTitle("tcmet/#sqrt{sumJetPt} (GeV^{1/2})");
            
      habcd_tprof_nopresel[i][j]->GetXaxis()->SetTitle("sumJetPt (GeV)");
      habcd_tprof_nopresel[i][j]->GetYaxis()->SetTitle("tcmet/#sqrt{sumJetPt} (GeV^{1/2})");
      
      hmt2core[i][j] = new TH1F(Form("%s_hmt2core_%s",prefix,suffix),
                                Form("%s_mt2core_%s" ,prefix,suffix),100,0,200);
  
      hmt2jcore[i][j] = new TH1F(Form("%s_hmt2jcore_%s",prefix,suffix),
                                 Form("%s_mt2jcore_%s" ,prefix,suffix),500,0,500);
            
      hmt2j[i][j] = new TH1F(Form("%s_hmt2j_%s",prefix,suffix),
                             Form("%s_mt2j_%s" ,prefix,suffix),500,0,500);

      hmet_dilpt_all[i][j] = new TH2F(Form("%s_met_dilpt_all_%s",prefix,suffix),
                                      Form("%s_met_dilpt_all_%s" ,prefix,suffix),500,0,500,500,0,500);

      hmet_dilpt_signal[i][j] = new TH2F(Form("%s_met_dilpt_signal_%s",prefix,suffix),
                                         Form("%s_met_dilpt_signal_%s" ,prefix,suffix),500,0,500,500,0,500);

      hmet_dilpt_control[i][j] = new TH2F(Form("%s_met_dilpt_control_%s",prefix,suffix),
                                          Form("%s_met_dilpt_control_%s" ,prefix,suffix),500,0,500,500,0,500);
 
      hmt2j_all[i][j] = new TH1F(Form("%s_hmt2j_all_%s",prefix,suffix),
                                 Form("%s_mt2j_all_%s" ,prefix,suffix),1000,0,1000);

      hmt2j_signal[i][j] = new TH1F(Form("%s_hmt2j_signal_%s",prefix,suffix),
                                    Form("%s_mt2j_signal_%s" ,prefix,suffix),1000,0,1000);

      hmt2j_control[i][j] = new TH1F(Form("%s_hmt2j_control_%s",prefix,suffix),
                                     Form("%s_mt2j_control_%s" ,prefix,suffix),1000,0,1000);
           
      hdilMass_tcmet[i][j] = new TH2F(Form("%s_dilMass_tcmet_%s",prefix,suffix),
                                      Form("%s_dilMass_tcmet_%s",prefix,suffix),500,0,500,500,0,500);
            
      hsumJetPt_tcmet[i][j] = new TH2F(Form("%s_sumJetPt_tcmet_%s",prefix,suffix),
                                       Form("%s_sumJetPt_tcmet_%s",prefix,suffix),200,0,2000,500,0,500);
            
      hsumJetPt_tcmetsqrtsumet[i][j] = new TH2F(Form("%s_sumJetPt_tcmetsqrtsumet_%s",prefix,suffix),
                                                Form("%s_sumJetPt_tcmetsqrtsumet_%s",prefix,suffix),200,0,2000,500,0,50);

      hsumJetPt_tcmetsqrtsumet_prof[i][j] = new TProfile(Form("%s_sumJetPt_tcmetsqrtsumet_prof_%s",prefix,suffix),
                                                         Form("%s_sumJetPt_tcmetsqrtsumet_prof_%s",prefix,suffix),200,0,2000,0,50);
      
      hsumJetPt_tcmetsqrtsumet_prof[i][j]->GetXaxis()->SetTitle("sumJetPt (GeV)");
      hsumJetPt_tcmetsqrtsumet_prof[i][j]->GetYaxis()->SetTitle("tcmet / #sqrt{tcsumet} (GeV^{1/2})");

      htcsumet_tcmet_prof[i][j] = new TProfile(Form("%s_tcsumet_tcmet_prof_%s",prefix,suffix),
                                               Form("%s_tcsumet_tcmet_prof_%s",prefix,suffix),200,0,2000,0,500);
      
      htcsumet_tcmet_prof[i][j]->GetXaxis()->SetTitle("tcsumet (GeV)");
      htcsumet_tcmet_prof[i][j]->GetYaxis()->SetTitle("tcmet (GeV)");
      
      hsumJetPt_tcmetsumet[i][j] = new TH2F(Form("%s_sumJetPt_tcmetsumet_%s",prefix,suffix),
                                            Form("%s_sumJetPt_tcmetsumet_%s",prefix,suffix),200,0,2000,500,0,1);
            
      hetaZ_tcmet[i][j] = new TH2F(Form("%s_etaZ_tcmet_%s",prefix,suffix),
                                   Form("%s_etaZ_tcmet_%s",prefix,suffix),200,0,5,500,0,500);
            
      hetaZ_tcmetsqrtsumet[i][j] = new TH2F(Form("%s_etaZ_tcmetsqrtsumet_%s",prefix,suffix),
                                            Form("%s_etaZ_tcmetsqrtsumet_%s",prefix,suffix),200,0,5,500,0,50);
            
      hetaZ_tcmetsumet[i][j] = new TH2F(Form("%s_etaZ_tcmetsumet_%s",prefix,suffix),
                                        Form("%s_etaZ_tcmetsumet_%s",prefix,suffix),200,0,5,500,0,1);

      hgensumet_genmet_prof[i][j] = new TProfile(Form("%s_gensumet_genmet_prof_%s",prefix,suffix),
                                                 Form("%s_gensumet_genmet_prof_%s",prefix,suffix),200,0,2000,0,500);
            
      //hsumJetPt[i][j] = new TH1F(Form("%s_hsumJetPt_%s",prefix,suffix),Form("%s_sumJetPt_%s",prefix,suffix),5,binedges1500);
      hsumJetPt[i][j] = new TH1F(Form("%s_hsumJetPt_%s",prefix,suffix),Form("%s_sumJetPt_%s",prefix,suffix),100,0,500);
      hsumJetPt_cut[i][j] = new TH1F(Form("%s_hsumJetPt_cut_%s",prefix,suffix),Form("%s_sumJetPt_cut_%s",prefix,suffix),100,0,500);
      //hmeffJet[i][j] = new TH1F(Form("%s_hmeffJet_%s",prefix,suffix),Form("%s_meffJet_%s",prefix,suffix),5,binedges1500);
      hmeffJet[i][j] = new TH1F(Form("%s_hmeffJet_%s",prefix,suffix),Form("%s_meffJet_%s",prefix,suffix),100,0,1000);
      hDtcmetgenmetVsumJetPt[i][j] = new TH2F(Form("%s_hDtcmetgenmetVsumJetPt_%s",prefix,suffix),Form("%s_DtcmetgenmetVsumJetPt_%s",prefix,suffix),5,binedges1500,40,-50.,50.);
      hDmetmuonjesgenmetVsumJetPt[i][j] = new TH2F(Form("%s_hDmetmuonjesgenmetVsumJetPt_%s",prefix,suffix),Form("%s_DmetmuonjesgenmetVsumJetPt_%s",prefix,suffix),5,binedges1500,40,-50.,50.);

      hsumJptPt[i][j] = new TH1F(Form("%s_hsumJptPt_%s",prefix,suffix),Form("%s_sumJptPt_%s",prefix,suffix),100,0,1000);
      hsumJptPt_cut[i][j] = new TH1F(Form("%s_hsumJptPt_cut_%s",prefix,suffix),Form("%s_sumJptPt_cut_%s",prefix,suffix),100,0,1000);
      hmeffJPT[i][j] = new TH1F(Form("%s_hmeffJPT_%s",prefix,suffix),Form("%s_meffJPT_%s",prefix,suffix),100,0,1000);

      hsumHypPt[i][j] = new TH1F(Form("%s_hsumHypPt_%s",prefix,suffix),Form("%s_sumHypPt_%s",prefix,suffix),5,binedges1500);
      hmeffHyp[i][j] = new TH1F(Form("%s_hmeffHyp_%s",prefix,suffix),Form("%s_meffJHyp_%s",prefix,suffix),5,binedges1500);

      //
      // HISTOGRAMS INHERITED FROM SLAVA'S TTDIL LOOPER
      //

      helePt[i][j] = new TH1F(Form("%s_helePt_%s",prefix,suffix),Form("%s_elePt_%s",prefix,suffix),60,0.,300.);
      helePt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hmuPt[i][j]  = new TH1F(Form("%s_hmuPt_%s",prefix,suffix),Form("%s_muPt_%s",prefix,suffix),60,0.,300.);
      hmuPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hminLepPt[i][j]  = new TH1F(Form("%s_hminLepPt_%s",prefix,suffix),Form("%s_minLepPt_%s",prefix,suffix),60,0.,300.);
      hminLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hmaxLepPt[i][j]  = new TH1F(Form("%s_hmaxLepPt_%s",prefix,suffix),Form("%s_maxLepPt_%s",prefix,suffix),60,0.,300.);
      hmaxLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hminLepEta[i][j]  = new TH1F(Form("%s_hminLepEta_%s",prefix,suffix),Form("%s_minLepEta_%s",prefix,suffix),60,-3,3);
      hminLepEta[i][j]->GetXaxis()->SetTitle("Eta (GeV)");

      hmaxLepEta[i][j]  = new TH1F(Form("%s_hmaxLepEta_%s",prefix,suffix),Form("%s_maxLepEta_%s",prefix,suffix),60,-3,3);
      hmaxLepEta[i][j]->GetXaxis()->SetTitle("Eta (GeV)");

      helePhi[i][j] = new TH1F(Form("%s_helePhi_%s",prefix,suffix),Form("%s_elePhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      helePhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmuPhi[i][j]  = new TH1F(Form("%s_hmuPhi_%s",prefix,suffix),Form("%s_muPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hmuPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hdphiLep[i][j]  = new TH1F(Form("%s_hdphiLep_%s",prefix,suffix),Form("%s_dphiLep_%s",prefix,suffix),50,0.,TMath::Pi());
      hdphiLep[i][j]->GetXaxis()->SetTitle("#delta#phi_{ll}");

      hdrLep[i][j]  = new TH1F(Form("%s_hdrLep_%s",prefix,suffix),Form("%s_drLep_%s",prefix,suffix),50,0.,5);
      hdrLep[i][j]->GetXaxis()->SetTitle("#DeltaR(ll)");

      hdrJ1J2[i][j]  = new TH1F(Form("%s_hdrJ1J2_%s",prefix,suffix),Form("%s_drJ1J2_%s",prefix,suffix),50,0.,5);
      hdrJ1J2[i][j]->GetXaxis()->SetTitle("#DeltaR(J1J2)");

      heleEta[i][j] = new TH1F(Form("%s_heleEta_%s",prefix,suffix),Form("%s_eleEta_%s",prefix,suffix),60,-3.,3.);
      heleEta[i][j]->GetXaxis()->SetTitle("#eta");

      hmuEta[i][j]  = new TH1F(Form("%s_hmuEta_%s",prefix,suffix),Form("%s_muEta_%s",prefix,suffix),60,-3.,3.);
      hmuEta[i][j]->GetXaxis()->SetTitle("#eta");

      hdilMass[i][j] = new TH1F(Form("%s_hdilMass_%s",prefix,suffix),Form("%s_dilMass_%s",prefix,suffix),60,0.,300.);
      hdilMass[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");

      hdilMass_cut[i][j] = new TH1F(Form("%s_hdilMass_cut_%s",prefix,suffix),Form("%s_dilMass_cut_%s",prefix,suffix),60,0.,300.);
      hdilMass_cut[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");

      hdilPt[i][j] = new TH1F(Form("%s_hdilPt_%s",prefix,suffix),Form("%s_dilPt_%s",prefix,suffix),60,0.,300.);
      hdilPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      htopMass[i][j] = new TH1F(Form("%s_htopMass_%s",prefix,suffix),Form("%s_topMass_%s",prefix,suffix),500,0.,500.);
      htopMass[i][j]->GetXaxis()->SetTitle("Top Mass Estimate (GeV)");

      htopMassAllComb[i][j] = new TH1F(Form("%s_htopMassAllComb_%s",prefix,suffix),Form("%s_topMassAllComb_%s",prefix,suffix),1000,0.,1000.);
      htopMassAllComb[i][j]->GetXaxis()->SetTitle("Top Mass Estimate for all jet combinations (GeV)");

      //dilepton pT with z-veto applied
      hdilPt_zveto[i][j] = new TH1F(Form("%s_hdilPt_zveto_%s",prefix,suffix),Form("%s_dilPt_zveto_%s",prefix,suffix),60,0.,300.);
      hdilPt_zveto[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hdilPtSmeared[i][j] = new TH1F(Form("%s_hdilPtSmeared_%s",prefix,suffix),Form("%s_dilPtSmeared_%s",prefix,suffix),60,0.,300.);
      hdilPtSmeared[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      // changed binning from 2 GeV to 10 GeV
      hgenmet[i][j] = new TH1F(Form("%s_hgenmet_%s",prefix,suffix),Form("%s_genmet_%s",prefix,suffix),60,0.,300.);
      hgenmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hgenmetPhi[i][j] = new TH1F(Form("%s_hgenmetPhi_%s",prefix,suffix),Form("%s_genmetPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hgenmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hpfmet[i][j] = new TH1F(Form("%s_hpfmet_%s",prefix,suffix),Form("%s_pfmet_%s",prefix,suffix),60,0.,300.);
      hpfmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hpfmetPhi[i][j] = new TH1F(Form("%s_hpfmetPhi_%s",prefix,suffix),Form("%s_pfmetPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hpfmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetmuon[i][j] = new TH1F(Form("%s_hmetmuon_%s",prefix,suffix),Form("%s_metmuon_%s",prefix,suffix),60,0.,300.);
      hmetmuon[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hmetmuonPhi[i][j] = new TH1F(Form("%s_hmetmuonPhi_%s",prefix,suffix),Form("%s_metmuonPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hmetmuonPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetmuonVsDilepPt[i][j] = new TH2F(Form("%s_hmetmuonVsDilepPt_%s",prefix,suffix),Form("%s_metmuonVsDilepPt_%s",prefix,suffix),60,0.,300.,60,0.,300.);
      hmetmuonVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hmetmuonVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");

      hmetmuonOverPtVsDphi[i][j] = new TH2F(Form("%s_hmetmuonOverPtVsDphi_%s",prefix,suffix),Form("%s_metmuonOverPtVsDphi_%s",prefix,suffix),30,0.,3.,25,0.,TMath::Pi());
      hmetmuonVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hmetmuonVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");

      // met with muon corr and jes
      hmetmuonjes[i][j] = new TH1F(Form("%s_hmetmuonjes_%s",prefix,suffix),Form("%s_metmuonjes_%s",prefix,suffix),60,0.,300.);
      hmetmuonjes[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hmetmuonjesPhi[i][j] = new TH1F(Form("%s_hmetmuonjesPhi_%s",prefix,suffix),Form("%s_metmuonjesPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hmetmuonjesPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetmuonjesVsDilepPt[i][j] = new TH2F(Form("%s_hmetmuonjesVsDilepPt_%s",prefix,suffix),Form("%s_metmuonjesVsDilepPt_%s",prefix,suffix),60,0.,300.,60,0.,300.);
      hmetmuonjesVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hmetmuonjesVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");

      hmetmuonjesOverPtVsDphi[i][j] = new TH2F(Form("%s_hmetmuonjesOverPtVsDphi_%s",prefix,suffix),Form("%s_metmuonjesOverPtVsDphi_%s",prefix,suffix),30,0.,3.,25,0.,TMath::Pi());
      hmetmuonjesVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hmetmuonjesVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");

      // tc
      htcmet[i][j] = new TH1F(Form("%s_htcmet_%s",prefix,suffix),Form("%s_tcmet_%s",prefix,suffix),60,0.,300.);
      htcmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      htcmet_sqrtht[i][j] = new TH1F(Form("%s_htcmet_sqrtht_%s",prefix,suffix),Form("%s_tcmet_sqrtht_%s",prefix,suffix),200,0.,20.);
      htcmet_sqrtht[i][j]->GetXaxis()->SetTitle("tcmet/#sqrt{sumJetPt} (GeV^{1/2})");

      htcmetPhi[i][j] = new TH1F(Form("%s_htcmetPhi_%s",prefix,suffix),Form("%s_tcmetPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      htcmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      htcmetVsDilepPt[i][j] = new TH2F(Form("%s_htcmetVsDilepPt_%s",prefix,suffix),Form("%s_tcmetVsDilepPt_%s",prefix,suffix),60,0.,300.,60,0.,300.);
      htcmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      htcmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");

      htcmetOverPtVsDphi[i][j] = new TH2F(Form("%s_htcmetOverPtVsDphi_%s",prefix,suffix),Form("%s_tcmetOverPtVsDphi_%s",prefix,suffix),30,0.,3.,25,0.,TMath::Pi());
      htcmetOverPtVsDphi[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      htcmetOverPtVsDphi[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");

      hdphillvsmll[i][j] = new TH2F(Form("%s_dphillvsmll_%s",prefix,suffix),Form("%s_dphillvsmll_%s",prefix,suffix),100,10.,210.,50,0.,TMath::Pi());
      hdphillvsmll[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      hdphillvsmll[i][j]->GetYaxis()->SetTitle("#delta#phi_{ll}");

      hptJet1[i][j] = new TH1F(Form("%s_hptJet1_%s",prefix,suffix),Form("%s_ptJet1_%s",prefix,suffix),60,0.,300.);
      hptJet1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet2[i][j] = new TH1F(Form("%s_hptJet2_%s",prefix,suffix),Form("%s_ptJet2_%s",prefix,suffix),60,0.,300.);
      hptJet2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet3[i][j] = new TH1F(Form("%s_hptJet3_%s",prefix,suffix),Form("%s_ptJet3_%s",prefix,suffix),60,0.,300.);
      hptJet3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet4[i][j] = new TH1F(Form("%s_hptJet4_%s",prefix,suffix),Form("%s_ptJet4_%s",prefix,suffix),60,0.,300.);
      hptJet4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hetaJet1[i][j] = new TH1F(Form("%s_hetaJet1_%s",prefix,suffix),Form("%s_etaJet1_%s",prefix,suffix),50,-4.,4.);
      hetaJet1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet2[i][j] = new TH1F(Form("%s_hetaJet2_%s",prefix,suffix),Form("%s_etaJet2_%s",prefix,suffix),50,-4.,4.);
      hetaJet2[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet3[i][j] = new TH1F(Form("%s_hetaJet3_%s",prefix,suffix),Form("%s_etaJet3_%s",prefix,suffix),50,-4.,4.);
      hetaJet3[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet4[i][j] = new TH1F(Form("%s_hetaJet4_%s",prefix,suffix),Form("%s_etaJet4_%s",prefix,suffix),50,-4.,4.);
      hetaJet4[i][j]->GetXaxis()->SetTitle("#eta");

      hptJpt1[i][j] = new TH1F(Form("%s_hptJpt1_%s",prefix,suffix),Form("%s_ptJpt1_%s",prefix,suffix),60,0.,300.);
      hptJpt1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJpt2[i][j] = new TH1F(Form("%s_hptJpt2_%s",prefix,suffix),Form("%s_ptJpt2_%s",prefix,suffix),60,0.,300.);
      hptJpt2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJpt3[i][j] = new TH1F(Form("%s_hptJpt3_%s",prefix,suffix),Form("%s_ptJpt3_%s",prefix,suffix),60,0.,300.);
      hptJpt3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJpt4[i][j] = new TH1F(Form("%s_hptJpt4_%s",prefix,suffix),Form("%s_ptJpt4_%s",prefix,suffix),60,0.,300.);
      hptJpt4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptBtagJpt1[i][j] = new TH1F(Form("%s_hptBtagJpt1_%s",prefix,suffix),Form("%s_ptBtagJpt1_%s",prefix,suffix),60,0.,300.);
      hptBtagJpt1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptBtagJpt2[i][j] = new TH1F(Form("%s_hptBtagJpt2_%s",prefix,suffix),Form("%s_ptBtagJpt2_%s",prefix,suffix),60,0.,300.);
      hptBtagJpt2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptBtagJpt3[i][j] = new TH1F(Form("%s_hptBtagJpt3_%s",prefix,suffix),Form("%s_ptBtagJpt3_%s",prefix,suffix),60,0.,300.);
      hptBtagJpt3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptBtagJpt4[i][j] = new TH1F(Form("%s_hptBtagJpt4_%s",prefix,suffix),Form("%s_ptBtagJpt4_%s",prefix,suffix),60,0.,300.);
      hptBtagJpt4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hetaJpt1[i][j] = new TH1F(Form("%s_hetaJpt1_%s",prefix,suffix),Form("%s_etaJpt1_%s",prefix,suffix),50,-4.,4.);
      hetaJpt1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJpt2[i][j] = new TH1F(Form("%s_hetaJpt2_%s",prefix,suffix),Form("%s_etaJpt2_%s",prefix,suffix),50,-4.,4.);
      hetaJpt2[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJpt3[i][j] = new TH1F(Form("%s_hetaJpt3_%s",prefix,suffix),Form("%s_etaJpt3_%s",prefix,suffix),50,-4.,4.);
      hetaJpt3[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJpt4[i][j] = new TH1F(Form("%s_hetaJpt4_%s",prefix,suffix),Form("%s_etaJpt4_%s",prefix,suffix),50,-4.,4.);
      hetaJpt4[i][j]->GetXaxis()->SetTitle("#eta");

      hptHypJet1[i][j] = new TH1F(Form("%s_hptHypJet1_%s",prefix,suffix),Form("%s_ptHypJet1_%s",prefix,suffix),60,0.,300.);
      hptHypJet1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptHypJet2[i][j] = new TH1F(Form("%s_hptHypJet2_%s",prefix,suffix),Form("%s_ptHypJet2_%s",prefix,suffix),60,0.,300.);
      hptHypJet2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptHypJet3[i][j] = new TH1F(Form("%s_hptHypJet3_%s",prefix,suffix),Form("%s_ptHypJet3_%s",prefix,suffix),60,0.,300.);
      hptHypJet3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptHypJet4[i][j] = new TH1F(Form("%s_hptHypJet4_%s",prefix,suffix),Form("%s_ptHypJet4_%s",prefix,suffix),60,0.,300.);
      hptHypJet4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hetaHypJet1[i][j] = new TH1F(Form("%s_hetaHypJet1_%s",prefix,suffix),Form("%s_etaHypJet1_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaHypJet2[i][j] = new TH1F(Form("%s_hetaHypJet2_%s",prefix,suffix),Form("%s_etaHypJet2_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet2[i][j]->GetXaxis()->SetTitle("#eta");

      hetaHypJet3[i][j] = new TH1F(Form("%s_hetaHypJet3_%s",prefix,suffix),Form("%s_etaHypJet3_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet3[i][j]->GetXaxis()->SetTitle("#eta");

      hetaHypJet4[i][j] = new TH1F(Form("%s_hetaHypJet4_%s",prefix,suffix),Form("%s_etaHypJet4_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet4[i][j]->GetXaxis()->SetTitle("#eta");

      /*************************/

      hsumJetPt[i][j]->Sumw2(); 
      hsumJetPt_cut[i][j]->Sumw2(); 
      hmeffJet[i][j]->Sumw2(); 

      hsumJptPt[i][j]->Sumw2(); 
      hsumJptPt_cut[i][j]->Sumw2(); 
      hmeffJPT[i][j]->Sumw2(); 

      hsumHypPt[i][j]->Sumw2(); 
      hmeffHyp[i][j]->Sumw2(); 

      if (j == 0)  {
        hnJet[i]->Sumw2();
        hnJpt[i]->Sumw2();
        hnHypJet[i]->Sumw2();
        hnBtagJpt[i]->Sumw2();
      }
      helePt[i][j]->Sumw2();
      hdrLep[i][j]->Sumw2();
      hdrJ1J2[i][j]->Sumw2();
      hmuPt[i][j]->Sumw2();
      hminLepPt[i][j]->Sumw2();
      hmaxLepPt[i][j]->Sumw2();
      hminLepEta[i][j]->Sumw2();
      hmaxLepEta[i][j]->Sumw2();
      helePhi[i][j]->Sumw2();
      hmuPhi[i][j]->Sumw2();
      hdphiLep[i][j]->Sumw2();
      heleEta[i][j]->Sumw2();
      hmuEta[i][j]->Sumw2();
      hdilMass[i][j]->Sumw2();
      hdilMass_cut[i][j]->Sumw2();
      hdilPt[i][j]->Sumw2();
      htopMass[i][j]->Sumw2();
      htopMassAllComb[i][j]->Sumw2();
      hdilPt_zveto[i][j]->Sumw2();
      hdilPtSmeared[i][j]->Sumw2();
      hgenmet[i][j]->Sumw2();
      hgenmetPhi[i][j]->Sumw2();
      hmetmuon[i][j]->Sumw2();
      hmetmuonPhi[i][j]->Sumw2();
      hmetmuonjes[i][j]->Sumw2();
      hmetmuonjesPhi[i][j]->Sumw2();
      htcmet[i][j]->Sumw2();
      htcmetPhi[i][j]->Sumw2();
      hpfmet[i][j]->Sumw2();
      hpfmetPhi[i][j]->Sumw2();
      hptJet1[i][j]->Sumw2();
      hptJet2[i][j]->Sumw2();
      hptJet3[i][j]->Sumw2();
      hptJet4[i][j]->Sumw2();
      hetaJet1[i][j]->Sumw2();
      hetaJet2[i][j]->Sumw2();
      hetaJet3[i][j]->Sumw2();
      hetaJet4[i][j]->Sumw2();
      hptJpt1[i][j]->Sumw2();
      hptJpt2[i][j]->Sumw2();
      hptJpt3[i][j]->Sumw2();
      hptJpt4[i][j]->Sumw2();
      hptBtagJpt1[i][j]->Sumw2();
      hptBtagJpt2[i][j]->Sumw2();
      hptBtagJpt3[i][j]->Sumw2();
      hptBtagJpt4[i][j]->Sumw2();
      hetaJpt1[i][j]->Sumw2();
      hetaJpt2[i][j]->Sumw2();
      hetaJpt3[i][j]->Sumw2();
      hetaJpt4[i][j]->Sumw2();
      hptHypJet1[i][j]->Sumw2();
      hptHypJet2[i][j]->Sumw2();
      hptHypJet3[i][j]->Sumw2();
      hptHypJet4[i][j]->Sumw2();
      hetaHypJet1[i][j]->Sumw2();
      hetaHypJet2[i][j]->Sumw2();
      hetaHypJet3[i][j]->Sumw2();
      hetaHypJet4[i][j]->Sumw2();
    } // njet loop
  } // channel loop

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

//function to add additional smearing to met. metscale = 1.5 --> additional 50% smearing
float ossusy_looper::smearMet( float met , float sumjetpt , float metscale ){

  if( metscale < 1. ){
    cout << "metscale must be >1!, quitting..." << endl;
    exit(0);
  }
  
  float res       = returnSigma( sumjetpt , e_tcmet );
  float addmetres = sqrt( metscale * metscale - 1 ) * res;
  float smear     = random3_->Gaus( 0 , addmetres );
  
  return met + smear;

}

//--------------------------------------------------------------------

float returnSigma (float sumJetPt, ossusy_looper::MetTypeEnum metType)
{
  const float xbins[]         = { 100   , 200   , 400   , 800           };
  const float tcmet_sigma[]   = { 13.02 , 13.44 , 14.64 , 17.19 , 21.95 };
  const float muonjes_sigma[] = { 21.13 , 20.3  , 21.08 , 24.06 , 29.08 };

  if (metType == ossusy_looper::e_tcmet) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return tcmet_sigma[j];
    }

    return tcmet_sigma[4];
  }
  else if (metType == ossusy_looper::e_muonjes) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return muonjes_sigma[j];
    }

    return muonjes_sigma[4];
  }
}

//--------------------------------------------------------------------

float returnBias(float sumJetPt, ossusy_looper::MetTypeEnum metType)
{
  const float xbins[]         = { 100  , 200  , 400  , 800         };
  const float tcmet_bias[]    = { 0.87 , 0.88 , 0.91 , 0.95 , 0.99 };
  const float muonjes_bias[]  = { 0.94 , 1.   , 1.03 , 1.03 , 1.07 };

  if (metType == ossusy_looper::e_tcmet) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return tcmet_bias[j];
    }

    return tcmet_bias[4];
  }
  else if (metType == ossusy_looper::e_muonjes) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return muonjes_bias[j];
    }

    return muonjes_bias[4];
  }
}

//--------------------------------------------------------------------
                                                                     
                                             
float ossusy_looper::getCosThetaStarWeight(){


  int nels = 0;
  int nmus  = 0;
  int ntaus = 0;
  int nleps = 0;
    
  nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
    
  ///////////////////////////
  // Generator Information //
  ///////////////////////////
  //
  int ntplus  = 0;
  int ntminus = 0;
  int nWplus  = 0;
  int nWminus = 0;
  int nbplus  = 0;
  int nbminus = 0;
  int nlplus  = 0;
  int nlminus = 0;
  int nnu     = 0;
  int nnubar  = 0;

  LorentzVector genp4_tplus_;
  LorentzVector genp4_tminus_;
  LorentzVector genp4_Wplus_;
  LorentzVector genp4_Wminus_;
  LorentzVector genp4_lplus_;
  LorentzVector genp4_lminus_;
  LorentzVector genp4_bplus_;
  LorentzVector genp4_bminus_;
  LorentzVector genp4_nu_;
  LorentzVector genp4_nubar_;
  LorentzVector genp4_Wplus_bminus_;
  LorentzVector genp4_Wminus_bplus_;
  LorentzVector genp4_lminus_nubar_;
  LorentzVector genp4_lplus_nu_;
  LorentzVector genp4_lminus_lplus_;
  LorentzVector genp4_nu_nubar_;
  LorentzVector genp4_Wminus_tCM_;
  LorentzVector genp4_Wplus_tCM_;
  LorentzVector genp4_lminus_tCM_;
  LorentzVector genp4_nubar_tCM_;
  LorentzVector genp4_lplus_tCM_;
  LorentzVector genp4_nu_tCM_;
  LorentzVector genp4_lminus_WCM_;
  LorentzVector genp4_nubar_WCM_;
  LorentzVector genp4_lplus_WCM_;
  LorentzVector genp4_nu_WCM_;
  LorentzVector genp4_lminus_nubar_WCM_;
  LorentzVector genp4_lplus_nu_WCM_;
  
  float theta_lminus_nubar_WCM_;
  float theta_lplus_nu_WCM_;
  int   genid_lplus_;
  int   genid_lminus_;

  // loop on generator particles
  for( unsigned int i=0; i < genps_id().size(); i++ ){

    // top
    if( genps_id().at(i) == 6 ){
      genp4_tplus_ = genps_p4().at(i);
      ntplus++;
    }
    if( genps_id().at(i) == -6 ){
      genp4_tminus_ = genps_p4().at(i);
      ntminus++;
    }

    // W's from top
    if( genps_id().at(i) == 24 && genps_id_mother().at(i) == 6 ){
      genp4_Wplus_ = genps_p4().at(i);
      nWplus++;
    }
    if( genps_id().at(i) == -24 && genps_id_mother().at(i) == -6 ){
      genp4_Wminus_ = genps_p4().at(i);
      nWminus++;
    }

    // leptons from W+
    if( genps_id_mother().at(i) == 24 ){
      if( genps_id().at(i) == -11 || genps_id().at(i) == -13 ){
        genp4_lplus_ = genps_p4().at(i);
        nlplus++;
        if( genps_id().at(i) == -11 ){
          genid_lplus_ = -11;
        } 
        else if( genps_id().at(i) == -13 ){
          genid_lplus_ = -13;
        }
      }
      if( genps_id().at(i) == 12 || genps_id().at(i) == 14 ){
        genp4_nu_ = genps_p4().at(i);
        nnu++;
      }
    }

    // leptons from W-
    if( genps_id_mother().at(i) == -24 ){
      if( genps_id().at(i) == 11 || genps_id().at(i) == 13 ){
        genp4_lminus_ = genps_p4().at(i);
        nlminus++;
        if( genps_id().at(i) == 11 ){
          genid_lminus_ = 11;;
        } 
        else if( genps_id().at(i) == 13 ){
          genid_lminus_ = 13;
        }
      }
      if( genps_id().at(i) == -12 || genps_id().at(i) == -14 ){
        genp4_nubar_ = genps_p4().at(i);
        nnubar++;
      }
    }

    // b's
    if( genps_id().at(i) == 5 && genps_id_mother().at(i) == 6){
      nbminus++;
      genp4_bminus_ = genps_p4().at(i);
    }
    if( genps_id().at(i) == -5 && genps_id_mother().at(i) == -6){
      nbplus++;
      genp4_bplus_ = genps_p4().at(i);
    }

  } // end loop on generator particles

  // Construct some composite 4-vectors
  genp4_Wplus_bminus_ = genp4_Wplus_  + genp4_bminus_;
  genp4_Wminus_bplus_ = genp4_Wminus_ + genp4_bplus_;
  genp4_lminus_nubar_ = genp4_lminus_ + genp4_nubar_;
  genp4_lplus_nu_     = genp4_lplus_  + genp4_nu_;
  genp4_lminus_lplus_ = genp4_lminus_ + genp4_lplus_;
  genp4_nu_nubar_     = genp4_nu_     + genp4_nubar_;

  // Sanity
  if( ntplus != 1 ){ 
    cout << "ERROR: did not find exactly 1 t" << endl;
    exit(1);
  }
  if( ntminus != 1 ){ 
    cout << "ERROR: did not find exactly 1 tbar" << endl;
    exit(1);
  }
  if( nWplus != 1 ){ 
    cout << "ERROR: did not find exactly 1 W+" << endl;
    exit(1);
  }
  if( nWminus != 1 ){ 
    cout << "ERROR: did not find exactly 1 W-" << endl;
    exit(1);
  }
  if( nbplus != 1 ){ 
    cout << "ERROR: did not find exactly 1 b+ from top" << endl;
    //exit(1);
  }
  if( nbminus != 1 ){ 
    cout << "ERROR: did not find exactly 1 b- from top" << endl;
    //exit(1);
  }
  // 
  if( ntaus == 0 ){
    if( nlplus != 1 ){ 
      cout << "ERROR: did not find exactly 1 e+ or mu+" << endl;
      exit(1);
    }
    if( nlminus != 1 ){ 
      cout << "ERROR: did not find exactly 1 e- or mu-" << endl;
      exit(1);
    }
    if( nnu != 1 ){
      cout << "ERROR: did not find exactly 1 nu" << endl;
      exit(1);
    }
    if( nnubar != 1 ){
      cout << "ERROR: did not find exactly 1 nubar" << endl;
      exit(1);
    }
  } 
  else {
    //cout << "Unhandled Case: N Taus = " << ntaus << endl << endl;
    return -2;
  }



      

  // Boost from the LAB to the top CM
  ROOT::Math::Boost boost_tplus_CM( genp4_tplus_.BoostToCM().x(), genp4_tplus_.BoostToCM().y(), genp4_tplus_.BoostToCM().z() );
  ROOT::Math::Boost boost_tminus_CM( genp4_tminus_.BoostToCM().x(), genp4_tminus_.BoostToCM().y(), genp4_tminus_.BoostToCM().z() );
  
  // Boost from the LAB to the W CM
  ROOT::Math::Boost boost_Wplus_CM( genp4_Wplus_.BoostToCM().x(), genp4_Wplus_.BoostToCM().y(), genp4_Wplus_.BoostToCM().z() );
  ROOT::Math::Boost boost_Wminus_CM( genp4_Wminus_.BoostToCM().x(), genp4_Wminus_.BoostToCM().y(), genp4_Wminus_.BoostToCM().z() );
  
  // Get the W, lepton, and neutrino 4-vectors in the top rest frame
  genp4_Wminus_tCM_       = boost_tminus_CM * genp4_Wminus_;
  genp4_Wplus_tCM_        = boost_tplus_CM  * genp4_Wplus_;
  genp4_lminus_tCM_       = boost_tminus_CM * genp4_lminus_;
  genp4_nubar_tCM_        = boost_tminus_CM * genp4_nubar_;
  genp4_lplus_tCM_        = boost_tplus_CM  * genp4_lplus_;
  genp4_nu_tCM_           = boost_tplus_CM  * genp4_nu_;

  // Boost from the top CM to the W CM ( tW notation )
  ROOT::Math::Boost boost_Wplus_tWCM( genp4_Wplus_tCM_.BoostToCM().x(), genp4_Wplus_tCM_.BoostToCM().y(), genp4_Wplus_tCM_.BoostToCM().z() );
  ROOT::Math::Boost boost_Wminus_tWCM( genp4_Wminus_tCM_.BoostToCM().x(), genp4_Wminus_tCM_.BoostToCM().y(), genp4_Wminus_tCM_.BoostToCM().z() );

  // Boost from the W CM to the top CM ( Wt notation )
  ROOT::Math::Boost boost_Wplus_WtCM  = boost_Wplus_tWCM.Inverse();
  ROOT::Math::Boost boost_Wminus_WtCM = boost_Wminus_tWCM.Inverse();
        
  // Get the lepton and neutrino 4-vectors in the W rest frame

  //--- LAB -> top CM -> W CM ---//
  genp4_lminus_WCM_       = boost_Wminus_tWCM * genp4_lminus_tCM_;
  genp4_nubar_WCM_        = boost_Wminus_tWCM * genp4_nubar_tCM_;
  genp4_lplus_WCM_        = boost_Wplus_tWCM  * genp4_lplus_tCM_;
  genp4_nu_WCM_           = boost_Wplus_tWCM  * genp4_nu_tCM_;

  /*
  //--- LAB -> W -> top ---//
  genp4_lminus_WCM_       = boost_Wminus_CM * genp4_lminus_;
  genp4_nubar_WCM_        = boost_Wminus_CM * genp4_nubar_;
  genp4_lplus_WCM_        = boost_Wplus_CM  * genp4_lplus_;
  genp4_nu_WCM_           = boost_Wplus_CM  * genp4_nu_;

  genp4_lminus_WCM_       = boost_Wminus_WtCM * genp4_lminus_WCM_;
  genp4_nubar_WCM_        = boost_Wminus_WtCM * genp4_nubar_WCM_;
  genp4_lplus_WCM_        = boost_Wplus_WtCM  * genp4_lplus_WCM_;
  genp4_nu_WCM_           = boost_Wplus_WtCM  * genp4_nu_WCM_;
  */

  // Composite 4-vectors
  genp4_lminus_nubar_WCM_ = genp4_lminus_WCM_ + genp4_nubar_WCM_;
  genp4_lplus_nu_WCM_     = genp4_lplus_WCM_  + genp4_nu_WCM_;

  // W +/- in the lab frame here
  // want W direction in the top rest frame
  // NB: this says theta, but it is real cosTheta

  theta_lminus_nubar_WCM_ = genp4_lminus_WCM_.Px()*genp4_Wminus_tCM_.Px() +
    genp4_lminus_WCM_.Py()*genp4_Wminus_tCM_.Py() +
    genp4_lminus_WCM_.Pz()*genp4_Wminus_tCM_.Pz();

  theta_lplus_nu_WCM_     = genp4_lplus_WCM_.Px()*genp4_Wplus_tCM_.Px() +
    genp4_lplus_WCM_.Py()*genp4_Wplus_tCM_.Py() +
    genp4_lplus_WCM_.Pz()*genp4_Wplus_tCM_.Pz();

  theta_lminus_nubar_WCM_ /= genp4_lminus_WCM_.P()*genp4_Wminus_tCM_.P();
  theta_lplus_nu_WCM_     /= genp4_lplus_WCM_.P()*genp4_Wplus_tCM_.P();

  float x1 = theta_lminus_nubar_WCM_;
  float x2 = theta_lplus_nu_WCM_;

  float f1 = 0.703 * (1 - x1 * x1 ) + 0.5 * 0.297 * ( 1 - x1 ) * ( 1 - x1 );
  float f2 = 0.703 * (1 - x2 * x2 ) + 0.5 * 0.297 * ( 1 - x2 ) * ( 1 - x2 );

  
  float weight = -1;

  if( f1 > 0 && f2 > 0 ) weight = 1. / ( f1 * f2 );

  return weight;

}



//*****************************************************************
// get the FR weight
//*****************************************************************

double ossusy_looper::getFRWeight(const int hypIdx, SimpleFakeRate* mufr, SimpleFakeRate * elfr, FREnum frmode, bool isData) {

  //std::cout<<"Called ossusy_looper::getFRWeight"<<std::endl;

  bool  estimateQCD   = false;
  bool  estimateWJets = false;

  if ( frmode == e_qcd ) {
    estimateQCD   = true;
    estimateWJets = false;
  } 
  else if( frmode == e_wjets ) {
    estimateQCD   = false;
    estimateWJets = true;
  }
  else {
    std::cout<<"ossusy_looper::getFRWeight: bad FR mode given, fix this!"<<std::endl;
    return -9999.;
  }

  if(hyp_type()[hypIdx] == 0) {

    bool isGoodMut = false;
    bool isGoodMul = false;
    bool isFOMut   = false;
    bool isFOMul   = false;
    
    unsigned int iMut = hyp_lt_index()[hypIdx];
    unsigned int iMul = hyp_ll_index()[hypIdx];
    
    if( muonId( iMut , OSGeneric_v3 ) ) {
      isGoodMut = true;
    }
    if( muonId( iMul , OSGeneric_v3 ) ) {
      isGoodMul = true;
    }
    if( muonId( iMut , OSGeneric_v3_FO ) ) {
      isFOMut = true;
    }
    if( muonId( iMul , OSGeneric_v3_FO ) ) {
      isFOMul = true;
    }

    //for both WJets and QCD, we need both to be FOs at least
    if(!isFOMut || !isFOMul)
      return -9999.;

    //if we want to estimate the fakes for QCD, then we ask that 
    //both are not num objects, and that both are FO
    if(estimateQCD) {
      
      //if at least one is a Numerator lepton, we return
      if( isGoodMut || isGoodMul) 
        return -9999.;
      
      double FRMut = mufr->getFR(mus_p4()[iMut].pt(), mus_p4()[iMut].eta());
      double FRMul = mufr->getFR(mus_p4()[iMul].pt(), mus_p4()[iMul].eta());
      return (FRMut/(1-FRMut))*(FRMul/(1-FRMul));
    } 
    else if(estimateWJets) {
      
      //need one to be a Numerator lepton, and the other to be FO but not num
      if( isGoodMut && !isGoodMul && isFOMul) {
        double FR = mufr->getFR(mus_p4()[iMul].pt(), mus_p4()[iMul].eta());
        //cout << "mm, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      
      //check the other muon
      if( isGoodMul && !isGoodMut && isFOMut) {
        double FR = mufr->getFR(mus_p4()[iMut].pt(), mus_p4()[iMut].eta());
        //cout << "mm, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
    }//estimate WJets
    return -9999.;
  }//mumu case
  

  //now we do the ee case
  if(hyp_type()[hypIdx] == 3) {
	  
    unsigned int iElt = hyp_lt_index()[hypIdx];
    unsigned int iEll = hyp_ll_index()[hypIdx];
	  
    bool isGoodElt = false;
    bool isGoodEll = false;
    bool isFOElt   = false;
    bool isFOEll   = false;

    if( pass_electronSelection( iElt , electronSelection_el_OSV3 ) ) {
      isGoodElt = true;
    }
    if( pass_electronSelection( iEll , electronSelection_el_OSV3 ) ) {
      isGoodEll = true;
    }
    if( pass_electronSelection( iElt , electronSelection_el_OSV3_FO ) ) {
      isFOElt   = true;
    }
    if( pass_electronSelection( iEll , electronSelection_el_OSV3_FO ) ) {
      isFOEll   = true;
    }
    
    //for both WJets and QCD, we need both to be FOs at least
    //if both are good, we continue
    if( !isFOElt || !isFOEll)
      return -9999.;

    if(estimateQCD) {
      
      //if at least one is a Numerator object, then we return -9999.
      if( isGoodElt || isGoodEll) 
        return -9999.;
      
      double FRElt = elfr->getFR(els_p4()[iElt].pt(), els_p4()[iElt].eta());
      double FREll = elfr->getFR(els_p4()[iEll].pt(), els_p4()[iEll].eta());
      //cout << "ee, FRlt, FRll, FR " << FRElt << " " << FREll << " " << (FRElt/(1-FRElt))*(FREll/(1-FREll)) << endl;
      return (FRElt/(1-FRElt))*(FREll/(1-FREll));
    } 
    else if(estimateWJets) {
      
      if(isGoodElt && !isGoodEll && isFOEll) {
        double FR = elfr->getFR(els_p4()[iEll].pt(), els_p4()[iEll].eta());
        //cout << "ee, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      //check the other electron 
      if(isGoodEll && !isGoodElt && isFOElt) {
        double FR = elfr->getFR(els_p4()[iElt].pt(), els_p4()[iElt].eta());
        //cout << "ee, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      return -9999.;
    }//estimateWJets
    
  }//ee case

  if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {

    int iEl = 0;
    int iMu = 0;
    
    if     ( abs(hyp_ll_id()[hypIdx])==11 && abs(hyp_lt_id()[hypIdx])==13 ){
      iEl = hyp_ll_index()[hypIdx];
      iMu = hyp_lt_index()[hypIdx];
    }
    else if( abs(hyp_ll_id()[hypIdx])==13 && abs(hyp_lt_id()[hypIdx])==11 ){
      iEl = hyp_lt_index()[hypIdx];
      iMu = hyp_ll_index()[hypIdx];
    }
    else{
      cout << "ID ll " << hyp_ll_id()[hypIdx] << endl;
      cout << "ID lt " << hyp_lt_id()[hypIdx] << endl;
      cout << "Error in getFRWeight, quitting!" << endl;
      exit(0); 
    }

    bool isGoodEl = false;
    bool isFOEl   = false;
    bool isGoodMu = false;
    bool isFOMu   = false;

    if( pass_electronSelection( iEl , electronSelection_el_OSV3 ) ){
      isGoodEl = true;
    }
    if( muonId( iMu , OSGeneric_v3 ) ) { 
      isGoodMu = true;
    }
    if( pass_electronSelection( iEl , electronSelection_el_OSV3_FO ) ){
      isFOEl = true;
    }
    if( muonId( iMu , OSGeneric_v3_FO ) ) { 
      isFOMu = true;
    }
    
    //if either fail FO, return!!!
    if(!isFOMu || !isFOEl)
      return -9999.;
    
    if(estimateQCD ) {
      
      //if at least one is a numerator, then we fail
      if(isGoodMu || isGoodEl)
        return -9999.;
      
      double FRMu = mufr->getFR(mus_p4()[iMu].pt(), mus_p4()[iMu].eta());
      double FREl = elfr->getFR(els_p4()[iEl].pt(), els_p4()[iEl].eta());
      return FRMu*FREl/(1-FRMu)/(1-FREl);
    } 
    else if(estimateWJets) {
      
      //need one to be a numerator lepton and the other to be a FO
      if(isGoodMu && !isGoodEl && isFOEl) {
        double FR = elfr->getFR(els_p4()[iEl].pt(), els_p4()[iEl].eta());
        //cout << "emu, el FR, FR/(1-FR): " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      
      if(isGoodEl && !isGoodMu && isFOMu) {
        double FR = mufr->getFR(mus_p4()[iMu].pt(), mus_p4()[iMu].eta());
        //cout << "emu, mu FR, FR/(1-FR): " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      return -9999.;
    }
  } //emu case

  return -9999.;
}

//--------------------------------------------------------------------

void ossusy_looper::makeTree(char *prefix, bool doFakeApp, FREnum frmode ){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();


  char* dir = "";
  if     ( g_trig == e_lowpt  ) dir = "lowpt";
  else if( g_trig == e_highpt ) dir = "highpt";

  //Super compressed ntuple here
  char* frsuffix = "";
  if( doFakeApp ){
    if ( frmode == e_qcd   ) frsuffix = "_doubleFake";
    if ( frmode == e_wjets ) frsuffix = "_singleFake";
  }

  outFile   = new TFile(Form("../output/%s/%s/%s_smallTree%s.root",g_version,dir,prefix,frsuffix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in ossusy_looper.h
  outTree->Branch("acc_2010",        &acc_2010_,         "acc_2010/I");
  outTree->Branch("acc_highmet",     &acc_highmet_,      "acc_highmet/I");
  outTree->Branch("acc_highht",      &acc_highht_,       "acc_highht/I");
  outTree->Branch("hbhe",            &hbhe_,             "hbhe/I");
  outTree->Branch("jetid",           &jetid_,            "jetid/I");
  outTree->Branch("jetid30",         &jetid30_,          "jetid30/I");
  outTree->Branch("json",            &json_,             "json/I");
  outTree->Branch("htoffset",        &htoffset_,         "htoffset/F");
  outTree->Branch("htuncor",         &htuncor_,          "htuncor/F");
  outTree->Branch("njetsoffset",     &njetsoffset_,      "njetsoffset/I");
  outTree->Branch("njetsuncor",      &njetsuncor_,       "njetsuncor/I");
  outTree->Branch("costhetaweight",  &costhetaweight_,   "costhetaweight/F");
  outTree->Branch("weight",          &weight_,           "weight/F");
  outTree->Branch("trgeff",          &trgeff_,           "trgeff/F");
  outTree->Branch("pthat",           &pthat_,            "pthat/F");
  outTree->Branch("qscale",          &qscale_,           "qscale/F");
  outTree->Branch("smeff",           &smeff_,            "smeff/F");
  outTree->Branch("k",               &k_,                "k/F");
  outTree->Branch("mllgen",          &mllgen_,           "mllgen/F");
  outTree->Branch("nlep",            &nlep_,             "nlep/I");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("ngoodmu",         &ngoodmu_,          "ngoodmu/I");
  outTree->Branch("mull",            &mull_,             "mull/I");
  outTree->Branch("mult",            &mult_,             "mult/I");
  outTree->Branch("mullgen",         &mullgen_,          "mullgen/I");
  outTree->Branch("multgen",         &multgen_,          "multgen/I");
  outTree->Branch("proc",            &proc_,             "proc/I");
  outTree->Branch("leptype",         &leptype_,          "leptype/I");
  outTree->Branch("topmass",         &topmass_,          "topmass/F");
  outTree->Branch("dilmass",         &dilmass_,          "dilmass/F");
  outTree->Branch("tcmet",           &tcmet_,            "tcmet/F");

  outTree->Branch("trkmet",          &trkmet_,           "trkmet/F");
  outTree->Branch("trkmetphi",       &trkmetphi_,        "trkmetphi/F");
  outTree->Branch("trkmetproj",      &trkmetproj_,       "trkmetproj/F");

  outTree->Branch("trkmet4",         &trkmet4_,          "trkmet4/F");
  outTree->Branch("trkmet4phi",      &trkmet4phi_,       "trkmet4phi/F");
  outTree->Branch("trkmet4proj",     &trkmet4proj_,      "trkmet4proj/F");

  outTree->Branch("trkmet8",         &trkmet8_,          "trkmet8/F");
  outTree->Branch("trkmet8phi",      &trkmet8phi_,       "trkmet8phi/F");
  outTree->Branch("trkmet8proj",     &trkmet8proj_,      "trkmet8proj/F");


  outTree->Branch("tcmet00",         &tcmet00_,          "tcmet00/F");
  outTree->Branch("tcmet10",         &tcmet10_,          "tcmet10/F");
  outTree->Branch("tcmet20",         &tcmet20_,          "tcmet20/F");
  outTree->Branch("tcmet30",         &tcmet30_,          "tcmet30/F");
  outTree->Branch("tcmet40",         &tcmet40_,          "tcmet40/F");
  outTree->Branch("tcmet50",         &tcmet50_,          "tcmet50/F");
  outTree->Branch("genmet",          &genmet_,           "genmet/F");
  outTree->Branch("pfmet",           &pfmet_,            "pfmet/F");
  outTree->Branch("pfmetsig",        &pfmetsig_,         "pfmetsig/F");
  outTree->Branch("pfmetphi",        &pfmetphi_,         "pfmetphi/F");
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
  outTree->Branch("mt2",             &mt2_,              "mt2/F");  
  outTree->Branch("mt2j",            &mt2j_,             "mt2j/F");  
  outTree->Branch("mt2jcore",        &mt2jcore_,         "mt2jcore/F");  
  outTree->Branch("sumjetpt",        &sumjetpt_,         "sumjetpt/F");
  outTree->Branch("dileta",          &dileta_,           "dileta/F");
  outTree->Branch("dilpt",           &dilpt_,            "dilpt/F");
  outTree->Branch("dildphi",         &dildphi_,          "dildphi/F");
  outTree->Branch("njets",           &njets_,            "njets/I");
  outTree->Branch("ngenjets",        &ngenjets_,         "ngenjets/I");
  outTree->Branch("npfjets",         &npfjets_,          "npfjets/I");
  outTree->Branch("njpt",            &njpt_,             "njpt/I");
  outTree->Branch("npfjets25",       &npfjets25_,        "npfjets25/I");
  outTree->Branch("npfjets40",       &npfjets40_,        "npfjets40/I");
  outTree->Branch("npfjetspv",       &npfjetspv_,        "npfjetspv/I");
  outTree->Branch("njetsUp",         &njetsUp_,          "njetsUp/I");
  outTree->Branch("njetsDown",       &njetsDown_,        "njetsDown/I");
  outTree->Branch("htUp",            &htUp_,             "htUp/F");
  outTree->Branch("htDown",          &htDown_,           "htDown/F");
  outTree->Branch("nvtx",            &nvtx_,             "nvtx/I");
  outTree->Branch("ndavtx",          &ndavtx_,           "ndavtx/I");
  outTree->Branch("ndavtxweight",    &ndavtxweight_,     "ndavtxweight/F");
  outTree->Branch("nbtags",          &nbtags_,           "nbtags/I");
  outTree->Branch("nbtagstcl",       &nbtagstcl_,        "nbtagstcl/I");
  outTree->Branch("nbtagstcm",       &nbtagstcm_,        "nbtagstcm/I");
  outTree->Branch("vecjetpt",        &vecjetpt_,         "vecjetpt/F");
  outTree->Branch("pass",            &pass_,             "pass/I");
  outTree->Branch("passz",           &passz_,            "passz/I");
  outTree->Branch("m0",              &m0_,               "m0/F");
  outTree->Branch("m12",             &m12_,              "m12/F");
  outTree->Branch("id1",             &id1_,              "id1/I");
  outTree->Branch("id2",             &id2_,              "id2/I");
  outTree->Branch("w1",              &w1_,               "w1/I");
  outTree->Branch("w2",              &w2_,               "w2/I");
  outTree->Branch("iso1",            &iso1_,             "iso1/F");
  outTree->Branch("isont1",          &isont1_,           "isont1/F");
  outTree->Branch("iso2",            &iso2_,             "iso2/F");
  outTree->Branch("isont2",          &isont2_,           "isont2/F");
  outTree->Branch("ptl1",            &ptl1_,             "ptl1/F");
  outTree->Branch("ptl2",            &ptl2_,             "ptl2/F");
  outTree->Branch("ptj1",            &ptj1_,             "ptj1/F");
  outTree->Branch("ptj2",            &ptj2_,             "ptj2/F");
  outTree->Branch("etal1",           &etal1_,            "etal1/F");
  outTree->Branch("etal2",           &etal2_,            "etal2/F");
  outTree->Branch("phil1",           &phil1_,            "phil1/F");
  outTree->Branch("phil2",           &phil2_,            "phil2/F");
  outTree->Branch("meff",            &meff_,             "meff/F");
  outTree->Branch("mt",              &mt_,               "mt/F");
  outTree->Branch("dataset",         &dataset_,          "dataset[200]/C");
  outTree->Branch("run",             &run_,              "run/I");
  outTree->Branch("lumi",            &lumi_,             "lumi/I");
  outTree->Branch("event",           &event_,            "event/I");
  outTree->Branch("y",               &y_,                "y/F");  
  outTree->Branch("ht",              &ht_,               "ht/F");  
  outTree->Branch("htgen",           &htgen_,            "htgen/F");  
  outTree->Branch("htpf",            &htpf_,             "htpf/F");  
  outTree->Branch("htjpt",           &htjpt_,            "htjpt/F");  
  outTree->Branch("htpf25",          &htpf25_,           "htpf25/F");  
  outTree->Branch("htpf40",          &htpf40_,           "htpf40/F");  
  outTree->Branch("htpfpv",          &htpfpv_,           "htpfpv/F");  
  outTree->Branch("nels",            &nels_,             "nels/I");  
  outTree->Branch("nmus",            &nmus_,             "nmus/I");  
  outTree->Branch("ntaus",           &ntaus_,            "ntaus/I");  
  outTree->Branch("dphijm",          &dphijm_,           "dphijm/F");  
  outTree->Branch("ptjetraw",        &ptjetraw_,         "ptjetraw/F");  
  outTree->Branch("ptjet23",         &ptjet23_,          "ptjet23/F");  
  outTree->Branch("ptjetF23",        &ptjetF23_,         "ptjetF23/F");  
  outTree->Branch("ptjetO23",        &ptjetO23_,         "ptjetO23/F");  
  outTree->Branch("cosphijz",        &cosphijz_,         "cosphijz/F");  
  outTree->Branch("njets15",         &njets15_,          "njets15/I");  
 
  outTree->Branch("dilep"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &dilep_	);
  outTree->Branch("lep1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_	);
  outTree->Branch("lep2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_	);
  outTree->Branch("jet"	    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_	        );


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

bool ossusy_looper::jetFromSignalPV( int ijet , int vtxIdx , int beta_exponent ){

  vector<int> myGoodVertices       = goodDAVertices();
  const unsigned int nGoodVertices = myGoodVertices.size();

  if( nGoodVertices == 0 ){
    cout << "Didn't find any good vertices!" << endl;
    return false;
  }

  float betamax = -1;
  int   imax    = -1;

  for (vector<int>::iterator ivtx = myGoodVertices.begin(); ivtx != myGoodVertices.end() ; ++ivtx ){
    
    float beta = beta_jet_vtx( ijet , *ivtx , beta_exponent );

    if( beta > betamax ){
      betamax = beta;
      imax    = *ivtx;
    }

  }
  
  if( imax == vtxIdx )  return true;
  return false;
}

//--------------------------------------------------------------------

float ossusy_looper::dz_trk_vtx(const unsigned int trkidx, const unsigned int vtxidx){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

//--------------------------------------------------------------------

float ossusy_looper::beta_jet_vtx( int ijet , int vtxIdx , int beta_exponent ){

  //------------------------------------------
  // loop over PFCandidates matched to pfjet
  //------------------------------------------
  
  vector<int> pfcands = pfjets_pfcandIndicies().at(ijet);
  
  LorentzVector v_pfcands(0,0,0,0);

  vector<int> matchedTracks;
  matchedTracks.clear();

  float sumPtTot = 0.;

  for( vector<int>::iterator ipf = pfcands.begin() ; ipf < pfcands.end() ; ++ipf ){
     
    v_pfcands += pfcands_p4().at(*ipf);

    //---------------------------------------------
    // find track index of charged PFCandidates
    //---------------------------------------------
     
    if( pfcands_charge().at(*ipf) == 0 ) continue;

    int itrk = cms2.pfcands_trkidx().at(*ipf);
     
    //note: this should only happen for electrons which do not have a matched track, currently ignoring these guys
    if( itrk >= trks_trk_p4().size() || itrk < 0 )  continue;
     
    //----------------------------------
    // store indices of matched tracks
    //----------------------------------
     
    matchedTracks.push_back( itrk );

    sumPtTot += pow( trks_trk_p4().at(itrk).pt() , beta_exponent );

  }

  //-------------
  //sanity check
  //-------------
  
  if( fabs( pfjets_p4().at(ijet).pt() - v_pfcands.pt() ) > 0.1 ){
    cout << "Warning: pfjet pt " << pfjets_p4().at(ijet).pt() 
         << " doesn't match sum of PFCandidates pt " << v_pfcands.pt() << endl;
  }

  //----------------------------------
  // find good vertices
  //----------------------------------

  vector<int> myGoodVertices       = goodDAVertices();
  const unsigned int nGoodVertices = myGoodVertices.size();

  if( nGoodVertices == 0 ){
    cout << "Didn't find any good vertices!" << endl;
    return -1.;
  }

  float beta = 0.;

  //------------------------------------------------------------------------------------------------
  // loop over tracks. if track_i is closest to signal PV, add pow(track_i pt,beta_exponent) to beta
  //------------------------------------------------------------------------------------------------

  for (vector<int>::iterator itrk = matchedTracks.begin(); itrk != matchedTracks.end(); ++itrk) {
     
    float mindz = 1000;
    int   vtxi  = -1;
     
    for (vector<int>::iterator ivtx = myGoodVertices.begin(); ivtx != myGoodVertices.end() ; ++ivtx ){
       
      float thisdz = dz_trk_vtx(*itrk,*ivtx);
       
      if ( fabs( thisdz ) < fabs( mindz ) ) {
        mindz = thisdz;
        vtxi  = *ivtx;
      }
    }
     
    if( vtxi == vtxIdx ){
      beta += pow( trks_trk_p4().at(*itrk).pt() , beta_exponent ) / sumPtTot;
    }
  }
 
  return beta;

}


std::pair<float,float> ossusy_looper::PFCandidateMET(const unsigned int vtxIdx, const unsigned int hypIdx, vector<int> goodjets, 
						     float dz_thresh, float pt_thresh, float etacut , bool usePFCandidatePt , bool correctJets ){
  
  float tmet_x = 0.;
  float tmet_y = 0.;
  
  //------------------------------------------------
  // start by adding hypothesis leptons to tmet
  //------------------------------------------------
  
  if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
    tmet_x -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].px();
    tmet_y -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].py();
  }
  else if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
    tmet_x -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].px();
    tmet_y -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].py();
  }
  
  if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
    tmet_x -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].px();
    tmet_y -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].py();
  }
  else if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
    tmet_x -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].px();
    tmet_y -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].py();
  }
  
  //------------------------------------------------
  // next, add selected jets to met
  //------------------------------------------------
  
  if( correctJets ){

    for( vector<int>::iterator igoodjet = goodjets.begin() ; igoodjet < goodjets.end() ; ++igoodjet ){
      
      LorentzVector vjet = pfjets_p4().at(*igoodjet) * pfjets_corL1FastL2L3().at(*igoodjet);
      
      tmet_x -= vjet.px();
      tmet_y -= vjet.py();
      
    }
    
  }
  
  //---------------------------------------------------
  // loop over PFCandidates
  //---------------------------------------------------

  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

    //----------------------------------------------------------------------
    // exclude particles belonging to jets that were already corrected for
    //----------------------------------------------------------------------

    if( correctJets ){
      
      bool skipPFCandidate = false;
      
      for( vector<int>::iterator jet_it = goodjets.begin() ; jet_it != goodjets.end() ; jet_it++ ){
        
        vector<int> jetConstituents = pfjets_pfcandIndicies().at(*jet_it);
        
        for( vector<int>::iterator pf_it = jetConstituents.begin() ; pf_it != jetConstituents.end() ; pf_it++ ){
          if( ipf == *pf_it ) skipPFCandidate = true;
        }
        
      }
    
      if( skipPFCandidate ) continue;
      
    }

    //--------------------
    //deal with neutrals
    //--------------------
  
    if( cms2.pfcands_charge().at(ipf) == 0 ){
      
      //--------------------
      // pt, eta cuts
      //--------------------
      
      if( pfcands_p4().at(ipf).pt()          < pt_thresh ) continue;
      if( fabs( pfcands_p4().at(ipf).eta() ) > etacut    ) continue;
      
      if( dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_lt_p4()[hypIdx]) < 0.1)  continue;
      if( dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_ll_p4()[hypIdx]) < 0.1)  continue;
      
      //-----------------------------------------
      // correct MET for neutral PFCandidate
      //-----------------------------------------
      
      tmet_x -= pfcands_p4().at(ipf).px();
      tmet_y -= pfcands_p4().at(ipf).py();
      
    }
    
    //-------------------------------
    // deal with charged particles
    //-------------------------------
    
    else{
      
      //------------------------------------
      // get track matched to PFCandidate
      //------------------------------------
    
      int itrk = cms2.pfcands_trkidx().at(ipf);
    
      if( itrk >= trks_trk_p4().size() || itrk < 0 ){
        //note: this should only happen for electrons which do not have a matched track
        //currently we are just ignoring these guys
        continue;
      }
    
      //----------------------------------------
      // find closest PV and dz w.r.t. that PV
      //----------------------------------------
    
      float mindz = 999.;
      int vtxi    = -1;
      
      for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
       
        float mydz = dz_trk_vtx(itrk,ivtx);
        //float mydz = cms2.vtxs_position()[ivtx].z() - cms2.trks_vertex_p4()[itrk].z();
         
        if (fabs(mydz) < fabs(mindz)) {
          mindz = mydz;
          vtxi = ivtx;
        }
         
      }
    
      //----------------------------------------------------------------------------
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
      //----------------------------------------------------------------------------
    
      if ( vtxi != vtxIdx )                                                               continue;
      if ( fabs(mindz) > dz_thresh )                                                      continue;
      if ( dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_lt_p4()[hypIdx]) < 0.1 )   continue;
      if ( dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_ll_p4()[hypIdx]) < 0.1 )   continue;
    

      //---------------------------------------
      // correct MET for charged PFCandidate
      //---------------------------------------

      if( usePFCandidatePt ){
        tmet_x -= cms2.pfcands_p4()[ipf].px();
        tmet_y -= cms2.pfcands_p4()[ipf].py();
      }
    
      else{
        tmet_x -= cms2.trks_trk_p4()[itrk].px();
        tmet_y -= cms2.trks_trk_p4()[itrk].py();
      }
    
    } 
  }// end loop over tracks
  
  float met = sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
  float metphi = atan2( tmet_y , tmet_x );
  
  return make_pair( met , metphi );
}
