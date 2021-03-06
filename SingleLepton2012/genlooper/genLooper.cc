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
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/LorentzVector.h"
#include "genLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"

#include "../CORE/CMS2.h"
#include "../CORE/utilities.h"
#include "../CORE/ssSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/MITConversionUtilities.h"
#include "../CORE/muonSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetcorr/FactorizedJetCorrector.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"
#include "../CORE/jetSelections.h"
#include "../CORE/photonSelections.h"
#include "../CORE/triggerUtils.h"
#include "../CORE/triggerSuperModel.h"
#include "../CORE/mcSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/mcSUSYkfactor.h"
#include "../CORE/SimpleFakeRate.h"
#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
#include "../Tools/msugraCrossSection.h"
//#include "BtagFuncs.h"
//#include "../Tools/bTagEff_BTV.h"

bool verbose              = false;
bool doTenPercent         = false;
bool vetoTransition       = true;
bool useOldIsolation      = false;

//#include "../CORE/topmass/getTopMassEstimate.icc" // REPLACETOPMASS
//#include "../CORE/triggerUtils.cc"

using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > P4;
typedef vector< P4 > VofP4;
//typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

struct indP4{
  LorentzVector p4obj;
  int p4ind;
};

typedef vector< indP4 > VofiP4;

inline bool sortIP4ByPt(indP4 iP41, indP4 iP42) {
  return iP41.p4obj.pt() > iP42.p4obj.pt();
}


//mSUGRA scan parameters-----------------------------

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
float getdltrigweight(int id1, int id2);
float getsltrigweight(int id1, float pt, float eta);
//--------------------------------------------------------------------

void genLooper::weight3D_init( std::string WeightFileName ) { 

  TFile *infile = new TFile(WeightFileName.c_str());
  TH1F *WHist = (TH1F*)infile->Get("WHist");

  // Check if the histogram exists           
  if (!WHist) {
    cout << "Error, could not find the histogram WHist in the file "
	 << WeightFileName << ", quitting" << endl;
    exit(0);
  }

  for (int i=0; i<50; i++) 
    for(int j=0; j<50; j++)
      for(int k=0; k<50; k++) {
	Weight3D[i][j][k] = WHist->GetBinContent(i+1,j+1,k+1);
      }

  cout << " 3D Weight Matrix initialized! " << endl;

  delete infile;

  return;


}

//--------------------------------------------------------------------

double genLooper::weight3D( int pv1, int pv2, int pv3 ) {

  using std::min;

  int npm1 = min(pv1,49);
  int np0 = min(pv2,49);
  int npp1 = min(pv3,49);

  return Weight3D[npm1][np0][npp1];

}

//--------------------------------------------------------------------

void checkElectron( int elidx ){

  cout << "Check electron" << endl;
  cout << "Pass all    " << pass_electronSelection( elidx , electronSelection_ssV5			) << endl;
  cout << "Pass ID     " << pass_electronSelection( elidx , electronSelection_ssV5_noIso		) << endl;
  cout << "Pass iso    " << pass_electronSelection( elidx , electronSelection_ssV5_iso	        	) << endl;
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

bool isGenBMatched ( LorentzVector p4, float dR ) {

  //For now only checking status 3 in dR
    for (unsigned int igen = 0; igen < cms2.genps_p4().size(); igen++) {
      
      int id = cms2.genps_id().at(igen);
      if( abs(id)!=5 ) continue;
      
      if( dRbetweenVectors( p4 , cms2.genps_p4().at(igen) ) < dR ) return true;

    } 
    return false;
}

//--------------------------------------------------------------------

int isGenQGMatched ( LorentzVector p4, float dR ) {
  //Start from the end that seems to have the decay products of the W first
  for (int igen = (cms2.genps_p4().size()-1); igen >-1; igen--) {
    float deltaR = dRbetweenVectors( p4 , cms2.genps_p4().at(igen) );
    if ( deltaR > dR ) continue;
    int id = cms2.genps_id().at(igen);
    int mothid = genps_id_mother().at(igen);
    // cout<<"status 3 particle ID "<<id<<" mother "<<mothid
    // 	<<" dR to jet "<<deltaR<<endl;
    if (abs(id)<6 && abs(mothid)==24) 
      return (mothid>0) ? 2 : -2;
    if (abs(id)==5 && abs(mothid)==6) 
      return (mothid>0) ? 1 : -1;
    if (abs(id)==21) return 3;
    if (abs(id)<6) return 4;
  }
  return -9;
}

//--------------------------------------------------------------------

float dRGenJet ( LorentzVector p4, float genminpt = 20.) {

  //return dR to closest gen-jet with pT > genminpt
  float mindeltaR = 9999.;
  for (unsigned int igen = 0; igen < cms2.genjets_p4().size(); igen++) {
    LorentzVector vgenj = cms2.genjets_p4().at(igen);
    if ( vgenj.Pt() < genminpt ) continue;
    float deltaR = dRbetweenVectors( p4 , vgenj );
    if ( deltaR< mindeltaR ) mindeltaR = deltaR;
  }
  return mindeltaR;

}

//--------------------------------------------------------------------

int isSSVMTagged ( int ijet ) {

  if ( pfjets_simpleSecondaryVertexHighEffBJetTag().at(ijet) > 1.74 )
    return 1;

  return 0;

}

//--------------------------------------------------------------------

int isCSVTagged ( int ijet ) {
  
  float discrim = pfjets_combinedSecondaryVertexBJetTag().at(ijet);
  if ( discrim > 0.679 ) return 1;//medium
  if ( discrim > 0.244 ) return 2;//loose

  return 0;

}

//--------------------------------------------------------------------

bool isBTagged ( LorentzVector p4, VofP4 bJets ) {

  for( int ijet = 0 ; ijet < (int)bJets.size() ; ijet++ ){
    if( dRbetweenVectors( p4 , bJets.at(ijet) ) < 0.4 ) return true;
  }

  return false;

}

//--------------------------------------------------------------------

int getLeptonMatchIndex ( LorentzVector *jet, LorentzVector *lep1, LorentzVector *lep2, float dR ) {
  if (lep1 && dRbetweenVectors( *lep1 , *jet ) < dR) return 1;
  if (lep2 && dRbetweenVectors( *lep2 , *jet ) < dR) return 2;
  return -1;
  
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

bool objectPassTrigger(const LorentzVector &obj, char* trigname, float drmax = 0.1 ){

  TString exact_trigname = triggerName( trigname );

  if( exact_trigname.Contains("TRIGGER_NOT_FOUND") ){
    cout << __FILE__ << " " << __LINE__ << " Error! couldn't find trigger name " << trigname << endl;
    return false;
  }

  std::vector<LorentzVector> trigp4 = cms2.hlt_trigObjs_p4()[findTriggerIndex(exact_trigname)];

  if( trigp4.size() == 0 ) return false;

  for (unsigned int i = 0; i < trigp4.size(); ++i){
    float dr = dRbetweenVectors(trigp4[i], obj);
    if( dr < drmax ) return true;
  }

  return false;
}

//--------------------------------------------------------------------

genLooper::genLooper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
  initialized = false;
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
  int index     = (int)((m0 - m0min) / binsize);
  return index;
    
}

//--------------------------------------------------------------------

int getIndexFromM12(float m12){
  
  float binsize = (m12max - m12min) / (float) nm12points;
  int index     = (int)((m12 - m12min) / binsize);
  return index;
    
}

//--------------------------------------------------------------------

void genLooper::InitBaby(){

  //pdf variables
  pdfid1_ = -999;
  pdfid2_ = -999;
  pdfQ_   = -99999.;
  pdfx1_  = -99999.;
  pdfx2_  = -99999.;

  mutrigweight_ = 1.;
  sltrigweight_ = 1.;
  dltrigweight_ = 1.;

  jetid_	= 1;
  jetid30_	= 1;

  // btag variables
  nbtagsssv_	= 0;
  nbtagstcl_	= 0;
  nbtagstcm_	= 0;
  nbtagscsvl_   = 0;
  nbtagscsvm_   = 0;
  nbtagscsvt_   = 0;
  nbtagsssvcorr_    = 0;
  nbtagstclcorr_    = 0;
  nbtagstcmcorr_    = 0;
  nbtagscsvlcorr_   = 0;
  nbtagscsvmcorr_   = 0;
  nbtagscsvtcorr_   = 0;

  // njets with JEC variation
  njetsUp_	= 0;
  njetsDown_	= 0;
  ht_           = -999.;
  htUp_		= -999.;
  htDown_	= -999.;

  // Type1 pfmet
  t1met10_	=-999.;
  t1met20_	=-999.;
  t1met30_	=-999.;
  t1met10phi_	=-999.;
  t1met20phi_	=-999.;
  t1met30phi_	=-999.;
  t1met10mt_	=-999.;
  t1met20mt_	=-999.;
  t1met30mt_	=-999.;
  lepmetpt_	=-999.;
  lept1met10pt_	=-999.;
  
  //trkmet
  trkmet_              =-999.;
  trkmetphi_           =-999.;
  trkmet_nolepcorr_    =-999.;
  trkmetphi_nolepcorr_ =-999.;

  //phi corrected type1 mets
  t1metphicorr_	      =-999.;
  t1metphicorrphi_    =-999.;
  t1metphicorrlep_    =-999.;
  t1metphicorrlepphi_ =-999.;
  t1metphicorrmt_     =-999.;
  t1metphicorrlepmt_  =-999.;
  
  // pfjet vars
  npfjets30_	= 0;
  npfjets35_	= 0;
  npfjets40_	= 0;
  npfjets45_	= 0;
  npfjets30lepcorr_ = 0;
  knjets_       = 1.;

  htpf30_	= 0.;
  htpf35_	= 0.;
  htpf40_	= 0.;
  htpf45_	= 0.;
  htpfres30_	= 0.;
  htpfres35_	= 0.;
  htpfres40_	= 0.;
  htpfres45_	= 0.;

  //iso trk vars
  trkpt5_ 	    = -999.;
  mleptrk5_ 	    = -999.;
  trkreliso5_ 	    = -999.;
  trkpt10_ 	    = -999.;
  mleptrk10_ 	    = -999.;
  trkreliso10_ 	    = -999.;
  trkpt5loose_ 	    = -999.;
  trkreliso5loose_  = -999.;
  trkpt10loose_     = -999.;
  trkreliso10loose_ = -999.;

  // MC truth info
  mcid1_	= -1;
  mcid2_	= -1;
  mclep1_	=  0;
  mclep2_	=  0;
  mctaud1_      =  0;
  mctaud2_      =  0;
  mctaudvis1_   =  0;
  mctaudvis2_   =  0;
  mcdecay1_	= -1;
  mcdecay2_	= -1;
  mcdr1_	= -1;
  mcdr2_	= -1;
  
  mlepid_       = -1;
  mlep_         =  0;
  mleppassid_   = -1;
  mleppassiso_  = -1;
  mlepiso_      = -1.0;

  mllgen_	= -1;
  pthat_	= -1;
  qscale_	= -1;
  genmet_       = -9999;
  gensumet_     = -9999;
  genmetphi_    = -9999;   
  m0_		= -9999;
  m12_		= -9999;
  mG_		= -9999;
  mL_		= -9999;
  x_		= -9999;
  ksusy_	= -999;
  ksusyup_	= -999;
  ksusydn_	= -999;
  xsecsusy_	= -999;
  xsecsusy2_	= -999;
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
  jet_		= 0;

  // jet p4's
  pfjet1_	= 0;
  pfjet2_	= 0;
  pfjet3_	= 0;
  pfjet4_	= 0;
  pfjet5_	= 0;
  pfjet6_	= 0;

  bjet1_      = -999; 
  bjet2_      = -999; 
  bjet3_      = -999; 
  bjet4_      = -999; 
  bjet5_      = -999; 
  bjet6_      = -999; 
  lepjet1_    = -999; 
  lepjet2_    = -999; 
  lepjet3_    = -999; 
  lepjet4_    = -999; 
  lepjet5_    = -999; 
  lepjet6_    = -999; 
  qgjet1_     = -999; 
  qgjet2_     = -999; 
  qgjet3_     = -999; 
  qgjet4_     = -999; 
  qgjet5_     = -999; 
  qgjet6_     = -999; 
  genjetdr1_  = 9999.; 
  genjetdr2_  = 9999.; 
  genjetdr3_  = 9999.; 
  genjetdr4_  = 9999.; 
  genjetdr5_  = 9999.; 
  genjetdr6_  = 9999.; 

  lep1chi2ndf_	= -9999.;
  lep2chi2ndf_	= -9999.;
  lep1dpt_	= -9999.;
  lep2dpt_	= -9999.;
  leptype1_	= -999;
  leptype2_	= -999;
  lep1_		= 0;
  lep2_		= 0;
  trklep1_	= 0;
  trklep2_	= 0;
  gfitlep1_	= 0;
  gfitlep2_	= 0;
  lepp_		= 0;
  lepm_		= 0;
  pflep1_	= 0;
  pflep2_	= 0;
  leppfjet1_	= 0;
  leppfjet2_	= 0;
  pflep_        = 0;
  pftaud_       = 0;
  mbb_		= -9999.;
  mcmln_	= -9999.;
  mcmtln_	= -9999.;
  pflepmindrj_   = 9999.;
  pftaudmindrj_  = 9999.;
  lep1pfjetdr_   = 9999.;
  lep2pfjetdr_   = 9999.;

  pfcand5_        = 0;
  pfcand10_       = 0;
  pfcandiso5_     = 9999.;     
  pfcandiso10_    = 9999.;     
  pfcandpt5_      = 9999.;
  pfcandpt10_     = 9999.;
  pfcandmindrj5_  = 9999.;
  pfcandmindrj10_ = 9999.;

  trkpt10pt0p1_	    = 9999.;
  trkreliso10pt0p1_ = 9999.;
  trkpt10pt0p2_	    = 9999.;
  trkreliso10pt0p2_ = 9999.;
  trkpt10pt0p3_	    = 9999.;
  trkreliso10pt0p3_ = 9999.;
  trkpt10pt0p4_	    = 9999.;
  trkreliso10pt0p4_ = 9999.;
  trkpt10pt0p5_	    = 9999.;
  trkreliso10pt0p5_ = 9999.;
  trkpt10pt0p6_	    = 9999.;
  trkreliso10pt0p6_ = 9999.;
  trkpt10pt0p7_	    = 9999.;
  trkreliso10pt0p7_ = 9999.;
  trkpt10pt0p8_	    = 9999.;
  trkreliso10pt0p8_ = 9999.;
  trkpt10pt0p9_	    = 9999.;
  trkreliso10pt0p9_ = 9999.;
  trkpt10pt1p0_	    = 9999.;
  trkreliso10pt1p0_ = 9999.;

  pfcandpt10pt0p1_  = 9999.;
  pfcandiso10pt0p1_ = 9999.;
  pfcandpt10pt0p2_  = 9999.;
  pfcandiso10pt0p2_ = 9999.;
  pfcandpt10pt0p3_  = 9999.;
  pfcandiso10pt0p3_ = 9999.;
  pfcandpt10pt0p4_  = 9999.;
  pfcandiso10pt0p4_ = 9999.;
  pfcandpt10pt0p5_  = 9999.;
  pfcandiso10pt0p5_ = 9999.;
  pfcandpt10pt0p6_  = 9999.;
  pfcandiso10pt0p6_ = 9999.;
  pfcandpt10pt0p7_  = 9999.;
  pfcandiso10pt0p7_ = 9999.;
  pfcandpt10pt0p8_  = 9999.;
  pfcandiso10pt0p8_ = 9999.;
  pfcandpt10pt0p9_  = 9999.;
  pfcandiso10pt0p9_ = 9999.;
  pfcandpt10pt1p0_  = 9999.;
  pfcandiso10pt1p0_ = 9999.;

  //lepton variables
  iso1_   = -9999; 
  isont1_ = -9999;
  isopfold1_ = -9999;
  isopf1_ = -9999;
  etasc1_ = -9999;
  eoverpin_  = -9999;
  eoverpout_ = -9999;
  dEtaIn_ = -9999;
  dPhiIn_ = -9999;
  sigmaIEtaIEta_ = -9999;
  hOverE_ = -9999;
  ooemoop_ = -9999;
  d0vtx_ = -9999;
  dzvtx_ = -9999;
  expinnerlayers_ = -9999;
  fbrem_ = -9999;
  pfisoch_ = -9999;
  pfisoem_ = -9999;
  pfisonh_ = -9999;
  eSC_ = -9999;
  phiSC_ = -9999;
  eSCRaw_ = -9999;
  eSCPresh_ = -9999;  
  ecalveto1_ = -9999;
  hcalveto1_ = -9999;

  iso2_   = -9999;
  isont2_ = -9999;
  isopf2_ = -9999;
  etasc2_ = -9999;
  eoverpin2_  = -9999;
  eoverpout2_ = -9999;
  dEtaIn2_ = -9999;
  dPhiIn2_ = -9999;
  sigmaIEtaIEta2_ = -9999;
  hOverE2_ = -9999;
  ooemoop2_ = -9999;
  d0vtx2_ = -9999;
  dzvtx2_ = -9999;
  expinnerlayers2_ = -9999;
  fbrem2_ = -9999;
  pfisoch2_ = -9999;
  pfisoem2_ = -9999;
  pfisonh2_ = -9999;
  eSC2_ = -9999;
  phiSC2_ = -9999;
  eSCRaw2_ = -9999;
  eSCPresh2_ = -9999;  
  ecalveto2_ = -9999;
  hcalveto2_ = -9999;
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

void genLooper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

//--------------------------------------------------------------------

float genLooper::stopPairCrossSection( float stopmass ){

  int   bin  = stop_xsec_hist->FindBin(stopmass);
  float xsec = stop_xsec_hist->GetBinContent(bin);
  return xsec;

}

//--------------------------------------------------------------------

pair<float,float> Type1PFMET( VofP4 jets_p4 , vector<float> cors , vector<float> l1cors , float minpt ){

  float metx = evt_pfmet() * cos( evt_pfmetPhi() );
  float mety = evt_pfmet() * sin( evt_pfmetPhi() );

  assert( jets_p4.size() == cors.size() );

  for( unsigned int i = 0 ; i < jets_p4.size() ; ++i ){
    float corrpt = jets_p4.at(i).pt() * cors.at(i);
    if( corrpt < minpt ) continue;
    float l1corr = (l1cors.size()==0) ? 1. : l1cors.at(i);
    metx += jets_p4.at(i).px() * l1corr - jets_p4.at(i).px() * cors.at(i);
    mety += jets_p4.at(i).py() * l1corr - jets_p4.at(i).py() * cors.at(i);
  }

  pair<float, float> type1met = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return type1met;
}

//--------------------------------------------------------------------

float getMT( float leppt , float lepphi , float met , float metphi ) {
  float dphi = fabs( lepphi - metphi );
      if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
      return sqrt( 2 * ( leppt * met * (1 - cos( dphi ) ) ) );
}

//--------------------------------------------------------------------

bool passSingleMuTrigger2011_pt30( bool isData , int lepType ) {
  
  //----------------------------
  // single muon triggers
  //----------------------------

  // no triggers required for MC
  //if( !isData ) return true;

  // false for electron channel
  if( lepType == 0 ){
    return false;
  }

  // muon channel
  else if( lepType == 1 ){    
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_v") )          return true; //  < 173212
    if( passUnprescaledHLTTriggerPattern("HLT_IsoMu30_eta2p1_v") )   return true; // >= 173212
  }

  else{
    cout << __FILE__ << " " << __LINE__ << " ERROR unrecognized lepType " << lepType << ", quitting" << endl;
    exit(0);
  }

  return false;
}

//--------------------------------------------------------------------

float getMuTriggerWeight( float pt, float eta ) {
  //Trigger efficiency for single muon triggers averaged over full 2011 dataset
  //From AN2011-456 Table 28
  float trigweights[2][3] = {{0.9002, 0.8352, 0.8266},
			     {0.9440, 0.8821, 0.8611}};
  int i_pt = -1;
  if ( pt > 30 && pt < 40 ) i_pt = 0;
  else if ( pt > 40 ) i_pt = 1;
  if ( i_pt < 0 ) return 1.;

  int i_eta = -1;
  if ( fabs(eta) < 0.8 ) i_eta = 0;
  else if ( fabs(eta) >= 0.8 && abs(eta) < 1.5 ) i_eta = 1;
  else if ( fabs(eta) >= 1.5 && abs(eta) < 2.1 ) i_eta = 2;
  if ( i_eta < 0 ) return 1.;

  return trigweights[i_pt][i_eta];

}


//--------------------------------------------------------------------

float getMuTriggerWeightNew( float pt, float eta ) {

  float aeta = fabs(eta);

  if( aeta < 0.8 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.89;
    if( pt >  40.0 && pt <   60.0 ) return 0.89;
    if( pt >  60.0 && pt <  100.0 ) return 0.88;
    if( pt > 100.0                ) return 0.87;
  }
  
  else if( aeta > 0.8 && aeta < 1.3 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.81;
    if( pt >  40.0 && pt <   60.0 ) return 0.82;
    if( pt >  60.0 && pt <  100.0 ) return 0.81;
    if( pt > 100.0                ) return 0.80;
  }

  else if( aeta > 1.3 && aeta < 1.8 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.82;
    if( pt >  40.0 && pt <   60.0 ) return 0.84;
    if( pt >  60.0 && pt <  100.0 ) return 0.82;
    if( pt > 100.0                ) return 0.80;
  }

  else if( aeta > 1.8 && aeta < 2.0 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.79;
    if( pt >  40.0 && pt <   60.0 ) return 0.81;
    if( pt >  60.0 && pt <  100.0 ) return 0.80;
    if( pt > 100.0                ) return 0.79;
  }

  else if( aeta > 2.0 && aeta < 2.1 ){
    if( pt >  30.0 && pt <   40.0 ) return 0.69;
    if( pt >  40.0 && pt <   60.0 ) return 0.71;
    if( pt >  60.0 && pt <  100.0 ) return 0.70;
    if( pt > 100.0                ) return 0.70;
  }


  return 1;
}

//--------------------------------------------------------------------

float getminjdr( VofiP4 jets, LorentzVector *particle ) {
  float mindr = 9999.;
  if (jets.size()==0 || particle==0) return mindr;
  for ( unsigned int ijet = 0; ijet<jets.size(); ++ijet ) {
    float partjdr = dRbetweenVectors(jets.at(ijet).p4obj,*particle);
    if ( partjdr<mindr ) mindr = partjdr;
  }
  return mindr;
}

//--------------------------------------------------------------------

int genLooper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale, float lumi,
				  FREnum frmode, bool doFakeApp)

{

  bool isLM = TString(prefix).Contains("LM");
  bool isData = false;
  if( TString(prefix).Contains("data") || TString(prefix).Contains("2012") 
      || TString(prefix).Contains("dimu") || TString(prefix).Contains("diel")
      || TString(prefix).Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
  }
  if( doTenPercent ) cout << "Processing 10% of MC" << endl;

  //------------------------------------------------------------------------------------------------------
  // set json, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){

    //set json
    //    cout << "setting json " << g_json << endl;
    //    set_goodrun_file( g_json );

    //    if( TString(prefix).Contains("ttall_massivebin") ) 
    //    set_vtxreweight_rootfile("../vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root",true);

    //    weight3D_init( "vtxreweight/Weight3D.root" );

    //set msugra cross section file
    //    set_msugra_file("../goodModelNames_tanbeta10.txt");


    initialized = true;
  }

  //------------------------------------------------
  // set stop cross section file
  //------------------------------------------------

  stop_xsec_file = TFile::Open("stop_xsec.root");
  
  if( !stop_xsec_file->IsOpen() ){
    cout << "Error, could not open stop cross section TFile, quitting" << endl;
    exit(0);
  }
  
  stop_xsec_hist        = (TH1D*) stop_xsec_file->Get("h_stop_xsec");
  
  if( stop_xsec_hist == 0 ){
    cout << "Error, could not retrieve stop cross section hist, quitting" << endl;
    exit(0);
  }

  // instanciate topmass solver REPLACETOPMASS
  //ttdilepsolve * d_llsol = new ttdilepsolve;

  //instantiate SimpleFakeRate class for electrons and muons
  //this is the default, can change it below if needed
  SimpleFakeRate* mufr = 0;
  SimpleFakeRate* elfr = 0;

  if(doFakeApp) {

    cout << "NOT CURRENTLY SET UP FOR FAKES!!!!! QUITTING!" << endl;
    exit(0);

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

  int nSkip_els_conv_dist = 0;

  float netot  = 0.;
  float nmtot  = 0.;
  float nepass = 0.;
  float nmpass = 0.;

  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    cout << currentFile->GetTitle() << endl;

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

      // PrintTriggers();
      // exit(0);

      if( cms2.evt_ww_rho_vor() != cms2.evt_ww_rho_vor() ){
	cout << "Skipping event with rho = nan!!!" << endl;
	continue;
      }

      InitBaby();

      isdata_ = isData ? 1 : 0;

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      TString datasetname(evt_dataset().at(0));
      bool isperiodA = datasetname.Contains("2011A") ? true : false;
      //      cout<<"dataset: "<<datasetname.Data()<<" isperiodA: "<<isperiodA<<endl;


      
      //--------------------------------
      // get MC quantities
      //--------------------------------
      
      int nels       =  0;
      int nmus       =  0;
      int ntaus      =  0;
      int nleps      =  0;
      float dilptgen = -1;

      ptwgen_   = -1;
      ptzgen_   = -1;
      ptttbar_  = -1;
      ptt_      = -1;
      pttbar_   = -1;
      mttbar_   = -1;
      etattbar_ = -999;
      t_        = 0;
      tbar_     = 0;
      ttbar_    = 0;
      mcStop1_ = 0;   
      mcStop2_ = 0;   
      neutralino_t_ = 0;   
      neutralino_tbar_ = 0; 
      W_t_ = 0;   
      W_tbar_ = 0; 
      b_     = 0;
      bbar_    = 0;
      mcid1_      = -1;
      mcid2_      = -1;
      mclep1_     =  0;
      mclep2_     =  0;

      npartons_    =  0;
      nwzpartons_  = -9;
      maxpartonpt_ = -1;

      mgcor_ = 1.0;
      wflav_ = -1;

      if( !isData ){

	bool foundwz = false;
	//	w1_     = leptonOrTauIsFromW( index1 , id1_ , isLM );
	pthat_  = genps_pthat();
	qscale_ = genps_qScale();
	
	//store W flavor history
	if (!TString(prefix).Contains("mcatnlo"))
	  wflav_ = (int)genps_flavorHistoryFilterResult();

	//splitting ttbar into ttdil/ttotr
	//nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
	nleps = leptonGenpCount(nels, nmus, ntaus);
	
	nels_  = nels;
	nmus_  = nmus;
	ntaus_ = ntaus;
	nleps_ = nleps;

	// this is a weight which corrects for the wrong MG W->lnu BF
	if( TString(prefix).Contains("ttall") ||
	    TString(prefix).Contains("tt_") ){
	  if( nleps == 0 ) mgcor_ = 1.028;
	  if( nleps == 1 ) mgcor_ = 0.986;
	  if( nleps == 2 ) mgcor_ = 0.945;
	}
	if( TString(prefix).Contains("powheg") ||
	    TString(prefix).Contains("sherpa") ) 
	  mgcor_ = 1.0;

	if( strcmp(prefix,"ttem")  == 0 && ( nels + nmus ) != 2 ) continue;
	if( strcmp(prefix,"ttdil") == 0 && nleps != 2           ) continue;
	if( strcmp(prefix,"ttotr") == 0 && nleps == 2           ) continue;
	
	LorentzVector vdilepton(0,0,0,0);
	LorentzVector vttbar(0,0,0,0);
	int ntops = 0;

	for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 
	  if ( abs( cms2.genps_id().at(igen) ) == 11) vdilepton += genps_p4().at(igen); 
	  if ( abs( cms2.genps_id().at(igen) ) == 13) vdilepton += genps_p4().at(igen); 

	  int id = cms2.genps_id().at(igen);
	  int mothid = genps_id_mother().at(igen);

	  if( id == 6 ){
	    t_         = &(genps_p4().at(igen));
	    ptt_       = genps_p4().at(igen).pt();
	    vttbar    += genps_p4().at(igen);
	    ntops++;
	  }
	  if( id == -6 ){
	    tbar_      = &(genps_p4().at(igen));
	    pttbar_    = genps_p4().at(igen).pt();
	    vttbar    += genps_p4().at(igen); 
	    ntops++;
	  }

	  //store stop
	  if ( id == 1000006)
	    mcStop1_ = &(genps_p4().at(igen));   
	  else if ( id == -1000006 )
	    mcStop2_ = &(genps_p4().at(igen));   

    //store neutralino
    if ( genps_id_mother().at(igen) == 1000006  && ( abs(id) == 1000022 ) ) {
      neutralino_t_ = &(genps_p4().at(igen));
    }
    if ( genps_id_mother().at(igen) == -1000006 && ( abs(id) == 1000022 ) ) {
      neutralino_tbar_ = &(genps_p4().at(igen));
    }


	  //store b
	  if( id == 5 ){
	    b_         = &(genps_p4().at(igen));
	  }

	  if( id == -5 ){
	    bbar_      = &(genps_p4().at(igen));
	  }



	  //store W
	  if( id == 24 ){
	    W_t_         = &(genps_p4().at(igen));
	  }

	  if( id == -24 ){
	    W_tbar_      = &(genps_p4().at(igen));
	  }


	  //store daughter lepton
	  if ( abs(mothid) == 24 && (abs(id) == 11 || abs(id) == 13 || abs(id) ==15)) {

	    if (genps_id_mother().at(igen)>0) {
	      // lept 1 is the particle 
	      mcid1_ = genps_id().at(igen);
	      mclep1_ = &(genps_p4().at(igen));
	    } else {
	      // lept 2 is the anti-particle
	      mcid2_ = genps_id().at(igen);
	      mclep2_ = &(genps_p4().at(igen));
	    }
	  }

	  // store W or Z pT 
	  // ignoring cases where have more than 1 boson for now
	  if ( abs(id) == 24 ) {
	    ptwgen_ = genps_p4().at(igen).pt();
	    foundwz = true;
	    nwzpartons_  = 0;
	  }
	  if ( abs(id) == 23 ) {
	    ptzgen_ = genps_p4().at(igen).pt();
	    foundwz = true;
	    nwzpartons_  = 0;
	  }
  
	  double pid=abs(id);

	  if (foundwz && ( pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5 || pid == 6 || pid == 21 ) )   
	    nwzpartons_++;

	  // skip lines up to t and tbar
	  if( igen < 8 ) continue;

	  // require particle is a quark or a gluon
	  if( !( pid==1 || pid==2 || pid==3 || pid==4 || pid==5 || pid==6 || pid == 21 ) ) continue;

	  // require mother is not a top or W
	  if( mothid == 6 || mothid == 24) continue;

	  // found additional parton
	  npartons_ ++;
	  if( genps_p4().at(igen).pt() > maxpartonpt_ ) maxpartonpt_ = genps_p4().at(igen).pt();
	  //	  cout << "found parton, igen " << igen << " id " << pid << " motherid " << mothid << " pt " << genps_p4().at(igen).pt() << endl;

	}

	// if( npartons_ > 0 ){
	//   cout << endl << endl;
	//   dumpDocLines();
	//   cout << endl << endl;
	//   cout << "number of partons " << npartons_    << endl;
	//   cout << "max parton pt     " << maxpartonpt_ << endl;
	// }

	//count tops and only get two
	//ttbar_    = &(vttbar);
	if (ntops==2) {
	  LorentzVector ttpair = *t_ + *tbar_;
	  ttbar_    = &ttpair;
	  ptttbar_  = ttbar_->pt();
	  mttbar_   = ttbar_->mass();
	  etattbar_ = ttbar_->eta();
	}
	
	if( nels + nmus == 2) dilptgen = vdilepton.pt();
        
	if ( strcmp(prefix , "DYee"     ) == 0 &&  nels  != 2  ) continue;
	if ( strcmp(prefix , "DYmm"     ) == 0 &&  nmus  != 2  ) continue;
	if ( strcmp(prefix , "DYtautau" ) == 0 &&  ntaus != 2  ) continue;
	
	//splice together the DY samples - if its madgraph, then we do nothing
	if(TString(prefix).Contains("DY") && TString(evt_dataset().at(0)).Contains("madgraph") == false) {	
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

      //---------------------------
      // set event weight
      //---------------------------

      weight_ = -1.;

      if( TString(prefix).Contains("T2") ){
	mG_ = -999; //sparm_mG();
        mL_ = -999; //sparm_mL();
        x_  = -999; //sparm_mf();

        for (int i=0; i<(int)sparm_values().size(); ++i) {
	  if( TString(prefix).Contains("T2tt") ){
	    if (sparm_names().at(i).Contains("mstop")) mG_ = sparm_values().at(i);
	    if (sparm_names().at(i).Contains("mlsp")) mL_ = sparm_values().at(i);
	  }
	  if( TString(prefix).Contains("T2bw") ){
	    if (sparm_names().at(i).Contains("x")) mG_ = sparm_values().at(i);
	    if (sparm_names().at(i).Contains("mstop")) mL_ = sparm_values().at(i);
	    if (sparm_names().at(i).Contains("mlsp")) mL_ = sparm_values().at(i);
	  }

	}
        xsecsusy_  = mG_ > 0. ? stopPairCrossSection(mG_) : -999;
        weight_ = xsecsusy_ > 0. ? lumi * xsecsusy_ * (1000./50000.) : -999.;

	if( doTenPercent )	  weight_ *= 10;
      }

      else if(strcmp(prefix,"LMscan") == 0){ 

	m0_  = -999; //sparm_m0();
	m12_ = -999; //sparm_m12();

	ksusy_     = kfactorSUSY(m0_,m12_,"tanbeta10");
	ksusyup_   = kfactorSUSY(m0_,m12_,"tanbeta10Scale20");
	ksusydn_   = kfactorSUSY(m0_,m12_,"tanbeta10Scale05");
	xsecsusy_  = cmssm_loxsec(m0_,m12_);
	xsecsusy2_ = getMsugraCrossSection(m0_,m12_,10);

	weight_ = lumi * ksusy_ * xsecsusy_ * (1000. / 10000.); // k * xsec / nevents

	if( doTenPercent )	  weight_ *= 10;
      }

      else if( isData ){
	weight_ = 1;
      }

      else{

	weight_ = kFactor * evt_scale1fb() * lumi;
        //do a signed weight for mcatnlo
        if ( TString(prefix).Contains("mcatnlo") && genps_weight()<0) weight_ *= -1.;

      }

 
	
      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << endl;
  cout << "-----------------------" << endl;
  cout << "| Lepton yields       |" << endl;
  cout << "-----------------------" << endl;
  cout << "ne  " << netot       << endl;
  cout << "nm  " << nmtot       << endl;
  cout << "tot " << netot+nmtot << endl;

  cout << endl;
  cout << "-----------------------" << endl;
  cout << "| Preselection yields |" << endl;
  cout << "-----------------------" << endl;
  cout << "ne  " << nepass        << endl;
  cout << "nm  " << nmpass        << endl;
  cout << "tot " << nepass+nmpass << endl;
  cout << endl;

  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  stop_xsec_file->Close();

  //delete d_llsol; //REPLACETOPMASS

  return 0;

}


//--------------------------------------------------------------------
 
void genLooper::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  h_PU_trkpt = new TH1F(Form("%s_PU_trkpt",prefix),"track pt from PU interactions",100,0,100);
  h_dz_vtx_trk = new TH1F(Form("%s_dz_vtx_trk",prefix),"dZ between vtx and tracks",200,-0.1,0.1);
 
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

void genLooper::makeTree(char *prefix, bool doFakeApp, FREnum frmode ){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();


  //char* dir = "";
  //if     ( g_trig == e_lowpt  ) dir = "lowpt";
  //else if( g_trig == e_highpt ) dir = "highpt";

  //Super compressed ntuple here
  char* frsuffix = (char*) "";
  if( doFakeApp ){
    if ( frmode == e_qcd   ) frsuffix = (char*) "_doubleFake";
    if ( frmode == e_wjets ) frsuffix = (char*) "_singleFake";
  }

  char* tpsuffix = (char*) "";
  if( doTenPercent ) tpsuffix = (char*) "_tenPercent";

  outFile   = new TFile(Form("baby.root",prefix,frsuffix,tpsuffix), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix,frsuffix,tpsuffix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in genLooper.h
  outTree->Branch("acc_2010",        &acc_2010_,         "acc_2010/I");
  outTree->Branch("acc_highmet",     &acc_highmet_,      "acc_highmet/I");
  outTree->Branch("acc_highht",      &acc_highht_,       "acc_highht/I");

  outTree->Branch("eldup"     ,  &eldup_     ,  "eldup/I");  
  outTree->Branch("csc"       ,  &csc_       ,  "csc/I");  
  outTree->Branch("hbhe"      ,  &hbhe_      ,  "hbhe/I");  
  outTree->Branch("hbhenew"   ,  &hbhenew_   ,  "hbhenew/I");  
  outTree->Branch("hcallaser" ,  &hcallaser_ ,  "hcallaser/I");  
  outTree->Branch("ecaltp"    ,  &ecaltp_    ,  "ecaltp/I");  
  outTree->Branch("trkfail"   ,  &trkfail_   ,  "trkfail/I");  
  outTree->Branch("eebadsc"   ,  &eebadsc_   ,  "eebadsc/I");  

  outTree->Branch("isdata",          &isdata_,           "isdata/I");
  outTree->Branch("jetid",           &jetid_,            "jetid/I");
  outTree->Branch("jetid30",         &jetid30_,          "jetid30/I");
  outTree->Branch("json",            &json_,             "json/I");
  outTree->Branch("htoffset",        &htoffset_,         "htoffset/F");
  outTree->Branch("htuncor",         &htuncor_,          "htuncor/F");
  outTree->Branch("ptt",             &ptt_,              "ptt/F");
  outTree->Branch("pttbar",          &pttbar_,           "pttbar/F");
  outTree->Branch("ptttbar",         &ptttbar_,          "ptttbar/F");
  outTree->Branch("mttbar",          &mttbar_,           "mttbar/F");
  outTree->Branch("npartons",        &npartons_,         "npartons/I");
  outTree->Branch("nwzpartons",      &nwzpartons_,       "nwzpartons/I");
  outTree->Branch("hyptype",         &hyptype_,          "hyptype/I");
  outTree->Branch("maxpartonpt",     &maxpartonpt_,      "maxpartonpt/F");
  outTree->Branch("etattbar",        &etattbar_,         "etatbar/F");
  outTree->Branch("njetsoffset",     &njetsoffset_,      "njetsoffset/I");
  outTree->Branch("njetsuncor",      &njetsuncor_,       "njetsuncor/I");
  outTree->Branch("costhetaweight",  &costhetaweight_,   "costhetaweight/F");
  outTree->Branch("weight",          &weight_,           "weight/F");
  outTree->Branch("mutrigweight",    &mutrigweight_,     "mutrigweight/F");
  outTree->Branch("mutrigweight2",   &mutrigweight2_,    "mutrigweight2/F");
  outTree->Branch("sltrigweight",    &sltrigweight_,     "sltrigweight/F");
  outTree->Branch("dltrigweight",    &dltrigweight_,     "dltrigweight/F");
  outTree->Branch("trgeff",          &trgeff_,           "trgeff/F");
  outTree->Branch("pthat",           &pthat_,            "pthat/F");
  outTree->Branch("qscale",          &qscale_,           "qscale/F");
  outTree->Branch("mgcor",           &mgcor_,            "mgcor/F");
  outTree->Branch("wflav",           &wflav_,            "wflav/I");
  outTree->Branch("ksusy",           &ksusy_,            "ksusy/F");
  outTree->Branch("ksusyup",         &ksusyup_,          "ksusyup/F");
  outTree->Branch("ksusydn",         &ksusydn_,          "ksusydn/F");
  outTree->Branch("xsecsusy",        &xsecsusy_,         "xsecsusy/F");
  outTree->Branch("xsecsusy2",       &xsecsusy2_,        "xsecsusy2/F");
  outTree->Branch("smeff",           &smeff_,            "smeff/F");
  outTree->Branch("k",               &k_,                "k/F");
  outTree->Branch("mllgen",          &mllgen_,           "mllgen/F");
  outTree->Branch("ptwgen",          &ptwgen_,           "ptwgen/F");
  outTree->Branch("ptzgen",          &ptzgen_,           "ptzgen/F");
  outTree->Branch("nlep",            &nlep_,             "nlep/I");
  outTree->Branch("nosel",           &nosel_,            "nosel/I");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("ngoodmu",         &ngoodmu_,          "ngoodmu/I");
  outTree->Branch("mull",            &mull_,             "mull/I");
  outTree->Branch("mult",            &mult_,             "mult/I");
  //outTree->Branch("eltrijet",        &eltrijet_,         "eltrijet/I");
  //outTree->Branch("mutrijet",        &mutrijet_,         "mutrijet/I");
  //outTree->Branch("ldi",             &ldi_,              "ldi/I");
  //outTree->Branch("ltri",            &ltri_,             "ltri/I");
  //outTree->Branch("smu",             &smu_,              "smu/I");
  //outTree->Branch("smu30",           &smu30_,            "smu30/I");
  //outTree->Branch("trgmu30",         &trgmu30_,          "trgmu30/I");
  //outTree->Branch("trg2mu30",        &trg2mu30_,         "trg2mu30/I");
  //outTree->Branch("dil",             &dil_,              "dil/I");
  outTree->Branch("mullgen",         &mullgen_,          "mullgen/I");
  outTree->Branch("multgen",         &multgen_,          "multgen/I");
  outTree->Branch("proc",            &proc_,             "proc/I");
  outTree->Branch("leptype",         &leptype_,          "leptype/I");
  outTree->Branch("topmass",         &topmass_,          "topmass/F");
  outTree->Branch("dilmass",         &dilmass_,          "dilmass/F");
  outTree->Branch("dilrecoil",       &dilrecoil_,        "dilrecoil/F");
  outTree->Branch("dilrecoilparl",   &dilrecoilparl_,    "dilrecoilparl/F");
  outTree->Branch("dilrecoilperp",   &dilrecoilperp_,    "dilrecoilperp/F");
  outTree->Branch("tcmet",           &tcmet_,            "tcmet/F");
  outTree->Branch("genmet",          &genmet_,           "genmet/F");
  outTree->Branch("gensumet",        &gensumet_,         "gensumet/F");
  outTree->Branch("genmetphi",       &genmetphi_,        "genmetphi/F");
  outTree->Branch("trkmet",          &trkmet_,           "trkmet/F");
  outTree->Branch("trkmetphi",       &trkmetphi_,        "trkmetphi/F");
  outTree->Branch("trkmet_nolepcorr",    &trkmet_nolepcorr_,    "trkmet_nolepcorr/F");
  outTree->Branch("trkmetphi_nolepcorr", &trkmetphi_nolepcorr_, "trkmetphi_nolepcorr/F");
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
  outTree->Branch("ngenjets",        &ngenjets_,         "ngenjets/I");
  outTree->Branch("njpt",            &njpt_,             "njpt/I");

  outTree->Branch("trgmu1"         ,  &trgmu1_          ,    "trgmu1/I"      );
  outTree->Branch("trgmu2"         ,  &trgmu2_          ,    "trgmu2/I"      );
  outTree->Branch("trgel1"         ,  &trgel1_          ,    "trgel1/I"      );
  outTree->Branch("trgel2"         ,  &trgel2_          ,    "trgel2/I"      );

  outTree->Branch("isomu24"        ,  &isomu24_         ,    "isomu24/I"     );
  outTree->Branch("ele27wp80"      ,  &ele27wp80_       ,    "ele27wp80/I"   );
  outTree->Branch("mm"             ,  &mm_              ,    "mm/I"          );
  outTree->Branch("mmtk"           ,  &mmtk_            ,    "mmtk/I"        );
  outTree->Branch("me"             ,  &me_              ,    "me/I"          );
  outTree->Branch("em"             ,  &em_              ,    "em/I"          );
  outTree->Branch("mu"             ,  &mu_              ,    "mu/I"          );
  outTree->Branch("ee"             ,  &ee_              ,    "ee/I"          );

  // pfjets L1FastL2L3Res
  outTree->Branch("npfjets30",        &npfjets30_,        "npfjets30/I");
  outTree->Branch("npfjets35",        &npfjets35_,        "npfjets35/I");
  outTree->Branch("npfjets40",        &npfjets40_,        "npfjets40/I");
  outTree->Branch("npfjets45",        &npfjets45_,        "npfjets45/I");
  outTree->Branch("npfjets30lepcorr", &npfjets30lepcorr_, "npfjets30lepcorr/I");
  outTree->Branch("knjets",           &knjets_,           "knjets/F");

  //rho correction
  outTree->Branch("rhovor",          &rhovor_,           "rhovor/F");

  outTree->Branch("htpf30",          &htpf30_,           "htpf30/F");
  outTree->Branch("htpf35",          &htpf35_,           "htpf35/F");
  outTree->Branch("htpf40",          &htpf40_,           "htpf40/F");
  outTree->Branch("htpf45",          &htpf45_,           "htpf45/F");

  // type1 met flavors
  outTree->Branch("t1met10",         &t1met10_,          "t1met10/F");
  outTree->Branch("t1met20",         &t1met20_,          "t1met20/F");
  outTree->Branch("t1met30",         &t1met30_,          "t1met30/F");
  outTree->Branch("t1met10phi",      &t1met10phi_,       "t1met10phi/F");
  outTree->Branch("t1met20phi",      &t1met20phi_,       "t1met20phi/F");
  outTree->Branch("t1met30phi",      &t1met30phi_,       "t1met30phi/F");
  outTree->Branch("t1met10mt",       &t1met10mt_,        "t1met10mt/F");
  outTree->Branch("t1met20mt",       &t1met20mt_,        "t1met20mt/F");
  outTree->Branch("t1met30mt",       &t1met30mt_,        "t1met30mt/F");
  outTree->Branch("lepmetpt",        &lepmetpt_,         "lepmetpt/F");
  outTree->Branch("lept1met10pt",    &lept1met10pt_,     "lept1met10pt/F");

  //met variables with phi correction
  outTree->Branch("t1metphicorr"       , &t1metphicorr_       , "t1metphicorr/F");
  outTree->Branch("t1metphicorrup"     , &t1metphicorrup_     , "t1metphicorrup/F");
  outTree->Branch("t1metphicorrdn"     , &t1metphicorrdn_     , "t1metphicorrdn/F");
  outTree->Branch("t1metphicorrphi"    , &t1metphicorrphi_    , "t1metphicorrphi/F");
  outTree->Branch("t1metphicorrphiup"  , &t1metphicorrphiup_  , "t1metphicorrphiup/F");
  outTree->Branch("t1metphicorrphidn"  , &t1metphicorrphidn_  , "t1metphicorrphidn/F");
  outTree->Branch("t1metphicorrlep"    , &t1metphicorrlep_    , "t1metphicorrlep/F");
  outTree->Branch("t1metphicorrlepphi" , &t1metphicorrlepphi_ , "t1metphicorrlepphi/F");
  outTree->Branch("t1metphicorrmt"     , &t1metphicorrmt_     , "t1metphicorrmt/F");
  outTree->Branch("t1metphicorrmtup"   , &t1metphicorrmtup_   , "t1metphicorrmtup/F");
  outTree->Branch("t1metphicorrmtdn"   , &t1metphicorrmtdn_   , "t1metphicorrmtdn/F");
  outTree->Branch("t1metphicorrlepmt"  , &t1metphicorrlepmt_  , "t1metphicorrlepmt/F");

  outTree->Branch("htpfres30",        &htpfres30_,        "htpfres30/F");
  outTree->Branch("htpfres35",        &htpfres35_,        "htpfres35/F");
  outTree->Branch("htpfres40",        &htpfres40_,        "htpfres40/F");
  outTree->Branch("htpfres45",        &htpfres45_,        "htpfres45/F");
				      
  // btag variables		      
  outTree->Branch("nbtagsssv",        &nbtagsssv_,        "nbtagsssv/I");
  outTree->Branch("nbtagstcl",        &nbtagstcl_,        "nbtagstcl/I");
  outTree->Branch("nbtagstcm",        &nbtagstcm_,        "nbtagstcm/I");
  outTree->Branch("nbtagscsvl",       &nbtagscsvl_,       "nbtagscsvl/I");
  outTree->Branch("nbtagscsvm",       &nbtagscsvm_,       "nbtagscsvm/I");
  outTree->Branch("nbtagscsvt",       &nbtagscsvt_,       "nbtagscsvt/I");
  outTree->Branch("nbtagsssvcorr",    &nbtagsssvcorr_,    "nbtagsssvcorr/I");
  outTree->Branch("nbtagstclcorr",    &nbtagstclcorr_,    "nbtagstclcorr/I");
  outTree->Branch("nbtagstcmcorr",    &nbtagstcmcorr_,    "nbtagstcmcorr/I");
  outTree->Branch("nbtagscsvlcorr",   &nbtagscsvlcorr_,   "nbtagscsvlcorr/I");
  outTree->Branch("nbtagscsvmcorr",   &nbtagscsvmcorr_,   "nbtagscsvmcorr/I");
  outTree->Branch("nbtagscsvtcott",   &nbtagscsvtcorr_,   "nbtagscsvtcorr/I");
  outTree->Branch("bjet1",            &bjet1_,            "bjet1/I");
  outTree->Branch("bjet2",            &bjet2_,            "bjet2/I");
  outTree->Branch("bjet3",            &bjet3_,            "bjet3/I");
  outTree->Branch("bjet4",            &bjet4_,            "bjet4/I");
  outTree->Branch("bjet5",            &bjet5_,            "bjet5/I");
  outTree->Branch("bjet6",            &bjet6_,            "bjet6/I");
  outTree->Branch("lepjet1",          &lepjet1_,          "lepjet1/I");
  outTree->Branch("lepjet2",          &lepjet2_,          "lepjet2/I");
  outTree->Branch("lepjet3",          &lepjet3_,          "lepjet3/I");
  outTree->Branch("lepjet4",          &lepjet4_,          "lepjet4/I");
  outTree->Branch("lepjet5",          &lepjet5_,          "lepjet5/I");
  outTree->Branch("lepjet6",          &lepjet6_,          "lepjet6/I");
  outTree->Branch("qgjet1",           &qgjet1_,           "qgjet1/I");
  outTree->Branch("qgjet2",           &qgjet2_,           "qgjet2/I");
  outTree->Branch("qgjet3",           &qgjet3_,           "qgjet3/I");
  outTree->Branch("qgjet4",           &qgjet4_,           "qgjet4/I");
  outTree->Branch("qgjet5",           &qgjet5_,           "qgjet5/I");
  outTree->Branch("qgjet6",           &qgjet6_,           "qgjet6/I");
  outTree->Branch("genjetdr1",        &genjetdr1_,        "genjetdr1/F");
  outTree->Branch("genjetdr2",        &genjetdr2_,        "genjetdr2/F");
  outTree->Branch("genjetdr3",        &genjetdr3_,        "genjetdr3/F");
  outTree->Branch("genjetdr4",        &genjetdr4_,        "genjetdr4/F");
  outTree->Branch("genjetdr5",        &genjetdr5_,        "genjetdr5/F");
  outTree->Branch("genjetdr6",        &genjetdr6_,        "genjetdr6/F");
  outTree->Branch("njetsUp",          &njetsUp_,          "njetsUp/I");
  outTree->Branch("njetsDown",        &njetsDown_,        "njetsDown/I");
  outTree->Branch("htUp",             &htUp_,             "htUp/F");
  outTree->Branch("htDown",           &htDown_,           "htDown/F");
  outTree->Branch("npu",              &npu_,              "npu/I");
  outTree->Branch("npuMinusOne",      &npuMinusOne_,      "npuMinusOne/I");
  outTree->Branch("npuPlusOne",       &npuPlusOne_,       "npuPlusOne/I");
  outTree->Branch("nvtx",             &nvtx_,             "nvtx/I");
  outTree->Branch("nvtxweight",       &nvtxweight_,       "nvtxweight/F");
  outTree->Branch("n3dvtxweight",     &n3dvtxweight_,     "n3dvtxweight/F");
  outTree->Branch("pdfid1",           &pdfid1_,           "pdfid1/I");
  outTree->Branch("pdfid2",           &pdfid2_,           "pdfid2/I");
  outTree->Branch("pdfx1",            &pdfx1_,            "pdfx1/F");
  outTree->Branch("pdfx2",            &pdfx2_,            "pdfx2/F");
  outTree->Branch("pdfQ",             &pdfQ_,             "pdfQ/F");
  outTree->Branch("vecjetpt",         &vecjetpt_,         "vecjetpt/F");
  outTree->Branch("pass",             &pass_,             "pass/I");
  outTree->Branch("passz",            &passz_,            "passz/I");
  outTree->Branch("m0",               &m0_,               "m0/F");
  outTree->Branch("mg",               &mG_,               "mg/F");
  outTree->Branch("ml",               &mL_,               "ml/F");
  outTree->Branch("x",                &x_,                "x/F");
  outTree->Branch("m12",              &m12_,              "m12/F");
  outTree->Branch("lep1chi2ndf",      &lep1chi2ndf_,      "lep1chi2ndf/F");
  outTree->Branch("lep2chi2ndf",      &lep2chi2ndf_,      "lep2chi2ndf/F");
  outTree->Branch("lep1dpt",          &lep1dpt_,          "lep1dpt/F");
  outTree->Branch("lep2dpt",          &lep2dpt_,          "lep2dpt/F");
  outTree->Branch("id1",              &id1_,              "id1/I");
  outTree->Branch("id2",              &id2_,              "id2/I");
  outTree->Branch("leptype1",         &leptype1_,         "leptype1/I");
  outTree->Branch("leptype2",         &leptype2_,         "leptype2/I");
  outTree->Branch("w1",               &w1_,               "w1/I");
  outTree->Branch("w2",               &w2_,               "w2/I");
  outTree->Branch("iso1",             &iso1_,             "iso1/F");
  outTree->Branch("isont1",           &isont1_,           "isont1/F");
  outTree->Branch("isopfold1",    &isopfold1_,     "isopfold1/F");
  outTree->Branch("isopf1",           &isopf1_,           "isopf1/F");
  outTree->Branch("etasc1",           &etasc1_,           "etasc1/F");
  outTree->Branch("etasc2",           &etasc2_,           "etasc2/F");
  outTree->Branch("eoverpin",         &eoverpin_,         "eoverpin/F");
  outTree->Branch("eoverpout",        &eoverpout_,        "eoverpout/F");
  outTree->Branch("dEtaIn", &dEtaIn_, "dEtaIn/F");
  outTree->Branch("dPhiIn", &dPhiIn_, "dPhiIn/F");
  outTree->Branch("sigmaIEtaIEta", &sigmaIEtaIEta_, "sigmaIEtaIEta/F");
  outTree->Branch("hOverE", &hOverE_, "hOverE/F");
  outTree->Branch("ooemoop", &ooemoop_, "ooemoop/F");
  outTree->Branch("d0vtx", &d0vtx_, "d0vtx/F");
  outTree->Branch("dzvtx", &dzvtx_, "dzvtx/F");
  outTree->Branch("expinnerlayers", &expinnerlayers_, "expinnerlayers/F");
  outTree->Branch("fbrem", &fbrem_, "fbrem/F");
  outTree->Branch("pfisoch", &pfisoch_, "pfisoch/F");
  outTree->Branch("pfisoem", &pfisoem_, "pfisoem/F");
  outTree->Branch("pfisonh", &pfisonh_, "pfisonh/F");
  outTree->Branch("eSC", & eSC_, "eSC/F");
  outTree->Branch("phiSC", & phiSC_, "phiSC/F");
  outTree->Branch("eSCRaw", & eSCRaw_, "eSCRaw/F");
  outTree->Branch("eSCPresh", & eSCPresh_, "eSCPresh/F");
  
  outTree->Branch("eoverpin2",         &eoverpin2_,         "eoverpin2/F");
  outTree->Branch("eoverpout2",        &eoverpout2_,        "eoverpout2/F");
  outTree->Branch("dEtaIn2", &dEtaIn2_, "dEtaIn2/F");
  outTree->Branch("dPhiIn2", &dPhiIn2_, "dPhiIn2/F");
  outTree->Branch("sigmaIEtaIEta2", &sigmaIEtaIEta2_, "sigmaIEtaIEta2/F");
  outTree->Branch("hOverE2", &hOverE2_, "hOverE2/F");
  outTree->Branch("ooemoop2", &ooemoop2_, "ooemoop2/F");
  outTree->Branch("d0vtx2", &d0vtx2_, "d0vtx2/F");
  outTree->Branch("dzvtx2", &dzvtx2_, "dzvtx2/F");
  outTree->Branch("expinnerlayers2", &expinnerlayers2_, "expinnerlayers2/F");
  outTree->Branch("fbrem2", &fbrem2_, "fbrem2/F");
  outTree->Branch("pfisoch2", &pfisoch2_, "pfisoch2/F");
  outTree->Branch("pfisoem2", &pfisoem2_, "pfisoem2/F");
  outTree->Branch("pfisonh2", &pfisonh2_, "pfisonh2/F");
  outTree->Branch("eSC2", & eSC2_, "eSC2/F");
  outTree->Branch("phiSC2", & phiSC2_, "phiSC2/F");
  outTree->Branch("eSCRaw2", & eSCRaw2_, "eSCRaw2/F");
  outTree->Branch("eSCPresh2", & eSCPresh2_, "eSCPresh2/F");
  
  outTree->Branch("iso2",             &iso2_,             "iso2/F");
  outTree->Branch("ecalveto1",        &ecalveto1_,        "ecalveto1/F");
  outTree->Branch("ecalveto2",        &ecalveto2_,        "ecalveto2/F");
  outTree->Branch("hcalveto1",        &hcalveto1_,        "hcalveto1/F");
  outTree->Branch("hcalveto2",        &hcalveto2_,        "hcalveto2/F");
  outTree->Branch("isont2",           &isont2_,           "isont2/F");
  outTree->Branch("isopf2",           &isopf2_,           "isopf2/F");
  outTree->Branch("ptl1",             &ptl1_,             "ptl1/F");
  outTree->Branch("ptl2",             &ptl2_,             "ptl2/F");
  outTree->Branch("etal1",            &etal1_,            "etal1/F");
  outTree->Branch("etal2",            &etal2_,            "etal2/F");
  outTree->Branch("phil1",            &phil1_,            "phil1/F");
  outTree->Branch("phil2",            &phil2_,            "phil2/F");
  outTree->Branch("meff",             &meff_,             "meff/F");
  outTree->Branch("mt",               &mt_,               "mt/F");
  outTree->Branch("dataset",          &dataset_,          "dataset[200]/C");
  outTree->Branch("run",              &run_,              "run/I");
  outTree->Branch("lumi",             &lumi_,             "lumi/I");
  outTree->Branch("event",            &event_,            "event/I");
  outTree->Branch("y",                &y_,                "y/F");  
  outTree->Branch("ht",               &ht_,               "ht/F");  
  outTree->Branch("htgen",            &htgen_,            "htgen/F");  
  outTree->Branch("htjpt",            &htjpt_,            "htjpt/F");  
  outTree->Branch("nels",             &nels_,             "nels/I");  
  outTree->Branch("nmus",             &nmus_,             "nmus/I");  
  outTree->Branch("ntaus",            &ntaus_,            "ntaus/I");  
  outTree->Branch("nleps",            &nleps_,            "nleps/I");  
  outTree->Branch("dphijm",           &dphijm_,           "dphijm/F");  
  outTree->Branch("ptjetraw",         &ptjetraw_,         "ptjetraw/F");  
  outTree->Branch("ptjet23",          &ptjet23_,          "ptjet23/F");  
  outTree->Branch("ptjetF23",         &ptjetF23_,         "ptjetF23/F");  
  outTree->Branch("ptjetO23",         &ptjetO23_,         "ptjetO23/F");  
  //outTree->Branch("cosphijz",         &cosphijz_,         "cosphijz/F");  
  outTree->Branch("mcid1",            &mcid1_,            "mcid1/I");  
  outTree->Branch("mcdr1",            &mcdr1_,            "mcdr1/F");  
  outTree->Branch("mcdecay1",         &mcdecay1_,         "mcdecay1/I");  
  outTree->Branch("mcndec1",          &mcndec1_,          "mcndec1/I");  
  outTree->Branch("mcndec2",          &mcndec2_,          "mcndec2/I");  
  outTree->Branch("mcndeckls1",       &mcndeckls1_,       "mcndeckls1/I");  
  outTree->Branch("mcndeckls2",       &mcndeckls2_,       "mcndeckls2/I");  
  outTree->Branch("mcndecem1",        &mcndecem1_,        "mcndecem1/I");  
  outTree->Branch("mcndecem2",        &mcndecem2_,        "mcndecem2/I");  
  outTree->Branch("mcid2",            &mcid2_,            "mcid2/I");  
  outTree->Branch("mcdr2",            &mcdr2_,            "mcdr2/F");  
  outTree->Branch("mcdecay2",         &mcdecay2_,         "mcdecay2/I");  
  outTree->Branch("mctaudpt1",        &mctaudpt1_,        "mctaudpt1/F");  
  outTree->Branch("mctaudpt2",        &mctaudpt2_,        "mctaudpt2/F");  
  outTree->Branch("mctaudid1",        &mctaudid1_,        "mctaudid1/I");  
  outTree->Branch("mctaudid2",        &mctaudid2_,        "mctaudid2/I");  
  outTree->Branch("mlepid",           &mlepid_,           "mlepid/I");  
  outTree->Branch("mleppassid",       &mleppassid_,       "mleppassid/I");  
  outTree->Branch("mleppassiso",      &mleppassiso_,      "mleppassiso/I");  
  outTree->Branch("mlepiso",          &mlepiso_,          "mlepiso/F");  
  outTree->Branch("mlepdr",           &mlepdr_,           "mlepdr/F");  
  outTree->Branch("pflepiso",         &pflepiso_,         "pflepiso/F");  
  outTree->Branch("pflepdr",          &pflepdr_,          "pflepdr/F");  
  outTree->Branch("pfleppt",          &pfleppt_,          "pfleppt/F");  
  outTree->Branch("pflepmindrj",      &pflepmindrj_,      "pflepmindrj/F");  
  outTree->Branch("pftaudiso",        &pftaudiso_,        "pftaudiso/F");  
  outTree->Branch("pftauddr",         &pftauddr_,         "pftauddr/F");  
  outTree->Branch("pftaudpt",         &pftaudpt_,         "pftaudpt/F");  
  outTree->Branch("pftaudmindrj",     &pftaudmindrj_,     "pftaudmindrj/F");  
  outTree->Branch("pfcandiso5",       &pfcandiso5_,       "pfcandiso5/F");  
  outTree->Branch("pfcandpt5",        &pfcandpt5_,        "pfcandpt5/F");  
  outTree->Branch("pfcandmindrj5",    &pfcandmindrj5_,    "pfcandmindrj5/F");  
  outTree->Branch("pfcandiso10",      &pfcandiso10_,      "pfcandiso10/F");  
  outTree->Branch("pfcandpt10",       &pfcandpt10_,       "pfcandpt10/F");  
  outTree->Branch("pfcandmindrj10",   &pfcandmindrj10_,   "pfcandmindrj10/F");  
  outTree->Branch("emjet10",          &emjet10_,          "emjet10/F");  
  outTree->Branch("mjj",              &mjj_,              "mjj/F");  
  outTree->Branch("emjet20",          &emjet20_,          "emjet20/F");  
  outTree->Branch("trkpt5",           &trkpt5_,           "trkpt5/F");  
  outTree->Branch("trkpt10",          &trkpt10_,          "trkpt10/F");  
  outTree->Branch("mleptrk5",         &mleptrk5_,         "mleptrk5/F");  
  outTree->Branch("mleptrk10",        &mleptrk10_,        "mleptrk10/F");  
  outTree->Branch("trkreliso5",       &trkreliso5_,       "trkreliso5/F");  
  outTree->Branch("trkreliso10",      &trkreliso10_,      "trkreliso10/F");  
  outTree->Branch("trkpt5loose",      &trkpt5loose_,      "trkpt5loose/F");  
  outTree->Branch("trkpt10loose",     &trkpt10loose_,     "trkpt10loose/F");  
  outTree->Branch("trkreliso5loose",  &trkreliso5loose_,  "trkreliso5loose/F");  
  outTree->Branch("trkreliso10loose", &trkreliso10loose_, "trkreliso10loose/F");  

  outTree->Branch("trkpt10pt0p1",     &trkpt10pt0p1_,      "trkpt10pt0p1/F");  
  outTree->Branch("trkpt10pt0p2",     &trkpt10pt0p2_,      "trkpt10pt0p2/F");  
  outTree->Branch("trkpt10pt0p3",     &trkpt10pt0p3_,      "trkpt10pt0p3/F");  
  outTree->Branch("trkpt10pt0p4",     &trkpt10pt0p4_,      "trkpt10pt0p4/F");  
  outTree->Branch("trkpt10pt0p5",     &trkpt10pt0p5_,      "trkpt10pt0p5/F");  
  outTree->Branch("trkpt10pt0p6",     &trkpt10pt0p6_,      "trkpt10pt0p6/F");  
  outTree->Branch("trkpt10pt0p7",     &trkpt10pt0p7_,      "trkpt10pt0p7/F");  
  outTree->Branch("trkpt10pt0p8",     &trkpt10pt0p8_,      "trkpt10pt0p8/F");  
  outTree->Branch("trkpt10pt0p9",     &trkpt10pt0p9_,      "trkpt10pt0p9/F");  
  outTree->Branch("trkpt10pt1p0",     &trkpt10pt1p0_,      "trkpt10pt1p0/F");  
  outTree->Branch("trkreliso10pt0p1", &trkreliso10pt0p1_,  "trkreliso10pt0p1/F");  
  outTree->Branch("trkreliso10pt0p2", &trkreliso10pt0p2_,  "trkreliso10pt0p2/F");  
  outTree->Branch("trkreliso10pt0p3", &trkreliso10pt0p3_,  "trkreliso10pt0p3/F");  
  outTree->Branch("trkreliso10pt0p4", &trkreliso10pt0p4_,  "trkreliso10pt0p4/F");  
  outTree->Branch("trkreliso10pt0p5", &trkreliso10pt0p5_,  "trkreliso10pt0p5/F");  
  outTree->Branch("trkreliso10pt0p6", &trkreliso10pt0p6_,  "trkreliso10pt0p6/F");  
  outTree->Branch("trkreliso10pt0p7", &trkreliso10pt0p7_,  "trkreliso10pt0p7/F");  
  outTree->Branch("trkreliso10pt0p8", &trkreliso10pt0p8_,  "trkreliso10pt0p8/F");  
  outTree->Branch("trkreliso10pt0p9", &trkreliso10pt0p9_,  "trkreliso10pt0p9/F");  
  outTree->Branch("trkreliso10pt1p0", &trkreliso10pt1p0_,  "trkreliso10pt1p0/F");  

  outTree->Branch("pfcandpt10pt0p1",  &pfcandpt10pt0p1_,   "pfcandpt10pt0p1/F");  
  outTree->Branch("pfcandpt10pt0p2",  &pfcandpt10pt0p2_,   "pfcandpt10pt0p2/F");  
  outTree->Branch("pfcandpt10pt0p3",  &pfcandpt10pt0p3_,   "pfcandpt10pt0p3/F");  
  outTree->Branch("pfcandpt10pt0p4",  &pfcandpt10pt0p4_,   "pfcandpt10pt0p4/F");  
  outTree->Branch("pfcandpt10pt0p5",  &pfcandpt10pt0p5_,   "pfcandpt10pt0p5/F");  
  outTree->Branch("pfcandpt10pt0p6",  &pfcandpt10pt0p6_,   "pfcandpt10pt0p6/F");  
  outTree->Branch("pfcandpt10pt0p7",  &pfcandpt10pt0p7_,   "pfcandpt10pt0p7/F");  
  outTree->Branch("pfcandpt10pt0p8",  &pfcandpt10pt0p8_,   "pfcandpt10pt0p8/F");  
  outTree->Branch("pfcandpt10pt0p9",  &pfcandpt10pt0p9_,   "pfcandpt10pt0p9/F");  
  outTree->Branch("pfcandpt10pt1p0",  &pfcandpt10pt1p0_,   "pfcandpt10pt1p0/F");  
  outTree->Branch("pfcandiso10pt0p1", &pfcandiso10pt0p1_,  "pfcandiso10pt0p1/F");  
  outTree->Branch("pfcandiso10pt0p2", &pfcandiso10pt0p2_,  "pfcandiso10pt0p2/F");  
  outTree->Branch("pfcandiso10pt0p3", &pfcandiso10pt0p3_,  "pfcandiso10pt0p3/F");  
  outTree->Branch("pfcandiso10pt0p4", &pfcandiso10pt0p4_,  "pfcandiso10pt0p4/F");  
  outTree->Branch("pfcandiso10pt0p5", &pfcandiso10pt0p5_,  "pfcandiso10pt0p5/F");  
  outTree->Branch("pfcandiso10pt0p6", &pfcandiso10pt0p6_,  "pfcandiso10pt0p6/F");  
  outTree->Branch("pfcandiso10pt0p7", &pfcandiso10pt0p7_,  "pfcandiso10pt0p7/F");  
  outTree->Branch("pfcandiso10pt0p8", &pfcandiso10pt0p8_,  "pfcandiso10pt0p8/F");  
  outTree->Branch("pfcandiso10pt0p9", &pfcandiso10pt0p9_,  "pfcandiso10pt0p9/F");  
  outTree->Branch("pfcandiso10pt1p0", &pfcandiso10pt1p0_,  "pfcandiso10pt1p0/F");  

  outTree->Branch("mbb",             &mbb_,              "mbb/F");
  outTree->Branch("lep1pfjetdr",     &lep1pfjetdr_,      "lep1pfjetdr/F");  
  outTree->Branch("lep2pfjetdr",     &lep2pfjetdr_,      "lep2pfjetdr/F");  

  outTree->Branch("mclep"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep_      );
  outTree->Branch("mcnu"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mcnu_       );
  outTree->Branch("mcmln",           &mcmln_,              "mcmln/F");
  outTree->Branch("mcmtln",          &mcmtln_,             "mcmtln/F");

  outTree->Branch("mlep"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mlep_	);
  outTree->Branch("lep1"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_	);
  outTree->Branch("lep2"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_	);
  outTree->Branch("trklep1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &trklep1_	);
  outTree->Branch("trklep2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &trklep2_	);
  outTree->Branch("gfitlep1"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &gfitlep1_	);
  outTree->Branch("gfitlep2"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &gfitlep2_	);
  outTree->Branch("lepp"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lepp_	);
  outTree->Branch("lepm"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lepm_	);
  outTree->Branch("pflep1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pflep1_	);
  outTree->Branch("pflep2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pflep2_	);
  outTree->Branch("leppfjet1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &leppfjet1_	);
  outTree->Branch("leppfjet2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &leppfjet2_	);
  outTree->Branch("mclep1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep1_	);
  outTree->Branch("mclep2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep2_	);
  outTree->Branch("mctaud1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaud1_	);
  outTree->Branch("mctaud2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaud2_	);
  outTree->Branch("mctaudvis1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaudvis1_	);
  outTree->Branch("mctaudvis2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaudvis2_	);
  outTree->Branch("pflep"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pflep_	);
  outTree->Branch("pftaud"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pftaud_	);
  outTree->Branch("pfcand5"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcand5_	);
  outTree->Branch("pfcand10"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfcand10_	);
  outTree->Branch("jet"	      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_	);

  outTree->Branch("pfjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet1_	);
  outTree->Branch("pfjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet2_	);
  outTree->Branch("pfjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet3_	);
  outTree->Branch("pfjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet4_	);
  outTree->Branch("pfjet5"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet5_	);
  outTree->Branch("pfjet6"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet6_	);

  outTree->Branch("nonisoel"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisoel_	);
  outTree->Branch("nonisomu"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisomu_	);
  outTree->Branch("t"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &t_   	);
  outTree->Branch("tbar"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &tbar_   	);
  outTree->Branch("ttbar"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &ttbar_   	);
  outTree->Branch("mcStop1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mcStop1_   );
  outTree->Branch("mcStop2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mcStop2_   );
  outTree->Branch("neutralino_t"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &neutralino_t_    );
  outTree->Branch("neutralino_tbar" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &neutralino_tbar_ );
  outTree->Branch("W_t"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &W_t_    );
  outTree->Branch("W_tbar" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &W_tbar_ );

  outTree->Branch("b"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &b_   	);
  outTree->Branch("bbar"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &bbar_   	);

}

//--------------------------------------------------------------------

float genLooper::dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

float genLooper::trackIso( int thisPf , float coneR , float dz_thresh , bool dovtxcut , float pt_thresh ){

  float iso = 0.0;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate
    if( cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals
    if( cms2.pfcands_p4().at(ipf).pt() < pt_thresh ) continue; // skip pfcands below pt threshold

    if( dRbetweenVectors( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) ) > coneR ) continue;

    int itrk = cms2.pfcands_trkidx().at(ipf);
    
    if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ){
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
	fillOverFlow( h_dz_vtx_trk , mydz );
	
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

    iso += cms2.pfcands_p4().at(ipf).pt();

  }

  return iso;
}

// std::vector<float> trackIsoPtRanges( int thisPf , float coneR , float dz_thresh ){

//   float iso[10];
//   for (int i=0; i<10; ++i) iso[i] = 0.0;
  
//   for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

//     if( ipf == thisPf                 ) continue; // skip this PFCandidate
//     if( cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals

//     if( dRbetweenVectors( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) ) > coneR ) continue;

//     int itrk = cms2.pfcands_trkidx().at(ipf);
    
//     if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ){
//       //note: this should only happen for electrons which do not have a matched track
//       //currently we are just ignoring these guys
//       continue;
//     }
    
//     float mindz = dz_trk_vtx(itrk,0);
//     if ( fabs(mindz) > dz_thresh )     continue;

//     //---------------------------------------
//     // passes cuts, add up isolation value
//     //---------------------------------------
    
//     //figure out which pT threshold this track passes
//     float pfcandpt = cms2.pfcands_p4().at(ipf).pt();
//     for int (i=0; i<10; ++i) {
// 	if ( pfcandpt > 0.1*i+0.1 ) 
// 	  iso[i] += pfcandpt;	  
//       }
//   } 

//   std::vector<float> isovec;
//   for (int i=0; i<10; ++i) isovec.push_back(iso[i]);
  
//   return isovec;

// }

std::vector<float> genLooper::totalIso( int thisPf , float coneR , float dz_thresh ){

  float iso = 0.; 

  float emiso = 0.;
  float nhiso = 0.;
  float chiso = 0.;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate
    float dR = dRbetweenVectors( pfcands_p4().at(ipf) , pfcands_p4().at(thisPf) );
    if( dR > coneR ) continue;
    float pfpt = cms2.pfcands_p4().at(ipf).pt();

    //----------------------------------------
    // neutrals
    //----------------------------------------

    if( cms2.pfcands_charge().at(ipf) == 0 ) {
      // skip neutrals with pT < 1 GeV to reduce pileup dependence
      if ( pfpt < 1. ) continue; 
      int pfid = abs(cms2.pfcands_particleId().at(ipf));
        // get isolation parameters
        
        float pfeta = cms2.pfcands_p4().at(ipf).eta();
        float deta = fabs(pfeta - cms2.pfcands_p4().at(thisPf).eta());
	//photons 
	if (pfid == 22) {
	  // to remove pi0s and possible radiation from their photons
	  if (deta <= 0.1) continue;
	  iso += pfpt;
	  emiso += pfpt;
	}
	else {
	  iso += pfpt;
	  nhiso += pfpt;
	}
    }

    //----------------------------------------
    // charged candidates
    //----------------------------------------

    else {
      
      int itrk = cms2.pfcands_trkidx().at(ipf);
    
      if( itrk >= (int)trks_trk_p4().size() || itrk < 0 ){
	//note: this should only happen for electrons which do not have a matched track
	//currently we are just ignoring these guys
	continue; 
      }
    
      //----------------------------------------
      // find closest PV and dz w.r.t. that PV
      //----------------------------------------
    
      // float mindz = 999.;
      // int vtxi    = -1;
      
      // for (unsigned int ivtx = 0; ivtx < cms2.davtxs_position().size(); ivtx++) {
	
      // 	if(!isGoodDAVertex(ivtx)) continue;

      // 	float mydz = dz_trk_vtx(itrk,ivtx);
      
      // 	if (fabs(mydz) < fabs(mindz)) {
      // 	  mindz = mydz;
      // 	  vtxi = ivtx;
      // 	}
         
      // }
    
      //----------------------------------------------------------------------------
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
      //----------------------------------------------------------------------------
    
      // if ( vtxi != 0               )     continue;
      float mindz = dz_trk_vtx(itrk,0);
      if ( fabs(mindz) > dz_thresh )     continue;

      //---------------------------------------
      // passes cuts, add up isolation value
      //---------------------------------------

      iso += pfpt;
      chiso += pfpt;
    }

  }//end loop over pfcands

  // cout<<"total iso: "<<iso<<" sum of emiso + nhiso + chiso: "
  //     <<emiso<<" + "<<nhiso<<" + "<<chiso<<" = "<<emiso+nhiso+chiso
  //     <<endl;
  //  return iso;
  std::vector<float> isos;
  isos.push_back(chiso);
  isos.push_back(emiso);
  isos.push_back(nhiso);
  return isos;
}

pair<float,float> genLooper::getPhiCorrMET( float met, float metphi, int nvtx, bool ismc ) {

  //using met phi corrections from C. Veelken (emails from Oct. 4th)
  //previous versions are available here:
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py

  // Data
  // ------
  // x :  "+2.87340e-01 + 3.29813e-01*Nvtx" 
  // y : "-2.27938e-01 - 1.71272e-01*Nvtx"
  // MC
  // ------            
  // x : "+8.72683e-02 - 1.66671e-02*Nvtx"
  // y :  "+1.86650e-01 - 1.21946e-01*Nvtx"
  

  float metx = met * cos( metphi );
  float mety = met * sin( metphi );

  float shiftx = 0.;
  float shifty = 0.;

  //use correction for data vs. mc 
  shiftx = ismc ? (+8.72683e-02 - 1.66671e-02*nvtx)
    : (+2.87340e-01 + 3.29813e-01*nvtx);
  shifty = ismc ? (+1.86650e-01 - 1.21946e-01*nvtx)
    : (-2.27938e-01 - 1.71272e-01*nvtx);
  
  metx -= shiftx;
  mety -= shifty;

  pair<float, float> phicorrmet = make_pair( sqrt( metx*metx + mety*mety ), atan2( mety , metx ) );
  return phicorrmet;
}

pair<float,float> genLooper::getTrackerMET( P4 *lep, double deltaZCut, bool dolepcorr )
{

  if ( cms2.vtxs_sumpt().empty() ) return make_pair(-999.,-999.);

  float pX = 0.;
  float pY = 0.;

  if (dolepcorr ) {
    pX -= lep->px();
    pY -= lep->py();
  }
  
  for (unsigned int i=0; i<cms2.pfcands_particleId().size(); ++i){
        if ( cms2.pfcands_charge().at(i)==0 ) continue;

	if ( dolepcorr && dRbetweenVectors( cms2.pfcands_p4().at(i) , *lep ) < 0.1 ) continue;

        int trkIndex = cms2.pfcands_trkidx().at(i);
        if (trkIndex<0) continue;
        double dzpv = dzPV(cms2.trks_vertex_p4()[trkIndex], cms2.trks_trk_p4()[trkIndex], cms2.vtxs_position().front());

        if ( fabs(dzpv) > deltaZCut) continue;

        pX -= cms2.pfcands_p4().at(i).px();
        pY -= cms2.pfcands_p4().at(i).py();
    }

  pair<float, float> trkmet = make_pair( sqrt( pX*pX + pY*pY ), atan2( pY , pX ) );
  return trkmet;
  
}

float getdltrigweight(int id1, int id2)
{ 
  if (abs(id1)==11 && abs(id2)==11) return 0.95;
  if (abs(id1)==13 && abs(id2)==13) return 0.88;
  if (abs(id1)!=abs(id2)) return 0.92;
  return -999.;

}

float getsltrigweight(int id1, float pt, float eta) 
{

  //electron efficiencies
  if ( abs(id1)==11 ) {
    if ( fabs(eta)<1.5) {
      if ( pt>=20 && pt<22 ) return 0.00;
      if ( pt>=22 && pt<24 ) return 0.00;
      if ( pt>=24 && pt<26 ) return 0.00;
      if ( pt>=26 && pt<28 ) return 0.08;
      if ( pt>=28 && pt<30 ) return 0.61;
      if ( pt>=30 && pt<32 ) return 0.86;
      if ( pt>=32 && pt<34 ) return 0.88;
      if ( pt>=34 && pt<36 ) return 0.90;
      if ( pt>=36 && pt<38 ) return 0.91;
      if ( pt>=38 && pt<40 ) return 0.92;
      if ( pt>=40 && pt<50 ) return 0.94;
      if ( pt>=50 && pt<60 ) return 0.95;
      if ( pt>=60 && pt<80 ) return 0.96;
      if ( pt>=80 && pt<100 ) return 0.96;
      if ( pt>=100 && pt<150 ) return 0.96;
      if ( pt>=150 && pt<200 ) return 0.97;
      if ( pt>=200 ) return 0.97;
    } else if ( fabs(eta)>=1.5 && fabs(eta)<2.1) {
      if ( pt>=20 && pt<22 ) return 0.00;
      if ( pt>=22 && pt<24 ) return 0.00;
      if ( pt>=24 && pt<26 ) return 0.02;
      if ( pt>=26 && pt<28 ) return 0.18;
      if ( pt>=28 && pt<30 ) return 0.50;
      if ( pt>=30 && pt<32 ) return 0.63;
      if ( pt>=32 && pt<34 ) return 0.68;
      if ( pt>=34 && pt<36 ) return 0.70;
      if ( pt>=36 && pt<38 ) return 0.72;
      if ( pt>=38 && pt<40 ) return 0.74;
      if ( pt>=40 && pt<50 ) return 0.76;
      if ( pt>=50 && pt<60 ) return 0.77;
      if ( pt>=60 && pt<80 ) return 0.78;
      if ( pt>=80 && pt<100 ) return 0.80;
      if ( pt>=100 && pt<150 ) return 0.79;
      if ( pt>=150 && pt<200 ) return 0.76;
      if ( pt>=200 ) return 0.81;
    }
  } else if ( abs(id1)==13 ) {//muon efficiencies

    if ( fabs(eta)<0.8 ) {
      if (pt>=20 && pt<22)  return  0.00;	 
      if (pt>=22 && pt<24)  return  0.03; 	 
      if (pt>=24 && pt<26)  return  0.87; 
      if (pt>=26 && pt<28)  return  0.90; 
      if (pt>=28 && pt<30)  return  0.91; 
      if (pt>=30 && pt<32)  return  0.91; 
      if (pt>=32 && pt<34)  return  0.92; 
      if (pt>=34 && pt<36)  return  0.93; 
      if (pt>=36 && pt<38)  return  0.93; 
      if (pt>=38 && pt<40)  return  0.93; 
      if (pt>=40 && pt<50)  return  0.94; 
      if (pt>=50 && pt<60)  return  0.95; 
      if (pt>=60 && pt<80)  return  0.95; 
      if (pt>=80 && pt<100) return 0.94; 
      if (pt>=100 && pt<150) return 0.94; 
      if (pt>=150 && pt<200) return 0.93; 
      if (pt>=200) return 0.92; 
    } else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) {
      if (pt>=20 && pt<22)  return  0.00;
      if (pt>=22 && pt<24)  return  0.05;
      if (pt>=24 && pt<26)  return  0.78;
      if (pt>=26 && pt<28)  return  0.81;
      if (pt>=28 && pt<30)  return  0.81;
      if (pt>=30 && pt<32)  return  0.81;
      if (pt>=32 && pt<34)  return  0.82;
      if (pt>=34 && pt<36)  return  0.82;
      if (pt>=36 && pt<38)  return  0.83;
      if (pt>=38 && pt<40)  return  0.83;
      if (pt>=40 && pt<50)  return  0.84;
      if (pt>=50 && pt<60)  return  0.84;
      if (pt>=60 && pt<80)  return  0.84;
      if (pt>=80 && pt<100) return 0.84; 
      if (pt>=100 && pt<150) return 0.84;
      if (pt>=150 && pt<200) return 0.84;
      if (pt>=200) return 0.82;
    } else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) {
      if (pt>=20 && pt<22)  return  0.00;
      if (pt>=22 && pt<24)  return  0.11;
      if (pt>=24 && pt<26)  return  0.76;
      if (pt>=26 && pt<28)  return  0.78;
      if (pt>=28 && pt<30)  return  0.79;
      if (pt>=30 && pt<32)  return  0.80;
      if (pt>=32 && pt<34)  return  0.80;
      if (pt>=34 && pt<36)  return  0.81;
      if (pt>=36 && pt<38)  return  0.81;
      if (pt>=38 && pt<40)  return  0.82;
      if (pt>=40 && pt<50)  return  0.82;
      if (pt>=50 && pt<60)  return  0.83;
      if (pt>=60 && pt<80)  return  0.83;
      if (pt>=80 && pt<100) return 0.83;
      if (pt>=100 && pt<150) return 0.83;
      if (pt>=150 && pt<200) return 0.82;
      if (pt>=200) return 0.82;
    }
  }//end check for muons

  return 1.;

}
