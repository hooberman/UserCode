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
#include "singleLeptonLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"

#include "../CORE/CMS2.h"
#ifndef __CINT__
#include "../CORE/utilities.h"
#include "../CORE/ssSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/MITConversionUtilities.h"
#include "../CORE/muonSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/metSelections.h"
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
//#include "../Tools/bTagEff_BTV.h"
#endif

bool verbose        = false;
bool doTenPercent   = false;

//#include "../CORE/topmass/getTopMassEstimate.icc" // REPLACETOPMASS
//#include "../CORE/triggerUtils.cc"

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

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

singleLeptonLooper::singleLeptonLooper()
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

void singleLeptonLooper::InitBaby(){

  mutrigweight_ = 1.;

  jetid_	= 1;
  jetid30_	= 1;

  // btag variables
  nbtagsssv_	= 0;
  nbtagstcl_	= 0;
  nbtagstcm_	= 0;
  nbtagscsvl_   = 0;
  nbtagscsvm_   = 0;
  nbtagscsvt_   = 0;

  // njets with JEC variation
  njetsUp_	= 0;
  njetsDown_	= 0;
  htUp_		= 0.;
  htDown_	= 0.;

  // pfjet vars
  npfjets30_	= 0;
  npfjets35_	= 0;
  npfjets40_	= 0;
  npfjets45_	= 0;
  npfresjets30_	= 0;
  npfresjets35_	= 0;
  npfresjets40_	= 0;
  npfresjets45_	= 0;

  htpf30_	= 0.;
  htpf35_	= 0.;
  htpf40_	= 0.;
  htpf45_	= 0.;
  htpfres30_	= 0.;
  htpfres35_	= 0.;
  htpfres40_	= 0.;
  htpfres45_	= 0.;

  // calojet vars
  ncjets30_	= 0;
  ncjets35_	= 0;
  ncjets40_	= 0;
  ncjets45_	= 0;
  ncresjets30_	= 0;
  ncresjets35_	= 0;
  ncresjets40_	= 0;
  ncresjets45_	= 0;

  htc30_	= 0.;
  htc35_	= 0.;
  htc40_	= 0.;
  htc45_	= 0.;
  htcres30_	= 0.;
  htcres35_	= 0.;
  htcres40_	= 0.;
  htcres45_	= 0.;

  // MC truth info
  mcid1_	= -1;
  mcid2_	= -1;
  mclep1_	=  0;
  mclep2_	=  0;
  mctaud1_      =  0;
  mctaud2_      =  0;
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
  cjet1_	= 0;
  cjet2_	= 0;
  cjet3_	= 0;
  cjet4_	= 0;

  pfjet1_	= 0;
  pfjet2_	= 0;
  pfjet3_	= 0;
  pfjet4_	= 0;

  cresjet1_	= 0;
  cresjet2_	= 0;
  cresjet3_	= 0;
  cresjet4_	= 0;

  pfresjet1_	= 0;
  pfresjet2_	= 0;
  pfresjet3_	= 0;
  pfresjet4_	= 0;

  lep1_		= 0;
  lep2_		= 0;
  mbb_		= -9999.;

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

void singleLeptonLooper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

//--------------------------------------------------------------------

float singleLeptonLooper::stopPairCrossSection( float stopmass ){

  // stop mass divisible by 10
  if( ((int)stopmass%10)<1 ){
    int   bin  = stop_xsec_hist->FindBin(stopmass);
    float xsec = stop_xsec_hist->GetBinContent(bin);
    return xsec;
  }

  // stop mass not divisible by 10
  else{
    int   bin   = stop_xsec_hist->FindBin(stopmass);
    float xsec1 = stop_xsec_hist->GetBinContent(bin);
    float xsec2 = stop_xsec_hist->GetBinContent(bin+1);
    float xsec  = 0.5 * ( xsec1 + xsec2 );
    return xsec;
  }
}

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

float getMT( float leppt , float lepphi , float met , float metphi ) {
  float dphi = fabs( lepphi - metphi );
      if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
      return sqrt( 2 * ( leppt * met * (1 - cos( dphi ) ) ) );
}

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
  if ( abs(eta) < 0.8 ) i_eta = 0;
  else if ( abs(eta) >= 0.8 && abs(eta) < 1.5 ) i_eta = 1;
  else if ( abs(eta) >= 1.5 && abs(eta) < 2.1 ) i_eta = 2;
  if ( i_eta < 0 ) return 1.;

  return trigweights[i_pt][i_eta];

}

//--------------------------------------------------------------------

int singleLeptonLooper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale, float lumi,
				  FREnum frmode, bool doFakeApp)

{

  bool isLM = TString(prefix).Contains("LM");
  bool isData = false;
  if( TString(prefix).Contains("data") || TString(prefix).Contains("2011") ){
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
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    //set vtx reweighting hist
    set_vtxreweight_rootfile("vtxreweight_Summer11MC_PUS4_4p7fb_Zselection.root",true);

    //set msugra cross section file
    set_msugra_file("goodModelNames_tanbeta10.txt");


    initialized = true;
  }

  //------------------------------------------------------------------------------------------------------
  // latest-and-greatest JEC
  //------------------------------------------------------------------------------------------------------

  std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

  std::vector<std::string> jetcorr_filenames_caloL1OffsetL2L3;
  FactorizedJetCorrector *jet_corrector_caloL1OffsetL2L3;

  jetcorr_filenames_pfL1FastJetL2L3.clear();
  jetcorr_filenames_caloL1OffsetL2L3.clear();
  
  //string pfUncertaintyFile;
  //string caloUncertaintyFile;

  if ( isData ) {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L1FastJet.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2Relative.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L3Absolute.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2L3Residual.txt");

    //pfUncertaintyFile = "jetCorrections/GR_R_42_V23_AK5PF_Uncertainty.txt";

    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/GR_R_42_V23_AK5Calo_L1Offset.txt");
    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/GR_R_42_V23_AK5Calo_L2Relative.txt");
    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/GR_R_42_V23_AK5Calo_L3Absolute.txt");
    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/GR_R_42_V23_AK5Calo_L2L3Residual.txt");

    //caloUncertaintyFile = "jetCorrections/GR_R_42_V23_AK5Calo_Uncertainty.txt";
  } 
  else {
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L1FastJet.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L2Relative.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L3Absolute.txt");
    
    //pfUncertaintyFile = "jetCorrections/DESIGN42_V17_AK5PF_Uncertainty.txt";

    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/DESIGN42_V17_AK5Calo_L1Offset.txt");
    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/DESIGN42_V17_AK5Calo_L2Relative.txt");
    jetcorr_filenames_caloL1OffsetL2L3.push_back ("jetCorrections/DESIGN42_V17_AK5Calo_L3Absolute.txt");

    //caloUncertaintyFile = "jetCorrections/DESIGN42_V17_AK5Calo_Uncertainty.txt";
  }

  jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);
  jet_corrector_caloL1OffsetL2L3 = makeJetCorrector(jetcorr_filenames_caloL1OffsetL2L3);

  //JetCorrectionUncertainty *pfUncertainty   = new JetCorrectionUncertainty( pfUncertaintyFile   );
  //JetCorrectionUncertainty *caloUncertainty = new JetCorrectionUncertainty( caloUncertaintyFile );

  //------------------------------------------------
  // set stop cross section file
  //------------------------------------------------

  stop_xsec_file = TFile::Open("data/reference_xSec_stop.root");
  
  if( !stop_xsec_file->IsOpen() ){
    cout << "Error, could not open stop cross section TFile, quitting" << endl;
    exit(0);
  }
  
  stop_xsec_hist        = (TH1D*) stop_xsec_file->Get("stop");
  
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

      // skip stop-pair events with m(stop) > 850 GeV
      if( TString(prefix).Contains("T2") ){
	if( sparm_mG() > 800.0 ) continue;
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

	//cout << "Found electron " << ngoodlep_ << " pt " << els_p4().at(iel).pt() << endl;
      }
          
      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 10 )           continue;
	if( !muonId( imu , OSGeneric_v3 ))         continue;
	goodLeptons.push_back( mus_p4().at(imu) );
	lepId.push_back( mus_charge().at(imu) * 13 );
	lepIndex.push_back(imu);
	ngoodmu_++;
	ngoodlep_++;

	//cout << "Found muon " << ngoodlep_ << " pt " << mus_p4().at(imu).pt() << endl;
      }  

      // REQUIRE AT LEAST 1 GOOD LEPTON!!!
      if( goodLeptons.size() < 1 ) continue;

      //---------------------------------------------
      // find leading lepton
      //---------------------------------------------

      float maxpt   = -1;
      int   imaxpt  = -1;
      int   imaxpt2 = -1;

      for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	if( goodLeptons.at(ilep).pt() > maxpt ){
	  maxpt  = goodLeptons.at(ilep).pt();
	  imaxpt = ilep;
	}
      }

      if( imaxpt < 0 ){
	cout << "ERROR! QUITTING imaxpt " << imaxpt << endl;
	exit(2);
      }

      // REQUIRE LEADING LEPTON PT > 20 GEV
      if( maxpt < 20 ) continue;

      id1_       = lepId.at(imaxpt);
      lep1_      = &goodLeptons.at(imaxpt);
      int index1 = lepIndex.at(imaxpt);

      //cout << "Leading lepton: pt " << lep1_->pt() << " id " << id1_ << endl;

      //---------------------------------------------
      // find 2nd leading lepton (if >=2 leptons)
      //---------------------------------------------

      id2_       = -999;
      lep2_      = 0;
      int index2 = -1;
      dilmass_   = -999;

      if( ngoodlep_ > 1 ){

	maxpt = -1;

	for( unsigned int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){

	  if( (int)ilep == imaxpt ) continue;

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
	dilmass_   = sqrt((goodLeptons.at(imaxpt) + goodLeptons.at(imaxpt2)).mass2());

	//cout << "2nd deading lepton: pt " << lep2_->pt() << " id " << id2_ << " mass " << dilmass_ << endl;
      }

      //--------------------------------
      // store dilepton type in myType
      //--------------------------------

      leptype_ = 99;
      if      ( abs(id1_) == 11 )                   leptype_ = 0; // e
      else if ( abs(id1_) == 13 )                   leptype_ = 1; // m
      else{
	cout << "Skipping unknown lepton type = " << id1_ << endl;
	continue;
      }

      //--------------------------------
      // require trigger
      //--------------------------------

      int hypType = (leptype_==1) ? 0 : 3;//lepton type for dilepton triggers
      if( !passSingleLepSUSYTrigger2011_v1( isData , leptype_ ) && !passSUSYTrigger2011_v1( isData , hypType , true ) ) continue;

      //-----------------------------------------------------------------
      // number of OS generic electrons *in addition to* primary lepton
      //-----------------------------------------------------------------

      nosel_ = 0;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                                 continue;
	if( !pass_electronSelection( iel , electronSelection_el_OSV3 , false , false ) ) continue;
	if( dRbetweenVectors( *lep1_ , els_p4().at(iel) ) < 0.1 )                        continue;
	nosel_++;
      }
      
      //--------------------------------
      // get MC quantities
      //--------------------------------
      
      int nels       =  0;
      int nmus       =  0;
      int ntaus      =  0;
      int nleps      =  0;
      float dilptgen = -1;

      ptttbar_  = -1;
      ptt_      = -1;
      pttbar_   = -1;
      mttbar_   = -1;
      etattbar_ = -999;
      t_        = 0;
      tbar_     = 0;
      ttbar_    = 0;

      npartons_    =  0;
      maxpartonpt_ = -1;

      mgcor_ = 1.0;
      wflav_ = -1;

      if( !isData ){

	w1_     = leptonOrTauIsFromW( index1 , id1_ , isLM );
	pthat_  = genps_pthat();
	qscale_ = genps_qScale();
	
	//store W flavor history
	wflav_ = (int)genps_flavorHistoryFilterResult();

	//splitting ttbar into ttdil/ttotr
	//nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
	nleps = leptonGenpCount(nels, nmus, ntaus);
	
	nels_  = nels;
	nmus_  = nmus;
	ntaus_ = ntaus;
	nleps_ = nleps;

	// this is a weight which corrects for the wrong MG W->lnu BF
	if( TString(prefix).Contains("ttall") ){
	  if( nleps == 0 ) mgcor_ = 1.029;
	  if( nleps == 1 ) mgcor_ = 0.986;
	  if( nleps == 2 ) mgcor_ = 0.947;
	}

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

	  // skip lines up to t and tbar
	  if( igen < 8 ) continue;

	  // require particle is a quark or a gluon
	  int pid = abs( genps_id().at(igen) );
	  if( !( pid==1 || pid==2 || pid==3 || pid==4 || pid==5 || pid==6 || pid == 21 ) ) continue;

	  // require mother is not a top or W
	  int mothid = abs(genps_id_mother().at(igen));
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
            
	//-----------------------------------------------------
	// track gen lepton information here
	//-----------------------------------------------------
	
	mcid1_     = -1;
	mcid2_     = -1;
	mclep1_    =  0;
	mclep2_    =  0;
	mcdecay1_  = -1;
	mcdecay2_  = -1;
	mcdr1_     = -1;
	mcdr2_     = -1;
	mcndec1_   =  0;
	mcndec2_   =  0;
	mctaudpt1_ = -1;
	mctaudpt2_ = -1;
	mctaud1_   =  0;
	mctaud2_   =  0;
	mctaudid1_ = -1;
	mctaudid2_ = -1;

	//-----------------------------------------------------
	// store single gen lepton info
	//-----------------------------------------------------

	if( nleps_ == 1 ){


	  int nfoundleps = 0;

	  for ( int igen = 0 ; igen < (int)cms2.genps_id().size() ; igen++ ) { 

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    nfoundleps++;
	    mcid1_   = id;
	    mclep1_  = &genps_p4().at(igen);
	    mcdr1_   = dRbetweenVectors( *lep1_ , *mclep1_ );

	    if( abs(id)==15 ){
	      mcdecay1_ = 1;

	      for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++) {
		int daughter = abs(cms2.genps_lepdaughter_id()[igen][kk]);

		// get tau daughter pt
		if( daughter == 11 || daughter == 13 || daughter == 211 || daughter == 321 ){
		  if( genps_lepdaughter_p4()[igen][kk].pt() > mctaudpt1_ ) {
		    mctaudpt1_ = genps_lepdaughter_p4()[igen][kk].pt();
		    mctaudid1_ = cms2.genps_lepdaughter_id()[igen][kk];
		    mctaud1_   = &genps_lepdaughter_p4()[igen][kk];
		  }
		}

		if( daughter == 211 || daughter == 321 ) mcndec1_  ++;  // count charged hadrons
		if( daughter ==  12 || daughter ==  14 ) mcdecay1_ = 2; // check for nu_e or nu_mu 
	      } 
	    } 
	  } 


	  if( nfoundleps != 1 ) cout << "ERROR! expected 1 lepton, found " << nfoundleps << endl;
	}

	//-----------------------------------------------------	
	// store both gen lepton info
	//-----------------------------------------------------
	
	else if( nleps_ == 2 ){
	  
	  float drmin      = 9999;
	  int   igenmin    =   -1;
	  int   nfoundleps =    0;
	  
	  //-----------------------------------------------------
	  // find gen lepton closest to reco lepton
	  //-----------------------------------------------------
	  
	  for ( int igen = 0 ; igen < (int)cms2.genps_id().size() ; igen++ ) { 

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    nfoundleps++;

	    if( dRbetweenVectors( *lep1_ , genps_p4().at(igen) ) < drmin ){
	      drmin   = dRbetweenVectors( *lep1_ , genps_p4().at(igen) );
	      igenmin = igen;
	    }
	  }

	  if( nfoundleps != 2 ) cout << "ERROR! expected 2 leptons, found " << nfoundleps << endl;

	  //-----------------------------------------------------
	  // store info for closest gen lepton
	  //-----------------------------------------------------

	  mcid1_   = genps_id().at(igenmin);
	  mclep1_  = &genps_p4().at(igenmin);
	  mcdr1_   = dRbetweenVectors( *lep1_ , *mclep1_ );

	  if( abs(mcid1_)==15 ){
	    mcdecay1_ = 1;
	    
	    for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igenmin).size(); kk++) {
	      int daughter = abs(cms2.genps_lepdaughter_id()[igenmin][kk]);

	      // get tau daughter pt
	      if( daughter == 11 || daughter == 13 || daughter == 211 || daughter == 321 ){
		if( genps_lepdaughter_p4()[igenmin][kk].pt() > mctaudpt1_ ) {
		  mctaudpt1_ = genps_lepdaughter_p4()[igenmin][kk].pt();
		  mctaudid1_ = cms2.genps_lepdaughter_id()[igenmin][kk];
		  mctaud1_   = &genps_lepdaughter_p4()[igenmin][kk];
		}
	      }

	      if( daughter == 211 || daughter == 321 ) mcndec1_  ++;  // count charged hadrons
	      if( daughter ==  12 || daughter ==  14 ) mcdecay1_ = 2; // check for nu_e or nu_mu 
	    } 
	  }

	  //-----------------------------------------------------
	  // find 2nd lepton
	  //-----------------------------------------------------

	  int igenmin2 = -1;

	  for ( int igen = 0 ; igen < (int)cms2.genps_id().size() ; igen++ ) { 

	    if( igen == igenmin ) continue; //skip closest lepton

	    int id = genps_id().at(igen);

	    if( !( abs(id)==11 || abs(id)==13 || abs(id)==15 ) ) continue;

	    igenmin2 = igen;

	    mcid2_   = id;
	    mclep2_  = &genps_p4().at(igen);
	    if( ngoodlep_ > 1 ) mcdr2_   = dRbetweenVectors( *lep2_ , *mclep2_ );

	    if( abs(id)==15 ){
	      mcdecay2_ = 1;

	      for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++) {
		int daughter = abs(cms2.genps_lepdaughter_id()[igen][kk]);
		
		// get tau daughter pt
		if( daughter == 11 || daughter == 13 || daughter == 211 || daughter == 321 ){
		  if( genps_lepdaughter_p4()[igen][kk].pt() > mctaudpt1_ ) {
		    mctaudpt2_ = genps_lepdaughter_p4()[igen][kk].pt();
		    mctaudid2_ = cms2.genps_lepdaughter_id()[igen][kk];
		    mctaud2_   = &genps_lepdaughter_p4()[igen][kk];
		  }
		}

		if( daughter == 211 || daughter == 321 ) mcndec2_  ++;  // count charged hadrons
		if( daughter == 12  || daughter == 14  ) mcdecay2_ = 2; // check for nu_e or nu_mu
	      } 
	    } 
	  } 

	  if( igenmin2 < 0 ) cout << __FILE__ << " " << __LINE__ << " Error! unable to find 2nd gen lepton" << endl;

	  //-----------------------------------------------------
	  // find reco lepton corresponding to 2nd gen lepton
	  // check if found, if pass ID and/or iso
	  //-----------------------------------------------------

	  int  nMatchLeptons =  0;
	  int  imatch        = -1;
	  int  ID            = -1;

	  float drminlep = 999;

	  for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	    int mc3idx = els_mc3idx().at(iel);

	    if( mc3idx != igenmin2 ) continue;
	    nMatchLeptons++;

	    float dr = els_mc3dr().at(iel);

	    if( dr < drminlep ){
	      drminlep   = dr;
	      imatch     = iel;
	      ID         = 1;
	    }
	  }
          
	  for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	    int mc3idx = mus_mc3idx().at(imu);

	    if( mc3idx != igenmin2 ) continue;
	    nMatchLeptons++;

	    float dr = mus_mc3dr().at(imu);

	    if( dr < drminlep ){
	      drminlep   = dr;
	      imatch     = imu;
	      ID         = 2;
	    }

	  }

	  mlepid_       = -1;
	  mlep_         =  0;
	  mleppassid_   = -1;
	  mleppassiso_  = -1;
	  mlepiso_      = -1.0;
	  mlepdr_       = -1.0;

	  if( nMatchLeptons > 0 ){

	    // found matched electron
	    if( ID == 1 ){
	      mlepid_       = 11 * els_charge().at(imatch);
	      mlep_         = &els_p4().at(imatch);
	      mleppassid_   = pass_electronSelection( imatch , electronSelection_ssV5_noIso ) ? 1 : 0;
	      mleppassiso_  = pass_electronSelection( imatch , electronSelection_ssV5_iso   ) ? 1 : 0;
	      mlepiso_      = electronIsolation_rel_v1(imatch, true );
	    }

	    // found matched muon
	    else if( ID == 2 ){
	      mlepid_       = 13 * mus_charge().at(imatch);
	      mlep_         = &mus_p4().at(imatch);
	      mleppassid_   = muonIdNotIsolated( imatch , OSGeneric_v3 ) ? 1 : 0;
	      mleppassiso_  = muonIsoValue(imatch,false) < 0.15 ? 1 : 0;
	      mlepiso_      = muonIsoValue(imatch,false);
	    }

	    mlepdr_ = dRbetweenVectors( *mlep_ , *mclep2_ );
	  }
	  
	}

	else if( nleps_ < 0 || nleps_ > 2 ){
	  cout << "ERROR nleptons = " << nleps_ << endl;
	}

      }

      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = cms2.pfcands_trkidx().at(ipf);
	
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > 0.2 ){
 	    fillOverFlow( h_PU_trkpt , pfcands_p4().at(ipf).pt() );
 	  }
 	}
      }
      //------------------------------------------------------
      // store pt and iso for most isolated track (pt>10 GeV)
      //------------------------------------------------------

      trkpt10_         = -1.0;
      trkreliso10_     = 1000.;
      mleptrk10_       = -1.0;
      float miniso10   = 999;

      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( pfcands_p4().at(ipf).pt() < 10  ) continue;
	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = cms2.pfcands_trkidx().at(ipf);
	
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > 0.2 ) continue;
 	}

	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) 
	    isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	float iso = trackIso(ipf) / pfcands_p4().at(ipf).pt();

	if( iso < miniso10 ){
	  miniso10       = iso;
	  trkpt10_       = pfcands_p4().at(ipf).pt();
	  mleptrk10_     = (*lep1_+pfcands_p4().at(ipf)).pt();
	  trkreliso10_   = iso;
	}
      }

      //------------------------------------------------------
      // store pt and iso for most isolated track (pt>5 GeV)
      //------------------------------------------------------

      trkpt5_          = -1.0;
      trkreliso5_      = 1000.;
      mleptrk5_        = -1.0;
      float miniso5    = 999;

      for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

	if( pfcands_p4().at(ipf).pt() < 5   ) continue;
	if( pfcands_charge().at(ipf) == 0   ) continue;

 	int itrk = cms2.pfcands_trkidx().at(ipf);
	
 	if( itrk < (int)trks_trk_p4().size() && itrk >= 0 ){
 	  if( fabs( dz_trk_vtx(itrk,0) ) > 0.2 ) continue;
 	}

	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( pfcands_p4().at(ipf) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	float iso = trackIso(ipf) / pfcands_p4().at(ipf).pt();

	if( iso < miniso5 ){
	  miniso5     = iso;
	  trkpt5_     = pfcands_p4().at(ipf).pt();
	  mleptrk5_   = (*lep1_+pfcands_p4().at(ipf)).pt();
	  trkreliso5_ = iso;
	  //itrk       = ipf;
	}
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

      //-------------------------------------
      // jet counting
      //-------------------------------------

      VofP4 mediumBJets;

      int   imaxjet   = -1;
      float maxjetpt  = -1.;

      VofP4 vpfjets_p4;
      VofP4 vpfrawjets_p4;
      VofP4 vpfresjets_p4;
      vpfjets_p4.clear();
      vpfrawjets_p4.clear();
      vpfresjets_p4.clear();

      vector<float> fullcors;
      vector<float> rescors;
      vector<float> l1cors;
      fullcors.clear();
      rescors.clear();
      l1cors.clear();

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

	// get L1FastL2L3Residual total correction
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();

	// get L1Fast, L2, L3, Residual individual corrections
	jet_corrector_pfL1FastJetL2L3->setRho   ( cms2.evt_ww_rho_vor()           );
	jet_corrector_pfL1FastJetL2L3->setJetA  ( cms2.pfjets_area().at(ijet)     );
	jet_corrector_pfL1FastJetL2L3->setJetPt ( cms2.pfjets_p4().at(ijet).pt()  );
	jet_corrector_pfL1FastJetL2L3->setJetEta( cms2.pfjets_p4().at(ijet).eta() );
	vector<float> factors = jet_corrector_pfL1FastJetL2L3->getSubCorrections();

	// get residual correction only
	float rescorr = 1;
	if( isData ){
	  if( factors.size() == 4 ) rescorr = factors.at(3) / factors.at(2);
	  else                      cout << "ERROR! " << factors.size() << " jetSubCorrections" << endl;
	}

	LorentzVector vjet      = corr    * pfjets_p4().at(ijet);
	LorentzVector vresjet   = rescorr * pfjets_p4().at(ijet);
	LorentzVector vjetUp    = corr    * pfjets_p4().at(ijet) * 1.075; // over-estimate...
	LorentzVector vjetDown  = corr    * pfjets_p4().at(ijet) * 0.925; // over-estimate...

	// lepton-jet overlap removal
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
	          
	// PFJetID
	if( !passesPFJetID(ijet) ){
	  jetid_ = 0;
	  if( vjet.pt() > 30 && fabs( vjet.eta() ) < 2.5 ) jetid30_ = 0;
	  continue;
	}

        // store raw pfjet p4's and corrections for type1 pfmet
        // using corr jet pT > 10 GeV and adding |eta|<4.7 as in AN2011/459 to avoid problems with
        // erroneously large corrections
        if( vjet.pt() > 10 && abs(vjet.eta()) < 4.7 ){
          vpfrawjets_p4.push_back( pfjets_p4().at(ijet) );
          fullcors.push_back( corr );
          rescors.push_back( rescorr );
          l1cors.push_back( factors.at(0) );
        }

	// store L1FastL2L3Residual jet p4's pt > 15 GeV
	if( vjet.pt() > 15 && fabs( vjet.eta() ) < 2.5 ){
	  vpfjets_p4.push_back( vjet );
	}

	// store Residual jet p4's pt > 15 GeV
	if( vresjet.pt() > 15 && fabs( vresjet.eta() ) < 2.5 ){
	  vpfresjets_p4.push_back( vresjet );
	}

	// njets JEC up
	if( vjetUp.pt() > 30. && fabs( vjetUp.eta() ) < 2.5 ){
	  njetsUp_++;
	  htUp_ += vjetUp.pt();
	}

	// njets JEC down
	if( vjetDown.pt() > 30. && fabs( vjetDown.eta() ) < 2.5 ){
	  njetsDown_++;
	  htDown_ += vjetDown.pt();
	}

	// njets: residual only, pt > 30 GeV
	if( vresjet.pt() > 30 && fabs( vresjet.eta() ) < 2.5 ){
	  npfresjets30_ ++;
	  htpfres30_ += vresjet.pt();
	}

	// njets: residual only, pt > 35 GeV
	if( vresjet.pt() > 35 && fabs( vresjet.eta() ) < 2.5 ){
	  npfresjets35_ ++;
	  htpfres35_ += vresjet.pt();
	}

	// njets: residual only, pt > 40 GeV
	if( vresjet.pt() > 40 && fabs( vresjet.eta() ) < 2.5 ){
	  npfresjets40_ ++;
	  htpfres40_ += vresjet.pt();
	}

	// njets: residual only, pt > 45 GeV
	if( vresjet.pt() > 45 && fabs( vresjet.eta() ) < 2.5 ){
	  npfresjets45_ ++;
	  htpfres45_ += vresjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 30 GeV
	if( vjet.pt() > 30 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets30_ ++;
	  htpf30_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 35 GeV
	if( vjet.pt() > 35 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets35_ ++;
	  htpf35_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 40 GeV
	if( vjet.pt() > 40 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets40_ ++;
	  htpf40_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 45 GeV
	if( vjet.pt() > 45 && fabs( vjet.eta() ) < 2.5 ){
	  npfjets45_ ++;
	  htpf45_ += vjet.pt();
	}

	if(       vjet.pt()    < 30. )           continue;
	if( fabs( vjet.eta() ) > 2.5 )           continue;

	// btag variables: SSV
	if( pfjets_simpleSecondaryVertexHighEffBJetTag().at(ijet) > 1.74 ){
	  nbtagsssv_++;
	}

	// btag variables: TCHEL
	if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 1.7 ){
	  nbtagstcl_++;
	}

	// btag variables: TCHEM
	if( pfjets_trackCountingHighEffBJetTag().at(ijet) > 3.3 ){
	  nbtagstcm_++;
	}

	// btag variables: CSVL
	if( pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.244 ){
	  nbtagscsvl_++;
	}

	// btag variables: CSVM
	if( pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.679 ){
	  nbtagscsvm_++;
	  mediumBJets.push_back(vjet);
	}

	// btag variables: CSVT
	if( pfjets_combinedSecondaryVertexBJetTag().at(ijet) > 0.898 ){
	  nbtagscsvt_++;
	}

	// store max jet pt
	if( vjet.pt() > maxjetpt ){
	  maxjetpt = vjet.pt();
	  imaxjet  = ijet;
	}
      }

      // type1 met's
      vector<float> empty(0);
      pair<float, float> p_t1met10    = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 10.0 );
      pair<float, float> p_t1met20    = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 20.0 );
      pair<float, float> p_t1met30    = Type1PFMET( vpfrawjets_p4 , fullcors , l1cors , 30.0 );
      pair<float, float> p_t1metres10 = Type1PFMET( vpfrawjets_p4 , rescors  , empty  , 10.0 );
      pair<float, float> p_t1metres20 = Type1PFMET( vpfrawjets_p4 , rescors  , empty  , 20.0 );
      pair<float, float> p_t1metres30 = Type1PFMET( vpfrawjets_p4 , rescors  , empty  , 30.0 );
      t1met10_       = p_t1met10.first;
      t1met20_       = p_t1met20.first;
      t1met30_       = p_t1met30.first;
      t1metres10_    = p_t1metres10.first;
      t1metres20_    = p_t1metres20.first;
      t1metres30_    = p_t1metres30.first;
      t1met10phi_    = p_t1met10.second;
      t1met20phi_    = p_t1met20.second;	  
      t1met30phi_    = p_t1met30.second;	  
      t1metres10phi_ = p_t1metres10.second;
      t1metres20phi_ = p_t1metres20.second;
      t1metres30phi_ = p_t1metres30.second;

      // store L1FastL2L3Residual pfjets
      sort(vpfjets_p4.begin(), vpfjets_p4.end(), sortByPt);
      if( vpfjets_p4.size() > 0 ) pfjet1_  = &vpfjets_p4.at(0);
      if( vpfjets_p4.size() > 1 ) pfjet2_  = &vpfjets_p4.at(1);
      if( vpfjets_p4.size() > 2 ) pfjet3_  = &vpfjets_p4.at(2);
      if( vpfjets_p4.size() > 3 ) pfjet4_  = &vpfjets_p4.at(3);

      // store Residual pfjets
      sort(vpfresjets_p4.begin(), vpfresjets_p4.end(), sortByPt);
      if( vpfresjets_p4.size() > 0 ) pfresjet1_  = &vpfresjets_p4.at(0);
      if( vpfresjets_p4.size() > 1 ) pfresjet2_  = &vpfresjets_p4.at(1);
      if( vpfresjets_p4.size() > 2 ) pfresjet3_  = &vpfresjets_p4.at(2);
      if( vpfresjets_p4.size() > 3 ) pfresjet4_  = &vpfresjets_p4.at(3);

      // max jet variables
      if( imaxjet > -1 ){ 
	LorentzVector vjetcorr = pfjets_corL1FastL2L3().at(imaxjet) * pfjets_p4().at(imaxjet);
	jet_ = &vjetcorr;

	LorentzVector vjetraw = pfjets_p4().at(imaxjet);

	ptjetraw_     = vjetraw.pt();
	ptjet23_      = pfjets_cor().at(imaxjet)           * vjetraw.pt();
	ptjetF23_     = pfjets_corL1FastL2L3().at(imaxjet) * vjetraw.pt();
	ptjetO23_     = pfjets_corL1L2L3().at(imaxjet)     * vjetraw.pt();
	//cosphijz_     = -1 * cos( vjetraw.phi() - hyp_p4()[hypIdx].phi() );
	
	LorentzVector vjet = pfjets_corL1FastL2L3().at(imaxjet) * pfjets_p4().at(imaxjet);
	dphijm_ = acos(cos(vjet.phi()-evt_pfmetPhi()));
      }

      emjet10_     = -1;
      emjet20_     = -1;

      //-------------------------------------------------------------
      // find jet with max EM fraction outside tracker acceptance
      //-------------------------------------------------------------

      for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
	
	LorentzVector vjet = jets_p4().at(ijet) * jets_corL1FastL2L3().at(ijet);

	if( fabs( vjet.eta() ) < 2.5 )         continue;

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;

	float emfrac = jets_emFrac().at(ijet);

	if( vjet.pt() < 10. ) continue;

	if( emfrac > emjet10_ ) emjet10_ = emfrac;

	if( vjet.pt() < 20. ) continue;

	if( emfrac > emjet20_ ) emjet20_ = emfrac;

      }

      //------------------------------------------
      // count calojets
      //------------------------------------------

      //VofP4 goodCaloJets;


      VofP4 vcjets_p4;     
      VofP4 vcresjets_p4;   

      vcjets_p4.clear();
      vcresjets_p4.clear();

      for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {

	// get total L1OffsetL2L3Residual corrections
	jet_corrector_caloL1OffsetL2L3->setNPV   ( ndavtx_                           );
	jet_corrector_caloL1OffsetL2L3->setJetE  ( cms2.jets_p4().at(ijet).energy()  );
	jet_corrector_caloL1OffsetL2L3->setJetPt ( cms2.jets_p4().at(ijet).pt()      );
	jet_corrector_caloL1OffsetL2L3->setJetEta( cms2.jets_p4().at(ijet).eta()     );
	double corr = jet_corrector_caloL1OffsetL2L3->getCorrection();

	// get individual L1Offset, L2, L3, Residual corrections
	jet_corrector_caloL1OffsetL2L3->setNPV   ( ndavtx_                           );
	jet_corrector_caloL1OffsetL2L3->setJetE  ( cms2.jets_p4().at(ijet).energy()  );
	jet_corrector_caloL1OffsetL2L3->setJetPt ( cms2.jets_p4().at(ijet).pt()      );
	jet_corrector_caloL1OffsetL2L3->setJetEta( cms2.jets_p4().at(ijet).eta()     );
	vector<float> factors = jet_corrector_caloL1OffsetL2L3->getSubCorrections();

	// get residual correction only
	float rescorr = 1;
	if( isData ){
	  if( factors.size() == 4 ) rescorr = factors.at(3) / factors.at(2);
	  else                      cout << "ERROR! " << factors.size() << " jetSubCorrections" << endl;
	}

	LorentzVector vjet    = jets_p4().at(ijet) * corr;
	LorentzVector vresjet = jets_p4().at(ijet) * rescorr;

	// jet-lepton overlap removal
	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;

	// jet ID
	if( !passesCaloJetID( vjet ) )         continue;	

	// store L1FastL2L3Residual jet p4's pt > 15 GeV
	if( vjet.pt() > 15 && fabs( vjet.eta() ) < 2.5 ){
	  vcjets_p4.push_back( vjet );
	}

	// store Residual jet p4's pt > 15 GeV
	if( vresjet.pt() > 15 && fabs( vresjet.eta() ) < 2.5 ){
	  vcresjets_p4.push_back( vresjet );
	}

	// njets: residual only, pt > 30 GeV
	if( vresjet.pt() > 30 && fabs( vresjet.eta() ) < 2.5 ){
	  ncresjets30_ ++;
	  htcres30_ += vresjet.pt();
	}

	// njets: residual only, pt > 35 GeV
	if( vresjet.pt() > 35 && fabs( vresjet.eta() ) < 2.5 ){
	  ncresjets35_ ++;
	  htcres35_ += vresjet.pt();
	}

	// njets: residual only, pt > 40 GeV
	if( vresjet.pt() > 40 && fabs( vresjet.eta() ) < 2.5 ){
	  ncresjets40_ ++;
	  htcres40_ += vresjet.pt();
	}

	// njets: residual only, pt > 45 GeV
	if( vresjet.pt() > 45 && fabs( vresjet.eta() ) < 2.5 ){
	  ncresjets45_ ++;
	  htcres45_ += vresjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 30 GeV
	if( vjet.pt() > 30 && fabs( vjet.eta() ) < 2.5 ){
	  ncjets30_ ++;
	  htc30_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 35 GeV
	if( vjet.pt() > 35 && fabs( vjet.eta() ) < 2.5 ){
	  ncjets35_ ++;
	  htc35_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 40 GeV
	if( vjet.pt() > 40 && fabs( vjet.eta() ) < 2.5 ){
	  ncjets40_ ++;
	  htc40_ += vjet.pt();
	}

	// njets: L1FastL2L3Residual, pt > 45 GeV
	if( vjet.pt() > 45 && fabs( vjet.eta() ) < 2.5 ){
	  ncjets45_ ++;
	  htc45_ += vjet.pt();
	}

	// // store non-b jets
	// if( jets_trackCountingHighEffBJetTag().at(ijet) < 3.3 ){
	//   goodCaloJets.push_back(vjet);
	// }
      }

      // store L1FastL2L3Residual cjets
      sort(vcjets_p4.begin(), vcjets_p4.end(), sortByPt);
      if( vcjets_p4.size() > 0 ) cjet1_  = &vcjets_p4.at(0);
      if( vcjets_p4.size() > 1 ) cjet2_  = &vcjets_p4.at(1);
      if( vcjets_p4.size() > 2 ) cjet3_  = &vcjets_p4.at(2);
      if( vcjets_p4.size() > 3 ) cjet4_  = &vcjets_p4.at(3);

      // store Residual cjets
      sort(vcresjets_p4.begin(), vcresjets_p4.end(), sortByPt);
      if( vcresjets_p4.size() > 0 ) cresjet1_  = &vcresjets_p4.at(0);
      if( vcresjets_p4.size() > 1 ) cresjet2_  = &vcresjets_p4.at(1);
      if( vcresjets_p4.size() > 2 ) cresjet3_  = &vcresjets_p4.at(2);
      if( vcresjets_p4.size() > 3 ) cresjet4_  = &vcresjets_p4.at(3);

      //--------------------------------
      // get non-isolated leptons
      //--------------------------------

      nonisoel_ = 0;
      nonisomu_ = 0;

      float maxelpt  = -1.;
      int   imaxelpt = -1;

      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){

	if( els_p4().at(iel).pt() < 10 )                                                    continue;
	if( !pass_electronSelection( iel , electronSelection_ssV5_noIso , false , false ) ) continue;

	// don't count the leptons that we already counted as good
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( els_p4().at(iel) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	// don't count leptons near b-jets (TCHEM)
	bool nearBJet = false;
	for( int ijet = 0 ; ijet < (int)mediumBJets.size() ; ijet++ ){
	  if( dRbetweenVectors( els_p4().at(iel) , mediumBJets.at(ijet) ) < 0.4 ) nearBJet = true;
	}
	if( nearBJet ) continue;
	
	if( els_p4().at(iel).pt() > maxelpt ){
	  maxelpt  = els_p4().at(iel).pt();
	  imaxelpt = iel;
	}
      }

      if( imaxelpt >= 0 ) nonisoel_ = &(els_p4().at(imaxelpt));

      float maxmupt  = -1.;
      int   imaxmupt = -1;

      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){

	if( mus_p4().at(imu).pt() < 10 )                   continue;
	if( !muonIdNotIsolated(imu , OSGeneric_v3 ) )      continue;

	// don't count the leptons that we already counted as good
	bool isGoodLepton = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( mus_p4().at(imu) , goodLeptons.at(ilep) ) < 0.1 ) isGoodLepton = true;  
	}
	if( isGoodLepton ) continue;

	// don't count leptons near b-jets (TCHEM)
	bool nearBJet = false;
	for( int ijet = 0 ; ijet < (int)mediumBJets.size() ; ijet++ ){
	  if( dRbetweenVectors( mus_p4().at(imu) , mediumBJets.at(ijet) ) < 0.4 ) nearBJet = true;
	}
	if( nearBJet ) continue;
	
	if( mus_p4().at(imu).pt() > maxmupt ){
	  maxmupt  = mus_p4().at(imu).pt();
	  imaxmupt = imu;
	}
      }

      if( imaxmupt >= 0 ) nonisomu_ = &(mus_p4().at(imaxmupt));

      //---------------------------------
      // jet mass variables
      //---------------------------------

      // mjj_ = -1;

      // if( nbctcm_ >= 2 && goodCaloJets.size() >= 2 ){
      // 	sort( goodCaloJets.begin(), goodCaloJets.end(), sortByPt);
      // 	mjj_ = ( goodCaloJets.at(0) + goodCaloJets.at(1) ).mass();
      // }

      //---------------------------------
      // L1offset jets
      //---------------------------------
      
      htoffset_    = 0.;
      njetsoffset_ = 0;

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
	LorentzVector vjet = pfjets_corL1L2L3().at(ijet) * pfjets_p4().at(ijet);

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( fabs( vjet.eta() ) > 2.5 )           continue;
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

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( fabs( vjet.eta() ) > 2.5 )           continue;
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

	  bool rejectJet = false;
	  for( int ilep = 0 ; ilep < (int)goodLeptons.size() ; ilep++ ){
	    if( dRbetweenVectors( vgjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	  }
	  if( rejectJet ) continue;
	    
	  if( vgjet.pt() < 30.                   )  continue;
	  if( fabs( vgjet.eta() ) > 2.5          )  continue;
	    
	  ngenjets_++;
	  htgen_ += vgjet.pt();
	}
      }
     
      if( !isData ){
	genmet_     = gen_met();
	gensumet_   = gen_sumEt();
	genmetphi_  = gen_metPhi();

      }

      //----------------------------
      // MET flavors
      //----------------------------

      pair<float, float> p_met; //met and met phi
      //p_met = getMet( "tcMET"    , hypIdx);
      p_met = make_pair( evt_tcmet() , evt_tcmetPhi() );

      tcmet_    = p_met.first;
      tcmetphi_ = p_met.second;
      tcsumet_  = evt_tcsumet();

      pfmet_    = evt_pfmet();
      pfmetphi_ = evt_pfmetPhi();
      pfsumet_  = evt_pfsumet();

      //p_met = getMet( "pfMET"    , hypIdx);
      p_met = make_pair( evt_pfmet() , evt_pfmetPhi() );      

      // pair<float, float> pfmetUp   = ScaleMET( p_met , hyp_p4().at(hypIdx) , 1.075 );
      // pair<float, float> pfmetDown = ScaleMET( p_met , hyp_p4().at(hypIdx) , 0.925 );
      // pair<float, float> pfmetTest = ScaleMET( p_met , hyp_p4().at(hypIdx) , 1.000 );

      // pfmetUp_      = pfmetUp.first;
      // pfmetDown_    = pfmetDown.first;
      // pfmetTest_    = pfmetTest.first;

      meff_ = ht_ + pfmet_ + lep1_->pt();

      m_events.insert(pair<int,int>(evt_event(), 1));

      //---------------------------
      // set event weight
      //---------------------------

      weight_ = -1.;

      if( TString(prefix).Contains("T2") ){
	mG_ = sparm_mG();
	mL_ = sparm_mL();
	x_  = sparm_mf();

	weight_ = lumi * stopPairCrossSection(mG_) * (1000./50000.);

	if( doTenPercent )	  weight_ *= 10;

 	if( TString(prefix).Contains("T2bw") ){
 	  if( fabs(x_-0.75) < 0.01 ) weight_ *= 5./4.;
 	}
      }

      else if(strcmp(prefix,"LMscan") == 0){

	m0_  = sparm_m0();
	m12_ = sparm_m12();

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

	if( doTenPercent )	  weight_ *= 10;

	if( TString(prefix).Contains("LM") ){
	  if( strcmp( prefix , "LM0" )  == 0 ) weight_ *= kfactorSUSY( "lm0" );
	  if( strcmp( prefix , "LM1" )  == 0 ) weight_ *= kfactorSUSY( "lm1" );
	  if( strcmp( prefix , "LM2" )  == 0 ) weight_ *= kfactorSUSY( "lm2" );
	  if( strcmp( prefix , "LM3" )  == 0 ) weight_ *= kfactorSUSY( "lm3" );
	  if( strcmp( prefix , "LM4" )  == 0 ) weight_ *= kfactorSUSY( "lm4" );
	  if( strcmp( prefix , "LM5" )  == 0 ) weight_ *= kfactorSUSY( "lm5" );
	  if( strcmp( prefix , "LM6" )  == 0 ) weight_ *= kfactorSUSY( "lm6" );
	  if( strcmp( prefix , "LM7" )  == 0 ) weight_ *= kfactorSUSY( "lm7" );
	  if( strcmp( prefix , "LM8" )  == 0 ) weight_ *= kfactorSUSY( "lm8" );
	  if( strcmp( prefix , "LM9" )  == 0 ) weight_ *= kfactorSUSY( "lm9" );
	  if( strcmp( prefix , "LM10" ) == 0 ) weight_ *= kfactorSUSY( "lm10");
	  if( strcmp( prefix , "LM11" ) == 0 ) weight_ *= kfactorSUSY( "lm11");
	  if( strcmp( prefix , "LM12" ) == 0 ) weight_ *= kfactorSUSY( "lm12");
	  if( strcmp( prefix , "LM13" ) == 0 ) weight_ *= kfactorSUSY( "lm13");
	}
      }


      //tranverse mass leading lepton & met
      dphilm_ = fabs( lep1_->phi() - pfmetphi_ );
      if( dphilm_ > TMath::Pi() ) dphilm_ = TMath::TwoPi() - dphilm_;

      mt_ = sqrt( 2 * ( lep1_->pt() * pfmet_ * (1 - cos( dphilm_ ) ) ) );

      //transverse mass for leading lepton & type1 mets
      t1met10mt_    = getMT( lep1_->pt() , lep1_->phi() , t1met10_ , t1met10phi_ );
      t1met20mt_    = getMT( lep1_->pt() , lep1_->phi() , t1met20_ , t1met20phi_ );
      t1met30mt_    = getMT( lep1_->pt() , lep1_->phi() , t1met30_ , t1met30phi_ );
      t1metres10mt_ = getMT( lep1_->pt() , lep1_->phi() , t1metres10_ , t1metres10phi_ );
      t1metres20mt_ = getMT( lep1_->pt() , lep1_->phi() , t1metres20_ , t1metres20phi_ );
      t1metres30mt_ = getMT( lep1_->pt() , lep1_->phi() , t1metres30_ , t1metres30phi_ );

      //dijet mass two bs highest pT b-tagged jets
      if (mediumBJets.size()>1) {
	mbb_ = (mediumBJets.at(0)+mediumBJets.at(1)).M();
      }

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

      //other vars for baby

      pass_         = ( npfjets30_ >= 3 && pfmet_ > 25. ) ? 1 : 0;
      proc_         = getProcessType(prefix);       //integer specifying sample
      topmass_      = -999;//topMass;                      //topepton mass //REPLACE TOPMASS
      y_	    = pfmet_ / sqrt( ht_ ); //y=MET/sqrt(HT)

      strcpy(dataset_, cms2.evt_dataset().Data());  //dataset name
      run_          = evt_run();                    //run
      lumi_         = evt_lumiBlock();              //lumi
      event_        = evt_event();                  //event
      ndavtxweight_ = vtxweight(isData,true);

      hbhe_         = evt_hbheFilter();

      k_ = 1;
      if     ( strcmp( prefix , "LM0"  )    == 0 ) k_ = kfactorSUSY( "lm0"  );
      else if( strcmp( prefix , "LM1"  )    == 0 ) k_ = kfactorSUSY( "lm1"  );
      else if( strcmp( prefix , "LM2"  )    == 0 ) k_ = kfactorSUSY( "lm2"  );
      else if( strcmp( prefix , "LM3"  )    == 0 ) k_ = kfactorSUSY( "lm3"  );
      else if( strcmp( prefix , "LM4"  )    == 0 ) k_ = kfactorSUSY( "lm4"  );
      else if( strcmp( prefix , "LM5"  )    == 0 ) k_ = kfactorSUSY( "lm5"  );
      else if( strcmp( prefix , "LM6"  )    == 0 ) k_ = kfactorSUSY( "lm6"  );
      else if( strcmp( prefix , "LM7"  )    == 0 ) k_ = kfactorSUSY( "lm7"  );
      else if( strcmp( prefix , "LM8"  )    == 0 ) k_ = kfactorSUSY( "lm8"  );
      else if( strcmp( prefix , "LM9"  )    == 0 ) k_ = kfactorSUSY( "lm9"  );
      else if( strcmp( prefix , "LM10" )    == 0 ) k_ = kfactorSUSY( "lm10" );
      else if( strcmp( prefix , "LM11" )    == 0 ) k_ = kfactorSUSY( "lm11" );
      else if( strcmp( prefix , "LM12" )    == 0 ) k_ = kfactorSUSY( "lm12" );
      else if( strcmp( prefix , "LMscan" )  == 0 ) k_ = kfactorSUSY(m0_,m12_,"tanbeta10");
      

      ecalveto1_ = -1.;
      ecalveto2_ = -1.;
      hcalveto1_ = -1.;
      hcalveto2_ = -1.;
      
      if( leptype_ == 0 ){
	iso1_   = electronIsolation_rel   ( index1 , true ); //truncated
	isont1_ = electronIsolation_rel_v1( index1 , true ); //non-truncated
	etasc1_ = els_etaSC()[index1];
      }
      else if( leptype_ == 1 ){
	iso1_   = muonIsoValue( index1 , true  ); //truncated 
	isont1_ = muonIsoValue( index1 , false ); //non-truncated
	etasc1_ = -999;
	
	ecalveto1_ = mus_iso_ecalvetoDep().at(index1);
	hcalveto1_ = mus_iso_hcalvetoDep().at(index1);
      }
            
      if     ( leptype_ == 0 ) netot += weight_;
      else if( leptype_ == 1 ) nmtot += weight_;

      if( pass_ == 1 ){
	if     ( leptype_ == 0 ) nepass += weight_;
	else if( leptype_ == 1 ) nmpass += weight_;
      }

      //-------------------------------------
      // triggers
      //-------------------------------------

      //eltrijet_ = passHLTTrigger("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v3") ? 1 : 0;
      //mutrijet_ = passHLTTrigger("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v3") ? 1 : 0;

      ldi_  = passSingleLep2JetSUSYTrigger2011( isData , leptype_ ) ? 1 : 0;
      ltri_ = passSingleLep3JetSUSYTrigger2011( isData , leptype_ ) ? 1 : 0;
      smu_  = passSingleMuTrigger2011(          isData , leptype_ ) ? 1 : 0;
      dil_  = passSUSYTrigger2011_v1(     isData , hypType , true ) ? 1 : 0;

      //set trigger weight
      mutrigweight_ = getMuTriggerWeight( lep1_->pt() , lep1_->eta() );
      
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
 
void singleLeptonLooper::BookHistos(char *prefix)
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

//*****************************************************************
// get the FR weight
//*****************************************************************

double singleLeptonLooper::getFRWeight(const int hypIdx, SimpleFakeRate* mufr, SimpleFakeRate * elfr, FREnum frmode, bool isData) {

  //std::cout<<"Called singleLeptonLooper::getFRWeight"<<std::endl;

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
    std::cout<<"singleLeptonLooper::getFRWeight: bad FR mode given, fix this!"<<std::endl;
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

    if( pass_electronSelection( iElt , electronSelection_ssV5 ) ) {
      isGoodElt = true;
    }
    if( pass_electronSelection( iEll , electronSelection_ssV5 ) ) {
      isGoodEll = true;
    }
    if( pass_electronSelection( iElt , electronSelectionFOV5_ssVBTF80_v1 ) ) {
      isFOElt   = true;
    }
    if( pass_electronSelection( iEll , electronSelectionFOV5_ssVBTF80_v1 ) ) {
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

    if( pass_electronSelection( iEl , electronSelection_ssV5 ) ){
      isGoodEl = true;
    }
    if( muonId( iMu , OSGeneric_v3 ) ) { 
      isGoodMu = true;
    }
    if( pass_electronSelection( iEl , electronSelectionFOV5_ssVBTF80_v1 ) ){
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

void singleLeptonLooper::makeTree(char *prefix, bool doFakeApp, FREnum frmode ){
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

  outFile   = new TFile(Form("output/%s_smallTree%s%s.root",prefix,frsuffix,tpsuffix), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix,frsuffix,tpsuffix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in singleLeptonLooper.h
  outTree->Branch("acc_2010",        &acc_2010_,         "acc_2010/I");
  outTree->Branch("acc_highmet",     &acc_highmet_,      "acc_highmet/I");
  outTree->Branch("acc_highht",      &acc_highht_,       "acc_highht/I");
  outTree->Branch("hbhe",            &hbhe_,             "hbhe/I");
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
  outTree->Branch("maxpartonpt",     &maxpartonpt_,      "maxpartonpt/F");
  outTree->Branch("etattbar",        &etattbar_,         "etatbar/F");
  outTree->Branch("njetsoffset",     &njetsoffset_,      "njetsoffset/I");
  outTree->Branch("njetsuncor",      &njetsuncor_,       "njetsuncor/I");
  outTree->Branch("costhetaweight",  &costhetaweight_,   "costhetaweight/F");
  outTree->Branch("weight",          &weight_,           "weight/F");
  outTree->Branch("mutrigweight",    &mutrigweight_,     "mutrigweight/F");
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
  outTree->Branch("nlep",            &nlep_,             "nlep/I");
  outTree->Branch("nosel",           &nosel_,            "nosel/I");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("ngoodmu",         &ngoodmu_,          "ngoodmu/I");
  outTree->Branch("mull",            &mull_,             "mull/I");
  outTree->Branch("mult",            &mult_,             "mult/I");
  //outTree->Branch("eltrijet",        &eltrijet_,         "eltrijet/I");
  //outTree->Branch("mutrijet",        &mutrijet_,         "mutrijet/I");
  outTree->Branch("ldi",             &ldi_,              "ldi/I");
  outTree->Branch("ltri",            &ltri_,             "ltri/I");
  outTree->Branch("smu",             &smu_,              "smu/I");
  outTree->Branch("dil",             &dil_,              "dil/I");
  outTree->Branch("mullgen",         &mullgen_,          "mullgen/I");
  outTree->Branch("multgen",         &multgen_,          "multgen/I");
  outTree->Branch("proc",            &proc_,             "proc/I");
  outTree->Branch("leptype",         &leptype_,          "leptype/I");
  outTree->Branch("topmass",         &topmass_,          "topmass/F");
  outTree->Branch("dilmass",         &dilmass_,          "dilmass/F");
  outTree->Branch("tcmet",           &tcmet_,            "tcmet/F");
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
  outTree->Branch("ngenjets",        &ngenjets_,         "ngenjets/I");
  outTree->Branch("njpt",            &njpt_,             "njpt/I");

  // pfjets L1FastL2L3Res
  outTree->Branch("npfjets30",       &npfjets30_,        "npfjets30/I");
  outTree->Branch("npfjets35",       &npfjets35_,        "npfjets35/I");
  outTree->Branch("npfjets40",       &npfjets40_,        "npfjets40/I");
  outTree->Branch("npfjets45",       &npfjets45_,        "npfjets45/I");

  outTree->Branch("htpf30",          &htpf30_,           "htpf30/F");
  outTree->Branch("htpf35",          &htpf35_,           "htpf35/F");
  outTree->Branch("htpf40",          &htpf40_,           "htpf40/F");
  outTree->Branch("htpf45",          &htpf45_,           "htpf45/F");

  // type1 met flavors
  outTree->Branch("t1met10",         &t1met10_,          "t1met10/F");
  outTree->Branch("t1met20",         &t1met20_,          "t1met20/F");
  outTree->Branch("t1met30",         &t1met30_,          "t1met30/F");
  outTree->Branch("t1metres10",      &t1metres10_,       "t1metres10/F");
  outTree->Branch("t1metres20",      &t1metres20_,       "t1metres20/F");
  outTree->Branch("t1metres30",      &t1metres30_,       "t1metres30/F");
  outTree->Branch("t1met10phi",      &t1met10phi_,       "t1met10phi/F");
  outTree->Branch("t1met20phi",      &t1met20phi_,       "t1met20phi/F");
  outTree->Branch("t1met30phi",      &t1met30phi_,       "t1met30phi/F");
  outTree->Branch("t1metres10phi",   &t1metres10phi_,    "t1metres10phi/F");
  outTree->Branch("t1metres20phi",   &t1metres20phi_,    "t1metres20phi/F");
  outTree->Branch("t1metres30phi",   &t1metres30phi_,    "t1metres30phi/F");
  outTree->Branch("t1met10mt",       &t1met10mt_,       "t1met10mt/F");
  outTree->Branch("t1met20mt",       &t1met20mt_,       "t1met20mt/F");
  outTree->Branch("t1met30mt",       &t1met30mt_,       "t1met30mt/F");
  outTree->Branch("t1metres10mt",    &t1metres10mt_,    "t1metres10mt/F");
  outTree->Branch("t1metres20mt",    &t1metres20mt_,    "t1metres20mt/F");
  outTree->Branch("t1metres30mt",    &t1metres30mt_,    "t1metres30mt/F");
  
  // pfjets Res
  outTree->Branch("npfresjets30",    &npfresjets30_,     "npfresjets30/I");
  outTree->Branch("npfresjets35",    &npfresjets35_,     "npfresjets35/I");
  outTree->Branch("npfresjets40",    &npfresjets40_,     "npfresjets40/I");
  outTree->Branch("npfresjets45",    &npfresjets45_,     "npfresjets45/I");

  outTree->Branch("htpfres30",       &htpfres30_,        "htpfres30/F");
  outTree->Branch("htpfres35",       &htpfres35_,        "htpfres35/F");
  outTree->Branch("htpfres40",       &htpfres40_,        "htpfres40/F");
  outTree->Branch("htpfres45",       &htpfres45_,        "htpfres45/F");

  /// calojets L1OffsetL2L3Res
  outTree->Branch("ncjets30",        &ncjets30_,         "ncjets30/I");
  outTree->Branch("ncjets35",        &ncjets35_,         "ncjets35/I");
  outTree->Branch("ncjets40",        &ncjets40_,         "ncjets40/I");
  outTree->Branch("ncjets45",        &ncjets45_,         "ncjets45/I");

  outTree->Branch("htc30",           &htc30_,            "htc30/F");
  outTree->Branch("htc35",           &htc35_,            "htc35/F");
  outTree->Branch("htc40",           &htc40_,            "htc40/F");
  outTree->Branch("htc45",           &htc45_,            "htc45/F");

  /// calojets Res
  outTree->Branch("ncresjets30",     &ncresjets30_,      "ncresjets30/I");
  outTree->Branch("ncresjets35",     &ncresjets35_,      "ncresjets35/I");
  outTree->Branch("ncresjets40",     &ncresjets40_,      "ncresjets40/I");
  outTree->Branch("ncresjets45",     &ncresjets45_,      "ncresjets45/I");

  outTree->Branch("htcres30",        &htcres30_,         "htcres30/F");
  outTree->Branch("htcres35",        &htcres35_,         "htcres35/F");
  outTree->Branch("htcres40",        &htcres40_,         "htcres40/F");
  outTree->Branch("htcres45",        &htcres45_,         "htcres45/F");

  // btag variables
  outTree->Branch("nbtagsssv",       &nbtagsssv_,        "nbtagsssv/I");
  outTree->Branch("nbtagstcl",       &nbtagstcl_,        "nbtagstcl/I");
  outTree->Branch("nbtagstcm",       &nbtagstcm_,        "nbtagstcm/I");
  outTree->Branch("nbtagscsvl",      &nbtagscsvl_,       "nbtagscsvl/I");
  outTree->Branch("nbtagscsvm",      &nbtagscsvm_,       "nbtagscsvm/I");
  outTree->Branch("nbtagscsvt",      &nbtagscsvt_,       "nbtagscsvt/I");

  outTree->Branch("njetsUp",         &njetsUp_,          "njetsUp/I");
  outTree->Branch("njetsDown",       &njetsDown_,        "njetsDown/I");
  outTree->Branch("htUp",            &htUp_,             "htUp/F");
  outTree->Branch("htDown",          &htDown_,           "htDown/F");
  outTree->Branch("nvtx",            &nvtx_,             "nvtx/I");
  outTree->Branch("ndavtx",          &ndavtx_,           "ndavtx/I");
  outTree->Branch("ndavtxweight",    &ndavtxweight_,     "ndavtxweight/F");
  outTree->Branch("vecjetpt",        &vecjetpt_,         "vecjetpt/F");
  outTree->Branch("pass",            &pass_,             "pass/I");
  outTree->Branch("passz",           &passz_,            "passz/I");
  outTree->Branch("m0",              &m0_,               "m0/F");
  outTree->Branch("mg",              &mG_,               "mg/F");
  outTree->Branch("ml",              &mL_,               "ml/F");
  outTree->Branch("x",               &x_,                "x/F");
  outTree->Branch("m12",             &m12_,              "m12/F");
  outTree->Branch("id1",             &id1_,              "id1/I");
  outTree->Branch("id2",             &id2_,              "id2/I");
  outTree->Branch("w1",              &w1_,               "w1/I");
  outTree->Branch("w2",              &w2_,               "w2/I");
  outTree->Branch("iso1",            &iso1_,             "iso1/F");
  outTree->Branch("isont1",          &isont1_,           "isont1/F");
  outTree->Branch("etasc1",          &etasc1_,           "etasc1/F");
  outTree->Branch("etasc2",          &etasc2_,           "etasc2/F");
  outTree->Branch("iso2",            &iso2_,             "iso2/F");
  outTree->Branch("ecalveto1",       &ecalveto1_,        "ecalveto1/F");
  outTree->Branch("ecalveto2",       &ecalveto2_,        "ecalveto2/F");
  outTree->Branch("hcalveto1",       &hcalveto1_,        "hcalveto1/F");
  outTree->Branch("hcalveto2",       &hcalveto2_,        "hcalveto2/F");
  outTree->Branch("isont2",          &isont2_,           "isont2/F");
  outTree->Branch("ptl1",            &ptl1_,             "ptl1/F");
  outTree->Branch("ptl2",            &ptl2_,             "ptl2/F");
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
  outTree->Branch("htjpt",           &htjpt_,            "htjpt/F");  
  outTree->Branch("nels",            &nels_,             "nels/I");  
  outTree->Branch("nmus",            &nmus_,             "nmus/I");  
  outTree->Branch("ntaus",           &ntaus_,            "ntaus/I");  
  outTree->Branch("nleps",           &nleps_,            "nleps/I");  
  outTree->Branch("dphijm",          &dphijm_,           "dphijm/F");  
  outTree->Branch("ptjetraw",        &ptjetraw_,         "ptjetraw/F");  
  outTree->Branch("ptjet23",         &ptjet23_,          "ptjet23/F");  
  outTree->Branch("ptjetF23",        &ptjetF23_,         "ptjetF23/F");  
  outTree->Branch("ptjetO23",        &ptjetO23_,         "ptjetO23/F");  
  //outTree->Branch("cosphijz",        &cosphijz_,         "cosphijz/F");  
  outTree->Branch("mcid1",           &mcid1_,            "mcid1/I");  
  outTree->Branch("mcdr1",           &mcdr1_,            "mcdr1/F");  
  outTree->Branch("mcdecay1",        &mcdecay1_,         "mcdecay1/I");  
  outTree->Branch("mcndec1",         &mcndec1_,          "mcndec1/I");  
  outTree->Branch("mcndec2",         &mcndec2_,          "mcndec2/I");  
  outTree->Branch("mcid2",           &mcid2_,            "mcid2/I");  
  outTree->Branch("mcdr2",           &mcdr2_,            "mcdr2/F");  
  outTree->Branch("mcdecay2",        &mcdecay2_,         "mcdecay2/I");  
  outTree->Branch("mctaudpt1",       &mctaudpt1_,        "mctaudpt1/F");  
  outTree->Branch("mctaudpt2",       &mctaudpt2_,        "mctaudpt2/F");  
  outTree->Branch("mctaudid1",       &mctaudid1_,        "mctaudid1/I");  
  outTree->Branch("mctaudid2",       &mctaudid2_,        "mctaudid2/I");  
  outTree->Branch("mlepid",          &mlepid_,           "mlepid/I");  
  outTree->Branch("mleppassid",      &mleppassid_,       "mleppassid/I");  
  outTree->Branch("mleppassiso",     &mleppassiso_,      "mleppassiso/I");  
  outTree->Branch("mlepiso",         &mlepiso_,          "mlepiso/F");  
  outTree->Branch("mlepdr",          &mlepdr_,           "mlepdr/F");  
  outTree->Branch("emjet10",         &emjet10_,          "emjet10/F");  
  outTree->Branch("mjj",             &mjj_,              "mjj/F");  
  outTree->Branch("emjet20",         &emjet20_,          "emjet20/F");  
  outTree->Branch("trkpt5",          &trkpt5_,           "trkpt5/F");  
  outTree->Branch("trkpt10",         &trkpt10_,          "trkpt10/F");  
  outTree->Branch("mleptrk5",        &mleptrk5_,         "mleptrk5/F");  
  outTree->Branch("mleptrk10",       &mleptrk10_,        "mleptrk10/F");  
  outTree->Branch("trkreliso5",      &trkreliso5_,       "trkreliso5/F");  
  outTree->Branch("trkreliso10",     &trkreliso10_,      "trkreliso10/F");  
  outTree->Branch("mbb",             &mbb_,              "mbb/F");

  outTree->Branch("mlep"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mlep_	);
  outTree->Branch("lep1"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_	);
  outTree->Branch("lep2"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_	);
  outTree->Branch("mclep1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep1_	);
  outTree->Branch("mclep2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mclep2_	);
  outTree->Branch("mctaud1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaud1_	);
  outTree->Branch("mctaud2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &mctaud2_	);
  outTree->Branch("jet"	     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_	);

  outTree->Branch("cjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet1_	);
  outTree->Branch("cjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet2_	);
  outTree->Branch("cjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet3_	);
  outTree->Branch("cjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet4_	);

  outTree->Branch("pfjet1"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet1_	);
  outTree->Branch("pfjet2"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet2_	);
  outTree->Branch("pfjet3"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet3_	);
  outTree->Branch("pfjet4"   , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfjet4_	);

  outTree->Branch("cresjet1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cresjet1_	);
  outTree->Branch("cresjet2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cresjet2_	);
  outTree->Branch("cresjet3" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cresjet3_	);
  outTree->Branch("cresjet4" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cresjet4_	);

  outTree->Branch("pfresjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfresjet1_	);
  outTree->Branch("pfresjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfresjet2_	);
  outTree->Branch("pfresjet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfresjet3_	);
  outTree->Branch("pfresjet4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pfresjet4_	);

  outTree->Branch("nonisoel" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisoel_	);
  outTree->Branch("nonisomu" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &nonisomu_	);
  outTree->Branch("t"        , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &t_   	);
  outTree->Branch("tbar"     , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &tbar_   	);
  outTree->Branch("ttbar"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &ttbar_   	);


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

float singleLeptonLooper::dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

float singleLeptonLooper::trackIso( int thisPf , float coneR , float dz_thresh ){

  float iso = 0.0;

  for (int ipf = 0; ipf < (int)cms2.pfcands_p4().size(); ipf++) {

    if( ipf == thisPf                 ) continue; // skip this PFCandidate
    if( cms2.pfcands_charge().at(ipf) == 0 ) continue; // skip neutrals

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
      
    for (unsigned int ivtx = 0; ivtx < cms2.davtxs_position().size(); ivtx++) {
	
      if(!isGoodDAVertex(ivtx)) continue;

      float mydz = dz_trk_vtx(itrk,ivtx);
      
      if (fabs(mydz) < fabs(mindz)) {
	mindz = mydz;
	vtxi = ivtx;
      }
         
    }
    
    //----------------------------------------------------------------------------
    // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
    //----------------------------------------------------------------------------
    
    if ( vtxi != 0               )     continue;
    if ( fabs(mindz) > dz_thresh )     continue;

    //---------------------------------------
    // passes cuts, add up isolation value
    //---------------------------------------

    iso += cms2.pfcands_p4().at(ipf).pt();

  }

  return iso;
}

std::vector<float> singleLeptonLooper::totalIso( int thisPf , float coneR , float dz_thresh ){

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
    
      float mindz = 999.;
      int vtxi    = -1;
      
      for (unsigned int ivtx = 0; ivtx < cms2.davtxs_position().size(); ivtx++) {
	
	if(!isGoodDAVertex(ivtx)) continue;

	float mydz = dz_trk_vtx(itrk,ivtx);
      
	if (fabs(mydz) < fabs(mindz)) {
	  mindz = mydz;
	  vtxi = ivtx;
	}
         
      }
    
      //----------------------------------------------------------------------------
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
      //----------------------------------------------------------------------------
    
      if ( vtxi != 0               )     continue;
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