#ifndef Z_looper_h
#define Z_looper_h

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TProfile.h"
#include <vector>
#include <map>
#include <fstream>
#include "Math/LorentzVector.h"
//#include "Math/PxPyPzE4D.h"

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class TChain;

class Z_looper
{
 public:
  Z_looper() {};
  ~Z_looper() {
    delete babyFile_;
    delete babyTree_;
  };

  void MakeBabyNtuple (const char *);
  void InitBabyNtuple ();
  void FillBabyNtuple ();
  void CloseBabyNtuple ();
  void ScanChain (TChain*, const char*, bool isData, bool calculateTCMET = false,
                  int my_nEvents = -1, float kFactor = 1.);
  void bookHistos();
  bool isGoodTrack(int, bool usePV = false);
  float deltaPhi( float phi1 , float phi2);
  void fillUnderOverFlow(TH1F *h1, float value, float weight);
  void fillHistos(TH1F *h1[4],    float value, float weight, int myType);
  void fillHistos(TH1F *h1[4][4], float value, float weight, int myType, int nJetsIdx);
  float PassGenSelection( bool isData );
  float getMetError(  vector<int> goodMuonIndices );
  float getMetError_claudio(  vector<int> goodMuonIndices );
        
 private:
                
  //ntuple, file
  TFile *babyFile_;
  TTree *babyTree_;
    
  //histos

  // event stuff
  char    dataset_[200];
  Int_t   run_;
  Int_t   goodrun_;
  Int_t   lumi_;
  Int_t   event_;
  Int_t   leptype_;
  Int_t   ecaltype_;
  Int_t   nGoodVertex_;
  Int_t   nGoodDAVertex_;
  Float_t weight_;
  Float_t pthat_;
  Float_t qscale_;
  Float_t mllgen_;
  Float_t vtxweight_;
  Float_t davtxweight_;
  Float_t maxemf_;
  Float_t dpdm_;
  Float_t metError_;
  Float_t metErrorC_;
  Int_t   id1_;
  Int_t   id2_;

  Float_t ml_;
  Float_t mg_;

  LorentzVector*  lep1_;
  LorentzVector*  lep2_;
  LorentzVector*  dilep_;
  LorentzVector*  jet_; 

  // genmet stuff
  Float_t genmet_;
  Float_t genmetphi_;
  Float_t gensumet_;

  // pfmet stuff
  Float_t pfmet_;
  Float_t pfmetphi_;
  Float_t pfsumet_;

  //pfmuon stuff
  Int_t   npfmuons_;
  Int_t   nmatchedpfmuons_;
  Float_t ptll_pf_;
  Float_t ptlt_pf_;

  // calomet stuff
  Float_t met_;
  Float_t metphi_;
  Float_t sumet_;

  Float_t ptlltrk_;
  Float_t ptllgen_;
  Float_t ptltgen_;
  Float_t ptlttrk_;
  Float_t ptllgfit_;
  Float_t ptltgfit_;

  // muon-corrected calomet stuff
  Float_t mumet_;
  Float_t mumetphi_;
  Float_t musumet_;

  // muon-corrected JES calomet stuff
  Float_t mujesmet_;
  Float_t mujesmetphi_;
  Float_t mujessumet_;

  // tcmet stuff
  Float_t dphixmet_;
  Float_t metPar_;
  Float_t metPerp_;
  Float_t dphijetmet_;
  Float_t tcmet_;
  Float_t tcmetphi_;
  Float_t tcsumet_;

  Float_t tcmetNew_;
  Float_t tcmetphiNew_;
  Float_t tcsumetNew_;
  Float_t tcmetcor_;
  Float_t pfmetcor_;
  Float_t dilmasscor_;

  // jet stuff
  Int_t   failjetid_;
  Int_t   nJets_;
  Int_t   nJPT_;
  Int_t   nJets40_;
  Float_t sumJetPt_;
  Float_t sumJetPt10_;
  Float_t vecJetPt_;
  Int_t   nbtags_;

  //electron-matched jet stuff
  Float_t drjet_ll_;
  Float_t jetpt_ll_;
  Int_t   pfjetid_ll_;
  Float_t drjet_lt_;
  Float_t jetpt_lt_;
  Int_t   pfjetid_lt_;

  //leading jet stuff
  Float_t jetmax_pt_;
  Float_t jetmax_dphimet_;

  //Z stuff
  Int_t   passz_;
  Int_t   pdgid_;
  Int_t   idll_;
  Int_t   idlt_;
  Float_t ptll_;
  Float_t ptlt_;
  Float_t pterrll_;
  Float_t pterrlt_;
  Float_t etall_;
  Float_t etalt_;
  Float_t phill_;
  Float_t philt_;
  Float_t dilmass_;
  Float_t dilmasspf_;
  Float_t dilpt_;
  Int_t   flagll_;
  Int_t   flaglt_;

  Int_t   bptx_;       
  Int_t   bsc_;        
  Int_t   beamhalo_;   
  Int_t   goodvtx_;    
  Int_t   goodtrks_;   

  TH1F* hgenmet_all;
  TH1F* hgenmet_pass;
  TProfile* hresponse;  
  TH1F*   hgenps_pthat;
  TH1F*   hphotonpt;


  TH1F* hptz[5];
  TH1F* htcmet[4][4];
  TH1F* htcmetNew[4][4];
  TH1F* hpfmet[4][4];
  TH1F* hdilMass[4];

  TH1F* metPredicted;
  TH1F* metObserved;
  TH1F* metPredicted_sf;
  TH1F* metObserved_sf;
  TH1F* metPredicted_df;
  TH1F* metObserved_df;
  TH1F* metPredicted_ee;
  TH1F* metObserved_ee;
  TH1F* metPredicted_mm;
  TH1F* metObserved_mm;
  TH1F* metParPredicted;
  TH1F* metParObserved;
  TH1F* metPerpPredicted;
  TH1F* metPerpObserved;

  TH1F* metPredicted_njets[11];
  TH1F* metObserved_njets[11];

  TH1F* metTemplate[11][23];
  TH1F* metParTemplate[11][23];
  TH1F* metPerpTemplate[11][23];


};





#endif
