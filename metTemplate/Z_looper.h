#ifndef Z_looper_h
#define Z_looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>


class TChain;

class Z_looper
{
 public:
  Z_looper() {};
  ~Z_looper() {
    delete babyFile_;
    delete babyTree_;
  };

  enum metAlgo   { e_makeTemplate = 0, e_photonSelection = 1, e_ZSelection = 2};


  void MakeBabyNtuple (const char *);
  void InitBabyNtuple ();
  void FillBabyNtuple ();
  void CloseBabyNtuple ();
  void ScanChain (TChain*, const char*, bool isData, bool calculateTCMET = false,
                  metAlgo algo = e_makeTemplate, int nEvents = -1, float kFactor = 1.);
  void bookHistos();
  bool isGoodTrack(int, bool usePV = false);
  float deltaPhi( float phi1 , float phi2);
  void fillUnderOverFlow(TH1F *h1, float value, float weight);
  void fillHistos(TH1F *h1[4],    float value, float weight, int myType);
  void fillHistos(TH1F *h1[4][4], float value, float weight, int myType, int nJetsIdx);
        
 private:
        
  metAlgo algo_;
        
  //ntuple, file
  TFile *babyFile_;
  TTree *babyTree_;
    
  //histos

  //triggers
  Int_t HLT_L1Jet6U_;
  Int_t HLT_L1Jet10U_;
  Int_t HLT_Jet15U_;
  Int_t HLT_Jet30U_;
  Int_t L1_SingleEG5_;
  Int_t HLT_Photon10_L1R_;
  Int_t HLT_Photon15_L1R_;
  Int_t HLT_Photon10_Cleaned_L1R_;
  Int_t HLT_Photon15_Cleaned_L1R_;
  Int_t HLT_Photon20_Cleaned_L1R_;
  Int_t HLT_Photon20_L1R_;

  // event stuff
  Int_t   run_;
  Int_t   lumi_;
  Int_t   event_;
  Int_t   leptype_;
  Int_t   nGoodVertex_;
  Float_t weight_;
  Float_t pthat_;

  // genmet stuff
  Float_t genmet_;
  Float_t genmetphi_;
  Float_t gensumet_;

  // pfmet stuff
  Float_t pfmet_;
  Float_t pfmetphi_;
  Float_t pfsumet_;

  // calomet stuff
  Float_t met_;
  Float_t metphi_;
  Float_t sumet_;

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
       
  // jet stuff
  Int_t   failjetid_;
  Int_t   nJets_;
  Int_t   nJets40_;
  Float_t sumJetPt_;
  Float_t sumJetPt10_;
  Float_t vecJetPt_;
  Int_t   nbtags_;

  //leading jet stuff
  Float_t jetmax_pt_;
  Float_t jetmax_dphimet_;

  //Z stuff
  Int_t   passz_;
  Int_t   passe_ll_ttbar_;
  Int_t   passe_ll_ttbarV1_;
  Int_t   passe_ll_ttbarV2_;
  Int_t   passe_ll_cand01_;
  Int_t   passm_ll_nomttbar_;
  Int_t   passm_ll_nomttbarV2_;
  Int_t   passm_ll_nom_;
  Int_t   passe_lt_ttbar_;
  Int_t   passe_lt_ttbarV1_;
  Int_t   passe_lt_ttbarV2_;
  Int_t   passe_lt_cand01_;
  Int_t   passm_lt_nomttbar_;
  Int_t   passm_lt_nomttbarV2_;
  Int_t   passm_lt_nom_;
  Int_t   pdgid_;
  Int_t   idll_;
  Int_t   idlt_;
  Float_t ptll_;
  Float_t ptlt_;
  Float_t etall_;
  Float_t etalt_;
  Float_t phill_;
  Float_t philt_;
  Float_t dilmass_;
  Float_t dilpt_;
  Int_t   flagll_;
  Int_t   flaglt_;

  Int_t   bptx_;       
  Int_t   bsc_;        
  Int_t   beamhalo_;   
  Int_t   goodvtx_;    
  Int_t   goodtrks_;   

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

  ofstream ofile_tcmet;
  ofstream ofile_events;

};





#endif
