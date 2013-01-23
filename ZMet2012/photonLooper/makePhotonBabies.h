#ifndef makephotonbabies_h
#define makephotonbabies_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

class TChain;
class FactorizedJetCorrector;

class makePhotonBabies
{
 public:
  makePhotonBabies() {};
  ~makePhotonBabies() {
    delete babyFile_;
    delete babyTree_;
  };

  void MakeBabyNtuple (const char *);
  void InitBabyNtuple ();
  void FillBabyNtuple ();
  void CloseBabyNtuple ();
  float dRGenJet ( LorentzVector p4, bool isData, float ptcut = 20.0 );
  int isGenQGLMatched ( LorentzVector p4, bool isData, float dR = 0.4 );
  int getJetIndex( LorentzVector thisJet , FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3 );
  void ScanChain (TChain*, const char*, bool isData, bool calculateTCMET = false,
                  int nEvents = -1, float kFactor = 1.);
  void  bookHistos();
  bool  isGoodTrack(int, bool usePV = false);
  float deltaPhi( float phi1 , float phi2);
  void  fillUnderOverFlow(TH1F *h1, float value, float weight);
  void  fillHistos(TH1F *h1[4],    float value, float weight, int myType);
  void  fillHistos(TH1F *h1[4][4], float value, float weight, int myType, int nJetsIdx);
  int   passThisHLTTrigger( char* hltname );
        
 private:
        
  //ntuple, file
  TFile *babyFile_;
  TTree *babyTree_;

  VofP4 pujets_;
  Int_t npujets_;
    
  //histos
  Float_t maxleppt_;
  Int_t   elveto_;

  //triggers
  Int_t hlt20_;
  Int_t hlt30_;
  Int_t hlt50_;
  Int_t hlt75_;
  Int_t hlt90_;
  Int_t hlt135_;
  Int_t hlt150_;
  Int_t hlt160_;

  Int_t hgg22_;
  Int_t hgg36_;
  Int_t hgg50_;
  Int_t hgg75_;
  Int_t hgg90_;

  Float_t rho_;
  Float_t ht30_;
  Float_t ht40_;
  Float_t jzb_;

  // event stuff
  char    dataset_[200];
  UInt_t  run_;
  UInt_t  lumi_;
  UInt_t  event_;
  Int_t   leptype_;
  Int_t   nGoodVertex_;
  Float_t weight_;
  Float_t pthat_;
  Int_t   failjetid_;
  Float_t maxemf_;
  
  // genmet stuff
  Float_t genmet_;
  Float_t genmetphi_;
  Float_t gensumet_;

  // pfmet stuff
  Float_t pfmet_;
  Float_t pfmetphi_;
  Float_t pfmett1_;
  Float_t pfmett1phi_;
  Float_t pfmett1new_;
  Float_t pfmett1newphi_;
  Float_t pfsumet_;
  Float_t pfmet_type1_pt30_;
  Float_t pfmet_type1_pt15_;

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
  Float_t tcmetNew_type1_pt30_;
  Float_t tcmetNew_type1_pt15_;

  // jet stuff
  Int_t   nJets_;
  Int_t   nJets40_;
  Int_t   nJets10_;
  Int_t   nJets15_;
  Int_t   nJets20_;
  Float_t ht_;
  Float_t ht10_;
  Float_t vecJetPt_;
  Int_t   nbl_;
  Int_t   nbm_;
  Int_t   nbt_;

  LorentzVector*  jet1_; 
  LorentzVector*  jet2_; 
  LorentzVector*  jet3_; 
  LorentzVector*  jet4_; 

  //leading jet stuff
  Float_t jetmax_pt_;
  Float_t jetmax_dphimet_;

  // photon stuff
  Int_t    nPhotons_;
  Float_t  etg_;
  Float_t  etag_;
  Float_t  phig_;
  Float_t  hoe_;
  Float_t  eciso_;
  Float_t  hciso_;
  Float_t  tkiso_;
  Float_t  swiss_;
  Float_t  seed_;
  Float_t  r4_;
  Float_t  s4_;
  Float_t  drel_;
  Float_t  jetidg_;

  //more photon stuff
  Int_t    photon_scidx_;        
  Int_t    photon_pixelseed_;        
  Float_t  photon_e15_;              
  Float_t  photon_e25max_;           
  Float_t  photon_e33_;              
  Float_t  photon_e55_;              
  Float_t  photon_ecalIso03_;        
  Float_t  photon_ecalIso04_;        
  Float_t  photon_hcalIso03_;        
  Float_t  photon_hcalIso04_;        
  Float_t  photon_ntkIsoHollow03_;   
  Float_t  photon_ntkIsoHollow04_;   
  Float_t  photon_ntkIsoSolid03_;    
  Float_t  photon_ntkIsoSolid04_;    
  Float_t  photon_sigmaEtaEta_;      
  Float_t  photon_sigmaIEtaIEta_;    
  Float_t  photon_sigmaIPhiIPhi_;    
  Float_t  photon_tkisoHollow03_;    
  Float_t  photon_tkisoHollow04_;    
  Float_t  photon_tkisoSolid03_;     
  Float_t  photon_tkisoSolid04_;     

  //photon-matched jet stuff
  Float_t jet_pt_;            
  Float_t calojet_pt_;            
  Float_t jet_dr_;            
  Float_t jet_eta_;            
  Float_t jet_energy_;            
  Float_t jet_chg_emfrac_;    
  Float_t jet_chg_hadfrac_;   
  Float_t jet_neu_emfrac_;    
  Float_t jet_neu_hadfrac_;   
  Int_t   jet_nchg_;              
  Int_t   jet_nmuon_;         
  Int_t   jet_nneu_;          
  Int_t   jet_pfjetid_;          
  Float_t jet_dphimet_;       
  Float_t jet_dpt_;           
  Float_t jet_drgen_;  

  Int_t   csc_;      
  Int_t   hbhe_;   
  Int_t   hbhenew_;   
  Int_t   hcallaser_;
  Int_t   ecaltp_;
  Int_t   trkfail_;
  Int_t   eebadsc_;
  
  Float_t jet1beta1_01_;
  Float_t jet2beta1_01_;
  Float_t jet3beta1_01_;
  Float_t jet4beta1_01_;

  Float_t jet1beta2_01_;
  Float_t jet2beta2_01_;
  Float_t jet3beta2_01_;
  Float_t jet4beta2_01_;

  Float_t jet1beta1_05_;
  Float_t jet2beta1_05_;
  Float_t jet3beta1_05_;
  Float_t jet4beta1_05_;

  Float_t jet1beta2_05_;
  Float_t jet2beta2_05_;
  Float_t jet3beta2_05_;
  Float_t jet4beta2_05_;

  Float_t jet1beta1_10_;
  Float_t jet2beta1_10_;
  Float_t jet3beta1_10_;
  Float_t jet4beta1_10_;

  Float_t jet1beta2_10_;
  Float_t jet2beta2_10_;
  Float_t jet3beta2_10_;
  Float_t jet4beta2_10_;

  Int_t   vtxidx_;

  Int_t   jet1flav_;
  Int_t   jet2flav_;
  Int_t   jet3flav_;
  Int_t   jet4flav_;

  Float_t jet1drgen_;
  Float_t jet2drgen_;
  Float_t jet3drgen_;
  Float_t jet4drgen_;

  TH1F* tcmetTemplate[3][7][4];
  TH1F* pfmetTemplate[3][7][4];
  TH1F* tcmetNewTemplate[3][7][4];

  TH1F* tcmetTemplate_combined[3][7];
  TH1F* pfmetTemplate_combined[3][7];
  TH1F* tcmetNewTemplate_combined[3][7];
  TH1F* tcmetNewType1Pt30Template_combined[3][7];
  TH1F* tcmetNewType1Pt15Template_combined[3][7];
  TH1F* pfmetType1Pt30Template_combined[3][7];
  TH1F* pfmetType1Pt15Template_combined[3][7];
  
  TH1F* hetg[5];
  TH1F* hdilMass[4];
  TH1F* htcmet[4][4];
  TH1F* hpfmet[4][4];
  TH1F* htcmetNew[4][4];

  TH1F*   hgenps_pthat;
  TH1F*   hphotonpt;
  
  ofstream ofile_tcmet;
  ofstream ofile_events;  
};





#endif
