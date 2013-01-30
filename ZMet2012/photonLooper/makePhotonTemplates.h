#ifndef makePhotonTemplates_h
#define makePhotonTemplates_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>
//#include "Math/LorentzVector.h"

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class TChain;

class makePhotonTemplates
{
 public:
  makePhotonTemplates() {};
  ~makePhotonTemplates() {
  };
  
  void setBranches (TTree *tree);
  void ScanChain (TChain*,char*,char*);
  void bookHistos();
  void fillUnderOverFlow(TH1F *h1, float value, float weight);
  
 private:

  Float_t   tcmet_;
  Float_t   pfmet_;
  Float_t   pfmetphi_;
  Int_t     nJets_;
  Int_t     run_;
  Int_t     nJets40_;
  Int_t     nbl_;
  Int_t     nbm_;
  Int_t     nbt_;
  Float_t   ht_;
  Float_t   ht40_;
  Int_t     nvtx_;
  Float_t   maxjetpt_;
  Float_t   prescale_;
  Float_t   jet_pt_;
  Float_t   calojet_pt_;
  Float_t   etg_;
  Float_t   etag_;
  Float_t   phig_;
  Int_t     pfjetid_;
  Float_t   hoe_;
  Int_t     photon_pixelseed_;
  Float_t   maxleppt_;
  Int_t     elveto_;
  Float_t   jetneutralemfrac_;
  Float_t   pfmett1new_;
  Float_t   pfmett1_;

  Int_t csc_;
  Int_t hbhe_;
  Int_t hcallaser_;
  Int_t ecaltp_;
  Int_t trkfail_;
  Int_t eebadsc_;
  Int_t hbhenew_;

  Int_t hlt20_;
  Int_t hlt30_;
  Int_t hlt50_;
  Int_t hlt75_;
  Int_t hlt90_;
  Int_t hlt125_;

  Int_t nbcsvl_;
  Int_t nbcsvm_;
  Int_t nbcsvt_;

  Int_t hgg22_;
  Int_t hgg36_;
  Int_t hgg50_;
  Int_t hgg75_;
  Int_t hgg90_;

  TH1F* hphotonPt20;
  TH1F* hphotonPt30;
  TH1F* hphotonPt50;
  TH1F* hphotonPt70;
  TH1F* hphotonPt90;
  TH1F* hphotonAll;

  TH1F* hphotonPt20_exc;
  TH1F* hphotonPt30_exc;
  TH1F* hphotonPt50_exc;
  TH1F* hphotonPt70_exc;
  TH1F* hphotonPt90_exc;

  TH1F* hnvtxPt20;
  TH1F* hnvtxPt30;
  TH1F* hnvtxPt50;
  TH1F* hnvtxPt70;
  TH1F* hnvtxPt90;
  TH1F* hnvtxAll;

  TH1F* tcmetTemplate[3][7][4];
  TH1F* pfmetTemplate[3][7][4];
  
  TH1F* tcmetTemplate_njets_ht_nvtx[3][7][3];
  TH1F* pfmetTemplate_njets_ht_nvtx[3][7][3];
  
  TH1F* tcmetTemplate_combined[3][7];
  TH1F* pfmetTemplate_combined[3][7];

  TH1F* tcmetTemplate_photon[5][3][7];
  TH1F* pfmetTemplate_photon[5][3][7];
  TH1F* t1pfmetTemplate_photon[5][3][7];
  TH1F* t1newpfmetTemplate_photon[5][3][7];

  // LorentzVector  jet1_;
  // LorentzVector  jet2_;
  
  // LorentzVector *jet1Ptr_; 
  // LorentzVector *jet2Ptr_; 
};





#endif

