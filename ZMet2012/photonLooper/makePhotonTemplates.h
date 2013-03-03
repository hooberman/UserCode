#ifndef makePhotonTemplates_h
#define makePhotonTemplates_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>


class TChain;

class makePhotonTemplates
{
 public:
  makePhotonTemplates() {};
  ~makePhotonTemplates() {
  };
  
  void setBranches (TTree *tree);
  void ScanChain (TChain*,char*);
  void bookHistos();
  void fillUnderOverFlow(TH1F *h1, float value, float weight);
  
 private:

  Float_t   tcmet_;
  Float_t   pfmet_;
  Int_t     nJets_;
  Float_t   sumJetPt_;
  Int_t     nvtx_;
  Float_t   maxjetpt_;
  Float_t   prescale_;
  Float_t   jet_pt_;
  Float_t   etg_;
  Float_t   etag_;
  Int_t     pfjetid_;

  Int_t     hlt20_;
  Int_t     hlt30_;
  Int_t     hlt50_;
  Int_t     hlt75_;
  Int_t     hlt125_;

  TH1F* hphotonPt20;
  TH1F* hphotonPt30;
  TH1F* hphotonPt50;
  TH1F* hphotonPt70;

  TH1F* tcmetTemplate[3][7][4];
  TH1F* pfmetTemplate[3][7][4];
  
  TH1F* tcmetTemplate_njets_ht_nvtx[3][7][3];
  TH1F* pfmetTemplate_njets_ht_nvtx[3][7][3];
  
  TH1F* tcmetTemplate_combined[3][7];
  TH1F* pfmetTemplate_combined[3][7];

  TH1F* tcmetTemplate_photon[4][3][7];
  TH1F* pfmetTemplate_photon[4][3][7];
  
};





#endif
