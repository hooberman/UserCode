#ifndef makeVictorTemplates_h
#define makeVictorTemplates_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>


class TChain;

class makeVictorTemplates
{
 public:
  makeVictorTemplates() {};
  ~makeVictorTemplates() {
  };
  
  void setBranches (TTree *tree);
  void ScanChain (TChain*);
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

  Float_t   jet15U_HLTPrescale_;
  Float_t   jet15U_L1Prescale_;
  Int_t     firedJet15U_;
  
  Float_t   jet30U_HLTPrescale_;
  Float_t   jet30U_L1Prescale_;
  Int_t     firedJet30U_;
  
  Float_t   jet50U_HLTPrescale_;
  Float_t   jet50U_L1Prescale_;
  Int_t     firedJet50U_;
  
  Float_t   jet100U_HLTPrescale_;
  Float_t   jet100U_L1Prescale_;
  Int_t     firedJet100U_;

  TH1F* hleadJetPt15;
  TH1F* hleadJetPt30;
  TH1F* hleadJetPt50;
  TH1F* hleadJetPt100;

  TH1F* tcmetTemplate[3][7][4];
  TH1F* pfmetTemplate[3][7][4];
  
  TH1F* tcmetTemplate_njets_ht_nvtx[3][7][3];
  TH1F* pfmetTemplate_njets_ht_nvtx[3][7][3];
  
  TH1F* tcmetTemplate_combined[3][7];
  TH1F* pfmetTemplate_combined[3][7];

  TH1F* tcmetTemplate_qcd[4][3][7];
  TH1F* pfmetTemplate_qcd[4][3][7];
  
};





#endif
