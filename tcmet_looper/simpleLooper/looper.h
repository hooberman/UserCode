#ifndef looper_h
#define looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "TLorentzVector.h"

using namespace std;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector; 


class TChain;

class looper
{
    public:
        looper() {};
        ~looper() {
            delete babyFile_;
            delete eventTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, const char*, bool isData, int nEvents = -1);
        void bookHistos();
	void fillUnderOverFlow(TH1F *h1, float value, float weight = 1);
        void printEvent(  ostream& ostr = cout );
        float deltaPhi( float phi1 , float phi2 );
        float dz_trk_vtx(const unsigned int trkidx, const unsigned int vtxidx);
        std::pair<float,float>  ZanettiMET(const unsigned int vtxIdx, const unsigned int hypIdx, float threshold);
        //float HooberMET(const unsigned int vtxIdx, const unsigned int hypIdx, 
        //                float dz_thresh, float pt_thresh, float etacut , bool usePFCandidatePt = false );
        std::pair<float,float> PFCandidateMET(const unsigned int vtxIdx, const unsigned int hypIdx, vector<int> goodjets, 
                                              float dz_thresh, float pt_thresh, float etacut , bool usePFCandidatePt , bool correctJets );

        bool jetFromSignalPV( int ijet , int vtxIdx , int beta_exponent = 2 );
        float beta_jet_vtx(   int ijet , int vtxIdx , int beta_exponent = 2 );
        std::pair<float,float> pfmetByHand( float ptcut , float etacut );
        vector<int> goodVertices();
        vector<int> goodDAVertices();
  
    private:

        LorentzVector*  lep1_;
        LorentzVector*  lep2_;
        LorentzVector*  dilep_;
        LorentzVector*  jet_;
        
        //ntuple, file
        TFile *babyFile_;
        TTree *eventTree_;

        //histos
	TH1F* htcmetNew;
	TH1F* htcmetNewPFC;
	TH1F* hpfmet;

        TH1F* hmumet;
        TH1F* hmujesmet;
        TH1F* htcmet;
        TH1F* hrawtcmet;

        TH1F* hdmumet;
        TH1F* hdmujesmet;
        TH1F* hdtcmet;
        TH1F* hdrawtcmet;

        TH1F* hdtcmet_mumet;
        TH1F* hdtcmet_mujesmet;

        TH1F* hneutralpt_0jet;
        TH1F* hneutralpt_1jet;

        TH1F* hmismatch_jetpt_tot;
        TH1F* hmismatch_jetpt_pass;
        TH1F* hmatch_jetpt_tot;
        TH1F* hmatch_jetpt_pass;     
        
        TH1F* hmismatch_nvtx_tot;    
        TH1F* hmismatch_nvtx_pass;   
        TH1F* hmatch_nvtx_tot;       
        TH1F* hmatch_nvtx_pass;      

        Float_t jetpt_;
        Float_t jeteta_;
        Float_t jetphi_;
        Float_t jetbeta_;
        Int_t   jetpv_;

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
        Int_t   nvtx_;
        Int_t   ndavtx_;
        Int_t   vtxIdx_;
        Int_t   leptype_;
        Int_t   njets25_;
        Int_t   njets30_;
        Float_t dilmass_;
        Float_t weight_;
        Float_t vtxweight_;

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

	// tcmet stuff
	Float_t tcmet_;
	Float_t tcmetphi_;
	Float_t tcsumet_;

	// pfcandmet stuff
	Float_t pfcandmet_;
	Float_t pfcandmetphi_;
	Float_t pfcandsumet_;

	// latest-and-greatest tcmet stuff
	Float_t tcmetNew_;
	Float_t tcmetphiNew_;
	Float_t tcsumetNew_;

        Float_t hmet_;
        Float_t hmetpfnotrks_;
        Float_t hmetpf_;
        Float_t hmetpf0_;
        Float_t hmetpf1_;
        Float_t hmetpf2_;
        Float_t hmetpf3_;
        Float_t hmetpf4_;
        Float_t hmetpf5_;
        Float_t hmetpf6_;
        Float_t hmetpf7_;
        Float_t hmetpf8_;
        Float_t hmetpf9_;
        Float_t hmetpf10_;
        Float_t zmet_;
        Float_t jetzmet_;
        Float_t jetzmetnotrks_;
        Float_t jetzmet4_;
        Float_t jetzmet8_;
        Float_t pfmet3_;

        Float_t hmetphi_;
        Float_t hmetphipfnotrks_;
        Float_t hmetphipf_;
        Float_t hmetphipf0_;
        Float_t hmetphipf1_;
        Float_t hmetphipf2_;
        Float_t hmetphipf3_;
        Float_t hmetphipf4_;
        Float_t hmetphipf5_;
        Float_t hmetphipf6_;
        Float_t hmetphipf7_;
        Float_t hmetphipf8_;
        Float_t hmetphipf9_;
        Float_t hmetphipf10_;
        Float_t zmetphi_;
        Float_t jetzmetphinotrks_;
        Float_t jetzmetphi_;
        Float_t jetzmetphi4_;
        Float_t jetzmetphi8_;
        Float_t pfmetphi3_;

        Float_t hmetproj_;
        Float_t hmetpfnotrksproj_;
        Float_t hmetpfproj_;
        Float_t hmetpf0proj_;
        Float_t hmetpf1proj_;
        Float_t hmetpf2proj_;
        Float_t hmetpf3proj_;
        Float_t hmetpf4proj_;
        Float_t hmetpf5proj_;
        Float_t hmetpf6proj_;
        Float_t hmetpf7proj_;
        Float_t hmetpf8proj_;
        Float_t hmetpf9proj_;
        Float_t hmetpf10proj_;
        Float_t zmetproj_;
        Float_t jetzmetnotrksproj_;
        Float_t jetzmetproj_;
        Float_t jetzmet4proj_;
        Float_t jetzmet8proj_;
        Float_t pfmet3proj_;

        ofstream ofile;
};





#endif
