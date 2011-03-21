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
        vector<int> goodVertices();

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

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
        Int_t   nvtx_;
        Int_t   leptype_;
        Int_t   njets_;
        Float_t weight_;


        
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

	// latest-and-greatest tcmet stuff
	Float_t tcmetNew_;
	Float_t tcmetphiNew_;
	Float_t tcsumetNew_;

        ofstream ofile;
};





#endif
