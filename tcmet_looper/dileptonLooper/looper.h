#ifndef looper_h
#define looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>
#include <iostream>

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
        void fillHistos(TH1F *h1[4],    float value, float weight, int myType);
        void fillHistos(TH1F *h1[4][4], float value, float weight, int myType, int nJetsIdx);
        bool isGoodTrack( int index );
        bool isMuon( int index);
        bool isElectron( int index);


    private:
        
        //ntuple, file
        TFile *babyFile_;
        TTree *eventTree_;

        //histos
        TH1F* htcmet[4][4];
        TH1F* htcmetNew_calo[4][4];
        TH1F* htcmetNew_pfc[4][4];
        TH1F* hpfmet[4][4];
        TH1F* hnjets[4];
        TH1F* hjetpt[4];
        TH1F* hdilMass[4];
        TH1F* hdphijetmet_metlt20[4][4];
        TH1F* hdphijetmet_metgt30[4][4];
        TH1F* heleta_metgt30;
        TH1F* heleta_metlt20;
        TH1F* hmueta_metgt30;
        TH1F* hmueta_metlt20;
        TH1F* htrkptnearel;

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
	Bool_t  hase_;

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

        // new tcmet stuff
        Float_t tcmetNew_calo_;
        Float_t tcmetNew_pfc_;

        // dilepton stuff
        Float_t dilMass_;

        //jet-met stuff
        Float_t dphijetmet_;

        Int_t leptype_;
        Int_t njets_;

        ofstream ofile;
};





#endif
