#ifndef tcmetLooperTemplate_h
#define tcmetLooperTemplate_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>


class TChain;

class tcmetLooperTemplate
{
    public:
        tcmetLooperTemplate() {};
        ~tcmetLooperTemplate() {
            delete babyFile_;
            delete eventTree_;
            delete trackTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, const char*, bool isData, int nEvents = -1);
	bool isGoodZee ();
	bool isGoodZmm ();
	bool isMuon ( int index );
	bool isElectron ( int index );
	bool isGoodDilepton ();
        bool isTruthZee ();
        bool isTruthZmm ();
        bool jetVeto ();
        void bookHistos();
	bool isGoodTrack(int, bool usePV = false);
	void fillUnderOverFlow(TH1F *h1, float value, float weight = 1);

    private:
        
        //ntuple, file
        TFile *babyFile_;
        TTree *eventTree_;
        TTree *trackTree_;

        //track tree variables
        Float_t trk_pt_;
        Float_t trk_d0vtx_;
        Float_t trk_d0corr_;
        Int_t   trk_nhits_;
        Float_t trk_chi2_;
        Int_t   trk_ndf_;
        Float_t trk_pterr_;
        Float_t trk_phi_;
        Float_t trk_eta_;
        Int_t   trk_qual_;
        Int_t   trk_pass_;

        //histos
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
	Bool_t  hase_;
	
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
	Float_t tcmet_;
	Float_t tcmetphi_;
	Float_t tcsumet_;

	// raw tcmet stuff
	Float_t rawtcmet_;
	Float_t rawtcmetphi_;
	Float_t rawtcsumet_;

	// raw tcmet stuff
	Float_t raw35Xtcmet_;
	Float_t raw35Xtcmetphi_;
	Float_t raw35Xtcsumet_;

	Float_t ectcmet_;
	Float_t ectcmetphi_;
	Float_t ectcsumet_;

	// electron stuff
	Float_t eldr_;
	Float_t elpt_;
	Float_t eleta_;

        // vars for dupTree
        Float_t dupdr_;
        Float_t dupdphi_;
        Float_t dupdcotth_;

        ofstream ofile;
};





#endif
