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
            delete babyTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, const char*, bool isData, int nEvents = -1);
	bool isGoodZee ();
	bool isGoodZmm ();
	bool isGoodDilepton ();
        bool isTruthZee ();
        bool isTruthZmm ();
        bool jetVeto ();
        void bookHistos();
	bool isGoodTrack(int, bool usePV = false);


    private:
        
        //ntuple, file
        TFile *babyFile_;
        TTree *babyTree_;
    
        //histos
        TH1F* hmumet;
        TH1F* hmujesmet;
        TH1F* htcmet;
        TH1F* hrawtcmet;

        TH1F* hdmumet;
        TH1F* hdmujesmet;
        TH1F* hdtcmet;
        TH1F* hdrawtcmet;


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
