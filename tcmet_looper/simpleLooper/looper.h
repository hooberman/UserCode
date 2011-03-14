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
        float dz(const unsigned int trkidx, const unsigned int vtxidx);
        float ZanettiMET(const unsigned int vtxIdx, const unsigned int hypIdx, float threshold);
        float HooberMET(const unsigned int vtxIdx, const unsigned int hypIdx, 
                        float dz_thresh, float pt_thresh, float etacut , bool usePFCandidatePt = false );

    private:
        
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

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
        Int_t   leptype_;
        Int_t   njets25_;

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

        ofstream ofile;
};





#endif
