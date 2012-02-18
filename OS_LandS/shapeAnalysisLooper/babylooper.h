#ifndef babylooper_h
#define babylooper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TProfile.h"
#include <vector>
#include <fstream>


class TChain;

class babylooper
{
    public:
        babylooper() {};
        ~babylooper() {
        };

        void setBranches (TTree *tree);
        void ScanChain (TChain*, const char*, const char*, bool isData);

    private:

        // branches in baby
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
        Int_t   leptype_;
        Int_t   passz_;
	Int_t   njets_;
	Int_t   njetsUp_;
	Int_t   njetsDown_;
        Float_t m0_;
        Float_t m12_;
	Float_t dilpt_;
	Float_t pfmet_;
	Float_t pfmetUp_;
	Float_t pfmetDown_;
	Float_t ht_;
	Float_t htUp_;
	Float_t htDown_;
	Float_t weight_;
	Float_t ndavtxweight_;
	Float_t trgeff_;
	Float_t lepscale_;
};





#endif
