#ifndef looper_h
#define looper_h

#include <vector>
#include <map>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

//#include "../CORE/topmass/ttdilepsolve.h" REPLACETOPMASS

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

class  TChain;
class  TH1F;
class  TH2F;
class  TRandom3;
class  TTree;
struct metStruct;

class looper
{
    public: 
        looper();
        ~looper() {}
	
        int  ScanChain(TChain *chain, char *prefix = "");
        void BookHistos (char *prefix);
	void InitBaby();
	float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx = 0 );
	float getMinDR( int type1, int type2, std::vector<int> hltid, std::vector<LorentzVector> hltobj );
	int getJetIndex(LorentzVector jet);
	
        // Set globals
        void set_susybaseline (bool  b)    { g_susybaseline = b; }
        void set_createTree   (bool  b)    { g_createTree   = b; }
        void set_useBitMask   (bool  b)    { g_useBitMask   = b; }
        void set_version      (char* v)    { g_version      = v; }
	void set_json         (char* v)    { g_json         = v; }        

        // Baby ntuple methods
        void makeTree (char *prefix );
        void closeTree ();

    private:

	TTree* outTree;
	TFile* outFile;

        // Globals
        bool  g_susybaseline;
        bool  g_createTree;
        bool  g_useBitMask;
        char* g_version;
	char* g_json;      
        TRandom3 *random3_;

        LorentzVector*  eledijet_hltele_; 
        LorentzVector*  eletrijet_hltele_; 
        LorentzVector*  mudijet_hltmu_; 
        LorentzVector*  mutrijet_hltmu_; 
        LorentzVector*  cjet1_; 
        LorentzVector*  cjet2_; 
        LorentzVector*  cjet3_; 
        LorentzVector*  cjet4_; 
        LorentzVector*  pjet1_; 
        LorentzVector*  pjet2_; 
        LorentzVector*  pjet3_; 
        LorentzVector*  pjet4_; 
	Float_t pjet1_L1Fast_;
	Float_t pjet2_L1Fast_;
	Float_t pjet3_L1Fast_;
	Float_t pjet4_L1Fast_;
	Float_t pjet1_L2L3_;
	Float_t pjet2_L2L3_;
	Float_t pjet3_L2L3_;
	Float_t pjet4_L2L3_;
        char    dataset_[200];
        UInt_t  run_;
        UInt_t  lumi_;
        UInt_t  event_;
	Int_t   njets_;
	Int_t   ncjets_;
	Float_t ht_;
	Float_t htc_;
	Float_t pfmet_;
	Float_t pfmetphi_;
	Float_t pfsumet_;
	Int_t   ngoodlep_;
	Int_t   ngoodmu_;
	Int_t   ngoodel_;
	Int_t   ndavtx_;
	Int_t   nvtx_;
	Int_t   eledijet_hg_;
	Int_t   eledijetmht15_;
	Int_t   eledijetmht25_;
	Int_t   eledijet_;
	Int_t   eletrijet_;
	Int_t   elequadjet_;
	Int_t   eledijet_n82_;
	Int_t   eletrijet_n82_;
	Int_t   eledijet_n85_;
	Int_t   eletrijet_n85_;
	Int_t   mudijet_n83_;
	Int_t   mutrijet_n83_;
	Int_t   mudijet_n85_;
	Int_t   mutrijet_n85_;
	Float_t eledijet_trigmindr_ejet_;
	Float_t eletrijet_trigmindr_ejet_;
	Float_t mudijet_trigmindr_mujet_;
	Float_t mutrijet_trigmindr_mujet_;
	Int_t   mudijet_;
	Int_t   mudijetmht15_;
	Int_t   mudijetmht25_;
	Int_t   mutrijet_;
	Int_t   muquadjet_;
	Float_t elptmatch_;
	Float_t eledijet_trigdr_pjet1_;
	Float_t eledijet_trigdr_pjet2_;
	Float_t eledijet_trigdr_pjet3_;
	Float_t eledijet_trigdr_pjet4_;
	Float_t eletrijet_trigdr_pjet1_;
	Float_t eletrijet_trigdr_pjet2_;
	Float_t eletrijet_trigdr_pjet3_;
	Float_t eletrijet_trigdr_pjet4_;
	Float_t mudijet_trigdr_pjet1_;
	Float_t mudijet_trigdr_pjet2_;
	Float_t mudijet_trigdr_pjet3_;
	Float_t mudijet_trigdr_pjet4_;
	Float_t mutrijet_trigdr_pjet1_;
	Float_t mutrijet_trigdr_pjet2_;
	Float_t mutrijet_trigdr_pjet3_;
	Float_t mutrijet_trigdr_pjet4_;
};

#endif
