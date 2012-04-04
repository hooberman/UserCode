#ifndef singleLeptonLooper_h
#define singleLeptonLooper_h

#include <vector>
#include <map>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
#include "../CORE/SimpleFakeRate.h" // will .h be ok? lets see.. 101007

//#include "../CORE/topmass/ttdilepsolve.h" REPLACETOPMASS

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
typedef map<unsigned int, unsigned int> m_uiui;

class  TChain;
class  TH1F;
class  TH2F;
class  TRandom3;
class  TTree;
struct metStruct;

class singleLeptonLooper
{
    public: 
        singleLeptonLooper();
        ~singleLeptonLooper() {}

        enum TrigEnum { e_highpt = 0, e_lowpt };  
        // e_highpt  :   high pt dilepton triggers, 20,10  
        // e_lowpt   :   dilepton-HT cross triggers, 10,5            
        enum JetTypeEnum { e_JPT = 0, e_calo , e_pfjet };
        // e_JPT     :   jpt jets
        // e_calo    :   l1 and l2 corrected calo jets
        // e_pfjet   :   corrected pfjets
        enum MetTypeEnum { e_tcmet = 0, e_muon, e_muonjes , e_pfmet };
        // e_tcmet   :   track corrected met
        // e_muon    :   calo met with muon corrections
        // e_muonjes :   calo met with muon and jet energy scale corrections
        // e_pfmet   :   particle-flow met
        enum ZVetoEnum   { e_standard = 0, e_allzveto, e_nozveto, e_selectz };
        // e_standard:   apply Z-veto to same-flavor pairs
        // e_allzveto:   apply Z-veto regardless of lepton flavor
        // e_nozveto :   no Z-veto
        // e_selectz :   select Z by requiring SF OS pair in Z mass window
        enum FREnum   { e_qcd = 0, e_wjets };
        // e_qcd     :   derive prediction for 2 fake leptons
        // e_wjets   :   derive prediction for 1 real and one fake lepton
	
        int  ScanChain(TChain *chain, char *prefix = "", float kFactor = 1., 
		       int prescale = 1., float lumi = 1.,
                       FREnum frmode  = e_wjets,
                       bool doFakeApp = false
                       );
        void BookHistos (char *prefix);
	void InitBaby();
	float dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx = 0 );
	
        // Set globals
        void set_susybaseline (bool  b)    { g_susybaseline = b; }
        void set_createTree   (bool  b)    { g_createTree   = b; }
        void set_useBitMask   (bool  b)    { g_useBitMask   = b; }
        void set_version      (char* v)    { g_version      = v; }
	void set_json         (char* v)    { g_json         = v; }        
        void set_trigger      (TrigEnum t) { g_trig         = t; } 

        // Baby ntuple methods
        void makeTree (char *prefix,bool doFakeApp, FREnum frmode );
	float stopPairCrossSection( float stopmass );
        void closeTree ();
	float trackIso( int thisPf , float coneR = 0.3 , float dz_thresh = 0.2 );
	std::vector<float> totalIso( int thisPf , float coneR = 0.3 , float dz_thresh = 0.2 );

	bool initialized;
	TH1D*   stop_xsec_hist;
	TFile*  stop_xsec_file;

    private:

        // Globals
        bool  g_susybaseline;
        bool  g_createTree;
        bool  g_useBitMask;
        char* g_version;
	char* g_json;      
	TrigEnum g_trig;
        TRandom3 *random3_;

	// MC truth lepton info
	Int_t   mcid1_;    
	Int_t   mcid2_;    
	Int_t	mcdecay1_; 
	Int_t   mcdecay2_; 
	Int_t	mcndec1_; 
	Int_t   mcndec2_; 
	Float_t mctaudpt1_;
	Float_t mctaudpt2_;
	Int_t   mctaudid1_;    
	Int_t   mctaudid2_;    
	Float_t mcdr1_;    
	Float_t mcdr2_;    

	// isolated track vars
	Float_t trkpt5_;
	Float_t mleptrk5_;
	Float_t trkreliso5_;
	Float_t trkpt10_;
	Float_t mleptrk10_;
	Float_t trkreliso10_;

	// btag variables
	Int_t   nbtagsssv_;     
	Int_t   nbtagstcl_;     
	Int_t   nbtagstcm_;     
	Int_t   nbtagscsvl_;    
	Int_t   nbtagscsvm_;    
	Int_t   nbtagscsvt_;    
 
	// pfjet counters
	Int_t   npfjets30_;
	Int_t   npfjets35_;
	Int_t   npfjets40_;
	Int_t   npfjets45_;
	Int_t   npfresjets30_;
	Int_t   npfresjets35_;
	Int_t   npfresjets40_;
	Int_t   npfresjets45_;

	// pfht vars
	Float_t htpf30_;
	Float_t htpf35_;
	Float_t htpf40_;
	Float_t htpf45_;
	Float_t htpfres30_;
	Float_t htpfres35_;
	Float_t htpfres40_;
	Float_t htpfres45_;

	// calojet counters
	Int_t   ncjets30_;
	Int_t   ncjets35_;
	Int_t   ncjets40_;
	Int_t   ncjets45_;
	Int_t   ncresjets30_;
	Int_t   ncresjets35_;
	Int_t   ncresjets40_;
	Int_t   ncresjets45_;

	// caloht vars
	Float_t htc30_;
	Float_t htc35_;
	Float_t htc40_;
	Float_t htc45_;
	Float_t htcres30_;
	Float_t htcres35_;
	Float_t htcres40_;
	Float_t htcres45_;

	// matched lepton vars
	Int_t   mlepid_;
	Int_t   mleppassid_;
	Int_t   mleppassiso_;
	Float_t mlepiso_;
	Float_t mlepdr_;

	// HLT variables
	Int_t   ldi_;
	Int_t   ltri_;
	Int_t   smu_;
	Int_t   dil_;

	// MC truth vars
	Int_t   npartons_;
	Float_t maxpartonpt_;
	Float_t ptt_;
	Float_t pttbar_;
	Float_t ptttbar_;
	Float_t mttbar_;
	Float_t etattbar_;
	Float_t mgcor_;
	Int_t   wflav_;

	// Type1 pfmet
	Float_t t1met10_;
	Float_t t1met20_;
	Float_t t1met30_;
	Float_t t1metres10_;
	Float_t t1metres20_;
	Float_t t1metres30_;
	Float_t t1met10phi_;
	Float_t t1met20phi_;
	Float_t t1met30phi_;
	Float_t t1metres10phi_;
	Float_t t1metres20phi_;
	Float_t t1metres30phi_;
	Float_t t1met10mt_;
	Float_t t1met20mt_;
	Float_t t1met30mt_;
	Float_t t1metres10mt_;
	Float_t t1metres20mt_;
	Float_t t1metres30mt_;

	// assorted p4's
	LorentzVector*  t_;   
	LorentzVector*  tbar_;   
	LorentzVector*  ttbar_;   
	LorentzVector*  mlep_;   
	LorentzVector*  mclep1_;   
	LorentzVector*  mclep2_;   
	LorentzVector*  mctaud1_;   
	LorentzVector*  mctaud2_;  
        LorentzVector*  lep1_;
        LorentzVector*  lep2_;
        LorentzVector*  dilep_;
        LorentzVector*  jet_; 

	// jet p4's
        LorentzVector*  cjet1_; 
        LorentzVector*  cjet2_; 
        LorentzVector*  cjet3_; 
        LorentzVector*  cjet4_; 
        LorentzVector*  pfjet1_; 
        LorentzVector*  pfjet2_; 
        LorentzVector*  pfjet3_; 
        LorentzVector*  pfjet4_; 
        LorentzVector*  cresjet1_; 
        LorentzVector*  cresjet2_; 
        LorentzVector*  cresjet3_; 
        LorentzVector*  cresjet4_; 
        LorentzVector*  pfresjet1_; 
        LorentzVector*  pfresjet2_; 
        LorentzVector*  pfresjet3_; 
        LorentzVector*  pfresjet4_; 

 	LorentzVector*  nonisoel_;   
 	LorentzVector*  nonisomu_;   

        // Baby ntuple variables
	Float_t mjj_;
	Float_t dphilm_;
	Float_t mG_;
	Float_t x_;
	Float_t mL_;
	Float_t ecalveto1_;
	Float_t ecalveto2_;
	Float_t hcalveto1_;
	Float_t hcalveto2_;
        TFile  *outFile;
        TTree  *outTree;
	Int_t   acc_2010_;
	Int_t   acc_highmet_;
	Int_t   acc_highht_;
	Float_t pthat_;
	Float_t qscale_;
        Float_t weight_;
        Float_t trgeff_;
        Float_t mutrigweight_;
        Float_t pfmetsig_;
        Float_t smeff_;
        Float_t k_;
        Float_t mllgen_;
        Float_t dphijm_;
        Float_t costhetaweight_;
	Int_t   njpt_;
	Float_t htjpt_;
        Int_t   hbhe_;
	Int_t   jetid_;
	Int_t   jetid30_;
        Int_t   mull_;
        Int_t   json_;
        Int_t   mult_;
        Int_t   mullgen_;
        Int_t   multgen_;
        Int_t   nlep_;
        Int_t   ngoodlep_;
        Int_t   ngoodel_;
        Int_t   nosel_;
        Int_t   ngoodmu_;
        Int_t   proc_;
        Int_t   leptype_;
        Int_t   ngenjets_;


        Int_t   njetsUp_;
        Int_t   npfjets25_;
        Int_t   njetsDown_;
	Float_t trkmet_;
	Float_t trkmetphi_;
	Float_t trkmetproj_;
	Float_t trkmet4_;
	Float_t trkmet4phi_;
	Float_t trkmet4proj_;
	Float_t trkmet8_;
	Float_t trkmet8phi_;
	Float_t trkmet8proj_;
	Float_t trkjetmet_;
	Float_t trkjetmetphi_;
	Float_t trkjetmetproj_;
        Float_t htUp_;
        Float_t htDown_;
        Int_t   nvtx_;
        Float_t dilmass_;
        Float_t topmass_;
        Float_t tcmet_;
        Float_t tcmet00_;
        Float_t tcmet10_;
        Float_t tcmet20_;
        Float_t tcmet30_;
        Float_t tcmet40_;
        Float_t tcmet50_;
        Float_t genmet_;
        Float_t gensumet_;
        Float_t genmetphi_;
        Float_t mucormet_;
        Float_t mucorjesmet_;
        Float_t pfmet_;
        Float_t pfsumet_;
        Float_t pfmetveto_;
        Float_t pfmetphi_;
        Float_t tcmet_35X_;
        Float_t tcmet_event_;
        Float_t tcmet_looper_;
        Float_t tcmetUp_;
        Float_t tcmetDown_;
        Float_t tcmetTest_;
        Float_t pfmetUp_;
        Float_t pfmetDown_;
        Float_t pfmetTest_;
        Float_t tcsumet_;
        Float_t tcmetphi_;
        Float_t mt2_;
        Float_t mt2j_;
        Float_t mt2jcore_;
        Float_t sumjetpt_;
        Float_t dileta_;
        Float_t dilpt_;
        Float_t dildphi_;
        Float_t vecjetpt_;
        Int_t   pass_;
        Int_t   passz_;
        Float_t m0_;
        Float_t m12_;
        Float_t ptl1_;
	Int_t   id1_;
	Int_t   id2_;
	Int_t   w1_;
	Int_t   w2_;
	Float_t iso1_;
	Float_t isont1_;
	Float_t iso2_;
	Float_t isont2_;
        Float_t ptl2_;
        Float_t etal1_;
        Float_t etal2_;
        Float_t phil1_;
        Float_t phil2_;
        Float_t meff_;
        Float_t mt_;
        char    dataset_[200];
        UInt_t  run_;
        UInt_t  lumi_;
        UInt_t  event_;
	Float_t y_;
	Float_t ht_;
	Float_t htoffset_;
	Float_t htuncor_;
	Int_t   njetsuncor_;
	Int_t   njetsoffset_;
	Float_t htgen_;
	Float_t htpf_;
	Float_t ptjetraw_;
	Float_t ptjet23_;
	Float_t ptjetF23_;
	Float_t ptjetO23_;
	Float_t cosphijz_;
	Int_t   ndavtx_;
	Int_t   nels_;
	Int_t   nmus_;
	Int_t   ntaus_;
	Int_t   nleps_;
	Float_t ndavtxweight_;
	Float_t etasc1_;
	Float_t etasc2_;
	Float_t emjet10_;
	Float_t emjet20_;

	Float_t ksusy_;
	Float_t ksusyup_;
	Float_t ksusydn_;
	Float_t xsecsusy_;
	Float_t xsecsusy2_;

	Float_t mbb_;

        // for fakeRates
        double getFRWeight(const int hypIdx, SimpleFakeRate *mufr, SimpleFakeRate *elfr, FREnum frmode, bool isData);

        // Lots and lots of histograms

        TH1F* h_PU_trkpt;
};

#endif
