#ifndef looper_h
#define looper_h

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

class looper
{
    public: 
        looper();
        ~looper() {}

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
                       JetTypeEnum jetType = e_JPT, 
                       MetTypeEnum metType = e_tcmet,
                       ZVetoEnum zveto = e_standard,
                       FREnum frmode  = e_wjets,
                       bool doFakeApp = false,
                       bool calculateTCMET = false
                       );
        void BookHistos (char *prefix);

	void InitBaby();
	
        // Set globals
        void set_susybaseline (bool  b)    { g_susybaseline = b; }
        void set_createTree   (bool  b)    { g_createTree   = b; }
        void set_useBitMask   (bool  b)    { g_useBitMask   = b; }
        void set_version      (char* v)    { g_version      = v; }
	void set_json         (char* v)    { g_json         = v; }        
        void set_trigger      (TrigEnum t) { g_trig         = t; } 

        // Baby ntuple methods
        void makeTree (char *prefix,bool doFakeApp, FREnum frmode );
        void closeTree ();
	bool objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt);

	bool initialized;

    private:

        // Globals
        bool  g_susybaseline;
        bool  g_createTree;
        bool  g_useBitMask;
        char* g_version;
	char* g_json;      
	TrigEnum g_trig;
        TRandom3 *random3_;

	Int_t   mcid1_;    
	Int_t   mcid2_;    
	Int_t	mcdecay1_; 
	Int_t   mcdecay2_; 
	Float_t mcdr1_;    
	Float_t mcdr2_;    

	Float_t trkpt_;
	Float_t trkreliso_;

	Int_t   hlt1_;
	Int_t   hlt2_;
	Int_t   mlepid_;
	Int_t   mleppassid_;
	Int_t   mleppassiso_;
	Float_t mlepiso_;
	Float_t mlepdr_;

	Float_t iso1_;
	Float_t iso2_;
	Int_t   passid1_;
	Int_t   passid2_;

	/* vector<int> pingruns; */
	/* vector<int> pinglumis; */
	/* vector<int> pingevents; */
	
	/* ifstream pingfile; */

        LorentzVector*  lep1trk_;
        LorentzVector*  lep1glb_;
        LorentzVector*  lep1sta_;
        LorentzVector*  lep2trk_;
        LorentzVector*  lep2glb_;
        LorentzVector*  lep2sta_;

        LorentzVector*  lep1_;
        LorentzVector*  lep2_;
        LorentzVector*  dilep_;
        LorentzVector*  jet1_; 
        LorentzVector*  jet2_; 
        LorentzVector*  jet3_; 
        LorentzVector*  jet4_; 

        // Baby ntuple variables
        TFile  *outFile;
        TTree  *outTree;
	Int_t   acc_2010_;
	Int_t   acc_highmet_;
	Int_t   acc_highht_;
	Float_t pingmass_;
	Float_t jetcor1_;
	Float_t jetcor2_;
	Float_t cor1_;
	Float_t cor2_;
	Float_t mmjj_;
	Float_t mmjjuncor_;
	Float_t mmjjc_;
	Float_t metnew_;
	Float_t mmjjdef_;
	Float_t mmjjtrk_;
	Float_t mmjjtrkdef_;
	Float_t mmjjtrkuncor_;
	Float_t mmjjglb_;
	Float_t mmjjsta_;
	Float_t mjjj_;
	Float_t pthat_;
	Float_t qscale_;
        Float_t weight_;
        Float_t trgeff_;
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
        Int_t   ngoodmu_;
        Int_t   proc_;
        Int_t   leptype_;
        Int_t   njets_;
        Int_t   njetsdef_;
        Int_t   njetsplain_;
        Int_t   ngenjets_;
        Int_t   npfjets_;
        Int_t   npfjets40_;
        Int_t   npfjetspv_;
        Int_t   njetsUp_;
        Int_t   npfjets25_;
        Int_t   njetsDown_;
        Int_t   nvtx_;
        Int_t   nbtags17_;
        Int_t   nbtags20_;
        Int_t   nbtags20_24_;
        Int_t   nbtags20c_;
        Int_t   nbtags20_24c_;
        Int_t   nbtags33_;
        Float_t dilmass_;
        Float_t topmass_;
        Float_t tcmet_;
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
	Int_t   diltrig_;
        Float_t ptl2_;
        Float_t ptj1_;
        Float_t ptj2_;
        Float_t etal1_;
        Float_t etal2_;
        Float_t phil1_;
        Float_t phil2_;
        Float_t meff_;
        Float_t mt_;
        char    dataset_[200];
        ULong64_t  run_;
        ULong64_t  lumi_;
        ULong64_t  event_;
	Float_t y_;
	Float_t ht_;
	Float_t htoffset_;
	Float_t htuncor_;
	Int_t   njetsuncor_;
	Int_t   njetsoffset_;
	Float_t htgen_;
	Float_t htpf_;
	Float_t htpf40_;
	Float_t htpf25_;
	Float_t htpfpv_;
	Float_t ptjetraw_;
	Float_t ptjet23_;
	Float_t ptjetF23_;
	Float_t ptjetO23_;
	Float_t cosphijz_;
	Int_t   njets15_;
	Int_t   ndavtx_;
	Int_t   nels_;
	Int_t   nmus_;
	Int_t   ntaus_;
	Int_t   nleps_;
	Float_t ndavtxweight_;
	Float_t etasc1_;
	Float_t etasc2_;
	Int_t   mtrg_;
	Int_t   eetrg_;
	Int_t   mmtrg_;
	Int_t   emtrg_;

	Float_t ksusy_;
	Float_t ksusyup_;
	Float_t ksusydn_;
	Float_t xsecsusy_;
	Float_t xsecsusy2_;

        // for fakeRates
        double getFRWeight(const int hypIdx, SimpleFakeRate *mufr, SimpleFakeRate *elfr, FREnum frmode, bool isData);

        // Lots and lots of histograms

        //Z histos
	TH1F* hping;
        TH1F* hdilMass_Z[4][4];
        TH1F* htcmet_event_Z[4][4];
        TH1F* htcmet_looper_Z[4][4];
        TH1F* hpfmet_Z[4][4];
        TH1F* hmucormet_Z[4][4];
        TH1F* hmucorjesmet_Z[4][4];

        TH1F* hyield;
        TH1F* hyield_weight;
        TH1F* hyieldsig;
        TH1F* hyield_unweighted;
        TH1F* hyieldsig_unweighted;

        TH2F*     hdtcmetevent_genmet[4][4];
        TProfile* tdtcmetevent_genmet[4][4];
        TH2F*     hdtcmetlooper_genmet[4][4];
        TProfile* tdtcmetlooper_genmet[4][4];
        TH2F*     hdpfmet_genmet[4][4];
        TProfile* tdpfmet_genmet[4][4];
        TH2F*     hdmucormet_genmet[4][4];
        TProfile* tdmucormet_genmet[4][4];
        TH2F*     hdmucorjesmet_genmet[4][4];
        TProfile* tdmucorjesmet_genmet[4][4];

        TH1F* hmt2j_signal[4][4];
        TH1F* hmt2j_control[4][4];
        TH1F* hmt2j_all[4][4];
        TH2F* hmet_dilpt_signal[4][4];
        TH2F* hmet_dilpt_control[4][4];
        TH2F* hmet_dilpt_all[4][4];
              
        TH1F* hmt[4][4];
        TH1F* hetaz[4][4];
        TProfile* htcsumet_tcmet_prof[4][4]; 
        TProfile* hsumJetPt_tcmetsqrtsumet_prof[4][4]; 
        TProfile* hgensumet_genmet_prof[4][4]; 
        //TProfile* hsumJetPt_tcmetpowtcsumet_prof[4][4][101]; 
        //TH2F* hsumJetPt_tcmetpowtcsumet_th2[4][4][101];
        //TH2F*     hsumJetPt_tcmetpowtcsumet_th2[101];
        //TProfile* hsumJetPt_tcmetpowtcsumet_prof[101];

	TH2F* msugra_highmet;
	TH2F* msugra_highht;
	TH2F* msugra_all;

	TH2F* msugra_highmet_kup;
	TH2F* msugra_highmet_kdn;
	TH2F* msugra_highmet_jup;
	TH2F* msugra_highmet_jdn;

	TH2F* msugra_highht_kup;
	TH2F* msugra_highht_kdn;
	TH2F* msugra_highht_jup;
	TH2F* msugra_highht_jdn;

        TH2F* hsumJetPt_tcmet[4][4];
        TH2F* hsumJetPt_tcmetsqrtsumet[4][4]; 
        TH2F* hsumJetPt_tcmetsumet[4][4]; 
        TH2F* hetaZ_tcmet[4][4]; 
        TH2F* hetaZ_tcmetsqrtsumet[4][4]; 
        TH2F* hetaZ_tcmetsumet[4][4]; 

        TH2F* hdilMass_tcmet[4][4];
        TH1F* hmt2core[4][4];               // MT2 from CORE
        TH1F* hmt2jcore[4][4];               // MT2J from CORE
        TH1F* hmt2j[4][4];                   // potentially custom MT2J
        TH1F* hsumJetPt[4][4];               // scalar sum jet Et
        TH1F* hsumJetPt_cut[4][4];               // scalar sum jet Et
        TH2F* hDtcmetgenmetVsumJetPt[4][4];
        TH2F* hDmetmuonjesgenmetVsumJetPt[4][4];
        TH1F* hmeffJet[4][4];                // scalar sum jet pt + scalar sum dil pt + (tc)met

        TH2F* habcd[4][4];                   
        TProfile* habcd_tprof[4][4];                   
        TH2F* habcd_nopresel[4][4];                   
        TProfile* habcd_tprof_nopresel[4][4];                   
        TH1F* hsumJptPt[4][4];               // scalar sum JPT jet Et
        TH1F* hsumJptPt_cut[4][4];               // scalar sum JPT jet Et
        TH1F* hmeffJPT[4][4];                // scalar sum JPT jet pt + scalar sum dil pt + (tc)met

        TH1F* hsumHypPt[4][4];               // scalar sum JPT jet Et
        TH1F* hmeffHyp[4][4];                // scalar sum JPT jet pt + scalar sum dil pt + (tc)met

        TH1F* hnJet[4];                      // Njet distributions
        TH1F* hnJpt[4];                      // Njpt distributions
        TH1F* hnBtagJpt[4];                  // N btag jpts
        TH1F* hnHypJet[4];                   // Hyp Njet distributions
        TH1F* helePt[4][4];                  // electron Pt
        TH1F* hmuPt[4][4];                   // muon Pt
        TH1F* hminLepPt[4][4];               // minimum lepton Pt
        TH1F* hmaxLepPt[4][4];               // maximum lepton Pt
        TH1F* hminLepEta[4][4];              // minimum lepton eta
        TH1F* hmaxLepEta[4][4];              // maximum lepton eta
        TH1F* helePhi[4][4];                 // electron phi
        TH1F* hmuPhi[4][4];                  // muon phi
        TH1F* hdphiLep[4][4];                // delta phi between leptons
        TH1F* hdrLep[4][4];                  // dR between leptons
        TH1F* hdrJ1J2[4][4];                 // dR between 2 leading jets
        TH1F* heleEta[4][4];                 // electron eta
        TH1F* hmuEta[4][4];                  // muon eta
        TH1F* htopMass[4][4];                // top mass estimate for 2 highest pt jets
        TH1F* htopMassAllComb[4][4];         // top mass estimate for all jets
        TH1F* hdilMass[4][4];                // dilepton mass
        TH1F* hdilMass_cut[4][4];            // dilepton mass
        TH1F* hdilPt[4][4];                  // dilepton Pt
        TH1F* hdilPt_zveto[4][4];            // dilepton Pt with z-veto applied
        TH1F* hdilPtSmeared[4][4];           // dilepton Pt with Gaussian smearing

        TH1F* hgenmet[4][4];                 // MET corrected for muons and JES
        TH1F* hgenmetPhi[4][4];              // MET corrected for muons and JES phi
        TH1F* hmetmuon[4][4];                // MET corrected for muons and JES
        TH1F* hmetmuonPhi[4][4];             // MET corrected for muons and JES phi
        TH1F* hmetmuonjes[4][4];             // MET corrected for muons and JES
        TH1F* hmetmuonjesPhi[4][4];          // MET corrected for muons and JES phi
        TH1F* htcmet[4][4];                  // tc MET
        TH1F* htcmet_sqrtht[4][4];           // tc MET
        TH1F* htcmetPhi[4][4];               // tc MET phi
        TH1F* hpfmet[4][4];                  // tc MET
        TH1F* hpfmetPhi[4][4];               // tc MET phi

        TH2F* hmetmuonVsDilepPt[4][4];       // MET vs dilepton Pt
        TH2F* hmetmuonOverPtVsDphi[4][4];    // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
        TH2F* hmetmuonjesVsDilepPt[4][4];    // MET vs dilepton Pt
        TH2F* hmetmuonjesOverPtVsDphi[4][4]; // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
        TH2F* htcmetVsDilepPt[4][4];         // tc MET vs dilepton Pt
        TH2F* htcmetOverPtVsDphi[4][4];      // tc MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

        TH2F* hdphillvsmll[4][4];            // delta phi between leptons vs dilepton mass
        TH1F* hptJet1[4][4];                 // Pt of 1st jet
        TH1F* hptJet2[4][4];                 // Pt of 2nd jet
        TH1F* hptJet3[4][4];                 // Pt of 3rd jet
        TH1F* hptJet4[4][4];                 // Pt of 4th jet
        TH1F* hetaJet1[4][4];                // eta of 1st jet
        TH1F* hetaJet2[4][4];                // eta of 2nd jet
        TH1F* hetaJet3[4][4];                // eta of 3rd jet
        TH1F* hetaJet4[4][4];                // eta of 4th jet
        TH1F* hptJpt1[4][4];                 // Pt of 1st JPT jet
        TH1F* hptJpt2[4][4];                 // Pt of 2nd JPT jet
        TH1F* hptJpt3[4][4];                 // Pt of 3rd JPT jet
        TH1F* hptJpt4[4][4];                 // Pt of 4th JPT jet
        TH1F* hptBtagJpt1[4][4];             // Pt of 1st JPT jet
        TH1F* hptBtagJpt2[4][4];             // Pt of 2nd JPT jet
        TH1F* hptBtagJpt3[4][4];             // Pt of 3rd JPT jet
        TH1F* hptBtagJpt4[4][4];             // Pt of 4th JPT jet
        TH1F* hetaJpt1[4][4];                // eta of 1st JPT jet
        TH1F* hetaJpt2[4][4];                // eta of 2nd JPT jet
        TH1F* hetaJpt3[4][4];                // eta of 3rd JPT jet
        TH1F* hetaJpt4[4][4];                // eta of 4th JPT jet=
        TH1F* hptHypJet1[4][4];              // Pt of 1st JPT jet
        TH1F* hptHypJet2[4][4];              // Pt of 2nd JPT jet
        TH1F* hptHypJet3[4][4];              // Pt of 3rd JPT jet
        TH1F* hptHypJet4[4][4];              // Pt of 4th JPT jet
        TH1F* hetaHypJet1[4][4];             // eta of 1st JPT jet
        TH1F* hetaHypJet2[4][4];             // eta of 2nd JPT jet
        TH1F* hetaHypJet3[4][4];             // eta of 3rd JPT jet
        TH1F* hetaHypJet4[4][4];             // eta of 4th JPT jet
};

#endif
