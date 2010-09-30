#ifndef looper_h
#define looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
#include <fstream>


class TChain;

class looper
{
    public:
        looper() {};
        ~looper() {
            delete babyFile_;
            delete babyTree_;
        };

        enum metAlgo   { e_makeTemplate = 0, e_photonSelection = 1, e_ZSelection = 2};

        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, const char*, bool isData, bool calculateTCMET = false, metAlgo algo = e_makeTemplate, int nEvents = -1, float kFactor = 1.);
        void bookHistos();
	bool isGoodTrack(int, bool usePV = false);
        float deltaPhi( float phi1 , float phi2);
        bool passZSelection();
        bool passMuon_Nominal();
        bool passMuon_NominalTTbar();
        bool passMuon_NominalTTbarV2();
        bool passElectron_ttbarV1( bool isData );
        bool passElectron_ttbar( bool isData );
        bool passElectron_cand01();
        void fillUnderOverFlow(TH1F *h1, float value, float weight);

    private:
        
        metAlgo algo_;
        
        //ntuple, file
        TFile *babyFile_;
        TTree *babyTree_;
    
        //histos

        //triggers
        Int_t HLT_L1Jet6U_;
        Int_t HLT_L1Jet10U_;
        Int_t HLT_Jet15U_;
        Int_t HLT_Jet30U_;
        Int_t L1_SingleEG5_;
        Int_t HLT_Photon10_L1R_;
        Int_t HLT_Photon15_L1R_;
        Int_t HLT_Photon10_Cleaned_L1R_;
        Int_t HLT_Photon15_Cleaned_L1R_;
        Int_t HLT_Photon20_Cleaned_L1R_;

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
        Int_t   nGoodVertex_;
	Float_t weight_;
        Float_t pthat_;

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
        Float_t dphixmet_;
        Float_t metPar_;
        Float_t metPerp_;

	Float_t tcmet_;
	Float_t tcmetphi_;
	Float_t tcsumet_;

	Float_t tcmetNew_;
	Float_t tcmetphiNew_;
	Float_t tcsumetNew_;

        // photon stuff
        Int_t    nPhotons_;
        Float_t  etg_;
        Float_t  etag_;
        Float_t  phig_;
        Float_t  hoe_;
        Float_t  eciso_;
        Float_t  hciso_;
        Float_t  tkiso_;
        Float_t  swiss_;
        Float_t  seed_;
        Float_t  r4_;
        Float_t  s4_;
        Float_t  drel_;

        Int_t    photon_pixelseed_;        
	Float_t  photon_e15_;              
	Float_t  photon_e25max_;           
	Float_t  photon_e33_;              
	Float_t  photon_e55_;              
	Float_t  photon_ecalIso03_;        
	Float_t  photon_ecalIso04_;        
	Float_t  photon_hcalIso03_;        
	Float_t  photon_hcalIso04_;        
	Float_t  photon_ntkIsoHollow03_;   
        Float_t  photon_ntkIsoHollow04_;   
	Float_t  photon_ntkIsoSolid03_;    
	Float_t  photon_ntkIsoSolid04_;    
	Float_t  photon_sigmaEtaEta_;      
	Float_t  photon_sigmaIEtaIEta_;    
	Float_t  photon_tkisoHollow03_;    
	Float_t  photon_tkisoHollow04_;    
	Float_t  photon_tkisoSolid03_;     
	Float_t  photon_tkisoSolid04_;     
       
        // jet stuff
        Int_t   nJets_;
        Int_t   nJets40_;
        Float_t sumJetPt_;
        Float_t sumJetPt10_;
        Float_t vecJetPt_;

        //photon-matched jet stuff
        Float_t jet_pt_;            
        Float_t jet_dr_;            
        Float_t jet_eta_;            
        Float_t jet_energy_;            
        Float_t jet_chg_emfrac_;    
        Float_t jet_chg_hadfrac_;   
        Float_t jet_neu_emfrac_;    
        Float_t jet_neu_hadfrac_;   
        Int_t   jet_nchg_;              
        Int_t   jet_nmuon_;         
        Int_t   jet_nneu_;          
        Float_t jet_dphimet_;       
        Float_t jet_dpt_;           
        Float_t jet_drgen_;         

        //leading jet stuff
        Float_t jetmax_pt_;
        Float_t jetmax_dphimet_;

        //Z stuff
        Int_t   passz_;
        Int_t   passe_ttbar_;
        Int_t   passe_ttbarV1_;
        Int_t   passe_cand01_;
        Int_t   passm_nomttbar_;
        Int_t   passm_nomttbarV2_;
        Int_t   passm_nom_;
        Int_t   pdgid_;
        Float_t ptll_;
        Float_t ptlt_;
        Float_t etall_;
        Float_t etalt_;
        Float_t dilmass_;
        Float_t dilpt_;
        Int_t   flagll_;
        Int_t   flaglt_;

        TH1F*   hgenps_pthat;
        TH1F*   hphotonpt;

        TH1F* metPredicted;
        TH1F* metObserved;
        TH1F* metParPredicted;
        TH1F* metParObserved;
        TH1F* metPerpPredicted;
        TH1F* metPerpObserved;

        TH1F* metPredicted_njets[11];
        TH1F* metObserved_njets[11];

        TH1F* metTemplate[11][23];
        TH1F* metParTemplate[11][23];
        TH1F* metPerpTemplate[11][23];
};





#endif
