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
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, const char*, bool isData, bool calculateTCMET = false, int nEvents = -1);
        void bookHistos();
	bool isGoodTrack(int, bool usePV = false);
        float deltaPhi( float phi1 , float phi2);

    private:
        
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

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
	
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

        // photon stuff
        Int_t   nPhotons_;
        Float_t maxPhotonPt_;

        // jet stuff
        Int_t   nJets_;
        Float_t sumJetPt_;

        //leading jet stuff
        Float_t jet_pt_;            
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
        Int_t   jet_overlapPhoton_;

};





#endif
