#ifndef babylooper_h
#define babylooper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TProfile.h"
#include <vector>
#include <fstream>
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class TChain;

class babylooper
{
    public:
        babylooper() {};
        ~babylooper() {
        };

        enum selectionType   { e_QCDSelection = 0, e_photonSelection = 1, e_ZSelection = 2};

        void setBranches (TTree *tree);
        void ScanChain (TChain*, const char*, const char*, const char*, bool isData, 
                        selectionType mySelectionType = e_photonSelection, 
                        bool makeTemplate = false, int nEvents = -1);
        void bookHistos();
	bool isGoodTrack(int, bool usePV = false);
        float deltaPhi( float phi1 , float phi2);
        void fillUnderOverFlow(TH1F *h1, float value, float weight);
        void fillHistos(TH1F *h1[4],    float value, float weight, int myType);
        void fillHistos(TH1F *h1[4][4], float value, float weight, int myType, int nJetsIdx);
        TH1F* getMetTemplate( TFile* file, int iTrigBin , int iJetBin , int iSumJetPtBin , 
                              int iBosonPtBin , int iVtxBin, float Zpt , float weight );
        void setErrors( TFile* file,  TH1F* hist , int n[3][7] );
        void setErrors( TFile* file,  TH1F* hist , int n[4][3][7] );
	TH1F* correctedMetTemplate( TH1F* h_metTemplate , float ptZ );

    private:
        
        //metAlgo algo_;
        bool makeTemplate_;
        selectionType selection_;

        //histos

        //triggers
        Int_t HLT_L1Jet6U_;
        Int_t HLT_L1Jet10U_;
        Int_t HLT_Jet15U_;
        Int_t HLT_Jet30U_;
        Int_t L1_SingleEG5_;
        Int_t HLT_Photon10_L1R_;
        Int_t HLT_Photon15_L1R_;
        Int_t HLT_Photon20_L1R_;
        Int_t HLT_Photon10_Cleaned_L1R_;
        Int_t HLT_Photon15_Cleaned_L1R_;
        Int_t HLT_Photon20_Cleaned_L1R_;

        // event stuff
        Int_t   run_;
        Int_t   lumi_;
        Int_t   event_;
        Int_t   nvtx_;
        Int_t   npfmuons_;
        Float_t weight_;
        Float_t vtxweight_;
        Float_t davtxweight_;
        Float_t pthat_;
	Float_t mjj_;

	Int_t   nbl_;
	Int_t   nbm_;
	Int_t   nbt_;

	// genmet stuff
	Float_t genmet_;
	Float_t genmetphi_;
	Float_t gensumet_;

	// pfmet stuff
	Float_t pfmet_;
	Float_t pfmett1new_;
	Float_t pfmetcor_;
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
	Float_t tcmetNew_;
	Float_t tcmetphiNew_;
	Float_t tcsumetNew_;
        Float_t dphixmet_;
        Float_t metPar_;
        Float_t metPerp_;


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
	Int_t    nlep_;

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
        Float_t vecJetPt_;
        Int_t   nbtags_;
        Float_t dphijetmet_;

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
        Int_t   pfjetid_;         
        Int_t   jet_nneu_;          
        Float_t jet_dphimet_;       
        Float_t jet_dpt_;           
        Float_t jet_drgen_;    

        Float_t drjet_ll_;
        Float_t jetpt_ll_;
        Int_t   pfjetid_ll_;
        Float_t drjet_lt_;
        Float_t jetpt_lt_;
        Int_t   pfjetid_lt_;
        Int_t   failjetid_;
        Float_t maxemf_;


        //Z stuff
        Int_t   passe_ll_ttbar_;
        Int_t   passe_ll_ttbarV1_;
        Int_t   passe_ll_ttbarV2_;
        Int_t   passe_ll_cand01_;
        Int_t   passm_ll_nomttbar_;
        Int_t   passm_ll_nomttbarV2_;
        Int_t   passm_ll_nom_;
        Int_t   passe_lt_ttbar_;
        Int_t   passe_lt_ttbarV1_;
        Int_t   passe_lt_ttbarV2_;
        Int_t   passe_lt_cand01_;
        Int_t   passm_lt_nomttbar_;
        Int_t   passm_lt_nomttbarV2_;
        Int_t   passm_lt_nom_;
        Int_t   TMLastStationTight_ll_;
        Int_t   TMLastStationTight_lt_;

        //Int_t   pdgid_;
        Int_t   leptype_;
        Float_t ptll_;
        Float_t ptlt_;
        Float_t ptllpf_;
        Float_t ptltpf_;
        Float_t etall_;
        Float_t etalt_;
        Float_t phill_;
        Float_t philt_;
        Int_t   flagll_;
        Int_t   flaglt_;
        Float_t dilmass_;
        Float_t dilmasspf_;
        Float_t dilmasscor_;
        Float_t dilpt_;

	Int_t   ee_;
	Int_t   mm_;
	Int_t   em_;
	Int_t   me_;

	LorentzVector*  lep1_;
	LorentzVector*  lep2_;
	LorentzVector*  lep3_;
	LorentzVector*  pflep1_;
	LorentzVector*  pflep2_;

	Int_t el1tv_;
	Int_t el2tv_;
	Int_t ngennu_;

        //leading jet stuff
        Float_t jetmax_pt_;
        Float_t jetmax_dphimet_;
      
	TH1F* hpfdilmassee;
	TH1F* hpfdilmassmm;
        TH1F* hyield_0j;
        TH1F* hyield_1j;
        TH1F* hyield_2j;
        TH1F* hyield_g2j;
        TH1F* hyield;
        TH1F* hyield_pfmet30;
        TH1F* hyield_pfmet60;
        TH1F* hyield_pfmet120;
        TH1F* hpthad[5];
        TH1F* hgenps_pthat;
        TH1F* hphotonpt;
        TH1F* hr4;
	TH1F* hnVtx;
	TH1F* hvecJetPt;

	TH1F* hgenmet_all;
	TH1F* hgenmet_pass;
	TProfile* hresponse;
        TH1F* metPredicted;
        int   n_metPredicted[3][7];
        int   n_metPredicted_ee[3][7];
        int   n_metPredicted_mm[3][7];
        int   nphoton_metPredicted[4][3][7];
        int   nphoton_metPredicted_ee[4][3][7];
        int   nphoton_metPredicted_mm[4][3][7];
        int   nqcd_metPredicted[4][3][7];
        int   nqcd_metPredicted_ee[4][3][7];
        int   nqcd_metPredicted_mm[4][3][7];
        TH1F* metObserved;
        TH1F* metPredicted_sf;
        TH1F* metObserved_sf;
        TH1F* metPredicted_df;
        TH1F* metObserved_df;
        TH1F* metObserved_df_nozveto;
        TH1F* metPredicted_ptlt40;
        TH1F* metObserved_ptlt40;
        TH1F* metPredicted_ptlt50;
        TH1F* metObserved_ptlt50;
        TH1F* metPredicted_ptgt50;
        TH1F* metObserved_ptgt50;
        TH1F* metPredicted_pt40_60;
        TH1F* metObserved_pt40_60;
        TH1F* metPredicted_ptgt60;
        TH1F* metObserved_ptgt60;
        TH1F* metPredicted_ee;
        TH1F* metObserved_ee;
        TH1F* metPredicted_mm;
        TH1F* metObserved_mm;
        TH1F* metParPredicted;
        TH1F* metParObserved;
        TH1F* metPerpPredicted;
        TH1F* metPerpObserved;

        TH1F* hdilmass[4][4];
        TH1F* hdilmass_pfmet60[4][4];
        TH1F* htcmet[4][4];
        TH1F* htcmetNew[4][4];
        TH1F* hpfmet[4][4];

	TH1F* hpfmet_ebeb;
	TH1F* hpfmet_ebee;
	TH1F* hpfmet_eeeep;
	TH1F* hpfmet_eeeem;
	TH1F* hpfmet_mm;

        TH1F* metPredicted_njets[11];
        TH1F* metObserved_njets[11];
        
        TH1F* tcmetTemplate[3][7][4];
        TH1F* tcmetNewTemplate[3][7][4];
        TH1F* pfmetTemplate[3][7][4];

        TH1F* tcmetTemplate_njets_ht_nvtx[3][7][3];
        TH1F* tcmetNewTemplate_njets_ht_nvtx[3][7][3];
        TH1F* pfmetTemplate_njets_ht_nvtx[3][7][3];

        TH1F* tcmetTemplate_combined[3][7];
        TH1F* tcmetNewTemplate_combined[3][7];
        TH1F* pfmetTemplate_combined[3][7];

        int nTemplate[3][7];

        ofstream ofile_tcmet;
        ofstream ofile_events;
};





#endif
