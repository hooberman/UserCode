#ifndef looper_h
#define looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TProfile.h"
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

        struct metStruct{
          
          metStruct() : 
            dmet_rms(-999.),
            dmet_mean(-999.),
            met_rms(-999.), 
            ebmet_rms(-999.), 
            eemet_rms(-999.), 
            hbmet_rms(-999.), 
            hemet_rms(-999.), 
            hfhmet_rms(-999.), 
            hfemet_rms(-999.)  {}
        
          float dmet_rms;
          float dmet_mean;
          float met_rms;
          float ebmet_rms;
          float eemet_rms;
          float hbmet_rms;
          float hemet_rms;
          float hfhmet_rms;
          float hfemet_rms;
      
        };

        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        metStruct ScanChain (TChain*, char*, const char* suffix = "", bool isData = false, int nEvents = -1,
                             float eb_threshold = 0., float ee_threshold = 0., float hb_threshold = 0.,
                             float he_threshold = 0., float hfh_threshold = 0., float hfe_threshold = 0.,
                             float hfshort_threshold = 0., float hflong_threshold = 0.);
        
        void bookHistos( char* prefix );
        void fillUnderOverFlow(TH1F *h1, float value, float weight = 1);


        
    private:
          
        //ntuple, file
        TFile *babyFile_;
        TTree *babyTree_;
    
        //histos
        TH1F* hgenmet;
        TH1F* hcalo_met;
        TH1F* hpfr_met;
        TH1F* hpfc_met;
        TH1F* hpfr_met_nothresh;

        TH1F* hcalo_dmet;
        TH1F* hpfc_dmet;

        TH1F* hgensumet;
        TH1F* hcalo_sumet;
        TH1F* hpfr_sumet;
        TH1F* hpfc_sumet;
        TH1F* hpfr_sumet_nothresh;
             
        TH1F* hcalo_ebsumet;   
        TH1F* hcalo_eesumet;   
        TH1F* hcalo_hbsumet;   
        TH1F* hcalo_hesumet;   
        TH1F* hcalo_hfhsumet;  
        TH1F* hcalo_hfesumet;  
                
        TH1F* hpfr_ebsumet;    
        TH1F* hpfr_eesumet;    
        TH1F* hpfr_hbsumet;    
        TH1F* hpfr_hesumet;    
        TH1F* hpfr_hfhsumet;   
        TH1F* hpfr_hfesumet;   

        TH1F* hcalo_ebmet;   
        TH1F* hcalo_eemet;   
        TH1F* hcalo_hbmet;   
        TH1F* hcalo_hemet;   
        TH1F* hcalo_hfhmet;  
        TH1F* hcalo_hfemet;  
                
        TH1F* hpfr_ebmet;    
        TH1F* hpfr_eemet;    
        TH1F* hpfr_hbmet;    
        TH1F* hpfr_hemet;    
        TH1F* hpfr_hfhmet;   
        TH1F* hpfr_hfemet;   

        TH1F* hebe;        
        TH1F* heee;        
        TH1F* hhbe;        
        TH1F* hhee;        
        TH1F* hhfhe;       
        TH1F* hhfee;       
  
        TH1F*   hpfr_metxy;       
        TH1F*   hpfr_ebmetxy;     
        TH1F*   hpfr_eemetxy;     
        TH1F*   hpfr_hbmetxy;     
        TH1F*   hpfr_hemetxy;     
        TH1F*   hpfr_hfhmetxy;    
        TH1F*   hpfr_hfemetxy;    
        TH1F*   hdpfr_met;        

        TH1F*   hpfc_ebsumet;     
        TH1F*   hpfc_eesumet;     
        TH1F*   hpfc_hbsumet;     
        TH1F*   hpfc_hesumet;     
        TH1F*   hpfc_hfhsumet;    
        TH1F*   hpfc_hfesumet;   

        TH1F*   hpfc_ebmet;     
        TH1F*   hpfc_eemet;     
        TH1F*   hpfc_hbmet;     
        TH1F*   hpfc_hemet;     
        TH1F*   hpfc_hfhmet;    
        TH1F*   hpfc_hfemet; 


        TProfile* tmet;

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
        Float_t dphixmet_;
        Float_t metPar_;
        Float_t metPerp_;
	Float_t tcmetphi_;
	Float_t tcsumet_;




};





#endif
