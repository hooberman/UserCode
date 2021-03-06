#ifndef looper_h
#define looper_h

#include <vector>
#include <map>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
#include "../CORE/jetSelections.h"
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
	int getJetIndex(LorentzVector jet, FactorizedJetCorrector* jet_corrector);
	
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

	Int_t   mmht150_;
	Int_t   emht150_;
	Int_t   eeht150_;

	Float_t dilmass_;
	Float_t elnoiso1_mt_;
	Float_t munoiso1_mt_;

	Float_t munoiso1_d0pv_;
	Float_t munoiso1_d0bs_;
	Float_t elnoiso1_d0pv_;
	Float_t elnoiso1_d0bs_;
	Float_t munoiso2_d0pv_;
	Float_t munoiso2_d0bs_;
	Float_t elnoiso2_d0pv_;
	Float_t elnoiso2_d0bs_;

	Float_t munoiso1_ev_;
	Float_t munoiso1_hv_;
	Float_t munoiso2_ev_;
	Float_t munoiso2_hv_;

	Int_t elnoiso1_conv_;
	Int_t elnoiso1_conv2_;
	Int_t elnoiso1_hitp_;
	Int_t elnoiso1_inner_;
	Int_t elnoiso1_mit_;
	Int_t elnoiso2_conv_;
	Int_t elnoiso2_conv2_;
	Int_t elnoiso2_hitp_;
	Int_t elnoiso2_inner_;
	Int_t elnoiso2_mit_;

	Int_t mu5_;
	Int_t mu8_;
	Int_t mu12_;
	Int_t mu15_;
	Int_t mu17_;
	Int_t elnoiso1_wp80_;
	Int_t elnoiso1_top_;
	Int_t elnoiso2_wp80_;
	Int_t elnoiso2_top_;
	Int_t elnoiso3_wp80_;
	Int_t elnoiso3_top_;
	Int_t elnoiso4_wp80_;
	Int_t elnoiso4_top_;
	Int_t munoiso1_mu24_;
	Int_t munoiso2_mu24_;
	Int_t munoiso3_mu24_;
	Int_t munoiso4_mu24_;
	Int_t munoiso1_mu30_;
	Int_t munoiso2_mu30_;
	Int_t munoiso3_mu30_;
	Int_t munoiso4_mu30_;

	Float_t elnoiso1_iso_;
	Float_t elnoiso1_isofj_;
	Float_t elnoiso1_isovtx_;
	Float_t elnoiso1_isopf_;
	Float_t elnoiso1_isopffj03_;
	Float_t elnoiso1_isopffj04_;

	Float_t elnoiso2_iso_;
	Float_t elnoiso2_isofj_;
	Float_t elnoiso2_isovtx_;
	Float_t elnoiso2_isopf_;
	Float_t elnoiso2_isopffj03_;
	Float_t elnoiso2_isopffj04_;

	Float_t elnoiso3_iso_;
	Float_t elnoiso3_isofj_;
	Float_t elnoiso3_isovtx_;
	Float_t elnoiso3_isopf_;
	Float_t elnoiso3_isopffj03_;
	Float_t elnoiso3_isopffj04_;

	Float_t elnoiso4_iso_;
	Float_t elnoiso4_isofj_;
	Float_t elnoiso4_isovtx_;
	Float_t elnoiso4_isopf_;
	Float_t elnoiso4_isopffj03_;
	Float_t elnoiso4_isopffj04_;

	Float_t munoiso1_iso_;
	Float_t munoiso1_isofj_;
	Float_t munoiso1_isovtx_;
	Float_t munoiso1_isopf_;
	Float_t munoiso1_isodb03_;
	Float_t munoiso1_isodb04_;
	Float_t munoiso2_iso_;
	Float_t munoiso2_isofj_;
	Float_t munoiso2_isovtx_;
	Float_t munoiso2_isopf_;
	Float_t munoiso2_isodb03_;
	Float_t munoiso2_isodb04_;
	Float_t munoiso3_iso_;
	Float_t munoiso3_isofj_;
	Float_t munoiso3_isovtx_;
	Float_t munoiso3_isopf_;
	Float_t munoiso4_iso_;
	Float_t munoiso4_isofj_;
	Float_t munoiso4_isovtx_;
	Float_t munoiso4_isopf_;

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
        LorentzVector*  pjetnoiso1_; 
        LorentzVector*  pjetnoiso2_; 
        LorentzVector*  pjetnoiso3_; 
        LorentzVector*  pjetnoiso4_; 
        LorentzVector*  lep1_; 
        LorentzVector*  lep2_; 
        LorentzVector*  lep3_; 
        LorentzVector*  lep4_; 
        LorentzVector*  elnoiso1_; 
        LorentzVector*  elnoiso2_; 
        LorentzVector*  elnoiso3_; 
        LorentzVector*  elnoiso4_; 
        LorentzVector*  munoiso1_; 
        LorentzVector*  munoiso2_; 
        LorentzVector*  munoiso3_; 
        LorentzVector*  munoiso4_; 
        Int_t   lep1_top_; 
	Float_t pjet1_res_;
	Float_t pjet2_res_;
	Float_t pjet3_res_;
	Float_t pjet4_res_;
	Float_t pjet1_L1Fast_;
	Float_t pjet2_L1Fast_;
	Float_t pjet3_L1Fast_;
	Float_t pjet4_L1Fast_;
	Float_t pjet1_L2L3_;
	Float_t pjet2_L2L3_;
	Float_t pjet3_L2L3_;
	Float_t pjet4_L2L3_;
        char    dataset_[500];
        UInt_t  run_;
        UInt_t  lumi_;
        UInt_t  event_;
	Int_t   njets_;
	Int_t   njetsnoiso_;
	Int_t   ncjets_;
	Float_t ht_;
	Float_t htnoiso_;
	Float_t htc_;
	Float_t pfmet_;
	Float_t pfmetphi_;
	Float_t pfsumet_;
	Int_t   ele8dijet30_;
	Int_t   ngoodlep_;
	Int_t   ngoodmu_;
	Int_t   ngoodel_;
	Int_t   nelnoiso_;
	Int_t   nmunoiso_;
	Int_t   nosel_;
	Int_t   ndavtx_;
	Int_t   nvtx_;
	Int_t   eledijet_hg_;
	Int_t   eledijetmht15_;
	Int_t   eledijetmht25_;
	Int_t   eledijet_;
	Int_t   ele27dijet25_;
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
	Float_t mindrej_;
	Float_t mindrmj_;
	Float_t eledijet_trigmindr_ejet_;
	Float_t eletrijet_trigmindr_ejet_;
	Float_t mudijet_trigmindr_mujet_;
	Float_t mutrijet_trigmindr_mujet_;
	Int_t   mudijetmht15_;
	Int_t   mudijetmht25_;
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

	// custom top triggers
	Int_t eltrijetcaloisovl_;
	Int_t eltrijettest_;

	// top electron+jets triggers
	Int_t eltrijet_;             
	Int_t eltrijetbackup_;       
	Int_t eldijet_;              
	Int_t eljet_;                
	Int_t elnoisotrijet_;        
	Int_t nelnoisotrijet_;        
	Int_t elnoisotrijetbackup_;  

	// top muon+jets triggers
	Int_t mutrijet_;             
	Int_t mutrijetbackup_;       
	Int_t mudijet_;              
	Int_t mujet_;                
	Int_t munoisotrijet_;        
	Int_t munoisotrijetbackup_;  

	// non-isolated dilepton-HT triggers
	Int_t eeht175_;              
	Int_t eeht225_;              
	Int_t mmht175_;              
	Int_t mmht225_;              
	Int_t emht175_;              
	Int_t emht225_;              

	// isolated dilepton-HT triggers
	Int_t mmisoht175_;           
	Int_t mmisoht225_;           
	Int_t emisoht175_;           
	Int_t emisoht225_;           

	// isolated single muon triggers
	Int_t isomu20_;              
	Int_t isomu24_;              
	Int_t isomu30_;              
	Int_t isomu34_;              
	Int_t isomu40_;              

	// custom isomu triggers
	Int_t isomu24test_;              
	Int_t iso20mu24_;              
	Int_t iso30mu24_;              
	Int_t iso40mu24_;              

	// non-isolated single muon triggers
	Int_t mu24_;                 
	Int_t nmu24_;                 
	Int_t mu30_;                 
	Int_t mu40_;                 
	Int_t mu50_;                 
	
	// single-electron triggers
	Int_t el27wp80_;             
	Int_t el27wp70_;             
	Int_t el27_;                 
	Int_t el32_;                 

	// multi-jet triggers
	Int_t quadjet70_;            
	Int_t quadjet80_;            
	Int_t quadjet90_;            

	// single photon triggers
	Int_t photon20_;             
	Int_t photon30_;             
	Int_t photon50_;             
	Int_t photon75_;             
	Int_t photon90_;             
	Int_t photon135_;            
	Int_t photon150_;            
	Int_t photon160_;            

	// higgs single photon triggers
	Int_t hphoton22_;            
	Int_t hphoton36_;            
	Int_t hphoton50_;            
	Int_t hphoton75_;            
	Int_t hphoton90_;            

	// single electron utility triggers
	Int_t el8_;                  
	Int_t el8jet30_;             
	Int_t el17_;                 
	Int_t el17jet30_;            
	Int_t el8vl_;                
	Int_t el17vl_;               
	Int_t el8noiso_;                  
	Int_t el8noisojet30_;                  
	Int_t el30noiso_;                  

	// Higgs dilepton triggers
	Int_t ee_;                   
	Int_t mmtrk_;                
	Int_t mm_;                   
	Int_t em_;                   
	Int_t me_;                   

	Int_t munoiso1_emht175_;
	Int_t munoiso1_mmht175_;
	Int_t munoiso2_emht175_;
	Int_t munoiso2_mmht175_;

};

#endif
