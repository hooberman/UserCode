
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "../CORE/CMS2.h"

using namespace tas;

class TPEff;

class MyScanChain {

    public:

        //
        // if even does not equal -1
        // then use it to determine if to take
        // odd or even event numbers only
        // - in this case, weight is adjusted as needed
        //
        int ScanChain(bool isData, std::string prefix, TChain* chain, float kFactor = 1.0, int useOdd=-1);
        void setGoodRunList(std::string fname);

    private:

        //
        // work out what passed for the electron and muon probes
        //

        void doElectrons(const float &weight, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);
        void doMuons(const float &weight);

        //
        //
        // utility functions

        // check if object passed a given trigger
        bool objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt = 0.0);
        // find the index of a given trigger with name trigName
        int findTriggerIndex(TString trigName);

        // init baby branches with dummy values
        void InitBaby();

        bool isData_;
        bool verbose_;

        //
        // baby ntuple variables
        // for debugging purposes
        //

        TFile *babyFile_;
        TTree *babyTree_;

        Int_t run_;
        Int_t ls_;
        Int_t evt_;
        Float_t rho_; // average energy per unit area from PU 
	Int_t goodvtx_; //isGoodVertex
	Int_t gooddavtx_; //isGoodDAVertex
	Float_t vtxweight_;

        Float_t mass_;
        Float_t weight_;

        Float_t pt1_;
        Float_t eta1_;
        Float_t phi1_;
        Float_t sieie1_;
        Float_t detain1_;
        Float_t hoe1_;
        Float_t dphiin1_;
        Float_t fbrem1_;
        Float_t eopin1_;
        Int_t q1_;
        Int_t q3agree1_;
        Int_t distdcot1_;
        Int_t hitpattern1_;
        Float_t d01_;
        Float_t ecaliso1_;
        Float_t hcaliso1_;
        Float_t tkiso1_;
        Float_t reliso1_;
        Float_t relisont1_;
        Float_t reliso1_fastjet_;
        Int_t smurfv11_;
        Int_t smurfv21_;

	Float_t pt2_;
        Float_t eta2_;
        Float_t phi2_;
        Float_t sieie2_;
        Float_t detain2_;
        Float_t hoe2_;
        Float_t dphiin2_;
        Float_t fbrem2_;
        Float_t eopin2_;
        Int_t q2_;
        Int_t q3agree2_;
        Int_t distdcot2_;
        Int_t hitpattern2_;
        Float_t d02_;
        Float_t ecaliso2_;
        Float_t hcaliso2_;
        Float_t tkiso2_;
        Float_t reliso2_;
        Float_t relisont2_;
        Float_t reliso2_fastjet_;
        Int_t smurfv12_;
        Int_t smurfv22_;

        Int_t type_;

        Int_t tt_;
        Int_t passPreselection_;
        Int_t passIso_;
        Int_t passIso_fastjet_;
        Int_t passId_;
        Int_t passTrigger_;

};

