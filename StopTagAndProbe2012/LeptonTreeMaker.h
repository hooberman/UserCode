#ifndef WW_LeptonTreeMaker_h
#define WW_LeptonTreeMaker_h

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "wwtypes.h"
#include "TChain.h"
#include <fstream>
#include <vector>

#include "../../../Smurf/Core/LeptonTree.h"
#include "../../../Smurf/Core/SmurfTree.h"

#include "../Tools/ElectronIDMVA.h"
#include "../CORE/jetcorr/FactorizedJetCorrector.h"
#include "analysisEnums.h"
#include "SmurfDataTypes.h"


class TChain;
//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorD;
//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
//typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

class LeptonTreeMaker {

    public:

        LeptonTreeMaker() {};
        LeptonTreeMaker(bool lockToCoreSelectors, bool useLHeleId, 
                bool useMVAeleId, bool doDYNNLOw, unsigned int prescale, bool realData);
        ~LeptonTreeMaker();

        void SetBatchMode(bool on) {
            batchMode_ = on;
        }

        void ScanChain(TString outfileid, 
                TChain* chain,
				SmurfTree::DataType sample,
                double integratedLumi,
                double xsec,
                int nProcessedEvents,
                bool identifyEvents = false,
                bool realData = false,
                TString cms2_json_file = "");


    private:

        bool loosefo(unsigned int index);

        //
        // efficiencies
        //

        void MakeElectronTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);
        void MakeMuonTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);

        // this one is a little bit special
        void MakeElectronRecoTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);

        //
        // fake rates
        //

        void MakeElectronFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection);
        void MakeMuonFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection);

        //
        // utilities
        //

        float GetAwayJetPt(LorentzVector lep1, LorentzVector lep2);

        //
        // common variables
        //

        void SetCommonTreeVariables(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);

	//VofP4 getStopJets();
	//VofP4 getGoodLeptons();

        //
        // data members 
        //

        bool lockToCoreSelectors_;
        bool useLHeleId_;
        bool useMVAeleId_;
        bool doDYNNLOw_;
        unsigned int prescale_;
        bool realData_;
        bool batchMode_;

        // jet correction
        FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3_;



        // ID MVAas
        ElectronIDMVA *electronIdMVA_;


	Float_t drprobe_;
	Float_t tkiso_old_;
	Float_t tkiso_new_;
	Float_t tkiso_new_noveto_;
	Float_t tkiso_new_pt1_;
	Float_t tkiso_new_pt2_;
	Float_t tkiso_new_pt3_;
	Float_t tkiso_new_pt4_;
	Float_t tkiso_new_pt5_;
	Float_t tkisov4_;
	Float_t vtxweight_;
	Float_t isoch_;

	Int_t   muispf_;
	Int_t   mustahits_;
	Int_t   munsegs_;
	Int_t   munpix_;
	Int_t   munlayers_;
	Int_t   mutrk_;
	Int_t   muglb_;
	Float_t mud0_;
	Float_t mudz_;
	Float_t muchi2ndf_;

	Int_t   probepassid_;
	Int_t   probepassiso_;

	Float_t probepfpt_;
	Float_t eoverpin_;

	unsigned int nbl_;
	unsigned int nbm_;
	
	/* unsigned int leptype_; */
	/* unsigned int elpassid_;    */
	/* unsigned int elpassiso_;   */
	/* unsigned int mupassid_;    */
	/* unsigned int mupassiso_; */

	unsigned int HLT_Ele27_WP80_tag_;
	unsigned int HLT_Ele27_WP80_probe_;
	unsigned int HLT_IsoMu24_tag_;
	unsigned int HLT_IsoMu24_probe_;

	Float_t tagiso_;
	Float_t probeiso_;
};

#endif
