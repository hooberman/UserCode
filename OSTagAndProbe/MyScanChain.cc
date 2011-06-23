#include <exception>

#include "MyScanChain.h"

#include "../CORE/utilities.h"
#include "../CORE/muonSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/eventSelections.h"
#include "../CORE/triggerUtils.h"
#include "../CORE/ttbarSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/CMS2.h"

#include "../Tools/goodrun.cc"
#include "../Tools/tools.cc"
#include "../Tools/vtxreweight.cc"

const float mu_rel_iso_cut = 0.15;

bool MyScanChain::objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt) 
{

    float drMin = 999.99;
    for (size_t i = 0; i < trigObjs.size(); ++i)
    {
        if (trigObjs[i].Pt() < pt) continue;
        float dr = dRbetweenVectors(trigObjs[i], obj);
        if (dr < drMin) drMin = dr;
    }

    if (drMin < 0.1) return true;
    return false;

}

int MyScanChain::findTriggerIndex(TString trigName)
{
    vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
    vector<TString>::const_iterator end_it = hlt_trigNames().end();
    vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
    if(found_it != end_it) return found_it - begin_it;
    return -1;
}

void MyScanChain::InitBaby()
{
    run_		= cms2.evt_run();
    ls_			= cms2.evt_lumiBlock();
    evt_		= cms2.evt_event();

    mass_		= -999.;
    weight_		= -999.;
    rho_		= -999.; // average energy per unit area from PU 
    goodvtx_		= -1; //isGoodVertex
    gooddavtx_		= -1;

    pt1_		= -999.;
    eta1_		= -999.;
    phi1_		= -999.;
    sieie1_		= -999.;
    detain1_		= -999.;
    hoe1_		= -999.; 
    dphiin1_		= -999.; 
    fbrem1_		= -999.; 
    eopin1_		= -999.;
    q1_			= -1; 
    q3agree1_		= -1; 
    distdcot1_		= -1;
    hitpattern1_	= -1; 
    d01_		= -999.; 
    ecaliso1_		= -999.; 
    hcaliso1_		= -999.; 
    tkiso1_		= -999.;
    reliso1_		= -999.;
    relisont1_		= -999.;
    reliso1_fastjet_    = -999.;
    smurfv11_		= -1; 
    smurfv21_		= -1; 
   
    pt2_		= -999.;
    eta2_		= -999.;
    phi2_		= -999.;
    sieie2_		= -999.;
    detain2_		= -999.;
    hoe2_		= -999.; 
    dphiin2_		= -999.; 
    fbrem2_		= -999.; 
    eopin2_		= -999.;
    q2_			= -1; 
    q3agree2_		= -1; 
    distdcot2_		= -1; 
    hitpattern2_	= -1; 
    d02_		= -999.; 
    ecaliso2_		= -999.; 
    hcaliso2_		= -999.; 
    tkiso2_		= -999.; 
    reliso2_		= -999.;
    relisont2_		= -999.;
    reliso2_fastjet_	= -999.;
    smurfv12_		= -1; 
    smurfv22_		= -1; 

    type_		= -1;

    tt_			= -1;
    passIso_		= -1;
    passIso_fastjet_	= -1;
    passId_		= -1;
    passTrigger_	= -1;

}

//
// set good run list
//

void MyScanChain::setGoodRunList(std::string fname)
{
    set_goodrun_file(fname.c_str());
}

//
// if even does not equal -1
// then use it to determine if to take
// odd or even event numbers only
// - in this case, weight is adjusted as needed
//

int MyScanChain::ScanChain(bool isData, std::string prefix, TChain* chain, float kFactor, int useOdd) {

    //
    // Configuration
    // 

    isData_ = isData;
    verbose_ = false;

    TObjArray *listOfFiles = chain->GetListOfFiles();
    unsigned int nEventsChain=0;
    nEventsChain = chain->GetEntries();
    unsigned int nEventsTotal = 0;

    //initialize
    bool verbose = true;
    char* vtxfile = "/tas/benhoob/vtxreweight/vtxreweight_Spring11MC_153pb_Zselection.root";
    set_vtxreweight_rootfile( vtxfile , verbose );
    
    //
    // Baby Ntuple
    // 
 
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();

    babyFile_ = new TFile(Form("output/baby_%s.root" , prefix.c_str()), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree"		, "A Baby Ntuple");

    babyTree_->Branch("run"		,       &run_,         "run/I"         );
    babyTree_->Branch("ls"		,       &ls_,          "ls/I"          );
    babyTree_->Branch("evt"		,       &evt_,         "evt/I"         );

    babyTree_->Branch("mass"		,       &mass_,        "mass/F"         );
    babyTree_->Branch("weight"		,       &weight_,      "weight/F"         );
    babyTree_->Branch("rho"		,	&rho_,      "rho/F"         );
    babyTree_->Branch("goodvtx"		,	&goodvtx_,      "goodvtx/I"         );
    babyTree_->Branch("gooddavtx"	,	&gooddavtx_,      "gooddavtx/I"         );
    babyTree_->Branch("vtxweight"	,	&vtxweight_,      "vtxweight/F"         );

    babyTree_->Branch("pt1"		,       &pt1_,          "pt1/F"         );
    babyTree_->Branch("eta1"		,       &eta1_,         "eta1/F"         );
    babyTree_->Branch("phi1"		,       &phi1_,         "phi1/F"         );
    babyTree_->Branch("sieie1"		,       &sieie1_,       "sieie1/F"         );
    babyTree_->Branch("detain1"		,       &detain1_,      "detain1/F"         );
    babyTree_->Branch("hoe1"		,       &hoe1_,       	"hoe1_/F"         );
    babyTree_->Branch("dphiin1"		,       &dphiin1_,      "dphiin1_/F"         );
    babyTree_->Branch("fbrem1"		,       &fbrem1_,       "fbrem1_/F"         );
    babyTree_->Branch("eopin1"		,       &eopin1_,       "eopin1_/F"         );
    babyTree_->Branch("q1"		,       &q1_,       	"q1/I"         );
    babyTree_->Branch("q3agree1"	,       &q3agree1_,       	"q3agree1/I"         );
    babyTree_->Branch("distdcot1"	,	&distdcot1_,    "distdcot1/I"         );
    babyTree_->Branch("hitpattern1"	,	&hitpattern1_,  "hitpattern1/I"         );
    babyTree_->Branch("d01"		,       &d01_,       	"d01/F"         );
    babyTree_->Branch("ecaliso1"	,       &ecaliso1_,     "ecaliso1/F"         );
    babyTree_->Branch("hcaliso1"	,       &hcaliso1_,     "hcaliso1/F"         );
    babyTree_->Branch("tkiso1"		,       &tkiso1_,     	"tkiso1/F"         );
    babyTree_->Branch("reliso1"		,       &reliso1_,      "reliso1/F"     );
    babyTree_->Branch("relisont1"	,       &relisont1_,    "relisont1/F"     );
    babyTree_->Branch("reliso1F"	,       &reliso1_fastjet_,      "reliso1F/F"     );
    babyTree_->Branch("smurfv11"	,       &smurfv11_,     "smurfv11/I"         );
    babyTree_->Branch("smurfv21"	,       &smurfv21_,     "smurfv21/I"         );

    babyTree_->Branch("pt2"		,       &pt2_,          "pt2/F"         );
    babyTree_->Branch("eta2"		,       &eta2_,         "eta2/F"         );
    babyTree_->Branch("phi2"		,       &phi2_,         "phi2/F"         );
    babyTree_->Branch("sieie2"		,       &sieie2_,       "sieie2/F"         );
    babyTree_->Branch("detain2"		,       &detain2_,       "detain2/F"         );
    babyTree_->Branch("hoe2"		,       &hoe2_,       	"hoe2_/F"         );
    babyTree_->Branch("dphiin2"		,       &dphiin2_,      "dphiin2_/F"         );
    babyTree_->Branch("fbrem2"		,       &fbrem2_,       "fbrem2_/F"         );
    babyTree_->Branch("eopin2"		,       &eopin2_,       "eopin2_/F"         );
    babyTree_->Branch("q2"		,       &q2_,       	"q2/I"         );
    babyTree_->Branch("q3agree2"	,       &q3agree2_,       	"q3agree2/I"         );
    babyTree_->Branch("distdcot2"	,	&distdcot2_,    "distdcot2/I"         );
    babyTree_->Branch("hitpattern2"	,	&hitpattern2_,  "hitpattern2/I"         );
    babyTree_->Branch("d02"		,       &d02_,       	"d02/F"         );
    babyTree_->Branch("ecaliso2"	,       &ecaliso2_,     "ecaliso2/F"         );
    babyTree_->Branch("hcaliso2"	,       &hcaliso2_,     "hcaliso2/F"         );
    babyTree_->Branch("tkiso2"		,       &tkiso2_,     	"tkiso2/F"         );
    babyTree_->Branch("reliso2"		,       &reliso2_,      "reliso2/F"     );
    babyTree_->Branch("relisont2"	,       &relisont2_,    "relisont2/F"     );
    babyTree_->Branch("reliso2F"	,       &reliso2_fastjet_,      "reliso2F/F"     );
    babyTree_->Branch("smurfv12"	,       &smurfv12_,     "smurfv12/I"         );
    babyTree_->Branch("smurfv22"	,       &smurfv22_,     "smurfv22/I"         );

    babyTree_->Branch("type"		,       &type_,                 "type/I"       );

    //
    // what did the probe pass
    //
    babyTree_->Branch("tt"		,       &tt_,        	"tt/I"   );
    babyTree_->Branch("passIso"		,       &passIso_,      "passIso/I"   );
    babyTree_->Branch("passIsoF"	,       &passIso_fastjet_,      "passIsoF/I"   );
    babyTree_->Branch("passId"		,       &passId_,       "passId/I"   );
    babyTree_->Branch("passTrigger"	,	&passTrigger_,  "passTrigger/I"   );                    

    //
    // Do looping
    //

    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile f(currentFile->GetTitle());
        TTree *tree = (TTree*)f.Get("Events");
        if(!tree) {
            std::cout<<"Could not extract tree from file, skipping."<<std::endl;
            tree = 0;
            continue;
        }
        cms2.Init(tree);

        //Event Loop
        unsigned int nEvents = tree->GetEntries();
        for( unsigned int event = 0; event < nEvents; ++event) {
            cms2.GetEntry(event);
            ++nEventsTotal;

            //do some event cleaning  - false because here we only do MC
            if (!cleaning_standardAugust2010(isData))
                continue;

            // get event weight and adjust for 1fb
            float weight = 1.0;
            if (!isData) weight = cms2.evt_scale1fb();

            //
            // do good run check
            //

            if (isData) {
                if (!goodrun(cms2.evt_run(), cms2.evt_lumiBlock())) continue;
            }
	    
            //
            // do duplicate check
            //

            if (isData) {
                DorkyEventIdentifier id(cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock());
                if (is_duplicate(id)) {
                    continue;
                }
            }

            if(verbose_) std::cout << "\t passed duplicate check" << std::endl;

            //
            // split DY MC into two parts if required
            //

            if (!isData) weight *= kFactor;

            if (useOdd != -1) {
                weight *= 2.0;
                // want to use odd numbers
                if ((useOdd == 1) && (event % 2 == 0)) continue;
                // want to use even numbers
                if ((useOdd == 0) && (event % 2 != 0)) continue;
            }

            // Progress feedback to the user
            if(nEventsTotal%10000 == 0) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
                    fflush(stdout);
                }
            }

            //
            // do tag and probe stuff
            //

            // for electrons need to decide if to apply 
            // alignment correction (first argument)
            // remove dEtaIn in endcap (second argument)
            // note for 38X neither are needed
            doElectrons(weight, false, false);

            // no such corrections needed for muons

            doMuons(weight);

        }//event loop

        delete tree;

    }//file loop

    //
    // close the baby file
    //

    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();

    return 0;

} 

//----------------------------------------------------------------//
// do Tag and Probe for Electrons
//----------------------------------------------------------------//
void MyScanChain::doElectrons(const float &weight, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {

    if (verbose_) std::cout << "[MyScanChain::doElectrons]" << std::endl;

    InitBaby();
    weight_ = weight; // set weight for filling tree

    //
    // Event must pass single electron HLT
    //

    //
    // Find the index of the trigger and
    // get the corresponding vector of objects
    //
    // tag and probe trigger for 2011
    std::vector<LorentzVector> trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30;

    if (verbose_) std::cout << evt_run() << std::endl;
    if (isData_) {

        // 2011
        if (evt_run() >= 160325) {

            // tag and probe triger changed version number
	  if (evt_run() >= 160325 && evt_run() <= 160877) { 
	    trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1")];
	  }
	  else if (evt_run() >= 160888 && evt_run() <= 163261){
	    trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2")];
	  }
	  else if (evt_run() >= 163269){
	    trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3")];
	  }
     
	}
	
    }

    if (verbose_) std::cout << "[MyScanChain::doElectrons] Finished looking for triggers" << std::endl;
          
    static const cuts_t electronSelection_probe =
      (1ll<<ELEPT_010)                         | // electron p_T > 10 GeV
      (1ll<<ELEETA_250);                         // |eta| < 2.5

    static const cuts_t electronSelection_id =
      (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL) | // VBTF90, tightened to match CaloIdT+TrkIdVL
      (1ll<<ELEIP_PV_OSV2)                     | // d0(PV) < 0.04 cm, dz(PV) < 1.0 cm
      (1ll<<ELENOMUON_010)                     | // no muon dR < 0.1
      (1ll<<ELENOTCONV_HITPATTERN)             | // <=1 missing hits
      (1ll<<ELENOTCONV_DISTDCOT002)            | // dist/dcot(theta) conversion rejection
      (1ll<<ELE_NOT_TRANSITION);                 // veto electrons with SC in transition region 
    
    static const cuts_t electronSelection_iso =       
      (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
      (1ll<<ELEISO_RELNT015);                    // reliso < 0.15, non-truncated, 1 GeV EB PS
    
    static const cuts_t electronSelection_iso_fastjet =       
      (1ll<<ELEISO_ECAL_RELNT020_NPS)          | // ecal/pt < 0.2 (matches HLT requirement)
      (1ll<<ELEISO_FASTJET_REL015);              // fastjet reliso < 0.15, truncated, 1 GeV EB PS
    
    // tag is the probe + all the selection criteria
    // trigger is checked separately
    static const cuts_t electronSelection_tag = 
      electronSelection_probe |
      electronSelection_id |
      electronSelection_iso;

    //
    // get the indices of electrons that are tags
    // index refers to SC branches
    //

    vector<unsigned int> v_tags;
    for (unsigned int iEl = 0; iEl < els_p4().size(); iEl++) {
        if (els_eSC()[iEl]/cosh(els_etaSC()[iEl]) < 17.0) continue;
        if (fabs(els_p4()[iEl].Pt()) < 20.0) continue;
        if (fabs(els_p4()[iEl].Eta()) > 2.5) continue;
        if (!pass_electronSelection(iEl, electronSelection_tag, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;

        //
        // beware that trigger on the tag is run dependent
        // for electrons
        //

        if(isData_) {

            // 2011
            if (evt_run() >= 160325 && evt_run() <= 999999) {
                if (!objectPassTrigger(els_p4()[iEl], trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)) continue;
            }

        }

        v_tags.push_back(iEl);
    }

    //
    // need at least one tag
    //

    if (v_tags.size() == 0) return;

    if (verbose_) std::cout << "\t found a tag" << std::endl;

    //
    // now make the probes
    //

    vector<unsigned int> v_probes;
    for (unsigned int iEl = 0; iEl < els_p4().size(); iEl++) {
        if (els_eSC()[iEl]/cosh(els_etaSC()[iEl]) < 8.0) continue;
        if (cms2.els_p4()[iEl].Pt() < 10.) continue;
        if (fabs(cms2.els_p4()[iEl].Eta()) > 2.5) continue;
        if (!pass_electronSelection(iEl, electronSelection_probe, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;


        // FIXME - I am not sure if this is necessary
        // 2011 - it's a di-object trigger so the probe must also be able to pass
        // ------ but assume this can generally be the loose leg...
        //if (isData_) {
        //    if (evt_run() >= 160325 && evt_run() <= 999999) {
        //        if (!objectPassTrigger(els_p4()[iEl], trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)) continue;
        //    }
        //}

        v_probes.push_back(iEl);
    }

    //
    // now work out possible tag-probe pairs
    // for this event
    //

    std::vector<std::pair<unsigned int, unsigned int> > vec_pair_tp_Z;
    std::vector<std::pair<unsigned int, unsigned int> > vec_pair_tp_all;

    std::sort(v_tags.begin(), v_tags.end());
    std::sort(v_probes.begin(), v_probes.end());

    for (size_t tt = 0; tt < v_tags.size(); ++tt) {

        unsigned int idx_tag = v_tags[tt];
        LorentzVector vec_tag = cms2.els_p4()[idx_tag];

        for (size_t pp = 0; pp < v_probes.size(); ++pp) {

            unsigned int idx_probe = v_probes[pp];

            // prevent double counting
            if ((pp < tt) || idx_probe == idx_tag) continue;

            // get mass of pair
            LorentzVector  vec_probe = cms2.els_p4()[idx_probe]; 

            LorentzVector vec_tp = vec_tag + vec_probe;
            float m_tp = sqrt(fabs(vec_tp.mass2()));

            // store index, noting if the pair is in the mass window or not
            vec_pair_tp_all.push_back(std::make_pair<unsigned int, unsigned int>(idx_tag, idx_probe));

            if (m_tp > 76.0 && m_tp < 106.0) {
                vec_pair_tp_Z.push_back(std::make_pair<unsigned int, unsigned int>(idx_tag, idx_probe));
            }

        }

    }

    //
    // now pick a tag-probe pair
    // for this event
    //

    unsigned int idx_probe;
    unsigned int idx_tag;

    // if there is exactly one pair in the mass window 
    // then use that pair 
    if (vec_pair_tp_Z.size() == 1) {
        idx_probe = vec_pair_tp_Z[0].second;
        idx_tag = vec_pair_tp_Z[0].first;
    }
    // if there are no pairs in the mass window
    // and more than zero pairs with any mass
    // then pick the first available pair
    // - note - the reason for doing this is to ensure
    // a mass spectrum that extends outside the mass window
    else if (vec_pair_tp_Z.size() == 0 && vec_pair_tp_all.size() > 0) {
        idx_probe = vec_pair_tp_all[0].second;
        idx_tag = vec_pair_tp_all[0].first;				
    }
    // if no pairs were found anywhere, or more than one
    // pair was found in the mass window then reject
    // this event
    else return;

    //
    // RECORD WHAT THE PROBE PASSED
    //

    // if the probe passed all tag criteria then this is a TT event
    // and each leg can be used twice, since either leg could have triggered
    bool probeIsATag = false;
    if (find(v_tags.begin(), v_tags.end(), idx_probe) != v_tags.end()) probeIsATag = true;

    //
    // find out what probe passed
    // and write it down in the tree for later analysis
    //

    cuts_t ele_cuts_passed   = electronSelection(idx_probe, applyAlignmentCorrection, removedEtaCutInEndcap);
    LorentzVector vec_TP = cms2.els_p4()[idx_tag] + cms2.els_p4()[idx_probe];

    //
    // what does the probe pass
    //
    tt_                             = probeIsATag;
    passIso_                        = pass_electronSelectionCompareMask(ele_cuts_passed, electronSelection_iso);
    passIso_fastjet_                = pass_electronSelectionCompareMask(ele_cuts_passed, electronSelection_iso_fastjet);
    passId_                         = pass_electronSelectionCompareMask(ele_cuts_passed, electronSelection_id);
    
    if (isData_) {
        passTrigger_                = objectPassTrigger(els_p4()[idx_probe], trigObjs_HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
    }
    else {
        passTrigger_                = 1;
    }

    // number of good vertices
    int ngoodvertex = 0;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v)
    {
      if(isGoodVertex(v)) ngoodvertex++;
    }

    // number of good DA vertices
    int ngooddavertex = 0;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v)
    {
      if(isGoodDAVertex(v)) ngooddavertex++;
    }

    // 
    // event quantities
    //
    type_ = 0;
    mass_ = sqrt(fabs(vec_TP.mass2()));

    // for fastjet correction: average energy per unit area from PU   
    rho_        = cms2.evt_rho();
    goodvtx_ 	= ngoodvertex; //isGoodVertex
    gooddavtx_ 	= ngooddavertex; //isGoodDAVertex

    //
    // tag variables
    //
    pt1_  	= cms2.els_p4()[idx_tag].Pt();
    eta1_ 	= cms2.els_etaSC()[idx_tag];
    phi1_ 	= cms2.els_phiSC()[idx_tag];
    sieie1_ 	= cms2.els_sigmaIEtaIEta()[idx_tag];
    detain1_ 	= cms2.els_dEtaIn()[idx_tag];
    hoe1_	= cms2.els_hOverE()[idx_tag];
    dphiin1_	= cms2.els_dPhiIn()[idx_tag];
    fbrem1_	= cms2.els_fbrem()[idx_tag];
    eopin1_	= cms2.els_eOverPIn()[idx_tag];
    q1_		= cms2.els_charge()[idx_tag];
    q3agree1_	= !isChargeFlip3agree(idx_tag);
    distdcot1_	= isFromConversionPartnerTrack(idx_tag);
    hitpattern1_= isFromConversionHitPattern(idx_tag);
    d01_	= fabs(cms2.els_d0corr()[idx_tag]);
    ecaliso1_	= cms2.els_ecalIso()[idx_tag];
    hcaliso1_	= cms2.els_hcalIso()[idx_tag];
    tkiso1_	= cms2.els_tkIso()[idx_tag];
    reliso1_    = electronIsolation_rel(idx_tag, true);
    relisont1_  = electronIsolation_rel_v1(idx_tag, true);
    reliso1_fastjet_ = electronIsolation_rel_FastJet(idx_tag, true);
   
    smurfv11_	= -1;//pass_electronSelection(idx_tag, electronSelection_smurfV1_id, applyAlignmentCorrection, removedEtaCutInEndcap);
    smurfv21_	= pass_electronSelection(idx_tag, electronSelection_smurfV2_id, applyAlignmentCorrection, removedEtaCutInEndcap);

    //
    // probe variables
    //
    pt2_		= cms2.els_p4()[idx_probe].Pt();
    eta2_		= cms2.els_etaSC()[idx_probe];
    phi2_		= cms2.els_phiSC()[idx_probe];
    sieie2_		= cms2.els_sigmaIEtaIEta()[idx_probe];
    detain2_		= cms2.els_dEtaIn()[idx_probe];
    hoe2_		= cms2.els_hOverE()[idx_probe];
    dphiin2_		= cms2.els_dPhiIn()[idx_probe];
    fbrem2_		= cms2.els_fbrem()[idx_probe];
    eopin2_		= cms2.els_eOverPIn()[idx_probe];
    q2_			= cms2.els_charge()[idx_probe];
    q3agree2_		= !isChargeFlip3agree(idx_probe);
    distdcot2_		= isFromConversionPartnerTrack(idx_probe);
    hitpattern2_	= isFromConversionHitPattern(idx_probe);
    d02_		= fabs(cms2.els_d0corr()[idx_probe]);
    ecaliso2_		= cms2.els_ecalIso()[idx_probe];
    hcaliso2_		= cms2.els_hcalIso()[idx_probe];
    tkiso2_		= cms2.els_tkIso()[idx_probe];
    reliso2_		= electronIsolation_rel(idx_probe, true);
    relisont2_          = electronIsolation_rel_v1(idx_probe, true);
    reliso2_fastjet_	= electronIsolation_rel_FastJet(idx_probe, true);
    smurfv12_		= -1;//pass_electronSelection(idx_probe, electronSelection_smurfV1_id, applyAlignmentCorrection, removedEtaCutInEndcap);
    smurfv22_		= pass_electronSelection(idx_probe, electronSelection_smurfV2_id, applyAlignmentCorrection, removedEtaCutInEndcap);
    vtxweight_          = vtxweight();

    //
    // fill the baby tree
    //
    babyTree_->Fill();

}

void MyScanChain::doMuons(const float &weight) {

    if (verbose_) std::cout << "[MyScanChain::doMuons]" << std::endl;

    InitBaby();
    weight_ = weight; // set weight for filling tree

    //
    // Event must pass single muon HLT
    //

    //
    // Find the index of the trigger and
    // get the corresponding vector of objects
    //
    std::vector<LorentzVector> trigObjs_HLT_DoubleMu7;

    if (isData_) {
        if (evt_run() >= 160325 && evt_run() <= 163261 ) {
	  trigObjs_HLT_DoubleMu7 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_DoubleMu7_v1")]; 
        }
        else if (evt_run() >= 163269 && evt_run() <= 164236 ) {
	  trigObjs_HLT_DoubleMu7 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_DoubleMu7_v2")]; 
        }
        else if (evt_run() >= 165088 && evt_run() <= 167043 ) {
	  trigObjs_HLT_DoubleMu7 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu13_Mu8_v2")]; 
        }
        else if (evt_run() >= 167078 ) {
	  trigObjs_HLT_DoubleMu7 = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu13_Mu8_v4")]; 
        }
    }

    //
    // get the indices of muons that are tags
    //

    vector<unsigned int> v_tags;
    for(unsigned int muIdx = 0; muIdx < mus_p4().size(); muIdx++) {
        if(mus_p4()[muIdx].Pt() < 20.)                 continue;
        if(!muonId(muIdx, OSGeneric_v3))               continue;
	if( muonIsoValue(muIdx,false) > mu_rel_iso_cut )     continue;
        if (isData_) {
            if (!objectPassTrigger(mus_p4()[muIdx], trigObjs_HLT_DoubleMu7)) continue;
        }

        v_tags.push_back(muIdx);
    }

    //
    // need at least one tag
    //

    if (v_tags.size() == 0) return;

    //
    // now make the probes
    //

    vector<unsigned int> v_probes;
    for(unsigned int muIdx = 0; muIdx < mus_p4().size(); muIdx++) {
        if (mus_p4()[muIdx].Pt() < 5.0) continue;
        if (TMath::Abs(cms2.mus_p4()[muIdx].eta()) > 2.4)   continue; // eta cut
        if (((cms2.mus_type().at(muIdx)) & (1<<1)) == 0)    continue; // global muon
        if (((cms2.mus_type().at(muIdx)) & (1<<2)) == 0)    continue; // tracker muon

        // FIXME - I am not sure if this is necessary
        // 2011 - it's a di-object trigger so the probe must also be able to pass
        // ------ but assume this can generally be the loose leg...
        //if (isData_) {
        //    if (evt_run() >= 160325 && evt_run() <= 999999) {
        //        if (!objectPassTrigger(mus_p4()[muIdx], trigObjs_HLT_DoubleMu7)) continue;
        //    }
        //}

        v_probes.push_back(muIdx);

    }

    //
    // now work out possible tag-probe pairs
    // for this event
    //

    std::vector<std::pair<unsigned int, unsigned int> > vec_pair_tp_Z;
    std::vector<std::pair<unsigned int, unsigned int> > vec_pair_tp_all;


    std::sort(v_tags.begin(), v_tags.end());
    std::sort(v_probes.begin(), v_probes.end());

    for (size_t tt = 0; tt < v_tags.size(); ++tt) {

        unsigned int idx_tag = v_tags[tt];
        LorentzVector vec_tag = cms2.mus_p4()[idx_tag];

        for (size_t pp = 0; pp < v_probes.size(); ++pp) {

            unsigned int idx_probe = v_probes[pp];

            // prevent double counting
            if ((pp<tt) || idx_probe == idx_tag) continue;

            // get mass of pair
            LorentzVector vec_probe = cms2.mus_p4()[idx_probe];
            LorentzVector vec_tp = vec_tag + vec_probe;
            float m_tp = sqrt(fabs(vec_tp.mass2()));

            // store index, noting if the pair is in the mass window or not
            vec_pair_tp_all.push_back(std::make_pair<unsigned int, unsigned int>(idx_tag, idx_probe));
            if (m_tp > 76.0 && m_tp < 106.0) {
                vec_pair_tp_Z.push_back(std::make_pair<unsigned int, unsigned int>(idx_tag, idx_probe));

            }

        }

    }

    //
    // now pick a tag-probe pair
    // for this event
    //

    unsigned int idx_probe;
    unsigned int idx_tag;

    // if there is exactly one pair in the mass window 
    // then use that pair 
    if (vec_pair_tp_Z.size() == 1) {
        idx_probe = vec_pair_tp_Z[0].second;
        idx_tag = vec_pair_tp_Z[0].first;
    }
    // if there are no pairs in the mass window
    // and more than zero pairs with any mass
    // then pick the first available pair
    // - note - the reason for doing this is to ensure
    // a mass spectrum that extends outside the mass window
    else if (vec_pair_tp_Z.size() == 0 && vec_pair_tp_all.size() > 0) {
        idx_probe = vec_pair_tp_all[0].second;
        idx_tag = vec_pair_tp_all[0].first;
    }
    // if no pairs were found anywhere, or more than one
    // pair was found in the mass window then reject
    // this event
    else return;

    // if the probe passed all tag criteria then this is a TT event
    // and each leg can be used twice, since either leg could have triggered
    bool probeIsATag = false;
    if(find(v_tags.begin(), v_tags.end(), idx_probe) != v_tags.end()) probeIsATag = true;

    //
    // find out what the probe passed
    // and record in the baby ntuple for later analysis
    //

    tt_              = probeIsATag;
    passIso_         = muonIsoValue(idx_probe,false) < mu_rel_iso_cut ? 1 : 0;
    passId_          = muonIdNotIsolated( idx_probe, OSGeneric_v3 );
    passIso_fastjet_ = muonIsoValue_FastJet( idx_probe ) < mu_rel_iso_cut ? 1 : 0;

    if (isData_) {
        passTrigger_    = objectPassTrigger(mus_p4()[idx_probe], trigObjs_HLT_DoubleMu7);
    }
    else {
        passTrigger_    = 1;
    }

    //
    // event quantities
    //
    type_ = 1;
    LorentzVector vec_TP = cms2.mus_p4()[idx_tag] + cms2.mus_p4()[idx_probe];
    mass_ = sqrt(fabs(vec_TP.mass2()));

    //
    // tag variables
    //
    pt1_		= cms2.mus_p4()[idx_tag].Pt();
    eta1_		= cms2.mus_p4()[idx_tag].Eta();
    phi1_		= cms2.mus_p4()[idx_tag].Phi();
    reliso1_		= muonIsoValue(idx_tag);
    relisont1_		= muonIsoValue(idx_tag,false);
    reliso1_fastjet_	= muonIsoValue_FastJet(idx_tag);

    //
    // probe variables
    //
    pt2_		= cms2.mus_p4()[idx_probe].Pt();
    eta2_		= cms2.mus_p4()[idx_probe].Eta();
    phi2_		= cms2.mus_p4()[idx_probe].Phi();
    reliso2_		= muonIsoValue(idx_probe);
    relisont2_		= muonIsoValue(idx_probe,false);
    reliso2_fastjet_	= muonIsoValue_FastJet(idx_probe);

    vtxweight_          = vtxweight();

    // number of good vertices
    int ngoodvertex = 0;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v)
    {
      if(isGoodVertex(v)) ngoodvertex++;
    }

    // number of good DA vertices
    int ngooddavertex = 0;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v)
    {
      if(isGoodDAVertex(v)) ngooddavertex++;
    }

    rho_        = cms2.evt_rho();
    goodvtx_ 	= ngoodvertex; //isGoodVertex
    gooddavtx_ 	= ngooddavertex; //isGoodDAVertex

    //
    // fill the baby tree
    //
    babyTree_->Fill();

}

