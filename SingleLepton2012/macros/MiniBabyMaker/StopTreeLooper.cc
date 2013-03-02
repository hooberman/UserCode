//#include "../../CORE/jetSmearingTools.h"
#include "../../CORE/Thrust.h"
#include "../../CORE/EventShape.h"

#include "StopTreeLooper.h"

//#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "TFitter.h"
#include "TRandom3.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace Stop;

std::set<DorkyEventIdentifier> already_seen; 
std::set<DorkyEventIdentifier> events_lasercalib; 
std::set<DorkyEventIdentifier> events_hcallasercalib; 

StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
  // t1metphicorr = -9999.;
  // t1metphicorrphi = -9999.;
  // t1metphicorrmt = -9999.;
  // min_mtpeak = -9999.;
  // max_mtpeak = -9999.; 
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;
  jets.clear();
  jets_btag.clear();
  jets_sigma.clear();
  metphi  = -99999.;
}

void StopTreeLooper::loop(TChain *chain, TString name)
{
  TRandom3 aleatorio;

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents("../Core/badlaser_events.txt", events_lasercalib);
  load_badlaserevents("../Core/badhcallaser_events.txt", events_hcallasercalib);

  //---------------------------------
  // check for valid chain
  //---------------------------------

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  if (listOfFiles->GetEntries() == 0) {
    cout << "[StopTreeLooper::loop] no files in chain" << endl;
    return;
  }

  //---------------------------------
  // set up histograms
  //---------------------------------

  gROOT->cd();

  makeTree(name.Data(), chain);

  TH2F *h_nsig, *h_nsig25, *h_nsig75 ;

  if( name.Contains("T2") ){
    char* h_nsig_filename = "";

    if( name.Contains("T2tt") )
      h_nsig_filename = "/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2tt_3/myMassDB_T2tt.root";

    if( name.Contains("T2bw") )
      h_nsig_filename = "/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_3/myMassDB_T2bw_fine.root";

    cout << "[StopTreeLooper::loop] opening mass TH2 file " << h_nsig_filename << endl;

    TFile *f_nsig = TFile::Open(h_nsig_filename);

    assert(f_nsig);

    h_nsig = (TH2F*) f_nsig->Get("masses");

    if( name.Contains("T2bw") ){
      h_nsig25 = (TH2F*) f_nsig->Get("masses25");
      h_nsig75 = (TH2F*) f_nsig->Get("masses75");
    }

  }

  const int NREG = 5;
  TMVA::Reader* reader[NREG];
  if ( __apply_mva ) {
    for (int i=1; i < NREG ; i++){
      reader[i] = new TMVA::Reader( "!Color:!Silent" );
      reader[i]->AddVariable("mini_met", &met_);
      reader[i]->AddVariable("mini_chi2minprob", &chi2_);
      reader[i]->AddVariable("mini_mt2wmin", &mt2w_);
      reader[i]->AddVariable("mini_htssm/(mini_htosm+mini_htssm)", &htratiom_);
      reader[i]->AddVariable("mini_dphimjmin", &dphimjmin_);
      //		  reader[i]->AddVariable("mini_pt_b", &pt_b_);
      //		  reader[i]->AddVariable("mini_lep1pt", &lep1pt_);
      //		  reader[i]->AddVariable("mini_nb", &nb_);

      TString dir, prefix;
      dir    = "/home/users/magania/stop/SingleLepton2012/MVA/weights/";
      prefix = "classification_T2tt_";
      prefix += i;
      prefix += "_BDT";

      TString weightfile = dir + prefix + TString(".weights.xml");
      reader[i]->BookMVA( "BDT" , weightfile );
    }
  }

  // TFile* vtx_file = TFile::Open("vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root");
  // if( vtx_file == 0 ){
  //   cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
  //   exit(0);
  // }

  // TH1F* h_vtx_wgt = (TH1F*)vtx_file->Get("hratio");
  // h_vtx_wgt->SetName("h_vtx_wgt");

  //
  // file loop
  //

  unsigned int nEventsPass=0;
  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  ULong64_t nEventsTotal = 0;
  //  int i_permille_old = 0;

  bool isData = name.Contains("data") ? true : false;

  cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

  cout << "[StopTreeLooper::loop] N jets requirement " << NJETS_CUT
       << " with pT>"<<JET_PT<<" GeV and |eta|<" << JET_ETA
       << endl;

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

    //---------------------------------
    // load the stop baby tree
    //---------------------------------
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    stopt.Init(tree);

    tree->CopyAddresses(outTree_);

    //---------------------------------
    // event loop
    //---------------------------------

    ULong64_t nEvents = tree->GetEntries();
    //nEvents = 1000;

    for(ULong64_t event = 0; event < nEvents; ++event) {
      stopt.GetEntry(event);
      tree->GetEntry(event);

      //---------------------------------
      // increment counters
      //---------------------------------

      ++nEventsTotal;
      if (nEventsTotal%10000==0) {
	ULong64_t i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	//if (i_permille != i_permille_old) {//this prints too often!
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	  fflush(stdout);
	}
	//	i_permille_old = i_permille;
      }

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {

	DorkyEventIdentifier id = {stopt.run(), stopt.event(), stopt.lumi() };

	if (is_duplicate(id, already_seen) ) continue;
	if (is_badLaserEvent(id,events_lasercalib) ) continue;
	if (is_badLaserEvent(id,events_hcallasercalib) ) continue;

      }

      //------------------------------------------ 
      // init baby
      //------------------------------------------ 

      initBaby(); // set all branches to -1

      //------------------------------------------ 
      // event weight
      //------------------------------------------ 
      weight_ = isData ? 1. : 
	( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );

      if( name.Contains("T2tt") ) {
	int bin = h_nsig->FindBin(stopt.mg(),stopt.ml());
	float nevents = h_nsig->GetBinContent(bin);
	//NOTE::need to add vtx. reweighting for the signal sample
	weight_  = stopt.xsecsusy() * 1000.0 / nevents * 19.5; 

	// cout << "mg " << tree->mg_ << " ml " << tree->ml_ 
	//      << " bin " << bin << " nevents " << nevents 
	//      << " xsec " << tree->xsecsusy_ 
	//      << " weight " << evtweight << endl;
      }

      if( name.Contains("T2bw") ) {
	float nevents = 0;
	if ( stopt.x() == 0.25 ){
	  int bin = h_nsig25->FindBin(stopt.mg(),stopt.ml());
	  nevents = h_nsig25->GetBinContent(bin);
	} else if ( stopt.x() == 0.50 ){
	  int bin = h_nsig->FindBin(stopt.mg(),stopt.ml());
	  nevents = h_nsig->GetBinContent(bin);
	} else if ( stopt.x() == 0.75 ){
	  int bin = h_nsig75->FindBin(stopt.mg(),stopt.ml());
	  nevents = h_nsig75->GetBinContent(bin);
	}

	assert ( nevents > 0 );

	//NOTE::need to add vtx. reweighting for the signal sample
	weight_ =  stopt.xsecsusy() * 1000.0 / nevents * 19.5; 

	// cout << "mg " << stopt.mg() << " ml " << stopt.ml() 
	//      << " x " << stopt.x() 
	//      << " nevents " << nevents 
	//      << " xsec " << stopt.xsecsusy() 
	//      << " weight " << evtweight << endl;
      }

      sltrigeff_   = isData ? 1. : 
	getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
      dltrigeff_ = isData ? 1. : 
	getdltrigweight(stopt.id1(), stopt.id2());
	    
      //------------------------------------------ 
      // selection criteria
      //------------------------------------------ 
      // >=1 lepton, rho cut, MET filters, veto 2 nearby leptons
      if ( !passEvtSelection(name) ) continue; 

      // MET > 100 GeV
      met_ = stopt.t1metphicorr();
      metphi = stopt.t1metphicorrphi();
      if ( met_ < 100.0 ) continue; 

      // jet selection
      jets.clear();
      jets_btag.clear();
      jets_sigma.clear();
      njets_ = 0;
      nb_ = 0;

      // kinematic variables
      htssl_ = 0.;
      htosl_ = 0.;
      htssm_ = 0.;
      htosm_ = 0.;

      for (unsigned int i =0; i<stopt.pfjets().size(); i++){
	if ( stopt.pfjets().at(i).pt() < JET_PT ) continue;
	if ( fabs(stopt.pfjets().at(i).eta()) > JET_ETA ) continue;
	if ( stopt.pfjets_beta2_0p5().at(i) < 0.2 ) continue;

	// jet information
	njets_++;
	if (njets_==1) pt_J1_ = stopt.pfjets().at(i).pt();
	if (njets_==2) pt_J2_ = stopt.pfjets().at(i).pt();
	jets.push_back(stopt.pfjets().at(i));
	jets_sigma.push_back(stopt.pfjets_sigma().at(i));

	// b-tagging information
	//NOTE::need to reweight b-tagging discriminator
	if (stopt.pfjets_csv().at(i) > 0.679) {
	  nb_++;
	  // leading b-jet variables
	  if (nb_==1) {
	    pt_b_ = stopt.pfjets().at(i).pt();
	    dRleptB1_ = deltaR(stopt.pfjets().at(i).eta() , 
			       stopt.pfjets().at(i).phi() , 
			       stopt.lep1().eta(), stopt.lep1().phi());
	  }
	}
	jets_btag.push_back(stopt.pfjets_csv().at(i));
		
	float dPhiL = getdphi( stopt.lep1().Phi(), stopt.pfjets().at(i).phi() );
	float dPhiM = getdphi(             metphi, stopt.pfjets().at(i).phi() );

	if( dPhiL  < (TMath::Pi()/2) ) htssl_=htssl_+stopt.pfjets().at(i).pt();
	if( dPhiL >= (TMath::Pi()/2) ) htosl_=htosl_+stopt.pfjets().at(i).pt();
	if( dPhiM  < (TMath::Pi()/2) ) htssm_=htssm_+stopt.pfjets().at(i).pt();
	if( dPhiM >= (TMath::Pi()/2) ) htosm_=htosm_+stopt.pfjets().at(i).pt();

      }

      if ( njets_ < NJETS_CUT ) continue; 

      // event shapes: HT in hemispheres
      htratiol_   = htssl_ / (htosl_ + htssl_);
      htratiom_   = htssm_ / (htosm_ + htssm_);

      //angles
      dphimj1_    = njets_>0 ? getdphi(metphi, jets.at(0).phi() ) : -99999.;
      dphimj2_    = njets_>1 ? getdphi(metphi, jets.at(1).phi() ) : -99999.;
      dphimjmin_  = njets_>1 ?  TMath::Min( dphimj1_ , dphimj2_ ) : -99999.;

      //------------------------------------------ 
      // kinematic variables
      //------------------------------------------ 

      mt2w_ = (float)calculateMT2w(jets, jets_btag, stopt.lep1(), met_, metphi);
      chi2_ = (float)calculateChi2SNT(jets, jets_sigma, jets_btag);

      //------------------------------------------ 
      // datasets bit
      //------------------------------------------ 

      bool dataset_1l=false;

      if((isData) && name.Contains("muo") 
	 && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
      if((isData) && name.Contains("ele") 
	 && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

      if(!isData) dataset_1l=true;

      bool dataset_CR4=false;

      if((isData) && name.Contains("dimu") 
	 && (abs(stopt.id1()) == 13 ) && (abs(stopt.id2())==13) 
	 && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
      if((isData) && name.Contains("diel") 
	 && (abs(stopt.id1()) == 11 ) && (abs(stopt.id2())==11) 
	 && (fabs( stopt.dilmass() - 91.) > 15.)) dataset_CR4=true;
      if((isData) && name.Contains("mueg") 
	 && abs(stopt.id1()) != abs(stopt.id2())) dataset_CR4=true;

      if(!isData) dataset_CR4=true;

      //------------------------------------------ 
      // variables to add to baby
      //------------------------------------------ 

      // isolated track veto selection
      passisotrk_   = passIsoTrkVeto_v4() ? 1 : 0; 

      // which selections are passed
      // pass signal region preselection
      sig_        = ( dataset_1l 
		      && passOneLeptonSelection(isData) 
		      && passisotrk_ == 1
		      && nb_>=1 ) ? 1 : 0; 
      // pass CR1 (b-veto) control region preselection
      cr1_        = ( dataset_1l 
		      && passOneLeptonSelection(isData) 
		      && passisotrk_ == 1
		      && nb_==0 ) ? 1 : 0; 
      // pass CR4 (dilepton) control region preselection
      cr4_        = ( dataset_CR4 
		      && passDileptonSelection(isData) 
		      && nb_>=1) ? 1 : 0; 
      // pass CR1 (lepton+isotrack) control region preselection
      cr5_        = ( dataset_1l 
		      && passOneLeptonSelection(isData) 
		      && passisotrk_ == 0
		      && nb_>=1) ? 1 : 0; 

      // Cut & Count analysis Selections
      t2ttLM_     = pass_T2tt_LM(isData);
      t2ttHM_     = pass_T2tt_HM(isData);

      // number of analysis selected leptons
      nlep_       = stopt.ngoodlep();              
      lep1pt_     = stopt.lep1().pt();             // 1st lepton pt
      lep1eta_    = fabs( stopt.lep1().eta() );    // 1st lepton eta

      if( nlep_ > 1 ){
	lep2pt_    = stopt.lep1().pt();            // 2nd lepton pt
	lep2eta_   = stopt.lep1().eta();           // 2nd lepton eta
	dilmass_   = stopt.dilmass();              // dilepton mass
      }

      rand_       = aleatorio.Uniform(0.0,1.0);

      for ( int i=0; i < NREG; i++){
	float bdtval = 0;
	if ( __apply_mva && i>0 )
	  bdtval = reader[i]->EvaluateMVA( "BDT" );
	bdt_.push_back(bdtval);
      }

      // fill me up
      nEventsPass++;
      outTree_->Fill();

    } // end event loop

    delete tree;
  } // end file loop

  //-------------------------
  // finish and clean up
  //-------------------------

  cout << "[StopTreeLooper::loop] saving mini-baby with total entries " << nEventsPass << endl;

  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
  delete outFile_;

  already_seen.clear();

  gROOT->cd();

}

//--------------------------------------------
// create the tree and set branch addresses
//--------------------------------------------

void StopTreeLooper::makeTree(const char *prefix, TChain* chain){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  string revision = "$Revision: 1.32 $";
  string revision_no = revision.substr(11, revision.length() - 13);
  outFile_   = new TFile(Form("output/%s_mini_%s.root",prefix,revision_no.c_str()), "RECREATE");
  outFile_->cd();

  outTree_ = chain->CloneTree(0);

  outTree_->Branch("mini_met"       , &met_       ,  "mini_met/F"	 );
  				      
  outTree_->Branch("mini_sig"       , &sig_       ,  "mini_sig/I"	 );
  outTree_->Branch("mini_cr1"       , &cr1_       ,  "mini_cr1/I"	 );
  outTree_->Branch("mini_cr4"       , &cr4_       ,  "mini_cr4/I"	 );
  outTree_->Branch("mini_cr5"       , &cr5_       ,  "mini_cr5/I"	 );
  				      
  outTree_->Branch("mini_chi2"      , &chi2_      ,  "mini_chi2/F"      );
  outTree_->Branch("mini_mt2w"      , &mt2w_      ,  "mini_mt2w/F"      );
  				      
  outTree_->Branch("mini_weight"    , &weight_    ,  "mini_weight/F"	 );
  outTree_->Branch("mini_sltrigeff" , &sltrigeff_ ,  "mini_sltrigeff/F" );
  outTree_->Branch("mini_dltrigeff" , &dltrigeff_ ,  "mini_dltrigeff/F" );
  				      
  outTree_->Branch("mini_nb"        , &nb_        ,  "mini_nb/I"	 );
  outTree_->Branch("mini_njets"     , &njets_     ,  "mini_njets/I"	 );
  				      
  outTree_->Branch("mini_passisotrk", &passisotrk_,  "mini_passisotrk/I");
  				      
  outTree_->Branch("mini_nlep"      , &nlep_      ,  "mini_nlep/I"	 );
  outTree_->Branch("mini_dRleptB1"  , &dRleptB1_  ,  "mini_dRleptB1/F"	 );
  outTree_->Branch("mini_lep1pt"    , &lep1pt_    ,  "mini_lep1pt/F"	 );
  outTree_->Branch("mini_lep1eta"   , &lep1eta_   ,  "mini_lep1eta/F"	 );
  outTree_->Branch("mini_lep2pt"    , &lep2pt_    ,  "mini_lep2pt/F"	 );
  outTree_->Branch("mini_lep2eta"   , &lep2eta_   ,  "mini_lep2eta/F"	 );
  outTree_->Branch("mini_dilmass"   , &dilmass_   ,  "mini_dilmass/F"	 );
  				      
  outTree_->Branch("mini_htssl"     , &htssl_     ,  "mini_htssl/F"     );
  outTree_->Branch("mini_htosl"     , &htosl_     ,  "mini_htosl/F"     );
  outTree_->Branch("mini_htratiol"  , &htratiol_  ,  "mini_htraiol/F"   );
  outTree_->Branch("mini_htssm"     , &htssm_     ,  "mini_htssm/F"     );
  outTree_->Branch("mini_htosm"     , &htosm_     ,  "mini_htosm/F"     );
  outTree_->Branch("mini_htratiom"  , &htratiom_  ,  "mini_htraiom/F"   );
  outTree_->Branch("mini_dphimj1"   , &dphimj1_   ,  "mini_dphimj1/F"   );
  outTree_->Branch("mini_dphimj2"   , &dphimj2_   ,  "mini_dphimj2/F"   );
  outTree_->Branch("mini_dphimjmin" , &dphimjmin_ ,  "mini_dphimjmin/F" );
  				      
  outTree_->Branch("mini_pt_b"      , &pt_b_      ,  "mini_pt_b/F"      );
  outTree_->Branch("mini_pt_J1"     , &pt_J1_     ,  "mini_pt_J1/F"     );
  outTree_->Branch("mini_pt_J2"     , &pt_J2_     ,  "mini_pt_J2/F"     );
  				      
  outTree_->Branch("mini_rand"      , &rand_      ,  "mini_rand/F"      );
  				      
  outTree_->Branch("mini_t2ttLM"    , &t2ttLM_    ,  "mini_t2ttLM/F"    );
  outTree_->Branch("mini_t2ttHM"    , &t2ttHM_    ,  "mini_t2ttHM/F"    );
  	
  outTree_->Branch("mini_bdt"       , "std::vector<float>" , &bdt_     );

}

//--------------------------------------------
// set all branches to -1
//--------------------------------------------

void StopTreeLooper::initBaby(){

  // which selections are passed
  sig_        = -1;
  cr1_        = -1;
  cr4_        = -1;
  cr5_        = -1; 

  // kinematic variables
  met_        = -1.0;
  chi2_       = -1.0;
  mt2w_       = -1.0;

  // event shapes
  htssl_      = -1.0;
  htosl_      = -1.0;
  htratiol_   = -1.0;
  htssm_      = -1.0;
  htosm_      = -1.0;
  htratiom_   = -1.0;

  // weights
  weight_     = -1.0;
  sltrigeff_  = -1.0;
  dltrigeff_  = -1.0;

  // jet counting
  nb_         = -1;
  njets_      = -1;

  // lepton variables
  passisotrk_ = -1;

  nlep_       = -1;
  dRleptB1_   = -1;
  lep1pt_     = -1.0;
  lep1eta_    = -1.0;
  lep2pt_     = -1.0;
  lep2eta_    = -1.0;
  dilmass_    = -1.0;

  // jet kinematics
  pt_b_  = -1.0;
  pt_J1_ = -1.0;
  pt_J2_ = -1.0;

  rand_  = -1.0;

  t2ttLM_ =-1;
  t2ttLM_ =-1;

  bdt_.clear();
}

