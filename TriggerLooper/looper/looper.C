#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/LorentzVector.h"
#include "histtools.h"
#include "looper.h"
#include "TTreeCache.h"
#include "../CORE/CMS2.h"
#include "../CORE/metSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/mcSelections.h"
#include "../CORE/muonSelections.h"
#include "../Tools/goodrun.cc"
#include "../CORE/utilities.cc"
#include "../CORE/ttbarSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/mcSUSYkfactor.h"
#include "../CORE/triggerSuperModel.h"
#include "../CORE/triggerUtils.h"
#include "../Tools/vtxreweight.cc"
#include "../Tools/msugraCrossSection.cc"

bool verbose        = false;
bool doTenPercent   = false;

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

//--------------------------------------------------------------------

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

//--------------------------------------------------------------------

int findTriggerIndex(TString trigName)
{
    vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
    vector<TString>::const_iterator end_it = hlt_trigNames().end();
    vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
    if(found_it != end_it) return found_it - begin_it;
    return -1;
}

//--------------------------------------------------------------------

bool objectPassTrigger(const LorentzVector &obj, const std::vector<LorentzVector> &trigObjs, float pt) 
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

//--------------------------------------------------------------------

looper::looper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
}

//--------------------------------------------------------------------

bool passesCaloJetID (const LorentzVector &jetp4)
{
  int jet_idx = -1;
  double minDR = 999;

  for (unsigned int i = 0; i < cms2.jets_p4().size(); i++)
    {
      double deltaR = ROOT::Math::VectorUtil::DeltaR(jetp4, cms2.jets_p4()[i]);

      if (deltaR < minDR)
	{
	  minDR = deltaR;
	  jet_idx = i;
	}
    }

  if (jet_idx < 0)
    return false;

  if (cms2.jets_emFrac()[jet_idx] < 0.01 || cms2.jets_fHPD()[jet_idx] > 0.98 || cms2.jets_n90Hits()[jet_idx] < 2)
    return false;

  return true;
}

//--------------------------------------------------------------------

bool passesPFJetID(unsigned int pfJetIdx) {

  float pfjet_chf_  = cms2.pfjets_chargedHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  float pfjet_cef_  = cms2.pfjets_chargedEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
  int   pfjet_cm_   = cms2.pfjets_chargedMultiplicity()[pfJetIdx];
  int   pfjet_mult_ = pfjet_cm_ + cms2.pfjets_neutralMultiplicity()[pfJetIdx] + cms2.pfjets_muonMultiplicity()[pfJetIdx];

  if (pfjet_nef_ >= 0.99)
    return false;
  if (pfjet_nhf_ >= 0.99)
    return false;
  if (pfjet_mult_ < 2)
    return false;

  if (fabs(cms2.pfjets_p4()[pfJetIdx].eta()) < 2.4)
    {
      if (pfjet_chf_ < 1e-6)
	return false;
      if (pfjet_cm_ < 1)
	return false;
      if (pfjet_cef_ >= 0.99)
	return false;
    }

  return true;
}  

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

void looper::InitBaby(){

  pjet1_ = 0;
  pjet2_ = 0;
  pjet3_ = 0;
  pjet4_ = 0;
  cjet1_ = 0;
  cjet2_ = 0;
  cjet3_ = 0;
  cjet4_ = 0;

  run_     = -999;
  event_   = -999;
  lumi_    = -999;
}

//--------------------------------------------------------------------

void looper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

int looper::ScanChain(TChain* chain, char *prefix){

  set_goodrun_file( g_json );

  bool isData = true;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  BookHistos(prefix);
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  int nSkip_els_conv_dist = 0;

  float netot  = 0.;
  float nmtot  = 0.;
  float nepass = 0.;
  float nmpass = 0.;

  if(g_createTree) makeTree(prefix);

  bool hasJptBtagBranch = true;

  char* thisFile = "blah";

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    if( !f || f->IsZombie() ) {
      cout << "Skipping bad input file: " << currentFile->GetTitle() << endl;
      continue; //exit(1);                                                                                             
    }

    if( strcmp(thisFile,currentFile->GetTitle()) != 0 ){
      thisFile = (char*) currentFile->GetTitle();
      cout << thisFile << endl;
    }

    TTree *tree = (TTree*)f->Get("Events");

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      if( doTenPercent ){
	if( !(nEventsTotal%10==0) ) continue;
      }

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //Matevz
      tree->LoadTree(z);

      cms2.GetEntry(z);

      InitBaby();

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset() << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      if( !cleaning_goodDAVertexApril2011() )                        continue;
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------
 
      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
          skipEvent = true;
        }
      }
             
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }
   
      //---------------------------------------------
      // find leptons passing analysis selection
      //---------------------------------------------

      VofP4 goodLeptons;
      vector<int> lepId;
      vector<int> lepIndex;

      ngoodlep_ = 0;
      ngoodel_  = 0;
      ngoodmu_  = 0;
            
      for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
	if( els_p4().at(iel).pt() < 10 )                                              continue;
	if( !pass_electronSelection( iel , electronSelection_ssV5 , false , false ) ) continue;
	goodLeptons.push_back( els_p4().at(iel) );
	lepId.push_back( els_charge().at(iel) * 11 );
	lepIndex.push_back(iel);
	ngoodel_++;
	ngoodlep_++;
      }

          
      for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
	if( mus_p4().at(imu).pt() < 10 )           continue;
	if( !muonId( imu , OSGeneric_v3 ))         continue;
	goodLeptons.push_back( mus_p4().at(imu) );
	lepId.push_back( mus_charge().at(imu) * 13 );
	lepIndex.push_back(imu);
	ngoodmu_++;
	ngoodlep_++;
      }  
      
      //------------------------------------------------
      // trigger study: turn-on curve for jet triggers
      //------------------------------------------------

      /*
      trgjet_  = 0;
      passtrg_ = -1;

      if( doTriggerStudy ){

	// only consider muon events
	if( abs(id1_) == 11 ) continue;

	// run range for which HLT_Mu17_CentralJet30 is un-prescaled
	if( evt_run() < 160329 || evt_run() > 163261 ) continue;

	// require event passes HLT_Mu20_v1
	if( !passUnprescaledHLTTrigger("HLT_Mu20_v1") ) continue;
 
	// require found trigger object
	std::vector<LorentzVector> muonObj;
	muonObj = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu20_v1")];

	if( muonObj.size() == 0 ) continue;

	// require offline muon matched to muon trigger object (dR < 0.1)
	if( dRbetweenVectors( *lep1_ , muonObj.at(0) ) > 0.1 ) continue;

	// now get trigger objects for HLT_Mu20_CentralJet30
	std::vector<LorentzVector> muonJetObj;
	if( evt_run() >= 160329 && evt_run() <= 161176 )
	  muonJetObj = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu17_CentralJet30_v1")];
	else if( evt_run() >= 161210 && evt_run() <= 163261 )
	  muonJetObj = cms2.hlt_trigObjs_p4()[findTriggerIndex("HLT_Mu17_CentralJet30_v2")];

	// did event pass HLT_Mu17_CentralJet30?
	passtrg_ = 0;
	if( evt_run() >= 160329 && evt_run() <= 161176 )
	  passtrg_ = passUnprescaledHLTTrigger("HLT_Mu17_CentralJet30_v1") ? 1 : 0;
	else if( evt_run() >= 161210 && evt_run() <= 163261 )
	  passtrg_ = passUnprescaledHLTTrigger("HLT_Mu17_CentralJet30_v2") ? 1 : 0;

	// now find highest pt offline jet matched to HLT jet
	int offlinejet = -1;
	int njets      =  0;
	float maxpt    = -1;

        for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
          
          LorentzVector vjet = jets_p4().at(ijet) * jets_corL1FastL2L3().at(ijet);

	  if( dRbetweenVectors( vjet , *lep1_ ) < 0.4 ) continue;
	  //if( vjet.pt() < 20. )                                          continue;

	  // if event passed trigger, require jet is matched to HLT object
	  if( passtrg_ == 1 ){
	    
	    bool HLTmatch = false;

	    for( unsigned int ihlt = 0 ; ihlt < muonJetObj.size() ; ++ihlt ){
	      
	      // exclude HLT muon object
	      if( dRbetweenVectors( vjet , muonObj.at(0) ) < 0.1 )       continue;
	    
	      // dr match to HLT jet object
	      if( dRbetweenVectors( vjet , muonJetObj.at(ihlt) ) > 0.4 ) continue;

	      HLTmatch = true;
	    }

	    if( !HLTmatch ) continue;
	  }

	  njets++;

	  if( vjet.pt() > maxpt ){
	    maxpt      = vjet.pt();
	    offlinejet = ijet;
	  }

        }

	//if( passtrg_ == 0 ) cout << "FAIL TRIGGER" << endl;
	if( offlinejet < 0 ) continue;
	//if( njets > 1      ) continue;

	// now store offline jet p4 and whether muon-jet trigger passed
	trgjet_  = &(jets_p4().at(offlinejet) * jets_corL1FastL2L3().at(offlinejet));

	outTree->Fill();
	continue;

	// cout << "muonJetObj.size() " << muonJetObj.size() << endl;
	// cout << "muon trigger       : pt, eta, phi    " << muonObj.at(0).pt() << ", " << muonObj.at(0).eta() << ", " << muonObj.at(0).phi() << endl;

	// if( muonJetObj.size() > 1 ){
	//   cout << "muon-jet trigger 1 : pt, eta, phi    " << muonJetObj.at(0).pt() << ", " << muonJetObj.at(0).eta() << ", " << muonJetObj.at(0).phi() << endl;
	//   cout << "muon-jet trigger 2 : pt, eta, phi    " << muonJetObj.at(1).pt() << ", " << muonJetObj.at(1).eta() << ", " << muonJetObj.at(1).phi() << endl;

	//
	//   float dr2 = dRbetweenVectors( muonObj.at(0) , muonJetObj.at(1) );

	//   cout << "dr1 " << dr1 << " dr2 " << dr2 << endl;
	// }



	//cout << endl << endl;
	//cout << "trigger objects " << muonObj.size() << endl;
	//if( muonObj.size() > 0 ) cout << "pt, eta, phi    " << muonObj.at(0).pt() << ", " << muonObj.at(0).eta() << ", " << muonObj.at(0).phi() << endl;
	//cout << "muon " << lep1_->pt() << ", " << lep1_->eta() << ", " << lep1_->phi() << endl;

      }
      */

      //-------------------------------------
      // jet counting
      //-------------------------------------

      VofP4 vpfjets_p4;

      njets_ = 0.;
      ht_    = 0.;

      int   imaxjet   = -1;
      float maxjetpt  = -1.;

      vector<int> goodjets;
      goodjets.clear();

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

	LorentzVector vjet      = pfjets_corL1FastL2L3().at(ijet) * pfjets_p4().at(ijet);

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;
          
	if( !passesPFJetID(ijet) )     continue;
	if( fabs( vjet.eta() ) > 3.0 ) continue;
	if( vjet.pt() < 30 )           continue;

	njets_++;
	ht_ += vjet.pt();

	vpfjets_p4.push_back( vjet );
      }

      sort( vpfjets_p4.begin(), vpfjets_p4.end(), sortByPt);

      if( njets_ > 0 ) 	pjet1_ = &( vpfjets_p4.at(0) );
      if( njets_ > 1 ) 	pjet2_ = &( vpfjets_p4.at(1) );
      if( njets_ > 2 ) 	pjet3_ = &( vpfjets_p4.at(2) );
      if( njets_ > 3 ) 	pjet4_ = &( vpfjets_p4.at(3) );

      //------------------------------------------
      // count calojets
      //------------------------------------------

      VofP4 vcalojets_p4;
      ncjets_ = 0;
      htc_    = 0.0;

      for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
	
	LorentzVector vjet = jets_p4().at(ijet) * jets_corL1FastL2L3().at(ijet);

	if( !passesCaloJetID( vjet ) )         continue;	
	if( fabs( vjet.eta() ) > 3.0 )         continue;
	if( vjet.pt() < 30           )         continue;

	bool rejectJet = false;
	for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
	  if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
	}
	if( rejectJet ) continue;

	htc_ += vjet.pt();
	ncjets_ ++;
	
	vcalojets_p4.push_back(vjet);
      }

      sort( vcalojets_p4.begin(), vcalojets_p4.end(), sortByPt);

      if( ncjets_ > 0 ) 	cjet1_ = &( vcalojets_p4.at(0) );
      if( ncjets_ > 1 ) 	cjet2_ = &( vcalojets_p4.at(1) );
      if( ncjets_ > 2 ) 	cjet3_ = &( vcalojets_p4.at(2) );
      if( ncjets_ > 3 ) 	cjet4_ = &( vcalojets_p4.at(3) );

      //----------------------------
      // MET flavors
      //----------------------------

      pfmet_    = evt_pfmet();
      pfmetphi_ = evt_pfmetPhi();
      pfsumet_  = evt_pfsumet();

      //----------------------------------------
      // nvertex variables
      //----------------------------------------

      nvtx_ = 0;
    
      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
	if(isGoodVertex(v)) ++nvtx_;
      }

      ndavtx_ = 0;
    
      for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
	if(isGoodDAVertex(v)) ++ndavtx_;
      }

      strcpy(dataset_, cms2.evt_dataset().Data());  //dataset name
      run_          = evt_run();                    //run
      lumi_         = evt_lumiBlock();              //lumi
      event_        = evt_event();                  //event

      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  return 0;

}

//--------------------------------------------------------------------
 
void looper::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()

//--------------------------------------------------------------------

void looper::makeTree(char *prefix ){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  outFile   = new TFile(Form("../output/%s/%s.root",g_version,prefix), "RECREATE");
  //outFile   = new TFile("temp.root","RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in looper.h
  outTree->Branch("njets",           &njets_,            "njets/I");
  outTree->Branch("ncjets",          &ncjets_,           "ncjets/I");
  outTree->Branch("ht",              &ht_,               "ht/F");
  outTree->Branch("htc",             &htc_,              "htc/F");
  outTree->Branch("pfmet",           &pfmet_,            "pfmet/F");
  outTree->Branch("pfmetphi",        &pfmetphi_,         "pfmetphi/F");
  outTree->Branch("pfsumet",         &pfsumet_,          "pfsumet/F");
  outTree->Branch("dataset",         &dataset_,          "dataset[200]/C");
  outTree->Branch("run",             &run_,              "run/I");
  outTree->Branch("lumi",            &lumi_,             "lumi/I");
  outTree->Branch("event",           &event_,            "event/I");
  outTree->Branch("ngoodlep",        &ngoodlep_,         "ngoodlep/I");
  outTree->Branch("ngoodel",         &ngoodel_,          "ngoodel/I");
  outTree->Branch("nvtx",            &nvtx_,             "nvtx/I");
  outTree->Branch("ndavtx",          &ndavtx_,           "ndavtx/I");
  outTree->Branch("cjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet1_	);
  outTree->Branch("cjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet2_	);
  outTree->Branch("cjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet3_	);
  outTree->Branch("cjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &cjet4_	);
  outTree->Branch("pjet1"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet1_	);
  outTree->Branch("pjet2"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet2_	);
  outTree->Branch("pjet3"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet3_	);
  outTree->Branch("pjet4"    , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &pjet4_	);

}

//--------------------------------------------------------------------

vector<int> goodDAVertices(){

  vector<int> myVertices;
  myVertices.clear();
  
  for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
    if( !isGoodDAVertex(v) ) continue;
    myVertices.push_back(v);
  }
  
  return myVertices;
}

//--------------------------------------------------------------------

float looper::dz_trk_vtx( const unsigned int trkidx, const unsigned int vtxidx ){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}
