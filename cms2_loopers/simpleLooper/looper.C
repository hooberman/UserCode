#include "looper.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <math.h>
#include <fstream>


#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>

#include "CORE/CMS2.h"
#include "CORE/trackSelections.h"
#include "CORE/metSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/ttbarSelections.cc"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"
#include "CORE/jetSelections.h"
#include "CORE/ttbarSelections.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

char* iter          = "default";
bool makebaby       = true;
bool debug          = false;
bool calculateTCMET = false;
float lumi          = 1.;

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

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
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

using namespace tas;
void looper::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  bookHistos();

  set_goodrun_file("Cert_TopAug25_Merged_135059-143336_goodruns.txt");
  ofile.open( Form( "output/%s_%s_events.txt" , prefix , iter) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  MakeBabyNtuple( Form( "output/%s_%s_baby.root" , prefix , iter) );

  if( debug ) cout << "Begin looping over files" << endl;

  float nee  = 0.;
  float nmm  = 0.;
  float nem  = 0.;
  float ntot = 0.;


  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next()))
    {
      TFile f(currentFile->GetTitle());
      TTree *tree = (TTree*)f.Get("Events");
      cms2.Init(tree);

      // event loop
      unsigned int nEvents = tree->GetEntries();

      for (unsigned int event = 0; event < nEvents; ++event)
        {
          if( debug ) cout << "Event " << event << endl;

          cms2.GetEntry(event);
          ++nEventsTotal;

          // progress feedback to user
          if (nEventsTotal % 1000 == 0)
            {
              // xterm magic from L. Vacavant and A. Cerri
              if (isatty(1))
                {
                  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                         "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
                  fflush(stdout);
                }
            }

          // skip duplicates
          if( isData ) {
            DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
            if (is_duplicate(id) ){
              continue;
            }
          }
         
          //---------------------------------------------
          // Event Selection
          //---------------------------------------------

          InitBabyNtuple();

	  //if ( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))  continue;
          //if( !cleaning_standardOctober2010() )                           continue;

          vector<unsigned int> v_goodHyps;
          v_goodHyps.clear();
          vector<unsigned int> v_vtxIndices;
          v_vtxIndices.clear();

          for(unsigned int i = 0; i < hyp_p4().size(); ++i) {
            
            //if( !passSUSYTrigger_v1( isData , hyp_type()[i] ) ) continue;
            
            //check that hyp leptons come from same vertex
            if( !hypsFromSameVtx( i ) )    continue;
            
            //OS, pt > (20,10) GeV, dilmass > 10 GeV
            if( hyp_lt_id()[i] * hyp_ll_id()[i] > 0 )  continue;
                      
            //pt > (20,10) GeV
            if( TMath::Max( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < 20. )   continue;
            if( TMath::Min( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < 10. )   continue;
            if( hyp_p4()[i].mass() < 10 )                                         continue;
            
            //cout << "Check muons" << endl;
            //muon ID
            if (abs(hyp_ll_id()[i]) == 13  && !( muonId(hyp_ll_index()[i] , NominalWWV1 ) ) )   continue;
            if (abs(hyp_lt_id()[i]) == 13  && !( muonId(hyp_lt_index()[i] , NominalWWV1 ) ) )   continue;
            
            //OSV1
            if (abs(hyp_ll_id()[i]) == 11  && !( pass_electronSelection( hyp_ll_index()[i] , electronSelection_wwV1 , false , false ))) continue;
            if (abs(hyp_lt_id()[i]) == 11  && !( pass_electronSelection( hyp_lt_index()[i] , electronSelection_wwV1 , false , false ))) continue;
            
            v_goodHyps.push_back( i );
            
          }
          
          //skip events with no good hyps
          if( v_goodHyps.size() == 0 ) continue;

          if( debug ) cout << "Pass event selection" << endl;

          //returns the index of the best hypothesis in the vector of hypotheses
          unsigned int hypIdx = selectHypByHighestSumPt(v_goodHyps);
          
          //---------------------------------------
          // jet counting
          //---------------------------------------

          njets_ = 0;

          vector<int> goodjets;
          goodjets.clear();

          int imaxjet    = -1;
          float maxjetpt = -1.;

          for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
            LorentzVector vjet = pfjets_corL1L2L3().at(ijet) * pfjets_p4().at(ijet);
            LorentzVector vlt  = hyp_lt_p4()[hypIdx];
            LorentzVector vll  = hyp_ll_p4()[hypIdx];
            
            if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
            if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
            if( fabs( vjet.eta() ) > 3.0 )           continue;
            if( !passesPFJetID(ijet) )               continue;
            if( vjet.pt() < 30. )                    continue;

            njets_++;
            goodjets.push_back(ijet);
            
            if( vjet.pt() > maxjetpt ){
              imaxjet = ijet;
              maxjetpt = vjet.pt();
            }
          }
          
          if( imaxjet > -1 ){
            jet_     = &(pfjets_corL1L2L3().at(imaxjet) * pfjets_p4().at(imaxjet));
            
          }
          
          //------------------------------------------
          // baby ntuple branches
          //------------------------------------------
      
          // event stuff
          run_    = cms2.evt_run();
          lumi_   = cms2.evt_lumiBlock();
          event_  = cms2.evt_event();

          nvtx_    = goodVertices().size();
          leptype_ = hyp_type().at(hypIdx);
          dilep_   = &hyp_p4().at(hypIdx);
          
          weight_ = 1;
          if( !isData ) weight_ = evt_scale1fb() * lumi;

          if( leptype_ == 0 ) nmm += weight_;
          if( leptype_ == 1 ) nem += weight_;
          if( leptype_ == 2 ) nem += weight_;
          if( leptype_ == 3 ) nee += weight_;
          ntot += weight_;

          
          if( hyp_ll_p4().at(hypIdx).pt() > hyp_lt_p4().at(hypIdx).pt() ){
            lep1_ = &hyp_ll_p4().at(hypIdx);
            lep2_ = &hyp_lt_p4().at(hypIdx);
          }else{
            lep1_ = &hyp_lt_p4().at(hypIdx);
            lep2_ = &hyp_ll_p4().at(hypIdx);
          }

          // pf met stuff
          pfmet_    = cms2.evt_pfmet();
          pfmetphi_ = cms2.evt_pfmetPhi();
          pfsumet_  = cms2.evt_pfsumet();

          // raw  tcmet stuff
          tcmet_    = cms2.evt_tcmet();
          tcmetphi_ = cms2.evt_tcmetPhi();
          tcsumet_  = cms2.evt_tcsumet();

          // genmet stuff
          if (!isData){
            genmet_     = cms2.gen_met();
            genmetphi_  = cms2.gen_metPhi();
            gensumet_   = cms2.gen_sumEt();
          }

          //calomet
          met_       = cms2.evt_met();
          metphi_    = cms2.evt_metPhi();
          sumet_     = cms2.evt_sumet();
      
          if( calculateTCMET ){
            
            metStruct myMetStruct = correctedTCMET( true, ofile );
            tcmetNew_    = myMetStruct.met;
            tcsumetNew_  = myMetStruct.sumet;
            tcmetphiNew_ = myMetStruct.metphi;

          }

          eventTree_->Fill();
    
        } // end loop over events
    } // end loop over files
  
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  cout << "Yields passing selection" << endl;
  cout << "ee    " << nee  << endl;
  cout << "mm    " << nmm  << endl;
  cout << "em    " << nem  << endl;
  cout << "to    " << ntot << endl;

  
  CloseBabyNtuple();

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form( "output/%s_%s_histos.root" , prefix , iter ) );
  deleteHistos();
  
} // end ScanChain

//--------------------------------------------------------------------

void looper::printEvent(  ostream& ostr ){
  ostr << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl; 
}

//--------------------------------------------------------------------

vector<int> looper::goodVertices(){

  vector<int> myVertices;
  myVertices.clear();
  
  for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
    if( !isGoodVertex(v) ) continue;
    myVertices.push_back(v);
  }
  
  return myVertices;
}

//--------------------------------------------------------------------

void looper::InitBabyNtuple ()
{
  // event stuff
  run_             = -999999;
  lumi_            = -999999;
  event_           = -999999;
  weight_          = -999999.;
  nvtx_            = -999999;
  leptype_         = -999999;
  njets_           = -999999;

  // genmet stuff
  genmet_          = -999999.;
  genmetphi_       = -999999.;
  gensumet_        = -999999.;

  // pfmet stuff
  pfmet_           = -999999.;
  pfmetphi_        = -999999.;
  pfsumet_         = -999999.;

  // calomet stuff
  met_             = -999999.;
  metphi_          = -999999.;
  sumet_           = -999999.;

  // tcmet stuff
  tcmet_           = -999999.;
  tcmetphi_        = -999999.;
  tcsumet_         = -999999.;

  // latest-and-greatest tcmet stuff
  tcmetNew_        = -999999.;
  tcmetphiNew_     = -999999.;
  tcsumetNew_      = -999999.;

  dilep_ = 0;
  jet_   = 0;
  lep1_  = 0;
  lep2_  = 0;

}

//--------------------------------------------------------------------

float looper::deltaPhi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//--------------------------------------------------------------------

void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hmumet    = new TH1F("hmumet",   "",100,0,100);
  hmujesmet = new TH1F("hmujesmet","",100,0,100);
  htcmet    = new TH1F("htcmet",   "",100,0,100);
  hrawtcmet = new TH1F("hrawtcmet","",100,0,100);

  hdmumet    = new TH1F("hdmumet",   "",400,-200,200);
  hdmujesmet = new TH1F("hdmujesmet","",400,-200,200);
  hdtcmet    = new TH1F("hdtcmet",   "",400,-200,200);
  hdrawtcmet = new TH1F("hdrawtcmet","",400,-200,200);

  hdtcmet_mumet    = new TH1F("hdtcmet_mumet","",   400,-200,200);
  hdtcmet_mujesmet = new TH1F("hdtcmet_mujesmet","",400,-200,200);

  htcmetNew    = new TH1F("htcmetNew",     "calo-tcmet",100,0,100);
  hpfmet       = new TH1F("hpfmet",        "pfmet",100,0,100);
  htcmetNewPFC = new TH1F("htcmetNewPFC",  "pfc-tcmet",100,0,100);


}

//--------------------------------------------------------------------

void looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  eventTree_ = new TTree("Events", "Events Tree");

  eventTree_->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &dilep_ );
  eventTree_->Branch("lep1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_  );
  eventTree_->Branch("lep2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_  );
  eventTree_->Branch("jet"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_   );

  eventTree_->Branch("weight"           , &weight_           , "weight/F"  );
  eventTree_->Branch("nvtx"             , &nvtx_             , "nvtx/I");
  eventTree_->Branch("leptype"          , &leptype_          , "leptype/I");
  eventTree_->Branch("njets"            , &njets_            , "njets/I");
  eventTree_->Branch("run"              , &run_              , "run/I"  );
  eventTree_->Branch("lumi"             , &lumi_             , "lumi/I" );
  eventTree_->Branch("event"            , &event_            , "event/I");

  eventTree_->Branch("genmet"           , &genmet_           , "genmet/F"   );
  eventTree_->Branch("genmetphi"        , &genmetphi_        , "genmetphi/F");
  eventTree_->Branch("gensumet"         , &gensumet_         , "gensumet/F" );
 
  eventTree_->Branch("pfmet"            , &pfmet_            , "pfmet/F"   );
  eventTree_->Branch("pfmetphi"         , &pfmetphi_         , "pfmetphi/F");
  eventTree_->Branch("pfsumet"          , &pfsumet_          , "pfsumet/F" );

  eventTree_->Branch("met"              , &met_              , "met/F"      );
  eventTree_->Branch("metphi"           , &metphi_           , "metphi/F"   );
  eventTree_->Branch("sumet"            , &sumet_            , "sumet/F"    );
 
  eventTree_->Branch("tcmet"            , &tcmet_            , "tcmet/F"      );
  eventTree_->Branch("tcmetphi"         , &tcmetphi_         , "tcmetphi/F"   );
  eventTree_->Branch("tcsumet"          , &tcsumet_          , "tcsumet/F"    );

  eventTree_->Branch("tcmetnew"         , &tcmetNew_         , "tcmetnew/F"      );
  eventTree_->Branch("tcmetphinew"      , &tcmetphiNew_      , "tcmetphinew/F"   );
  eventTree_->Branch("tcsumetnew"       , &tcsumetNew_       , "tcsumetnew/F"    );

}

//--------------------------------------------------------------------

void looper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

//--------------------------------------------------------------------

void looper::CloseBabyNtuple (){

  babyFile_->cd();
  eventTree_->Write();
  babyFile_->Close();
  
}

//--------------------------------------------------------------------
