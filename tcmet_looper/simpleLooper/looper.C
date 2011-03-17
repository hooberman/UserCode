#include "looper.h"
#include <algorithm>
#include <iostream>
#include <vector>
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
#include "CORE/susySelections.h"
#include "CORE/ttbarSelections.cc"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"
#include "CORE/jetSelections.cc"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

char* iter          = "default";
bool makebaby       = true;
bool debug          = false;
bool calculateTCMET = false;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

// inline double fround(double n, double d){
//   return floor(n * pow(10., d) + .5) / pow(10., d);
// }


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

bool jetMatchedToGenJet( LorentzVector vjet , int hypIdx){

    
  for( int ijet = 0 ; ijet < genjets_p4().size() ; ijet++ ){

    LorentzVector vgenjet = genjets_p4().at(ijet);
    LorentzVector vlt     = hyp_lt_p4()[hypIdx];
    LorentzVector vll     = hyp_ll_p4()[hypIdx];
    
    if( dRbetweenVectors(vgenjet, vll) < 0.4 )   continue;
    if( dRbetweenVectors(vgenjet, vlt) < 0.4 )   continue;
    if( fabs( vgenjet.eta() ) > 3.0 )            continue;
    if( vgenjet.pt() < 30. )                     continue;
    if( dRbetweenVectors(vgenjet, vjet) > 0.4 )  continue;

    return true;
  }

  return false;

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

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next()))
    {
      TFile f(currentFile->GetTitle());
      TTree *tree = (TTree*)f.Get("Events");
      cms2.Init(tree);

      // event loop
      unsigned int nEvents = tree->GetEntries() / 10;

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
         
          //--------------------------------------------------------------
          // Apply basic event selections
          //--------------------------------------------------------------

	  //if (!isData || goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))  continue;
          //if( !cleaning_standardOctober2010() )                          continue;
          
          //------------------------------------------
          // hyp selection
          //------------------------------------------

          vector<unsigned int> v_goodHyps;
          v_goodHyps.clear();
          vector<unsigned int> v_vtxIndices;
          v_vtxIndices.clear();

          for(unsigned int i = 0; i < hyp_p4().size(); ++i) {
            
            if( !passSUSYTrigger_v1( isData , hyp_type()[i] ) ) continue;
        
            //check that hyp leptons come from same vertex
            if( !hypsFromSameVtx( i ) )    continue;

            //OS, pt > (20,10) GeV, dilmass > 10 GeV
            if( hyp_lt_id()[i] * hyp_ll_id()[i] > 0 )  continue;
          
            //pt > (20,10) GeV
            if( TMath::Max( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < 20. )   continue;
            if( TMath::Min( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < 10. )   continue;
            if( hyp_p4()[i].mass() < 10 )                                         continue;
                    
            //muon ID
            if (abs(hyp_ll_id()[i]) == 13  && !( muonId(hyp_ll_index()[i] , OSGeneric_v1 ) ) )   continue;
            if (abs(hyp_lt_id()[i]) == 13  && !( muonId(hyp_lt_index()[i] , OSGeneric_v1 ) ) )   continue;
            
            //OSV1
            if (abs(hyp_ll_id()[i]) == 11  && !( pass_electronSelection( hyp_ll_index()[i] , electronSelection_el_OSV1 , false , false ))) continue;
            if (abs(hyp_lt_id()[i]) == 11  && !( pass_electronSelection( hyp_lt_index()[i] , electronSelection_el_OSV1 , false , false ))) continue;
            
            v_goodHyps.push_back( i );
            
          }
          
          //skip events with no good hyps
          if( v_goodHyps.size() == 0 ) continue;

          //returns the index of the best hypothesis in the vector of hypotheses
          unsigned int hypIdx = selectHypByHighestSumPt(v_goodHyps);
          int          vtxIdx = hypsFromSameVtx_int( hypIdx );

          InitBabyNtuple();

//           bool skipevent = true;
          
//           if( evt_run() == 1 && evt_lumiBlock() == 3 && evt_event() == 1393221 ) skipevent = false;
//           if( evt_run() == 1 && evt_lumiBlock() == 2 && evt_event() == 722068  ) skipevent = false;
          
//           if( skipevent ) continue;
          
//           cout << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
          
          //---------------------------------------
          // jet counting
          //---------------------------------------

          njets25_ = 0;

          vector<int> goodjets;
          goodjets.clear();

          for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
            LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
            LorentzVector vlt  = hyp_lt_p4()[hypIdx];
            LorentzVector vll  = hyp_ll_p4()[hypIdx];
            
            if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
            if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
            if( fabs( vjet.eta() ) > 3.0 )           continue;
            if( !passesPFJetID(ijet) )               continue;
            if( vjet.pt() < 30. )                    continue;

            if( !jetFromSignalPV( ijet , vtxIdx ) ) continue;
           
            njets25_++;
            goodjets.push_back(ijet);
            
//             LorentzVector pfcands_tot(0,0,0,0);

//             vector<int> pfcands = pfjets_pfcandIndicies().at(ijet);

//             float dzmin =  100;
//             float dzmax = -100;

//             for( vector<int>::iterator ipf = pfcands.begin() ; ipf < pfcands.end() ; ++ipf ){
              
//               pfcands_tot += pfcands_p4().at(*ipf);
              
//               if( pfcands_charge().at(*ipf) == 0 ) continue;
              
//               int itrk = pfcands_trkidx().at(*ipf);

//               if( itrk > trks_trk_p4().size() ) continue;

//               float thisdz = dz_trk_vtx(itrk,vtxIdx);
//               cout << "dz " << thisdz << endl;

//               if( thisdz < dzmin )    dzmin = thisdz;
//               if( thisdz > dzmax )    dzmax = thisdz;
                
//             }
          
//             cout << endl;
//             cout << "pt " << pfjets_p4().at(ijet).pt() << " " << pfcands_tot.pt() << endl;
//             cout << "spread " << dzmax - dzmin << endl;
          }


          for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
            LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
            LorentzVector vlt  = hyp_lt_p4()[hypIdx];
            LorentzVector vll  = hyp_ll_p4()[hypIdx];
            
            if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
            if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
            if( fabs( vjet.eta() ) > 3.0 )           continue;
            if( !passesPFJetID(ijet) )               continue;
            if( vjet.pt() < 30. )                    continue;

            if( jetFromSignalPV( ijet , vtxIdx ) ){
             
              fillUnderOverFlow( hmismatch_jetpt_tot , vjet.pt() );
              fillUnderOverFlow( hmismatch_nvtx_tot  , goodVertices().size()    );
              
              if( !jetMatchedToGenJet( vjet , hypIdx) ){
                fillUnderOverFlow( hmismatch_jetpt_pass , vjet.pt() ); 
                fillUnderOverFlow( hmismatch_nvtx_pass  , goodVertices().size() ); 
              }
            }

            if( jetMatchedToGenJet( vjet , hypIdx) ){
              
              fillUnderOverFlow( hmatch_jetpt_tot , vjet.pt() );
              fillUnderOverFlow( hmatch_nvtx_tot  , goodVertices().size()     );
              
              if( jetFromSignalPV( ijet , vtxIdx ) ){
                fillUnderOverFlow( hmatch_jetpt_pass , vjet.pt() ); 
                fillUnderOverFlow( hmatch_nvtx_pass  , goodVertices().size() ); 
              }
            }
          }

          if( debug ) cout << "Filling ntuple" << endl;

          // event stuff
          run_     = cms2.evt_run();
          lumi_    = cms2.evt_lumiBlock();
          event_   = cms2.evt_event();
          leptype_ = cms2.hyp_type().at(hypIdx);
          nvtx_    = goodVertices().size();

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
          
          float pfcandmet_x = 0;
          float pfcandmet_y = 0;
          pfcandsumet_ = 0;

          int nChargedPFCandidates = 0;

          for( unsigned int i = 0 ; i < pfcands_p4().size() ; ++i ){
            pfcandmet_x  -= pfcands_p4().at(i).px();
            pfcandmet_y  -= pfcands_p4().at(i).py();
            pfcandsumet_ += pfcands_p4().at(i).pt();
            
            if( pfcands_charge().at(i) != 0 )  ++nChargedPFCandidates;

            if( pfcands_charge().at(i) != 0 )            continue;
            if( fabs( pfcands_p4().at(i).eta() ) > 2.5 ) continue;
            if( dRbetweenVectors(pfcands_p4().at(i) , hyp_lt_p4().at(hypIdx)) < 0.1)  continue;
            if( dRbetweenVectors(pfcands_p4().at(i) , hyp_ll_p4().at(hypIdx)) < 0.1)  continue;

            if( njets25_ == 0 )        fillUnderOverFlow( hneutralpt_0jet , pfcands_p4().at(i).pt()    );
            if( njets25_ == 1 )        fillUnderOverFlow( hneutralpt_1jet , pfcands_p4().at(i).pt()    );
            
            /*  
            if( njets25_ == 0 ){
              fillUnderOverFlow( hneutralpt_0jet , pfcands_p4().at(i).pt()    );

                   
              if( pfcands_p4().at(i).pt() > 25 ){
                
                cout << endl            << endl 
                     << evt_run()       << " " 
                     << evt_lumiBlock() << " " 
                     << evt_event()     << endl;
                
                cout << " pt " << pfcands_p4().at(i).pt() 
                     << " eta " << pfcands_p4().at(i).eta() 
                     << " phi " << pfcands_p4().at(i).phi() 
                     << " ID " << pfcands_particleId().at(i) << endl;

                cout << "drll " << dRbetweenVectors( pfcands_p4().at(i) , hyp_ll_p4().at(hypIdx) ) << endl;
                cout << "drlt " << dRbetweenVectors( pfcands_p4().at(i) , hyp_lt_p4().at(hypIdx) ) << endl;
        
                float drmin = 999;
                float imin  = -1;
                
                for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
                  
                  
                  float dr = dRbetweenVectors(pfcands_p4().at(i),pfjets_p4().at(ijet));

                  if( dr < drmin ){
                    drmin = dr;
                    imin = ijet;
                  }
                }
                
                cout << "drmin " << drmin << endl;
                
                if( imin > -1 ){
                  for (unsigned int ijet = imin ; ijet < imin+1 ; ijet++) {
                  
                    LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
                    LorentzVector vlt  = hyp_lt_p4()[hypIdx];
                    LorentzVector vll  = hyp_ll_p4()[hypIdx];
                    
                    cout << "jet pt " << vjet.pt() << " eta " << vjet.eta() << " phi " << vjet.phi() << endl;
                    
                    if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
                    cout << "pass ll veto" << endl;
                    if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
                    cout << "pass lt veto" << endl;
                    if( fabs( vjet.eta() ) > 3.0 )           continue;
                    cout << "pass eta cut" << endl;
                    if( !passesPFJetID(ijet) )               continue;
                    cout << "pass pfjetid" << endl;
                    if( vjet.pt() < 25. )                    continue;
                    cout << "pass ptcut" << endl;
                    
                  }
                }
              }
            }
            */

           

            
            if( abs( pfcands_particleId().at(i) ) == 11 ){
              int iel = pfcands_pfelsidx().at(i);
              float elpt = pfels_p4().at(iel).pt();
              float cpt  = pfcands_p4().at(i).pt();
              if( fabs( elpt - cpt ) > 0.001 ) cout << "ELECTRON ERROR ---------------------------" << endl;
            }
            
            if( abs( pfcands_particleId().at(i) ) == 13 ){
              int imu = pfcands_pfmusidx().at(i);
              float mupt = pfmus_p4().at(imu).pt();
              float cpt  = pfcands_p4().at(i).pt();
              if( fabs( mupt - cpt ) > 0.001 ) cout << "MUON ERROR ---------------------------" << endl;
            }
            

          }
          
          //cout << "nChargedPFCandidates " << nChargedPFCandidates << endl;
          //cout << "ntrks                " << trks_trk_p4().size() << endl << endl;
          
          pfcandmet_    = sqrt( pow( pfcandmet_x , 2 ) + pow( pfcandmet_y , 2 ) );
          pfcandmetphi_ = atan2( pfcandmet_y , pfcandmet_x );
          
          //-------------------------------------
          // Zanetti/HooberMET
          //-------------------------------------

          float dzcut  = 0.1; // dz(trk,vtx) requirement
          float etacut = 3.0; // neutral PFCandidate eta requirement
          
          zmet_     = ZanettiMET ( vtxIdx,  hypIdx, dzcut                        );
          hmet_     = HooberMET  ( vtxIdx,  hypIdx, dzcut, 1.e10 , etacut        );
          hmetpf_   = HooberMET  ( vtxIdx,  hypIdx, dzcut, 1.e10 , etacut , true );
          hmetpf0_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 0.    , etacut , true );
          hmetpf1_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 1.    , etacut , true );
          hmetpf2_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 2.    , etacut , true );
          hmetpf3_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 3.    , etacut , true );
          hmetpf4_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 4.    , etacut , true );
          hmetpf5_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 5.    , etacut , true );
          hmetpf6_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 6.    , etacut , true );
          hmetpf7_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 7.    , etacut , true );
          hmetpf8_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 8.    , etacut , true );
          hmetpf9_  = HooberMET  ( vtxIdx,  hypIdx, dzcut, 9.    , etacut , true );
          hmetpf10_ = HooberMET  ( vtxIdx,  hypIdx, dzcut, 10.   , etacut , true );
       
          //hmet_     = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut, 1.e10 , etacut , false , false );
          //hmetpf_   = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut, 1.e10 , etacut ,  true , false );
          //hmetpf4_  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,    4. , etacut ,  true , false );
          //jetzmet_  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,    4. , etacut ,  true ,  true );

          //if( fabs( hmet_2 - hmet_ )       > 0.001 ) cout << "Error hmet "    << hmet_    << " " << hmet_2 << endl;
          //if( fabs( hmetpf_2 - hmetpf_ )   > 0.001 ) cout << "Error hmetpf "  << hmetpf_  << " " << hmetpf_2 << endl;
          //if( fabs( hmetpf4_2 - hmetpf4_ ) > 0.001 ) cout << "Error hmetpf4 " << hmetpf4_ << " " << hmetpf4_2 << endl;

          if( calculateTCMET ){
            
            metStruct myMetStruct = correctedTCMET( true, ofile );
            tcmetNew_    = myMetStruct.met;
            tcsumetNew_  = myMetStruct.sumet;
            tcmetphiNew_ = myMetStruct.metphi;

          }

 	  fillUnderOverFlow( htcmet ,       tcmet_    );
	  fillUnderOverFlow( htcmetNew ,    tcmetNew_ );
 	  fillUnderOverFlow( hpfmet ,       pfmet_    );
          
          eventTree_->Fill();
          
        } // end loop over events
    } // end loop over files
  
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
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

void looper::InitBabyNtuple ()
{
  // event stuff
  run_             = -999999;
  lumi_            = -999999;
  event_           = -999999;
  nvtx_            = -999999;
  leptype_         = -999999;
  njets25_         = -999999;

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

  // met from pfcandidates
  pfcandmet_       = -999999.;
  pfcandmetphi_    = -999999.;
  pfcandsumet_     = -999999.;

  hmet_            = -999999.;
  hmetpf_          = -999999.;
  hmetpf0_         = -999999.;
  hmetpf1_         = -999999.;
  hmetpf2_         = -999999.;
  hmetpf3_         = -999999.;
  hmetpf4_         = -999999.;
  hmetpf5_         = -999999.;
  hmetpf6_         = -999999.;
  hmetpf7_         = -999999.;
  hmetpf8_         = -999999.;
  hmetpf9_         = -999999.;
  hmetpf10_        = -999999.;
  zmet_            = -999999.;
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

  hneutralpt_0jet = new TH1F("hneutralpt_0jet","",100,0,30);
  hneutralpt_1jet = new TH1F("hneutralpt_1jet","",100,0,30);

  hmismatch_jetpt_tot    =    new TH1F("hmismatch_jetpt_tot",    "", 100,0,100);
  hmismatch_jetpt_pass   =    new TH1F("hmismatch_jetpt_pass",   "", 100,0,100);
  hmatch_jetpt_tot       =    new TH1F("hmatch_jetpt_tot",       "", 100,0,100);
  hmatch_jetpt_pass      =    new TH1F("hmatch_jetpt_pass",      "", 100,0,100);

  hmismatch_jetpt_tot->Sumw2();
  hmismatch_jetpt_pass->Sumw2();
  hmatch_jetpt_tot->Sumw2();
  hmatch_jetpt_pass->Sumw2();

  hmismatch_nvtx_tot     =    new TH1F("hmismatch_nvtx_tot",     "", 15,0.5,15.5);
  hmismatch_nvtx_pass    =    new TH1F("hmismatch_nvtx_pass",    "", 15,0.5,15.5);
  hmatch_nvtx_tot        =    new TH1F("hmatch_nvtx_tot",        "", 15,0.5,15.5);
  hmatch_nvtx_pass       =    new TH1F("hmatch_nvtx_pass",       "", 15,0.5,15.5);

  hmismatch_nvtx_tot->Sumw2();
  hmismatch_nvtx_pass->Sumw2();
  hmatch_nvtx_tot->Sumw2();
  hmatch_nvtx_pass->Sumw2();
  
}

//--------------------------------------------------------------------

void looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  eventTree_ = new TTree("Events", "Events Tree");

  eventTree_->Branch("run"              , &run_              , "run/I"  );
  eventTree_->Branch("lumi"             , &lumi_             , "lumi/I" );
  eventTree_->Branch("event"            , &event_            , "event/I");
  eventTree_->Branch("nvtx"             , &nvtx_             , "nvtx/I");
  eventTree_->Branch("leptype"          , &leptype_          , "leptype/I");
  eventTree_->Branch("njets25"          , &njets25_          , "njets25/I");

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

  eventTree_->Branch("pfcandmet"        , &pfcandmet_        , "pfcandmet/F"      );
  eventTree_->Branch("pfcandmetphi"     , &pfcandmetphi_     , "pfcandmetphi/F"   );
  eventTree_->Branch("pfcandsumet"      , &pfcandsumet_      , "pfcandsumet/F"    );

  eventTree_->Branch("tcmetnew"         , &tcmetNew_         , "tcmetnew/F"      );
  eventTree_->Branch("tcmetphinew"      , &tcmetphiNew_      , "tcmetphinew/F"   );
  eventTree_->Branch("tcsumetnew"       , &tcsumetNew_       , "tcsumetnew/F"    );

  eventTree_->Branch("hmet"             , &hmet_             , "hmet/F"       );
  eventTree_->Branch("hmetpf"           , &hmetpf_           , "hmetpf/F"     );
  eventTree_->Branch("hmetpf0"          , &hmetpf0_          , "hmetpf0/F"    );
  eventTree_->Branch("hmetpf1"          , &hmetpf1_          , "hmetpf1/F"    );
  eventTree_->Branch("hmetpf2"          , &hmetpf2_          , "hmetpf2/F"    );
  eventTree_->Branch("hmetpf3"          , &hmetpf3_          , "hmetpf3/F"    );
  eventTree_->Branch("hmetpf4"          , &hmetpf4_          , "hmetpf4/F"    );
  eventTree_->Branch("hmetpf5"          , &hmetpf5_          , "hmetpf5/F"    );
  eventTree_->Branch("hmetpf6"          , &hmetpf6_          , "hmetpf6/F"    );
  eventTree_->Branch("hmetpf7"          , &hmetpf7_          , "hmetpf7/F"    );
  eventTree_->Branch("hmetpf8"          , &hmetpf8_          , "hmetpf8/F"    );
  eventTree_->Branch("hmetpf9"          , &hmetpf9_          , "hmetpf9/F"    );
  eventTree_->Branch("hmetpf10"         , &hmetpf10_         , "hmetpf10/F"   );
  eventTree_->Branch("jetzmet"          , &jetzmet_          , "jetzmet/F"    );

  eventTree_->Branch("zmet"             , &zmet_             , "zmet/F"      );

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

float looper::ZanettiMET(const unsigned int vtxIdx, const unsigned int hypIdx, float threshold)
{
    float tmet_x = 0.;
    float tmet_y = 0.;

    // start by adding hypothesis leptons to tmet
    if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
        tmet_x -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].px();
        tmet_y -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].py();
    }
    else if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
        tmet_x -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].px();
        tmet_y -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].py();
    }

    if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
        tmet_x -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].px();
        tmet_y -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].py();
    }
    else if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
        tmet_x -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].px();
        tmet_y -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].py();
    }

    // now, loop over tracks and for each track determine the closest
    // vertex;  if this vertex is the same one to which the hypothesis
    // leptons are matched AND the track is within a DeltaZ of 500 um
    // include the track in the tmet calculation
    for (unsigned int itrk = 0; itrk < cms2.trks_trk_p4().size(); itrk++) {
        float trkz = cms2.trks_vertex_p4()[itrk].z();
        float mindz = 999.;
        int vtxi = -1;
        for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
            float thisdz = cms2.vtxs_position()[ivtx].z() - trkz;
            if (fabs(thisdz) < fabs(mindz)) {
                mindz = thisdz;
                vtxi = ivtx;
            }
        }

        if (vtxi != vtxIdx)
            continue;
        if (fabs(mindz) > threshold)
            continue;
        if (dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_lt_p4()[hypIdx]) < 0.03)
            continue;
        if (dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_ll_p4()[hypIdx]) < 0.03)
            continue;

        tmet_x -= cms2.trks_trk_p4()[itrk].px();
        tmet_y -= cms2.trks_trk_p4()[itrk].py();
    } // end loop over tracks

    return sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
}

//--------------------------------------------------------------------

float looper::HooberMET(const unsigned int vtxIdx, const unsigned int hypIdx, float dz_thresh, 
                        float pt_thresh, float etacut , bool usePFCandidatePt ){
  
  float tmet_x = 0.;
  float tmet_y = 0.;
  
  //------------------------------------------------
  // start by adding hypothesis leptons to tmet
  //------------------------------------------------
  
  if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
    tmet_x -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].px();
    tmet_y -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].py();

    //cout << "Correcting for electron (pt,phi) (" << cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].pt() 
    //     << " , " << cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].phi() << endl; 
  }
  else if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
    tmet_x -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].px();
    tmet_y -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].py();
  }
  
  if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
    tmet_x -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].px();
    tmet_y -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].py();

    //cout << "Correcting for electron (pt,phi) (" << cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].pt() 
    //     << " , " << cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].phi() << endl; 

  }
  else if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
    tmet_x -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].px();
    tmet_y -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].py();
  }

  //---------------------------------------------------
  // loop over PFCandidates
  //---------------------------------------------------
  
  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {
  
    //--------------------
    //deal with neutrals
    //--------------------
    
    if( cms2.pfcands_charge().at(ipf) == 0 ){

      //--------------------
      // pt, eta cuts
      //--------------------
      
      if( pfcands_p4().at(ipf).pt()  < pt_thresh         ) continue;
      if( fabs( pfcands_p4().at(ipf).eta() ) > etacut    ) continue;
      
      if (dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_lt_p4()[hypIdx]) < 0.1)  continue;
      if (dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_ll_p4()[hypIdx]) < 0.1)  continue;

      //--------------------------------------------------
      // for photons, require dr(photon,electron) > 0.1
      //--------------------------------------------------
      /*
      float dr = 100.;
       
      if( pfcands_particleId().at(ipf) == 22 ){
        
        if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
          float dr_lt = dRbetweenVectors( pfcands_p4().at(ipf) , cms2.hyp_lt_p4()[hypIdx] );
          if( dr_lt < dr ) dr = dr_lt;
        }
        
        if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
          float dr_ll = dRbetweenVectors( pfcands_p4().at(ipf) , cms2.hyp_ll_p4()[hypIdx] );
          if( dr_ll < dr ) dr = dr_ll;
        }
        
        if( dr < 0.1 ) continue;
      }
      */
      //-----------------------------------------
      // correct MET for neutral PFCandidate
      //-----------------------------------------
      
      tmet_x -= pfcands_p4().at(ipf).px();
      tmet_y -= pfcands_p4().at(ipf).py();
      
      //cout << "Correcting for neutral PFCandidate (pt,phi) (" << cms2.pfcands_p4()[ipf].pt() 
      //     << " , " << cms2.pfcands_p4()[ipf].phi() << " id " << pfcands_particleId()[ipf] << endl; 
      
    }
    
    //-------------------------------
    // deal with charged particles
    //-------------------------------

    else{

      //------------------------------------
      // get track matched to PFCandidate
      //------------------------------------
      
      int itrk = cms2.pfcands_trkidx().at(ipf);
      
      if( itrk >= trks_trk_p4().size() || itrk < 0 ){
        //note: this should only happen for electrons which do not have a matched track
        //currently we are just ignoring these guys
        continue;
      }
      
      //----------------------------------------
      // find closest PV and dz w.r.t. that PV
      //----------------------------------------
      
      float mindz = 999.;
      int vtxi    = -1;
        
      for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
        float mydz = dz_trk_vtx(itrk,ivtx);
        //float mydz = cms2.vtxs_position()[ivtx].z() - cms2.trks_vertex_p4()[itrk].z();
        if (fabs(mydz) < fabs(mindz)) {
          mindz = mydz;
          vtxi = ivtx;
        }
      }
      
      //----------------------------------------------------------------------------
      // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
      //----------------------------------------------------------------------------
      
      if (vtxi != vtxIdx)
        continue;
      if (fabs(mindz) > dz_thresh)
        continue;
      if (dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_lt_p4()[hypIdx]) < 0.03)
        continue;
      if (dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_ll_p4()[hypIdx]) < 0.03)
        continue;
      

      //---------------------------------------
      // correct MET for charged PFCandidate
      //---------------------------------------

      if( usePFCandidatePt ){
        tmet_x -= cms2.pfcands_p4()[ipf].px();
        tmet_y -= cms2.pfcands_p4()[ipf].py();

        //cout << "Correcting for charged PFCandidate (pt,phi) (" << cms2.pfcands_p4()[ipf].pt() 
        //     << " , " << cms2.pfcands_p4()[ipf].phi() << endl; 
      }
      
      else{
        tmet_x -= cms2.trks_trk_p4()[itrk].px();
        tmet_y -= cms2.trks_trk_p4()[itrk].py();
      }
      
    } 
  }// end loop over tracks

  return sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
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

float looper::PFCandidateMET(const unsigned int vtxIdx, const unsigned int hypIdx, vector<int> goodjets, 
                             float dz_thresh, float pt_thresh, float etacut , bool usePFCandidatePt , bool correctJets ){
  
  float tmet_x = 0.;
  float tmet_y = 0.;
  
  //------------------------------------------------
  // start by adding hypothesis leptons to tmet
  //------------------------------------------------
  
  if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
    tmet_x -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].px();
    tmet_y -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].py();
  }
  else if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
    tmet_x -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].px();
    tmet_y -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].py();
  }
  
  if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
    tmet_x -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].px();
    tmet_y -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].py();
  }
  else if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
    tmet_x -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].px();
    tmet_y -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].py();
  }
  
  //------------------------------------------------
  // next, add selected jets to met
  //------------------------------------------------
  
  if( correctJets ){

    for( vector<int>::iterator igoodjet = goodjets.begin() ; igoodjet < goodjets.end() ; ++igoodjet ){
      
      LorentzVector vjet = pfjets_p4().at(*igoodjet) * pfjets_cor().at(*igoodjet);
      
      tmet_x -= vjet.px();
      tmet_y -= vjet.py();
      
    }
    
  }
  
  //---------------------------------------------------
  // loop over PFCandidates
  //---------------------------------------------------

  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {

    //----------------------------------------------------------------------
    // exclude particles belonging to jets that were already corrected for
    //----------------------------------------------------------------------

    if( correctJets ){
      
      bool skipPFCandidate = false;
      
      for( vector<int>::iterator jet_it = goodjets.begin() ; jet_it != goodjets.end() ; jet_it++ ){
        
        vector<int> jetConstituents = pfjets_pfcandIndicies().at(*jet_it);
        
        for( vector<int>::iterator pf_it = jetConstituents.begin() ; pf_it != jetConstituents.end() ; pf_it++ ){
          if( ipf == *pf_it ) skipPFCandidate = true;
        }
        
      }
    
      if( skipPFCandidate ) continue;
      
    }

    //--------------------
    //deal with neutrals
    //--------------------
  
    if( cms2.pfcands_charge().at(ipf) == 0 ){
      
      //--------------------
      // pt, eta cuts
      //--------------------
      
      if( pfcands_p4().at(ipf).pt()          < pt_thresh ) continue;
      if( fabs( pfcands_p4().at(ipf).eta() ) > etacut    ) continue;
      
      if (dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_lt_p4()[hypIdx]) < 0.1)  continue;
      if (dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_ll_p4()[hypIdx]) < 0.1)  continue;
      
      //-----------------------------------------
      // correct MET for neutral PFCandidate
      //-----------------------------------------
      
      tmet_x -= pfcands_p4().at(ipf).px();
      tmet_y -= pfcands_p4().at(ipf).py();
      
    }
    
    //-------------------------------
    // deal with charged particles
    //-------------------------------
    
    else{
      
       //------------------------------------
       // get track matched to PFCandidate
       //------------------------------------
    
       int itrk = cms2.pfcands_trkidx().at(ipf);
    
       if( itrk >= trks_trk_p4().size() || itrk < 0 ){
         //note: this should only happen for electrons which do not have a matched track
         //currently we are just ignoring these guys
         continue;
       }
    
       //----------------------------------------
       // find closest PV and dz w.r.t. that PV
       //----------------------------------------
    
       float mindz = 999.;
       int vtxi    = -1;
      
       for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
       
         float mydz = dz_trk_vtx(itrk,ivtx);
         //float mydz = cms2.vtxs_position()[ivtx].z() - cms2.trks_vertex_p4()[itrk].z();
         
         if (fabs(mydz) < fabs(mindz)) {
           mindz = mydz;
           vtxi = ivtx;
         }
         
       }
    
       //----------------------------------------------------------------------------
       // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
       //----------------------------------------------------------------------------
    
       if ( vtxi != vtxIdx )                                                               continue;
       if ( fabs(mindz) > dz_thresh )                                                      continue;
       if ( dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_lt_p4()[hypIdx]) < 0.03 )  continue;
       if ( dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_ll_p4()[hypIdx]) < 0.03 )  continue;
    

       //---------------------------------------
       // correct MET for charged PFCandidate
       //---------------------------------------

       if( usePFCandidatePt ){
         tmet_x -= cms2.pfcands_p4()[ipf].px();
         tmet_y -= cms2.pfcands_p4()[ipf].py();
       }
    
       else{
         tmet_x -= cms2.trks_trk_p4()[itrk].px();
         tmet_y -= cms2.trks_trk_p4()[itrk].py();
       }
    
     } 
   }// end loop over tracks
  

   return sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
 }



 bool looper::jetFromSignalPV( int ijet , int vtxIdx ){

   //cout << endl << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;

   //------------------------------------------
   // loop over PFCandidates matched to pfjet
   //------------------------------------------
   vector<int> pfcands = pfjets_pfcandIndicies().at(ijet);

   //cout << pfcands.size() << " matched PFCandidates" << endl;

   LorentzVector v_pfcands(0,0,0,0);

   vector<int> matchedTracks;
   matchedTracks.clear();

   float sumPtTot = 0.;

   for( unsigned int i = 0 ; i < pfcands.size() ; ++i ){
     
     int ipf = pfcands.at(i);

     v_pfcands += pfcands_p4().at(ipf);

     //---------------------------------------------
     // find track index of charged PFCandidates
     //---------------------------------------------
     
     if( pfcands_charge().at(ipf) == 0 ) continue;

     int itrk = cms2.pfcands_trkidx().at(ipf);
     
     //note: this should only happen for electrons which do not have a matched track, currently ignoring these guys
     if( itrk >= trks_trk_p4().size() || itrk < 0 )  continue;
     
     //----------------------------------
     // store indices of matched tracks
     //----------------------------------
     
     //cout << "Storing matched track " << itrk << " pt " << trks_trk_p4().at(itrk).pt() << endl;
     matchedTracks.push_back( itrk );

     sumPtTot += trks_trk_p4().at(itrk).pt();

   }

   //cout << "Found " << matchedTracks.size() << " matched tracks, sumpt " << sumPtTot << endl;

   if( fabs( pfjets_p4().at(ijet).pt() - v_pfcands.pt() ) > 0.1 ){
     cout << "Warning: pfjet pt " << pfjets_p4().at(ijet).pt() 
          << " doesn't match sum of PFCandidates pt " << v_pfcands.pt() << endl;
   }

   vector<int> myGoodVertices = goodVertices();
   const unsigned int nGoodVertices = myGoodVertices.size();

   if( nGoodVertices == 0 ){
     cout << "Didn't find any good vertices!" << endl;
     return false;
   }

   float beta[nGoodVertices];
   for( unsigned int ivtx = 0 ; ivtx < nGoodVertices ; ++ivtx ) beta[ivtx] = 0.0;

   for (vector<int>::iterator itrk = matchedTracks.begin(); itrk != matchedTracks.end(); ++itrk) {
     
     float mindz = 1000;
     int   vtxi  = -1;
     
     for (vector<int>::iterator ivtx = myGoodVertices.begin(); ivtx != myGoodVertices.end() ; ++ivtx ){
       
       float thisdz = dz_trk_vtx(*itrk,*ivtx);
       
       if (fabs(thisdz) < fabs(mindz)) {
         mindz = thisdz;
         vtxi  = *ivtx;
       }
     }
     
     beta[vtxi] += trks_trk_p4().at(*itrk).pt() / sumPtTot;
     
   }
   
   float betamax = -1.;
   int   imax    = -1;
   for( unsigned int ivtx = 0 ; ivtx < nGoodVertices ; ++ivtx ){
     //cout << "beta[" << ivtx << "] " << beta[ivtx] << endl;

     if( beta[ivtx] > betamax ){
       betamax = beta[ivtx];
       imax    = ivtx;
     }
   }

   //cout << "imax " << imax << " betamax " << betamax << endl; 

   if( imax == vtxIdx ){
     //cout << "It's good!" << endl;
     return true;
   }

   return false;
 }

//--------------------------------------------------------------------

float looper::dz_trk_vtx(const unsigned int trkidx, const unsigned int vtxidx){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

//--------------------------------------------------------------------
