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

const float lumi         = 0.023;
const char* path         = "output/PVT/promptreco/dcsonly";
const bool  doTenPercent = true;

char* iter          = "default";
bool makebaby       = true;
bool debug          = false;
bool calculateTCMET = false;
bool doReweight     = true;

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
    if( fabs( vgenjet.eta() ) > 5.0 )            continue;
    if( vgenjet.pt() < 30. )                     continue;
    if( dRbetweenVectors(vgenjet, vjet) > 0.4 )  continue;

    return true;
  }

  return false;

}

//--------------------------------------------------------------------

using namespace tas;
void looper::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  if( doTenPercent ){
    cout << endl;
    cout << "-----------------------------------------------" << endl;
    cout << "| PROCESSING TEN PERCENT OF SAMPLE!!!!!!!!!!!!|" << endl;
    cout << "-----------------------------------------------" << endl;
    cout << endl;
  }


  bookHistos();

  set_goodrun_file("json_DCSONLY.txt_160404-161312.goodruns");
  //set_goodrun_file("Cert_160404-161216_7TeV_PromptReco_Collisions11_JSON_goodruns.txt");

  ofile.open( Form( "%s/%s_%s_events.txt" , path , prefix , iter) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  MakeBabyNtuple( Form( "%s/%s_%s_baby.root" , path , prefix , iter) );

  if( debug ) cout << "Begin looping over files" << endl;

  TH1F* h_reweight    = new TH1F();
  TH1F* h_da_reweight = new TH1F();

  if( doReweight ){
    TFile* f_reweight    = TFile::Open("vtx_reweight.root");
    TFile* f_da_reweight = TFile::Open("vtx_DA_reweight.root");

    h_reweight    = (TH1F*) f_reweight->Get("hratio");
    h_da_reweight = (TH1F*) f_da_reweight->Get("hratio");

    cout << "Doing reweighting (standard)" << endl;
    for( unsigned int ibin = 1 ; ibin <= h_reweight->GetNbinsX() ; ibin++ ){
      cout << ibin << " " << h_reweight->GetBinContent(ibin) << endl;
    }

    cout << "Doing reweighting (DA)" << endl;
    for( unsigned int ibin = 1 ; ibin <= h_da_reweight->GetNbinsX() ; ibin++ ){
      cout << ibin << " " << h_da_reweight->GetBinContent(ibin) << endl;
    }

  }

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
      unsigned int nEvents = tree->GetEntries()/100;

      for (unsigned int event = 0; event < nEvents; ++event)
        {
          if( debug ) cout << "Event " << event << endl;

          ++nEventsTotal;

	  // //SKIP 9/10 EVENTS!!!
	  // if( doTenPercent ){
	  //   if( nEventsTotal%100 != 0 ) continue;
	  // }

          cms2.GetEntry(event);


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

	  if ( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))  continue;
          if( !cleaning_standardOctober2010() )                           continue;
          
          //------------------------------------------
          // hyp selection
          //------------------------------------------

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
            
	    //if leading lepton is electron, require pt > 27 GeV
	    //int id = -1;
	    //if( hyp_ll_p4()[i].pt() > hyp_lt_p4()[i].pt() ) id = hyp_ll_id()[i];
	    //else                                            id = hyp_lt_id()[i];
	    //if( abs(id) == 11 && TMath::Max( hyp_ll_p4()[i].pt() , hyp_lt_p4()[i].pt() ) < 27. )   continue;

            // OSZ lepton selection
            //if (abs(hyp_ll_id()[i]) == 13  && !( muonId(hyp_ll_index()[i] , OSZ_v1 ) ) )   continue;
            //if (abs(hyp_lt_id()[i]) == 13  && !( muonId(hyp_lt_index()[i] , OSZ_v1 ) ) )   continue;
            //if (abs(hyp_ll_id()[i]) == 11  && !( pass_electronSelection( hyp_ll_index()[i] , electronSelection_el_OSV1 ))) continue;
            //if (abs(hyp_lt_id()[i]) == 11  && !( pass_electronSelection( hyp_lt_index()[i] , electronSelection_el_OSV1 ))) continue;

            //WW lepton selection
            if (abs(hyp_ll_id()[i]) == 13  && !( muonId(hyp_ll_index()[i] , NominalWWV1 ) ) )   continue;
            if (abs(hyp_lt_id()[i]) == 13  && !( muonId(hyp_lt_index()[i] , NominalWWV1 ) ) )   continue;
            if (abs(hyp_ll_id()[i]) == 11  && !( pass_electronSelection( hyp_ll_index()[i] , electronSelection_wwV1 ))) continue;
            if (abs(hyp_lt_id()[i]) == 11  && !( pass_electronSelection( hyp_lt_index()[i] , electronSelection_wwV1 ))) continue;
            
            v_goodHyps.push_back( i );
            
	    }
	  


          //skip events with no good hyps
          if( v_goodHyps.size() == 0 ) continue;

          if( debug ) cout << "Pass event selection" << endl;

          //returns the index of the best hypothesis in the vector of hypotheses
          unsigned int hypIdx = selectHypByHighestSumPt(v_goodHyps);
          int          vtxIdx = hypsFromSameVtx_int( hypIdx );

          InitBabyNtuple();
          
          //---------------------------------------
          // jet counting
          //---------------------------------------

          njets25_ = 0;
          njets30_ = 0;

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
            if( fabs( vjet.eta() ) > 5.0 )           continue;
            if( !passesPFJetID(ijet) )               continue;

            //if( !jetFromSignalPV( ijet , vtxIdx ) )  continue;

            if( vjet.pt() < 25. )                    continue;
           
            njets25_++;

            if( vjet.pt() < 30. )                    continue;

            njets30_++;
            goodjets.push_back(ijet);
            
            if( vjet.pt() > maxjetpt ){
              imaxjet = ijet;
              maxjetpt = vjet.pt();
            }
          }
          
          if( imaxjet > -1 ){
            LorentzVector vjet = pfjets_corL1L2L3().at(imaxjet) * pfjets_p4().at(imaxjet);

            jetpt_   = vjet.pt();
            jeteta_  = vjet.eta();
            jetphi_  = vjet.phi();
            jetpv_   = jetFromSignalPV( imaxjet , vtxIdx , 2 ) ? 1 : 0;
            jetbeta_ = beta_jet_vtx( imaxjet , vtxIdx , 2 );
            jet_     = &(pfjets_corL1L2L3().at(imaxjet) * pfjets_p4().at(imaxjet));
	    jetgen_  = jetMatchedToGenJet( vjet , hypIdx) ? 1 : 0;
          }
        

          if( !isData ){
            for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
              
              LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
              LorentzVector vlt  = hyp_lt_p4()[hypIdx];
              LorentzVector vll  = hyp_ll_p4()[hypIdx];
              
              if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
              if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
              if( fabs( vjet.eta() ) > 5.0 )           continue;
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
          }
          

          if( debug ) cout << "Filling ntuple" << endl;

          // event stuff
          run_     = cms2.evt_run();
          lumi_    = cms2.evt_lumiBlock();
          event_   = cms2.evt_event();
          leptype_ = cms2.hyp_type().at(hypIdx);
          nvtx_    = goodVertices().size();
	  ndavtx_  = goodDAVertices().size();
          dilmass_ = hyp_p4().at(hypIdx).mass();
          vtxIdx_  = vtxIdx;
          dilep_   = &hyp_p4().at(hypIdx);


          weight_      = 1;
          davtxweight_ = 1;
	  vtxweight_   = 1;

          if( !isData ){

            if( doReweight ){
              vtxweight_   = h_reweight->GetBinContent(    nvtx_   + 1 );
              davtxweight_ = h_da_reweight->GetBinContent( ndavtx_ + 1 );
            }

	    //if( TString(prefix).Contains("dymm_spring11") ){
	    //  weight_ = ( 0.84 ) * lumi;
	    // }
	    //else{
	    weight_ = evt_scale1fb() * lumi;
	    //}
	  }

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
	  pair<float,float> p_tcmet = getMet( "tcMET"    , hypIdx);
	  tcmet_    = p_tcmet.first;
	  tcmetphi_ = p_tcmet.second;

          //tcmet_    = cms2.evt_tcmet();
          //tcmetphi_ = cms2.evt_tcmetPhi();
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

            if( njets30_ == 0 )        fillUnderOverFlow( hneutralpt_0jet , pfcands_p4().at(i).pt()    );
            if( njets30_ == 1 )        fillUnderOverFlow( hneutralpt_1jet , pfcands_p4().at(i).pt()    );
            
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
          

          pair<float, float> p_test = pfmetByHand( -1 , 1e20 );
          float pfmet_test = p_test.first;
          if( fabs( pfmet_test - evt_pfmet() ) > 1 ){
            cout << "pfmet " << evt_pfmet() << " by hand " << pfmet_test << endl;
          }

          //-------------------------------------
          // custom MET types
          //-------------------------------------
         
          float dzcut  = 0.1; // dz(trk,vtx) requirement
          float etacut = 3.0; // neutral PFCandidate eta requirement

          //met built from tracks
          pair<float, float> p_pfmetByHand = pfmetByHand( -1 , 3.0 );
          pfmet3_      = p_pfmetByHand.first;
          pfmetphi3_   = p_pfmetByHand.second;
          pfmet3proj_  = projectedMET( pfmet3_ , pfmetphi3_ , hypIdx );

          //met built from tracks
          pair<float, float> p_zmet = ZanettiMET( vtxIdx,  hypIdx, dzcut );
          zmet_      = p_zmet.first;
          zmetphi_   = p_zmet.second;
          zmetproj_  = projectedMET( zmet_ , zmetphi_ , hypIdx );

          //met built from tracks matched to PFCandidates
          pair<float, float> p_hmet    = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut, 1.e10 , etacut , false , false );
          hmet_      = p_hmet.first;
          hmetphi_   = p_hmet.second;
          hmetproj_  = projectedMET( hmet_ , hmetphi_ , hypIdx );

          //met built from leptons only
          pair<float, float> p_hmetpfnotrks    = PFCandidateMET( vtxIdx,  hypIdx, goodjets, -1. , 1.e10 , etacut , true , false );
          hmetpfnotrks_        = p_hmetpfnotrks.first;
          hmetphipfnotrks_     = p_hmetpfnotrks.second;
          hmetpfnotrksproj_    = projectedMET( hmetpfnotrks_ , hmetphipfnotrks_ , hypIdx );

          //met built from charged PFCandidates
          pair<float, float> p_hmetpf  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut, 1.e10 , etacut ,  true , false );
          hmetpf_      = p_hmetpf.first;
          hmetphipf_   = p_hmetpf.second;
          hmetpfproj_  = projectedMET( hmetpf_ , hmetphipf_ , hypIdx );

          //met built from charged PFCandidates and neutral PFCandidates pt > 4 GeV
          pair<float, float> p_hmetpf4  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,    4. , etacut ,  true , false );
          hmetpf4_      = p_hmetpf4.first;
          hmetphipf4_   = p_hmetpf4.second;
          hmetpf4proj_  = projectedMET( hmetpf4_ , hmetphipf4_ , hypIdx );

          //met built from charged PFCandidates and neutral PFCandidates pt > 8 GeV
          pair<float, float> p_hmetpf8  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,    8. , etacut ,  true , false );
          hmetpf8_      = p_hmetpf8.first;
          hmetphipf8_   = p_hmetpf8.second;
          hmetpf8proj_  = projectedMET( hmetpf8_ , hmetphipf8_ , hypIdx );

          //met built from jets and leptons
          pair<float, float> p_jetzmetnotrks  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, -1. ,  1.e10 , etacut ,  true ,  true );
          jetzmetnotrks_         = p_jetzmetnotrks.first;
          jetzmetphinotrks_      = p_jetzmetnotrks.second;
          jetzmetnotrksproj_     = projectedMET( jetzmetnotrks_ , jetzmetphinotrks_ , hypIdx );
          
          //met built from jets and charged PFCandidates
          pair<float, float> p_jetzmet  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,  1.e10 , etacut ,  true ,  true );
          jetzmet_      = p_jetzmet.first;
          jetzmetphi_   = p_jetzmet.second;
          jetzmetproj_  = projectedMET( jetzmet_ , jetzmetphi_ , hypIdx );          

          //met built from jets and charged PFCandidates and neutral PFCandidates pt > 4 GeV
          pair<float, float> p_jetzmet4  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,     4. , etacut ,  true ,  true );
          jetzmet4_      = p_jetzmet4.first;
          jetzmetphi4_   = p_jetzmet4.second;
          jetzmet4proj_  = projectedMET( jetzmet4_ , jetzmetphi4_ , hypIdx );          
          
          //met built from jets and charged PFCandidates
          pair<float, float> p_jetzmet8  = PFCandidateMET( vtxIdx,  hypIdx, goodjets, dzcut,     8. , etacut ,  true ,  true );
          jetzmet8_      = p_jetzmet8.first;
          jetzmetphi8_   = p_jetzmet8.second;
          jetzmet8proj_  = projectedMET( jetzmet8_ , jetzmetphi8_ , hypIdx );          

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

  cout << "Yields passing selection" << endl;
  cout << "ee    " << nee  << endl;
  cout << "mm    " << nmm  << endl;
  cout << "em    " << nem  << endl;
  cout << "to    " << ntot << endl;

  
  CloseBabyNtuple();

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form( "%s/%s_%s_histos.root" , path , prefix , iter ) );
  deleteHistos();

  already_seen.clear();
} // end ScanChain

//--------------------------------------------------------------------

void looper::printEvent(  ostream& ostr ){
  ostr << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl; 
}

//--------------------------------------------------------------------

void looper::InitBabyNtuple ()
{
  // event stuff
  weight_          = -999999.;
  vtxweight_       = -999999.;
  davtxweight_     = -999999.;
  run_             = -999999;
  lumi_            = -999999;
  event_           = -999999;
  nvtx_            = -999999;
  ndavtx_          = -999999;
  vtxIdx_          = -999999;
  dilmass_         = -999999.;
  leptype_         = -999999;
  njets25_         = -999999;
  njets30_         = -999999;

  //leading jet stuff
  jetpt_           = -999999.;
  jeteta_          = -999999.;
  jetphi_          = -999999.;
  jetpv_           = -999999;
  jetbeta_         = -999999.;
  jetgen_          = -999999.;
  
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
  hmetpfnotrks_    = -999999.;
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
  jetzmetnotrks_   = -999999.;
  jetzmet_         = -999999.;
  jetzmet4_        = -999999.;
  jetzmet8_        = -999999.;
  pfmet3_          = -999999.;


  hmetphi_            = -999999.;
  hmetphipf_          = -999999.;
  hmetphipfnotrks_    = -999999.;
  hmetphipf0_         = -999999.;
  hmetphipf1_         = -999999.;
  hmetphipf2_         = -999999.;
  hmetphipf3_         = -999999.;
  hmetphipf4_         = -999999.;
  hmetphipf5_         = -999999.;
  hmetphipf6_         = -999999.;
  hmetphipf7_         = -999999.;
  hmetphipf8_         = -999999.;
  hmetphipf9_         = -999999.;
  hmetphipf10_        = -999999.;
  zmetphi_            = -999999.;
  jetzmetphi_         = -999999.;
  jetzmetphinotrks_   = -999999.;
  jetzmetphi4_        = -999999.;
  jetzmetphi8_        = -999999.;
  pfmetphi3_          = -999999.;

  hmetproj_            = -999999.;
  hmetpfproj_          = -999999.;
  hmetpfnotrksproj_    = -999999.;
  hmetpf0proj_         = -999999.;
  hmetpf1proj_         = -999999.;
  hmetpf2proj_         = -999999.;
  hmetpf3proj_         = -999999.;
  hmetpf4proj_         = -999999.;
  hmetpf5proj_         = -999999.;
  hmetpf6proj_         = -999999.;
  hmetpf7proj_         = -999999.;
  hmetpf8proj_         = -999999.;
  hmetpf9proj_         = -999999.;
  hmetpf10proj_        = -999999.;
  zmetproj_            = -999999.;
  jetzmetproj_         = -999999.;
  jetzmetnotrksproj_   = -999999.;
  jetzmet4proj_        = -999999.;
  jetzmet8proj_        = -999999.;
  pfmet3proj_          = -999999.;

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

  eventTree_->Branch("dilep", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &dilep_ );
  eventTree_->Branch("lep1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep1_  );
  eventTree_->Branch("lep2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &lep2_  );
  eventTree_->Branch("jet"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", &jet_   );

  eventTree_->Branch("weight"           , &weight_           , "weight/F"  );
  eventTree_->Branch("vtxweight"        , &vtxweight_        , "vtxweight/F"  );
  eventTree_->Branch("davtxweight"      , &davtxweight_      , "davtxweight/F"  );
  eventTree_->Branch("run"              , &run_              , "run/I"  );
  eventTree_->Branch("lumi"             , &lumi_             , "lumi/I" );
  eventTree_->Branch("event"            , &event_            , "event/I");
  eventTree_->Branch("nvtx"             , &nvtx_             , "nvtx/I");
  eventTree_->Branch("ndavtx"           , &ndavtx_           , "ndavtx/I");
  eventTree_->Branch("vtxIdx"           , &vtxIdx_           , "vtxIdx/I");
  eventTree_->Branch("dilmass"          , &dilmass_          , "dilmass/F");
  eventTree_->Branch("leptype"          , &leptype_          , "leptype/I");
  eventTree_->Branch("njets25"          , &njets25_          , "njets25/I");
  eventTree_->Branch("njets30"          , &njets30_          , "njets30/I");

  eventTree_->Branch("jetpt"            , &jetpt_            , "jetpt/F"   );
  eventTree_->Branch("jetgen"           , &jetgen_           , "jetgen/F"  );
  eventTree_->Branch("jeteta"           , &jeteta_           , "jeteta/F"  );
  eventTree_->Branch("jetphi"           , &jetphi_           , "jetphi/F"  );
  eventTree_->Branch("jetpv"            , &jetpv_            , "jetpv/I"   );
  eventTree_->Branch("jetbeta"          , &jetbeta_          , "jetbeta/F" );


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

  eventTree_->Branch("hmet"             , &hmet_             , "hmet/F"           );
  eventTree_->Branch("hmetpfnotrks"     , &hmetpfnotrks_     , "hmetpfnotrks/F"   );
  eventTree_->Branch("hmetpf"           , &hmetpf_           , "hmetpf/F"         );
  eventTree_->Branch("hmetpf0"          , &hmetpf0_          , "hmetpf0/F"        );
  eventTree_->Branch("hmetpf1"          , &hmetpf1_          , "hmetpf1/F"        );
  eventTree_->Branch("hmetpf2"          , &hmetpf2_          , "hmetpf2/F"        );
  eventTree_->Branch("hmetpf3"          , &hmetpf3_          , "hmetpf3/F"        );
  eventTree_->Branch("hmetpf4"          , &hmetpf4_          , "hmetpf4/F"        );
  eventTree_->Branch("hmetpf5"          , &hmetpf5_          , "hmetpf5/F"        );
  eventTree_->Branch("hmetpf6"          , &hmetpf6_          , "hmetpf6/F"        );
  eventTree_->Branch("hmetpf7"          , &hmetpf7_          , "hmetpf7/F"        );
  eventTree_->Branch("hmetpf8"          , &hmetpf8_          , "hmetpf8/F"        );
  eventTree_->Branch("hmetpf9"          , &hmetpf9_          , "hmetpf9/F"        );
  eventTree_->Branch("hmetpf10"         , &hmetpf10_         , "hmetpf10/F"       );
  eventTree_->Branch("jetzmetnotrks"    , &jetzmetnotrks_    , "jetzmetnotrks/F"  );
  eventTree_->Branch("jetzmet"          , &jetzmet_          , "jetzmet/F"        );
  eventTree_->Branch("jetzmet4"         , &jetzmet4_         , "jetzmet4/F"       );
  eventTree_->Branch("jetzmet8"         , &jetzmet8_         , "jetzmet8/F"       );
  eventTree_->Branch("zmet"             , &zmet_             , "zmet/F"           );
  eventTree_->Branch("pfmet3"           , &pfmet3_           , "pfmet3/F"         );

  eventTree_->Branch("hmetphi"          , &hmetphi_          , "hmetphi/F"            );
  eventTree_->Branch("hmetphipfnotrks"  , &hmetphipfnotrks_  , "hmetphipfnotrks/F"    );
  eventTree_->Branch("hmetphipf"        , &hmetphipf_        , "hmetphipf/F"          );
  eventTree_->Branch("hmetphipf0"       , &hmetphipf0_       , "hmetphipf0/F"         );
  eventTree_->Branch("hmetphipf1"       , &hmetphipf1_       , "hmetphipf1/F"         );
  eventTree_->Branch("hmetphipf2"       , &hmetphipf2_       , "hmetphipf2/F"         );
  eventTree_->Branch("hmetphipf3"       , &hmetphipf3_       , "hmetphipf3/F"         );
  eventTree_->Branch("hmetphipf4"       , &hmetphipf4_       , "hmetphipf4/F"         );
  eventTree_->Branch("hmetphipf5"       , &hmetphipf5_       , "hmetphipf5/F"         );
  eventTree_->Branch("hmetphipf6"       , &hmetphipf6_       , "hmetphipf6/F"         );
  eventTree_->Branch("hmetphipf7"       , &hmetphipf7_       , "hmetphipf7/F"         );
  eventTree_->Branch("hmetphipf8"       , &hmetphipf8_       , "hmetphipf8/F"         );
  eventTree_->Branch("hmetphipf9"       , &hmetphipf9_       , "hmetphipf9/F"         );
  eventTree_->Branch("hmetphipf10"      , &hmetphipf10_      , "hmetphipf10/F"        );
  eventTree_->Branch("jetzmetphinotrks" , &jetzmetphinotrks_ , "jetzmetphinotrks/F"   );
  eventTree_->Branch("jetzmetphi"       , &jetzmetphi_       , "jetzmetphi/F"         );
  eventTree_->Branch("jetzmetphi4"      , &jetzmetphi4_      , "jetzmetphi4/F"        );
  eventTree_->Branch("jetzmetphi8"      , &jetzmetphi8_      , "jetzmetphi8/F"        );
  eventTree_->Branch("zmetphi"          , &zmetphi_          , "zmetphi/F"            );
  eventTree_->Branch("pfmetphi3"        , &pfmetphi3_        , "pfmetphi3/F"         );

  eventTree_->Branch("hmetproj"          , &hmetproj_          , "hmetproj/F"            );
  eventTree_->Branch("hmetpfnotrksproj"  , &hmetpfnotrksproj_  , "hmetpfnotrksproj/F"    );
  eventTree_->Branch("hmetpfproj"        , &hmetpfproj_        , "hmetpfproj/F"          );
  eventTree_->Branch("hmetpf0proj"       , &hmetpf0proj_       , "hmetpf0proj/F"         );
  eventTree_->Branch("hmetpf1proj"       , &hmetpf1proj_       , "hmetpf1proj/F"         );
  eventTree_->Branch("hmetpf2proj"       , &hmetpf2proj_       , "hmetpf2proj/F"         );
  eventTree_->Branch("hmetpf3proj"       , &hmetpf3proj_       , "hmetpf3proj/F"         );
  eventTree_->Branch("hmetpf4proj"       , &hmetpf4proj_       , "hmetpf4proj/F"         );
  eventTree_->Branch("hmetpf5proj"       , &hmetpf5proj_       , "hmetpf5proj/F"         );
  eventTree_->Branch("hmetpf6proj"       , &hmetpf6proj_       , "hmetpf6proj/F"         );
  eventTree_->Branch("hmetpf7proj"       , &hmetpf7proj_       , "hmetpf7proj/F"         );
  eventTree_->Branch("hmetpf8proj"       , &hmetpf8proj_       , "hmetpf8proj/F"         );
  eventTree_->Branch("hmetpf9proj"       , &hmetpf9proj_       , "hmetpf9proj/F"         );
  eventTree_->Branch("hmetpf10proj"      , &hmetpf10proj_      , "hmetpf10proj/F"        );
  eventTree_->Branch("jetzmetnotrksproj" , &jetzmetnotrksproj_ , "jetzmetnotrksproj/F"   );
  eventTree_->Branch("jetzmetproj"       , &jetzmetproj_       , "jetzmetproj/F"         );
  eventTree_->Branch("jetzmet4proj"      , &jetzmet4proj_      , "jetzmet4proj/F"        );
  eventTree_->Branch("jetzmet8proj"      , &jetzmet8proj_      , "jetzmet8proj/F"        );
  eventTree_->Branch("zmetproj"          , &zmetproj_          , "zmetproj/F"            );
  eventTree_->Branch("pfmet3proj"        , &pfmet3proj_        , "pfmet3proj/F"         );
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

std::pair<float,float> looper::pfmetByHand( float ptcut , float etacut ){

  float met_x = 0;
  float met_y = 0;

  for( unsigned int i = 0 ; i < pfcands_p4().size() ; ++i ){

    if( pfcands_p4().at(i).pt()          < ptcut  ) continue;
    if( fabs( pfcands_p4().at(i).eta() ) > etacut ) continue;

    met_x  -= pfcands_p4().at(i).px();
    met_y  -= pfcands_p4().at(i).py();
    
  }

  float met = sqrt( met_x * met_x + met_y * met_y);
  float metphi = atan2( met_y , met_x );
  
  return make_pair( met , metphi );

}

//--------------------------------------------------------------------

std::pair<float,float> looper::ZanettiMET(const unsigned int vtxIdx, const unsigned int hypIdx, float threshold)
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

  float met = sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
  float metphi = atan2( tmet_y , tmet_x );
  
  return make_pair( met , metphi );

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

vector<int> looper::goodDAVertices(){

  vector<int> myVertices;
  myVertices.clear();
  
  for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
    if( !isGoodDAVertex(v) ) continue;
    myVertices.push_back(v);
  }
  
  return myVertices;
}

//--------------------------------------------------------------------

std::pair<float,float> looper::PFCandidateMET(const unsigned int vtxIdx, const unsigned int hypIdx, vector<int> goodjets, 
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
      
      LorentzVector vjet = pfjets_p4().at(*igoodjet) * pfjets_corL1L2L3().at(*igoodjet);
      
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
  
  float met = sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
  float metphi = atan2( tmet_y , tmet_x );
  
  return make_pair( met , metphi );
}

float looper::beta_jet_vtx( int ijet , int vtxIdx , int beta_exponent ){

  //------------------------------------------
  // loop over PFCandidates matched to pfjet
  //------------------------------------------
  
  vector<int> pfcands = pfjets_pfcandIndicies().at(ijet);
  
  LorentzVector v_pfcands(0,0,0,0);

  vector<int> matchedTracks;
  matchedTracks.clear();

  float sumPtTot = 0.;

  for( vector<int>::iterator ipf = pfcands.begin() ; ipf < pfcands.end() ; ++ipf ){
     
    v_pfcands += pfcands_p4().at(*ipf);

    //---------------------------------------------
    // find track index of charged PFCandidates
    //---------------------------------------------
     
    if( pfcands_charge().at(*ipf) == 0 ) continue;

    int itrk = cms2.pfcands_trkidx().at(*ipf);
     
    //note: this should only happen for electrons which do not have a matched track, currently ignoring these guys
    if( itrk >= trks_trk_p4().size() || itrk < 0 )  continue;
     
    //----------------------------------
    // store indices of matched tracks
    //----------------------------------
     
    matchedTracks.push_back( itrk );

    sumPtTot += pow( trks_trk_p4().at(itrk).pt() , beta_exponent );

  }

  //-------------
  //sanity check
  //-------------
  
  if( fabs( pfjets_p4().at(ijet).pt() - v_pfcands.pt() ) > 0.1 ){
    cout << "Warning: pfjet pt " << pfjets_p4().at(ijet).pt() 
         << " doesn't match sum of PFCandidates pt " << v_pfcands.pt() << endl;
  }

  //----------------------------------
  // find good vertices
  //----------------------------------

  vector<int> myGoodVertices = goodVertices();
  const unsigned int nGoodVertices = myGoodVertices.size();

  if( nGoodVertices == 0 ){
    cout << "Didn't find any good vertices!" << endl;
    return -1.;
  }

  float beta = 0.;

  //------------------------------------------------------------------------------------------------
  // loop over tracks. if track_i is closest to signal PV, add pow(track_i pt,beta_exponent) to beta
  //------------------------------------------------------------------------------------------------

  for (vector<int>::iterator itrk = matchedTracks.begin(); itrk != matchedTracks.end(); ++itrk) {
     
    float mindz = 1000;
    int   vtxi  = -1;
     
    for (vector<int>::iterator ivtx = myGoodVertices.begin(); ivtx != myGoodVertices.end() ; ++ivtx ){
       
      float thisdz = dz_trk_vtx(*itrk,*ivtx);
       
      if ( fabs( thisdz ) < fabs( mindz ) ) {
        mindz = thisdz;
        vtxi  = *ivtx;
      }
    }
     
    if( vtxi == vtxIdx ){
      beta += pow( trks_trk_p4().at(*itrk).pt() , beta_exponent ) / sumPtTot;
    }
  }
 
  return beta;

}

//--------------------------------------------------------------------

bool looper::jetFromSignalPV( int ijet , int vtxIdx , int beta_exponent ){

  vector<int> myGoodVertices       = goodVertices();
  const unsigned int nGoodVertices = myGoodVertices.size();

  if( nGoodVertices == 0 ){
    cout << "Didn't find any good vertices!" << endl;
    return false;
  }

  float betamax = -1;
  int   imax    = -1;

  for (vector<int>::iterator ivtx = myGoodVertices.begin(); ivtx != myGoodVertices.end() ; ++ivtx ){
    
    float beta = beta_jet_vtx( ijet , *ivtx , beta_exponent );

    if( beta > betamax ){
      betamax = beta;
      imax    = *ivtx;
    }

  }
  
  if( imax == vtxIdx )  return true;
  return false;
}

//--------------------------------------------------------------------

float looper::dz_trk_vtx(const unsigned int trkidx, const unsigned int vtxidx){
  
  return ((cms2.trks_vertex_p4()[trkidx].z()-cms2.vtxs_position()[vtxidx].z()) - ((cms2.trks_vertex_p4()[trkidx].x()-cms2.vtxs_position()[vtxidx].x()) * cms2.trks_trk_p4()[trkidx].px() + (cms2.trks_vertex_p4()[trkidx].y() - cms2.vtxs_position()[vtxidx].y()) * cms2.trks_trk_p4()[trkidx].py())/cms2.trks_trk_p4()[trkidx].pt() * cms2.trks_trk_p4()[trkidx].pz()/cms2.trks_trk_p4()[trkidx].pt());
  
}

//--------------------------------------------------------------------































//--------------------------------------------------------------------

// std::pair<float,float> looper::HooberMET(const unsigned int vtxIdx, const unsigned int hypIdx, float dz_thresh, 
//                                          float pt_thresh, float etacut , bool usePFCandidatePt ){
  
//   float tmet_x = 0.;
//   float tmet_y = 0.;
  
//   //------------------------------------------------
//   // start by adding hypothesis leptons to tmet
//   //------------------------------------------------
  
//   if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
//     tmet_x -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].px();
//     tmet_y -= cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].py();

//     //cout << "Correcting for electron (pt,phi) (" << cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].pt() 
//     //     << " , " << cms2.els_p4()[cms2.hyp_lt_index()[hypIdx]].phi() << endl; 
//   }
//   else if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
//     tmet_x -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].px();
//     tmet_y -= cms2.mus_p4()[cms2.hyp_lt_index()[hypIdx]].py();
//   }
  
//   if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
//     tmet_x -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].px();
//     tmet_y -= cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].py();

//     //cout << "Correcting for electron (pt,phi) (" << cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].pt() 
//     //     << " , " << cms2.els_p4()[cms2.hyp_ll_index()[hypIdx]].phi() << endl; 

//   }
//   else if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
//     tmet_x -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].px();
//     tmet_y -= cms2.mus_p4()[cms2.hyp_ll_index()[hypIdx]].py();
//   }

//   //---------------------------------------------------
//   // loop over PFCandidates
//   //---------------------------------------------------
  
//   for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ipf++) {
  
//     //--------------------
//     //deal with neutrals
//     //--------------------
    
//     if( cms2.pfcands_charge().at(ipf) == 0 ){

//       //--------------------
//       // pt, eta cuts
//       //--------------------
      
//       if( pfcands_p4().at(ipf).pt()  < pt_thresh         ) continue;
//       if( fabs( pfcands_p4().at(ipf).eta() ) > etacut    ) continue;
      
//       if (dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_lt_p4()[hypIdx]) < 0.1)  continue;
//       if (dRbetweenVectors(pfcands_p4().at(ipf) , cms2.hyp_ll_p4()[hypIdx]) < 0.1)  continue;

//       //--------------------------------------------------
//       // for photons, require dr(photon,electron) > 0.1
//       //--------------------------------------------------
//       /*
//         float dr = 100.;
       
//         if( pfcands_particleId().at(ipf) == 22 ){
        
//         if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
//         float dr_lt = dRbetweenVectors( pfcands_p4().at(ipf) , cms2.hyp_lt_p4()[hypIdx] );
//         if( dr_lt < dr ) dr = dr_lt;
//         }
        
//         if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
//         float dr_ll = dRbetweenVectors( pfcands_p4().at(ipf) , cms2.hyp_ll_p4()[hypIdx] );
//         if( dr_ll < dr ) dr = dr_ll;
//         }
        
//         if( dr < 0.1 ) continue;
//         }
//       */
//       //-----------------------------------------
//       // correct MET for neutral PFCandidate
//       //-----------------------------------------
      
//       tmet_x -= pfcands_p4().at(ipf).px();
//       tmet_y -= pfcands_p4().at(ipf).py();
      
//       //cout << "Correcting for neutral PFCandidate (pt,phi) (" << cms2.pfcands_p4()[ipf].pt() 
//       //     << " , " << cms2.pfcands_p4()[ipf].phi() << " id " << pfcands_particleId()[ipf] << endl; 
      
//     }
    
//     //-------------------------------
//     // deal with charged particles
//     //-------------------------------

//     else{

//       //------------------------------------
//       // get track matched to PFCandidate
//       //------------------------------------
      
//       int itrk = cms2.pfcands_trkidx().at(ipf);
      
//       if( itrk >= trks_trk_p4().size() || itrk < 0 ){
//         //note: this should only happen for electrons which do not have a matched track
//         //currently we are just ignoring these guys
//         continue;
//       }
      
//       //----------------------------------------
//       // find closest PV and dz w.r.t. that PV
//       //----------------------------------------
      
//       float mindz = 999.;
//       int vtxi    = -1;
        
//       for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
//         float mydz = dz_trk_vtx(itrk,ivtx);
//         //float mydz = cms2.vtxs_position()[ivtx].z() - cms2.trks_vertex_p4()[itrk].z();
//         if (fabs(mydz) < fabs(mindz)) {
//           mindz = mydz;
//           vtxi = ivtx;
//         }
//       }
      
//       //----------------------------------------------------------------------------
//       // require closest PV is signal PV, dz cut, exclude tracks near hyp leptons
//       //----------------------------------------------------------------------------
      
//       if (vtxi != vtxIdx)
//         continue;
//       if (fabs(mindz) > dz_thresh)
//         continue;
//       if (dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_lt_p4()[hypIdx]) < 0.03)
//         continue;
//       if (dRbetweenVectors(cms2.trks_trk_p4()[itrk], cms2.hyp_ll_p4()[hypIdx]) < 0.03)
//         continue;
      

//       //---------------------------------------
//       // correct MET for charged PFCandidate
//       //---------------------------------------

//       if( usePFCandidatePt ){
//         tmet_x -= cms2.pfcands_p4()[ipf].px();
//         tmet_y -= cms2.pfcands_p4()[ipf].py();

//         //cout << "Correcting for charged PFCandidate (pt,phi) (" << cms2.pfcands_p4()[ipf].pt() 
//         //     << " , " << cms2.pfcands_p4()[ipf].phi() << endl; 
//       }
      
//       else{
//         tmet_x -= cms2.trks_trk_p4()[itrk].px();
//         tmet_y -= cms2.trks_trk_p4()[itrk].py();
//       }
      
//     } 
//   }// end loop over tracks

  
//   float met = sqrt(tmet_x * tmet_x + tmet_y * tmet_y);
//   float metphi = atan2( tmet_y , tmet_x );
  
//   return make_pair( met , metphi );

// }
