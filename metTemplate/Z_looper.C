#include "Z_looper.h"
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
#include "CORE/metSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/muonSelections.h"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"


bool debug = false;
float lumi = 0.0027945;
char* iter = "2p8pb_zlooper_v2";

//float lumi = 0.001932;
//char* iter = "1p9pb_zlooper_v2";

const int nJetBins      = 11;
const int nSumJetPtBins = 23;

using namespace tas;

int getJetBin( int njets ){
  
  int bin = njets;
  if( bin >= nJetBins ) bin = nJetBins - 1;
  
  return bin;
}

string jetString( int bin ){

  stringstream s;
  s << "nJets = " << bin ;

  return s.str();

}

int getSumJetPtBin( float x ){

  //bins array defines the sumJetPt binning
  float bins[nSumJetPtBins+1]={0,25,50,75,100,125,150,175,200,250,300,350,400,450,
                               500,600,700,800,900,1000,2000,3000,4000,5000};
  
  if( x < bins[0] )              return 0;
  if( x >= bins[nSumJetPtBins] ) return nSumJetPtBins - 1;

  int ptbin = -1;

  for( int ibin = 0 ; ibin < nSumJetPtBins+1 ; ibin++){
    if( x >= bins[ibin] && x< bins[ibin+1] ){
      ptbin = ibin;
      break;
    }
  }

  if( ptbin == -1 ) 
    cout << "ERROR CANNOT FIND BIN FOR SUMJETPT " << x << endl;

  return ptbin;
}

string sumJetPtString( int bin ){
  
  float bins[nSumJetPtBins+1]={0,25,50,75,100,125,150,175,200,250,300,350,400,450,
                               500,600,700,800,900,1000,2000,3000,4000,5000};

  stringstream s;
  s << bins[bin] << " < sumJetPt < " << bins[bin+1] << " GeV";

  return s.str();
}

void Z_looper::ScanChain (TChain* chain, const char* prefix, bool isData, bool calculateTCMET, metAlgo algo, int nEvents, float kFactor){

  algo_ = algo;

  if( algo_ == e_makeTemplate )    cout << "metAlgo makeTemplate" << endl;
  if( algo_ == e_photonSelection ) cout << "metAlgo photonSelection" << endl;
  if( algo_ == e_ZSelection )      cout << "metAlgo ZSelection" << endl;

  TFile *metTemplateFile;
  string metTemplateString = "";
  
  DorkyEventIdentifier dei;

  //select templates
  if( algo_ != e_makeTemplate){
    if( isData ){
      metTemplateString = "_dataTemplate";
      metTemplateFile = TFile::Open("root/JetMETTau_250nb.root");
    }else{
      metTemplateString = "_mcTemplate";
      metTemplateFile = TFile::Open("root/QCD_Pt15_1p9pb.root");
    }
  }

  //set_goodrun_file("Cert_TopAug13_Merged_135059-142664_goodruns.txt");
  //set_goodrun_file("Cert_TopAug25_Merged_135059-143336_goodruns.txt");
  //set_goodrun_file("Cert_TopAug26_Merged_135059-143835_recover_noESDCS_goodruns.txt");
  set_goodrun_file("Cert_TopAug30_Merged_135059-144114_recover_noESDCS_goodruns.txt");

  bookHistos();

  // make a baby ntuple
  stringstream babyfilename;
  babyfilename << prefix << "_baby.root";
  MakeBabyNtuple( Form("root/%s_%s_baby.root", prefix , iter ) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  //pass fail counters
  int nPassGoodRun    = 0;
  int nPassBPTX       = 0;
  int nPassBSC        = 0;
  int nPassBeamHalo   = 0;
  int nPassGoodTracks = 0;
  int nPassGoodVertex = 0;
  int nSkip_els_conv_dist = 0;
  if(debug) cout << "Begin file loop" << endl;

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())){
    
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();
 
    for (unsigned int event = 0 ; event < nEvents; ++event){
      
      cms2.GetEntry(event);
      ++nEventsTotal;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
            
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //duplicate event veto CHECK THIS
      bool myduplicate = dei.is_duplicate(DorkyEvent());
      if(myduplicate) continue;
      
      //skip events with bad els_conv_dist 
      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
      }
      
      if( skipEvent ){
        cout << "SKIPPING EVENT WITH BAD ELS_CONV_DIST" << endl;       
        nSkip_els_conv_dist++;
        continue;
      }

      //basic event selection
      if (!isData || goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))
        nPassGoodRun++;
      else continue;
          
      bptx_        =  cleaning_BPTX(isData)             ? 1 : 0;
      bsc_         =  cleaning_BSC()                    ? 1 : 0;
      beamhalo_    =  cleaning_beamHalo()               ? 1 : 0;
      goodvtx_     =  cleaning_goodVertexAugust2010()   ? 1 : 0;
      goodtrks_    =  cleaning_goodTracks()             ? 1 : 0;

//       // determine if current event passes BPTX triggers
//       if (cleaning_BPTX( isData ))
//         nPassBPTX++;
//       else continue;

//       // determine if current event passes BSC triggers
//       if (cleaning_BSC())
//         nPassBSC++;
//       else continue;

//       // determine if current event passes beam halo triggers
//       if (cleaning_beamHalo())
//         nPassBeamHalo++;
//       else continue;

//       // determine if current event is a beam scraping event
//       if (cleaning_goodTracks())
//         nPassGoodTracks++;
//       else continue;

//       // determine if current event has a good vertex
//       if (cleaning_goodVertexAugust2010())
//         nPassGoodVertex++;
//       else continue;
      
      if(debug) cout << "Pass event selection" << endl;

      InitBabyNtuple();

      // event stuff
      run_    = cms2.evt_run();
      lumi_   = cms2.evt_lumiBlock();
      event_  = cms2.evt_event();

      weight_ = 1.;
      pthat_  = -1;
      if( !isData ){
        weight_ = cms2.evt_scale1fb() * kFactor * lumi;
        pthat_  = cms2.genps_pthat();
      }

      //access HLT triggers------------------------------------------------------------------

      for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_L1Jet6U" ) == 0 ){
          if( passHLTTrigger("HLT_L1Jet6U") )                  HLT_L1Jet6U_ = 1;
          else                                                 HLT_L1Jet6U_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_L1Jet10U" ) == 0 ){
          if( passHLTTrigger("HLT_L1Jet10U") )                 HLT_L1Jet10U_ = 1;
          else                                                 HLT_L1Jet10U_ = 0;
        }
        
        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Jet15U" ) == 0 ){
          if( passHLTTrigger("HLT_Jet15U") )                   HLT_Jet15U_ = 1;
          else                                                 HLT_Jet15U_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Jet30U" ) == 0 ){
          if( passHLTTrigger("HLT_Jet30U") )                   HLT_Jet30U_ = 1;
          else                                                 HLT_Jet30U_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "L1_SingleEG5" ) == 0 ){
          if( passHLTTrigger("L1_SingleEG5") )                 L1_SingleEG5_ = 1;
          else                                                 L1_SingleEG5_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon10_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon10_L1R") )             HLT_Photon10_L1R_ = 1;
          else                                                 HLT_Photon10_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon15_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon15_L1R") )             HLT_Photon15_L1R_ = 1;
          else                                                 HLT_Photon15_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon10_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon10_Cleaned_L1R") )     HLT_Photon10_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon10_Cleaned_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon15_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon15_Cleaned_L1R") )     HLT_Photon15_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon15_Cleaned_L1R_ = 0;
        }

        if( strcmp( hlt_trigNames().at(itrig) , "HLT_Photon20_Cleaned_L1R" ) == 0 ){
          if( passHLTTrigger("HLT_Photon20_Cleaned_L1R") )     HLT_Photon20_Cleaned_L1R_ = 1;
          else                                                 HLT_Photon20_Cleaned_L1R_ = 0;
        }
      }      
      //met quantities--------------------------------------------------------------------------
      
      //calomet
      met_       = cms2.evt_met();
      metphi_    = cms2.evt_metPhi();
      sumet_     = cms2.evt_sumet();

      // pf met stuff
      pfmet_    = cms2.evt_pfmet();
      pfmetphi_ = cms2.evt_pfmetPhi();
      pfsumet_  = cms2.evt_pfsumet();

      //muon-corrected met stuff
      mumet_    = cms2.evt_metMuonCorr();
      mumetphi_ = cms2.evt_metMuonCorrPhi();
      musumet_  = cms2.evt_sumetMuonCorr();

      //muon-corrected JES met stuff
      mujesmet_    = cms2.evt_metMuonJESCorr();
      mujesmetphi_ = cms2.evt_metMuonJESCorrPhi();
      mujessumet_  = -9999.; //cms2.evt_sumetMuonJESCorr(); //branch doesn't exist!!!

      // genmet stuff
      if (!isData){
        genmet_     = cms2.gen_met();
        genmetphi_  = cms2.gen_metPhi();
        gensumet_   = cms2.gen_sumEt();
      }

      //Z_looper-level tcmet
      if( calculateTCMET ){
        metStruct tcmetNewStruct = correctedTCMET();
        tcmetNew_     = tcmetNewStruct.met;
        tcmetphiNew_  = tcmetNewStruct.metphi;
        tcsumetNew_   = tcmetNewStruct.sumet;
      }

      //just take tcmet from event
      tcmet_     = evt_tcmet();
      tcmetphi_  = evt_tcmetPhi();
      tcsumet_   = evt_tcsumet();

      if( hyp_p4().size() != hyp_type().size() ){
        cout << "ERROR hyp_p4().size() " << hyp_p4().size() << " hyp_type().size() " << hyp_type().size() << endl;
      }

      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); ++hypIdx) {
   
        //loop over pfjets, find nJets and sumJetPt
        nJets_        = 0;
        sumJetPt_     = 0.;
        nJets40_      = 0;
        sumJetPt10_   = 0.;
        
        LorentzVector jetSystem(0.,0.,0.,0.);
        int imaxjet = -1;
        float maxpt = -1;
        
        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          
          if (dRbetweenVectors(vjet, vll) < 0.4)  continue;
          if (dRbetweenVectors(vjet, vlt) < 0.4)  continue;
          if( fabs( vjet.eta() ) > 5.)            continue;

          if( pfjets_p4().at(ijet).pt() > maxpt ){
            maxpt = pfjets_p4().at(ijet).pt();
            imaxjet = ijet;
          }
            
          if ( vjet.pt() > 30. ){
            nJets_++;
          }        
          if ( vjet.pt() > 15. ){
            sumJetPt_ += vjet.pt();
            jetSystem += vjet;
          }
          
          if ( vjet.pt() > 40. ){
            nJets40_++;
          }
          if ( vjet.pt() > 10. ){
            sumJetPt10_ += vjet.pt();
          }
        }
      
        jetmax_pt_ = -1;

        if( imaxjet > -1 ){
          jetmax_pt_       = pfjets_p4().at(imaxjet).pt();
          jetmax_dphimet_  = deltaPhi( pfjets_p4().at(imaxjet).phi() , tcmetphi_);
        }

        vecJetPt_ = jetSystem.pt();

        //passz_ = passZSelection( hypIdx ) ? 1 : 0;
        //if( passz_ == 0 ) continue;

        
        if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                   continue;
        if( hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  continue;
        if( hyp_ll_p4()[hypIdx].pt() < 10. )                                  continue;
        if( hyp_lt_p4()[hypIdx].pt() < 10. )                                  continue;
                   
        leptype_ = 99;
        if (hyp_type()[hypIdx] == 3) leptype_ = 0;                           // ee
        if (hyp_type()[hypIdx] == 0) leptype_ = 1;                           // mm
        if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) leptype_=2;  // em
        if (leptype_ == 99) {
          cout << "Skipping unknown dilepton type = " << hyp_type()[hypIdx] << endl;
          continue;
        }
        
        if( abs( hyp_ll_id()[hypIdx] ) == 11 ){
          passe_ll_ttbarV1_      = pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV1 , isData ) ? 1 : 0;
          passe_ll_ttbar_        = pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbar   , isData ) ? 1 : 0;
          passe_ll_cand01_       = pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_cand01 )           ? 1 : 0;
        }else if( abs( hyp_ll_id()[hypIdx] ) == 13 ){
          flagll_                = mus_tcmet_flag().at(hyp_ll_index()[hypIdx]);
          passm_ll_nom_          = muonId(hyp_ll_index()[hypIdx])                ? 1 : 0;
          passm_ll_nomttbar_     = muonId(hyp_ll_index()[hypIdx],NominalTTbar)   ? 1 : 0;
          passm_ll_nomttbarV2_   = muonId(hyp_ll_index()[hypIdx],NominalTTbarV2) ? 1 : 0;
        }else{
          cout << "ERROR UNRECOGNIZED LL ID " << hyp_ll_id()[hypIdx] << endl;
          continue;
        }
        
        if( abs( hyp_lt_id()[hypIdx] ) == 11 ){
          passe_lt_ttbarV1_      = pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV1 , isData ) ? 1 : 0;
          passe_lt_ttbar_        = pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbar   , isData ) ? 1 : 0;
          passe_lt_cand01_       = pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_cand01 )           ? 1 : 0;
        }else if( abs( hyp_lt_id()[hypIdx] ) == 13 ){
          flaglt_                = mus_tcmet_flag().at(hyp_lt_index()[hypIdx]);
          passm_lt_nom_          = muonId(hyp_lt_index()[hypIdx])                ? 1 : 0;
          passm_lt_nomttbar_     = muonId(hyp_lt_index()[hypIdx],NominalTTbar)   ? 1 : 0;
          passm_lt_nomttbarV2_   = muonId(hyp_lt_index()[hypIdx],NominalTTbarV2) ? 1 : 0;
        }else{
          cout << "ERROR UNRECOGNIZED LT ID " << hyp_lt_id()[hypIdx] << endl;
          continue;
        }

        
        ptll_             = hyp_ll_p4()[hypIdx].pt();
        ptlt_             = hyp_lt_p4()[hypIdx].pt();
        etall_            = hyp_ll_p4()[hypIdx].eta();
        etalt_            = hyp_lt_p4()[hypIdx].eta();
        dilmass_          = hyp_p4()[hypIdx].mass(); 
        dilpt_            = hyp_p4()[hypIdx].pt(); 

        metStruct tcmetStruct = correctTCMETforHypMuons( hypIdx ,  
                                                         evt_tcmet() * cos( evt_tcmetPhi() ), 
                                                         evt_tcmet() * sin( evt_tcmetPhi() ), 
                                                         evt_tcsumet() );
        
        // out-of-the-box  tcmet stuff (corrected for hyp muons)
        tcmet_    = tcmetStruct.met;
        tcmetphi_ = tcmetStruct.metphi;
        tcsumet_  = tcmetStruct.sumet;   

        metStruct tcmetNewStruct = correctTCMETforHypMuons( hypIdx ,  
                                                            tcmetNew_ * cos( tcmetphiNew_ ), 
                                                            tcmetNew_ * sin( tcmetphiNew_ ), 
                                                            tcsumetNew_ );
            
        // Z_looper-level  tcmet stuff (corrected for hyp muons)
        tcmetNew_    = tcmetNewStruct.met;
        tcmetphiNew_ = tcmetNewStruct.metphi;
        tcsumetNew_  = tcmetNewStruct.sumet;   

        nGoodVertex_ = 0;
        for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
          if(isGoodVertex(v)) nGoodVertex_++;
        }
        
        //fill baby ntuple
        if ( nJets_ < 1 )     continue;
        if( jetmax_pt_ < 30 ) continue;
        FillBabyNtuple();

        //selection
        if ( nJets_ < 2 )                                                     continue;
        if( jetmax_pt_ < 50 )                                                 continue;
        //if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                   continue;
        //if( hyp_type()[hypIdx]==1 || hyp_type()[hypIdx]==2)                   continue;
        //if( hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  continue;
        //if( hyp_ll_p4()[hypIdx].pt() < 10 )                                   continue;
        //if( hyp_lt_p4()[hypIdx].pt() < 10 )                                   continue;
     
        //ttbar muon ID
        if (abs(hyp_ll_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx],NominalTTbarV2))))   continue;
        if (abs(hyp_lt_id()[hypIdx]) == 13  && (!(fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx],NominalTTbarV2))))   continue;
        
        //ttbarV1 electron ID
        if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV1 , isData ))) continue;
        if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV1 , isData ))) continue;
 
        dphixmet_  = deltaPhi( tcmetphi_ , hyp_p4()[hypIdx].phi() );
        metPar_    = tcmet_ * cos( dphixmet_ );
        metPerp_   = tcmet_ * sin( dphixmet_ );
        
        //fill predicted and observed met histos
        int iJetBin      = getJetBin( nJets_ );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        TH1F* hmet     = (TH1F*) metTemplateFile->Get(Form("metTemplate_%i_%i",iJetBin,iSumJetPtBin));
        TH1F* hmetPar  = (TH1F*) metTemplateFile->Get(Form("metParTemplate_%i_%i",iJetBin,iSumJetPtBin));
        TH1F* hmetPerp = (TH1F*) metTemplateFile->Get(Form("metPerpTemplate_%i_%i",iJetBin,iSumJetPtBin));
        
        fillUnderOverFlow( metObserved_njets[iJetBin]  ,  tcmet_ , weight_ );
        metPredicted_njets[iJetBin]->Add( hmet );
        
        fillUnderOverFlow( metObserved , tcmet_ , weight_  );
        metPredicted->Add( hmet );

        if( leptype_ == 0 || leptype_ == 1 ){
          fillUnderOverFlow( metObserved_sf , tcmet_ , weight_  );
          metPredicted_sf->Add( hmet );
        }else if( leptype_ == 2 ){
          fillUnderOverFlow( metObserved_df , tcmet_ , weight_  );
          metPredicted_df->Add( hmet );
        }else{
          cout << "UNRECOGNIZED LEPTYPE " << leptype_ << endl;
        }
        
        fillUnderOverFlow( metParObserved , metPar_ , weight_  );
        metParPredicted->Add( hmetPar );
        
        fillUnderOverFlow( metPerpObserved ,  metPerp_ , weight_ );
        metPerpPredicted->Add( hmetPerp );
        
        delete hmet;
        delete hmetPar;
        delete hmetPerp;

      }// end loop over hypIdx
    } // end loop over events
  } // end loop over files
  
  cout << "\n\n********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events that pass good run/lumi: " << nPassGoodRun 
       << " (" << 100*(double)nPassGoodRun/nEventsTotal << "% of total)" << endl;
  cout << "Total number of events that pass BPTX trigger: " << nPassBPTX
       << " (" << 100*(double)nPassBPTX/nPassGoodRun << "%)" << endl;
  cout << "Total number of events that pass BSC trigger: " << nPassBSC
       << " (" << 100*(double)nPassBSC/nPassBPTX << "%)" << endl;
  cout << "Total number of events that pass BeamHalo trigger: " << nPassBeamHalo
       << " (" << 100*(double)nPassBeamHalo/nPassBSC << "%)" << endl;
  cout << "Total number of events that pass tracking cuts: " << nPassGoodTracks
       << " (" << 100*(double)nPassGoodTracks/nPassBeamHalo << "%)" << endl;
  cout << "Total number of events that pass vertex cuts: " << nPassGoodVertex
       << " (" << 100*(double)nPassGoodVertex/nPassGoodTracks << "%)" << endl;
  cout << endl << endl;

  if( nSkip_els_conv_dist > 0 ){
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist branch" << endl;
  }

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  CloseBabyNtuple();

  // make histos rootfile
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form("root/%s_%s%s.root", prefix , iter , metTemplateString.c_str() ) );
  deleteHistos();
  
} // end ScanChain

float Z_looper::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

bool Z_looper::passMuon_NominalTTbar( int hypIdx ){

  //ttbar muon ID
  if (abs(hyp_ll_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx],NominalTTbar))))   return false;
  if (abs(hyp_lt_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx],NominalTTbar))))   return false;
  return true;

}

bool Z_looper::passMuon_NominalTTbarV2( int hypIdx ){

  //ttbar muon ID
  if (abs(hyp_ll_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx],NominalTTbarV2))))   return false;
  if (abs(hyp_lt_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx],NominalTTbarV2))))   return false;
  return true;

}

bool Z_looper::passMuon_Nominal( int hypIdx ){

  //nominal muon ID
  if (abs(hyp_ll_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx]))))   return false;
  if (abs(hyp_lt_id()[hypIdx]) == 13  && (!(fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx]))))   return false;
  return true;

}

bool Z_looper::passElectron_ttbarV1( bool isData ,  int hypIdx ){

  //ttbar electron ID
  if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV1 , isData ))) return false;
  if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV1 , isData ))) return false;
  return true;

}

bool Z_looper::passElectron_ttbar( bool isData , int hypIdx ){

  //ttbar electron ID
  if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbar , isData ))) return false;
  if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbar , isData ))) return false;
  return true;

}

bool Z_looper::passElectron_cand01( int hypIdx ){

  //cand01 electron ID
  if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_cand01 ))) return false;
  if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_cand01 ))) return false;
  return true;

}

bool Z_looper::passZSelection( int hypIdx ){
  
  if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                   return false;
  if( hyp_type()[hypIdx]==1 || hyp_type()[hypIdx]==2)                   return false;
  if( hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  return false;

  return true;

}

void Z_looper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

void Z_looper::InitBabyNtuple (){

  //triggers
  HLT_L1Jet6U_          = -1;
  HLT_L1Jet10U_         = -1;
  HLT_Jet15U_           = -1;
  HLT_Jet30U_           = -1;
  L1_SingleEG5_         = -1;
  HLT_Photon10_L1R_     = -1;
  HLT_Photon15_L1R_     = -1;
  HLT_Photon10_Cleaned_L1R_     = -1;
  HLT_Photon15_Cleaned_L1R_     = -1;
  HLT_Photon20_Cleaned_L1R_     = -1;

  // event stuff
  run_          = -999999;
  lumi_         = -999999;
  event_        = -999999;
  weight_       = -999999.;
  pthat_        = -999999.;
  nGoodVertex_  = -999999;
  leptype_      = -999999;

  // genmet stuff
  genmet_     = -999999.;
  genmetphi_  = -999999.;
  gensumet_   = -999999.;

  // pfmet stuff
  pfmet_     = -999999.;
  pfmetphi_  = -999999.;
  pfsumet_   = -999999.;

  // calomet stuff
  met_          = -999999.;
  metphi_       = -999999.;
  sumet_        = -999999.;

  // muon-corrected calomet stuff
  mumet_        = -999999.;
  mumetphi_     = -999999.;
  musumet_      = -999999.;

  // calomet stuff
  mujesmet_     = -999999.;
  mujesmetphi_  = -999999.;
  mujessumet_   = -999999.;

  // tcmet stuff
  dphixmet_     = -999999.;
  metPar_       = -999999.;
  metPerp_      = -999999.;

  tcmet_        = -999999.;
  tcmetphi_     = -999999.;
  tcsumet_      = -999999.;

  tcmetNew_     = -999999.;
  tcsumetNew_   = -999999.;
  tcmetphiNew_  = -999999.;

  nJets_        = -999999;
  sumJetPt_     = -999999;
  vecJetPt_     = -999999;
  nJets40_      = -999999;
  sumJetPt10_   = -999999;

  //photon stuff
  nPhotons_ = -999999;
  etg_      = -999999.;
  etag_     = -999999.;
  phig_     = -999999.;
  hoe_      = -999999.;
  eciso_    = -999999.;
  hciso_    = -999999.;
  tkiso_    = -999999.;
  swiss_    = -999999.;
  seed_     = -999999.;
  s4_       = -999999.;
  r4_       = -999999.;
  drel_     = -999999.;

  photon_pixelseed_        = -999999;
  photon_e15_              = -999999.;
  photon_e25max_           = -999999.;
  photon_e33_              = -999999.;
  photon_e55_              = -999999.;
  photon_ecalIso03_        = -999999.;
  photon_ecalIso04_        = -999999.;
  photon_hcalIso03_        = -999999.;
  photon_hcalIso04_        = -999999.;
  photon_ntkIsoHollow03_   = -999999.;
  photon_ntkIsoHollow04_   = -999999.;
  photon_ntkIsoSolid03_    = -999999.;
  photon_ntkIsoSolid04_    = -999999.;
  photon_sigmaEtaEta_      = -999999.;
  photon_sigmaIEtaIEta_    = -999999.;
  photon_tkisoHollow03_    = -999999.;
  photon_tkisoHollow04_    = -999999.;
  photon_tkisoSolid03_     = -999999.;
  photon_tkisoSolid04_     = -999999.;
                                  
  //photon-matched jet stuff
  jet_eta_          = -999999.;  
  jet_energy_       = -999999.;  
  jet_pt_           = -999999.;  
  jet_chg_emfrac_   = -999999.;  
  jet_chg_hadfrac_  = -999999.;  
  jet_neu_emfrac_   = -999999.;  
  jet_neu_hadfrac_  = -999999.;  
  jet_nchg_         = -999999;      
  jet_nmuon_        = -999999;  
  jet_nneu_         = -999999;  
  jet_dphimet_      = -999999.;  
  jet_dpt_          = -999999.;  
  jet_drgen_        = -999999.;  

  //leading jet stuff
  jetmax_pt_        = -999999;
  jetmax_dphimet_   = -999999;

  //Z stuff
  passz_           = -999999;
  passe_ll_ttbar_     = -999999;
  passe_ll_ttbarV1_   = -999999;
  passe_ll_cand01_    = -999999;
  passm_ll_nomttbar_  = -999999;
  passm_ll_nomttbarV2_= -999999;
  passm_ll_nom_       = -999999;
  passe_lt_ttbar_     = -999999;
  passe_lt_ttbarV1_   = -999999;
  passe_lt_cand01_    = -999999;
  passm_lt_nomttbar_  = -999999;
  passm_lt_nomttbarV2_= -999999;
  passm_lt_nom_       = -999999;
  pdgid_           = -999999;
  ptll_            = -999999;
  ptlt_            = -999999;
  etall_           = -999999;
  etalt_           = -999999;
  dilmass_         = -999999.;
  dilpt_           = -999999.;
  flagll_          = -999999;
  flaglt_          = -999999;

  bptx_        =  -999999;
  bsc_         =  -999999;
  beamhalo_    =  -999999;
  goodvtx_     =  -999999;
  goodtrks_    =  -999999;

}

void Z_looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hgenps_pthat = new TH1F("hgenps_pthat","",100,0,100);
  hphotonpt    = new TH1F("hphotonpt","",100,0,100);

  hgenps_pthat->GetXaxis()->SetTitle("gen p_{T}(hat) (GeV)");
  hphotonpt->GetXaxis()->SetTitle("max photon p_{T} (GeV)");

  metObserved  = new TH1F("metObserved", "Observed MET",500,0,500);
  metPredicted = new TH1F("metPredicted","Predicted MET",500,0,500);

  metObserved_sf  = new TH1F("metObserved_sf", "Observed MET (SF)",500,0,500);
  metPredicted_sf = new TH1F("metPredicted_sf","Predicted MET (SF)",500,0,500);

  metObserved_df  = new TH1F("metObserved_df", "Observed MET (DF)",500,0,500);
  metPredicted_df = new TH1F("metPredicted_df","Predicted MET (DF)",500,0,500);

  metParObserved  = new TH1F("metParObserved", "Observed MET (Parallel)",1000,-500,500);
  metParPredicted = new TH1F("metParPredicted","Predicted MET (Parallel)",1000,-500,500);

  metPerpObserved  = new TH1F("metPerpObserved", "Observed MET (Perpendicular)",500,0,500);
  metPerpPredicted = new TH1F("metPerpPredicted","Predicted MET (Perpendicular)",500,0,500);
 
  metObserved->Sumw2();
  metPredicted->Sumw2();

  metParObserved->Sumw2();
  metParPredicted->Sumw2();

  metPerpObserved->Sumw2();
  metPerpPredicted->Sumw2();

  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){

    metObserved_njets[iJetBin]  = new TH1F(Form("metObserved_njets%i",iJetBin), Form("Observed MET NJets %i", iJetBin),500,0,500);
    metPredicted_njets[iJetBin] = new TH1F(Form("metPredicted_njets%i",iJetBin),Form("Predicted MET NJets %i",iJetBin),500,0,500);
    
    metObserved_njets[iJetBin] ->Sumw2();
    metPredicted_njets[iJetBin]->Sumw2();
  }
  

  if( algo_ == e_makeTemplate ) {
    
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
        metTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("metTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                          Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);


        metParTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("metParTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                             Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),1000,-500,500);

        metPerpTemplate[ iJetBin ][ iSumJetPtBin ] = new TH1F(Form("metPerpTemplate_%i_%i",iJetBin,iSumJetPtBin),
                                                              Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
       
        metTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
        metParTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 
        metPerpTemplate[ iJetBin ][ iSumJetPtBin ]->Sumw2(); 

        metTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("met");
        metParTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("met_{#parallel}");
        metPerpTemplate[ iJetBin ][ iSumJetPtBin ]->GetXaxis()->SetTitle("met_{#perp}");
      }
    }
  }
  
}


void Z_looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("T1", "A Baby Ntuple");

  babyTree_->Branch("run",          &run_,          "run/I"  );
  babyTree_->Branch("lumi",         &lumi_,         "lumi/I" );
  babyTree_->Branch("event",        &event_,        "event/I");
  babyTree_->Branch("leptype",      &leptype_,      "leptype/I");
  babyTree_->Branch("nvtx",         &nGoodVertex_,  "nvtx/I");
  babyTree_->Branch("weight",       &weight_,       "weight/F");
  babyTree_->Branch("pthat",        &pthat_,        "pthat/F");
  babyTree_->Branch("pfmet",        &pfmet_,        "pfmet/F"   );
  babyTree_->Branch("pfmetphi",     &pfmetphi_,     "pfmetphi/F");
  babyTree_->Branch("pfsumet",      &pfsumet_,      "pfsumet/F" );
  babyTree_->Branch("met",          &met_,          "met/F"      );
  babyTree_->Branch("metphi",       &metphi_,       "metphi/F"   );
  babyTree_->Branch("sumet",        &sumet_,        "sumet/F"    );
  babyTree_->Branch("mumet",        &mumet_,        "mumet/F"      );
  babyTree_->Branch("mumetphi",     &mumetphi_,     "mumetphi/F"   );
  babyTree_->Branch("musumet",      &musumet_,      "musumet/F"    );
  babyTree_->Branch("mujesmet",     &mujesmet_,     "mujesmet/F"      );
  babyTree_->Branch("mujesmetphi",  &mujesmetphi_,  "mujesmetphi/F"   );
  babyTree_->Branch("mujessumet",   &mujessumet_,   "mujessumet/F"    );
  babyTree_->Branch("genmet",       &genmet_,       "genmet/F"   );
  babyTree_->Branch("genmetphi",    &genmetphi_,    "genmetphi/F");
  babyTree_->Branch("gensumet",     &gensumet_,     "gensumet/F" );
  babyTree_->Branch("dphixmet",     &dphixmet_,      "dphixmet/F"    );
  babyTree_->Branch("metpar",       &metPar_,        "metpar/F"      );
  babyTree_->Branch("metperp",      &metPerp_,       "metperp/F"     );

  babyTree_->Branch("tcmet",        &tcmet_,        "tcmet/F"      );
  babyTree_->Branch("tcmetphi",     &tcmetphi_,     "tcmetphi/F"   );
  babyTree_->Branch("tcsumet",      &tcsumet_,      "tcsumet/F"    );

  babyTree_->Branch("tcmetNew",     &tcmetNew_,     "tcmetNew/F"      );
  babyTree_->Branch("tcmetphiNew",  &tcmetphiNew_,  "tcmetphiNew/F"   );
  babyTree_->Branch("tcsumetNew",   &tcsumetNew_,   "tcsumetNew/F"    );

  babyTree_->Branch("njets",        &nJets_,        "njets/I"       );
  babyTree_->Branch("njets40",      &nJets40_,      "njets40/I"     );
  babyTree_->Branch("sumjetpt",     &sumJetPt_,     "sumjetpt/F"    );
  babyTree_->Branch("sumjetpt10",   &sumJetPt10_,   "sumjetpt10/F"    );
  babyTree_->Branch("vecjetpt",     &vecJetPt_,     "vecjetpt/F"    );

  //photon stuff
  babyTree_->Branch("ng",      &nPhotons_, "ng/I"); 
  babyTree_->Branch("etg",     &etg_,      "etg/F");	   
  babyTree_->Branch("etag",    &etag_,     "etag/F");	   
  babyTree_->Branch("phig",    &phig_,     "phig/F");
  babyTree_->Branch("hoe",     &hoe_,      "hoe/F");	   
  babyTree_->Branch("eciso",   &eciso_,    "eciso/F");	   
  babyTree_->Branch("hciso",   &hciso_,    "hciso/F");	   
  babyTree_->Branch("tkiso",   &tkiso_,    "tkiso/F");
  babyTree_->Branch("swiss",   &swiss_,    "swiss/F");
  babyTree_->Branch("seed",    &seed_,     "seed/F");
  babyTree_->Branch("s4",      &s4_,       "s4/F");
  babyTree_->Branch("r4",      &r4_,       "r4/F");

  babyTree_->Branch("photon_pixelseed",         &photon_pixelseed_,         "photon_pixelseed/I");         
  babyTree_->Branch("photon_e15",               &photon_e15_,               "photon_e15/F");                
  babyTree_->Branch("photon_e25max",            &photon_e25max_,            "photon_e25max/F");             
  babyTree_->Branch("photon_e33",               &photon_e33_,               "photon_e33/F");                
  babyTree_->Branch("photon_e55",               &photon_e55_,               "photon_e55/F");                
  babyTree_->Branch("photon_ecalIso03",         &photon_ecalIso03_,         "photon_ecalIso03/F");          
  babyTree_->Branch("photon_ecalIso04",         &photon_ecalIso04_,         "photon_ecalIso04/F");          
  babyTree_->Branch("photon_hcalIso03",         &photon_hcalIso03_,         "photon_hcalIso03/F");          
  babyTree_->Branch("photon_hcalIso04",         &photon_hcalIso04_,         "photon_hcalIso04/F");          
  babyTree_->Branch("photon_ntkIsoHollow03",    &photon_ntkIsoHollow03_,    "photon_ntkIsoHollow03/F");     
  babyTree_->Branch("photon_ntkIsoHollow04",    &photon_ntkIsoHollow04_,    "photon_ntkIsoHollow04/F");     
  babyTree_->Branch("photon_ntkIsoSolid03",     &photon_ntkIsoSolid03_,     "photon_ntkIsoSolid03/F");      
  babyTree_->Branch("photon_ntkIsoSolid04",     &photon_ntkIsoSolid04_,     "photon_ntkIsoSolid04/F");      
  babyTree_->Branch("photon_sigmaEtaEta",       &photon_sigmaEtaEta_,       "photon_sigmaEtaEta/F");        
  babyTree_->Branch("photon_sigmaIEtaIEta",     &photon_sigmaIEtaIEta_,     "photon_sigmaIEtaIEta/F");      
  babyTree_->Branch("photon_tkisoHollow03",     &photon_tkisoHollow03_,     "photon_tkisoHollow03/F");      
  babyTree_->Branch("photon_tkisoHollow04",     &photon_tkisoHollow04_,     "photon_tkisoHollow04/F");      
  babyTree_->Branch("photon_tkisoSolid03",      &photon_tkisoSolid03_,      "photon_tkisoSolid03/F");      
  babyTree_->Branch("photon_tkisoSolid04",      &photon_tkisoSolid04_,      "photon_tkisoSolid04/F");           
                                                                            
  //photon-matched jet stuff
  babyTree_->Branch("jetdr",                 &jet_dr_,               "jetdr/F");
  babyTree_->Branch("jetpt",                 &jet_pt_,               "jetpt/F");
  babyTree_->Branch("jeteta",                &jet_eta_,              "jeteta/F");
  babyTree_->Branch("jetenergy",             &jet_energy_,           "jetenergy/F");
  babyTree_->Branch("jetchargedemfrac",      &jet_chg_emfrac_,       "jetchargedemfrac/F");
  babyTree_->Branch("jetchargedhadfrac",     &jet_chg_hadfrac_,      "jetchargedhadfrac/F");
  babyTree_->Branch("jetneutralemfrac",      &jet_neu_emfrac_,       "jetneutralemfrac/F");
  babyTree_->Branch("jetneutralhadfrac",     &jet_neu_hadfrac_,      "jetneutralhadfrac/F");
  babyTree_->Branch("jetncharged",           &jet_nchg_,             "jetncharged/I");
  babyTree_->Branch("jetnmuon",              &jet_nmuon_,            "jetnmuon/I");
  babyTree_->Branch("jetnneutral",           &jet_nneu_,             "jetnneutral/I");
  babyTree_->Branch("jetdphimet",            &jet_dphimet_,          "jetdphimet/F");
  babyTree_->Branch("jetdpt",                &jet_dpt_,              "jetdpt/F");
  babyTree_->Branch("jetdrgen",              &jet_drgen_,            "jetdrgen/F");
  babyTree_->Branch("drel",                  &drel_,                 "drel/F");

  babyTree_->Branch("maxjetpt",              &jetmax_pt_,            "maxjetpt/F");
  babyTree_->Branch("maxjetdphimet",         &jetmax_dphimet_,       "maxjetdphimet/F");

  //trigger
  //babyTree_->Branch("HLT_Jet6U",             &HLT_Jet6U_,            "HLT_Jet6U/I");
  //babyTree_->Branch("HLT_Jet10U",            &HLT_Jet10U_,           "HLT_Jet10U/I");
  babyTree_->Branch("HLT_Jet15U",            &HLT_Jet15U_,           "HLT_Jet15U/I");
  babyTree_->Branch("HLT_Jet30U",            &HLT_Jet30U_,           "HLT_Jet30U/I");

  babyTree_->Branch("HLT_Photon10_L1R",      &HLT_Photon10_L1R_,     "HLT_Photon10_L1R/I");
  babyTree_->Branch("HLT_Photon15_L1R",      &HLT_Photon15_L1R_,     "HLT_Photon15_L1R/I");

  babyTree_->Branch("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_,     "HLT_Photon10_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_,     "HLT_Photon15_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_,     "HLT_Photon20_Cleaned_L1R/I");  

  //Z stuff
  babyTree_->Branch("passz",              &passz_,              "passz/I");  
  babyTree_->Branch("pdgid",              &pdgid_,              "pdgid/I");  
  babyTree_->Branch("passm_ll_nom",          &passm_ll_nom_,          "passm_ll_nom/I");  
  babyTree_->Branch("passm_ll_nomttbar",     &passm_ll_nomttbar_,     "passm_ll_nomttbar/I");  
  babyTree_->Branch("passm_ll_nomttbarV2",   &passm_ll_nomttbarV2_,   "passm_ll_nomttbarV2/I");  
  babyTree_->Branch("passe_ll_ttbar",        &passe_ll_ttbar_,        "passe_ll_ttbar/I");  
  babyTree_->Branch("passe_ll_ttbarV1",      &passe_ll_ttbarV1_,      "passe_ll_ttbarV1/I");  
  babyTree_->Branch("passe_ll_cand01",       &passe_ll_cand01_,       "passe_ll_cand01/I");  
  babyTree_->Branch("passm_lt_nom",          &passm_lt_nom_,          "passm_lt_nom/I");  
  babyTree_->Branch("passm_lt_nomttbar",     &passm_lt_nomttbar_,     "passm_lt_nomttbar/I");  
  babyTree_->Branch("passm_lt_nomttbarV2",   &passm_lt_nomttbarV2_,   "passm_lt_nomttbarV2/I");  
  babyTree_->Branch("passe_lt_ttbar",        &passe_lt_ttbar_,        "passe_lt_ttbar/I");  
  babyTree_->Branch("passe_lt_ttbarV1",      &passe_lt_ttbarV1_,      "passe_lt_ttbarV1/I");  
  babyTree_->Branch("passe_lt_cand01",       &passe_lt_cand01_,       "passe_lt_cand01/I");  
  babyTree_->Branch("ptll",               &ptll_,               "ptll/F");  
  babyTree_->Branch("ptlt",               &ptlt_,               "ptlt/F");  
  babyTree_->Branch("etall",              &etall_,              "etall/F");  
  babyTree_->Branch("etalt",              &etalt_,              "etalt/F");  
  babyTree_->Branch("dilmass",            &dilmass_,            "dilmass/F");  
  babyTree_->Branch("dilpt",              &dilpt_,              "dilpt/F");  
  babyTree_->Branch("flagll",             &flagll_,             "flagll/I");  
  babyTree_->Branch("flaglt",             &flaglt_,             "flaglt/I");  
  babyTree_->Branch("bptx"             ,  &bptx_            ,   "bptx/I"         );
  babyTree_->Branch("bsc"              ,  &bsc_             ,   "bsc/I"          );
  babyTree_->Branch("beamhalo"         ,  &beamhalo_        ,   "beamhalo/I"     );
  babyTree_->Branch("goodvtx"          ,  &goodvtx_         ,   "goodvtx/I"      );
  babyTree_->Branch("goodtrks"         ,  &goodtrks_        ,   "goodtrks/I"     );

}


void Z_looper::FillBabyNtuple ()
{
  babyTree_->Fill();
}

void Z_looper::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}


       
//ttbar electron ID
//if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbar , isData ))) continue;
//if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbar , isData ))) continue;

//nominal muon ID
//if (abs(hyp_ll_id()[hypIdx]) == 13  && (! (fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx]))))   continue;
//if (abs(hyp_lt_id()[hypIdx]) == 13  && (! (fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx]))))   continue;

//ttbar electron ID
//if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_cand01 ))) continue;
//if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_cand01 ))) continue;
