
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
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

char* iter          = "e_ttbarV1align";
//char* iter          = "e_ttbar";
bool usePV          = false;
bool makebaby       = false;
bool debug          = false;
bool calculateTCMET = true;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

using namespace tas;
void looper::ScanChain (TChain* chain, const char* prefix, bool isData, int nEvents){

  bookHistos();

  set_goodrun_file("Cert_TopAug13_Merged_135059-142664_goodruns.txt");
  ofile.open( Form( "output/%s_%s_events.txt" , prefix , iter ) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  MakeBabyNtuple( Form( "output/%s_%s_baby.root" , prefix , iter) );

  DorkyEventIdentifier dei;

  //pass fail counters
  int nPassDuplicate  = 0;
  int nPassGoodRun    = 0;

  float nGoodMu = 0;
  float nGoodEl = 0;

  if( debug ) cout << "Begin looping over files" << endl;

  if( isData ){
    cout << "|" << setw(8)  << "run"          << setw(4) 
         << "|" << setw(6)  << "lumi"         << setw(4) 
         << "|" << setw(12) << "event"        << setw(4) 
         << "|" << setw(6)  << "type"         << setw(4) 
         << "|" << setw(6)  << "njets"        << setw(4) 
         << "|" << setw(6)  << "nbtags"       << setw(4) 
         << "|" << setw(8)  << "tcmet"        << setw(4) 
         << "|" << setw(8)  << "cltcmet"      << setw(4) 
         << "|" << setw(8)  << "pftcmet"      << setw(4) 
         << "|" << setw(8)  << "pfmet"        << setw(4) 
         << "|" << setw(8)  << "dphi"         << setw(4) << "|" << endl; 
  }

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
	  

          //APPLY BASIC EVENT SELECTIONS
          //TRIGGER, TRACKING AND VERTEX
          if( dei.is_duplicate( DorkyEvent() ) ) continue;
          nPassDuplicate++;

          if (!isData || goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))
            nPassGoodRun++;
	  else continue;
        
          float weight = 1.;
          if     ( strcmp(prefix,"data") == 0  ) weight = 1;
          else if( strcmp(prefix,"ZJets") == 0 ) weight = 1.27 * 2.2121427 * 0.84e-3;
          else if( strcmp(prefix,"TTBar") == 0 ) weight = 0.1112306 * 0.84e-3;
          else if( strcmp(prefix,"Vqq") == 0   ) weight = 0.0408562 * 0.84e-3;
          else{
            cout << "I DON'T RECOGNIZE " << prefix << endl;
            exit(0);
          }

          for( unsigned int hypIdx = 0 ; hypIdx < hyp_type().size() ; hypIdx++ ){
      
            InitBabyNtuple();
            
            int myType = 99;
            if (hyp_type()[hypIdx] == 3) myType = 0;                          // ee
            if (hyp_type()[hypIdx] == 0) myType = 1;                          // mm
            if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) myType=2; // em
            if (myType == 99) {
              cout << "Skipping unknown dilepton type = " << hyp_type()[hypIdx] << endl;
              continue;
            }

            //----------SELECTION----------------------------------------------------

            //pT > (20,20) GeV
            if( hyp_ll_p4()[hypIdx].pt() < 20 ) continue;
            if( hyp_lt_p4()[hypIdx].pt() < 20 ) continue;
 
            //ttbar muon ID
            if (abs(hyp_ll_id()[hypIdx]) == 13  && (! (fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx],NominalTTbar))))   continue;
            if (abs(hyp_lt_id()[hypIdx]) == 13  && (! (fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx],NominalTTbar))))   continue;
         
            //ttbarV1 electron ID
            if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV1 , isData ))) continue;
            if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV1 , isData ))) continue;

            //ttbar electron ID
            //if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbar , isData ))) continue;
            //if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbar , isData ))) continue;

            //nominal muon ID
            //if (abs(hyp_ll_id()[hypIdx]) == 13  && (! (fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx]))))   continue;
            //if (abs(hyp_lt_id()[hypIdx]) == 13  && (! (fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx]))))   continue;
         
            //ttbar electron ID
            //if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_cand01 ))) continue;
            //if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_cand01 ))) continue;

            //same flavor
            if( hyp_type()[hypIdx] == 1 ) continue;
            if( hyp_type()[hypIdx] == 2 ) continue;
            
            //opposite sign
            if( hyp_lt_id()[hypIdx] *  hyp_ll_id()[hypIdx] > 0 ) continue;

            //fill dilmass histo before cutting on dilmass
            dilMass_ = hyp_p4()[hypIdx].mass();
            fillHistos( hdilMass          , dilMass_  , weight , myType );

            if( dilMass_ < 76. || dilMass_ > 106. ) continue;

            if( myType == 0 ) nGoodEl+=weight;
            if( myType == 1 ) nGoodMu+=weight;

            //----------CORRECT TCMET FOR HYP MUONS---------------------------------------------

            metStruct tcmetStruct = correctTCMETforHypMuons( hypIdx , 
                                                             evt_tcmet() * cos( evt_tcmetPhi() ), 
                                                             evt_tcmet() * sin( evt_tcmetPhi() ),
                                                             evt_tcsumet() );
              
 

            // out-of-the-box  tcmet stuff (corrected for hyp muons)
            tcmet_    = tcmetStruct.met;
            tcmetphi_ = tcmetStruct.metphi;
            tcsumet_  = tcmetStruct.sumet;

            if( isData ){
              
              metStruct tcmetStructNewCalo = correctTCMETforHypMuons( hypIdx , 
                                                                      evtNew_tcmet() * cos( evtNew_tcmetPhi() ), 
                                                                      evtNew_tcmet() * sin( evtNew_tcmetPhi() ),
                                                                      evtNew_tcsumet() );
              tcmetNew_calo_ = tcmetStructNewCalo.met;

              metStruct tcmetStructNewPFC = correctTCMETforHypMuons( hypIdx , 
                                                                     evtNewPFC_tcmet() * cos( evtNewPFC_tcmetPhi() ), 
                                                                     evtNewPFC_tcmet() * sin( evtNewPFC_tcmetPhi() ),
                                                                     evtNewPFC_tcsumet() );
              tcmetNew_pfc_ = tcmetStructNewPFC.met;

            }else{
              tcmetNew_calo_ = tcmet_;
              tcmetNew_pfc_  = tcmet_;
            }


            //---------JET STUFF--------------------------------------------------------------

            int nJets        = 0;
            float sumJetPt   = 0;
            int nBTags       = 0;
            float maxcosdphi = -99;
            int imax         = -1;
            
            for (unsigned int ijet = 0; ijet < pfjets_p4().size(); ijet++) {

              LorentzVector vjet = pfjets_p4().at(ijet);
              LorentzVector vlt  = hyp_lt_p4()[hypIdx];
              LorentzVector vll  = hyp_ll_p4()[hypIdx];
              if (dRbetweenVectors(vjet, vll) < 0.4) continue;
              if (dRbetweenVectors(vjet, vlt) < 0.4) continue;
              
              if ( vjet.pt() > 30. && fabs(vjet.eta()) < 2.5 ) {
                nJets++;
                sumJetPt+=vjet.pt();
                fillHistos( hjetpt , vjet.pt() , weight , myType );

                float dRmin = 100;
                int imin = -1;
                
                if( fabs( cos( tcmetphi_ - vjet.phi() ) ) > maxcosdphi ){
                  maxcosdphi = fabs( cos( tcmetphi_ - vjet.phi() ) );
                  imax = ijet;
                  dphijetmet_ = fabs( tcmetphi_ - vjet.phi() );
                  if( dphijetmet_ > TMath::Pi() ) dphijetmet_ = TMath::TwoPi() - dphijetmet_;
                }

                //find closest calojet to use btagging info
                for( unsigned int icalojet = 0 ; icalojet < jets_p4().size() ; icalojet++ ){

                  LorentzVector vcalojet = jets_p4().at(icalojet);
                  if( vcalojet.pt() * jets_cor().at(icalojet) < 10 ) continue;
                  
                  float dR = dRbetweenVectors(vjet, vcalojet);
                  if( dR < dRmin ){
                    dRmin = dR;
                    imin = icalojet;
                  }
                }

                if( imin > -1 ){
                  if( jets_simpleSecondaryVertexHighEffBJetTag().at(imin) > 1.74 ) nBTags++;
                  //if( jets_trackCountingHighEffBJetTag().at(imin) > 1.7 ) nBTags++;
                }

              }
            }

          
            
      

            // event stuff
            run_    = cms2.evt_run();
            lumi_   = cms2.evt_lumiBlock();
            event_  = cms2.evt_event();

            // pf met stuff
            pfmet_    = cms2.evt_pfmet();
            pfmetphi_ = cms2.evt_pfmetPhi();
            pfsumet_  = cms2.evt_pfsumet();

            //calomet
            met_       = cms2.evt_met();
            metphi_    = cms2.evt_metPhi();
            sumet_     = cms2.evt_sumet();
     
            //electron eta
            if( myType == 0 ){
              if( tcmet_ < 20 ){
                fillUnderOverFlow( heleta_metlt20 , hyp_lt_p4()[hypIdx].eta() , weight );
                fillUnderOverFlow( heleta_metlt20 , hyp_ll_p4()[hypIdx].eta() , weight );
              }
              else if( tcmet_ > 30 ){
                fillUnderOverFlow( heleta_metgt30 , hyp_lt_p4()[hypIdx].eta() , weight );
                fillUnderOverFlow( heleta_metgt30 , hyp_ll_p4()[hypIdx].eta() , weight );
              }
            }

            //muon eta
            if( myType == 1 ){
              if( tcmet_ < 20 ){
                fillUnderOverFlow( hmueta_metlt20 , hyp_lt_p4()[hypIdx].eta() , weight );
                fillUnderOverFlow( hmueta_metlt20 , hyp_ll_p4()[hypIdx].eta() , weight );
              }
              else if( tcmet_ > 30 ){
                fillUnderOverFlow( hmueta_metgt30 , hyp_lt_p4()[hypIdx].eta() , weight );
                fillUnderOverFlow( hmueta_metgt30 , hyp_ll_p4()[hypIdx].eta() , weight );
              }
            }

            //look at track pT's near electron
            if( myType == 0 && tcmet_ > 30 ){
              for( unsigned int itrk = 0 ; itrk < trks_trk_p4().size() ; ++itrk ){
               
                if( isMuon( itrk ) )            continue;
                if( isElectron( itrk ) )        continue;
                if( !isGoodTrack( itrk ) )      continue;
             
                LorentzVector vtrk = trks_trk_p4().at(itrk);
                LorentzVector vlt  = hyp_lt_p4()[hypIdx];
                LorentzVector vll  = hyp_ll_p4()[hypIdx];

                float drll = dRbetweenVectors(vtrk, vll);
                float drlt = dRbetweenVectors(vtrk, vlt);

                if( drll < 0.3 || drlt < 0.3 ){
                  fillUnderOverFlow( htrkptnearel , trks_trk_p4().at(itrk).pt() , weight );
                  //printEvent();
                  //cout << "pt eta phi " << vtrk.pt() << " " << vtrk.eta() << " " << vtrk.phi() << endl;
                }
              }
            }

            string leptype[2]={"ee","mm"};

            if( isData && ( tcmet_ > 30. || pfmet_ > 30. ) ){

              if( calculateTCMET ){
                
                // calculate tcmet on-the-fly
                metStruct structMET = correctedTCMET( true, ofile );
                
              }

              cout << "|" << setw(8)  << evt_run()                   << setw(4) 
                   << "|" << setw(6)  << evt_lumiBlock()             << setw(4) 
                   << "|" << setw(12) << evt_event()                 << setw(4) 
                   << "|" << setw(6)  << leptype[myType]             << setw(4) 
                   << "|" << setw(6)  << nJets                       << setw(4) 
                   << "|" << setw(6)  << nBTags                      << setw(4) 
                   << "|" << setw(8)  << fround(tcmet_,1)            << setw(4) 
                   << "|" << setw(8)  << fround(tcmetNew_calo_,1)    << setw(4) 
                   << "|" << setw(8)  << fround(tcmetNew_pfc_,1)     << setw(4) 
                   << "|" << setw(8)  << fround(pfmet_,1)            << setw(4) 
                   << "|" << setw(8)  << fround(dphijetmet_,2)       << setw(4) << "|" << endl; 
            }

            if( imax > -1 ){

              if( tcmet_ < 20 )
                fillHistos( hdphijetmet_metlt20  , dphijetmet_    , weight , myType , nJets  );

              if( tcmet_ > 30 )
                fillHistos( hdphijetmet_metgt30  , dphijetmet_    , weight , myType , nJets  );

            }

            fillHistos( hnjets            , nJets            , weight , myType );
            fillHistos( htcmet            , tcmet_           , weight , myType , nJets );
            fillHistos( htcmetNew_calo    , tcmetNew_calo_   , weight , myType , nJets  );
            fillHistos( htcmetNew_pfc     , tcmetNew_pfc_    , weight , myType , nJets  );
            fillHistos( hpfmet            , pfmet_           , weight , myType , nJets  );
            
            njets_   = nJets;
            leptype_ = myType;

            eventTree_->Fill();
          }// end loop over hypotheses
        } // end loop over events
    } // end loop over files
  
  cout << "\n\n********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events that pass dup veto: " << nPassDuplicate 
       << " (" << 100*(double)nPassDuplicate/nEventsTotal << "% of total)" << endl;
  cout << "Total number of events that pass good run/lumi: " << nPassGoodRun 
       << " (" << 100*(double)nPassGoodRun/nEventsTotal << "% of total)" << endl;

  cout << nGoodEl << " ee events in Z mass window" << endl;
  cout << nGoodMu << " mm events in Z mass window" << endl;

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;

  CloseBabyNtuple();

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form( "output/%s_%s_histos.root" , prefix , iter ) );
  deleteHistos();

  
} // end ScanChain

void looper::printEvent(  ostream& ostr ){
  ostr << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl; 
}

void looper::InitBabyNtuple ()
{
  // event stuff
  run_    = -999999;
  lumi_   = -999999;
  event_  = -999999;

  // pfmet stuff
  pfmet_         = -999999.;
  pfmetphi_      = -999999.;
  pfsumet_       = -999999.;

  // calomet stuff
  met_           = -999999.;
  metphi_        = -999999.;
  sumet_         = -999999.;

  // out-of-the-box
  tcmet_         = -999999.;
  tcmetphi_      = -999999.;
  tcsumet_       = -999999.;

  // new tcmet stuff
  tcmetNew_calo_ = -999999.;
  tcmetNew_pfc_  = -999999.;

  // dilepton stuff
  dilMass_       = -999999.;

  //jet-met stuff
  dphijetmet_    = -999999.;

  leptype_       = -999999;
  njets_         = -999999;
}


float looper::deltaPhi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  heleta_metgt30           = new TH1F("heleta_metgt30","tcmet > 30 GeV",1000,-5,5);
  heleta_metlt20           = new TH1F("heleta_metlt20","tcmet < 20 GeV",1000,-5,5);
  hmueta_metgt30           = new TH1F("hmueta_metgt30","tcmet > 30 GeV",1000,-5,5);
  hmueta_metlt20           = new TH1F("hmueta_metlt20","tcmet < 20 GeV",1000,-5,5);
  heleta_metgt30->GetXaxis()->SetTitle("electron #eta");
  heleta_metlt20->GetXaxis()->SetTitle("electron #eta");
  hmueta_metgt30->GetXaxis()->SetTitle("muon #eta");
  hmueta_metlt20->GetXaxis()->SetTitle("muon #eta");

  htrkptnearel             = new TH1F("htrkptnearel","tcmet > 30, #DeltaR(trk,el) < 0.3",100,0,100);
  htrkptnearel->GetXaxis()->SetTitle("track p_{T} (GeV)");

  char* leptype[4]   = {"ee", "mm", "em", "all"};
  char* jetbin[4]    = {"0j", "1j", "geq2j", "allj"};

  char* leptype_title[4]   = {"ee", "#mu#mu", "e#mu", "all leptons"};
  char* jetbin_title[4]    = {"0 jets", "1 jet", "#geq 2 jets", "all jets"};

  for (int i = 0; i < 4; i++) {
  
    hnjets[i] = new TH1F(Form("hnjets_%s",leptype[i]),      leptype_title[i],   5,-0.5,4.5);
    hnjets[i]->GetXaxis()->SetTitle("nJets");

    hjetpt[i] = new TH1F(Form("hjetpt_%s",leptype[i]),      leptype_title[i],   100,0,200);
    hjetpt[i]->GetXaxis()->SetTitle("jet p_{T} (GeV)");

    hdilMass[i] = new TH1F(Form("hdilMass_%s",leptype[i]),  leptype_title[i],   150,0,300);
    hdilMass[i]->GetXaxis()->SetTitle("M(ll) (GeV)");

    for (int j = 0; j < 4; j++) {

      char* suffix       = Form("%s_%s",leptype[i],jetbin[j]);
      char* suffix_title = Form("%s %s",leptype_title[i],jetbin_title[j]);
    
      htcmet[i][j]             = new TH1F(Form("htcmet_%s",suffix),         suffix_title,   100,0,100);
      hpfmet[i][j]             = new TH1F(Form("hpfmet_%s",suffix),         suffix_title,   100,0,100);
      htcmetNew_calo[i][j]     = new TH1F(Form("htcmetNew_calo_%s",suffix), suffix_title,   100,0,100);
      htcmetNew_pfc[i][j]      = new TH1F(Form("htcmetNew_pfc_%s",suffix),  suffix_title,   100,0,100);
    
        
      hdphijetmet_metlt20[i][j]    = new TH1F(Form("hdphijetmet_metlt20_%s",suffix),
                                              Form("tcmet < 20 GeV %s",suffix_title), 100,0,3.15);
      hdphijetmet_metgt30[i][j]    = new TH1F(Form("hdphijetmet_metgt30_%s",suffix),
                                              Form("tcmet > 30 GeV %s",suffix_title), 100,0,3.15);
            
      htcmet[i][j]->               GetXaxis()->SetTitle("tcmet (GeV)");
      hpfmet[i][j]->               GetXaxis()->SetTitle("pfmet (GeV)");
      htcmetNew_calo[i][j]->       GetXaxis()->SetTitle("tcmet (GeV)");
      htcmetNew_pfc[i][j]->        GetXaxis()->SetTitle("tcmet (GeV)");
   
      hdphijetmet_metlt20[i][j]->  GetXaxis()->SetTitle("#Delta#phi(jet,tcmet)");
      hdphijetmet_metgt30[i][j]->  GetXaxis()->SetTitle("#Delta#phi(jet,tcmet)");
    }
  }
}

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
   eventTree_->Branch("hase"             , &hase_             , "hase/O" );
   
   eventTree_->Branch("pfmet"            , &pfmet_            , "pfmet/F"   );
   eventTree_->Branch("pfmetphi"         , &pfmetphi_         , "pfmetphi/F");
   eventTree_->Branch("pfsumet"          , &pfsumet_          , "pfsumet/F" );
   
   eventTree_->Branch("met"              , &met_              , "met/F"      );
   eventTree_->Branch("metphi"           , &metphi_           , "metphi/F"   );
   eventTree_->Branch("sumet"            , &sumet_            , "sumet/F"    );
   
   eventTree_->Branch("tcmet"            , &tcmet_            , "tcmet/F"      );
   eventTree_->Branch("tcmetphi"         , &tcmetphi_         , "tcmetphi/F"   );
   eventTree_->Branch("tcsumet"          , &tcsumet_          , "tcsumet/F"    );
   
   eventTree_->Branch("tcmetcalo"        , &tcmetNew_calo_    , "tcmetcalo/F"      );
   eventTree_->Branch("tcmetpfc"         , &tcmetNew_pfc_     , "tcmetpfc/F"      );
   eventTree_->Branch("dilmass"          , &dilMass_          , "dilMass/F"      );

   eventTree_->Branch("leptype"          , &leptype_          , "leptype/I"      );
   eventTree_->Branch("njets"            , &njets_            , "njets/I"      );

   
}

//--------------------------------------------------------------------

void looper::fillHistos(TH1F *h1[4],float value, float weight, int myType)
{

  fillUnderOverFlow(h1[myType], value, weight);      
  fillUnderOverFlow(h1[3],      value, weight);      
}

//--------------------------------------------------------------------

void looper::fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{

  if( nJetsIdx > 2 ) nJetsIdx = 2;

  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      
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

bool looper::isGoodTrack( int index ) {
  
     float corrected_d0 = trks_d0corr().at(index);
   
     if( trks_algo().at(index) < 8 ) {

       float d0cut = sqrt( pow(0.015,2) + pow(0.5/trks_trk_p4().at(index).pt(),2) );
       if( d0cut > 0.3 ) d0cut = 0.3;

       if( fabs( corrected_d0 ) > d0cut )                return false;
       //if( trks_nlayers().at(index) < nlayerscut_4567_ ) return false;
     }
     else {
       //if( trks_nlayers().at(index) < nlayerscut_89_ )                     return false;
       if( trks_validHits().at(index) < 9 )                                return false;
       if( trks_chi2().at(index) / trks_ndof().at(index) > 5. )            return false;
       if( trks_ptErr().at(index) / trks_trk_p4().at(index).pt() > 0.20 )  return false;
     }

     if( trks_validHits().at(index) < 6 )                               return false;
     if( trks_chi2().at(index) / trks_ndof().at(index) > 5 )            return false;
     if( fabs( trks_trk_p4().at(index).eta() ) > 2.65 )                 return false;
     if( trks_trk_p4().at(index).pt() > 100 )                           return false;
     if( trks_ptErr().at(index) / trks_trk_p4().at(index).pt() > 0.20 ) return false;
     if( !isTrackQuality( index, (1 << highPurity) ) )                  return false;

     if( trks_trk_p4().at(index).pt() > 0 && fabs(trks_trk_p4().at(index).eta()) > 2.5 ) return false;
     
     return true;
}

//--------------------------------------------------------------------

bool looper::isMuon( int index ) {

     for( unsigned int i = 0; i < mus_p4().size(); i++ ) {

	  if( mus_trkidx().at(i) == index ) return true;
     }

     return false;
}

//--------------------------------------------------------------------

bool looper::isElectron( int index ) {

     for( unsigned int i = 0; i < els_p4().size(); i++ ) {

	  if( els_trkidx().at(i) == index && els_hOverE().at(i) < 0.1 ) return true;
     }

     return false;
}

//--------------------------------------------------------------------
