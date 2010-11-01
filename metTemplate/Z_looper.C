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
#include "CORE/ttbarSelections.cc"
#include "CORE/susySelections.cc"
#include "CORE/jetSelections.cc"
//#include "CORE/triggerUtils.cc"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

using namespace tas;
//inline double fround(double n, double d){
//  return floor(n * pow(10., d) + .5) / pow(10., d);
//}


enum metType   { e_tcmet = 0, e_tcmetNew = 1, e_pfmet = 2};
enum templateSource { e_QCD = 0, e_PhotonJet = 1 };

//--------------------------------------------------------------------

bool debug = false;
const int nJetBins        = 3;
const int nSumJetPtBins   = 7;
const int nBosonPtBins    = 4;

float lumi         = 0.01106;
char* iter         = "oct15th_v2";
char* jsonfilename = "Cert_TopOct15_Merged_135821-147454_allPVT_V2_goodruns.txt";

metType myMetType               = e_tcmet;
templateSource myTemplateSource = e_PhotonJet;

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

int getBosonPtBin( float x ){

  if     ( x < 40.           ) return 0;
  else if( x > 40 && x < 60  ) return 1;
  else if( x > 60 && x < 100 ) return 2;
  else if( x > 100           ) return 3;
  else { 
    cout << "Error could not find boson pt bin" << endl;
    exit(0);
  }
}

//--------------------------------------------------------------------

string bosonPtString( int bin ){

  if     ( bin == 0 ) return "boson p_{T} < 40 GeV";
  else if( bin == 1 ) return "40 < boson p_{T} < 60 GeV";
  else if( bin == 2 ) return "60 < boson p_{T} < 100 GeV";
  else if( bin == 3 ) return "boson p_{T} > 100 GeV";
  else{
    cout << "Error unrecognized boson pt bin " << bin << endl;
    exit(0);
  }

}

//--------------------------------------------------------------------

int getJetBin( int njets ){

  int bin = -1;

  if     ( njets == 1 ) bin = 0;
  else if( njets == 2 ) bin = 1;
  else if( njets >= 3 ) bin = 2;
  else{
    cout << "Error " << njets << " jets passed as argument to met templates" << endl;
    exit(0);
  }

  return bin;
}

//--------------------------------------------------------------------

string jetString( int bin ){

  string js;

  if     ( bin == 0 ) js = "njets = 1";
  else if( bin == 1 ) js = "njets = 2";
  else if( bin == 2 ) js = "njets >= 3";
  else{
    cout << "Error invalid bin passed as argument to met templates" << endl;
    exit(0);
  }

  return js;

}

//--------------------------------------------------------------------

int getSumJetPtBin( float x ){

  //bins array defines the sumJetPt binning
  //float bins[nSumJetPtBins+1]={0,25,50,75,100,150,250,5000};
  float bins[nSumJetPtBins+1]={0,30,60,90,120,150,250,5000};
  
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

//--------------------------------------------------------------------

string sumJetPtString( int bin ){

  //float bins[nSumJetPtBins+1]={0,25,50,75,100,150,250,5000};
  float bins[nSumJetPtBins+1]={0,30,60,90,120,150,250,5000};

  stringstream s;
  s << bins[bin] << " < sumJetPt < " << bins[bin+1] << " GeV";

  return s.str();
}

//--------------------------------------------------------------------

void Z_looper::ScanChain (TChain* chain, const char* prefix, bool isData,
                          bool calculateTCMET, metAlgo algo, int nEvents, float kFactor){

  algo_ = algo;

  if     ( algo_ == e_makeTemplate )    cout << "metAlgo makeTemplate" << endl;
  else if( algo_ == e_photonSelection ) cout << "metAlgo photonSelection" << endl;
  else if( algo_ == e_ZSelection )      cout << "metAlgo ZSelection" << endl;
  else { cout << "Unrecognized algo" << endl; exit(0); }

  if     ( myMetType == e_tcmet    ) cout << "metType tcmet" << endl;
  else if( myMetType == e_pfmet    ) cout << "metType pfmet" << endl;
  else if( myMetType == e_tcmetNew ) cout << "metType tcmetNew" << endl;
  else { cout << "Unrecognized metType" << endl; exit(0); }

  if     ( myTemplateSource == e_QCD       ) cout << "QCD templates" << endl;
  else if( myTemplateSource == e_PhotonJet ) cout << "photon+jets templates" << endl;
  else { cout << "Unrecognized templateSource" << endl; exit(0); }

  TFile *metTemplateFile   = new TFile();
  string metTemplateString = "";
  char* templateFileName   = "";

  //select templates
  if( isData ){
    
    if( myTemplateSource == e_QCD ){
      cout << "QCD templates are deprecated. If you want to use this you"
           << "need to make QCD templates using makeTemplates.C." << endl;
      exit(0);
      //       templateFileName = Form("output/%s/JetMETTau_templates.root",iter);
      //       cout << "Using template file " << templateFileName << endl;
      //       metTemplateString = "_JetMETTauTemplate";
      //       metTemplateFile = TFile::Open( templateFileName );
    }
    
    else if( myTemplateSource == e_PhotonJet ){
      templateFileName = Form("output/%s/EG_templates.root",iter);
      cout << "Using template file " << templateFileName << endl;
      metTemplateString = "_EGTemplate";
      metTemplateFile = TFile::Open( templateFileName );
    }
    
  }else{
    
    if( myTemplateSource == e_QCD ){
      cout << "QCD templates are deprecated. If you want to use this you"
           << "need to make QCD templates using makeTemplates.C." << endl;
      exit(0);
      //       templateFileName = Form("output/%s/QCD_Pt15_templates.root",iter);
      //       cout << "Using template file " << templateFileName << endl;
      //       metTemplateString = "_QCD_Pt15Template";
      //       metTemplateFile = TFile::Open( templateFileName ); 
    }
        
    else if( myTemplateSource == e_PhotonJet ){
      templateFileName = Form("output/%s/PhotonJet_templates.root",iter);
      cout << "Using template file " << templateFileName << endl;
      metTemplateString = "_PhotonJetTemplate";
      metTemplateFile = TFile::Open( templateFileName );
    }
  }
  
  if( metTemplateFile == 0 ){
    cout << "Error couldn't find " << templateFileName << endl;
    exit(0);
  }

  set_goodrun_file( jsonfilename );
 
  if( isData ){
    
    ofile_tcmet.open(  Form( "output/%s/%s_tcmetprintout.txt" , iter , prefix  ) );
    ofile_events.open( Form( "output/%s/%s_highmetevents.txt" , iter , prefix  ) );
    
    ofile_events << "|" << setw(8)  << "run"          << setw(4) 
                 << "|" << setw(6)  << "lumi"         << setw(4) 
                 << "|" << setw(12) << "event"        << setw(4) 
                 << "|" << setw(6)  << "type"         << setw(4) 
                 << "|" << setw(6)  << "njets"        << setw(4) 
                 << "|" << setw(6)  << "nbtags"       << setw(4) 
                 << "|" << setw(8)  << "tcmet"        << setw(4) 
                 << "|" << setw(8)  << "pfmet"        << setw(4) 
                 << "|" << setw(8)  << "dphi"         << setw(4) << "|" << endl; 
  }


  bookHistos();

  // make a baby ntuple
  stringstream babyfilename;
  babyfilename << prefix << "_baby.root";
  MakeBabyNtuple( Form("output/%s/%s_baby.root", iter , prefix ) );

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  //pass fail counters
  float nGoodMu = 0;
  float nGoodEl = 0;
  float nGoodEM = 0;
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

      if( isData ) {
	DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
	if (is_duplicate(id) )
	  continue;
      }
     
      //skip events with bad els_conv_dist 
      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
      }
      
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }

      //good run+event selection-----------------------------------------------------------
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      if( !cleaning_standardAugust2010( isData) )                    continue;
      
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

      // calomet, pfmet, genmet
      met_       = cms2.evt_met();
      metphi_    = cms2.evt_metPhi();
      sumet_     = cms2.evt_sumet();

      pfmet_    = cms2.evt_pfmet();
      pfmetphi_ = cms2.evt_pfmetPhi();
      pfsumet_  = cms2.evt_pfsumet();

      if (!isData){
        genmet_     = cms2.gen_met();
        genmetphi_  = cms2.gen_metPhi();
        gensumet_   = cms2.gen_sumEt();
      }
      
      vector<unsigned int> v_goodHyps;
      v_goodHyps.clear();

      int nHypPass = 0;

      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); ++hypIdx) {

        if( !passSUSYTrigger_v1( isData , hyp_type()[hypIdx] ) ) continue;

        //check that hyp leptons come from same vertex
        if(!hypsFromSameVtx(hypIdx))   continue;
        
        //selection--------------------------------------------------------------------- 

        //OS, pt > (20,20) GeV
        if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                             continue;
        //if( hyp_ll_p4()[hypIdx].pt() < 20. )                                          continue;
        //if( hyp_lt_p4()[hypIdx].pt() < 20. )                                          continue;
        if( TMath::Max( hyp_ll_p4()[hypIdx].pt() , hyp_lt_p4()[hypIdx].pt() ) < 20. )   continue;
        if( TMath::Min( hyp_ll_p4()[hypIdx].pt() , hyp_lt_p4()[hypIdx].pt() ) < 20. )   continue;
        if( hyp_p4()[hypIdx].mass() < 10 )                                              continue;

        //ttbarV2 muon ID
        //if (abs(hyp_ll_id()[hypIdx]) == 13  && (! muonId(hyp_ll_index()[hypIdx],NominalTTbarV2)))   continue;
        //if (abs(hyp_lt_id()[hypIdx]) == 13  && (! muonId(hyp_lt_index()[hypIdx],NominalTTbarV2)))   continue;
 
        //nominal muon ID
        if (abs(hyp_ll_id()[hypIdx]) == 13  && !( fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx],Nominal)))   continue;
        if (abs(hyp_lt_id()[hypIdx]) == 13  && !( fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx],Nominal)))   continue;
        
        //OSV1 electron ID
        if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_el_OSV1 , false , false ))) continue;
        if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_el_OSV1 , false , false ))) continue;
        
        //Z-mass constraint
     
        if( hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  continue;

        nHypPass++;
      
        v_goodHyps.push_back( hypIdx );
      
      }

      if( v_goodHyps.size() == 0 ) continue;

      unsigned int hypIdx = selectBestZHyp(v_goodHyps);

      leptype_ = 99;
      if (hyp_type()[hypIdx] == 3) leptype_ = 0;                           // ee
      if (hyp_type()[hypIdx] == 0) leptype_ = 1;                           // mm
      if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) leptype_=2;  // em
      if (leptype_ == 99) {
        cout << "Skipping unknown dilepton type = " << hyp_type()[hypIdx] << endl;
        continue;
      }

      dilmass_ = hyp_p4()[hypIdx].mass();
      fillHistos( hdilMass          , dilmass_  , weight_ , leptype_ );

      if( leptype_ == 0 ) nGoodEl+=weight_;
      if( leptype_ == 1 ) nGoodMu+=weight_;
      if( leptype_ == 2 ) nGoodEM+=weight_;
  
      //dilepton stuff---------------------------------------------------------------- 
        
      idll_             = hyp_ll_id()[hypIdx];   
      idlt_             = hyp_lt_id()[hypIdx];        
      ptll_             = hyp_ll_p4()[hypIdx].pt();
      ptlt_             = hyp_lt_p4()[hypIdx].pt();
      etall_            = hyp_ll_p4()[hypIdx].eta();
      etalt_            = hyp_lt_p4()[hypIdx].eta();
      phill_            = hyp_ll_p4()[hypIdx].phi();
      philt_            = hyp_lt_p4()[hypIdx].phi();
      dilmass_          = hyp_p4()[hypIdx].mass(); 
      dilpt_            = hyp_p4()[hypIdx].pt(); 

      if( abs( hyp_ll_id()[hypIdx] ) == 11 ){
        passe_ll_ttbarV1_      = pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV1 , isData , true ) ? 1 : 0;
        passe_ll_ttbarV2_      = pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV2 , isData , true ) ? 1 : 0;
        passe_ll_ttbar_        = pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbar   , isData , true ) ? 1 : 0;
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
        passe_lt_ttbarV1_      = pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV1 , isData , true ) ? 1 : 0;
        passe_lt_ttbarV2_      = pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV2 , isData , true ) ? 1 : 0;
        passe_lt_ttbar_        = pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbar   , isData , true ) ? 1 : 0;
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

      //tcmet stuff------------------------------------------------------------------- 

      metStruct tcmetStruct = correctTCMETforHypMuons( hypIdx ,  
                                                       evt_tcmet() * cos( evt_tcmetPhi() ), 
                                                       evt_tcmet() * sin( evt_tcmetPhi() ), 
                                                       evt_tcsumet() );
        
      // out-of-the-box  tcmet stuff (corrected for hyp muons)
      tcmet_    = tcmetStruct.met;
      tcmetphi_ = tcmetStruct.metphi;
      tcsumet_  = tcmetStruct.sumet;   

      if( calculateTCMET ){
        
          
        metStruct tcmetNewStruct = correctedTCMET();
        tcmetNew_     = tcmetNewStruct.met;
        tcmetphiNew_  = tcmetNewStruct.metphi;
        tcsumetNew_   = tcmetNewStruct.sumet;
          

        metStruct tcmetNewStruct_corr = correctTCMETforHypMuons( hypIdx ,  
                                                                 tcmetNew_ * cos( tcmetphiNew_ ), 
                                                                 tcmetNew_ * sin( tcmetphiNew_ ), 
                                                                 tcsumetNew_ );
          
        // Z_looper-level  tcmet stuff (corrected for hyp muons)
        tcmetNew_    = tcmetNewStruct_corr.met;
        tcmetphiNew_ = tcmetNewStruct_corr.metphi;
        tcsumetNew_  = tcmetNewStruct_corr.sumet;   

      }else{
        
        tcmetNew_    = -9999;
        tcmetphiNew_ = -9999;
        tcsumetNew_  = -9999;
       
      }

      nGoodVertex_ = 0;
      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
        if(isGoodVertex(v)) nGoodVertex_++;
      }
         
      //jet stuff--------------------------------------------------------------------- 
        
      nJets_        = 0;
      sumJetPt_     = 0.;
      nJets40_      = 0;
      sumJetPt10_   = 0.;
      nbtags_       = 0;

      LorentzVector jetSystem(0.,0.,0.,0.);        
      float maxcosdphi  = -99;
      //int   imaxcosdphi = -1;
      int   imaxjet     = -1;
      float maxpt       = -1;

      failjetid_ = 0;
        
      //loop over pfjets pt > 30 GeV |eta| < 2.5
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
        LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
        LorentzVector vlt  = hyp_lt_p4()[hypIdx];
        LorentzVector vll  = hyp_ll_p4()[hypIdx];
          
        if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
        if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
        if( fabs( vjet.eta() ) > 2.5 )           continue;
        if( !passesPFJetID(ijet) ){
          failjetid_ = 1;
          continue;
        }

        if ( vjet.pt() > 10. ){
          sumJetPt10_ += vjet.pt();
        }
        if ( vjet.pt() > 15. ){
          sumJetPt_ += vjet.pt();
          jetSystem += vjet;
        }

        if( vjet.pt() < 30. )                    continue;
          
        //find max jet pt
        if( vjet.pt() > maxpt ){
          maxpt   = vjet.pt();
          imaxjet = ijet;
        }
          
        //find jet (anti-)aligned with tcmet
        if( fabs( cos( tcmetphi_ - vjet.phi() ) ) > maxcosdphi ){
          maxcosdphi  = fabs( cos( tcmetphi_ - vjet.phi() ) );
          //imaxcosphi  = ijet;
          dphijetmet_ = fabs( tcmetphi_ - vjet.phi() );
          if( dphijetmet_ > TMath::Pi() ) dphijetmet_ = TMath::TwoPi() - dphijetmet_;
        }

        //find closest calojet to use btagging info
        float dRmin    = 100;
        int   iCaloJet = -1;
          
        for( unsigned int iC = 0 ; iC < jets_p4().size() ; iC++ ){
            
          LorentzVector vcalojet = jets_p4().at(iC);
          if( vcalojet.pt() * jets_cor().at(iC) < 10 ) continue;
            
          float dR = dRbetweenVectors(vjet, vcalojet);
          if( dR < dRmin ){
            dRmin = dR;
            iCaloJet = iC;
          }
        }
          
        if( iCaloJet > -1 ){
          if( jets_simpleSecondaryVertexHighEffBJetTag().at(iCaloJet) > 1.74 ) ++nbtags_;
          //if( jets_trackCountingHighEffBJetTag().at(iCaloJet) > 1.7 ) ++nbtags_;
        }

        if ( vjet.pt() > 30. ) nJets_++;
        if ( vjet.pt() > 40. ) nJets40_++;
          
      }
       
      jetmax_pt_ = -1;

      if( imaxjet > -1 ){
        jetmax_pt_       = pfjets_cor().at(imaxjet) * pfjets_p4().at(imaxjet).pt();
        jetmax_dphimet_  = deltaPhi( pfjets_p4().at(imaxjet).phi() , tcmetphi_);
      }

      vecJetPt_ = jetSystem.pt();
      
      //fill histos and ntuple----------------------------------------------------------- 

      FillBabyNtuple();
      fillHistos( htcmet            , tcmet_           , weight_ , leptype_ , nJets_ );
      fillHistos( htcmetNew         , tcmetNew_        , weight_ , leptype_ , nJets_ );
      fillHistos( hpfmet            , pfmet_           , weight_ , leptype_ , nJets_  );

      if( isData && ( tcmet_ > 30 || pfmet_ > 30 ) ){
        metStruct dummyStruct = correctedTCMET( true, ofile_tcmet );
          
        string lepstring[3]={"ee","mm","em"};

        ofile_events << "|" << setw(8)  << evt_run()                   << setw(4) 
                     << "|" << setw(6)  << evt_lumiBlock()             << setw(4) 
                     << "|" << setw(12) << evt_event()                 << setw(4) 
                     << "|" << setw(6)  << lepstring[leptype_]         << setw(4) 
                     << "|" << setw(6)  << nJets_                      << setw(4) 
                     << "|" << setw(6)  << nbtags_                     << setw(4) 
                     << "|" << setw(8)  << fround(tcmet_,1)            << setw(4) 
                     << "|" << setw(8)  << fround(pfmet_,1)            << setw(4) 
                     << "|" << setw(8)  << fround(dphijetmet_,2)       << setw(4) << "|" << endl; 
          
      }

      //met templates-------------------------------------------------------------------- 
 
      //if ( nJets_ < 2 )                                                     continue;
      //if( jetmax_pt_ < 50 )                                                 continue;
      if ( nJets_ < 1 )                                                     continue;
      if( jetmax_pt_ < 30 )                                                 continue;

      int iJ = nJets_;
      if( iJ > 4 ) iJ = 4;

      fillUnderOverFlow( hptz[iJ] , dilpt_ , weight_ );
      fillUnderOverFlow( hptz[0]  , dilpt_ , weight_ );
 
      float theMet    = -1;
      float theMetPhi = -1;

      if( myMetType == e_tcmet ){
        theMet    = tcmet_;
        theMetPhi = tcmetphi_;
      }
      else if( myMetType == e_tcmetNew ){
        theMet    = tcmetNew_;
        theMetPhi = tcmetphiNew_;
      }
      else if( myMetType == e_pfmet ){
        theMet    = pfmet_;
        theMetPhi = pfmetphi_;
      }

      dphixmet_  = deltaPhi( theMetPhi , hyp_p4()[hypIdx].phi() );
      metPar_    = theMet * cos( dphixmet_ );
      metPerp_   = theMet * sin( dphixmet_ );
        
      //fill predicted and observed met histos
      int iJetBin      = getJetBin( nJets_ );
      int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
      int iBosonPtBin  = getBosonPtBin( dilpt_ );

      TH1F* hmet = new TH1F();
        
      if( !isData && nJets_ > 1 ){
        if( myMetType == e_tcmet ){
          hmet     = (TH1F*) metTemplateFile->Get(Form("tcmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin));
        }
        else if( myMetType == e_tcmetNew ){
          hmet     = (TH1F*) metTemplateFile->Get(Form("tcmetNewTemplate_combined_%i_%i",iJetBin,iSumJetPtBin));
        }
        else if( myMetType == e_pfmet ){
          hmet     = (TH1F*) metTemplateFile->Get(Form("pfmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin));
        }
      }
        
      else{
        if     ( myMetType == e_tcmet ){
          hmet     = (TH1F*) metTemplateFile->Get(Form("tcmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin));
        }
        else if( myMetType == e_tcmetNew ){
          hmet     = (TH1F*) metTemplateFile->Get(Form("tcmetNewTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin));
        }
        else if( myMetType == e_pfmet ){
          hmet     = (TH1F*) metTemplateFile->Get(Form("pfmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin));
        }
      }

      hmet->Scale( weight_ );

      iJ = iJetBin;
      if( iJ > 3 ) iJ = 3;
      fillUnderOverFlow( metObserved_njets[iJ]  ,  theMet , weight_ );
      metPredicted_njets[iJ]->Add( hmet );
        
      fillUnderOverFlow( metObserved , theMet , weight_  );
      metPredicted->Add( hmet );

      // SF vs. DF
      if( leptype_ == 0 || leptype_ == 1 ){
        fillUnderOverFlow( metObserved_sf , theMet , weight_  );
        metPredicted_sf->Add( hmet );
      }
      else if( leptype_ == 2 ){
        fillUnderOverFlow( metObserved_df , theMet , weight_  );
        metPredicted_df->Add( hmet );
      }

      // ee vs. mumu
      if( leptype_ == 0 ){
        fillUnderOverFlow( metObserved_ee , theMet , weight_  );
        metPredicted_ee->Add( hmet );
      }
      else if( leptype_ == 1 ){
        fillUnderOverFlow( metObserved_mm , theMet , weight_  );
        metPredicted_mm->Add( hmet );
      }
        
      delete hmet;

      //}// end loop over hypIdx
 
      if( nHypPass > 1 && isData ) 
        cout << "Found " << nHypPass << " hypotheses passing selection" << endl;
    } // end loop over events
  } // end loop over files

  if( nSkip_els_conv_dist > 0 ){
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist branch" << endl;
  }

  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  cout << nGoodEl << " ee events in Z mass window" << endl;
  cout << nGoodMu << " mm events in Z mass window" << endl;
  cout << nGoodEM << " em events in Z mass window" << endl;

  CloseBabyNtuple();

  // make histos rootfile
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form("output/%s/%s%s.root", iter , prefix , metTemplateString.c_str() ) );
  deleteHistos();
  
} // end ScanChain

//--------------------------------------------------------------------

float Z_looper::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//--------------------------------------------------------------------

void Z_looper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

//--------------------------------------------------------------------

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
  HLT_Photon20_L1R_     = -1;

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

  nbtags_       = -999999;
  dphijetmet_   = -999999;

  //leading jet stuff
  jetmax_pt_        = -999999;
  jetmax_dphimet_   = -999999;

  //Z stuff
  passz_           = -999999;
  passe_ll_ttbar_     = -999999;
  passe_ll_ttbarV1_   = -999999;
  passe_ll_ttbarV2_   = -999999;
  passe_ll_cand01_    = -999999;
  passm_ll_nomttbar_  = -999999;
  passm_ll_nomttbarV2_= -999999;
  passm_ll_nom_       = -999999;
  passe_lt_ttbar_     = -999999;
  passe_lt_ttbarV1_   = -999999;
  passe_lt_ttbarV2_   = -999999;
  passe_lt_cand01_    = -999999;
  passm_lt_nomttbar_  = -999999;
  passm_lt_nomttbarV2_= -999999;
  passm_lt_nom_       = -999999;
  pdgid_           = -999999;
  ptll_            = -999999;
  ptlt_            = -999999;
  idll_            = -999999;
  idlt_            = -999999;
  etall_           = -999999;
  etalt_           = -999999;
  phill_           = -999999;
  philt_           = -999999;
  dilmass_         = -999999.;
  dilpt_           = -999999.;
  flagll_          = -999999;
  flaglt_          = -999999;
  failjetid_       = -999999;

  bptx_        =  -999999;
  bsc_         =  -999999;
  beamhalo_    =  -999999;
  goodvtx_     =  -999999;
  goodtrks_    =  -999999;

}

//--------------------------------------------------------------------

void Z_looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  char* pttitle[5]={"all jets","1 jet","2 jet","3 jet","#geq 4 jet"};

  for( int iJ = 0 ; iJ < 5 ; iJ++ ){
    hptz[iJ] = new TH1F(Form("hptz_%i",iJ),pttitle[iJ],200,0,200);
    hptz[iJ]->GetXaxis()->SetTitle("Z p_{T} (GeV)");
  }

  Double_t maxmet = 200;

  metObserved     = new TH1F("metObserved", "Observed MET",maxmet,0,maxmet);
  metPredicted    = new TH1F("metPredicted","Predicted MET",maxmet,0,maxmet);
  metObserved->Sumw2();
  metPredicted->Sumw2();

  metObserved_sf  = new TH1F("metObserved_sf", "Observed MET (SF)",maxmet,0,maxmet);
  metPredicted_sf = new TH1F("metPredicted_sf","Predicted MET (SF)",maxmet,0,maxmet);
  metObserved_sf->Sumw2();
  metPredicted_sf->Sumw2();

  metObserved_df  = new TH1F("metObserved_df", "Observed MET (DF)",maxmet,0,maxmet);
  metPredicted_df = new TH1F("metPredicted_df","Predicted MET (DF)",maxmet,0,maxmet);
  metObserved_df->Sumw2();
  metPredicted_df->Sumw2();

  metObserved_ee  = new TH1F("metObserved_ee", "Observed MET (ee)",maxmet,0,maxmet);
  metPredicted_ee = new TH1F("metPredicted_ee","Predicted MET (ee)",maxmet,0,maxmet);
  metObserved_ee->Sumw2();
  metPredicted_ee->Sumw2();

  metObserved_mm  = new TH1F("metObserved_mm", "Observed MET (#mu#mu)",maxmet,0,maxmet);
  metPredicted_mm = new TH1F("metPredicted_mm","Predicted MET (#mu#mu)",maxmet,0,maxmet);
  metObserved_mm->Sumw2();
  metPredicted_mm->Sumw2();

  metParObserved  = new TH1F("metParObserved", "Observed MET (Parallel)",1000,-maxmet,maxmet);
  metParPredicted = new TH1F("metParPredicted","Predicted MET (Parallel)",1000,-maxmet,maxmet);
  metParObserved->Sumw2();
  metParPredicted->Sumw2();

  metPerpObserved  = new TH1F("metPerpObserved", "Observed MET (Perpendicular)",maxmet,0,maxmet);
  metPerpPredicted = new TH1F("metPerpPredicted","Predicted MET (Perpendicular)",maxmet,0,maxmet);
  metPerpObserved->Sumw2();
  metPerpPredicted->Sumw2();

  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){

    metObserved_njets[iJetBin]  = new TH1F(Form("metObserved_njets%i",iJetBin), Form("Observed MET NJets %i", iJetBin),maxmet,0,maxmet);
    metPredicted_njets[iJetBin] = new TH1F(Form("metPredicted_njets%i",iJetBin),Form("Predicted MET NJets %i",iJetBin),maxmet,0,maxmet);
    
    metObserved_njets[iJetBin] ->Sumw2();
    metPredicted_njets[iJetBin]->Sumw2();
  }
  
  char* leptype[4]   = {"ee", "mm", "em", "all"};
  char* jetbin[4]    = {"0j", "1j", "geq2j", "allj"};

  char* leptype_title[4]   = {"ee", "#mu#mu", "e#mu", "all leptons"};
  char* jetbin_title[4]    = {"0 jets", "1 jet", "#geq 2 jets", "all jets"};

  for (int i = 0; i < 4; i++) {
   
    hdilMass[i] = new TH1F(Form("hdilMass_%s",leptype[i]),  leptype_title[i],   150,0,300);
    hdilMass[i]->GetXaxis()->SetTitle("M(ll) (GeV)");
 
    for (int j = 0; j < 4; j++) {

      char* suffix       = Form("%s_%s",leptype[i],jetbin[j]);
      char* suffix_title = Form("%s %s",leptype_title[i],jetbin_title[j]);
    
      htcmet[i][j]    = new TH1F(Form("htcmet_%s",suffix),    suffix_title, 100,0,100);
      htcmetNew[i][j] = new TH1F(Form("htcmetNew_%s",suffix), suffix_title, 100,0,100);
      hpfmet[i][j]    = new TH1F(Form("hpfmet_%s",suffix),    suffix_title, 100,0,100);
      htcmet[i][j]->GetXaxis()->SetTitle("tcmet (GeV)");
      htcmetNew[i][j]->GetXaxis()->SetTitle("tcmetNew (GeV)");
      hpfmet[i][j]->GetXaxis()->SetTitle("pfmet (GeV)");
    }
  }
  
}
 
//--------------------------------------------------------------------

void Z_looper::MakeBabyNtuple (const char* babyFileName)
{

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("T1", "A Baby Ntuple");

  //event stuff
  babyTree_->Branch("run",          &run_,          "run/I"  );
  babyTree_->Branch("lumi",         &lumi_,         "lumi/I" );
  babyTree_->Branch("event",        &event_,        "event/I");
  babyTree_->Branch("failjetid",    &failjetid_,    "failjetid/I");
  babyTree_->Branch("nvtx",         &nGoodVertex_,  "nvtx/I");
  babyTree_->Branch("weight",       &weight_,       "weight/F");
  babyTree_->Branch("pthat",        &pthat_,        "pthat/F");

  //met stuff
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

  //jet stuff
  babyTree_->Branch("njets",          &nJets_,            "njets/I"       );
  babyTree_->Branch("njets40",        &nJets40_,          "njets40/I"     );
  babyTree_->Branch("sumjetpt",       &sumJetPt_,         "sumjetpt/F"    );
  babyTree_->Branch("sumjetpt10",     &sumJetPt10_,       "sumjetpt10/F"    );
  babyTree_->Branch("vecjetpt",       &vecJetPt_,         "vecjetpt/F"    );
  babyTree_->Branch("nbtags",         &nbtags_,           "nbtags/I");
  babyTree_->Branch("ndphijetmet",    &dphijetmet_,       "dphijetmet/F");
  babyTree_->Branch("maxjetpt",       &jetmax_pt_,        "maxjetpt/F");
  babyTree_->Branch("maxjetdphimet",  &jetmax_dphimet_,   "maxjetdphimet/F");
                                                    
  //trigger stuff
  babyTree_->Branch("HLT_Jet15U",                    &HLT_Jet15U_,                   "HLT_Jet15U/I");
  babyTree_->Branch("HLT_Jet30U",                    &HLT_Jet30U_,                   "HLT_Jet30U/I");
  babyTree_->Branch("HLT_Photon10_L1R",              &HLT_Photon10_L1R_,             "HLT_Photon10_L1R/I");
  babyTree_->Branch("HLT_Photon15_L1R",              &HLT_Photon15_L1R_,             "HLT_Photon15_L1R/I");
  babyTree_->Branch("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_,     "HLT_Photon10_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_,     "HLT_Photon15_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_,     "HLT_Photon20_Cleaned_L1R/I");  
  babyTree_->Branch("HLT_Photon20_L1R",              &HLT_Photon20_L1R_,             "HLT_Photon20__L1R/I");  

  //Z stuff
  babyTree_->Branch("leptype",               &leptype_,               "leptype/I");
  babyTree_->Branch("passz",                 &passz_,                 "passz/I");  
  babyTree_->Branch("pdgid",                 &pdgid_,                 "pdgid/I");  
  babyTree_->Branch("passm_ll_nom",          &passm_ll_nom_,          "passm_ll_nom/I");  
  babyTree_->Branch("passm_ll_nomttbar",     &passm_ll_nomttbar_,     "passm_ll_nomttbar/I");  
  babyTree_->Branch("passm_ll_nomttbarV2",   &passm_ll_nomttbarV2_,   "passm_ll_nomttbarV2/I");  
  babyTree_->Branch("passe_ll_ttbar",        &passe_ll_ttbar_,        "passe_ll_ttbar/I");  
  babyTree_->Branch("passe_ll_ttbarV1",      &passe_ll_ttbarV1_,      "passe_ll_ttbarV1/I");  
  babyTree_->Branch("passe_ll_ttbarV2",      &passe_ll_ttbarV2_,      "passe_ll_ttbarV2/I");  
  babyTree_->Branch("passe_ll_cand01",       &passe_ll_cand01_,       "passe_ll_cand01/I");  
  babyTree_->Branch("passm_lt_nom",          &passm_lt_nom_,          "passm_lt_nom/I");  
  babyTree_->Branch("passm_lt_nomttbar",     &passm_lt_nomttbar_,     "passm_lt_nomttbar/I");  
  babyTree_->Branch("passm_lt_nomttbarV2",   &passm_lt_nomttbarV2_,   "passm_lt_nomttbarV2/I");  
  babyTree_->Branch("passe_lt_ttbar",        &passe_lt_ttbar_,        "passe_lt_ttbar/I");  
  babyTree_->Branch("passe_lt_ttbarV1",      &passe_lt_ttbarV1_,      "passe_lt_ttbarV1/I");  
  babyTree_->Branch("passe_lt_ttbarV2",      &passe_lt_ttbarV2_,      "passe_lt_ttbarV2/I");  
  babyTree_->Branch("passe_lt_cand01",       &passe_lt_cand01_,       "passe_lt_cand01/I");  
  babyTree_->Branch("ptll",                  &ptll_,                  "ptll/F");  
  babyTree_->Branch("ptlt",                  &ptlt_,                  "ptlt/F");  
  babyTree_->Branch("idll",                  &idll_,                  "idll/F");  
  babyTree_->Branch("idlt",                  &idlt_,                  "idlt/F");  
  babyTree_->Branch("etall",                 &etall_,                 "etall/F");  
  babyTree_->Branch("etalt",                 &etalt_,                 "etalt/F");  
  babyTree_->Branch("phill",                 &phill_,                 "phill/F");  
  babyTree_->Branch("philt",                 &philt_,                 "philt/F");  
  babyTree_->Branch("dilmass",               &dilmass_,               "dilmass/F");  
  babyTree_->Branch("dilpt",                 &dilpt_,                 "dilpt/F");  
  babyTree_->Branch("flagll",                &flagll_,                "flagll/I");  
  babyTree_->Branch("flaglt",                &flaglt_,                "flaglt/I");  

}

//--------------------------------------------------------------------

void Z_looper::FillBabyNtuple ()
{
  babyTree_->Fill();
}

//--------------------------------------------------------------------

void Z_looper::CloseBabyNtuple ()
{
  babyFile_->cd();
  babyTree_->Write();
  babyFile_->Close();
}

//--------------------------------------------------------------------

void Z_looper::fillHistos(TH1F *h1[4],float value, float weight, int myType)
{

  fillUnderOverFlow(h1[myType], value, weight);      
  fillUnderOverFlow(h1[3],      value, weight);      
}

//--------------------------------------------------------------------

void Z_looper::fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{

  if( nJetsIdx > 2 ) nJetsIdx = 2;

  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      
}

//--------------------------------------------------------------------
