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

#include "../CORE/CMS2.h"
#include "../CORE/metSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/muonSelections.h"
#include "../Tools/goodrun.cc"
#include "../CORE/utilities.cc"
#include "histtools.h"
#include "../CORE/ttbarSelections.cc"
#include "../CORE/susySelections.cc"
#include "../CORE/jetSelections.cc"
//#include "../CORE/triggerUtils.cc"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
//#include "../CORE/mcSelections.cc"
#include "../CORE/mcSUSYkfactor.cc"

using namespace tas;
inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

enum metType   { e_tcmet = 0, e_tcmetNew = 1, e_pfmet = 2};
enum templateSource { e_QCD = 0, e_PhotonJet = 1 };

//--------------------------------------------------------------------

bool debug = false;
const int nJetBins              = 3;
const int nSumJetPtBins         = 7;
const int nBosonPtBins          = 4;
const bool generalLeptonVeto    = true;
const bool doTemplatePrediction = false;
metType myMetType               = e_tcmet;
templateSource myTemplateSource = e_PhotonJet;

float lumi         = 33.96e-3; 
char* iter         = "V00-00-00";
char* jsonfilename = "Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3_goodrun.txt";
//char* jsonfilename = "Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON_goodruns.txt";

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
    cout << "Error could not find boson pt bin x = " << x << endl;
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

bool isMuMuEvent(){

 
  if( evt_run() ==  146430 && evt_lumiBlock() ==        3  && evt_event() ==    460014 ) return true;
  if( evt_run() ==  147216 && evt_lumiBlock() ==       48  && evt_event() ==  35885648 ) return true;
  if( evt_run() ==  147217 && evt_lumiBlock() ==       75  && evt_event() ==  55188718 ) return true;
  if( evt_run() ==  148031 && evt_lumiBlock() ==      765  && evt_event() == 595250802 ) return true;
  if( evt_run() ==  147450 && evt_lumiBlock() ==       82  && evt_event() ==  29253181 ) return true;
  if( evt_run() ==  148029 && evt_lumiBlock() ==      534  && evt_event() == 414899947 ) return true;
  if( evt_run() ==  148862 && evt_lumiBlock() ==      350  && evt_event() == 522383338 ) return true;
  if( evt_run() ==  149181 && evt_lumiBlock() ==     1769  && evt_event() == 1675896175) return true;
  if( evt_run() ==  149182 && evt_lumiBlock() ==      167  && evt_event() == 145682218)  return true;
  if( evt_run() ==  149291 && evt_lumiBlock() ==      205  && evt_event() == 199787369)  return true;
  if( evt_run() ==  149291 && evt_lumiBlock() ==      232  && evt_event() == 235101408)  return true;
  if( evt_run() ==  149291 && evt_lumiBlock() ==      616  && evt_event() == 641847074)  return true;
  
  return false;
}

//--------------------------------------------------------------------

void printEvent(){

  cout << evt_dataset() << endl;
  cout << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;

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

  if( doTemplatePrediction ){

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

  }

  set_goodrun_file( jsonfilename );
 
  if( isData ){
    
    ofile_tcmet.open(  Form( "../output/%s/%s_tcmetprintout.txt" , iter , prefix  ) );
    ofile_events.open( Form( "../output/%s/%s_highmetevents.txt" , iter , prefix  ) );
    
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
  MakeBabyNtuple( Form("../output/%s/%s_baby.root", iter , prefix ) );

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

  int  nTot               = 0;  //total number of events
  
  int nGenPass60         = 0;  //number of events generated in sig window
  float nGenPass60_K     = 0;  //number of events generated in sig window
  int nRecoPassGenPass60 = 0;  //number of events reconstructed in sig window which pass gen
  int nRecoPassGenFail60 = 0;  //number of events reconstructed in sig window which fail gen

  int nGenPass120         = 0;  //number of events generated in sig window
  float nGenPass120_K     = 0;  //number of events generated in sig window
  int nRecoPassGenPass120 = 0;  //number of events reconstructed in sig window which pass gen
  int nRecoPassGenFail120 = 0;  //number of events reconstructed in sig window which fail gen

  float sigma      = 1;
  int   nTotEvents = 1;

  const int ncuts = 10;
  int nRecoPass_cut[ncuts];
  for( unsigned int icut = 0 ; icut < ncuts ; ++icut )
    nRecoPass_cut[icut] = 0;

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

      if( !isData ) sigma = cms2.evt_xsec_incl();

      nTot++;

      
      //if( !isMuMuEvent() ) continue;

      //if( !( evt_run()==142557 && evt_lumiBlock()==192 && evt_event()==80895186) ) continue;

      //cout << "Found event " << currentFile->GetTitle() <<  endl;

      float ksusy = 1;
      if( strcmp( prefix , "LM4" ) == 0 ) ksusy = kfactorSUSY( "lm4" );
      if( strcmp( prefix , "LM8" ) == 0 ) ksusy = kfactorSUSY( "lm8" );

      if( PassGenSelection( isData ) > 60. ){
	nGenPass60++;
	nRecoPass_cut[0]++;
	nGenPass60_K += ksusy;
      }
      if( PassGenSelection( isData ) > 120. ){
	nGenPass120++;
	nRecoPass_cut[0]++;
	nGenPass120_K += ksusy;
      }

      if( !isData ){
	hresponse->Fill( gen_met() , evt_pfmet() / gen_met() );
      	hgenmet_all->Fill( gen_met() );
	if( evt_pfmet() > 60 ) hgenmet_pass->Fill( gen_met() );
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
      
      if( !isData ){

	//splice together the DY samples - if its madgraph, then we do nothing
	if(TString(prefix).Contains("DY") && TString(evt_dataset()).Contains("madgraph") == false) {	
	  bool doNotContinue = false;
	  for(unsigned int i = 0; i < genps_p4().size(); i++){
	    if(abs(genps_id()[i]) == 23 && genps_p4()[i].M() > 50.)
	      doNotContinue = true;
	  }
	  if(doNotContinue)
	    continue;	
	}
	
	//extract pthat
	if(TString(prefix).Contains("DY")){
	  int nz = 0;
	  for(unsigned int i = 0; i < genps_p4().size(); i++){
	    if(abs(genps_id()[i]) == 23){
	      mllgen_ = genps_p4()[i].M();
	      nz++;
	    }
	  }
	  if(nz != 1 ) cout << "ERROR NZ " << nz << endl;
	}
      }


      //good run+event selection-----------------------------------------------------------
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      if( !cleaning_standardAugust2010( isData) )                    continue;

      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[1]++;
      
      if(debug) cout << "Pass event selection" << endl;

      InitBabyNtuple();

      // event stuff
      strcpy(dataset_, cms2.evt_dataset().Data());
      run_    = cms2.evt_run();
      lumi_   = cms2.evt_lumiBlock();
      event_  = cms2.evt_event();

      weight_ = 1.;
      pthat_  = -1;
      if( !isData ){
        weight_ = cms2.evt_scale1fb() * kFactor * lumi;
	if( TString(prefix).Contains("LM") ){
	  if( strcmp( prefix , "LM0" ) == 0 ) weight_ *= kfactorSUSY( "lm0" );
	  if( strcmp( prefix , "LM1" ) == 0 ) weight_ *= kfactorSUSY( "lm1" );
	  if( strcmp( prefix , "LM2" ) == 0 ) weight_ *= kfactorSUSY( "lm2" );
	  if( strcmp( prefix , "LM3" ) == 0 ) weight_ *= kfactorSUSY( "lm3" );
	  if( strcmp( prefix , "LM4" ) == 0 ) weight_ *= kfactorSUSY( "lm4" );
	  if( strcmp( prefix , "LM5" ) == 0 ) weight_ *= kfactorSUSY( "lm5" );
	  if( strcmp( prefix , "LM6" ) == 0 ) weight_ *= kfactorSUSY( "lm6" );
	  if( strcmp( prefix , "LM7" ) == 0 ) weight_ *= kfactorSUSY( "lm7" );
	  if( strcmp( prefix , "LM8" ) == 0 ) weight_ *= kfactorSUSY( "lm8" );
	  if( strcmp( prefix , "LM9" ) == 0 ) weight_ *= kfactorSUSY( "lm9" );
	}
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

      VofP4 goodLeptons;
      vector<int> goodMuonIndices;
      vector<int> goodPFMuonIndices;
      vector<bool>  killedJet;
      goodLeptons.clear();
      killedJet.clear();
      goodMuonIndices.clear();
      goodPFMuonIndices.clear();

      if( generalLeptonVeto ){
        
        for( unsigned int iel = 0 ; iel < els_p4().size(); ++iel ){
          if( els_p4().at(iel).pt() < 20 )                                                 continue;
          if( !pass_electronSelection( iel , electronSelection_el_OSV1 , false , false ) ) continue;
          goodLeptons.push_back( els_p4().at(iel) );
          killedJet.push_back( false );
        }
        
        for( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu ){
          if( mus_p4().at(imu).pt() < 20 )           continue;
          if( !muonId( imu , OSZ_v1 ))               continue;
          goodLeptons.push_back( mus_p4().at(imu) );
          killedJet.push_back( false );
	  goodMuonIndices.push_back( imu );
	  int ipf = mus_pfmusidx().at(imu);
	  if( ipf < pfmus_p4().size() && ipf >= 0 ){
	    goodPFMuonIndices.push_back( imu );
	  }
        }

      }


      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); ++hypIdx) {

        if( !passSUSYTrigger_v1( isData , hyp_type()[hypIdx] ) ) continue;

        //check that hyp leptons come from same vertex
        if(!hypsFromSameVtx(hypIdx))   continue;
        
        //selection--------------------------------------------------------------------- 

        //OS, pt > (20,20) GeV, dilmass > 10 GeV
        if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                             continue;
        if( TMath::Max( hyp_ll_p4()[hypIdx].pt() , hyp_lt_p4()[hypIdx].pt() ) < 20. )   continue;
        if( TMath::Min( hyp_ll_p4()[hypIdx].pt() , hyp_lt_p4()[hypIdx].pt() ) < 20. )   continue;
        if( hyp_p4()[hypIdx].mass() < 10 )                                              continue;

        //nominal muon ID
        if (abs(hyp_ll_id()[hypIdx]) == 13  && !( muonId( hyp_ll_index()[hypIdx] , OSZ_v1 )))   continue;
        if (abs(hyp_lt_id()[hypIdx]) == 13  && !( muonId( hyp_lt_index()[hypIdx] , OSZ_v1 )))   continue;
        
        //OSV1 electron ID
        if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_el_OSV1 , false , false ))) continue;
        if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_el_OSV1 , false , false ))) continue;
        
        //Z-mass constraint
        //if( hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  continue;

        nHypPass++;
      
        v_goodHyps.push_back( hypIdx );
      
      }

      if( v_goodHyps.size() == 0 ) continue;

      unsigned int hypIdx = selectBestZHyp(v_goodHyps);

      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[2]++;

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
      metError_         = getMetError( goodMuonIndices );
      metErrorC_        = getMetError_claudio( goodMuonIndices );


      float mom1 = hyp_ll_p4()[hypIdx].P();
      float mom2 = hyp_lt_p4()[hypIdx].P();
      if( mom2 > mom1 ){
        mom1 = mom2;
        mom2 = mom1;
      }

      float one_minus_cosalpha = (dilmass_ * dilmass_) / (2 * mom1 * mom2 );
      
      dpdm_ = dilmass_ /( mom2 * one_minus_cosalpha );


      npfmuons_ = 0;
      nmatchedpfmuons_ = 0;

      for( unsigned int ipf = 0 ; ipf < pfmus_p4().size() ; ipf++ ){
	if( pfmus_p4().at(ipf).pt()  < 20. )         continue;
	if( fabs( pfmus_p4().at(ipf).eta() ) > 2.5 ) continue;
	npfmuons_++;

      }
      
      if( abs( hyp_ll_id()[hypIdx] ) == 13 ){
	int ipf_ll = mus_pfmusidx().at(hyp_ll_index().at(hypIdx));
	if( ipf_ll >= pfmus_p4().size() || ipf_ll < 0 ){
	  //cout << "Error, pfmuon ll index out of range " << ipf_ll << endl;
          //printEvent();
	}else{
	  ptll_pf_ = pfmus_p4().at(ipf_ll).pt();
	  nmatchedpfmuons_ ++;
	  if( fabs( ptll_ - ptll_pf_ ) > 0.01 ){
	    //cout << "ERROR: " << ptll_ << " " << ptll_pf_ << endl;
	  }
	}

	int muidx  = hyp_ll_index().at(hypIdx);
	ptlltrk_   = mus_trk_p4().at(muidx).pt();
	ptllgfit_  = mus_gfit_p4().at(muidx).pt();
	pterrll_   = mus_ptErr().at(muidx);

        if( !isData ){
          int mcid   = mus_mc3idx().at(muidx);
          //cout << "ll mcid " << mcid << endl;
          if( mcid >=0 && mcid < genps_p4().size() ){
            ptllgen_   = genps_p4().at(mcid).pt();
          }else{
            //cout << "Error, ll MC index " << mcid << " size " << genps_p4().size() << endl;
          }
        }
      }

      if( abs( hyp_lt_id()[hypIdx] ) == 13 ){
	int ipf_lt = mus_pfmusidx()[hyp_lt_index()[hypIdx]];
	if( ipf_lt >= pfmus_p4().size() || ipf_lt < 0 ){
	  //cout << "Error, pfmuon lt index out of range " << ipf_lt << endl;
          //printEvent();
	}else{
	  ptlt_pf_ = pfmus_p4().at(ipf_lt).pt();
	  nmatchedpfmuons_ ++;
	  if( fabs( ptlt_ - ptlt_pf_ ) > 0.01 ){
	    //cout << "ERROR: " << ptlt_ << " " << ptlt_pf_ << endl;
	  }
	}

	int muidx  = hyp_lt_index().at(hypIdx);
	ptlttrk_   = mus_trk_p4().at(muidx).pt();
	ptltgfit_  = mus_gfit_p4().at(muidx).pt();
	pterrlt_   = mus_ptErr().at(muidx);

        if( !isData ){
          int mcid   = mus_mc3idx().at(muidx);
          //cout << "lt mcid " << mcid << endl;
          if( mcid >=0 && mcid < genps_p4().size() ){
            ptltgen_   = genps_p4().at(mcid).pt();
          }else{
            //cout << "Error, lt MC index " << mcid << " size " << genps_p4().size() <<endl;
          }
        }

      }

      if( leptype_ == 1 && npfmuons_ >= 2 ){

	int ipf_ll = mus_pfmusidx().at(hyp_ll_index().at(hypIdx));
	int ipf_lt = mus_pfmusidx().at(hyp_lt_index().at(hypIdx));
	
	dilmasspf_ = -9999;
	
	if( ipf_lt >= pfmus_p4().size() || ipf_ll >= pfmus_p4().size() ){
	  cout << "Error, pfmuon out of range: SHOULDN'T GET HERE!!" << endl;
	}else{
	  dilmasspf_ = ( pfmus_p4().at(ipf_ll) + pfmus_p4().at(ipf_lt) ).mass();
	}
      }


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

      //pfjet matched to electron (ll)

      if( abs( hyp_ll_id()[hypIdx] ) == 11 ){

        float drjet_ll = 100;
        int   ijet_ll  = -1;

        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          
          if( vjet.pt()  < 10  )             continue;
          if( fabs(vjet.eta()) > 2.5 )       continue;
          
          float dr = dRbetweenVectors(vjet, vll);
          
          if( dr < drjet_ll ){
            drjet_ll    = dr;
            ijet_ll     = ijet;
          }
        }
      
        if( ijet_ll >= 0 ){
          drjet_ll_   = drjet_ll;
          jetpt_ll_   = pfjets_cor().at(ijet_ll) * pfjets_p4().at(ijet_ll).pt();
          pfjetid_ll_ = passesPFJetID( ijet_ll ) ? 1 : 0;
        }
      }

      //pfjet matched to electron (lt)

      if( abs( hyp_lt_id()[hypIdx] ) == 11 ){

        float drjet_lt = 100;
        int   ijet_lt  = -1;

        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          
          if( vjet.pt()  < 10  )             continue;
          if( fabs(vjet.eta()) > 2.5 )       continue;
          
          float dr = dRbetweenVectors(vjet, vlt);
          
          if( dr < drjet_lt ){
            drjet_lt    = dr;
            ijet_lt     = ijet;
          }
        }

        if( ijet_lt >= 0 ){
          drjet_lt_   = drjet_lt;
          jetpt_lt_   = pfjets_cor().at(ijet_lt) * pfjets_p4().at(ijet_lt).pt();
          pfjetid_lt_ = passesPFJetID( ijet_lt ) ? 1 : 0;
        }
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

      // electron energy scale stuff

      dilmasscor_  = dilmass_;
      tcmetcor_    = tcmet_;
      pfmetcor_    = pfmet_;

      if( leptype_ == 0 ){

        if     ( fabs( etall_ ) < 1.474 && fabs( etalt_ ) < 1.479 ) ecaltype_ = 1; //EB-EB
        else if( fabs( etall_ ) > 1.474 && fabs( etalt_ ) > 1.479 ) ecaltype_ = 2; //EE-EE
        else                                                        ecaltype_ = 3; //EB-EE
        
        if( isData ){

          LorentzVector vllcor;
          LorentzVector vltcor;

          float tcmetcor_x = tcmet_ * cos( tcmetphi_ );
          float tcmetcor_y = tcmet_ * sin( tcmetphi_ );

          float pfmetcor_x = pfmet_ * cos( pfmetphi_ );
          float pfmetcor_y = pfmet_ * sin( pfmetphi_ );
       
          //ll correction
          if( fabs( etall_ ) < 1.474 ){ //correct for EB electron
            vllcor      = 1.01 * hyp_ll_p4()[hypIdx];  
            
            tcmetcor_x -= 0.01 * hyp_ll_p4()[hypIdx].x();
            tcmetcor_y -= 0.01 * hyp_ll_p4()[hypIdx].y();

            pfmetcor_x -= 0.01 * hyp_ll_p4()[hypIdx].x();
            pfmetcor_y -= 0.01 * hyp_ll_p4()[hypIdx].y();
          }
          else{                         //correct for EE electron
            vllcor      = 1.03 * hyp_ll_p4()[hypIdx]; 

            tcmetcor_x -= 0.03 * hyp_ll_p4()[hypIdx].x();
            tcmetcor_y -= 0.03 * hyp_ll_p4()[hypIdx].y();

            pfmetcor_x -= 0.03 * hyp_ll_p4()[hypIdx].x();
            pfmetcor_y -= 0.03 * hyp_ll_p4()[hypIdx].y();
          }

          //lt correction
          if( fabs( etalt_ ) < 1.474 ){ //correct for EB electron
            vltcor      = 1.01 * hyp_lt_p4()[hypIdx];  
            
            tcmetcor_x -= 0.01 * hyp_lt_p4()[hypIdx].x();
            tcmetcor_y -= 0.01 * hyp_lt_p4()[hypIdx].y();

            pfmetcor_x -= 0.01 * hyp_lt_p4()[hypIdx].x();
            pfmetcor_y -= 0.01 * hyp_lt_p4()[hypIdx].y();
          }
          else{                         //correct for EE electron
            vltcor      = 1.03 * hyp_lt_p4()[hypIdx]; 

            tcmetcor_x -= 0.03 * hyp_lt_p4()[hypIdx].x();
            tcmetcor_y -= 0.03 * hyp_lt_p4()[hypIdx].y();

            pfmetcor_x -= 0.03 * hyp_lt_p4()[hypIdx].x();
            pfmetcor_y -= 0.03 * hyp_lt_p4()[hypIdx].y();
          }
          
          dilmasscor_  = ( vllcor + vltcor ).mass();
          tcmetcor_    = sqrt( tcmetcor_x * tcmetcor_x + tcmetcor_y * tcmetcor_y );
          pfmetcor_    = sqrt( pfmetcor_x * pfmetcor_x + pfmetcor_y * pfmetcor_y );
        }
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

      failjetid_ =  0;
      maxemf_    = -1;

      /*
      //store in goodpfjetsXX all pfjets with pt > XX GeV, eta < 2.5, passesPFJetID 
      VofP4 goodpfjets30;
      goodpfjets30.clear();
      VofP4 goodpfjets15;
      goodpfjets15.clear();

      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {

        LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
      
        if( fabs( vjet.eta() ) > 2.5 )   continue;
        if( !passesPFJetID(ijet) )       continue;
        
        if( vjet.pt() < 15 )             continue;
        goodpfjets15.push_back( vjet );
        if( vjet.pt() < 30 )             continue;
        goodpfjets30.push_back( vjet );
      }


      //veto the closest good pfjet to all good leptons if dr < 0.4
      vector<unsigned int> vetoJet;
      vetoJet.clear();

      if( generalLeptonVeto ){
       
        for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){

          float drmin = 100;
          int   imin  = -1;

          for (unsigned int ijet = 0 ; ijet < goodpfjets30.size() ; ijet++){
            
            float dr = dRbetweenVectors( goodpfjets30.at(ijet) , goodLeptons.at(ilep) );
            
            if( dr < drmin ){
              drmin = dr;
              imin  = ijet;
            }
          }
          
          if( drmin < 0.4 ) vetoJet.push_back( imin );
          
        }
      }

      drll_ = 100;

      int   iminll    = -1;
      float mindrll   = 100.;

      for ( unsigned int ijet = 0 ; ijet < goodpfjets30.size() ; ijet++ ){
        float dr = dRbetweenVectors( goodpfjets30.at(ijet) , hyp_ll_p4()[hypIdx] );
        if( dr < mindrll ){
          iminll  = ijet;
          mindrll = dr;
        }
      }

      int iminll2    = -1;
      float mindrll2 = 100.;

      for ( unsigned int ijet = 0 ; ijet < goodpfjets30.size() ; ijet++ ){
        if( ijet == iminll ) continue;
        float dr = dRbetweenVectors( goodpfjets30.at(ijet) , hyp_ll_p4()[hypIdx] );
        if( dr < mindrll2 ){
          iminll2  = ijet;
          mindrll2 = dr;
        }
      }

      drll_ = mindrll2;
      */


      //count JPT's
      nJPT_ = 0;

      for ( unsigned int ijet = 0 ; ijet < jpts_p4().size() ; ijet++ ) {
          
        LorentzVector vjet = jpts_cor().at(ijet) * jpts_p4().at(ijet);
        LorentzVector vlt  = hyp_lt_p4()[hypIdx];
        LorentzVector vll  = hyp_ll_p4()[hypIdx];

        if( fabs( vjet.eta() ) > 2.5 )           continue;
        if( vjet.pt()  < 30.         )           continue;
        if( !passesCaloJetID( vjet ) )           continue;

        if( generalLeptonVeto ){
          bool rejectJet = false;
          for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
            if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
          }
          if( rejectJet ) continue;
        }

        if( dRbetweenVectors(vjet, vll) < 0.4 ) continue;
        if( dRbetweenVectors(vjet, vlt) < 0.4 ) continue;

        /*
        if( generalLeptonVeto ){
          
          bool rejectJet = false;
          for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
            if( killedJet.at( ilep ) ) continue;
            if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ){
              rejectJet            = true;  
              killedJet.at( ilep ) = true;
            }
          }
          
          if( rejectJet ) continue;
       
        }

        else{
          
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          
          if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
       
        }
        */

        nJPT_++;
      } 
        
      //reset killedJet vector
      for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
        killedJet.at( ilep ) = false;
      }
        
      //loop over pfjets pt > 30 GeV |eta| < 2.5
      for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
        LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
        LorentzVector vlt  = hyp_lt_p4()[hypIdx];
        LorentzVector vll  = hyp_ll_p4()[hypIdx];

        if( fabs( vjet.eta() ) > 2.5 )           continue;
     
        /*
        if( generalLeptonVeto ){
        
          bool rejectJet = false;
          for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
            if( killedJet.at( ilep ) ) continue;
            if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ){
              rejectJet            = true;  
              killedJet.at( ilep ) = true;
            }
          }
        
          if( rejectJet ) continue;
       
        }

        else{
          
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          
          if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
       
        }
        */

        if( generalLeptonVeto ){
          bool rejectJet = false;
          for( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ ){
            if( dRbetweenVectors( vjet , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;  
          }
          if( rejectJet ) continue;
        }

        if( dRbetweenVectors(vjet, vll) < 0.4 ) continue;
        if( dRbetweenVectors(vjet, vlt) < 0.4 ) continue;


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

          float emfrac = pfjets_neutralEmE().at(ijet) / pfjets_p4().at(ijet).energy();
          if( emfrac > maxemf_ ) maxemf_ = emfrac;
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

      //-----------------------
      //signal region selection
      //-----------------------

      if( dilmass_ < 81. || dilmass_ > 101. )    continue;
      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[3]++;

      if( leptype_ == 2 )                        continue;
      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[4]++;

      if( nJets_ < 2 )                           continue;
      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[5]++;

      if( pfmet_ < 60 )                          continue;
      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[6]++;

      if( leptype_ == 0 ){
	if( jetpt_ll_ - ptll_ < -5  ) continue; 
	if( jetpt_lt_ - ptlt_ < -5  ) continue; 
      }
      if( PassGenSelection( isData ) > 60. )   nRecoPass_cut[7]++;
      

      if( pfmet_ > 60 ){
	if( PassGenSelection( isData ) > 60. ) nRecoPassGenPass60++;
	else                           nRecoPassGenFail60++;
      }
      if( pfmet_ > 120 ){
	if( PassGenSelection( isData ) > 120. ) nRecoPassGenPass120++;
	else                            nRecoPassGenFail120++;
      }


      continue;
      /*
      //met templates-------------------------------------------------------------------- 
 
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
        
      TH1F* hmet = new TH1F();
      int iJetBin      = getJetBin( nJets_ );

      if( doTemplatePrediction ){

        //fill predicted and observed met histos
  
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        int iBosonPtBin  = getBosonPtBin( dilpt_ );
        
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
        
      }
      else{
        hmet = (TH1F*) metObserved->Clone("hdummy");
        hmet->Reset();
      }

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
      */
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

  cout << endl;
  cout << "nTot                " << nTot               << endl;

  cout << "Loose signal region-----------------------" << endl;
  cout << "nGenPass            " << nGenPass60         << endl;
  cout << "nGenPass X K        " << nGenPass60_K       << endl;
  cout << "nRecoPassGenPass    " << nRecoPassGenPass60 << endl;
  cout << "nRecoPassGenFail    " << nRecoPassGenFail60 << endl;
  cout << "Efficiency          " << nRecoPassGenPass60 / (float) nGenPass60 << endl;

  cout << "Tight signal region-----------------------"  << endl;
  cout << "nGenPass            " << nGenPass120         << endl;
  cout << "nGenPass X K        " << nGenPass120_K       << endl;
  cout << "nRecoPassGenPass    " << nRecoPassGenPass120 << endl;
  cout << "nRecoPassGenFail    " << nRecoPassGenFail120 << endl;
  cout << "Efficiency          " << nRecoPassGenPass120 / (float) nGenPass120 << endl;

  cout << "Cut flow----------------------------------" << endl;
  for( int icut = 0 ; icut < 8 ; icut++ ){
    float eff = 1;
    if( icut > 0 ) eff = nRecoPass_cut[icut] / (float) nRecoPass_cut[icut-1];
    cout << "nRecoPass cut " << icut << " " << nRecoPass_cut[icut] << " " << Form("%.2f",eff) << endl;
  }

  CloseBabyNtuple();

  already_seen.clear();

  // make histos rootfile
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist( Form("../output/%s/%s%s.root", iter , prefix , metTemplateString.c_str() ) );
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

  //electron-matched jet stuff
  drjet_ll_       = -999999.;
  jetpt_ll_       = -999999.;
  pfjetid_ll_     = -999999;
  drjet_lt_       = -999999.;
  jetpt_lt_       = -999999.;
  pfjetid_lt_     = -999999;

  // event stuff
  run_          = -999999;
  memset(dataset_, '\0', 200);
  lumi_         = -999999;
  event_        = -999999;
  weight_       = -999999.;
  pthat_        = -999999.;
  mllgen_       = -999999.;
  nGoodVertex_  = -999999;
  leptype_      = -999999;
  ecaltype_     = -999999;

  // genmet stuff
  genmet_     = -999999.;
  genmetphi_  = -999999.;
  gensumet_   = -999999.;

  // pfmet stuff
  pfmet_     = -999999.;
  pfmetphi_  = -999999.;
  pfsumet_   = -999999.;

  //pfmuon stuff
  npfmuons_         = -999999;
  nmatchedpfmuons_  = -999999;
  ptll_pf_   =  999999.;
  ptlt_pf_   =  999999.;

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
  tcmetcor_     = -999999.;
  pfmetcor_     = -999999.;
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

  metError_         = -999999;
  metErrorC_        = -999999;
  dpdm_             = -999999;

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
  pterrll_         = -999999;
  pterrlt_         = -999999;
  ptlltrk_         = -999999;
  ptllgen_         = -999999;
  ptltgen_         = -999999;
  ptlttrk_         = -999999;
  ptllgfit_        = -999999;
  ptltgfit_        = -999999;
  idll_            = -999999;
  idlt_            = -999999;
  etall_           = -999999;
  etalt_           = -999999;
  phill_           = -999999;
  philt_           = -999999;
  dilmass_         = -999999.;
  dilmasspf_       = -999999.;
  dilmasscor_      = -999999.;
  dilpt_           = -999999.;
  flagll_          = -999999;
  flaglt_          = -999999;
  failjetid_       = -999999;
  maxemf_          = -999999.;

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

  hgenmet_all  = new TH1F("hgenmet_all","",100,0,200);
  hgenmet_pass = new TH1F("hgenmet_pass","",100,0,200);
  hresponse    = new TProfile("hresponse","",100,0,200,0,2);

  hgenmet_all->Sumw2();
  hgenmet_pass->Sumw2();

  Int_t maxmet = 200;

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
  babyTree_->Branch("dataset",      &dataset_,      "dataset[200]/C");
  babyTree_->Branch("run",          &run_,          "run/I"  );
  babyTree_->Branch("lumi",         &lumi_,         "lumi/I" );
  babyTree_->Branch("event",        &event_,        "event/I");
  babyTree_->Branch("failjetid",    &failjetid_,    "failjetid/I");
  babyTree_->Branch("maxemf",       &maxemf_,       "maxemf/F");
  babyTree_->Branch("nvtx",         &nGoodVertex_,  "nvtx/I");
  babyTree_->Branch("weight",       &weight_,       "weight/F");
  babyTree_->Branch("pthat",        &pthat_,        "pthat/F");
  babyTree_->Branch("mllgen",       &mllgen_,       "mllgen/F");

  //electron-matched jet stuff
  babyTree_->Branch("drjetll",      &drjet_ll_,     "drjetll/F"     );
  babyTree_->Branch("jetptll",      &jetpt_ll_,     "jetptll/F"     );
  babyTree_->Branch("pfjetidll",    &pfjetid_ll_,   "pfjetidll/I"   );
  babyTree_->Branch("drjetlt",      &drjet_lt_,     "drjetlt/F"     );
  babyTree_->Branch("jetptlt",      &jetpt_lt_,     "jetptlt/F"     );
  babyTree_->Branch("pfjetidlt",    &pfjetid_lt_,   "pfjetidlt/I"   );


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
  babyTree_->Branch("dphixmet",     &dphixmet_,     "dphixmet/F"    );
  babyTree_->Branch("metpar",       &metPar_,       "metpar/F"      );
  babyTree_->Branch("metperp",      &metPerp_,      "metperp/F"     );
  babyTree_->Branch("tcmet",        &tcmet_,        "tcmet/F"      );
  babyTree_->Branch("tcmetphi",     &tcmetphi_,     "tcmetphi/F"   );
  babyTree_->Branch("tcsumet",      &tcsumet_,      "tcsumet/F"    );
  babyTree_->Branch("tcmetNew",     &tcmetNew_,     "tcmetNew/F"      );
  babyTree_->Branch("tcmetphiNew",  &tcmetphiNew_,  "tcmetphiNew/F"   );
  babyTree_->Branch("tcsumetNew",   &tcsumetNew_,   "tcsumetNew/F"    );
  babyTree_->Branch("tcmetcor",     &tcmetcor_,     "tcmetcor/F");  
  babyTree_->Branch("pfmetcor",     &pfmetcor_,     "pfmetcor/F");  

  //jet stuff
  babyTree_->Branch("njets",          &nJets_,            "njets/I"       );
  babyTree_->Branch("njpt",           &nJPT_,             "njpt/I"        );
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
  babyTree_->Branch("ecaltype",              &ecaltype_,              "ecaltype/I");
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
  babyTree_->Branch("dpdm",                  &dpdm_,                  "dpdm/F");  
  babyTree_->Branch("meterror",              &metError_,              "metError/F");  
  babyTree_->Branch("meterrorc",             &metErrorC_,             "metErrorc/F");  
  babyTree_->Branch("ptll",                  &ptll_,                  "ptll/F");  
  babyTree_->Branch("ptlt",                  &ptlt_,                  "ptlt/F");  
  babyTree_->Branch("pterrll",               &pterrll_,               "pterrll/F");  
  babyTree_->Branch("pterrlt",               &pterrlt_,               "pterrlt/F");  
  babyTree_->Branch("ptlltrk",               &ptlltrk_,               "ptlltrk/F");  
  babyTree_->Branch("ptlttrk",               &ptlttrk_,               "ptlttrk/F");  
  babyTree_->Branch("ptllgfit",              &ptllgfit_,              "ptllgfit/F");  
  babyTree_->Branch("ptltgfit",              &ptltgfit_,              "ptltgfit/F");  
  babyTree_->Branch("ptllpf",                &ptll_pf_,               "ptllpf/F");  
  babyTree_->Branch("ptltpf",                &ptlt_pf_,               "ptltpf/F");  
  babyTree_->Branch("ptllgen",               &ptllgen_,               "ptllgen/F");  
  babyTree_->Branch("ptltgen",               &ptltgen_,               "ptltgen/F");  
  babyTree_->Branch("npfmuons",              &npfmuons_,              "npfmuons/I");  
  babyTree_->Branch("nmatchedpfmuons",       &nmatchedpfmuons_,       "nmatchedpfmuons/I");  
  babyTree_->Branch("idll",                  &idll_,                  "idll/I");  
  babyTree_->Branch("idlt",                  &idlt_,                  "idlt/I");  
  babyTree_->Branch("etall",                 &etall_,                 "etall/F");  
  babyTree_->Branch("etalt",                 &etalt_,                 "etalt/F");  
  babyTree_->Branch("phill",                 &phill_,                 "phill/F");  
  babyTree_->Branch("philt",                 &philt_,                 "philt/F");  
  babyTree_->Branch("dilmass",               &dilmass_,               "dilmass/F");  
  babyTree_->Branch("dilmasspf",             &dilmasspf_,             "dilmasspf/F");  
  babyTree_->Branch("dilmasscor",            &dilmasscor_,            "dilmasscor/F");  
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

float Z_looper::getMetError(  vector<int> goodMuonIndices ){

  //error is zero if no muons
  if( goodMuonIndices.size() == 0 ) return 0.;

  float metError = 0.;

  float met_x = evt_pfmet() * cos( evt_pfmetPhi() );
  float met_y = evt_pfmet() * sin( evt_pfmetPhi() );
  float met   = evt_pfmet();

  //loop over all muons in goodMuonIndices (pt > 20 GeV, |eta| < 2.4, OSZ_v1)
  for( unsigned int i = 0 ; i < goodMuonIndices.size() ; ++i ){
    int imu = goodMuonIndices.at(i);

    //calculate uncertainties in muon x,y momentum components
    //I am assuming that dpx = dpt * cos(phi)
    float dpx = cos( mus_p4().at(imu).phi() ) * mus_ptErr().at(imu);
    float dpy = sin( mus_p4().at(imu).phi() ) * mus_ptErr().at(imu);

    //check how much met changes when I alter the muon pt by (dpx,dpy)
    float metxprime = met_x + dpx;
    float metyprime = met_y + dpy;
    float metprime  = sqrt( metxprime * metxprime + metyprime * metyprime );
    float dmet      = metprime - met;

    //sum the errors from each muon in quadrature
    metError += dmet * dmet;
    
  }

  metError = sqrt( metError );

  return metError;
}

//--------------------------------------------------------------------

float Z_looper::getMetError_claudio(  vector<int> goodMuonIndices ){
   
  //error is zero if no muons
  if( goodMuonIndices.size() == 0 ) return 0.;
  

  float met_x = evt_pfmet() * cos( evt_pfmetPhi() );
  float met_y = evt_pfmet() * sin( evt_pfmetPhi() );
  float met   = evt_pfmet();

  float metError = 0.;

  //loop over all muons in goodMuonIndices (pt > 20 GeV, |eta| < 2.4, OSZ_v1)
  for( unsigned int i = 0 ; i < goodMuonIndices.size() ; ++i ){
    int imu = goodMuonIndices.at(i);
 
    float phi  = mus_p4().at(imu).phi();
    float dpt  = mus_ptErr().at(imu);
    float dmet = ( ( cos(phi) * met_x + sin(phi) * met_y ) * dpt ) / met;
    metError += dmet * dmet;
  
  }
  
  metError = sqrt(metError);
  
  return metError;
  
}

//--------------------------------------------------------------------

float Z_looper::PassGenSelection( bool isData ){
  
  if( isData ) return -999.;

  //---------------------------------------------
  // does this event pass the analysis selection?
  //---------------------------------------------
  
  // mc leptons
  std::vector<unsigned int> mcLeptonIndices;
  int nGoodLep = 0;
  for (size_t i = 0; i < cms2.genps_id().size(); ++i){
    
    //electron or muon
    if (!(abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13))      continue;

    //pt > 20 GeV, |eta| < 2.5
    if ( cms2.genps_p4()[i].Pt() < 20.0 || abs(cms2.genps_p4()[i].Eta()) > 2.5) continue;

    nGoodLep++;
    mcLeptonIndices.push_back(i);
  }

  if( nGoodLep < 2 ) return -1.;

  //look for OS pt > 20,10 GeV pair Z mass veto
  bool foundPair = false;

  for( unsigned int i = 0 ; i < mcLeptonIndices.size() ; ++i ){
    unsigned int ilep = mcLeptonIndices.at(i);
    for( unsigned int j = i + 1 ; j < mcLeptonIndices.size() ; ++j ){
      unsigned int jlep = mcLeptonIndices.at(j);

	//OS
	if ( cms2.genps_id()[ilep] * cms2.genps_id()[jlep] > 0 )                            continue;

	//SF
	if ( abs( cms2.genps_id()[ilep] ) != abs( cms2.genps_id()[jlep] ) )                 continue;

	//Z mass 81-101 GeV
	float dilmass = ( cms2.genps_p4()[ilep] + cms2.genps_p4()[jlep] ).mass();
	if( dilmass < 81.0 || dilmass > 101. ) continue;
	
	//found OS pair!
	foundPair = true;
	     
      }
    }

    if( !foundPair ) return -2.;
   
    // mc jets
    
    int nGoodJet   = 0;
    float sumJetPt = 0.;
    for (size_t j = 0; j < cms2.evt_ngenjets(); ++j) 
    {
        if (cms2.genjets_p4()[j].Pt() < 30.0)       continue;
        if (fabs(cms2.genjets_p4()[j].Eta()) > 2.5) continue;
        bool clean = true;
        for ( size_t i = 0; i < mcLeptonIndices.size(); ++i) 
        {
            if (ROOT::Math::VectorUtil::DeltaR(cms2.genjets_p4()[j], cms2.genps_p4()[mcLeptonIndices[i]]) < 0.4) {
                clean = false;
                break;
            }
        }
        if (clean){
	  nGoodJet ++;
	  sumJetPt += genjets_p4()[j].Pt();
	}
    }
    
    if( nGoodJet < 2          ) return -3.;

    //hooray! return met
    return cms2.gen_met();
    
}
