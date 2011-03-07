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

char* iter          = "default";
bool makebaby       = true;
bool debug          = false;
bool calculateTCMET = true;

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
         
          //APPLY BASIC EVENT SELECTIONS
          //TRIGGER, TRACKING AND VERTEX
      
	  if (!isData || goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))  continue;
          if( !cleaning_standardAugust2010( isData) )                    continue;
          
          // N.B. BABY NTUPLE IS FILLED
          // FOR EACH EVENT

          InitBabyNtuple();

          // event stuff
          run_    = cms2.evt_run();
          lumi_   = cms2.evt_lumiBlock();
          event_  = cms2.evt_event();

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
