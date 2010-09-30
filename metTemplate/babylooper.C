#include "babylooper.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>


#include "TChain.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

enum metType   { e_tcmet = 0, e_tcmetNew = 1, e_pfmet = 2};
enum templateSource { e_QCD = 0, e_PhotonJet = 1 };

//------------USER PARAMS-------------------------
bool debug                       = true;
bool reweight                    = false;
metType myMetType                = e_pfmet;
templateSource myTemplateSource  = e_PhotonJet;

const int nJetBins               = 3;
const int nSumJetPtBins          = 7;
const int nBosonPtBins           = 4;

char* iter                       = "V01-02";
bool useCombinedTemplates        = true;
const int nMetEntries            = 50;
//------------------------------------------------

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

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

string jetString( int bin ){

  string js;

  if     ( bin == 0 ) js = "njets = 1";
  else if( bin == 1 ) js = "njets = 2";
  else if( bin == 2 ) js = "njets #geq 3";
  else{
    cout << "Error invalid bin passed as argument to met templates" << endl;
    exit(0);
  }

  return js;

}

int getSumJetPtBin( float x ){

  //bins array defines the sumJetPt binning
  float bins[nSumJetPtBins+1]={0,25,50,75,100,150,250,5000};
  
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

  float bins[nSumJetPtBins+1]={0,25,50,75,100,150,250,5000};

  stringstream s;
  s << bins[bin] << " < sumJetPt < " << bins[bin+1] << " GeV";

  return s.str();
}

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

/*
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
  
  int bin = njets;
  if( bin >= nJetBins ) bin = nJetBins - 1;
  
  return bin;
}

//--------------------------------------------------------------------

string jetString( int bin ){

  stringstream s;
  s << "nJets = " << bin ;

  return s.str();

}

//--------------------------------------------------------------------

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

//--------------------------------------------------------------------

string sumJetPtString( int bin ){
  
  float bins[nSumJetPtBins+1]={0,25,50,75,100,125,150,175,200,250,300,350,400,450,
                               500,600,700,800,900,1000,2000,3000,4000,5000};

  stringstream s;
  s << bins[bin] << " < sumJetPt < " << bins[bin+1] << " GeV";

  return s.str();
}

//--------------------------------------------------------------------
*/

void babylooper::ScanChain (TChain* chain, const char* prefix, bool isData, 
                            selectionType mySelectionType, bool makeTemplate, int nEvents){

  int npass = 0;
  selection_    = mySelectionType;
  makeTemplate_ = makeTemplate;

  if     ( selection_ == e_QCDSelection )     cout << "QCD selection" << endl;
  else if( selection_ == e_photonSelection )  cout << "photon selection" << endl;
  else if( selection_ == e_ZSelection )       cout << "Z selection" << endl;
  else { cout << "Unrecognized selection type, quitting" << endl ; exit(0) ; }

  if( makeTemplate_ ) cout << "Making templates" << endl;

  TFile*  metTemplateFile   = new TFile();
  string  metTemplateString = "";
  char*   templateFileName  = "";

  TH1F* reweightHist[3];
  TRandom3 rand;

  if( reweight ){
  
    TFile *reweightFile;
    reweightFile = TFile::Open("reweight.root");
    reweightHist[0] = (TH1F*) reweightFile->Get("ratio_1jet");
    reweightHist[1] = (TH1F*) reweightFile->Get("ratio_2jet");
    reweightHist[2] = (TH1F*) reweightFile->Get("ratio_3jet");
 
//     TCanvas *c1=new TCanvas("c1","",1200,400);
//     c1->Divide(3,1);
//     c1->cd(1);
//     reweightHist[0]->Draw();
//     c1->cd(2);
//     reweightHist[1]->Draw();
//     c1->cd(3);
//     reweightHist[2]->Draw();
  }

  if( !makeTemplate ){

    //select templates
    if( isData ){
      
      if( myTemplateSource == e_QCD ){
        cout << "QCD templates are deprecated. If you want to use this you"
             << "need to make QCD templates using makeTemplates.C." << endl;
        exit(0);
        //templateFileName = Form("output/%s/JetMETTau_templates.root",iter);
        //cout << "Using template file " << templateFileName << endl;
        //metTemplateString = "_JetMETTauTemplate";
        //metTemplateFile = TFile::Open( templateFileName );
      }
      
      else if( myTemplateSource == e_PhotonJet ){
        //templateFileName = Form("output/%s/EG_templates.root",iter);
        templateFileName = Form("output/%s/babylooper_EG_templates.root",iter);
        cout << "Using template file " << templateFileName << endl;
        metTemplateString = "_EGTemplate";
        metTemplateFile = TFile::Open( templateFileName );
      }
      
    }else{
      
      if( myTemplateSource == e_QCD ){
        cout << "QCD templates are deprecated. If you want to use this you"
             << "need to make QCD templates using makeTemplates.C." << endl;
        exit(0);
        //templateFileName = Form("output/%s/QCD_Pt15_templates.root",iter);
        //cout << "Using template file " << templateFileName << endl;
        //metTemplateString = "_QCD_Pt15Template";
        //metTemplateFile = TFile::Open( templateFileName ); 
      }
      
      else if( myTemplateSource == e_PhotonJet ){
        //templateFileName = Form("output/%s/PhotonJet_templates.root",iter);
        templateFileName = Form("output/%s/babylooper_PhotonJet_templates.root",iter);
        cout << "Using template file " << templateFileName << endl;
        metTemplateString = "_PhotonJetTemplate";
        metTemplateFile = TFile::Open( templateFileName );
      }
    }
  }
  
  if( isData ){
  
    ofile_tcmet.open(  Form( "output/%s/babylooper_%s_tcmetprintout.txt" , iter , prefix ) );
    ofile_events.open( Form( "output/%s/babylooper_%s_highmetevents.txt" , iter , prefix ) );
    
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
  

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEvents == -1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  if(debug) cout << "Begin file loop" << endl;

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())){
    
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("T1");
    //cms2.Init(tree);
    setBranches(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();
 
    for (unsigned int event = 0 ; event < nEvents; ++event){
        
      tree->GetEntry(event);
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

      //if( run_ > 140387 ) continue;
      weight_ = 1;

      fillUnderOverFlow( hgenps_pthat  , pthat_ , weight_ );
      fillUnderOverFlow( hphotonpt     , etg_   , weight_ );

      float theMet = -1;
      if     ( myMetType == e_tcmet    ) theMet = tcmet_;
      else if( myMetType == e_tcmetNew ) theMet = tcmetNew_;
      else if( myMetType == e_pfmet    ) theMet = pfmet_;
 
     
      //apply good photon selection-------------------------------------------------------------
  
      if( selection_ == e_photonSelection ) {
        
        if( isData ){
          if( !( HLT_Photon20_Cleaned_L1R_ == 1)  ) continue;
        }else{
          if( !( HLT_Photon20_L1R_ == 1)  )         continue;
        }

        if( nJets_ < 1      ) continue;
        if( jetmax_pt_ < 30 ) continue;

        //Ben selection
        if ( etag_ > 1 )                                 continue;
        if ( etg_ < 22 )                                 continue;
        if ( (1.-r4_) < 0.05 )                           continue;
        if ( hoe_ > 0.1 )                                continue;
        if ( jet_neu_emfrac_ < 0.95 )                    continue; 
        //if ( nvtx_ != 1 )                                continue;

        /*
        //Warren selection
        if ( etag_ > 1 )                                 continue;
        if ( etg_ < 20 )                                 continue;
        if ( (1.-r4_) < 0.05 )                           continue;
        if ( hoe_ > 0.05 )                               continue;
        if ( photon_hcalIso03_ > 2.4 && photon_hcalIso03_ / etg_ > 0.05)        continue;              
        if ( photon_ecalIso03_ > 1.7 && photon_ecalIso03_ / etg_ > 0.05)        continue;       
        if ( photon_ntkIsoSolid03_ > 2 || photon_tkisoSolid03_ / etg_ > 0.1 )   continue;
        */
    
        /*
        //EWK PAS selection
        if( photon_ecalIso04_ > 4.2 * 0.004 * etg_ )       continue;
        if( photon_ecalIso04_ > 4.2 * 0.004 * etg_ )       continue;
        */

        //if ( jet_eta_ > 1 )                              continue;
        //if ( jet_dr_ > 0.5 )                             continue;
        //if ( jet_neu_emfrac_ + jet_chg_emfrac_< 0.95 )   continue; 
        //if ( drel_ < 0.5 )                               continue;
        //if ( ( 1 - jet_neu_emfrac_ ) * etg_ > 1 )        continue;
        //if ( sumJetPt_ < 200. )                          continue;

        //dphixmet_  = deltaPhi( tcmetphi_ , phig_ );
        //metPar_    = theMet * cos( dphixmet_ );
        //metPerp_   = theMet * sin( dphixmet_ );
    
      }
      
      
      //apply Z selection-----------------------------------------------------------------------

      if( selection_ == e_ZSelection ) {

        //if( leptype_ == 2 ) continue;
        if( nJets_ < 1 ) continue;
        
        if( leptype_ == 1 ){
          if( !passm_ll_nom_ == 1) continue;
          if( !passm_lt_nom_ == 1) continue;
        }

        if( reweight ){
          int iJetBin   = getJetBin( nJets_ );
          int bin       = reweightHist[iJetBin]->FindBin( dilpt_ );
          float ratio   = reweightHist[iJetBin]->GetBinContent( bin );
          float randnum = rand.Uniform(0,1);

          //cout << iJetBin << " " << dilpt_ << " " << bin << " " << ratio << " " << randnum << endl;
          if( randnum > ratio ) continue;
          //cout << "pass" << endl;
        }

        //if( nvtx_    != 1 ) continue;
        //if( ptll_ < 20 || ptlt_ < 20 )        continue;
        
        //if( leptype_ == 0 ){
        //  if( !passe_ll_ttbarV2_ == 1 )    continue;
        //  if( !passe_lt_ttbarV2_ == 1 )    continue;
        //}
        //else if( leptype_ == 1 ){
        //  if( !passm_ll_nomttbarV2_ == 1 ) continue;
        //  if( !passm_lt_nomttbarV2_ == 1 ) continue;
        //}
        //else if( leptype_ == 2 ){
        //  if( !( (passe_ll_ttbarV2_ == 1 && passm_lt_nomttbarV2_ == 1) || 
        //         (passe_lt_ttbarV2_ == 1 && passm_ll_nomttbarV2_ == 1) ) ) continue;
        //}

        fillHistos( htcmet            , tcmet_           , weight_ , leptype_ , nJets_ );
        fillHistos( htcmetNew         , tcmetNew_        , weight_ , leptype_ , nJets_ );
        fillHistos( hpfmet            , pfmet_           , weight_ , leptype_ , nJets_  );
        
        if( isData && ( tcmet_ > 30 || pfmet_ > 30 ) ){
          string lepstring[3]={"ee","mm","em"};
          
          
          ofile_events << "|" << setw(8)  << run_                        << setw(4) 
                       << "|" << setw(6)  << lumi_                       << setw(4) 
                       << "|" << setw(12) << event_                      << setw(4) 
                       << "|" << setw(6)  << lepstring[leptype_]         << setw(4) 
                       << "|" << setw(6)  << nJets_                      << setw(4) 
                       << "|" << setw(6)  << nbtags_                     << setw(4) 
                       << "|" << setw(8)  << fround(tcmet_,1)            << setw(4) 
                       << "|" << setw(8)  << fround(pfmet_,1)            << setw(4) 
                       << "|" << setw(8)  << fround(dphijetmet_,2)       << setw(4) << "|" << endl; 
          
        }
      }
      
      
      //fill histos and ntuple----------------------------------------------------------- 
      
          


      if( myTemplateSource == e_QCD ){
        if( nJets_ < 2      ) continue;
        if( jetmax_pt_ < 50 ) continue;
      }
      else if ( myTemplateSource == e_PhotonJet ){
        if( nJets_ < 1      ) continue;
        if( jetmax_pt_ < 30 ) continue;
      }
      

      npass++;

      float pthad = -1;
      if( selection_ == e_photonSelection ) pthad = etg_;
      if( selection_ == e_ZSelection )      pthad = dilpt_;

      //fill met template-----------------------------------------------------------------------
            

      if( makeTemplate_ ) {

        int iJetBin          = getJetBin( nJets_ );
        int iSumJetPtBin     = getSumJetPtBin( sumJetPt_ );
        int iBosonPtBin      = getBosonPtBin( etg_ );
        
        //fill templates binned by njets, sumjetpt, boson pt        
        fillUnderOverFlow( tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  tcmet_    , weight_ );
        fillUnderOverFlow( pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  pfmet_    , weight_ );
        fillUnderOverFlow( tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] ,  tcmetNew_ , weight_ );

        //fill templates binned by njets, sumjetpt
        fillUnderOverFlow( tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , weight_ );
        fillUnderOverFlow( pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , weight_ );
        fillUnderOverFlow( tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] ,  tcmetNew_ , weight_ );
        
      }

          

      int iJ = nJets_;
      if( iJ > 4 ) iJ = 4;
       
      fillUnderOverFlow( hpthad[iJ] ,  pthad , weight_ );
      fillUnderOverFlow( hpthad[0]  ,  pthad , weight_ );

      //fill predicted and observed met histos--------------------------------------------------

      if( !makeTemplate_ ){
        
        int iJetBin      = getJetBin( nJets_ );
        int iBosonPtBin  = getBosonPtBin( pthad );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
                
        //TH1F* hmet = getMetTemplate( metTemplateFile , iJetBin , iSumJetPtBin , iBosonPtBin , weight_ );

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
        
        //int bin = hmet->FindBin(60);
        //float int60 = hmet->Integral(bin,100000);
        //if( int60 > 0 ) cout << "Found event " << int60 * hmet->GetEntries() << endl;
        
        if( selection_ == e_ZSelection ){
          if( leptype_ == 2 ){
            fillUnderOverFlow( metObserved_df , tcmet_ , weight_  );
            metPredicted_df->Add( hmet );
            continue; //exclude emu final state
          }
        }
        
        fillUnderOverFlow( metObserved_njets[iJetBin] , theMet , weight_ );
        metPredicted_njets[iJetBin]->Add( hmet );
        
        fillUnderOverFlow( metObserved , theMet , weight_ );
        metPredicted->Add( hmet );
        

        if( selection_ == e_ZSelection ) {
                
          if( leptype_ == 0 ){
            fillUnderOverFlow( metObserved_ee , theMet , weight_ );
            metPredicted_ee->Add( hmet );
          }
          if( leptype_ == 1 ){
            fillUnderOverFlow( metObserved_mm , theMet , weight_ );
            metPredicted_mm->Add( hmet );
          }
        
          if( leptype_ == 0 || leptype_ == 1 ){
            fillUnderOverFlow( metObserved_sf , tcmet_ , weight_  );
            metPredicted_sf->Add( hmet );
          }else if( leptype_ == 2 ){
            fillUnderOverFlow( metObserved_df , tcmet_ , weight_  );
            metPredicted_df->Add( hmet );
          }else{
            cout << "UNRECOGNIZED LEPTYPE " << leptype_ << endl;
          }
        }
                        
        if( selection_ == e_photonSelection ||  selection_ == e_ZSelection){

          if( pthad < 40 ){
            fillUnderOverFlow( metObserved_ptlt40 , theMet , weight_ );
            metPredicted_ptlt40->Add( hmet );
          }
          if( pthad > 40 && pthad < 60 ){
            fillUnderOverFlow( metObserved_pt40_60 , theMet , weight_ );
            metPredicted_pt40_60->Add( hmet );
          }
          if( pthad > 60 ){
            fillUnderOverFlow( metObserved_ptgt60 , theMet , weight_ );
            metPredicted_ptgt60->Add( hmet );
          }
          if( pthad < 50 ){
            fillUnderOverFlow( metObserved_ptlt50 , theMet , weight_ );
            metPredicted_ptlt50->Add( hmet );
          }else{
            fillUnderOverFlow( metObserved_ptgt50 , theMet , weight_ );
            metPredicted_ptgt50->Add( hmet );
          }
        }
      
        delete hmet;
        
      }
    } // end loop over events
  } // end loop over files
  
  cout << npass << " events passing selection" << endl;
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  
  //normalize met templates
  if( makeTemplate_ ) {
    
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        for( int iBosonPtBin = 0 ; iBosonPtBin < nBosonPtBins ; iBosonPtBin++ ){
          
          float scale = tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
          if( scale > 0 )  tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
          
          scale = tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
          if( scale > 0 )  tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
          
          scale = pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Integral();
          if( scale > 0 )  pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] -> Scale ( 1. / scale );
          
        }
      }
    }

    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        
        float scale = tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
        
        scale = tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
        
        scale = pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Integral();
        if( scale > 0 )  pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ] -> Scale ( 1. / scale );
        
        
      }
    }
  }

  // make histos rootfile
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  if( makeTemplate ) saveHist( Form("output/%s/babylooper_%s_templates.root" , iter , prefix ) );
  else               saveHist( Form("output/%s/babylooper_%s%s.root"         , iter , prefix , metTemplateString.c_str() ) );
  deleteHistos();
  
} // end ScanChain

TH1F* babylooper::getMetTemplate( TFile* file, int iJetBin , int iSumJetPtBin , int iBosonPtBin , float weight ){

  char* metstring = "";

  if     ( myMetType == e_tcmet    ) metstring = "tcmet";
  else if( myMetType == e_tcmetNew ) metstring = "tcmetNew";
  else if( myMetType == e_pfmet    ) metstring = "pfmet";
    
  TH1F* hmet = new TH1F();

  if( useCombinedTemplates ){
    hmet     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin));
  }else{
    hmet     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin,iBosonPtBin));
  }
    
  //if there are at least nMetEntries in the template, return template
  if( hmet->GetEntries() >= nMetEntries ){
    hmet->Scale( weight );
    return hmet;
  }
  
  if( hmet->GetEntries() > 0 )
    hmet->Scale( hmet->GetEntries() / hmet->Integral() );
  int counter = 1;
 
  //cout << endl << "Found template with " << hmet->GetEntries() << " entries" << endl;
  //cout << iJetBin << " " << iSumJetPtBin << endl;

  //add histograms in sumjetpt bins adjacent to current bin to get enough entries
  while( hmet->GetEntries() < nMetEntries && counter < nSumJetPtBins ){

    //add template 1 bin lower in sumJetPt, if it exists
    if( iSumJetPtBin - counter >= 0 ){
     
      TH1F *hmetLo = new TH1F();

      if( useCombinedTemplates ){
        hmetLo     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin-counter));
      }else{
        hmetLo     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin-counter,iBosonPtBin));
      }
 
      //cout << "Adding lo entries " << hmetLo->GetEntries() << endl;
      if( hmetLo->GetEntries() > 0 ){
        hmetLo->Scale( hmetLo->GetEntries() / hmetLo->Integral() );
        hmet->Add(hmetLo);
      }

      delete hmetLo;
    }
     
    //add template 1 bin higher in sumJetPt, if it exists
    if( iSumJetPtBin + counter < nSumJetPtBins ){
      
      TH1F *hmetHi = new TH1F();
 
      if( useCombinedTemplates ){
        hmetHi     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin+counter));
      }else{
        hmetHi     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin+counter,iBosonPtBin));
      }
      
      //cout << "Adding hi entries " << hmetHi->GetEntries() << endl;
      if( hmetHi->GetEntries() > 0 ){
        hmetHi->Scale( hmetHi->GetEntries() / hmetHi->Integral() );
        hmet->Add(hmetHi);
      }
 
      delete hmetHi;
    }

    counter++;
  }
  
  if( hmet->GetEntries() > 0 )
  hmet->Scale( weight / hmet->Integral() );
  return hmet;
}

float babylooper::deltaPhi( float phi1 , float phi2){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}


void babylooper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);

}

void babylooper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  char* leptype[4]   = {"ee", "mm", "em", "all"};
  char* jetbin[4]    = {"0j", "1j", "geq2j", "allj"};

  char* leptype_title[4]   = {"ee", "#mu#mu", "e#mu", "all leptons"};
  char* jetbin_title[4]    = {"0 jets", "1 jet", "#geq 2 jets", "all jets"};

  for (int i = 0; i < 4; i++) {
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

  hr4 = new TH1F("hr4","hr4",200,0,2);

  char* pttitle[5]={"all jets","1 jet","2 jet","3 jet","#geq 4 jet"};
  
  for( int i = 0 ; i < 5 ; ++i ){
    hpthad[i] = new TH1F(Form("pthad_%i",i),pttitle[i],200,0,200);
    hpthad[i]->GetXaxis()->SetTitle("p_{T} (GeV)");
  }

  hgenps_pthat = new TH1F("hgenps_pthat","",100,0,100);
  hphotonpt    = new TH1F("hphotonpt","",100,0,100);

  if( !makeTemplate_ ){
    
    metObserved  = new TH1F("metObserved", "Observed MET",500,0,500);
    metPredicted = new TH1F("metPredicted","Predicted MET",500,0,500);
    
    metObserved_sf  = new TH1F("metObserved_sf", "Observed MET (SF)",500,0,500);
    metPredicted_sf = new TH1F("metPredicted_sf","Predicted MET (SF)",500,0,500);
    
    metObserved_df  = new TH1F("metObserved_df", "Observed MET (DF)",500,0,500);
    metPredicted_df = new TH1F("metPredicted_df","Predicted MET (DF)",500,0,500);
    
    metObserved_ptlt40     = new TH1F("metObserved_ptlt40", "Observed MET (p_{T}<40 GeV)",500,0,500);
    metPredicted_ptlt40    = new TH1F("metPredicted_ptlt40","Predicted MET  (p_{T}<40 GeV)",500,0,500);
    
    metObserved_ptlt50     = new TH1F("metObserved_ptlt50", "Observed MET (p_{T}<50 GeV)",500,0,500);
    metPredicted_ptlt50    = new TH1F("metPredicted_ptlt50","Predicted MET  (p_{T}<50 GeV)",500,0,500);
    
    metObserved_ptgt50     = new TH1F("metObserved_ptgt50", "Observed MET (p_{T}>50 GeV)",500,0,500);
    metPredicted_ptgt50    = new TH1F("metPredicted_ptgt50","Predicted MET  (p_{T}>50 GeV)",500,0,500);
    
    metObserved_pt40_60    = new TH1F("metObserved_pt40_60", "Observed MET (40<p_{T}<60 GeV)",500,0,500);
    metPredicted_pt40_60   = new TH1F("metPredicted_pt40_60","Predicted MET  (40<p_{T}<60 GeV)",500,0,500);
    
    metObserved_ptgt60     = new TH1F("metObserved_ptgt60", "Observed MET (p_{T}>60 GeV)",500,0,500);
    metPredicted_ptgt60    = new TH1F("metPredicted_ptgt60","Predicted MET  (p_{T}>60 GeV)",500,0,500);
    
    metObserved_ee  = new TH1F("metObserved_ee", "Observed MET (ee)",500,0,500);
    metPredicted_ee = new TH1F("metPredicted_ee","Predicted MET (ee)",500,0,500);
    
    metObserved_mm  = new TH1F("metObserved_mm", "Observed MET (#mu#mu)",500,0,500);
    metPredicted_mm = new TH1F("metPredicted_mm","Predicted MET (#mu#mu)",500,0,500);
    
    metParObserved  = new TH1F("metParObserved", "Observed MET (Parallel)",1000,-500,500);
    metParPredicted = new TH1F("metParPredicted","Predicted MET (Parallel)",1000,-500,500);
    
    metPerpObserved  = new TH1F("metPerpObserved", "Observed MET (Perpendicular)",500,0,500);
    metPerpPredicted = new TH1F("metPerpPredicted","Predicted MET (Perpendicular)",500,0,500);
    
    metObserved->Sumw2();
    metPredicted->Sumw2();
    
    metObserved_ptlt40->Sumw2();
    metPredicted_ptlt40->Sumw2();
    
    metObserved_pt40_60->Sumw2();
    metPredicted_pt40_60->Sumw2();
    
    metObserved_ptgt60->Sumw2();
    metPredicted_ptgt60->Sumw2();
    
    metObserved_ptgt50->Sumw2();
    metPredicted_ptgt50->Sumw2();
    
    metObserved_ptlt50->Sumw2();
    metPredicted_ptlt50->Sumw2();
    
    metObserved_ee->Sumw2();
    metPredicted_ee->Sumw2();
    
    metObserved_mm->Sumw2();
    metPredicted_mm->Sumw2();
    
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
  }
   
  if( makeTemplate_ ){

    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        for( int iBosonPtBin = 0 ; iBosonPtBin < nBosonPtBins ; iBosonPtBin++ ){

          tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("tcmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                       Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                            bosonPtString(iBosonPtBin).c_str()),500,0,500);

          tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("tcmetNewTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                          Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                               bosonPtString(iBosonPtBin).c_str()),500,0,500);

          pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("pfmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                       Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                            bosonPtString(iBosonPtBin).c_str()),500,0,500);
        
          tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();
          tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();
          pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->Sumw2();

          tcmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
          tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("tcmetNew (GeV)");
          pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          

        }
      }
    }

  
   for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      
  
        tcmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                 Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
        
        tcmetNewTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetNewTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                    Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
        
        pfmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("pfmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                 Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),500,0,500);
        
        tcmetTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
        tcmetNewTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
        pfmetTemplate_combined[iJetBin][iSumJetPtBin]->Sumw2();
         
        tcmetTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmet (GeV)");
        tcmetNewTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("tcmetNew (GeV)");
        pfmetTemplate_combined[iJetBin][iSumJetPtBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
         
      }
    }
  
  }
}


void babylooper::setBranches (TTree* tree){
  tree->SetBranchAddress("run",          &run_          );
  tree->SetBranchAddress("weight",       &weight_       );
  tree->SetBranchAddress("pthat",        &pthat_        );
  tree->SetBranchAddress("lumi",         &lumi_         );
  tree->SetBranchAddress("event",        &event_        );

  tree->SetBranchAddress("nvtx",         &nvtx_         );
  tree->SetBranchAddress("pfmet",        &pfmet_        );
  tree->SetBranchAddress("pfmetphi",     &pfmetphi_     );
  tree->SetBranchAddress("pfsumet",      &pfsumet_      );
  tree->SetBranchAddress("met",          &met_          );
  tree->SetBranchAddress("metphi",       &metphi_       );
  tree->SetBranchAddress("sumet",        &sumet_        );
  tree->SetBranchAddress("mumet",        &mumet_        );
  tree->SetBranchAddress("mumetphi",     &mumetphi_     );
  tree->SetBranchAddress("musumet",      &musumet_      );
  tree->SetBranchAddress("mujesmet",     &mujesmet_     );
  tree->SetBranchAddress("mujesmetphi",  &mujesmetphi_  );
  tree->SetBranchAddress("mujessumet",   &mujessumet_   );
  tree->SetBranchAddress("genmet",       &genmet_       );
  tree->SetBranchAddress("genmetphi",    &genmetphi_    );
  tree->SetBranchAddress("gensumet",     &gensumet_     );
  tree->SetBranchAddress("tcmet",        &tcmet_        );
  tree->SetBranchAddress("tcmetphi",     &tcmetphi_     );     
  tree->SetBranchAddress("tcsumet",      &tcsumet_      );  
  tree->SetBranchAddress("tcmetNew",     &tcmetNew_        );
  tree->SetBranchAddress("tcmetphiNew",  &tcmetphiNew_     );     
  tree->SetBranchAddress("tcsumetNew",   &tcsumetNew_      );  
  tree->SetBranchAddress("dphixmet",     &dphixmet_     );
  tree->SetBranchAddress("metpar",       &metPar_       );
  tree->SetBranchAddress("metperp",      &metPerp_      );
  tree->SetBranchAddress("njets",        &nJets_        );     
  tree->SetBranchAddress("sumjetpt",     &sumJetPt_     );     
  tree->SetBranchAddress("vecjetpt",     &vecJetPt_     );     
  tree->SetBranchAddress("nbtags",       &nbtags_       );
  tree->SetBranchAddress("ndphijetmet",  &dphijetmet_   );

  if( selection_ == e_photonSelection ){

    //photon stuff
    tree->SetBranchAddress("ng",      &nPhotons_          ); 
    tree->SetBranchAddress("etg",     &etg_               ); 
    tree->SetBranchAddress("etag",    &etag_              ); 
    tree->SetBranchAddress("phig",    &phig_              ); 
    tree->SetBranchAddress("hoe",     &hoe_               ); 
    tree->SetBranchAddress("eciso",   &eciso_             ); 
    tree->SetBranchAddress("hciso",   &hciso_             ); 
    tree->SetBranchAddress("tkiso",   &tkiso_             );   
    tree->SetBranchAddress("swiss",   &swiss_             ); 
    tree->SetBranchAddress("seed",    &seed_              ); 
    tree->SetBranchAddress("s4",      &s4_                ); 
    tree->SetBranchAddress("r4",      &r4_                ); 
    
    //photon-matched jet stuff
    tree->SetBranchAddress("jetdr",                 &jet_dr_                );           
    tree->SetBranchAddress("jetpt",                 &jet_pt_                );   
    tree->SetBranchAddress("jeteta",                &jet_eta_               ); 
    tree->SetBranchAddress("jetenergy",             &jet_energy_            ); 
    tree->SetBranchAddress("jetchargedemfrac",      &jet_chg_emfrac_        ); 
    tree->SetBranchAddress("jetchargedhadfrac",     &jet_chg_hadfrac_       ); 
    tree->SetBranchAddress("jetneutralemfrac",      &jet_neu_emfrac_        ); 
    tree->SetBranchAddress("jetneutralhadfrac",     &jet_neu_hadfrac_       ); 
    //tree->SetBranchAddress("jetncharged",           &jet_nchg_              ); 
    //tree->SetBranchAddress("jetnmuon",              &jet_nmuon_             ); 
    //tree->SetBranchAddress("jetnneutral",           &jet_nneu_              ); 
    tree->SetBranchAddress("jetdphimet",            &jet_dphimet_           ); 
    tree->SetBranchAddress("jetdpt",                &jet_dpt_               ); 
    //tree->SetBranchAddress("jetdrgen",              &jet_drgen_             ); 
    //tree->SetBranchAddress("drel",                  &drel_                  ); 
    
    tree->SetBranchAddress("photon_pixelseed",         &photon_pixelseed_);       
    tree->SetBranchAddress("photon_e15",               &photon_e15_);             
    tree->SetBranchAddress("photon_e25max",            &photon_e25max_);          
    tree->SetBranchAddress("photon_e33",               &photon_e33_);             
    tree->SetBranchAddress("photon_e55",               &photon_e55_);             
    tree->SetBranchAddress("photon_ecalIso03",         &photon_ecalIso03_);       
    tree->SetBranchAddress("photon_ecalIso04",         &photon_ecalIso04_);       
    tree->SetBranchAddress("photon_hcalIso03",         &photon_hcalIso03_);       
    tree->SetBranchAddress("photon_hcalIso04",         &photon_hcalIso04_);       
    tree->SetBranchAddress("photon_ntkIsoHollow03",    &photon_ntkIsoHollow03_);  
    tree->SetBranchAddress("photon_ntkIsoHollow04",    &photon_ntkIsoHollow04_);  
    tree->SetBranchAddress("photon_ntkIsoSolid03",     &photon_ntkIsoSolid03_);   
    tree->SetBranchAddress("photon_ntkIsoSolid04",     &photon_ntkIsoSolid04_);   
    tree->SetBranchAddress("photon_sigmaEtaEta",       &photon_sigmaEtaEta_);     
    tree->SetBranchAddress("photon_sigmaIEtaIEta",     &photon_sigmaIEtaIEta_);   
    tree->SetBranchAddress("photon_tkisoHollow03",     &photon_tkisoHollow03_);   
    tree->SetBranchAddress("photon_tkisoHollow04",     &photon_tkisoHollow04_);   
    tree->SetBranchAddress("photon_tkisoSolid03",      &photon_tkisoSolid03_);    
    tree->SetBranchAddress("photon_tkisoSolid04",      &photon_tkisoSolid04_);     
  }
  
  tree->SetBranchAddress("maxjetpt",              &jetmax_pt_             ); 
  tree->SetBranchAddress("maxjetdphimet",         &jetmax_dphimet_        ); 
  
  //trigger
  tree->SetBranchAddress("HLT_Jet15U",            &HLT_Jet15U_            ); 
  tree->SetBranchAddress("HLT_Jet30U",            &HLT_Jet30U_            ); 

  tree->SetBranchAddress("HLT_Photon10_L1R",      &HLT_Photon10_L1R_      ); 
  tree->SetBranchAddress("HLT_Photon15_L1R",      &HLT_Photon15_L1R_      ); 
  tree->SetBranchAddress("HLT_Photon20_L1R",      &HLT_Photon20_L1R_      ); 
  tree->SetBranchAddress("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_  ); 
  tree->SetBranchAddress("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_  ); 
  tree->SetBranchAddress("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_  ); 

  if( selection_ == e_ZSelection ){
    //Z stuff
    //tree->SetBranchAddress("passz",                  &passz_                  );
    //tree->SetBranchAddress("pdgid",                  &pdgid_                  );
    tree->SetBranchAddress("passm_ll_nom",           &passm_ll_nom_           );
    tree->SetBranchAddress("passm_ll_nomttbar",      &passm_ll_nomttbar_      );  
    tree->SetBranchAddress("passm_ll_nomttbarV2",    &passm_ll_nomttbarV2_    );  
    tree->SetBranchAddress("passe_ll_ttbar",         &passe_ll_ttbar_         );
    tree->SetBranchAddress("passe_ll_ttbarV1",       &passe_ll_ttbarV1_       );
    tree->SetBranchAddress("passe_ll_ttbarV2",       &passe_ll_ttbarV2_       );
    tree->SetBranchAddress("passe_ll_cand01",        &passe_ll_cand01_        );  
    tree->SetBranchAddress("passm_lt_nom",           &passm_lt_nom_           );
    tree->SetBranchAddress("passm_lt_nomttbar",      &passm_lt_nomttbar_      );  
    tree->SetBranchAddress("passm_lt_nomttbarV2",    &passm_lt_nomttbarV2_    );  
    tree->SetBranchAddress("passe_lt_ttbar",         &passe_lt_ttbar_         );
    tree->SetBranchAddress("passe_lt_ttbarV1",       &passe_lt_ttbarV1_       );
    tree->SetBranchAddress("passe_lt_ttbarV2",       &passe_lt_ttbarV2_       );
    tree->SetBranchAddress("passe_lt_cand01",        &passe_lt_cand01_        );
    //tree->SetBranchAddress("TMLastStationTight_ll",  &TMLastStationTight_ll_  );
    //tree->SetBranchAddress("TMLastStationTight_lt",  &TMLastStationTight_lt_  );
  
    tree->SetBranchAddress("ptll",                 &ptll_                 );  
    tree->SetBranchAddress("ptlt",                 &ptlt_                 ); 
    tree->SetBranchAddress("dilmass",              &dilmass_              ); 
    tree->SetBranchAddress("dilpt",                &dilpt_                ); 
    tree->SetBranchAddress("flagll",               &flagll_               );  
    tree->SetBranchAddress("flaglt",               &flaglt_               ); 
    tree->SetBranchAddress("leptype",              &leptype_              );


  }

}

//--------------------------------------------------------------------

void babylooper::fillHistos(TH1F *h1[4],float value, float weight, int myType)
{

  fillUnderOverFlow(h1[myType], value, weight);      
  fillUnderOverFlow(h1[3],      value, weight);      
}

//--------------------------------------------------------------------

void babylooper::fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{

  if( nJetsIdx > 2 ) nJetsIdx = 2;

  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      
}

//--------------------------------------------------------------------


