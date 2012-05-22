#include "babylooper.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "../CORE/metTemplatesSelections.cc"
#include "TChain.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

enum metType   { e_tcmet = 0, e_tcmetNew = 1, e_pfmet = 2};
enum templateSource { e_QCD = 0, e_PhotonJet = 1, e_PhotonJetStitched };
enum templateType { e_njets_ht = 0, e_njets_ht_nvtx = 1, e_njets_ht_vecjetpt };

//------------------------------------------------
//USER PARAMS
//------------------------------------------------
bool           debug              = false;
bool           reweight           = false;   //reweight for photon vs. Z pt
bool           doVtxReweight      = false;   //reweight templates for nVertices
bool           setTemplateErrors  = true; 
metType        myMetType          = e_pfmet;
//templateSource myTemplateSource   = e_PhotonJet;
templateSource myTemplateSource   = e_PhotonJetStitched;
//templateSource myTemplateSource   = e_QCD;
templateType   myTemplateType     = e_njets_ht;
//templateType   myTemplateType     = e_njets_ht_nvtx;
//templateType   myTemplateType     = e_njets_ht_vecjetpt;
char*          iter               = "njetsgeq3";

//------------------------------------------------



using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void babylooper::ScanChain (TChain* chain, const char* Z_version, const char* template_version, const char* prefix, 
			    bool isData, selectionType mySelectionType, bool makeTemplate, int nEvents){

  //cout << "Setting min entries to 50!" << endl;
  
  int npass = 0;
  selection_    = mySelectionType;
  makeTemplate_ = makeTemplate;

  if( !isData && myTemplateSource == e_PhotonJetStitched ){
    myTemplateSource = e_PhotonJet;
    cout << "Switching MC template to PhotonJet" << endl;
  }

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
  TH1F* vtxReweightHist;

  if( doVtxReweight && makeTemplate_ ){
    TFile *vtxReweightFile = TFile::Open("vtxFile.root");
    vtxReweightHist = (TH1F*) vtxReweightFile->Get("vtxReweightHist");
    cout << "Using nVertex reweighting" << endl;
    cout << "weight (1 vtx) " << vtxReweightHist->GetBinContent(2) << endl;
    cout << "weight (2 vtx) " << vtxReweightHist->GetBinContent(3) << endl;
    cout << "weight (3 vtx) " << vtxReweightHist->GetBinContent(4) << endl;
    cout << "weight (4 vtx) " << vtxReweightHist->GetBinContent(5) << endl;
    cout << "weight (5 vtx) " << vtxReweightHist->GetBinContent(6) << endl;
  }

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
	cout << "NO QCD TEMPLATES" << endl;
	exit(0);
        templateFileName = Form("qcd-mini-ntuple/victor_templates_overlap.root");
        cout << "Using template file " << templateFileName << endl;
        metTemplateString = "_victorTemplate";
        metTemplateFile = TFile::Open( templateFileName );
      }
      
      else if( myTemplateSource == e_PhotonJet ){
	templateFileName = Form("../templates/%s/photon_templates.root",template_version);
        cout << "Using template file " << templateFileName << endl;
        metTemplateString = "_EGTemplate";
        metTemplateFile = TFile::Open( templateFileName );
      }

      else if( myTemplateSource == e_PhotonJetStitched ){
        templateFileName = Form("../templates/%s/photon_templates.root",template_version);
        cout << "Using template file " << templateFileName << endl;
        metTemplateString = "_PhotonStitchedTemplate";
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
      
      else if( myTemplateSource == e_PhotonJet || myTemplateSource == e_PhotonJetStitched ){
	cout << "WARNING! USING DATA TEMPLATES FOR MC" << endl;
	templateFileName = Form("../templates/%s/photon_templates.root",template_version);
	//templateFileName = "/tas03/home/benhoob/metTemplate/output/nov5th/babylooper_PhotonJet_templates.root";
	cout << "Using template file " << templateFileName << endl;
        metTemplateString = "_PhotonJetTemplate";
        metTemplateFile = TFile::Open( templateFileName );
      }
    }
  }
  
  if( isData ){
  
    ofile_tcmet.open(  Form( "../output/%s/babylooper_%s_tcmetprintout.txt" , Z_version , prefix ) );
    ofile_events.open( Form( "../output/%s/babylooper_%s_highmetevents.txt" , Z_version , prefix ) );
    
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


  //set template counters to zero
  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
    
      n_metPredicted[iJetBin][iSumJetPtBin] = 0;
      n_metPredicted_ee[iJetBin][iSumJetPtBin] = 0;
      n_metPredicted_mm[iJetBin][iSumJetPtBin] = 0;
      
      for( int i = 0 ; i < 4 ; ++i ){

        nphoton_metPredicted[i][iJetBin][iSumJetPtBin] = 0;
        nphoton_metPredicted_ee[i][iJetBin][iSumJetPtBin] = 0;
        nphoton_metPredicted_mm[i][iJetBin][iSumJetPtBin] = 0;
        
        nqcd_metPredicted[i][iJetBin][iSumJetPtBin] = 0;
        nqcd_metPredicted_ee[i][iJetBin][iSumJetPtBin] = 0;
        nqcd_metPredicted_mm[i][iJetBin][iSumJetPtBin] = 0;
     
      }
    }
  }

  ofstream eefile(Form("../output/%s/%s_ee.txt",Z_version,prefix));
  ofstream mmfile(Form("../output/%s/%s_mm.txt",Z_version,prefix));

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

      //weight_ = 1;
      float mcweight = 1;
      if( !isData ) mcweight = weight_ * davtxweight_;

      fillUnderOverFlow( hgenps_pthat  , pthat_ , mcweight );
      fillUnderOverFlow( hphotonpt     , etg_   , mcweight );

      float theMet = -1;
      if     ( myMetType == e_tcmet    ) theMet = tcmet_;
      else if( myMetType == e_tcmetNew ) theMet = tcmetNew_;
      else if( myMetType == e_pfmet    ) theMet = pfmet_;
      //theMet = pfmetcor_;
     
      //apply good photon selection-------------------------------------------------------------
  
      if( selection_ == e_photonSelection ) {
        
        //-------------------------------------------------------
        //trigger + jet selection
        //-------------------------------------------------------
        if( isData ){
          if( !( HLT_Photon20_Cleaned_L1R_ == 1)  ) continue;
        }else{
          if( !( HLT_Photon20_L1R_ == 1)  )         continue;
        }

        if( nJets_ < 2      ) continue;
        if( jetmax_pt_ < 30 ) continue;
        //if( nvtx_ < 3 )       continue;

        //-------------------------------------------------------
        //selection applied in addition to isGoodEMObject
        //-------------------------------------------------------
        if( jet_pt_ - etg_ < -5 )                                 continue; //cleaning
        if( fabs( etag_ ) > 2 )                                   continue; //photon eta < 2
        if( jet_neu_emfrac_ < 0.7 )                               continue; //0.8 < emf < 0.95
        if( pfjetid_ != 1 )                                       continue; //pass PFJetID
        if( isData && run_ > 144114 )                             continue; //restrict run range to Run2010A
        
      }
      
      //apply Z selection-----------------------------------------------------------------------

      if( selection_ == e_ZSelection ) {

        if( leptype_ == 0 ) eefile << dilmasscor_ << endl;
        if( leptype_ == 1 ) mmfile << dilmasscor_ << endl;
     
        if( dilmass_ > 76. && dilmass_ < 106. ){

          if( nJets_ == 0 ){
            hyield_0j->Fill(0.5,          mcweight);
            hyield_0j->Fill(1.5+leptype_, mcweight);
          }
          if( nJets_ == 1 ){
            hyield_1j->Fill(0.5,          mcweight);
            hyield_1j->Fill(1.5+leptype_, mcweight);
          }
          if( nJets_ == 2 ){
            hyield_2j->Fill(0.5,          mcweight);
            hyield_2j->Fill(1.5+leptype_, mcweight);
          }
          if( nJets_ > 2 ){
            hyield_g2j->Fill(0.5,          mcweight);
            hyield_g2j->Fill(1.5+leptype_, mcweight);
          }
        }


        //------------------------------------------------------------
        //event selection
        //------------------------------------------------------------

	if( leptype_ == 0 ){
	  if( jetpt_ll_ - ptll_ < -5  ) continue; 
	  if( jetpt_lt_ - ptlt_ < -5  ) continue; 
	}
	if( leptype_ == 1 && npfmuons_ < 2 ) continue;
	if( leptype_ == 2 && npfmuons_ < 1 ) continue;

	//jetselection
        if( nJets_ < 3 )                                                     continue; //>=2 jets

	//------------------------------
	//fill histos before Z-veto
	//------------------------------

        fillHistos( hdilmass , dilmass_ , mcweight , leptype_ , nJets_ );
	if( pfmet_ > 60. ){
	  fillHistos( hdilmass_pfmet60 , dilmass_ , mcweight , leptype_ , nJets_ );
	}
	if( leptype_ == 2) fillUnderOverFlow( metObserved_df_nozveto , theMet , mcweight  );

	//Zmass
	if( leptype_ == 1 ){
	  if( dilmasspf_ < 81. || dilmasspf_ > 101. )                          continue; 
	}else{
	  if( dilmass_ < 81. || dilmass_ > 101. )                              continue; 
	}
 
	hresponse->Fill( genmet_ , pfmet_ / genmet_ );
	hgenmet_all->Fill( genmet_ );
	if( pfmet_ > 60 ) hgenmet_pass->Fill( genmet_ );

        //------------------------------------------------------------
        //reweighting procedure for photon vs. Z pt
        //------------------------------------------------------------
        if( reweight ){
          int iJetBin   = getJetBin( nJets_ );
          int bin       = reweightHist[iJetBin]->FindBin( dilpt_ );
          float ratio   = reweightHist[iJetBin]->GetBinContent( bin );
          float randnum = rand.Uniform(0,1);

          //cout << iJetBin << " " << dilpt_ << " " << bin << " " << ratio << " " << randnum << endl;
          if( randnum > ratio ) continue;
          //cout << "pass" << endl;
        }
            
        fillHistos( htcmet            , tcmet_           , mcweight , leptype_ , nJets_ );
        fillHistos( htcmetNew         , tcmetNew_        , mcweight , leptype_ , nJets_ );
        fillHistos( hpfmet            , pfmet_           , mcweight , leptype_ , nJets_  );
        
        if( isData && ( pfmet_ > 60 ) ){
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
      npass++;

      hyield->Fill(0.5,          mcweight);
      hyield->Fill(1.5+leptype_, mcweight);

      if( pfmet_ > 30. ){
        hyield_pfmet30->Fill(0.5,          mcweight);
        hyield_pfmet30->Fill(1.5+leptype_, mcweight);
      }
      if( pfmet_ > 60. ){
        hyield_pfmet60->Fill(0.5,          mcweight);
        hyield_pfmet60->Fill(1.5+leptype_, mcweight);
      }
      if( pfmet_ > 120. ){
        hyield_pfmet120->Fill(0.5,          mcweight);
        hyield_pfmet120->Fill(1.5+leptype_, mcweight);
      }

      hnVtx->Fill( nvtx_ , mcweight );
      hvecJetPt->Fill( vecJetPt_ , mcweight );

      float pthad = -1;
      if( selection_ == e_photonSelection ) pthad = etg_;
      if( selection_ == e_ZSelection )      pthad = dilpt_;

      //fill met template-----------------------------------------------------------------------

      if( makeTemplate_ ) {

        int iJetBin          = getJetBin( nJets_ );
        int iSumJetPtBin     = getSumJetPtBin( sumJetPt_ );
        int iBosonPtBin      = getBosonPtBin( etg_ );
        int iVtxBin          = getVtxBin( nvtx_ );

        float templateWeight = mcweight;
        if( doVtxReweight ){
          int vtxBin = nvtx_ + 1;
          if( vtxBin > 6 ) vtxBin = 6;
          templateWeight *= vtxReweightHist->GetBinContent( vtxBin );
          //cout << "nVtx " << nvtx_ << " weight " << vtxReweightHist->GetBinContent( vtxBin ) << endl;
        }

        //fill templates binned by njets, sumjetpt, boson pt        
        fillUnderOverFlow( tcmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ]    ,  pfmet_    , templateWeight );
        fillUnderOverFlow( tcmetNewTemplate[ iJetBin ][ iSumJetPtBin ][ iBosonPtBin ] ,  tcmetNew_ , templateWeight );

        //fill templates binned by njets, sumjetpt, nVtx
        fillUnderOverFlow( tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ]    ,  pfmet_    , templateWeight );
        fillUnderOverFlow( tcmetNewTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] ,  tcmetNew_ , templateWeight );

        //fill templates binned by njets, sumjetpt
        fillUnderOverFlow( tcmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_combined[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_    , templateWeight );
        fillUnderOverFlow( tcmetNewTemplate_combined[ iJetBin ][ iSumJetPtBin ] ,  tcmetNew_ , templateWeight );

        /*
        float tcmet_ee    = smearMet( tcmet_    , tcmetPhi_    , "ee" );
        float tcmetnew_ee = smearMet( tcmetNew_ , tcmetPhiNew_ , "ee" );
        float pfmet_ee    = smearMet( pfmet_    , pfmetPhi_    , "ee" );

        float tcmet_mm    = smearMet( tcmet_    , tcmetPhi_    , "mm" );
        float tcmetnew_mm = smearMet( tcmetNew_ , tcmetPhiNew_ , "mm" );
        float pfmet_mm    = smearMet( pfmet_    , pfmetPhi_    , "mm" );

        //fill templates binned by njets, sumjetpt
        fillUnderOverFlow( tcmetTemplate_combined_ee[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_ee    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_combined_ee[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_ee    , templateWeight );
        fillUnderOverFlow( tcmetNewTemplate_combined_ee[ iJetBin ][ iSumJetPtBin ] ,  tcmetNew_ee , templateWeight );

        //fill templates binned by njets, sumjetpt
        fillUnderOverFlow( tcmetTemplate_combined_mm[ iJetBin ][ iSumJetPtBin ]    ,  tcmet_mm    , templateWeight );
        fillUnderOverFlow( pfmetTemplate_combined_mm[ iJetBin ][ iSumJetPtBin ]    ,  pfmet_mm    , templateWeight );
        fillUnderOverFlow( tcmetNewTemplate_combined_mm[ iJetBin ][ iSumJetPtBin ] ,  tcmetNew_mm , templateWeight );
        */
      }

      int iJ = nJets_;
      if( iJ > 4 ) iJ = 4;
       
      fillUnderOverFlow( hpthad[iJ] ,  pthad , mcweight );
      fillUnderOverFlow( hpthad[0]  ,  pthad , mcweight );

      //fill predicted and observed met histos--------------------------------------------------

      if( !makeTemplate_ ){
        
        int iJetBin      = getJetBin( nJets_ );
        //int iBosonPtBin  = getBosonPtBin( dilpt_ );
        int iBosonPtBin  = getBosonPtBin( vecJetPt_ );
        int iSumJetPtBin = getSumJetPtBin( sumJetPt_ );
        int iVtxBin      = getVtxBin( nvtx_ );
        int iTrigBin     = getTrigBin( jetmax_pt_ );
        
        //cout << run_ << " " << lumi_ << " " << event_ << endl ;
        //cout << jetmax_pt_ << " " << nJets_ << " " << sumJetPt_ << endl;
        //cout << iTrigBin << " " << iJetBin << " " << iSumJetPtBin << endl;

        TH1F* hmet = getMetTemplate( metTemplateFile , iTrigBin , iJetBin , iSumJetPtBin , 
                                     iBosonPtBin , iVtxBin, dilpt_ , mcweight );

        hmet->Scale( mcweight );
        
        if( selection_ == e_ZSelection ){
          if( leptype_ == 2 ){
            fillUnderOverFlow( metObserved_df , theMet , mcweight  );
            metPredicted_df->Add( hmet );
            continue; //exclude emu final state
          }
        }
        
        fillUnderOverFlow( metObserved_njets[iJetBin] , theMet , mcweight );
        metPredicted_njets[iJetBin]->Add( hmet );
        
        fillUnderOverFlow( metObserved , theMet , mcweight );
        metPredicted->Add( hmet );
        
        //counters to track how many times each template is picked
        if( myTemplateSource   == e_PhotonJetStitched ){
          nphoton_metPredicted[ iBosonPtBin ][ iJetBin ][ iSumJetPtBin ]++;
          if( leptype_ == 0 ) nphoton_metPredicted_ee[ iBosonPtBin ][ iJetBin ][ iSumJetPtBin ]++;
          if( leptype_ == 1 ) nphoton_metPredicted_mm[ iBosonPtBin ][ iJetBin ][ iSumJetPtBin ]++;
        }
        else if( myTemplateSource   == e_QCD ){
          nqcd_metPredicted[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]++;
          if( leptype_ == 0 ) nqcd_metPredicted_ee[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]++;
          if( leptype_ == 1 ) nqcd_metPredicted_mm[ iTrigBin ][ iJetBin ][ iSumJetPtBin ]++;
        }
        else if( myTemplateSource   == e_PhotonJet ){
          n_metPredicted[ iJetBin ][ iSumJetPtBin ]++;
          if( leptype_ == 0 ) n_metPredicted_ee[ iJetBin ][ iSumJetPtBin ]++;
          if( leptype_ == 1 ) n_metPredicted_mm[ iJetBin ][ iSumJetPtBin ]++;
        }
        else{
          cout << "Error unrecognized template source! " << myTemplateSource << endl;
          exit(0);
        }


        if( selection_ == e_ZSelection ) {
                
          if( leptype_ == 0 ){
            fillUnderOverFlow( metObserved_ee , theMet , mcweight );
            metPredicted_ee->Add( hmet );
          }
          if( leptype_ == 1 ){
            fillUnderOverFlow( metObserved_mm , theMet , mcweight );
            metPredicted_mm->Add( hmet );
          }
        
          if( leptype_ == 0 || leptype_ == 1 ){
            fillUnderOverFlow( metObserved_sf , theMet , mcweight  );
            metPredicted_sf->Add( hmet );
          }
          else if( leptype_ == 2 ){
            //fillUnderOverFlow( metObserved_df , tcmet_ , mcweight  );
            //metPredicted_df->Add( hmet );
            cout << "SHOULD NOT GET HERE!!!!!!" << endl;
          }
          else{
            cout << "UNRECOGNIZED LEPTYPE " << leptype_ << endl;
          }
        }
                        
        if( selection_ == e_photonSelection ||  selection_ == e_ZSelection){

          if( pthad < 40 ){
            fillUnderOverFlow( metObserved_ptlt40 , theMet , mcweight );
            metPredicted_ptlt40->Add( hmet );
          }
          if( pthad > 40 && pthad < 60 ){
            fillUnderOverFlow( metObserved_pt40_60 , theMet , mcweight );
            metPredicted_pt40_60->Add( hmet );
          }
          if( pthad > 60 ){
            fillUnderOverFlow( metObserved_ptgt60 , theMet , mcweight );
            metPredicted_ptgt60->Add( hmet );
          }
          if( pthad < 50 ){
            fillUnderOverFlow( metObserved_ptlt50 , theMet , mcweight );
            metPredicted_ptlt50->Add( hmet );
          }else{
            fillUnderOverFlow( metObserved_ptgt50 , theMet , mcweight );
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
  
  eefile.close();
  mmfile.close();

  if( !makeTemplate_ && isData && setTemplateErrors ){

    if( myTemplateSource   == e_QCD ){
      setErrors( metTemplateFile , metPredicted    , nqcd_metPredicted    );
      setErrors( metTemplateFile , metPredicted_ee , nqcd_metPredicted_ee );
      setErrors( metTemplateFile , metPredicted_mm , nqcd_metPredicted_mm );
    }
    
    else if( myTemplateSource   == e_PhotonJetStitched ){
      setErrors( metTemplateFile , metPredicted    , nphoton_metPredicted    );
      setErrors( metTemplateFile , metPredicted_ee , nphoton_metPredicted_ee );
      setErrors( metTemplateFile , metPredicted_mm , nphoton_metPredicted_mm );
    }
    
    else if( myTemplateSource   == e_PhotonJet ){
      setErrors( metTemplateFile , metPredicted    , n_metPredicted    );
      setErrors( metTemplateFile , metPredicted_ee , n_metPredicted_ee );
      setErrors( metTemplateFile , metPredicted_mm , n_metPredicted_mm );
    }
  }


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
        for( int iVtxBin = 0 ; iVtxBin < nVtxBins ; iVtxBin++ ){
          
          float scale = tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Integral();
          if( scale > 0 )  tcmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Scale ( 1. / scale );
          
          scale = tcmetNewTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Integral();
          if( scale > 0 )  tcmetNewTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Scale ( 1. / scale );
          
          scale = pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Integral();
          if( scale > 0 )  pfmetTemplate_njets_ht_nvtx[ iJetBin ][ iSumJetPtBin ][ iVtxBin ] -> Scale ( 1. / scale );
          
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
  if( makeTemplate ) saveHist( Form("../output/%s/babylooper_%s_templates.root" , template_version , prefix ) );
  else               saveHist( Form("../output/%s/babylooper_%s%s%s.root"         , Z_version , prefix , metTemplateString.c_str() , iter ) );
  deleteHistos();
  
} // end ScanChain

void babylooper::setErrors( TFile* file,  TH1F* hist , int n[3][7] ){

  char* metstring = "";
  if     ( myMetType == e_tcmet    ) metstring = "tcmet";
  else if( myMetType == e_tcmetNew ) metstring = "tcmetNew";
  else if( myMetType == e_pfmet    ) metstring = "pfmet";

  cout << hist->Integral() << endl;

  int ntemplates = 0;
  
  for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
    for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      cout << jetString(iJetBin) << " " << sumJetPtString(iSumJetPtBin) << " " << n[iJetBin][iSumJetPtBin] << endl; 
      ntemplates += n[iJetBin][iSumJetPtBin];
    }
  }

  cout << "Num entries in templates " << ntemplates << endl;

  if( fabs( hist->Integral() - ntemplates ) > 1.e-3 ){
    cout << "Error: " << hist->Integral() << " hist entries does not match " << ntemplates << endl;
    //exit(0);
  }

  TH1F* hmet = new TH1F();

  for( int ibin = 1 ; ibin <= hist->GetXaxis()->GetNbins() ; ibin++ ){
    
    float err2 = 0;

    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){

        hmet = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin));
        
        err2 += pow( n[iJetBin][iSumJetPtBin] * hmet->GetBinError( ibin ) , 2);
        
      }
    }

    hist->SetBinError( ibin , sqrt( err2 ) );
  }

}

void babylooper::setErrors( TFile* file,  TH1F* hist , int n[4][3][7] ){

  cout << "setErrors: " << hist->GetName() << endl;

  char* metstring = "";
  if     ( myMetType == e_tcmet    ) metstring = "tcmet";
  else if( myMetType == e_tcmetNew ) metstring = "tcmetNew";
  else if( myMetType == e_pfmet    ) metstring = "pfmet";

  cout << hist->Integral() << endl;

  int ntemplates = 0;

  for( int i = 0 ; i < 4 ; i++ ){
    for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
        cout << i << " " << jetString(iJetBin) << " " << sumJetPtString(iSumJetPtBin) << " " << n[i][iJetBin][iSumJetPtBin] << endl; 
        ntemplates += n[i][iJetBin][iSumJetPtBin];
      }
    }
  }

  cout << "Num entries in templates " << ntemplates << endl;

  if( fabs( hist->Integral() - ntemplates ) > 1.e-3 ){
    cout << "Error: " << hist->Integral() << " hist entries does not match " << ntemplates << endl;
    //exit(0);
  }

  TH1F* hmet = new TH1F();

  for( int ibin = 1 ; ibin <= hist->GetXaxis()->GetNbins() ; ibin++ ){
    
    float err2 = 0;

    for( int i = 0 ; i < 4 ; i++ ){
      for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
        for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
          
          if( myTemplateSource   == e_QCD ){
            hmet     = (TH1F*) file->Get(Form("%sTemplate_qcd_%i_%i_%i",metstring,i,iJetBin,iSumJetPtBin));
          }
          
          else if( myTemplateSource   == e_PhotonJetStitched ){
            hmet     = (TH1F*) file->Get(Form("%sTemplate_photon_%i_%i_%i",metstring,i,iJetBin,iSumJetPtBin));
          }
          
          else {
            cout << "Invalid template source! " << myTemplateSource << endl;
            exit(0);
          }
          
          err2 += pow( n[i][iJetBin][iSumJetPtBin] * hmet->GetBinError( ibin ) , 2);
          
        }
      }
    }

    hist->SetBinError( ibin , sqrt( err2 ) );
  }

}


TH1F* babylooper::getMetTemplate( TFile* file, int iTrigBin , int iJetBin , 
                                  int iSumJetPtBin , int iBosonPtBin , int iVtxBin, float Zpt , float weight ){
  
  char* metstring = "";
  
  if     ( myMetType == e_tcmet    ) metstring = "tcmet";
  else if( myMetType == e_tcmetNew ) metstring = "tcmetNew";
  else if( myMetType == e_pfmet    ) metstring = "pfmet";
  
  TH1F* hmet = new TH1F();
  
  if( myTemplateType == e_njets_ht ){
    
    if( myTemplateSource   == e_QCD ){
      TH1F* htemp     = (TH1F*) file->Get(Form("%sTemplate_qcd_%i_%i_%i",metstring,iTrigBin,iJetBin,iSumJetPtBin));
      hmet = correctedMetTemplate( htemp , Zpt );
      //hmet     = (TH1F*) file->Get(Form("%sTemplate_qcd_%i_%i_%i",metstring,iTrigBin,iJetBin,iSumJetPtBin));
      //hmet     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin));
    }
    
    else if( myTemplateSource   == e_PhotonJetStitched ){
      hmet     = (TH1F*) file->Get(Form("%sTemplate_photon_%i_%i_%i",metstring,iBosonPtBin,iJetBin,iSumJetPtBin));
    }
    
    else if( myTemplateSource   == e_PhotonJet ){
      hmet     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin));
    }
    
    else{
      cout << "Error unrecognized template source! " << myTemplateSource << endl;
      exit(0);
    }
    
  }
  
  else{
    cout << "Error unrecognized template source! " << myTemplateSource << endl;
    exit(0);
  }
  
  if( hmet->GetEntries() < nMetEntries ){
    cout << "Error found template with " << hmet->GetEntries() << endl;
    //exit(0);
  }
  
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

  hnVtx     = new TH1F("hnVtx","",10,0,10);
  hvecJetPt = new TH1F("hvecJetPt","",50,0,500);
  hnVtx->Sumw2();
  hvecJetPt->Sumw2();

  hgenmet_all  = new TH1F("hgenmet_all","",100,0,200);
  hgenmet_pass = new TH1F("hgenmet_pass","",100,0,200);
  hresponse    = new TProfile("hresponse","",100,0,200,0,2);

  hgenmet_all->Sumw2();
  hgenmet_pass->Sumw2();
  
  hyield_0j = new TH1F("yield_0j","Event Yields (0 jets)",4,0,4);
  hyield_0j->GetXaxis()->SetTitle("dil type");
  hyield_0j->GetXaxis()->SetBinLabel(1,"all");
  hyield_0j->GetXaxis()->SetBinLabel(2,"ee");
  hyield_0j->GetXaxis()->SetBinLabel(3,"mm");
  hyield_0j->GetXaxis()->SetBinLabel(4,"em");
  hyield_0j->Sumw2();

  hyield_1j = new TH1F("yield_1j","Event Yields (1 jet)",4,0,4);
  hyield_1j->GetXaxis()->SetTitle("dil type");
  hyield_1j->GetXaxis()->SetBinLabel(1,"all");
  hyield_1j->GetXaxis()->SetBinLabel(2,"ee");
  hyield_1j->GetXaxis()->SetBinLabel(3,"mm");
  hyield_1j->GetXaxis()->SetBinLabel(4,"em");
  hyield_1j->Sumw2();

  hyield_2j = new TH1F("yield_2j","Event Yields (2 jets)",4,0,4);
  hyield_2j->GetXaxis()->SetTitle("dil type");
  hyield_2j->GetXaxis()->SetBinLabel(1,"all");
  hyield_2j->GetXaxis()->SetBinLabel(2,"ee");
  hyield_2j->GetXaxis()->SetBinLabel(3,"mm");
  hyield_2j->GetXaxis()->SetBinLabel(4,"em");
  hyield_2j->Sumw2();

  hyield_g2j = new TH1F("yield_g2j","Event Yields (>2 jets)",4,0,4);
  hyield_g2j->GetXaxis()->SetTitle("dil type");
  hyield_g2j->GetXaxis()->SetBinLabel(1,"all");
  hyield_g2j->GetXaxis()->SetBinLabel(2,"ee");
  hyield_g2j->GetXaxis()->SetBinLabel(3,"mm");
  hyield_g2j->GetXaxis()->SetBinLabel(4,"em");
  hyield_g2j->Sumw2();

  hyield = new TH1F("yield","Event Yields",4,0,4);
  hyield->GetXaxis()->SetTitle("dil type");
  hyield->GetXaxis()->SetBinLabel(1,"all");
  hyield->GetXaxis()->SetBinLabel(2,"ee");
  hyield->GetXaxis()->SetBinLabel(3,"mm");
  hyield->GetXaxis()->SetBinLabel(4,"em");
  hyield->Sumw2();

  hyield_pfmet30 = new TH1F("yield_pfmet30","Event Yields (pfmet > 30 GeV)",4,0,4);
  hyield_pfmet30->GetXaxis()->SetTitle("dil type");
  hyield_pfmet30->GetXaxis()->SetBinLabel(1,"all");
  hyield_pfmet30->GetXaxis()->SetBinLabel(2,"ee");
  hyield_pfmet30->GetXaxis()->SetBinLabel(3,"mm");
  hyield_pfmet30->GetXaxis()->SetBinLabel(4,"em");
  hyield_pfmet30->Sumw2();

  hyield_pfmet60 = new TH1F("yield_pfmet60","Event Yields (pfmet > 60 GeV)",4,0,4);
  hyield_pfmet60->GetXaxis()->SetTitle("dil type");
  hyield_pfmet60->GetXaxis()->SetBinLabel(1,"all");
  hyield_pfmet60->GetXaxis()->SetBinLabel(2,"ee");
  hyield_pfmet60->GetXaxis()->SetBinLabel(3,"mm");
  hyield_pfmet60->GetXaxis()->SetBinLabel(4,"em");
  hyield_pfmet60->Sumw2();

  hyield_pfmet120 = new TH1F("yield_pfmet120","Event Yields (pfmet > 120 GeV)",4,0,4);
  hyield_pfmet120->GetXaxis()->SetTitle("dil type");
  hyield_pfmet120->GetXaxis()->SetBinLabel(1,"all");
  hyield_pfmet120->GetXaxis()->SetBinLabel(2,"ee");
  hyield_pfmet120->GetXaxis()->SetBinLabel(3,"mm");
  hyield_pfmet120->GetXaxis()->SetBinLabel(4,"em");
  hyield_pfmet120->Sumw2();

  char* leptype[4]   = {"ee", "mm", "em", "all"};
  char* jetbin[4]    = {"0j", "1j", "geq2j", "allj"};

  char* leptype_title[4]   = {"ee", "#mu#mu", "e#mu", "all leptons"};
  char* jetbin_title[4]    = {"0 jets", "1 jet", "#geq 2 jets", "all jets"};

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {

      char* suffix       = Form("%s_%s",leptype[i],jetbin[j]);
      char* suffix_title = Form("%s %s",leptype_title[i],jetbin_title[j]);

      hdilmass[i][j]          = new TH1F(Form("hdilmass_%s",suffix),          suffix_title, 200,0,200);
      hdilmass_pfmet60[i][j]  = new TH1F(Form("hdilmass_pfmet60_%s",suffix),  suffix_title, 200,0,200);
      
      htcmet[i][j]    = new TH1F(Form("htcmet_%s",suffix),    suffix_title, 100,0,100);
      htcmetNew[i][j] = new TH1F(Form("htcmetNew_%s",suffix), suffix_title, 100,0,100);
      hpfmet[i][j]    = new TH1F(Form("hpfmet_%s",suffix),    suffix_title, 100,0,100);

      hdilmass[i][j]->GetXaxis()->SetTitle("dilepton mass (GeV)");
      hdilmass_pfmet60[i][j]->GetXaxis()->SetTitle("dilepton mass (GeV)");
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

  int maxmet = 250;

  if( !makeTemplate_ ){
    
    metObserved  = new TH1F("metObserved", "Observed MET",maxmet,0,maxmet);
    metPredicted = new TH1F("metPredicted","Predicted MET",maxmet,0,maxmet);
    
    metObserved_sf  = new TH1F("metObserved_sf", "Observed MET (SF)",maxmet,0,maxmet);
    metPredicted_sf = new TH1F("metPredicted_sf","Predicted MET (SF)",maxmet,0,maxmet);
    
    metObserved_df  = new TH1F("metObserved_df", "Observed MET (DF)",maxmet,0,maxmet);
    metObserved_df_nozveto  = new TH1F("metObserved_df_nozveto", "Observed MET (DF) No Z-veto",maxmet,0,maxmet);
    metPredicted_df = new TH1F("metPredicted_df","Predicted MET (DF)",maxmet,0,maxmet);
    
    metObserved_ptlt40     = new TH1F("metObserved_ptlt40", "Observed MET (p_{T}<40 GeV)",maxmet,0,maxmet);
    metPredicted_ptlt40    = new TH1F("metPredicted_ptlt40","Predicted MET  (p_{T}<40 GeV)",maxmet,0,maxmet);
    
    metObserved_ptlt50     = new TH1F("metObserved_ptlt50", "Observed MET (p_{T}<50 GeV)",maxmet,0,maxmet);
    metPredicted_ptlt50    = new TH1F("metPredicted_ptlt50","Predicted MET  (p_{T}<50 GeV)",maxmet,0,maxmet);
    
    metObserved_ptgt50     = new TH1F("metObserved_ptgt50", "Observed MET (p_{T}>50 GeV)",maxmet,0,maxmet);
    metPredicted_ptgt50    = new TH1F("metPredicted_ptgt50","Predicted MET  (p_{T}>50 GeV)",maxmet,0,maxmet);
    
    metObserved_pt40_60    = new TH1F("metObserved_pt40_60", "Observed MET (40<p_{T}<60 GeV)",maxmet,0,maxmet);
    metPredicted_pt40_60   = new TH1F("metPredicted_pt40_60","Predicted MET  (40<p_{T}<60 GeV)",maxmet,0,maxmet);
    
    metObserved_ptgt60     = new TH1F("metObserved_ptgt60", "Observed MET (p_{T}>60 GeV)",maxmet,0,maxmet);
    metPredicted_ptgt60    = new TH1F("metPredicted_ptgt60","Predicted MET  (p_{T}>60 GeV)",maxmet,0,maxmet);
    
    metObserved_ee  = new TH1F("metObserved_ee", "Observed MET (ee)",maxmet,0,maxmet);
    metPredicted_ee = new TH1F("metPredicted_ee","Predicted MET (ee)",maxmet,0,maxmet);
    
    metObserved_mm  = new TH1F("metObserved_mm", "Observed MET (#mu#mu)",maxmet,0,maxmet);
    metPredicted_mm = new TH1F("metPredicted_mm","Predicted MET (#mu#mu)",maxmet,0,maxmet);
    
    metParObserved  = new TH1F("metParObserved", "Observed MET (Parallel)",1000,-maxmet,maxmet);
    metParPredicted = new TH1F("metParPredicted","Predicted MET (Parallel)",1000,-maxmet,maxmet);
    
    metPerpObserved  = new TH1F("metPerpObserved", "Observed MET (Perpendicular)",maxmet,0,maxmet);
    metPerpPredicted = new TH1F("metPerpPredicted","Predicted MET (Perpendicular)",maxmet,0,maxmet);
    
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
      
      metObserved_njets[iJetBin]  = new TH1F(Form("metObserved_njets%i",iJetBin), Form("Observed MET NJets %i", iJetBin),maxmet,0,maxmet);
      metPredicted_njets[iJetBin] = new TH1F(Form("metPredicted_njets%i",iJetBin),Form("Predicted MET NJets %i",iJetBin),maxmet,0,maxmet);
      
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
                                                                            bosonPtString(iBosonPtBin).c_str()),maxmet,0,maxmet);

          tcmetNewTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("tcmetNewTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                          Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                               bosonPtString(iBosonPtBin).c_str()),maxmet,0,maxmet);

          pfmetTemplate[iJetBin][iSumJetPtBin][iBosonPtBin] = new TH1F(Form("pfmetTemplate_%i_%i_%i",iJetBin,iSumJetPtBin,iBosonPtBin),
                                                                       Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                            bosonPtString(iBosonPtBin).c_str()),maxmet,0,maxmet);
        
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
        for( int iVtxBin = 0 ; iVtxBin < nVtxBins ; iVtxBin++ ){

          tcmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin] = new TH1F(Form("tcmetTemplate_njets_ht_nvtx_%i_%i_%i",iJetBin,iSumJetPtBin,iVtxBin),
                                                                                 Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                                      nVTXString(iVtxBin).c_str()),maxmet,0,maxmet);
          
          tcmetNewTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin] = new TH1F(Form("tcmetNewTemplate_njets_ht_nvtx_%i_%i_%i",iJetBin,iSumJetPtBin,iVtxBin),
                                                                                    Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                                         nVTXString(iVtxBin).c_str()),maxmet,0,maxmet);
          
          pfmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin] = new TH1F(Form("pfmetTemplate_njets_ht_nvtx_%i_%i_%i",iJetBin,iSumJetPtBin,iVtxBin),
                                                                                 Form("%s, %s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str(),
                                                                                      nVTXString(iVtxBin).c_str()),maxmet,0,maxmet);
          
          tcmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->Sumw2();
          tcmetNewTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->Sumw2();
          pfmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->Sumw2();
          
          tcmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->GetXaxis()->SetTitle("tcmet (GeV)");
          tcmetNewTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->GetXaxis()->SetTitle("tcmetNew (GeV)");
          pfmetTemplate_njets_ht_nvtx[iJetBin][iSumJetPtBin][iVtxBin]->GetXaxis()->SetTitle("pfmet (GeV)");          
          
        }
      }
    }

  
   for( int iJetBin = 0 ; iJetBin < nJetBins ; iJetBin++ ){
      for( int iSumJetPtBin = 0 ; iSumJetPtBin < nSumJetPtBins ; iSumJetPtBin++ ){
      
        tcmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                 Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
        
        tcmetNewTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("tcmetNewTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                    Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
        
        pfmetTemplate_combined[iJetBin][iSumJetPtBin] = new TH1F(Form("pfmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin),
                                                                 Form("%s, %s",jetString(iJetBin).c_str(),sumJetPtString(iSumJetPtBin).c_str()),maxmet,0,maxmet);
        
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
  tree->SetBranchAddress("vtxweight",    &vtxweight_    );
  tree->SetBranchAddress("davtxweight",  &davtxweight_  );
  tree->SetBranchAddress("pthat",        &pthat_        );
  tree->SetBranchAddress("lumi",         &lumi_         );
  tree->SetBranchAddress("event",        &event_        );

  tree->SetBranchAddress("nvtx",         &nvtx_         );
  tree->SetBranchAddress("npfmuons",     &npfmuons_     );
  tree->SetBranchAddress("pfmet",        &pfmet_        );
  tree->SetBranchAddress("pfmetcor",     &pfmetcor_     );
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
  tree->SetBranchAddress("tcmetNew",     &tcmetNew_     );
  tree->SetBranchAddress("tcmetphiNew",  &tcmetphiNew_  );     
  tree->SetBranchAddress("tcsumetNew",   &tcsumetNew_   );  
  tree->SetBranchAddress("dphixmet",     &dphixmet_     );
  tree->SetBranchAddress("metpar",       &metPar_       );
  tree->SetBranchAddress("metperp",      &metPerp_      );
  tree->SetBranchAddress("njets",        &nJets_        );     
  tree->SetBranchAddress("njets40",      &nJets40_      );     
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
    tree->SetBranchAddress("pfjetid",               &pfjetid_               );           
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

    tree->SetBranchAddress("drjetll",      &drjet_ll_     ); 
    tree->SetBranchAddress("failjetid",    &failjetid_    ); 
    tree->SetBranchAddress("jetptll",      &jetpt_ll_     ); 
    tree->SetBranchAddress("pfjetidll",    &pfjetid_ll_   ); 
    tree->SetBranchAddress("drjetlt",      &drjet_lt_     ); 
    tree->SetBranchAddress("jetptlt",      &jetpt_lt_     ); 
    tree->SetBranchAddress("pfjetidlt",    &pfjetid_lt_   ); 

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
    tree->SetBranchAddress("etall",                &etall_                );  
    tree->SetBranchAddress("etalt",                &etalt_                ); 
    tree->SetBranchAddress("phill",                &phill_                );  
    tree->SetBranchAddress("philt",                &philt_                ); 
    tree->SetBranchAddress("dilmass",              &dilmass_              ); 
    tree->SetBranchAddress("dilmasspf",            &dilmasspf_            ); 
    tree->SetBranchAddress("dilmasscor",           &dilmasscor_           ); 
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

TH1F* babylooper::correctedMetTemplate( TH1F* h_metTemplate , float ptZ ){

  TH1F* h_metEstimate = (TH1F*) h_metTemplate->Clone();
  h_metEstimate->Reset();

  //the MET templates give the fake MET. some fraction of Z-boson PT's
  //gives the MET scale shift in each Z+jets event.
  
  //The approach is to emulate the MET scale shift (the ET under-measurements of 
  //the jets recoiling against the Z boson) by convoluting the 
  //MET template with a fraction of the Z PT. 
  
  // ============ Here is how this is done ================
  
  //Below, h_metEstimate holds the prediction for a given Z+jets event.
  //h_metTemplate is a MET template for a given Z+jets event;
  //this is its definition: TH1D(histName,  histName, 1000, 0, 1000);
  //it should be normalized to 1.0.
  
  const double   MyPi = 3.141592654; 
  
  //ptZ is the Z-boson PT in the event. In data the MET scale 
  //shift is roughly 5%, so the absolute MET scale shift from the 
  //system of jets in this event is: 
  
  double metFromZ = ptZ * 0.075;
  
  //This is a loop over the angle (myphi) btw the fake MET PT-vector
  //and the Z-boson PT vector. Assuming they interfere at a random angle, 
  //I add them vectorially averaging over 2Pi:
  
  for( int phiCounter = 0; phiCounter < 20; phiCounter++ )
    {
      float myphi = phiCounter * MyPi * 0.05;

      //loop over the met bins in h_metTemplate

      for( int metCounter = 1 ; metCounter <= h_metEstimate->GetNbinsX() ; metCounter++ )
	{		     		
	  double met      = metCounter-0.5;
	  
	  double metStarX = metFromZ + met*cos(myphi);
	  double metStarY =            met*sin(myphi);
	  
	  double metStar  =  sqrt(pow(metStarX,2)+pow(metStarY,2));
	  
	  double newBin   = int( metStar );
	  
	  //add this to the prediction for this Z+jets event: 
	  h_metEstimate->Fill( newBin, h_metTemplate->GetBinContent(metCounter) );
	  
	}
    }
  
  //normalize h_metEstimate to 1.0
  h_metEstimate->Scale( 1.0 / h_metEstimate->Integral() );
  
  return h_metEstimate;
  
}




/*
  TH1F* hmet = new TH1F();
  
  //if( !isData && nJets_ > 1 ){
  if( myMetType == e_tcmet ){
  hmet     = (TH1F*) metTemplateFile->Get(Form("tcmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin));
  }
  else if( myMetType == e_tcmetNew ){
  hmet     = (TH1F*) metTemplateFile->Get(Form("tcmetNewTemplate_combined_%i_%i",iJetBin,iSumJetPtBin));
  }
  else if( myMetType == e_pfmet ){
  hmet     = (TH1F*) metTemplateFile->Get(Form("pfmetTemplate_combined_%i_%i",iJetBin,iSumJetPtBin));
  }
*/
//}
/*
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
*/










/*
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

string sumJetPtString( int bin ){

  //float bins[nSumJetPtBins+1]={0,25,50,75,100,150,250,5000};
  float bins[nSumJetPtBins+1]={0,30,60,90,120,150,250,5000};

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
*/

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
*/














        
        //cout << __LINE__ << endl;
        //-------------------------------------------------------
        //selection which is applied upstream in isGoodEMObject
        //-------------------------------------------------------
        //if ( etg_ < 22 )                                 continue;
        //if ( (1.-r4_) < 0.05 )                           continue;
        //if ( hoe_ > 0.1 )                                continue;
        //if ( jet_neu_emfrac_ < 0.95 )                    continue; 
        //if ( jet_dr_ > 0.3 )                             continue;
        //if ( nvtx_ != 1 )                                continue;
        
        //-------------------------------------------------------
        //Warren's photon selection
        //-------------------------------------------------------
        //if ( etag_ > 1 )                                 continue;
        //if ( etg_ < 20 )                                 continue;
        //if ( (1.-r4_) < 0.05 )                           continue;
        //if ( hoe_ > 0.05 )                               continue;
        //if ( photon_hcalIso03_ > 2.4 && photon_hcalIso03_ / etg_ > 0.05)        continue;              
        //if ( photon_ecalIso03_ > 1.7 && photon_ecalIso03_ / etg_ > 0.05)        continue;       
        //if ( photon_ntkIsoSolid03_ > 2 || photon_tkisoSolid03_ / etg_ > 0.1 )   continue;
            
        //-------------------------------------------------------
        //EWK PAS selection
        //-------------------------------------------------------
        //if( photon_ecalIso04_ > 4.2 * 0.004 * etg_ )       continue;
        //if( photon_ecalIso04_ > 4.2 * 0.004 * etg_ )       continue;
        
        //-------------------------------------------------------
        //misc selection
        //-------------------------------------------------------
        //if ( jet_eta_ > 1 )                              continue;
        //if ( jet_dr_ > 0.5 )                             continue;
        //if ( jet_neu_emfrac_ + jet_chg_emfrac_< 0.95 )   continue; 
        //if ( drel_ < 0.5 )                               continue;
        //if ( ( 1 - jet_neu_emfrac_ ) * etg_ > 1 )        continue;
        //if ( sumJetPt_ < 200. )                          continue;
        

        //dphixmet_  = deltaPhi( tcmetphi_ , phig_ );
        //metPar_    = theMet * cos( dphixmet_ );
        //metPerp_   = theMet * sin( dphixmet_ );

  
        //------------------------------------------------------------
        //misc selections
        //------------------------------------------------------------
        //if( ptll_ < 20 || ptlt_ < 20 )                  continue;
        //if( leptype_ == 1 ){
        //  if( !passm_ll_nom_ == 1) continue;
        //  if( !passm_lt_nom_ == 1) continue;
        //}
        //if( nJets40_ < 2 )                              continue;
        //if( ptll_ < 20 || ptlt_ < 20 )                  continue;
        //if( fabs(etall_) > 2.4 || fabs(etalt_) > 2.4 )  continue;
        //if( nvtx_ < 4 )                                 continue;
        
        //if( leptype_ == 1 ){
        //  if( !passm_ll_nom_ == 1) continue;
        //  if( !passm_lt_nom_ == 1) continue;
        //}

        //if( leptype_ == 2 )                   continue;
        //if( nvtx_    != 1 )                   continue;
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
 
     

 /*
TH1F* babylooper::getMetTemplate( TFile* file, int iTrigBin , int iJetBin , 
                                  int iSumJetPtBin , int iBosonPtBin , int iVtxBin, float weight ){

  char* metstring = "";

  if     ( myMetType == e_tcmet    ) metstring = "tcmet";
  else if( myMetType == e_tcmetNew ) metstring = "tcmetNew";
  else if( myMetType == e_pfmet    ) metstring = "pfmet";
    
  TH1F* hmet = new TH1F();



  
  if( myTemplateType == e_njets_ht ){

    
    if( myTemplateSource   == e_QCD ){
      //cout << "Taking QCD template" << endl;
      //cout << Form("%sTemplate_qcd_%i_%i_%i",metstring,iTrigBin,iJetBin,iSumJetPtBin) << endl;
      hmet     = (TH1F*) file->Get(Form("%sTemplate_qcd_%i_%i_%i",metstring,iTrigBin,iJetBin,iSumJetPtBin));
      //cout << hmet->GetTitle() << endl;
    }

    else if( myTemplateSource   == e_PhotonJetStitched ){
      //cout << "Taking photon template" << endl;
      //iBosonPtBin = 0;
      hmet     = (TH1F*) file->Get(Form("%sTemplate_photon_%i_%i_%i",metstring,iBosonPtBin,iJetBin,iSumJetPtBin));
    }

    else if( myTemplateSource   == e_PhotonJet ){
      hmet     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin));
    }
    else{
      cout << "Error unrecognized template source!" << endl;
      exit(0);
    }
    
  }
  else if( myTemplateType == e_njets_ht_nvtx ){
    hmet     = (TH1F*) file->Get(Form("%sTemplate_njets_ht_nvtx_%i_%i_%i",metstring,iJetBin,iSumJetPtBin,iVtxBin));
  }
  else if( myTemplateType == e_njets_ht_vecjetpt ){
    hmet     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin,iBosonPtBin));
  }
  else{
    cout << "Error unrecognized templateType " << myTemplateType << endl;
    exit(0);
  }
  

  //cout << endl << "Found template with " << hmet->GetEntries() << " entries" << endl;
  //cout << jetString(iJetBin) << " " << sumJetPtString(iSumJetPtBin) << " " << nVTXString(iVtxBin) << endl;
  
  //int nentries = hmet->GetEntries();
   
  //cout << iJetBin << " " << iSumJetPtBin << " " << nentries << endl;
  //cout << endl << "Found template with " << hmet->GetEntries() << " entries" << endl;
  //cout << jetString(iJetBin) << " " << sumJetPtString(iSumJetPtBin) << endl;

  //if there are at least nMetEntries in the template, return template
  if( hmet->GetEntries() >= nMetEntries ){
    hmet->Scale( weight );

    return hmet;
  }
  
  cout << "FOUND TEMPLATE WITH " << hmet->GetEntries() << " QUITTING!!!!!!" << endl;
  exit(0);

  if( hmet->GetEntries() > 0 )
    hmet->Scale( hmet->GetEntries() / hmet->Integral() );
  int counter = 1;
 
  //cout << "Less than " << nMetEntries << " entries!!!!!!!!!!" << endl;
  //cout << "Found template with " << hmet->GetEntries() << " entries" << endl;
  //cout << iJetBin << " " << iSumJetPtBin << endl;
  //cout << jetString(iJetBin) << " " << sumJetPtString(iSumJetPtBin) << endl;

  //add histograms in sumjetpt bins adjacent to current bin to get enough entries
  while( hmet->GetEntries() < nMetEntries && counter < nSumJetPtBins ){

    //add template 1 bin lower in sumJetPt, if it exists
    if( iSumJetPtBin - counter >= 0 ){
     
      TH1F *hmetLo = new TH1F();

      if( myTemplateType == e_njets_ht ){
        hmetLo     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin-counter));
      }
      else if( myTemplateType == e_njets_ht_nvtx ){
        hmetLo     = (TH1F*) file->Get(Form("%sTemplate_njets_ht_nvtx_%i_%i_%i",metstring,iJetBin,iSumJetPtBin-counter,iVtxBin));
      }
      else if( myTemplateType == e_njets_ht_vecjetpt ){
        hmetLo     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin-counter,iBosonPtBin));
      }
      else{
        cout << "Error unrecognized templateType " << myTemplateType << endl;
        exit(0);
      }

      //if( myTemplateType == e_njets_ht ){
      //  hmetLo     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin-counter));
      //}else{
      //  hmetLo     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin-counter,iBosonPtBin));
      //}
 
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

      if( myTemplateType == e_njets_ht ){
        hmetHi     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin+counter));
      }
      else if( myTemplateType == e_njets_ht_nvtx ){
        hmetHi     = (TH1F*) file->Get(Form("%sTemplate_njets_ht_nvtx_%i_%i_%i",metstring,iJetBin,iSumJetPtBin+counter,iVtxBin));
      }
      else if( myTemplateType == e_njets_ht_vecjetpt ){
        hmetHi     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin+counter,iBosonPtBin));
      }
      else{
        cout << "Error unrecognized templateType " << myTemplateType << endl;
        exit(0);
      }

//       if( myTemplateType == e_njets_ht ){
//         hmetHi     = (TH1F*) file->Get(Form("%sTemplate_combined_%i_%i",metstring,iJetBin,iSumJetPtBin+counter));
//       }else{
//         hmetHi     = (TH1F*) file->Get(Form("%sTemplate_%i_%i_%i",metstring,iJetBin,iSumJetPtBin+counter,iBosonPtBin));
//       }
      
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
 */


