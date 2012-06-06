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

enum metType        { e_tcmet    = 0 , e_tcmetNew      = 1 , e_pfmet             = 2 };
enum templateSource { e_QCD      = 0 , e_PhotonJet     = 1 , e_PhotonJetStitched = 2 };
enum templateType   { e_njets_ht = 0 , e_njets_ht_nvtx = 1 , e_njets_ht_vecjetpt = 2 };

//------------------------------------------------
//USER PARAMS
//------------------------------------------------
bool           debug              = false;                  // debug printout statements
bool           doVtxReweight      = true;                   // reweight templates for nVertices
bool           setTemplateErrors  = true;                   // calculate template errors
metType        myMetType          = e_pfmet;                // MET type
templateSource myTemplateSource   = e_PhotonJetStitched;    // source of templates
templateType   myTemplateType     = e_njets_ht;             // bin templates in njets and HT
bool           reweight           = false;                  // reweight for photon vs. Z pt
char*          iter               = "njetsgeq2";            // label for output file
float          lumi               = 2.4;                    // luminosity
//------------------------------------------------

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void babylooper::ScanChain (TChain* chain, const char* Z_version, const char* template_version, const char* prefix, 
			    bool isData, selectionType mySelectionType, bool makeTemplate, int nEvents){

  int npass = 0;
  selection_    = mySelectionType;
  makeTemplate_ = makeTemplate;

  // if( !isData && myTemplateSource == e_PhotonJetStitched ){
  //   myTemplateSource = e_PhotonJet;
  //   cout << "Switching MC template to PhotonJet" << endl;
  // }

  if( selection_ == e_ZSelection )       cout << "Z selection" << endl;
  else { cout << "Unrecognized selection type, quitting" << endl ; exit(0) ; }

  // if     ( selection_ == e_QCDSelection )     cout << "QCD selection" << endl;
  // else if( selection_ == e_photonSelection )  cout << "photon selection" << endl;
  // else if( selection_ == e_ZSelection )       cout << "Z selection" << endl;
  // else { cout << "Unrecognized selection type, quitting" << endl ; exit(0) ; }

  TFile*  metTemplateFile;//   = new TFile();
  string  metTemplateString = "";
  char*   templateFileName  = "";

  TRandom3 rand;
  //  if( isData ){

  if( myTemplateSource == e_PhotonJetStitched ){
    if( doVtxReweight ) templateFileName = Form("../photon_output/%s/DoubleElectron_templates_vtxreweight.root",template_version);
    else                templateFileName = Form("../photon_output/%s/DoubleElectron_templates.root",template_version);
    cout << "Using template file " << templateFileName << endl;
    metTemplateString = "_PhotonStitchedTemplate";
    metTemplateFile = TFile::Open( templateFileName );
  }else{
    cout << __FILE__ << " " << __LINE__ << " ERROR UNRECOGNIZED TEMPLATES" << endl;
    exit(0);
  }
      
  // }else{
      
  //   if( myTemplateSource == e_QCD ){
  //     cout << "QCD templates are deprecated. If you want to use this you"
  // 	   << "need to make QCD templates using makeTemplates.C." << endl;
  //     exit(0);
  //     //templateFileName = Form("output/%s/QCD_Pt15_templates.root",iter);
  //     //cout << "Using template file " << templateFileName << endl;
  //     //metTemplateString = "_QCD_Pt15Template";
  //     //metTemplateFile = TFile::Open( templateFileName ); 
  //   }
      
  //   else if( myTemplateSource == e_PhotonJet || myTemplateSource == e_PhotonJetStitched ){
  //     cout << "WARNING! USING DATA TEMPLATES FOR MC" << endl;
  //     templateFileName = Form("../templates/%s/photon_templates.root",template_version);
  //     //templateFileName = "/tas03/home/benhoob/metTemplate/output/nov5th/babylooper_PhotonJet_templates.root";
  //     cout << "Using template file " << templateFileName << endl;
  //     metTemplateString = "_PhotonJetTemplate";
  //     metTemplateFile = TFile::Open( templateFileName );
  //   }
  // }
  
  
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

  int neeall  = 0;
  int nmmall  = 0;
  int nee0jet = 0;
  int nmm0jet = 0;
  int nee1jet = 0;
  int nmm1jet = 0;
  int nee2jet = 0;
  int nmm2jet = 0;

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
      if( !isData ) mcweight = lumi * weight_ * vtxweight_;

      fillUnderOverFlow( hgenps_pthat  , pthat_ , mcweight );
      fillUnderOverFlow( hphotonpt     , etg_   , mcweight );
 
      float theMet = -1;
      if     ( myMetType == e_tcmet    ) theMet = tcmet_;
      else if( myMetType == e_tcmetNew ) theMet = tcmetNew_;
      else if( myMetType == e_pfmet    ) theMet = pfmet_;

      //---------------------------------------------------
      // apply event selection
      //---------------------------------------------------
      
      if( selection_ == e_ZSelection ) {

        // if( leptype_ == 0 ) eefile << dilmasscor_ << endl;
        // if( leptype_ == 1 ) mmfile << dilmasscor_ << endl;
     
        // if( dilmass_ > 76. && dilmass_ < 106. ){

        //   if( nJets_ == 0 ){
        //     hyield_0j->Fill(0.5,          mcweight);
        //     hyield_0j->Fill(1.5+leptype_, mcweight);
        //   }
        //   if( nJets_ == 1 ){
        //     hyield_1j->Fill(0.5,          mcweight);
        //     hyield_1j->Fill(1.5+leptype_, mcweight);
        //   }
        //   if( nJets_ == 2 ){
        //     hyield_2j->Fill(0.5,          mcweight);
        //     hyield_2j->Fill(1.5+leptype_, mcweight);
        //   }
        //   if( nJets_ > 2 ){
        //     hyield_g2j->Fill(0.5,          mcweight);
        //     hyield_g2j->Fill(1.5+leptype_, mcweight);
        //   }
        // }


        //------------------------------------------------------------
        //event selection
        //------------------------------------------------------------

	if( pflep1_->pt() < 20 )           continue; // PF lepton 1 pt > 20 GeV
	if( pflep2_->pt() < 20 )           continue; // PF lepton 2 pt > 20 GeV
	if( el1tv_ == 1 || el2tv_ == 1 )   continue; // veto transition region electrons
	//if( lep3_->pt() > 10 )             continue; // 3rd lepton veto

	// if( leptype_ == 0 ){
	//   if( ee_ == 0 ) continue;
	// }

	// if( leptype_ == 1 ){
	//   if( mm_ == 0 ) continue;
	// }

	// if( leptype_ == 2 ){
	//   if( em_ == 0 && me_ == 0 ) continue;
	// }

	//if( fabs( lep1_->pt() - pflep1_->pt() ) > 5.0 ) continue;
	//if( fabs( lep2_->pt() - pflep2_->pt() ) > 5.0 ) continue;

	//if( pflep1_->pt() > 100.0 )        continue;
	//if( pflep2_->pt() > 100.0 )        continue;

	//float pfdilmass = (*pflep1_ + *pflep2_).mass();

	// if( fabs( pfdilmass - dilmasspf_ ) > 0.01 ){
	//   cout << pfdilmass << " " << dilmasspf_ << " " << pfdilmass-dilmasspf_ << endl;
	// }

	if( leptype_ == 0 ){
	  if( jetpt_ll_ - ptll_ < -5  ) continue; 
	  if( jetpt_lt_ - ptlt_ < -5  ) continue; 
	}

	if( leptype_ == 0 ) fillUnderOverFlow( hpfdilmassee , dilmasspf_ , 1 );
	if( leptype_ == 1 ) fillUnderOverFlow( hpfdilmassmm , dilmasspf_ , 1 );

	if( dilmasspf_ > 81 && dilmasspf_ < 101 ){

	  //if( nJets_ >= 2 ){

	    if( leptype_ == 0 ){
	      float eta1 = lep1_->eta();
	      float eta2 = lep2_->eta();

	      if( fabs(eta1) < 1.4442 && fabs(eta2) < 1.4442 ){
		fillUnderOverFlow( hpfmet_ebeb , pfmet_ , 1 );
	      }

	      else if( ( fabs(eta1) < 1.4442 && fabs(eta2) > 1.556 ) || ( fabs(eta1) > 1.556 && fabs(eta2) < 1.4442 ) ){
		fillUnderOverFlow( hpfmet_ebee , pfmet_ , 1 );
	      }

	      else if( (eta1 >  1.556 && eta2 >  1.556) || (eta1 < -1.556 && eta2 < -1.556) ){
		fillUnderOverFlow( hpfmet_eeeep , pfmet_ , 1 );
	      }

	      else if( (eta1 >  1.556 && eta2 <  -1.556) || (eta1 < -1.556 && eta2 > 1.556) ){
		fillUnderOverFlow( hpfmet_eeeem , pfmet_ , 1 );
		cout << "eeeem" << endl;

	      }
	    }

	    if( leptype_ == 1 ){
	      fillUnderOverFlow( hpfmet_mm , pfmet_ , 1 );
	    }
	    //}

	  if( leptype_ == 0 ) neeall++;
	  if( leptype_ == 1 ) nmmall++;

	  if( nJets_ == 0 ){
	    if( leptype_ == 0 ) nee0jet++;
	    if( leptype_ == 1 ) nmm0jet++;
	  }
	  if( nJets_ == 1 ){
	    if( leptype_ == 0 ) nee1jet++;
	    if( leptype_ == 1 ) nmm1jet++;
	  }
	  if( nJets_ >= 2 ){
	    if( leptype_ == 0 ) nee2jet++;
	    if( leptype_ == 1 ) nmm2jet++;
	  }
	}

	if( nJets_ < 2 )                 continue; // >=2 jets

	//if( fabs( lep1_->eta() ) > 1.479 ) continue;
	//if( fabs( lep2_->eta() ) > 1.479 ) continue;
	
	// if( leptype_ == 1 ){
	//   if( fabs( ptllpf_ - ptll_ ) > 1  ) continue; 
	//   if( fabs( ptltpf_ - ptlt_ ) > 1  ) continue; 
	// }

	//if( fabs( lep1_->pt() - pflep1_->pt() ) > 1.0 ) continue;
	//if( fabs( lep2_->pt() - pflep2_->pt() ) > 1.0 ) continue;

	//if( leptype_ == 1 && npfmuons_ < 2 ) continue;
	//if( leptype_ == 2 && npfmuons_ < 1 ) continue;

	// //jetselection
        // if( nJets_ < 2 )                                                     continue; //>=2 jets

	//------------------------------
	//fill histos before Z-veto
	//------------------------------

        fillHistos( hdilmass , dilmass_ , mcweight , leptype_ , nJets_ );
	if( pfmet_ > 60. ){
	  fillHistos( hdilmass_pfmet60 , dilmass_ , mcweight , leptype_ , nJets_ );
	}
	if( leptype_ == 2) fillUnderOverFlow( metObserved_df_nozveto , theMet , mcweight  );

	// //Zmass
	// if( leptype_ == 1 ){
	//   if( dilmasspf_ < 81. || dilmasspf_ > 101. )                          continue; 
	// }else{
	//   if( dilmass_ < 81. || dilmass_ > 101. )                              continue; 
	// }
 


	//cout << dilmass_ << " " << pfdilmass << endl;
	//if( dilmass_ < 81. || dilmass_ > 101. )                              continue; 
	if( dilmasspf_ < 81. || dilmasspf_ > 101. )                              continue; 
	//if( pfdilmass < 86. || pfdilmass > 96. )                              continue; 

	hresponse->Fill( genmet_ , pfmet_ / genmet_ );
	hgenmet_all->Fill( genmet_ );
	if( pfmet_ > 60 ) hgenmet_pass->Fill( genmet_ );

        //------------------------------------------------------------
        //reweighting procedure for photon vs. Z pt
        //------------------------------------------------------------
        // if( reweight ){
        //   int iJetBin   = getJetBin( nJets_ );
        //   int bin       = reweightHist[iJetBin]->FindBin( dilpt_ );
        //   float ratio   = reweightHist[iJetBin]->GetBinContent( bin );
        //   float randnum = rand.Uniform(0,1);

        //   //cout << iJetBin << " " << dilpt_ << " " << bin << " " << ratio << " " << randnum << endl;
        //   if( randnum > ratio ) continue;
        //   //cout << "pass" << endl;
        // }
            
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

      int iJ = nJets_;
      if( iJ > 4 ) iJ = 4;
       
      fillUnderOverFlow( hpthad[iJ] ,  pthad , mcweight );
      fillUnderOverFlow( hpthad[0]  ,  pthad , mcweight );

      //fill predicted and observed met histos--------------------------------------------------
        
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
        
      
    } // end loop over events
  } // end loop over files
  
  cout << npass << " events passing selection" << endl;
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  eefile.close();
  mmfile.close();

  if( isData && setTemplateErrors ){

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

  float R     = sqrt((float) nmmall / (float) neeall);
  float Rerr  = R * sqrt( (1.0/(float)nmmall) + (1.0/(float)neeall));

  float R0    = sqrt((float) nmm0jet / (float) nee0jet);
  float R0err = R0 * sqrt( (1.0/(float)nmm0jet) + (1.0/(float)nee0jet));

  float R1    = sqrt((float) nmm1jet / (float) nee1jet);
  float R1err = R1 * sqrt( (1.0/(float)nmm1jet) + (1.0/(float)nee1jet));

  float R2    = sqrt((float) nmm2jet / (float) nee2jet);
  float R2err = R2 * sqrt( (1.0/(float)nmm2jet) + (1.0/(float)nee2jet));

  cout << endl;
  cout << "all" << endl;
  cout << "nee " << neeall << endl;
  cout << "nmm " << nmmall << endl;
  cout << "R   " << Form("%.2f +/- %.3f",R,Rerr) << endl;

  cout << endl;
  cout << "0 jets" << endl;
  cout << "nee " << nee0jet << endl;
  cout << "nmm " << nmm0jet << endl;
  cout << "R   " << Form("%.2f +/- %.3f",R0,R0err) << endl;
  
  cout << endl;
  cout << "1 jet" << endl;
  cout << "nee " << nee1jet << endl;
  cout << "nmm " << nmm1jet << endl;
  cout << "R   " << Form("%.2f +/- %.3f",R1,R1err) << endl;
  
  cout << endl;
  cout << "2 jets" << endl;
  cout << "nee " << nee2jet << endl;
  cout << "nmm " << nmm2jet << endl;
  cout << "R   " << Form("%.2f +/- %.3f",R2,R2err) << endl;
  


  // make histos rootfile
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  
  saveHist( Form("../output/%s/babylooper_%s%s%s.root"         , Z_version , prefix , metTemplateString.c_str() , iter ) );
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

  hpfdilmassee = new TH1F("hpfdilmassee","",200,0,200);
  hpfdilmassmm = new TH1F("hpfdilmassmm","",200,0,200);

  hpfmet_ebeb   = new TH1F("hpfmet_ebeb" ,"",100,0,200);
  hpfmet_ebee   = new TH1F("hpfmet_ebee" ,"",100,0,200);
  hpfmet_eeeep  = new TH1F("hpfmet_eeeep","",100,0,200);
  hpfmet_eeeem  = new TH1F("hpfmet_eeeem","",100,0,200);
  hpfmet_mm     = new TH1F("hpfmet_mm"   ,"",100,0,200);

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

  int maxmet = 300;


    
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


void babylooper::setBranches (TTree* tree){
  tree->SetBranchAddress("run",          &run_          );
  tree->SetBranchAddress("weight",       &weight_       );
  tree->SetBranchAddress("vtxweight",    &vtxweight_    );
  //tree->SetBranchAddress("davtxweight",  &davtxweight_  );
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
  // tree->SetBranchAddress("HLT_Jet15U",            &HLT_Jet15U_            ); 
  // tree->SetBranchAddress("HLT_Jet30U",            &HLT_Jet30U_            ); 

  // tree->SetBranchAddress("HLT_Photon10_L1R",      &HLT_Photon10_L1R_      ); 
  // tree->SetBranchAddress("HLT_Photon15_L1R",      &HLT_Photon15_L1R_      ); 
  // tree->SetBranchAddress("HLT_Photon20_L1R",      &HLT_Photon20_L1R_      ); 
  // tree->SetBranchAddress("HLT_Photon10_Cleaned_L1R",      &HLT_Photon10_Cleaned_L1R_  ); 
  // tree->SetBranchAddress("HLT_Photon15_Cleaned_L1R",      &HLT_Photon15_Cleaned_L1R_  ); 
  // tree->SetBranchAddress("HLT_Photon20_Cleaned_L1R",      &HLT_Photon20_Cleaned_L1R_  ); 

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
    // tree->SetBranchAddress("passm_ll_nom",           &passm_ll_nom_           );
    // tree->SetBranchAddress("passm_ll_nomttbar",      &passm_ll_nomttbar_      );  
    // tree->SetBranchAddress("passm_ll_nomttbarV2",    &passm_ll_nomttbarV2_    );  
    // tree->SetBranchAddress("passe_ll_ttbar",         &passe_ll_ttbar_         );
    // tree->SetBranchAddress("passe_ll_ttbarV1",       &passe_ll_ttbarV1_       );
    // tree->SetBranchAddress("passe_ll_ttbarV2",       &passe_ll_ttbarV2_       );
    // tree->SetBranchAddress("passe_ll_cand01",        &passe_ll_cand01_        );  
    // tree->SetBranchAddress("passm_lt_nom",           &passm_lt_nom_           );
    // tree->SetBranchAddress("passm_lt_nomttbar",      &passm_lt_nomttbar_      );  
    // tree->SetBranchAddress("passm_lt_nomttbarV2",    &passm_lt_nomttbarV2_    );  
    // tree->SetBranchAddress("passe_lt_ttbar",         &passe_lt_ttbar_         );
    // tree->SetBranchAddress("passe_lt_ttbarV1",       &passe_lt_ttbarV1_       );
    // tree->SetBranchAddress("passe_lt_ttbarV2",       &passe_lt_ttbarV2_       );
    // tree->SetBranchAddress("passe_lt_cand01",        &passe_lt_cand01_        );
    //tree->SetBranchAddress("TMLastStationTight_ll",  &TMLastStationTight_ll_  );
    //tree->SetBranchAddress("TMLastStationTight_lt",  &TMLastStationTight_lt_  );
  
    tree->SetBranchAddress("ptll",                 &ptll_                 );  
    tree->SetBranchAddress("ptlt",                 &ptlt_                 ); 
    tree->SetBranchAddress("ptllpf",               &ptllpf_               );  
    tree->SetBranchAddress("ptltpf",               &ptltpf_               ); 
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
    tree->SetBranchAddress("lep1",                 &lep1_                 );
    tree->SetBranchAddress("lep2",                 &lep2_                 );
    tree->SetBranchAddress("lep3",                 &lep3_                 );
    tree->SetBranchAddress("pflep1",               &pflep1_               );
    tree->SetBranchAddress("pflep2",               &pflep2_               );
    tree->SetBranchAddress("el1tv",                &el1tv_                );
    tree->SetBranchAddress("el2tv",                &el2tv_                );

    tree->SetBranchAddress("ee",                   &ee_                   );
    tree->SetBranchAddress("mm",                   &mm_                   );
    tree->SetBranchAddress("em",                   &em_                   );
    tree->SetBranchAddress("me",                   &me_                   );

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


