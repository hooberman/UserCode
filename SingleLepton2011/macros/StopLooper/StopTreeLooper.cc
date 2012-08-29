
#include "StopTreeLooper.h"
#include "Core/StopTree.h"
#include "Plotting/PlotUtilities.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"

#include <algorithm>
#include <utility>
#include <map>

float StopTreeLooper::dRbetweenVectors(LorentzVector vec1, 
				       LorentzVector vec2 ){ 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}


StopTreeLooper::StopTreeLooper()
{
    m_outfilename_ = "histos.root";
}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}

void StopTreeLooper::loop(TChain *chain, TString name)
{

    printf("[StopTreeLooper::loop] %s\n", name.Data());

    //
    // check for valid chain
    //

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        cout << "[StopTreeLooper::loop] no files in chain" << endl;
        return;
    }

    //
    // set up histograms
    //

    gROOT->cd();

    cout << "[StopTreeLooper::loop] setting up histos" << endl;

    //plotting map
    std::map<std::string, TH1F*> h_1d;

    //
    // file loop
    //

    unsigned int nEventsChain=0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    bool isData = name.Contains("data") ? true : false;

    cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

        //
        // load the stop baby tree
        //

        StopTree *tree = new StopTree();
        tree->LoadTree(currentFile->GetTitle());
        tree->InitTree();

        //
        // event loop
        //

        ULong64_t nEvents = tree->tree_->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) {
            tree->tree_->GetEntry(event);

            //
            // increment counters
            //

            ++nEventsTotal;
	    if (nEventsTotal%10000==0) {
	      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	      //if (i_permille != i_permille_old) {//this prints too often!
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            // 
            // event weight
            // 

	    float evtweight = isData ? 1. : ( tree->weight_ * 4.98 * tree->ndavtxweight_ * tree->mgcor_ );

            // 
            // selection criteria
            // 

	    //preselection
	    if ( !passSingleMuonSelection(tree) ) continue;

	    //tighten MET and MT for SRA
	    // if ( tree->t1metphicorr_   < 100 ) continue;
	    // if ( tree->t1metphicorrmt_ < 150 ) continue;
	    
	    //horrible code to attempt to reconstruct the hadronic top following the ATLAS procedure

	    LorentzVector jets[6];
	    jets[0] = tree->pfjet1_;
	    jets[1] = tree->pfjet2_;
	    jets[2] = tree->pfjet3_;
	    jets[3] = tree->pfjet4_;
	    jets[4] = tree->pfjet5_;
	    jets[5] = tree->pfjet6_;
	    
	    //loop over all jet pairs and find the one with the smallest dR with inv mass > 60 GeV

	    float mhadtop = 0.;
	    int i_j1 = -9;
	    int i_j2 = -9;
	    LorentzVector hadW = LorentzVector();
	    float mindRjj = 9999.;
	    for ( int i_jet=0; i_jet<6; ++i_jet ) {
	      if ( jets[i_jet].Pt()<20. ) continue;
	      for ( int j_jet=i_jet+1; j_jet<6; ++j_jet ) {
		if ( jets[j_jet].Pt()<20. ) continue;
		float jjmass = (jets[i_jet]+jets[j_jet]).mass();
		if (jjmass<60.) continue;
		float dRjj = dRbetweenVectors(jets[i_jet],jets[j_jet]);
		if (dRjj<mindRjj) {
		  hadW = jets[i_jet]+jets[j_jet];
		  mindRjj = dRjj;
		  i_j1 = i_jet;
		  i_j2 = j_jet;
		}
	      }
	    }

	    //look for third jet closest to hadronic W, with m(jjj) < 130 GeV

	    if (mindRjj<9998.) {
	      LorentzVector hadtop = LorentzVector();
	      float mindRWj = 9999.;
	      for ( int i_jet=0; i_jet<6; ++i_jet ) {
	    	if ( jets[i_jet].Pt()<20. ) continue;
	    	if (i_jet==i_j1 || i_jet==i_j2) continue;
	    	float Wjmass = (jets[i_jet]+hadW).mass();
	    	if (Wjmass<130.) continue;
	    	float dRWj = dRbetweenVectors(hadW, jets[i_jet]);
	    	if (dRWj<mindRWj) {
	    	  hadtop = jets[i_jet]+hadW;
	    	  mindRWj = dRWj;
	    	  mhadtop = hadtop.mass();
	    	}
	      }
	    } 

	    //store
	    plot1D("h_mhadtop",min(mhadtop, (float)499.99), evtweight, h_1d, 50, 0, 500.);

        } // end event loop

	// delete tree;

    } // end file loop

    //
    // finish
    //

    TFile outfile(m_outfilename_.c_str(),"RECREATE") ; 
    printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());

    std::map<std::string, TH1F*>::iterator it1d;
    for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
     it1d->second->Write(); 
     delete it1d->second;
    }

    outfile.Write();
    outfile.Close();
    
    gROOT->cd();

}

bool StopTreeLooper::passSingleMuonSelection(const StopTree *sTree) 
{
  //single lepton muon selection for 7 TeV muon analysis

  //rho requirement
  if ( sTree->rhovor_<0. || sTree->rhovor_>=40. ) return false;
  
  //4-jet requirement
  if ( sTree->npfjets30_<4) return false;

  //pass b-tag requirement
  if ( sTree->nbtagsssv_<1) return false;
  
  //at least one lepton
  if ( sTree->ngoodlep_ < 1 ) return false;

  //lepton flavor - ask for muon
  if ( sTree->leptype_ != 1 ) return false;

  //pass trigger if data - single muon
  if ( sTree->smu30_ != 1 )   return false;
  if ( sTree->trgmu30_ != 1 ) return false;

  //passes pt and eta requirements
  if ( sTree->lep1_.Pt() < 30 )          return false;
  if ( fabs(sTree->lep1_.Eta() ) > 2.1)  return false;
  
  //passes muon id selection criteria
  if ( sTree->lep1chi2ndf_>=10) return false;
  if ( sTree->lep1dpt_>=0.1)    return false;

  //passes met requirement
  if ( sTree->t1metphicorr_<50) return false;

  //pass isolated track veto
  //unfortunately changed default value to 9999.
  if ( sTree->pfcandpt10_ <9998. && sTree->pfcandiso10_ < 0.1 ) return false;

  return true;

}
