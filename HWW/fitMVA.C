#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TMath.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TStyle.h"
#include "TLegend.h"
#include "RooGaussian.h"
#include "RooNumConvPdf.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"

using namespace std;
using namespace RooFit;


void fitMVA(){

  //-------------------------------
  // user parameters
  //-------------------------------
  
  const int  nTrials  = 1000;  // number of trials in fit validation
  const bool display  = false; // display fit? Set to false for large nTrials
  const bool printgif = true; // print canvases to gif files?
  
  //---------------------------------------
  // name of MVA branch in ntuple
  //---------------------------------------
  
  //const char*   mvaname = "nn_hww160_ww";
  const char*   mvaname = "bdt_hww160_ww";

  //-----------------------------------------
  // selection to apply to sig, bkg samples
  //-----------------------------------------
  
  TCut jet("njets == 0");  
  TCut mll100("dilep->mass() < 100");
  TCut OS("lq1*lq2 < 0");
  TCut pt2010("lep1->pt() >= 20 && lep2->pt() >= 10");
  TCut pmet1("pmet >= 20 && (type == 1 || type == 2)");	 
  TCut pmet2("pmet >= 35 && (type == 0 || type == 3)");
  TCut pmet = pmet1 || pmet2;
  TCut lid("lid3 == 0");
  TCut lowptb("jetLowBtag <= 2.1");
  TCut softmu("nSoftMuons == 0");  
  TCut sel160("lep1.pt()>30 && lep2.pt()>25 && dilep.mass()<50 && dPhi<1.05");
 
  TCut sel       = jet + mll100 + OS + pt2010 + pmet;
  TCut weight    = "0.5 * scale1fb";
  TCut selweight = sel*weight;

  //--------------------------------------------------
  // get sig and bkg samples
  //--------------------------------------------------

  char* babyPath = "/smurf/benhoob/MVA/SmurfTraining/hww160_ww/output";

  TChain* sig = new TChain("tree");
  sig->Add(Form("%s/hww160.root",babyPath));

  TChain* bkg = new TChain("tree");
  bkg->Add(Form("%s/ww.root",babyPath));

  //------------------------------------
  // set binning
  //------------------------------------

  int   nbins = 0;
  float xmin  = 0.;
  float xmax  = 0.;

  if( strcmp( mvaname , "nn_hww160_ww" ) == 0 ){
    nbins   =   20;
    xmin    = -0.5;
    xmax    =  1.5;
  }
  else if( strcmp( mvaname , "bdt_hww160_ww" ) == 0 ){
    nbins   =   30;
    xmin    = -1.0;
    xmax    =  0.5;
  }else{
    cout << "Unrecognized MVA " << mvaname << ", quitting" << endl;
    exit(0);
  }
  
  //------------------------
  // declare variables
  //------------------------

  RooRealVar nsig ("nsig"   ,  "Signal Yield"       ,   50 ,    0 , 1000 );
  RooRealVar nbkg ("nbkg"   ,  "Background Yield"   ,  300 ,    0 , 1000 );
  RooRealVar mva  ("mva"    ,  "MVA Output"                , xmin , xmax );

  nsig.setVal(50);
  nbkg.setVal(300);
  mva.setBins(nbins);

  //---------------------------------------------------
  // get MVA output distributions from sig, bkg samp
  //---------------------------------------------------
  
  TH1F* mva_sig = new TH1F("mva_sig","MVA Output for Signal"     , nbins , xmin , xmax );
  TH1F* mva_bkg = new TH1F("mva_bkg","MVA Output for Background" , nbins , xmin , xmax );

  mva_sig->Sumw2();
  mva_bkg->Sumw2();

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  sig->Draw(Form("%s >> mva_sig",mvaname),selweight);
  bkg->Draw(Form("%s >> mva_bkg",mvaname),selweight);
  delete ctemp;

  float nsigtrue = mva_sig->Integral();
  float nbkgtrue = mva_bkg->Integral();

  //---------------------------------------------------
  // convert MVA distributions to RooHistPdf objects
  //---------------------------------------------------

  RooDataHist sigpdfhist("sigpdfhist","Signal PDF Hist"     , RooArgSet(mva),mva_sig);
  RooDataHist bkgpdfhist("bkgpdfhist","Background PDF Hist" , RooArgSet(mva),mva_bkg);

  RooHistPdf sigpdf("sigpdf", "Signal PDF"     , RooArgSet(mva), sigpdfhist);
  RooHistPdf bkgpdf("bkgpdf", "Background PDF" , RooArgSet(mva), bkgpdfhist);

  nsig.setVal( nsigtrue );
  nbkg.setVal( nbkgtrue );

  RooAddPdf datapdf("datapdf", "Data PDF", RooArgList(sigpdf,bkgpdf), RooArgList(nsig,nbkg));

  //-------------------------------------------------
  // perform fit (nTrials iterations)
  //-------------------------------------------------

  float nsigfit[nTrials];
  float nbkgfit[nTrials];
  float nsigerrfit[nTrials];
  float nbkgerrfit[nTrials];
  float significance[nTrials];
  TCanvas *can[nTrials];

  TH1F* hsig		= new TH1F("hsig","Sig Yield",100,0,100);
  TH1F* hbkg		= new TH1F("hbkg","Bkg Yield",100,0,200);
  TH1F* hsigpull	= new TH1F("hsigpull","Sig Pull",100,-5,5);
  TH1F* hbkgpull	= new TH1F("hbkgpull","Bkg Pull",100,-5,5);
  TH1F* hsignificance	= new TH1F("hsignificance","Significance",100,0,10);

  for( int i = 0 ; i < nTrials ; ++i ){

    nsig.setVal( nsigtrue );
    nbkg.setVal( nbkgtrue );

    //-------------------------------------
    // generate pseudo-dataset from PDF
    //-------------------------------------

    RooDataSet *gendata = datapdf.generate(RooArgList(mva),nsigtrue+nbkgtrue,Extended(kTRUE));

    //-------------------------------------
    // perform fit with nsig fixed to 0
    //-------------------------------------

    nsig.setVal(0);
    nbkg.setVal(gendata->sumEntries());
    nsig.setConstant();

    RooAbsReal* mynll_bkg = bkgpdf.createNLL( *gendata );
    RooAbsReal* mynll_tot = datapdf.createNLL( *gendata );
    
    cout << "nll_bkg  " << mynll_bkg->getVal() << endl;
    cout << "nll_tot  " << mynll_tot->getVal() << endl;

    cout << endl;
    cout << "-------------------------------------------" << endl;
    cout << "Performing fit with signal yield fixed to 0" << endl;
    cout << "-------------------------------------------" << endl;
    cout << endl;

    RooFitResult *bkg_result = datapdf.fitTo( *gendata , Save() , Extended(kTRUE) );

    /*
    TCanvas *bkgcan = new TCanvas("bkgcan","bkgcan",600,600);
    bkgcan->cd();

    RooPlot* bkgframe = mva.frame();
    bkgframe->SetXTitle(mvaname);
    gendata->plotOn(bkgframe);
    datapdf.plotOn(bkgframe,Components(sigpdf));
    datapdf.plotOn(bkgframe,Components(bkgpdf),LineColor(kRed));
    datapdf.plotOn(bkgframe,LineColor(kOrange));
    bkgframe->Draw();
    */

    //-------------------------------------
    // perform fit with nsig floated in fit
    //-------------------------------------
    
    //nsig.setVal(50);
    nsig.setConstant(kFALSE);

    cout << endl;
    cout << "----------------------------------------" << endl;
    cout << "Performing fit with signal yield floated" << endl;
    cout << "----------------------------------------" << endl;
    cout << endl;

    RooFitResult *result = datapdf.fitTo( *gendata , Save() , Extended(kTRUE) );
    
    if( display ){

      can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),600,600);
      can[i]->cd();
    
      RooPlot* frame = mva.frame();
      frame->SetXTitle(mvaname);
      gendata->plotOn(frame);
      datapdf.plotOn(frame,Components(sigpdf));
      datapdf.plotOn(frame,Components(bkgpdf),LineColor(kRed));
      datapdf.plotOn(frame,LineColor(kOrange));
      frame->Draw();

      if( printgif ) can[i]->Print(Form("plots/%s_%i.gif",mvaname,i));
    }

    //---------------------------------------------
    // print output to screen
    //---------------------------------------------

    float nll     = result->minNll();
    float bkg_nll = bkg_result->minNll();
    float signif  = sqrt( -2 * ( nll - bkg_nll ) );

    cout << endl << endl;
    cout << "Fit Results-------------------------------------------"       << endl;
    cout << "gen events = " << gendata->sumEntries()                       << endl;
    cout << "nll (BKG)  = " << bkg_nll                                     << endl;
    cout << "nll        = " << nll                                         << endl;
    cout << "sig        = " << signif                                      << endl;
    cout << "nsig       = " << nsig.getVal() << " +/- " << nsig.getError() << endl;
    cout << "nsig(true) = " << nsigtrue                                    << endl;
    cout << "nbkg       = " << nbkg.getVal() << " +/- " << nbkg.getError() << endl;
    cout << "nbkg(true) = " << nbkgtrue                                    << endl;
    cout << "------------------------------------------------------"       << endl;
    cout << endl << endl;
    
    
    nsigfit[i]       = nsig.getVal();
    nsigerrfit[i]    = nsig.getError();
    nbkgfit[i]       = nbkg.getVal();
    nbkgerrfit[i]    = nbkg.getError();
    significance[i]  = signif;

    hsig->Fill( nsigfit[i] );
    hsigpull->Fill( ( nsigfit[i] - nsigtrue ) / nsigerrfit[i] );
    hbkg->Fill( nbkgfit[i] );
    hbkgpull->Fill( ( nbkgfit[i] - nbkgtrue ) / nbkgerrfit[i] );
    hsignificance->Fill( significance[i] );
  }

  gStyle->SetOptFit(0111);

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(3,2);

  c1->cd(1);
  hsig->GetXaxis()->SetTitle("Sig Yield");
  hsig->Draw();
  hsig->Fit("gaus");

  c1->cd(2);
  hsigpull->GetXaxis()->SetTitle("Sig Pull");
  hsigpull->Draw();
  hsigpull->Fit("gaus");

  c1->cd(4);
  hbkg->GetXaxis()->SetTitle("Bkg Yield");
  hbkg->Draw();
  hbkg->Fit("gaus");

  c1->cd(5);
  hbkgpull->GetXaxis()->SetTitle("Bkg Pull");
  hbkgpull->Draw();
  hbkgpull->Fit("gaus");

  c1->cd(3);
  hsignificance->GetXaxis()->SetTitle("Significance");
  hsignificance->Draw();
  hsignificance->Fit("gaus");

  if( printgif ) c1->Print(Form("plots/%s_fitval.gif",mvaname));


}
