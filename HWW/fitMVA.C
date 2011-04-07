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

  
  //---------------------------------------
  // name of MVA branch in ntuple
  //---------------------------------------
  
  const char*   mvaname = "nn_hww160_ww";
  const int     nbins   =   20;
  const float   xmin    = -0.5;
  const float   xmax    =  1.5;

  //const char*   mvaname = "bdt_hww160_ww";
  //const int     nbins   =   30;
  //const float   xmin    = -1.0;
  //const float   xmax    =  0.5;


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
  //TCut lid("lid3 == 0");
  //TCut lowptb("jetLowBtag <= 2.1");
  //TCut softmu("nSoftMuons == 0");  
  //TCut sel160 = "lep1.pt()>30 && lep2.pt()>25 && dilep.mass()<50 && dPhi<1.05";
 
  TCut sel = jet + mll100 + OS + pt2010 + pmet;
  TCut weight = "0.5 * scale1fb";

  TCut selweight = sel*weight;

  //--------------------------------------------------
  // Get sig and bkg samples.
  // These samples will be used to __MAKE__ the PDFs.
  //--------------------------------------------------

  char* babyPath = "/smurf/benhoob/MVA/SmurfTraining/hww160_ww/output";

  TChain* sig = new TChain("tree");
  sig->Add(Form("%s/hww160.root",babyPath));

  TChain* bkg = new TChain("tree");
  bkg->Add(Form("%s/ww.root",babyPath));

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
  sig->Draw(Form("%s >> mva_sig",mvaname),selweight);
  bkg->Draw(Form("%s >> mva_bkg",mvaname),selweight);

  float nsigtrue = mva_sig->Integral();
  float nbkgtrue = mva_bkg->Integral();

  /*
  TCanvas *c1 = new TCanvas();
  c1->cd();
  mva_bkg->GetXaxis()->SetTitle(mvaname);
  mva_bkg->Draw("hist");
  mva_sig->SetLineColor(2);
  mva_sig->SetMarkerColor(2);
  mva_sig->Draw("sameE1");
  
  TLegend *leg = new TLegend(0.65,0.6,0.85,0.8);
  leg->AddEntry(mva_sig,"Signal","p");
  leg->AddEntry(mva_bkg,"Background");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();
  */

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

  const unsigned int nTrials = 1;
  
  float nsigfit[nTrials];
  float nbkgfit[nTrials];
  float nsigerrfit[nTrials];
  float nbkgerrfit[nTrials];
  float significance[nTrials];

  TH1F* hsig		= new TH1F("hsig","Sig Yield",100,0,100);
  TH1F* hbkg		= new TH1F("hbkg","Bkg Yield",100,0,500);
  TH1F* hsigpull	= new TH1F("hsigpull","Sig Pull",100,-5,5);
  TH1F* hbkgpull	= new TH1F("hbkgpull","Bkg Pull",100,-5,5);
  TH1F* hsignificance	= new TH1F("hsignificance","Significance",100,0,10);

  for( unsigned int i = 0 ; i < nTrials ; ++i ){

    nsig.setVal( nsigtrue );
    nbkg.setVal( nbkgtrue );

    RooDataSet *gendata = datapdf.generate(RooArgList(mva),nsigtrue+nbkgtrue,Extended(kTRUE));

    cout << "nsigtrue " << nsigtrue << endl;
    cout << "nbkgtrue " << nbkgtrue << endl;
    cout << "ntot     " << nsigtrue + nbkgtrue << endl;
    cout << "data     " << gendata->sumEntries() << endl;
    
    
    TCanvas *c4 = new TCanvas();
    c4->cd();
    
    RooPlot* frame4 = mva.frame();
    frame4->SetXTitle(mvaname);
    gendata->plotOn(frame4);
    //sigpdf.plotOn(frame4,Normalization(nsig.getVal()/gendata->sumEntries()));
    //bkgpdf.plotOn(frame4,LineColor(kRed),Normalization(nbkg.getVal()/gendata->sumEntries()));
    datapdf.plotOn(frame4,Components(sigpdf));
    datapdf.plotOn(frame4,Components(bkgpdf),LineColor(kRed));
    datapdf.plotOn(frame4,LineColor(kOrange));

    //datapdf.paramOn(frame4,Format("NEU",AutoPrecision(1)));
    frame4->Draw();
        
    nsig.setVal(0);
    nsig.setConstant();
    
    RooFitResult *bkg_result = datapdf.fitTo( *gendata , Save() , Extended(kTRUE) );
    
    nsig.setVal(50);
    nsig.setConstant(kFALSE);
    
    RooFitResult *result = datapdf.fitTo( *gendata , Save() , Extended(kTRUE) );
    
    /*
      TCanvas *c3 = new TCanvas();
      c3->cd();
      
      RooPlot* dataframe = mva.frame();
      gendata->plotOn(dataframe);
      sigpdf.plotOn(dataframe,LineStyle(kDashed),LineColor(kRed),Normalization(nsig.getVal()/gendata->sumEntries()));
      bkgpdf.plotOn(dataframe,LineStyle(kDashed),LineColor(kOrange),Normalization(nbkg.getVal()/gendata->sumEntries()));
      datapdf.plotOn(dataframe);
      dataframe->Draw();
    */
    
    //int nfloatpars = result->floatParsFinal().getSize();
    //int ndf        = nbins - nfloatpars;
    //float chi2     = dataframe->chiSquare(nfloatpars)*ndf;
    //float prob     = TMath::Prob(chi2,ndf);
    
    float nll     = result->minNll();
    float bkg_nll = bkg_result->minNll();
    float signif  = sqrt( -2 * ( nll - bkg_nll ) );

    cout << endl << endl;
    cout << "Fit Results-------------------------------------------"       << endl;
    //cout << "chi2/ndf   = " << chi2 << " / " << ndf                        << endl;
    //cout << "prob       = " << prob                                        << endl;  
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

  gStyle->SetOptFit(1111);

  TCanvas *c5 = new TCanvas("c5","",1200,800);
  c5->Divide(3,2);

  c5->cd(1);
  hsig->GetXaxis()->SetTitle("Sig Yield");
  hsig->Draw();
  hsig->Fit("gaus");

  c5->cd(2);
  hsigpull->GetXaxis()->SetTitle("Sig Pull");
  hsigpull->Draw();
  hsigpull->Fit("gaus");

  c5->cd(4);
  hbkg->GetXaxis()->SetTitle("Bkg Yield");
  hbkg->Draw();
  hbkg->Fit("gaus");

  c5->cd(5);
  hbkgpull->GetXaxis()->SetTitle("Bkg Pull");
  hbkgpull->Draw();
  hbkgpull->Fit("gaus");

  c5->cd(3);
  hsignificance->GetXaxis()->SetTitle("Significance");
  hsignificance->Draw();
  hsignificance->Fit("gaus");

  //c5->Print("plots/fitMVA.gif");
















  /*
  //--------------------------------------------------
  // now get "data" MVA output distributions
  // for now, using sig+bkg MC as "pseudo-data"
  //--------------------------------------------------

  TChain* data = new TChain("tree");
  data->Add(Form("%s/hww160.root" , babyPath));
  data->Add(Form("%s/ww.root"     , babyPath));

  TH1F* mva_data = new TH1F("mva_data","MVA Output for Data" , nbins , xmin , xmax );
  mva_data->Sumw2();

  data->Draw(Form("%s >> mva_data",mvaname),selodd);
  delete ctemp;

  RooDataHist datahist("datahist","Data PDF Hist" , RooArgSet(mva),mva_data);

  TCanvas *c2 = new TCanvas();
  c2->cd();

  RooPlot* frame = mva.frame();
  frame->SetXTitle(mvaname);
  sigpdf.plotOn(frame);
  bkgpdf.plotOn(frame,LineColor(kRed));
  datahist.plotOn(frame,LineColor(kOrange));
  frame->Draw();
  //gPad->SetLogy();
  */

  //-----------------------------------------------------------
  // construct PDF from sum of sig, bkg PDFs and fit to data
  //-----------------------------------------------------------


  /*
  RooFitResult *result = datapdf.fitTo( datahist , Save() , SumW2Error(kFALSE) , Extended(kTRUE) );

  TCanvas *c3 = new TCanvas();
  c3->cd();

  RooPlot* dataframe = mva.frame();
  datahist.plotOn(dataframe);
  datapdf.plotOn(dataframe);
  sigpdf.plotOn(dataframe,LineStyle(kDashed),LineColor(kRed),Normalization(nsig.getVal()/mva_data->Integral()));
  bkgpdf.plotOn(dataframe,LineStyle(kDashed),LineColor(kOrange),Normalization(nbkg.getVal()/mva_data->Integral()));
  datahist.plotOn(dataframe);
  datapdf.plotOn(dataframe);
  dataframe->Draw();

  int nfloatpars = result->floatParsFinal().getSize();
  int ndf        = nbins - nfloatpars;
  float chi2     = dataframe->chiSquare(nfloatpars)*ndf;
  float prob     = TMath::Prob(chi2,ndf);
 
  cout << endl << endl;
  cout << "Fit Results-------------------------------------------"       << endl;
  cout << "chi2/ndf   = " << chi2 << " / " << ndf                        << endl;
  cout << "prob       = " << prob                                        << endl;  
  cout << "nll        = " << result->minNll()                            << endl;
  cout << "nsig       = " << nsig.getVal() << " +/- " << nsig.getError() << endl;
  cout << "nsig(true) = " << nsigtrue                                    << endl;
  cout << "nbkg       = " << nbkg.getVal() << " +/- " << nbkg.getError() << endl;
  cout << "nbkg(true) = " << nbkgtrue                                    << endl;
  cout << "------------------------------------------------------"       << endl;
  cout << endl << endl;
  
  */








}
