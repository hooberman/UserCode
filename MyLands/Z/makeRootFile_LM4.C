#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <sstream>

using namespace std;

void makeRootFile_LM4(){

  TFile *f = TFile::Open("Z_shape_LM4.root","RECREATE");
  f->cd();

  //-----------------------------------------------------------
  // set yields and uncertainties
  //-----------------------------------------------------------

  // number of bins for shape analysis
  const unsigned int nbins                  = 4; 

  //------------------------------------------------------------------------------------------------------
  // MET bins                                    30-60            60-100          100-200          >200  
  //------------------------------------------------------------------------------------------------------

  // yields: NOTE- ALL YIELD DEFINED IN EXCLUSIVE BINS!!!!!
  int   data_yield[nbins]          = {        2287-206  ,         206-57  ,          57-4  ,          4 };
  float LM4_yield[nbins]           = {       25.4-22.9  ,      22.9-20.1  ,     20.1-12.3  ,       12.3 };
  float OF_yield[nbins]            = {     246.6-152.5  ,   152.5 - 50.6  ,    50.6 - 3.2  ,        3.2 };
  float Z_yield[nbins]             = {   2060.3 - 60.8  ,     60.8 - 5.1  ,    5.1 - 0.09  ,       0.09 };

  // uncertainties
  float LM4_JES[nbins]             = {           0.000  ,          0.000  ,         0.040  ,      0.120 };
  float OF_stat[nbins]             = {           0.026  ,          0.032  ,         0.055  ,      0.220 };
  float OF_syst[nbins]             = {           0.090  ,          0.090  ,         0.090  ,      0.090 };
  float  Z_stat[nbins]             = {           0.014  ,          0.067  ,         0.200  ,      0.440 };
  float  Z_syst[nbins]             = {           0.150  ,          0.150  ,         0.150  ,      0.150 };

  //---------------------------------------------------------------------------------------------
  // EVERYTHING BELOW HERE IS JUST MAKING HISTOGRAMS FROM THE ABOVE YIELDS
  // NO NEED TO ALTER ANYTHING BELOW, UNLESS YOU WANT TO ALTER THE UNCERTAINTIES FOR EACH SOURCE
  //---------------------------------------------------------------------------------------------

  //------------------------------------------------------
  // histogram containing observed yields in data
  //------------------------------------------------------
  
  TH1F* histo_Data = new TH1F("histo_Data","histo_Data",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_Data->SetBinContent( ibin+1 , data_yield[ibin] );
  }

  //-------------------------------------------------------------------------
  // histogram containing LM4 expected yield and shapes for JES up, JES down
  //-------------------------------------------------------------------------

  TH1F* histo_LM4               = new TH1F("histo_LM4","histo_LM4",nbins,0,nbins);
  TH1F* histo_LM4_JES_shapeUp   = new TH1F("histo_LM4_JES_shapeUp"  ,"histo_LM4_JES_shapeUp"  ,nbins,0,nbins);
  TH1F* histo_LM4_JES_shapeDown = new TH1F("histo_LM4_JES_shapeDown","histo_LM4_JES_shapeDown",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_LM4               -> SetBinContent(ibin+1,LM4_yield[ibin]);
    histo_LM4_JES_shapeUp   -> SetBinContent(ibin+1,(1 + LM4_JES[ibin]) * LM4_yield[ibin]);
    histo_LM4_JES_shapeDown -> SetBinContent(ibin+1,(1 - LM4_JES[ibin]) * LM4_yield[ibin]);
  }

  //-------------------------------------------------------------------------
  // histogram containing OF expected yield and shapes for JES up, JES down
  //-------------------------------------------------------------------------

  TH1F* histo_OF               = new TH1F("histo_OF","histo_OF",nbins,0,nbins);
  TH1F* histo_OF_statUp        = new TH1F("histo_OF_statUp"  ,"histo_OF_statUp"  ,nbins,0,nbins);
  TH1F* histo_OF_statDown      = new TH1F("histo_OF_statDown","histo_OF_statDown",nbins,0,nbins);
  TH1F* histo_OF_systUp        = new TH1F("histo_OF_systUp"  ,"histo_OF_systUp"  ,nbins,0,nbins);
  TH1F* histo_OF_systDown      = new TH1F("histo_OF_systDown","histo_OF_systDown",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_OF               -> SetBinContent(ibin+1,OF_yield[ibin]);
    histo_OF_statUp        -> SetBinContent(ibin+1, (1 + OF_stat[ibin]) * OF_yield[ibin] );
    histo_OF_statDown      -> SetBinContent(ibin+1, (1 - OF_stat[ibin]) * OF_yield[ibin] );
    histo_OF_systUp        -> SetBinContent(ibin+1, (1 + OF_syst[ibin]) * OF_yield[ibin] );
    histo_OF_systDown      -> SetBinContent(ibin+1, (1 - OF_syst[ibin]) * OF_yield[ibin] );
  }

  //-------------------------------------------------------------------------
  // histogram containing Z expected yield and shapes for JES up, JES down
  //-------------------------------------------------------------------------

  TH1F* histo_Z               = new TH1F("histo_Z","histo_Z",nbins,0,nbins);
  TH1F* histo_Z_statUp        = new TH1F("histo_Z_statUp"  ,"histo_Z_statUp"  ,nbins,0,nbins);
  TH1F* histo_Z_statDown      = new TH1F("histo_Z_statDown","histo_Z_statDown",nbins,0,nbins);
  TH1F* histo_Z_systUp        = new TH1F("histo_Z_systUp"  ,"histo_Z_systUp"  ,nbins,0,nbins);
  TH1F* histo_Z_systDown      = new TH1F("histo_Z_systDown","histo_Z_systDown",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_Z               -> SetBinContent(ibin+1,Z_yield[ibin]);
    histo_Z_statUp        -> SetBinContent(ibin+1, (1 + Z_stat[ibin]) * Z_yield[ibin] );
    histo_Z_statDown      -> SetBinContent(ibin+1, (1 - Z_stat[ibin]) * Z_yield[ibin] );
    histo_Z_systUp        -> SetBinContent(ibin+1, (1 + Z_syst[ibin]) * Z_yield[ibin] );
    histo_Z_systDown      -> SetBinContent(ibin+1, (1 - Z_syst[ibin]) * Z_yield[ibin] );
  }

  //---------------------------------------------
  // write histos to root file, then we're done
  //---------------------------------------------
  
  histo_Data->Write();
  histo_LM4->Write();
  histo_LM4_JES_shapeUp->Write();
  histo_LM4_JES_shapeDown->Write();
  histo_OF->Write();
  histo_OF_statUp->Write();
  histo_OF_statDown->Write();
  histo_OF_systUp->Write();
  histo_OF_systDown->Write();
  histo_Z->Write();
  histo_Z_statUp->Write();
  histo_Z_statDown->Write();
  histo_Z_systUp->Write();
  histo_Z_systDown->Write();

  f->Close();

}
