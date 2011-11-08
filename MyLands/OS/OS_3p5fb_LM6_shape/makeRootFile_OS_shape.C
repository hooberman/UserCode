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

void makeRootFile_OS_shape(){

  TFile *f = TFile::Open("OS_shape.root","RECREATE");
  f->cd();

  //-----------------------------------------------------------
  // set yields and uncertainties
  //-----------------------------------------------------------

  // number of bins for shape analysis
  const unsigned int nbins                  = 6; 

  // yields: NOTE- all yields defined in EXCLUSIVE bins
  
  // signal regions                       R1(SF)          R1(OF)    R2(SF)     R2(OF)   R3(SF)       R3(OF)
  int     data_yield[nbins]           = {      4        ,   9     ,  6      ,    4     ,  3        ,   9      };
  float   bkg_yield[nbins]            = {    6.83/2.    , 6.83/2. , 8.24/2. ,  8.24/2. , 10.43/2.  , 10.43/2. };
  float   LM6_yield[nbins]            = {     2.9       ,   2.8   ,   13.4  ,   8.1    ,   2.5     ,   0.4    };
  float   LM6_yield_JESup[nbins]      = {  3.4   ,   2.5   ,   14.7   ,   9.3   ,   2.1   ,   0.7  };
  float   LM6_yield_JESdn[nbins]      = {  2.7   ,   3.0   ,   11.7   ,   6.3   ,   2.7   ,   1.2  };
  float   bkg_syst[nbins]             = {    3.65/2.    ,       3.65/2. ,      3.72/2. ,      3.72/2. ,    4.04/2.   ,     4.04/2.  };
  float   bkg_stat[nbins]             = { 5.95/sqrt(2.) , 5.95/sqrt(2.) , 4.92/sqrt(2) , 4.92/sqrt(2) , 4.89/sqrt(2) , 4.89/sqrt(2) };

  //---------------------------------------------------------------------------------------------
  // EVERYTHING BELOW HERE IS JUST MAKING HISTOGRAMS FROM THE ABOVE YIELDS
  // NO NEED TO ALTER ANYTHING BELOW
  //---------------------------------------------------------------------------------------------

  //------------------------------------------------------
  // histogram containing observed yields in data
  //------------------------------------------------------
  
  TH1F* histo_Data = new TH1F("histo_Data","histo_Data",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_Data->SetBinContent( ibin+1 , data_yield[ibin] );
  }

  //-------------------------------------------------------------------------
  // histogram containing LM6 expected yield
  //-------------------------------------------------------------------------

  TH1F* histo_LM6               = new TH1F("histo_LM6","histo_LM6",nbins,0,nbins);
  TH1F* histo_LM6_JES_shapeUp   = new TH1F("histo_LM6_JES_shapeUp"  ,"histo_LM6_JES_shapeUp"  ,nbins,0,nbins);
  TH1F* histo_LM6_JES_shapeDown = new TH1F("histo_LM6_JES_shapeDown","histo_LM6_JES_shapeDown",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_LM6                  -> SetBinContent(ibin+1,LM6_yield[ibin]);
    histo_LM6_JES_shapeUp      -> SetBinContent(ibin+1,LM6_yield_JESup[ibin]);
    histo_LM6_JES_shapeDown    -> SetBinContent(ibin+1,LM6_yield_JESdn[ibin]);
  }

  //-------------------------------------------------------------------------
  // histogram containing LM6 expected yield
  //-------------------------------------------------------------------------

  TH1F* histo_bkg               = new TH1F("histo_bkg","histo_bkg",nbins,0,nbins);
  TH1F* histo_bkg_statUp        = new TH1F("histo_bkg_statUp"  ,"histo_bkg_statUp"  ,nbins,0,nbins);
  TH1F* histo_bkg_statDown      = new TH1F("histo_bkg_statDown","histo_bkg_statDown",nbins,0,nbins);
  TH1F* histo_bkg_systUp        = new TH1F("histo_bkg_systUp"  ,"histo_bkg_systUp"  ,nbins,0,nbins);
  TH1F* histo_bkg_systDown      = new TH1F("histo_bkg_systDown","histo_bkg_systDown",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_bkg               -> SetBinContent(ibin+1,bkg_yield[ibin]);
    histo_bkg_statUp        -> SetBinContent(ibin+1, bkg_yield[ibin] + bkg_stat[ibin] );
    histo_bkg_statDown      -> SetBinContent(ibin+1, bkg_yield[ibin] - bkg_stat[ibin] );
    histo_bkg_systUp        -> SetBinContent(ibin+1, bkg_yield[ibin] + bkg_syst[ibin] );
    histo_bkg_systDown      -> SetBinContent(ibin+1, bkg_yield[ibin] - bkg_syst[ibin] );
  }

  histo_Data->Write();
  histo_LM6->Write();
  histo_LM6_JES_shapeUp->Write();
  histo_LM6_JES_shapeDown->Write();
  histo_bkg->Write();
  histo_bkg_statUp->Write();
  histo_bkg_statDown->Write();
  histo_bkg_systUp->Write();
  histo_bkg_systDown->Write();

  f->Close();

}
