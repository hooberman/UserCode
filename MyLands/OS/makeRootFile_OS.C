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

void makeRootFile_OS(){

  TFile *f = TFile::Open("OS_shape.root","RECREATE");
  f->cd();

  //-----------------------------------------------------------
  // set yields and uncertainties
  //-----------------------------------------------------------

  // number of bins for shape analysis
  const unsigned int nbins                  = 3; 

  // yields: NOTE- all yields defined in EXCLUSIVE bins
  
  // signal regions                       R1                  R2              R3
  int   data_yield[nbins]          = {    11          ,        6       ,       5   };
  int   bkg_yield[nbins]           = {    7.4         ,      5.1       ,     7.2   };
  int   LM6_yield[nbins]           = {    7.          ,       6.       ,      7.   };

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

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_LM6               -> SetBinContent(ibin+1,LM6_yield[ibin]);
  }

  //-------------------------------------------------------------------------
  // histogram containing LM6 expected yield
  //-------------------------------------------------------------------------

  TH1F* histo_bkg               = new TH1F("histo_bkg","histo_bkg",nbins,0,nbins);

  for( unsigned int ibin = 0 ; ibin < nbins ; ibin++){
    histo_bkg               -> SetBinContent(ibin+1,bkg_yield[ibin]);
  }
  
  histo_Data->Write();
  histo_LM6->Write();
  histo_bkg->Write();

  f->Close();

}
