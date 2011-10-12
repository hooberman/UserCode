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

  // yields: NOTE- all yields defined in EXCLUSIVE bins
  
  // MET bins                                30-60         60-100              100-200          >200
  int   data_yield[nbins]          = {    2287-206 ,       206-57 ,               57-4 ,           4 };

  float LM4_yield[nbins]           = {   25.4-22.9 ,    22.9-20.1 ,         20.1-12.3  ,         12.3 };
  float LM4_yield_JESup[nbins]     = {   25.4-22.9 ,    22.9-20.1 , 1.04 * (20.1-12.3) ,  1.12 * 12.3 }; // assume 0% uncertainty for MET30 & MET60 (needs to be updated)
  float LM4_yield_JESdown[nbins]   = {   25.4-22.9 ,    22.9-20.1 , 0.96 * (20.1-12.3) ,  0.88 * 12.3 }; // 4% and 12% for MET100 and MET200

  float OF_yield[nbins]            = {            246.6-152.5  ,            152.5 - 50.6  ,            50.6 - 3.2  ,          3.2 };
  float OF_yield_statup[nbins]     = { (1+0.026)*(246.6-152.5) , (1+0.032)*(152.5 - 50.6) , (1+0.055)*(50.6 - 3.2) , (1+0.22)*3.2 }; // OF histo stat up
  float OF_yield_statdown[nbins]   = { (1-0.026)*(246.6-152.5) , (1-0.032)*(152.5 - 50.6) , (1-0.055)*(50.6 - 3.2) , (1-0.22)*3.2 }; // OF histo stat down
  float OF_yield_systup[nbins]     = { (1+0.090)*(246.6-152.5) , (1+0.090)*(152.5 - 50.6) , (1+0.090)*(50.6 - 3.2) , (1+0.09)*3.2 }; // OF histo syst up
  float OF_yield_systdown[nbins]   = { (1-0.090)*(246.6-152.5) , (1-0.090)*(152.5 - 50.6) , (1-0.090)*(50.6 - 3.2) , (1-0.09)*3.2 }; // OF histo syst down
  
  float Z_yield[nbins]             = {            2060.3 - 60.8  ,            60.8 - 5.1  ,           5.1 - 0.09  ,          0.09 };
  float Z_yield_statup[nbins]      = { (1+0.014)*(2060.3 - 60.8) , (1+0.067)*(60.8 - 5.1) , (1+0.20)*(5.1 - 0.09) , (1+0.44)*0.09 };
  float Z_yield_statdown[nbins]    = { (1-0.014)*(2060.3 - 60.8) , (1-0.067)*(60.8 - 5.1) , (1-0.20)*(5.1 - 0.09) , (1-0.44)*0.09 };
  float Z_yield_systup[nbins]      = {  (1+0.15)*(2060.3 - 60.8) ,  (1+0.15)*(60.8 - 5.1) , (1+0.15)*(5.1 - 0.09) , (1+0.15)*0.09 };
  float Z_yield_systdown[nbins]    = {  (1-0.15)*(2060.3 - 60.8) ,  (1-0.15)*(60.8 - 5.1) , (1-0.15)*(5.1 - 0.09) , (1-0.15)*0.09 };

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
    histo_LM4_JES_shapeUp   -> SetBinContent(ibin+1,LM4_yield_JESup[ibin]);
    histo_LM4_JES_shapeDown -> SetBinContent(ibin+1,LM4_yield_JESdown[ibin]);
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
    histo_OF_statUp        -> SetBinContent(ibin+1,OF_yield_statup[ibin]);
    histo_OF_statDown      -> SetBinContent(ibin+1,OF_yield_statdown[ibin]);
    histo_OF_systUp        -> SetBinContent(ibin+1,OF_yield_systup[ibin]);
    histo_OF_systDown      -> SetBinContent(ibin+1,OF_yield_systdown[ibin]);
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
    histo_Z_statUp        -> SetBinContent(ibin+1,Z_yield_statup[ibin]);
    histo_Z_statDown      -> SetBinContent(ibin+1,Z_yield_statdown[ibin]);
    histo_Z_systUp        -> SetBinContent(ibin+1,Z_yield_systup[ibin]);
    histo_Z_systDown      -> SetBinContent(ibin+1,Z_yield_systdown[ibin]);
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
