#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
//#include "histtools.h"
#include <iostream>
#include <fstream>
#endif

using namespace std;

TH2* suppressHist(TH2* hist,int iclone,float xmin,float xmax){
  
  TH2F* h=new TH2F(Form("%s_%s%i",hist->GetName(),"clone",iclone),hist->GetTitle(),
		   hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),
		   hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
  
  for(int ibinx=1;ibinx<=hist->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=hist->GetNbinsY();ibiny++){
      if(hist->GetBinCenter(ibinx)>xmin && hist->GetBinCenter(ibinx)<xmax){
	h->SetBinContent(ibinx,ibiny,hist->GetBinContent(ibinx,ibiny));
      }else{
	h->SetBinContent(ibinx,ibiny,0);
      }
    }
  }
  return h;
}

TGraphErrors *getTGraphFromTH2(TH2F* h,vector<float> xbins, int method){

  gROOT->LoadMacro("scripts/suppressHist.C");

  static const int nbins = (int)xbins.size()-1;
  TH2* htemp[nbins];
  TH1* hx[nbins];
  TH1* hy[nbins];
    
  float x[nbins];
  float y[nbins];
  float xerr[nbins];
  float yerr[nbins];
  
  for(int ibin=0;ibin<nbins;ibin++){
    
    htemp[ibin] = suppressHist(h,ibin,xbins.at(ibin),xbins.at(ibin+1));
    
    //float width = xbins.at(ibin+1) - xbins.at(ibin);
    hx[ibin]    = htemp[ibin]->ProjectionX();
    hy[ibin]    = htemp[ibin]->ProjectionY();

    hx[ibin] -> StatOverflows(kFALSE);
    hy[ibin] -> StatOverflows(kFALSE);

    
    x[ibin]     = hx[ibin]->GetMean(1);
    //xerr[ibin]  = width/2.;
    xerr[ibin]  = (hx[ibin]->GetEntries()>0) ? 
      hx[ibin]->GetRMS(1)/sqrt(hx[ibin]->GetEntries()) : 0.;
    //xerr[ibin]  = hx[ibin]->GetRMS(1);
    
    if(method == 0){
      y[ibin]     = hy[ibin]->GetMean(1);
      yerr[ibin]  = (hy[ibin]->GetEntries()>0) ? 
	hy[ibin]->GetRMS(1)/sqrt(hy[ibin]->GetEntries()) : 0.;
      //yerr[ibin]  = hy[ibin]->GetRMS(1);
    }
    if(method == 1){
      y[ibin]     = hy[ibin]->GetRMS(1);
      yerr[ibin]  = 0.;
    }

    //cout<<"bin "<<ibin<<" x "<<x[ibin]<<" y "<<y[ibin]<<" xerr "<<xerr[ibin]<<" yerr "<<yerr[ibin]<<endl;
  }
 

  TGraphErrors *g=new TGraphErrors(nbins,x,y,xerr,yerr);
  g->GetXaxis()->SetLimits(xbins.at(0),xbins.at(nbins));
  
  return g;
}


