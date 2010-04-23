#ifndef __CINT__
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TProfile.h"
#include <sstream>
#include <vector>
#endif

TGraphErrors *getTGraphFromTH2(TH2F* h,vector<float> xbins,int color=2);
TH2* suppressHist(TH2* hist,int iclone,float xmin,float xmax);
inline double fround(double n, unsigned d);
TGraphErrors* diffTGraph(TGraphErrors* g1, TGraphErrors *g2, 
                         char* title = "DECO - PEAK", char* xtitle = "Layer", char* ytitle = "");

void corrections(bool printgif = false){
  
  //load macros and set style-----------------------------------
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  int colors[]       = {2,4,6,1};
  
  //files to process
  vector<char*> filenames;
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_ALLPEAK/merged.root");
  filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_ALLPEAK/merged.root");
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_copy/merged_2.root");
  filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8/merged.root");
  filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_toblatebp06/merged.root");
  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom_mintrkmom1/merged.root");
  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom/merged.root");
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_toblatebp063/merged.root");
  
  //labels for files
  vector<char*> types;
  types.push_back("peak");
  types.push_back("deco");
  types.push_back("deco (BP)");  
  //types.push_back("MC");
  
  assert ( filenames.size() == types.size() );

  //list of subdetectors
  vector<char*> subdets;
  subdets.push_back("TOB");
  subdets.push_back("TIB");

  const unsigned int nfiles = filenames.size();
  const unsigned int ndet = subdets.size();

  TFile*   file[nfiles];
  TH2F*    hduthp[nfiles][ndet];
  TH2F*    hduthm[nfiles][ndet];
  TCanvas* can[ndet];

  //dummy legend
  TLegend *leg1=new TLegend(0.15,0.65,0.35,0.85);

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
    TH1F *hdummy=new TH1F(Form("hdummy_%i",ifile),Form("hdummy_%i",ifile),1,0,1);
    hdummy->SetLineColor(colors[ifile]);
    hdummy->SetMarkerColor(colors[ifile]);
    leg1->AddEntry(hdummy,types.at(ifile));
  }

  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);

  cout << "Getting histos... " << endl;

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
    cout << "Opening " << types.at(ifile) <<" file   " << filenames.at(ifile) << endl;
    file[ifile] = TFile::Open(filenames.at(ifile));
    
    for(unsigned int idet = 0 ; idet < ndet ; idet++){
    
      if( strcmp(types.at(ifile),"deco (BP)") == 0){

        hduthp[ifile][idet] = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vp_wp_all",subdets.at(idet)));
        TH2F* htp     = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vp_wm_all",subdets.at(idet)));
        hduthp[ifile][idet] -> Add(htp);
        
        hduthm[ifile][idet] = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vm_wp_all",subdets.at(idet)));
        TH2F* htm     = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vm_wm_all",subdets.at(idet)));
        hduthm[ifile][idet] -> Add(htm);

      }else{
        
        hduthp[ifile][idet] = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vp_wp",subdets.at(idet)));
        TH2F* htp     = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vp_wm",subdets.at(idet)));
        hduthp[ifile][idet] -> Add(htp);
        
        hduthm[ifile][idet] = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vm_wp",subdets.at(idet)));
        TH2F* htm     = (TH2F*) file[ifile]->Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vm_wm",subdets.at(idet)));
        hduthm[ifile][idet] -> Add(htm);
        
      }
    }
  }
  
  vector<float> thetabinsp;
  vector<float> thetabinsm;



  for(unsigned int idet = 0 ; idet < ndet ; idet++){

    cout << "Plotting " << subdets.at(idet) << "..." << endl;

    can[idet]=new TCanvas(Form("%s_can",subdets.at(idet)),Form("%s canvas",subdets.at(idet)),1200,900);
    can[idet]->Divide(2,2);
  
    thetabinsp.clear();
    thetabinsm.clear();

    if( strcmp( subdets.at(idet) , "TOB" ) == 0 ){
      for(int ibin=14;ibin<=30;ibin++)    thetabinsp.push_back(ibin*0.05-1);
      for(int ibin=10;ibin<=26;ibin++)    thetabinsm.push_back(ibin*0.05-1);
    }   
    if( strcmp( subdets.at(idet) , "TIB" ) == 0 ){
      for(int ibin=14;ibin<=26;ibin++)    thetabinsp.push_back(ibin*0.05-1);
      for(int ibin=14;ibin<=26;ibin++)    thetabinsm.push_back(ibin*0.05-1);
    }

    TLatex *t=new TLatex();
    t->SetNDC();

    TGraphErrors *gduthp[nfiles];
    TGraphErrors *gduthm[nfiles];
        
    TF1 *fduthp[nfiles];
    TF1 *fduthm[nfiles];

    stringstream sduthp1[nfiles];
    stringstream sduthp2[nfiles];
    
    stringstream sduthm1[nfiles];
    stringstream sduthm2[nfiles];

    bool fit = true;

    for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++ ){

      cout << "Plotting " << types.at(ifile) << "..." << endl;
      
      can[idet]->cd(1);
      
      hduthp[ifile][idet]->SetName(Form("hduthp_%i",ifile));
      gduthp[ifile] = getTGraphFromTH2(hduthp[ifile][idet],thetabinsp, colors[ifile] );
      fduthp[ifile]=new TF1(Form("fduthp_%i",ifile),"pol1");
      fduthp[ifile]->SetLineColor( colors[ifile] );
      if(fit) gduthp[ifile]->Fit(fduthp[ifile]);
         
      float dwp     = fduthp[ifile]->GetParameter(1);
      float dwperr  = fduthp[ifile]->GetParError(1);
      float bp      = fduthp[ifile]->GetParameter(0);
      float bperr   = fduthp[ifile]->GetParError(0);
      
      sduthp1[ifile] << "#DeltaW = " << fround(dwp,3) << " #pm " << fround(dwperr,3) << " #mum" << endl;
      sduthp2[ifile] << "#Deltatan(LA) = " << fround(bp/(235.-dwp),3) << " #pm " << fround(bperr/(235.-dwp),4) << endl;
      
      gduthp[ifile]->SetTitle(Form("%s (v+ modules)",subdets.at(idet)));
      gduthp[ifile]->Draw("AP");
      gduthp[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      gduthp[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      //gduthp[ifile]->GetXaxis()->SetLimits(-1,1);
      gduthp[ifile]->SetMinimum(-10);
      gduthp[ifile]->SetMaximum(20);
      
      
      can[idet]->cd(2);

      hduthm[ifile][idet]->SetName(Form("hduthm_%i",ifile));
      gduthm[ifile] = getTGraphFromTH2(hduthm[ifile][idet],thetabinsm, colors[ifile] );
      fduthm[ifile]=new TF1(Form("fduthm_%i",ifile),"pol1");
      fduthm[ifile]->SetLineColor( colors[ifile] );
      if(fit) gduthm[ifile]->Fit(fduthm[ifile]);
         
      float dwm     = fduthm[ifile]->GetParameter(1);
      float dwmerr  = fduthm[ifile]->GetParError(1);
      float bm      = fduthm[ifile]->GetParameter(0);
      float bmerr   = fduthm[ifile]->GetParError(0);
      
      sduthm1[ifile] << "#DeltaW = " << fround(dwm,3) << " #pm " << fround(dwmerr,3) << " #mum" << endl;
      sduthm2[ifile] << "#Deltatan(LA) = " << fround(bm/(235.-dwm),3) << " #pm " << fround(bmerr/(235.-dwm),4) << endl;
      
      gduthm[ifile]->SetTitle(Form("%s (v- modules)",subdets.at(idet)));
      gduthm[ifile]->Draw("AP");
      gduthm[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      gduthm[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      //gduthm[ifile]->GetXaxis()->SetLimits(-1,1);
      gduthm[ifile]->SetMinimum(-20);
      gduthm[ifile]->SetMaximum(20);
      //leg1->Draw();
    }

    for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++ ){
      
      can[idet]->cd(1);
      if( ifile == 0 )   gduthp[ifile]->Draw("AP");
      else               gduthp[ifile]->Draw("sameP");
      //leg1->Draw();

      if(fit){
        t->SetTextColor( colors[ifile] );
        t->DrawLatex(0.15,0.85-ifile*0.2,sduthp1[ifile].str().c_str());
        t->DrawLatex(0.15,0.75-ifile*0.2,sduthp2[ifile].str().c_str());
      }

      can[idet]->cd(2);
      if( ifile == 0 )   gduthm[ifile]->Draw("AP");
      else               gduthm[ifile]->Draw("sameP");
      //leg1->Draw();

      if(fit){
        t->SetTextColor( colors[ifile] );
        t->DrawLatex(0.15,0.85-ifile*0.2,sduthm1[ifile].str().c_str());
        t->DrawLatex(0.15,0.75-ifile*0.2,sduthm2[ifile].str().c_str());
      }
    
    }

    can[idet]->cd(3);

    TGraphErrors *gduthpdiff = diffTGraph(gduthp[0],gduthp[2],Form("%s DECO - PEAK (v+ modules)",subdets.at(idet)),
                                          "<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
    TF1* fduthpdiff=new TF1("fduthpdiff","pol1",-0.5,0.5);
    gduthpdiff->Fit(fduthpdiff,"R");
    gduthpdiff->Draw("AP");
    
    float dwp     = fduthpdiff->GetParameter(1);
    float dwperr  = fduthpdiff->GetParError(1);
    float bp      = fduthpdiff->GetParameter(0);
    float bperr   = fduthpdiff->GetParError(0);
    
    stringstream sduthp1diff;
    stringstream sduthp2diff;

    sduthp1diff << "#DeltaW = " << fround(dwp,3) << " #pm " << fround(dwperr,3) << " #mum" << endl;
    sduthp2diff << "#Deltatan(LA) = " << fround(bp/(235.-dwp),3) << " #pm " << fround(bperr/(235.-dwp),4) << endl;
  
    t->SetTextColor(1);
    t->DrawLatex(0.15,0.85,sduthp1diff.str().c_str());
    t->DrawLatex(0.15,0.75,sduthp2diff.str().c_str());

    can[idet]->cd(4);

    TGraphErrors *gduthmdiff = diffTGraph(gduthm[0],gduthm[2],Form("%s DECO - PEAK (v- modules)",subdets.at(idet)),
                                          "<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
    TF1* fduthmdiff=new TF1("fduthmdiff","pol1",-0.5,0.5);
    gduthmdiff->Fit(fduthmdiff,"R");
    gduthmdiff->Draw("AP");

    float dwm     = fduthmdiff->GetParameter(1);
    float dwmerr  = fduthmdiff->GetParError(1);
    float bm      = fduthmdiff->GetParameter(0);
    float bmerr   = fduthmdiff->GetParError(0);
    
    stringstream sduthm1diff;
    stringstream sduthm2diff;

    sduthm1diff << "#DeltaW = " << fround(dwm,3) << " #pm " << fround(dwmerr,3) << " #mum" << endl;
    sduthm2diff << "#Deltatan(LA) = " << fround(bm/(235.-dwm),3) << " #pm " << fround(bmerr/(235.-dwm),4) << endl;
  
    t->DrawLatex(0.15,0.85,sduthm1diff.str().c_str());
    t->DrawLatex(0.15,0.75,sduthm2diff.str().c_str());

  }
  

  
  for(unsigned int ican = 0 ; ican < ndet ; ican++ ){
    can[ican]->Modified();
    can[ican]->Update();
    if(printgif) can[ican]->Print(Form("plots/duvsdtantheta_%s.gif",subdets.at(ican)));
  }
  
}




inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}


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

TGraphErrors *getTGraphFromTH2(TH2F* h,vector<float> xbins, int color){

  const unsigned int nbins = xbins.size()-1;
  TH2* htemp[nbins];
  TH1* hx[nbins];
  TH1* hy[nbins];
    
  float x[nbins];
  float y[nbins];
  float xerr[nbins];
  float yerr[nbins];
  
  for(unsigned int ibin=0;ibin<nbins;ibin++){
    
    htemp[ibin] = suppressHist(h,ibin,xbins.at(ibin),xbins.at(ibin+1));
    
    hx[ibin]    = htemp[ibin]->ProjectionX();
    hy[ibin]    = htemp[ibin]->ProjectionY();

    hx[ibin] -> StatOverflows(kFALSE);
    hy[ibin] -> StatOverflows(kFALSE);
    
    x[ibin]     = hx[ibin]->GetMean(1);
    xerr[ibin]  = (hx[ibin]->GetEntries()>0) ? 
      hx[ibin]->GetRMS(1)/sqrt(hx[ibin]->GetEntries()) : 0.;
    
    y[ibin]     = hy[ibin]->GetMean(1);
      
    yerr[ibin]  = (hy[ibin]->GetEntries()>0) ? 
      hy[ibin]->GetRMS(1)/sqrt(hy[ibin]->GetEntries()) : 0.;
  
  }
  
  TGraphErrors *g=new TGraphErrors(nbins,x,y,xerr,yerr);
  g->GetXaxis()->SetLimits(xbins.at(0),xbins.at(nbins));

  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(0.5);
  
  return g;
}



TGraphErrors* diffTGraph(TGraphErrors* g1, TGraphErrors *g2, char* title, char* xtitle, char* ytitle){

  Double_t* x1  = g1->GetX();
  Double_t* y1  = g1->GetY();
  Double_t* ex1 = g1->GetEX();
  Double_t* ey1 = g1->GetEY();
  Int_t n1      = g1->GetN();
  
  Double_t* x2  = g2->GetX();
  Double_t* y2  = g2->GetY();
  Double_t* ex2 = g2->GetEX();
  Double_t* ey2 = g2->GetEY();
  Int_t n2      = g2->GetN();

  assert(n1 == n2);
  
  Double_t* x  = new Double_t[n1];
  Double_t* y  = new Double_t[n1];
  Double_t* ex = new Double_t[n1];
  Double_t* ey = new Double_t[n1];

  for(int i = 0 ; i < n1 ; i++){

    x[i]  = 0.5 * ( x1[i] + x2[i] );
    ex[i] = 0.5 * sqrt( pow(ex1[i],2) + pow(ex2[i],2) );

    ey[i] = sqrt( pow(ey1[i],2) + pow(ey2[i],2) );
    y[i]  = y2[i] - y1[i];

  }



  TGraphErrors *g = new TGraphErrors(n1,x,y,ex,ey);
  g->SetTitle(title);
  g->GetYaxis()->SetTitle(ytitle);
  g->GetXaxis()->SetTitle(xtitle);
  g -> SetMarkerStyle(8);
  g -> SetMarkerSize(0.5);
  return g;
}
