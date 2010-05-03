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
                         //char* title = "DECO - PEAK", char* xtitle = "Layer", char* ytitle = "");
                         string title = "DECO - PEAK", string xtitle = "Layer", string ytitle = "");

void corrections(bool printgif = false){
  
  //load macros and set style-----------------------------------
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  int colors[]       = {2,4,6,1};
  
  //files to process
  vector<string> filenames;
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_ALLPEAK/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8_ALLPEAK/res/merged.root");
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_copy/merged_2.root");
  //filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8_newcfg_ALLPEAK/res/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8_newcfg_ALLPEAK/res/merged.root");
  filenames.push_back("crabjobs/Commissioning10-Apr20Skim_GOODCOLL-v1_ALLPEAK/res/merged.root");
  filenames.push_back("crabjobs/Commissioning10-Apr20Skim_GOODCOLL-v1/res/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-Apr20Skim_GOODCOLL-v1_lateBP/res/merged.root");  
  //filenames.push_back("crabjobs/Commissioning10-Apr20Skim_GOODCOLL-v1_tanla/res/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-Apr20Skim_GOODCOLL-v1_tanla_lateBP/res/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-Apr20Skim_GOODCOLL-v1_copy/res/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8_Geometry/res/merged.root");
  //filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8_toblatebp06/merged.root");
  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom_mintrkmom1/merged.root");
  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom/merged.root");
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_toblatebp063/merged.root");
 
  //string corrtitle = "DECO (#Deltatan(#theta_{LA})+BP-corrected) - PEAK";
  //string corrtitle = "DECO (#Deltatan(#theta_{LA})-corrected) - PEAK";
  //string corrtitle = "DECO (BP-corrected) - PEAK";
  string corrtitle = "DECO - PEAK";
  int ncorr = 1;


  //labels for files
  vector<string> types;
  types.push_back("peak");
  types.push_back("deco");
  //types.push_back("deco (BP)");
  //types.push_back("deco (#Deltatan(#theta_{LA}))");  
  //types.push_back("deco (#Deltatan(#theta_{LA}) + BP)");
  
  assert ( filenames.size() == types.size() );

  //list of subdetectors
  vector<string> subdets;
  subdets.push_back("TOB");
  subdets.push_back("TIB");
  //subdets.push_back("TECthick");
  //subdets.push_back("TEC");

  const unsigned int nfiles = filenames.size();
  const unsigned int ndet = subdets.size();

  TFile*   file[nfiles];
  TH2F*    hduthp[nfiles][ndet];
  TH2F*    hduthm[nfiles][ndet];
  TCanvas* can[ndet];
  
  float dwp_[ndet];
  float dwm_[ndet];
  float dtanlap_[ndet];
  float dtanlam_[ndet];

  //dummy legend
  TLegend *leg1=new TLegend(0.15,0.65,0.35,0.85);

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
    TH1F *hdummy=new TH1F(Form("hdummy_%i",ifile),Form("hdummy_%i",ifile),1,0,1);
    hdummy->SetLineColor(colors[ifile]);
    hdummy->SetMarkerColor(colors[ifile]);
    leg1->AddEntry(hdummy,types.at(ifile).c_str());
  }

  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);

  cout << "Getting histos... " << endl;

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
    cout << "Opening " << types.at(ifile) <<" file   " << filenames.at(ifile) << endl;
    file[ifile] = TFile::Open(filenames.at(ifile).c_str());
    
    for(unsigned int idet = 0 ; idet < ndet ; idet++){
    
      hduthp[ifile][idet] = (TH2F*) file[ifile]->
        Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vp_wp_all",subdets.at(idet).c_str()));
      
      TH2F* htp     = (TH2F*) file[ifile]->
        Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vp_wm_all",subdets.at(idet).c_str()));
      
      hduthp[ifile][idet] -> Add(htp);
        
      hduthm[ifile][idet] = (TH2F*) file[ifile]->
        Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vm_wp_all",subdets.at(idet).c_str()));
        
      TH2F* htm     = (TH2F*) file[ifile]->
        Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_%s_vm_wm_all",subdets.at(idet).c_str()));
      
      hduthm[ifile][idet] -> Add(htm);
      
    }
  }
  
  vector<float> thetabinsp;
  vector<float> thetabinsm;
  float thickness = 0;

  for(unsigned int idet = 0 ; idet < ndet ; idet++){

    cout << "Plotting " << subdets.at(idet) << "..." << endl;

    can[idet]=new TCanvas(Form("%s_can",subdets.at(idet).c_str()),Form("%s canvas",subdets.at(idet).c_str()),1200,900);
    can[idet]->Divide(2,2);
  
    thetabinsp.clear();
    thetabinsm.clear();

    if( strcmp( subdets.at(idet).c_str() , "TOB" ) == 0 ){
      for(int ibin=14;ibin<=30;ibin++)    thetabinsp.push_back(ibin*0.05-1);
      for(int ibin=10;ibin<=26;ibin++)    thetabinsm.push_back(ibin*0.05-1);
      thickness = 470;
    }   
    else if( strcmp( subdets.at(idet).c_str() , "TIB" ) == 0 ){
      for(int ibin=14;ibin<=26;ibin++)    thetabinsp.push_back(ibin*0.05-1);
      for(int ibin=14;ibin<=26;ibin++)    thetabinsm.push_back(ibin*0.05-1);
      thickness = 290;
    }
    else if( strcmp( subdets.at(idet).c_str() , "TECthick" ) == 0 ){
      for(int ibin=70;ibin<=130;ibin++)    thetabinsp.push_back(ibin*0.01-1);
      for(int ibin=70;ibin<=130;ibin++)    thetabinsm.push_back(ibin*0.01-1);
      thickness = 470;
    }
    else{
      cout << "ERROR SET TAN(THETA) RANGE" <<endl;
      exit(0);
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

    bool fit = false;

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
      
      gduthp[ifile]->SetTitle(Form("%s (v+ modules)",subdets.at(idet).c_str()));
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
      
      gduthm[ifile]->SetTitle(Form("%s (v- modules)",subdets.at(idet).c_str()));
      gduthm[ifile]->Draw("AP");
      gduthm[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      gduthm[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      //gduthm[ifile]->GetXaxis()->SetLimits(-1,1);
      gduthm[ifile]->SetMinimum(-20);
      gduthm[ifile]->SetMaximum(10);
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

    if( nfiles > 1){
      
      can[idet]->cd(3);
      
      TGraphErrors *gduthpdiff = diffTGraph(gduthp[0],gduthp[ncorr],Form("%s %s (v+ modules)",subdets.at(idet).c_str(),corrtitle.c_str()),
                                            "<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
      TF1* fduthpdiff=new TF1("fduthpdiff","pol1",-0.5,0.5);
      gduthpdiff->Fit(fduthpdiff,"R");
      gduthpdiff->Draw("AP");
      gduthpdiff->SetMinimum(-10);
      gduthpdiff->SetMaximum(10);

      float dwp     = fduthpdiff->GetParameter(1);
      float dwperr  = fduthpdiff->GetParError(1);
      float bp      = fduthpdiff->GetParameter(0);
      float bperr   = fduthpdiff->GetParError(0);
      
      stringstream sduthp1diff;
      stringstream sduthp2diff;

      sduthp1diff << "#DeltaW = " << fround(dwp,3) << " #pm " << fround(dwperr,3) << " #mum" << endl;
      sduthp2diff << "#Deltatan(LA) = " << fround(bp/(thickness/2-dwp),3) << " #pm " << fround(bperr/(thickness/2-dwp),4) << endl;
    
      t->SetTextColor(1);
      t->DrawLatex(0.15,0.85,sduthp1diff.str().c_str());
      t->DrawLatex(0.15,0.75,sduthp2diff.str().c_str());

      dwp_[idet]     = dwp;
      dtanlap_[idet] = bp / ( thickness / 2. - dwp ); 

      can[idet]->cd(4);
      
      TGraphErrors *gduthmdiff = diffTGraph(gduthm[0],gduthm[ncorr],Form("%s %s (v- modules)",subdets.at(idet).c_str(),corrtitle.c_str()),
                                            "<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
      TF1* fduthmdiff=new TF1("fduthmdiff","pol1",-0.5,0.5);
      gduthmdiff->Fit(fduthmdiff,"R");
      gduthmdiff->Draw("AP");
      gduthmdiff->SetMinimum(-10);
      gduthmdiff->SetMaximum(10);

      
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

      dwm_[idet]     = dwm;
      dtanlam_[idet] = bm / ( thickness / 2. - dwm ); 

      
    }
  }

  
  for(unsigned int ican = 0 ; ican < ndet ; ican++ ){
    can[ican]->Modified();
    can[ican]->Update();
    if(printgif) can[ican]->Print(Form("plots/duvsdtantheta_%s.gif",subdets.at(ican).c_str()));
  }

  for( unsigned int idet = 0 ; idet < ndet ; idet++ ){
    
    if( strcmp( subdets.at(idet).c_str() , "TOB" ) == 0 )               thickness = 470;
    else if( strcmp( subdets.at(idet).c_str() , "TIB" ) == 0 )          thickness = 290;
    else if( strcmp( subdets.at(idet).c_str() , "TECthick" ) == 0 )     thickness = 470;
    else{
      cout << "ERROR SET TAN(THETA) RANGE" <<endl;
      exit(0);
    }


    cout << endl << "-------------------------------------------------" << endl << endl;
    
    cout << "Corrections for " << subdets.at(idet) << " thickness " << thickness << endl;
    cout << endl;
    cout << "dW (v+) " << dwp_[idet] << " um, dW (v-) " << dwm_[idet] << " um" << endl;
    cout << "dW (avg) " << 0.5 * ( dwp_[idet] + dwm_[idet] ) << endl;
    cout << "frac " << ( dwp_[idet] + dwm_[idet] ) / thickness << endl;
    cout << endl;
    cout << "dtan(LA) (v+) " << dtanlap_[idet] << " , dtan(LA) (v-) " << dtanlam_[idet] << endl;
    cout << "dtan(LA) (avg) " << 0.5 * ( dtanlap_[idet] - dtanlam_[idet] ) << endl;
    
    cout << endl << "-------------------------------------------------" << endl << endl;
  }
}


//15.919 19.246

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



//TGraphErrors* diffTGraph(TGraphErrors* g1, TGraphErrors *g2, char* title, char* xtitle, char* ytitle){
TGraphErrors* diffTGraph(TGraphErrors* g1, TGraphErrors *g2, string title, string xtitle, string ytitle){

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
  g->SetTitle(title.c_str());
  g->GetYaxis()->SetTitle(ytitle.c_str());
  g->GetXaxis()->SetTitle(xtitle.c_str());
  g -> SetMarkerStyle(8);
  g -> SetMarkerSize(0.5);
  return g;
}
