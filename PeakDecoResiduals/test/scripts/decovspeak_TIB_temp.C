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

void plotHists(TH1* h1, TH1* h2,char* title, char* xtitle,float xmin, float xmax,int rebin);
TGraphErrors *getTGraphFromTH2(TH2F* h,vector<float> xbins,int color=2);
TH2* suppressHist(TH2* hist,int iclone,float xmin,float xmax);
inline double fround(double n, unsigned d);
bool draw(char* var, vector<char*> cantitles);
void plotHistsTH2(TH2* h,char* title, char* xtitle, char* ytitle, 
                  float xmin, float xmax, float ymin, float ymax,int fit,TCanvas* c, int ipad);
void plotHists(TH1F** h, const unsigned int nhist, char* title, char* xtitle, 
               float xmin, float xmax, int rebin, int fit);
void setStats(TH1F** h, const int nhist, double startingY, double startingX, double height);
TGraphErrors* diffTGraph(TGraphErrors* g1, TGraphErrors *g2, 
                         char* title = "DECO - PEAK", char* xtitle = "Layer", char* ytitle = "");

void decovspeak(bool printgif = false){
  
  //load macros and set style-----------------------------------
  //gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  //gStyle->SetOptStat("mr");
  //gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  int colors[]       = {2,4,6,1};
  
  vector<char*> filenames;
  vector<char*> types;
  
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_ALLPEAK/merged.root");
  filenames.push_back("Commissioning10-GOODCOLL-v8_ALLPEAK/res/merged.root");
  
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_copy/merged_2.root");
  filenames.push_back("Commissioning10-GOODCOLL-v8/res/merged.root");

  filenames.push_back("Commissioning10-GOODCOLL-v8_toblatebp06/res/merged.root");

  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom_mintrkmom1/merged.root");
  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom/merged.root");
  
  
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_toblatebp063/merged.root");
  
  types.push_back("peak");
  types.push_back("deco");
  types.push_back("deco (BP)");  
  //types.push_back("MC");

  vector<char*> cantitles;
  //cantitles.push_back("du_dw");         //delta u, delta w TH1s
  //cantitles.push_back("duvsdth");         //v+/- split du vs. delta tan(theta) TH2s
  cantitles.push_back("duvsdth_tgraph");  //v+/- split du vs. delta tan(theta) TGraphs
  
  const unsigned int nfiles = filenames.size();
  const unsigned int ncan   = cantitles.size();

  TCanvas *can[ncan];
  int idx = 0;

  TFile* file[nfiles];
  TH1F*  hdu[nfiles];
  TH1F*  hdw[nfiles];
  TH2F*  hduthp[nfiles];
  TH2F*  hduthm[nfiles];

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
   
    hdu[ifile] = (TH1F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/du_TIB");
    hdw[ifile] = (TH1F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/dw_TIB");

    //    if( strcmp(types.at(ifile),"deco (BP)") == 0){
      hduthp[ifile] = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vp_wp_all");
      TH2F* htp     = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vp_wm_all");
      hduthp[ifile] -> Add(htp);
      
      hduthm[ifile] = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vm_wp_all");
      TH2F* htm     = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vm_wm_all");
      hduthm[ifile] -> Add(htm);
    // }else{
//       hduthp[ifile] = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vp_wp");
//       TH2F* htp     = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vp_wm");
//       hduthp[ifile] -> Add(htp);
      
//       hduthm[ifile] = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vm_wp");
//       TH2F* htm     = (TH2F*) file[ifile]->Get("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_TIB_vm_wm");
//       hduthm[ifile] -> Add(htm);
//     }
     
  }
  
  
  cout << "Plotting du, dw histos... " << endl;

  if(draw("du_dw",cantitles)){
    can[idx]=new TCanvas(Form("can_%i",idx),"du_dw",1200,450);
    can[idx]->Divide(2,1);
    
    can[idx]->cd(1);
    setStats(hdu,nfiles,0.8,0.65,0.15);  
    plotHists(hdu,nfiles,"TIB","#Delta u [#mum]",-1000,1000,10,3);
    leg1->Draw();

    can[idx]->cd(2);
    setStats(hdw,nfiles,0.8,0.65,0.15);  
    plotHists(hdw,nfiles,"TIB","#Delta w [#mum]",-2000,2000,10,3);
    leg1->Draw();
        
    idx++;
  }

  cout << "Plotting du vs. delta tan theta split v+/v- histos... " << endl;
  //du vs. delta tan(theta) split v+/v- histos-----------------------------------
  
  if(draw("duvsdth",cantitles)){

    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450*nfiles);
    can[idx]->Divide(nfiles,2);
    
    for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
      
      plotHistsTH2(hduthp[ifile],
                   Form("TIB %s (v+)",types.at(ifile)),"tan(#theta_{trk})-tan(#theta_{L})",
                   "#Delta u [#mum]",-0.3,0.5,-200,200,0,can[idx],ifile+1);
      
      plotHistsTH2(hduthm[ifile],
                   Form("TIB %s (v-)",types.at(ifile)),"tan(#theta_{trk})-tan(#theta_{L})",
                        "#Delta u [#mum]",-0.5,0.3,-200,200,0,can[idx],ifile+2);
      
    }
  
    idx++;

  }
  
  cout << "Plotting du vs. delta tan theta split v+/v- TGraphs... " << endl;

  //TGraph du vs. delta tan theta----------------------------------------

  if(draw("duvsdth_tgraph",cantitles)){

    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
  
    vector<float> thetabinsp;
    vector<float> thetabinsm;

    for(int ibin=14;ibin<=30;ibin++)    thetabinsp.push_back(ibin*0.05-1);
    for(int ibin=10;ibin<=26;ibin++)    thetabinsm.push_back(ibin*0.05-1);
    
 
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
      
      can[idx]->cd(1);
      
      hduthp[ifile]->SetName(Form("hduthp_%i",ifile));
      gduthp[ifile] = getTGraphFromTH2(hduthp[ifile],thetabinsp, colors[ifile] );
      fduthp[ifile]=new TF1(Form("fduthp_%i",ifile),"pol1");
      fduthp[ifile]->SetLineColor( colors[ifile] );
      if(fit) gduthp[ifile]->Fit(fduthp[ifile]);
         
      float dwp     = fduthp[ifile]->GetParameter(1);
      float dwperr  = fduthp[ifile]->GetParError(1);
      float bp      = fduthp[ifile]->GetParameter(0);
      float bperr   = fduthp[ifile]->GetParError(0);
      
      sduthp1[ifile] << "#DeltaW = " << fround(dwp,3) << " #pm " << fround(dwperr,3) << " #mum" << endl;
      sduthp2[ifile] << "#Deltatan(LA) = " << fround(bp/(235.-dwp),3) << " #pm " << fround(bperr/(235.-dwp),4) << endl;
      
      gduthp[ifile]->SetTitle("TIB (V+)");
      gduthp[ifile]->Draw("AP");
      gduthp[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      gduthp[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      //gduthp[ifile]->GetXaxis()->SetLimits(-1,1);
      gduthp[ifile]->SetMinimum(-25);
      gduthp[ifile]->SetMaximum(20);
      
      
      can[idx]->cd(2);

      hduthm[ifile]->SetName(Form("hduthm_%i",ifile));
      gduthm[ifile] = getTGraphFromTH2(hduthm[ifile],thetabinsm, colors[ifile] );
      fduthm[ifile]=new TF1(Form("fduthm_%i",ifile),"pol1");
      fduthm[ifile]->SetLineColor( colors[ifile] );
      if(fit) gduthm[ifile]->Fit(fduthm[ifile]);
         
      float dwm     = fduthm[ifile]->GetParameter(1);
      float dwmerr  = fduthm[ifile]->GetParError(1);
      float bm      = fduthm[ifile]->GetParameter(0);
      float bmerr   = fduthm[ifile]->GetParError(0);
      
      sduthm1[ifile] << "#DeltaW = " << fround(dwm,3) << " #pm " << fround(dwmerr,3) << " #mum" << endl;
      sduthm2[ifile] << "#Deltatan(LA) = " << fround(bm/(235.-dwm),3) << " #pm " << fround(bmerr/(235.-dwm),4) << endl;
      
      gduthm[ifile]->SetTitle("TIB (V-)");
      gduthm[ifile]->Draw("AP");
      gduthm[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      gduthm[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      //gduthm[ifile]->GetXaxis()->SetLimits(-1,1);
      gduthm[ifile]->SetMinimum(-35);
      gduthm[ifile]->SetMaximum(15);
      //leg1->Draw();
    }

    for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++ ){
      
      can[idx]->cd(1);
      if( ifile == 0 )   gduthp[ifile]->Draw("AP");
      else               gduthp[ifile]->Draw("sameP");
      leg1->Draw();

      if(fit){
        t->SetTextColor( colors[ifile] );
        t->DrawLatex(0.15,0.85-ifile*0.2,sduthp1[ifile].str().c_str());
        t->DrawLatex(0.15,0.75-ifile*0.2,sduthp2[ifile].str().c_str());
      }

      can[idx]->cd(2);
      if( ifile == 0 )   gduthm[ifile]->Draw("AP");
      else               gduthm[ifile]->Draw("sameP");
      //leg1->Draw();

      if(fit){
        t->SetTextColor( colors[ifile] );
        t->DrawLatex(0.15,0.85-ifile*0.2,sduthm1[ifile].str().c_str());
        t->DrawLatex(0.15,0.75-ifile*0.2,sduthm2[ifile].str().c_str());
      }
    
    }

    TCanvas *cdiff = new TCanvas("cdiff","duvsdtanthetadiff",1200,450);
    cdiff->Divide(2,1);

    cdiff->cd(1);


    TGraphErrors *gduthpdiff = diffTGraph(gduthp[0],gduthp[2],"TIB DECO - PEAK (v+)","<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
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

    cdiff->cd(2);

    TGraphErrors *gduthmdiff = diffTGraph(gduthm[0],gduthm[2],"TIB DECO - PEAK (v-)","<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
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

    idx++;
  }

  
  for(unsigned int ican=0;ican<ncan;ican++){
    can[ican]->Modified();
    can[ican]->Update();
    if(printgif) can[ican]->Print(Form("plots/%s_TIB.gif",cantitles.at(ican)));
  }
  
}



void plotHists(TH1* h1, TH1* h2, 
	       char* title, char* xtitle,
	       float xmin, float xmax,int rebin){
   
  if(h1==0 && h2==0){
    cout << " ERROR TWO NULL HISTOS " << endl;

  }

  TLatex *l1=new TLatex();
  TLatex *l2=new TLatex();
  stringstream s1;
  stringstream s2;


  if(h1!=0){

    if(rebin>1)   h1->Rebin(rebin);
    
    if(h1->Integral()>0) h1->Scale(1./h1->Integral());

    TF1 *f1=new TF1("f1","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",xmin,xmax);
    f1->SetParameters(h1->GetMaximum()/3.,0,200,h1->GetMaximum()/3.,400,h1->GetMaximum()/3.,600);
    f1->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
    
    f1->SetLineWidth(1);
    f1->SetLineColor(2);
    h1->Fit(f1);
    h1->SetTitle(title);
    h1->GetXaxis()->SetTitle(xtitle);
    h1->SetLineColor(2);
    h1->Draw();
   
    s1<<"#bar{x} = "<<fround(f1->GetParameter(1),2)<<" #mum"<<endl;
    l1->SetTextSize(0.06);
    l1->SetNDC();
    l1->SetTextColor(2);
    l1->DrawLatex(0.15,0.8,s1.str().c_str());
  }


  if(h2!=0){

    if(rebin>1)   h2->Rebin(rebin);
    
    if(h2->Integral()>0) h2->Scale(1./h2->Integral());

    TF1 *f2=new TF1("f2","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",xmin,xmax);
    f2->SetParameters(h2->GetMaximum()/3.,0,200,h2->GetMaximum()/3.,400,h2->GetMaximum()/3.,600);
    f2->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
    
    f2->SetLineWidth(1);
    f2->SetLineColor(4);
    h2->Fit(f2);
    h2->SetTitle(title);
    h2->GetXaxis()->SetTitle(xtitle);
    h2->SetLineColor(4);
    h2->Draw();
    
    s2<<"#bar{x} = "<<fround(f2->GetParameter(1),2)<<" #mum"<<endl;
    l2->SetTextSize(0.06);
    l2->SetNDC();
    l2->SetTextColor(4);
    l2->DrawLatex(0.15,0.8,s2.str().c_str());
  }

  if(h1!=0 && h2!=0){
    h1->Draw();
    h2->Draw("same");

    l1->DrawLatex(0.15,0.8,s1.str().c_str());
    l2->DrawLatex(0.15,0.7,s2.str().c_str());

    float max=h1->GetMaximum();
    if(h2->GetMaximum()>max)max=h2->GetMaximum();
    h1->SetMaximum(1.05*max);
  }

}

inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

bool draw(char* var, vector<char*> cantitles){
  
  for(unsigned int i=0;i<cantitles.size();i++) {
    if(strcmp(var,cantitles.at(i)) == 0) return true;
  }
  return false;
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

void plotHistsTH2(TH2* h,char* title, char* xtitle, char* ytitle, 
                  float xmin, float xmax, float ymin, float ymax,int fit,TCanvas* c, int ipad){

  c->cd(ipad);
  TPad *plotpad=new TPad("plotpad","",0.1,0.1,1.,0.9);
  plotpad->Draw();
  plotpad->cd();

  if(fit>0){
    TF1 *f=new TF1("f","[0]+[1]*x",xmin,xmax);
    f->SetParNames("y-int","slope");
    f->SetLineWidth(2);
    f->SetLineColor(2);
    h->Fit(f);
  }
  
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->SetTitle("");
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("");
  h->Draw("colz");



  c->cd(ipad);
  TPad *xpad=new TPad("xpad","",0.,0.,1.,0.1);
  xpad->Draw();
  xpad->cd();
  TLatex *x=new TLatex();
  x->SetNDC();
  x->SetTextSize(0.8);
  x->DrawLatex(0.5,0.3,xtitle);

  c->cd(ipad);
  TPad *ypad=new TPad("ypad","",0.,0.,0.1,1.);
  ypad->Draw();
  ypad->cd();
  TLatex *y=new TLatex();
  y->SetNDC();
  y->SetTextSize(0.8);
  y->SetTextAngle(90);
  y->DrawLatex(0.6,0.5,ytitle);

  c->cd(ipad);
  TPad *tpad=new TPad("tpad","",0.,0.9,1.,1.);
  tpad->Draw();
  tpad->cd();
  TLatex *t=new TLatex();
  t->SetNDC();
  t->SetTextSize(0.6);
  t->DrawLatex(0.2,0.1,title);

  TLatex *l=new TLatex();
  l->SetNDC();
  l->SetTextSize(0.6);
  l->SetTextColor(2);
  float cor=h->GetCorrelationFactor(1,2);
  stringstream s;
  s<<"C="<<(int)(1000*cor)<<"#times10^{-3}"<<endl;
  l->DrawLatex(0.7,0.1,s.str().c_str());
}


void setStats(TH1F** h, const int nhist, double startingY, double startingX, double height){

  int colors[]={2,4,6,1}; 
  TPaveStats* st[nhist];

  if (startingY<0){
    
    for(int i=0;i<nhist;i++)  h[i]->SetStats(0);
    
  } else {
    
    gStyle->SetOptStat("mr");

    for(int i=0 ; i < nhist ; i++){
      
      if(i == 0) h[i]->Draw();
      else       h[i]->Draw("sames");

      gPad->Update(); 
      st[i] = (TPaveStats*) h[i]->GetListOfFunctions()->FindObject("stats");
      st[i]->SetX1NDC(startingX);
      st[i]->SetX2NDC(startingX+0.30);
      st[i]->SetY1NDC(startingY-(i+1)*height+0.05);
      st[i]->SetY2NDC(startingY-i*height+0.05);
      st[i]->SetTextColor(colors[i]);
    }
  }
}



void plotHists(TH1F** h, const unsigned int nhist, char* title, char* xtitle, 
    float xmin, float xmax, int rebin, int fit){
   
  int colors[]={2,4,6,1};
  TF1* f[nhist];
  
  for(unsigned int i = 0 ; i < nhist ; ++i){
  
    if(rebin>1)   h[i]->Rebin(rebin);

    if(h[i] -> Integral() > 0) h[i] -> Scale( 1. / h[i]->Integral() );
    
    if(fit>0){
      if(fit==1){
        f[i] = new TF1("f","[0]*exp(-0.5*pow((x-[1])/[2],2))",xmin,xmax);
        f[i] -> SetParameters( h[i]->GetMaximum() , 0 , h[i]->GetRMS(1));
      }
      if(fit==2){
        f[i] = new TF1("f","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))",xmin,xmax);
        f[i]->SetParameters(h[i]->GetMaximum()/2.,0,200,h[i]->GetMaximum()/2.,500);
        f[i]->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}");
      }
      if(fit==3){
        f[i] = new TF1("f","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",xmin,xmax);
        f[i]->SetParameters(h[i]->GetMaximum()/3.,0,200,h[i]->GetMaximum()/3.,400,h[i]->GetMaximum()/3.,600);
        f[i]->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
      }
      
      f[i]->SetLineWidth(1);
      f[i]->SetLineColor(colors[i]);
      h[i]->Fit(f[i]);
      
    }

    h[i]->SetLineColor( colors[i] );
    h[i]->SetTitle( title );
    h[i]->GetXaxis()->SetTitle( xtitle );
  }

  for(unsigned int i = 0 ; i < nhist ; ++i){
    
    if(i == 0)   h[i]->Draw( "" );
    else         h[i]->Draw("same");
    
    if(fit>0&&fit<4){
      stringstream s;
      s<<"#bar{x} = "<<fround(f[i]->GetParameter(1),1)<<" #mum"<<endl;
      TLatex *l=new TLatex();
      l->SetTextSize(0.06);
      l->SetNDC();
      l->SetTextColor( colors[i] );
      l->DrawLatex(0.15,0.45-i*0.1,s.str().c_str());
    }
  }
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
