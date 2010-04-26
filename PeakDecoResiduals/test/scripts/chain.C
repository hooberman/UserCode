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
#endif

TH2* suppressHist(TH2* hist,int iclone,float xmin,float xmax);
void setStats(TH1F** h, const int nhist, double startingY, double startingX = .1, double height = 0.15);
void plotHists(TH1F** h, const unsigned int nhist, string title, string xtitle, float xmin, float xmax, int rebin = 1, int fit = 0);
TGraphErrors *getTGraphFromTH2(TH2F* h,vector<float> xbins, int method = 0 , int invert = 0, int color = 1);
TGraphErrors* diffTGraph(TGraphErrors* g1, TGraphErrors *g2, string title = "DECO - PEAK", string xtitle = "Layer", string ytitle = "");
inline double fround(double n, unsigned d){ return floor(n * pow(10., d) + .5) / pow(10., d); }
bool draw(string var, vector<string> cantitles);
int getDTSlice(float t);
string getStringFromTF1( TF1* f);

void chain(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
   
  enum subdetenum      { PixelBarrel = 0, PixelEndcap = 1, TIB = 2, TID = 3,
			 TOB = 4, TECthick = 5, TECthin = 6};

  string detnames[]   = {"PixBa","PixEC","TIB","TID","TOB","TECthick","TECthin"};

  int detlayers[]    = {0,0,4,0,6,0,0};
  int colors[]       = {2,4,6,1};

  //User input----------------------------------------------------------------------
  
  bool writeTFile   = false;
  bool printgif     = true;
  int  mysubdet     = TOB; 
  int  ndiv         = 100;
  const int nlayers = detlayers[mysubdet];
  
  vector<string> cantitles;
  //cantitles.push_back("du_dw");
  //cantitles.push_back("nstripsvstantrk");
  //cantitles.push_back("nstripsvstantrktgraph");
  //cantitles.push_back("nstripsvstantrktgraphdt");
  //cantitles.push_back("tanladt");
  cantitles.push_back("duvsdtantheta_tgraph");
  //cantitles.push_back("duvsdtantheta_layers_tgraph");
  
  vector<string> filenames; 
  vector<string> filetypes;

  filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8_ALLPEAK/res/merged.root");
  filetypes.push_back("PEAK");
  
  filenames.push_back("crabjobs/Commissioning10-GOODCOLL-v8/res/merged.root");
  filetypes.push_back("DECO");

  //   filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom/merged.root");
  //   filetypes.push_back("MC");

  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_ALLPEAK/merged.root");
  //filetypes.push_back("PEAK");
  
  //filenames.push_back("crabjobs/lpc/Commissioning10-GOODCOLL-v8_mintrkmom1_copy/merged_2.root");
  //filetypes.push_back("DECO");
  
  //filenames.push_back("crabjobs/lpc/Spring10-START3X_V26A_356ReReco-v1_standard_geom_mintrkmom1/merged.root");
  //filetypes.push_back("MC");

  //filenames.push_back("crabjobs/trial3/root/peak.root");
  //filetypes.push_back("PEAK");
  
  //filenames.push_back("crabjobs/trial3/root/deco.root");
  //filetypes.push_back("DECO");
  
  const unsigned int nfiles = filenames.size();
  
  //----------------------------------------------------------------------------------

  TCanvas *can[20];
  int idx = 0;

  TFile *ofile;
  if(writeTFile){
    ofile=TFile::Open("histos.root","RECREATE");
    ofile->cd();
  }

  //dummy legend
  TLegend *leg1=new TLegend(0.15,0.65,0.35,0.85);

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
    TH1F *hdummy=new TH1F(Form("hdummy_%i",ifile),Form("hdummy_%i",ifile),1,0,1);
    hdummy->SetLineColor(colors[ifile]);
    hdummy->SetMarkerColor(colors[ifile]);
    leg1->AddEntry(hdummy,filetypes.at(ifile).c_str());
  }

  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);
  
  //Declare/initialize histos
  TH1F* hdu[nfiles];
  TH1F* hdw[nfiles];
  
  TH2F* hduthetap[nfiles];
  TH2F* hduthetam[nfiles];
  
  TProfile* hduthetap_tprof[nfiles];
  TProfile* hduthetam_tprof[nfiles];
  
  TH2F* duthetap_layer[nfiles][nlayers];
  TH2F* duthetam_layer[nfiles][nlayers];
  
  TProfile* duthetap_layer_tprof[nfiles][nlayers];
  TProfile* duthetam_layer_tprof[nfiles][nlayers];
  
  TH1F* tanla_layer[nfiles][nlayers];
  TH2F* nstripstantrk[nfiles];
  TH2F* nstripstantrk_dt[nfiles][8];
  TH1F* dttime_dt[nfiles][8];
  
  cout << "Booking histos..." << endl;

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
  
    nstripstantrk[ifile]  = new TH2F(Form("nstripstantrk%s",filetypes.at(ifile).c_str()),
                                     Form("nstripstantrk%s",filetypes.at(ifile).c_str()),400,-2,2,21,-0.5,20.5);

    hdu[ifile]       = new TH1F(Form("du_%s",filetypes.at(ifile).c_str()),
                                Form("du_%s",filetypes.at(ifile).c_str()),1000,-1000,1000);
    
    hdw[ifile]       = new TH1F(Form("dw_%s",filetypes.at(ifile).c_str()),
                                Form("dw_%s",filetypes.at(ifile).c_str()),1000,-2000,2000);
    
    hduthetap[ifile] = new TH2F(Form("duthetap_%s",filetypes.at(ifile).c_str()),
                                Form("duthetap_%s",filetypes.at(ifile).c_str()),100,-2,2,200,-500,500);
    

    hduthetam[ifile] = new TH2F(Form("duthetam_%s",filetypes.at(ifile).c_str()),
                                Form("duthetam_%s",filetypes.at(ifile).c_str()),100,-2,2,200,-500,500);

    hduthetap_tprof[ifile] = new TProfile(Form("duthetap_tprof_%s",filetypes.at(ifile).c_str()),
                                          Form("duthetap_tprof_%s",filetypes.at(ifile).c_str()),100,-2,2,-500,500);
    
    hduthetam_tprof[ifile] = new TProfile(Form("duthetam_tprof_%s",filetypes.at(ifile).c_str()),
                                          Form("duthetam_tprof_%s",filetypes.at(ifile).c_str()),100,-2,2,-500,500);

    for(int ih = 0 ; ih < nlayers ; ih++){

      duthetap_layer[ifile][ih] = new TH2F(Form("dutheta%sp_%i",filetypes.at(ifile).c_str(),ih),
                                           Form("dutheta%sp_%i",filetypes.at(ifile).c_str(),ih),100,-2,2,200,-500,500);
      
      duthetam_layer[ifile][ih] = new TH2F(Form("dutheta%sm_%i",filetypes.at(ifile).c_str(),ih),
                                           Form("dutheta%sm_%i",filetypes.at(ifile).c_str(),ih),100,-2,2,200,-500,500);
      
//       duthetap_layer_tprof[ifile][ih] = new TProfile(Form("dutheta%sp_tprof_%i",filetypes.at(ifile).c_str(),ih),
//                                                      Form("dutheta%sp_tprof_%i",filetypes.at(ifile).c_str(),ih),100,-2,2,-500,500);
      
//       duthetam_layer_tprof[ifile][ih] = new TProfile(Form("dutheta%sm_tprof_%i",filetypes.at(ifile).c_str(),ih),
//                                                      Form("dutheta%sm_tprof_%i",filetypes.at(ifile).c_str(),ih),100,-2,2,-500,500);

      duthetap_layer_tprof[ifile][ih] = new TProfile(Form("dutheta%sp_tprof_%i",filetypes.at(ifile).c_str(),ih),
                                                     Form("dutheta%sp_tprof_%i",filetypes.at(ifile).c_str(),ih),20,-0.5,0.5,-500,500);
      
      duthetam_layer_tprof[ifile][ih] = new TProfile(Form("dutheta%sm_tprof_%i",filetypes.at(ifile).c_str(),ih),
                                                     Form("dutheta%sm_tprof_%i",filetypes.at(ifile).c_str(),ih),20,-0.5,0.5,-500,500);
      
      tanla_layer[ifile][ih]    = new TH1F(Form("tanla%s_%i",filetypes.at(ifile).c_str(),ih), 
                                           Form("tanla%s_%i",filetypes.at(ifile).c_str(),ih), 10000,0.0,0.2);
      
    }

    for(int ih = 0 ; ih < 8 ; ih++){
      nstripstantrk_dt[ifile][ih]    = new TH2F(Form("nstripstantrk%s_dt_%i",filetypes.at(ifile).c_str(),ih),
                                                Form("nstripstantrk%s_dt_%i",filetypes.at(ifile).c_str(),ih),400,-2,2,21,-0.5,20.5);

      dttime_dt[ifile][ih]           = new TH1F(Form("dttime%s_dt_%i",filetypes.at(ifile).c_str(),ih),
                                                Form("dttime%s_dt_%i",filetypes.at(ifile).c_str(),ih),100,-50,50);

    }
    
  }

  TChain* chain[nfiles];

  Float_t du;
  Float_t dw;
  Float_t dtanth;
  Float_t tantrk;
  Int_t   subdet;
  Int_t   v;
  Int_t   w;
  Int_t   nstrips;
  Float_t dttime;
  Int_t   layer;
  Float_t tanla;
  
  cout << "Looping over files..." << endl;

  for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++){
    
    chain[ifile] = new TChain("PeakDecoResiduals/t","Tree");
    chain[ifile] -> Add( filenames.at(ifile).c_str() );
    
    int nentriestot = chain[ifile]->GetEntries();
    int nentries = 0;

    TObjArray *listOfFiles = chain[ifile]->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TChainElement* currentFile = 0;


    while((currentFile = (TChainElement*)fileIter.Next())) {
      
      TFile file(currentFile->GetTitle());
      
      TTree *tree = (TTree*)file.Get("PeakDecoResiduals/t");
      tree->SetBranchAddress("du",       &du);
      tree->SetBranchAddress("dw",       &dw);
      tree->SetBranchAddress("dtanth",   &dtanth);
      tree->SetBranchAddress("subdet",   &subdet);
      tree->SetBranchAddress("v",        &v);
      tree->SetBranchAddress("w",        &w);
      tree->SetBranchAddress("nstrips",  &nstrips);
      tree->SetBranchAddress("tantrk" ,  &tantrk);
      tree->SetBranchAddress("dttime" ,  &dttime);
      tree->SetBranchAddress("layer"  ,  &layer);
      tree->SetBranchAddress("tanla"  ,  &tanla);
      
      for(unsigned int ientry = 0; ientry < tree->GetEntries()/ndiv  ; ++ientry) {
        
        if(nentries % 100000 == 0) cout<<filetypes[ifile]<<" event "<<nentries<<" / "<<nentriestot<<endl;
        nentries++;
        
        tree->GetEntry(ientry);
        
        if(subdet == mysubdet ){
          
          hdu[ifile]->Fill(dtanth > 0 ? du : -du);
          hdw[ifile]->Fill(dw);
          
          if(getDTSlice(dttime)>-1)
            dttime_dt[ifile][getDTSlice(dttime)]->Fill(dttime);
            
          if(mysubdet == TIB || mysubdet == TOB)
            tanla_layer[ifile][layer-1]->Fill(v > 0 ? -tanla : tanla);
            
          if(v>0){
           
            hduthetap[ifile]->Fill(dtanth,du);
            
            if(mysubdet == TIB || mysubdet == TOB){
              duthetap_layer[ifile][layer-1]->Fill(dtanth,du);
              duthetap_layer_tprof[ifile][layer-1]->Fill(dtanth,du);
            }
            
            nstripstantrk[ifile]->Fill(tantrk,nstrips);
            
            if(getDTSlice(dttime)>-1)
              nstripstantrk_dt[ifile][getDTSlice(dttime)]->Fill(tantrk,nstrips);
            
          }
          else{
          
            hduthetam[ifile]->Fill(dtanth,du);
            
            if(mysubdet == TIB || mysubdet == TOB){
              duthetam_layer[ifile][layer-1]->Fill(dtanth,du);
              duthetam_layer_tprof[ifile][layer-1]->Fill(dtanth,du);
            }

            nstripstantrk[ifile]->Fill(-tantrk,nstrips);
            
            if(getDTSlice(dttime)>-1)
              nstripstantrk_dt[ifile][getDTSlice(dttime)]->Fill(-tantrk,nstrips);
          }          
        }
      }
    }
  }


  

  //delta u, delta w TH1s
  if(draw("du_dw",cantitles)){
    can[idx]=new TCanvas(Form("can_%i",idx),"du_dw",1200,450);
    can[idx]->Divide(2,1);
    
    can[idx]->cd(1);
    setStats(hdu,nfiles,0.8,0.65,0.15);  
    plotHists(hdu,nfiles,"TOB","#Delta u [#mum]",-1000,1000,10,3);
    leg1->Draw();
    can[idx]->cd(2);
    setStats(hdw,nfiles,0.8,0.65,0.15);  
    plotHists(hdw,nfiles,"TOB","#Delta w [#mum]",-2000,2000,10,3);
    leg1->Draw();
        
    idx++;
  }


  if(draw("duvsdtantheta_tgraph",cantitles)){

    can[idx]=new TCanvas("duvsdtantheta_tgraph_can","duvsdtantheta_tgraph_can",1200,900);
    can[idx]->Divide(2,2);

      vector<float> thetabinsp;
      vector<float> thetabinsm;
      
      thetabinsp.clear();
      thetabinsm.clear();
      
      for(int ibin=14;ibin<=30;ibin++)    thetabinsp.push_back(ibin*0.05-1);
      for(int ibin=10;ibin<=26;ibin++)    thetabinsm.push_back(ibin*0.05-1);
            
      TLatex *t=new TLatex();
      t->SetNDC();
      
      TGraphErrors *gduthetap[nfiles];
      TGraphErrors *gduthetam[nfiles];
      
      TF1 *fduthetap[nfiles];
      TF1 *fduthetam[nfiles];
      
      stringstream sduthetap1[nfiles];
      stringstream sduthetap2[nfiles];
      
      stringstream sduthetam1[nfiles];
      stringstream sduthetam2[nfiles];
      
      bool fit = false;
      
      for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++ ){
        
        cout << "Plotting " << filetypes.at(ifile) << "..." << endl;
        
        can[idx]->cd(1);
        
        hduthetap[ifile]->SetName(Form("hduthetap_%i",ifile));
        gduthetap[ifile] = getTGraphFromTH2(hduthetap[ifile],thetabinsp, colors[ifile] );
        fduthetap[ifile]=new TF1(Form("fduthetap_%i",ifile),"pol1");
        fduthetap[ifile]->SetLineColor( colors[ifile] );
        if(fit) gduthetap[ifile]->Fit(fduthetap[ifile]);
         
        float dwp     = fduthetap[ifile]->GetParameter(1);
        float dwperr  = fduthetap[ifile]->GetParError(1);
        float bp      = fduthetap[ifile]->GetParameter(0);
        float bperr   = fduthetap[ifile]->GetParError(0);
      
        sduthetap1[ifile] << "#DeltaW = " << fround(dwp,3) << " #pm " << fround(dwperr,3) << " #mum" << endl;
        sduthetap2[ifile] << "#Deltatan(LA) = " << fround(bp/(235.-dwp),3) << " #pm " << fround(bperr/(235.-dwp),4) << endl;
      
        gduthetap[ifile]->SetTitle(Form("%s (v+ modules)",detnames[mysubdet].c_str()));
        gduthetap[ifile]->Draw("AP");
        gduthetap[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
        gduthetap[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
        //gduthetap[ifile]->GetXaxis()->SetLimits(-1,1);
        //gduthetap[ifile]->SetMinimum(-10);
        //gduthetap[ifile]->SetMaximum(20);
      
      
        can[idx]->cd(2);

        hduthetam[ifile]->SetName(Form("hduthetam_%i",ifile));
        gduthetam[ifile] = getTGraphFromTH2(hduthetam[ifile],thetabinsm, colors[ifile] );
        fduthetam[ifile]=new TF1(Form("fduthetam_%i",ifile),"pol1");
        fduthetam[ifile]->SetLineColor( colors[ifile] );
        if(fit) gduthetam[ifile]->Fit(fduthetam[ifile]);
         
        float dwm     = fduthetam[ifile]->GetParameter(1);
        float dwmerr  = fduthetam[ifile]->GetParError(1);
        float bm      = fduthetam[ifile]->GetParameter(0);
        float bmerr   = fduthetam[ifile]->GetParError(0);
      
        sduthetam1[ifile] << "#DeltaW = " << fround(dwm,3) << " #pm " << fround(dwmerr,3) << " #mum" << endl;
        sduthetam2[ifile] << "#Deltatan(LA) = " << fround(bm/(235.-dwm),3) << " #pm " << fround(bmerr/(235.-dwm),4) << endl;
      
        gduthetam[ifile]->SetTitle(Form("%s (v- modules)",detnames[mysubdet].c_str()));
        gduthetam[ifile]->Draw("AP");
        gduthetam[ifile]->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
        gduthetam[ifile]->GetYaxis()->SetTitle("<#Deltau> [#mum]");
        //gduthetam[ifile]->GetXaxis()->SetLimits(-1,1);
        //gduthetam[ifile]->SetMinimum(-20);
        //gduthetam[ifile]->SetMaximum(20);
        //leg1->Draw();
      }

      for(unsigned int ifile = 0 ; ifile < nfiles ; ifile++ ){
      
        can[idx]->cd(1);
        if( ifile == 0 )   gduthetap[ifile]->Draw("AP");
        else               gduthetap[ifile]->Draw("sameP");
        //leg1->Draw();

        if(fit){
          t->SetTextColor( colors[ifile] );
          t->DrawLatex(0.15,0.85-ifile*0.2,sduthetap1[ifile].str().c_str());
          t->DrawLatex(0.15,0.75-ifile*0.2,sduthetap2[ifile].str().c_str());
        }

        can[idx]->cd(2);
        if( ifile == 0 )   gduthetam[ifile]->Draw("AP");
        else               gduthetam[ifile]->Draw("sameP");
        //leg1->Draw();

        if(fit){
          t->SetTextColor( colors[ifile] );
          t->DrawLatex(0.15,0.85-ifile*0.2,sduthetam1[ifile].str().c_str());
          t->DrawLatex(0.15,0.75-ifile*0.2,sduthetam2[ifile].str().c_str());
        }
    
      }

      if( nfiles > 1){
      
        can[idx]->cd(3);
      
        TGraphErrors *gduthetapdiff = diffTGraph(gduthetap[0],gduthetap[1],Form("%s DECO - PEAK (v+ modules)",detnames[mysubdet].c_str()),
                                              "<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
        TF1* fduthetapdiff=new TF1("fduthetapdiff","pol1",-0.5,0.5);
        gduthetapdiff->Fit(fduthetapdiff,"R");
        gduthetapdiff->Draw("AP");
      
        float dwp     = fduthetapdiff->GetParameter(1);
        float dwperr  = fduthetapdiff->GetParError(1);
        float bp      = fduthetapdiff->GetParameter(0);
        float bperr   = fduthetapdiff->GetParError(0);
      
        stringstream sduthetap1diff;
        stringstream sduthetap2diff;
      
        sduthetap1diff << "#DeltaW = " << fround(dwp,3) << " #pm " << fround(dwperr,3) << " #mum" << endl;
        sduthetap2diff << "#Deltatan(LA) = " << fround(bp/(235.-dwp),3) << " #pm " << fround(bperr/(235.-dwp),4) << endl;
      
        t->SetTextColor(1);
        t->DrawLatex(0.15,0.85,sduthetap1diff.str().c_str());
        t->DrawLatex(0.15,0.75,sduthetap2diff.str().c_str());
      
        can[idx]->cd(4);
      
        TGraphErrors *gduthetamdiff = diffTGraph(gduthetam[0],gduthetam[1],Form("%s DECO - PEAK (v- modules)",detnames[mysubdet].c_str()),
                                                 "<tan(#theta_{trk})-tan(#theta_{LA})>","<#Deltau> (#mum)");
        TF1* fduthetamdiff=new TF1("fduthetamdiff","pol1",-0.5,0.5);
        gduthetamdiff->Fit(fduthetamdiff,"R");
        gduthetamdiff->Draw("AP");
      
        float dwm     = fduthetamdiff->GetParameter(1);
        float dwmerr  = fduthetamdiff->GetParError(1);
        float bm      = fduthetamdiff->GetParameter(0);
        float bmerr   = fduthetamdiff->GetParError(0);
      
        stringstream sduthetam1diff;
        stringstream sduthetam2diff;
      
        sduthetam1diff << "#DeltaW = " << fround(dwm,3) << " #pm " << fround(dwmerr,3) << " #mum" << endl;
        sduthetam2diff << "#Deltatan(LA) = " << fround(bm/(235.-dwm),3) << " #pm " << fround(bmerr/(235.-dwm),4) << endl;
      
        t->DrawLatex(0.15,0.85,sduthetam1diff.str().c_str());
        t->DrawLatex(0.15,0.75,sduthetam2diff.str().c_str());
      
      }
    }


              
  if(draw("duvsdtantheta_layers_tgraph",cantitles)){
             
    //declare stuff
    vector<float> thetabinsp;
    thetabinsp.clear();
    thetabinsp.push_back(-0.30);
    thetabinsp.push_back(-0.25);
    thetabinsp.push_back(-0.20);
    thetabinsp.push_back(-0.15);
    thetabinsp.push_back(-0.10);
    thetabinsp.push_back(-0.05);
    thetabinsp.push_back(0.00);
    thetabinsp.push_back(0.05);
    thetabinsp.push_back(0.10);
    thetabinsp.push_back(0.15);
    thetabinsp.push_back(0.20);
    thetabinsp.push_back(0.25);
    thetabinsp.push_back(0.30);
    thetabinsp.push_back(0.35);
    thetabinsp.push_back(0.40);
    thetabinsp.push_back(0.45);
    thetabinsp.push_back(0.50);
    
    vector<float> thetabinsm;
    thetabinsm.clear();
    thetabinsm.push_back(-0.50);
    thetabinsm.push_back(-0.45);
    thetabinsm.push_back(-0.40);
    thetabinsm.push_back(-0.35);
    thetabinsm.push_back(-0.30);
    thetabinsm.push_back(-0.25);
    thetabinsm.push_back(-0.20);
    thetabinsm.push_back(-0.15);
    thetabinsm.push_back(-0.10);
    thetabinsm.push_back(-0.05);
    thetabinsm.push_back(0.00);
    thetabinsm.push_back(0.05);
    thetabinsm.push_back(0.10);
    thetabinsm.push_back(0.15);
    thetabinsm.push_back(0.20);
    thetabinsm.push_back(0.25);
    thetabinsm.push_back(0.30);
    
    float thetapmin = thetabinsp.at(0);
    float thetapmax = thetabinsp.at(thetabinsp.size()-1);
    float thetammin = thetabinsm.at(0);
    float thetammax = thetabinsm.at(thetabinsm.size()-1);

    //vector<float> thetabins;
    //thetabins.clear();
    //thetabins.push_back(-1.);
    //thetabins.push_back(-0.75);
    //thetabins.push_back(-0.5);
    //thetabins.push_back(-0.25);
    //thetabins.push_back(0.);
    //thetabins.push_back(0.25);
    //thetabins.push_back(0.5);
    //thetabins.push_back(0.75);
    //thetabins.push_back(1.);
    
    TGraphErrors *gduthetap_layer[nfiles][nlayers];
    TGraphErrors *gduthetam_layer[nfiles][nlayers];
    
    TF1 *fduthetap_layer[nfiles][nlayers];
    TF1 *fduthetam_layer[nfiles][nlayers];
    
    string sduthetap_layer[nfiles][nlayers];
    string sduthetam_layer[nfiles][nlayers];

    TF1 *fduthetap_layer_tprof[nfiles][nlayers];
    TF1 *fduthetam_layer_tprof[nfiles][nlayers];
    
    string sduthetap_layer_tprof[nfiles][nlayers];
    string sduthetam_layer_tprof[nfiles][nlayers];

    TLatex *t=new TLatex();
    t->SetNDC();

    TCanvas *tanlacan[nlayers];
   
    float tanla[nfiles][nlayers];
    float tanlaerr[nfiles][nlayers];
   
    float deltaw[nfiles][nlayers];
    float deltawerr[nfiles][nlayers];
    float deltawp[nfiles][nlayers];
    float deltawperr[nfiles][nlayers];
    float deltawm[nfiles][nlayers];
    float deltawmerr[nfiles][nlayers];
    
    float bp;
    float bm;
    float mp;
    float mm;
    float berrp;
    float berrm;
    float merrp;
    float merrm;
    float b;
    float m;
    float berr; 
    float merr; 

    //get du vs. dtantheta for each layer
    bool drawlayers = true;
    
    for( int ih = 0 ; ih < nlayers ; ih++ ){

      if(drawlayers){
        tanlacan[ih]=new TCanvas(Form("tanlacan_%i",ih),Form("tanlacan_%i",ih),1200,900);
        tanlacan[ih]->Divide(2,2);

      }
      
      for( unsigned int ifile = 0 ; ifile < nfiles ; ifile++){

        tanla[ifile][ih]    = tanla_layer[ifile][ih]->GetMean(1);
        tanlaerr[ifile][ih] = pow(tanla_layer[ifile][ih]->GetRMS(1),2);
        
        //get du vs. dtantheta TGraph from TH2
        tanlacan[ih]->cd(1);
        duthetap_layer[ifile][ih]  -> SetName(Form("duthetap_layer_%i_%i",ifile,ih));
        gduthetap_layer[ifile][ih] =  getTGraphFromTH2(duthetap_layer[ifile][ih],thetabinsp,0,0,colors[ifile]);
        //gduthetap_layer[ifile][ih] -> SetMarkerColor(colors[ifile]);
        //gduthetap_layer[ifile][ih] -> SetLineColor(colors[ifile]);
        fduthetap_layer[ifile][ih] =  new TF1(Form("fduthetap_layer_%i_%i",ifile,ih),"pol1");
        fduthetap_layer[ifile][ih] -> SetLineColor(colors[ifile]);
        gduthetap_layer[ifile][ih] -> Fit(fduthetap_layer[ifile][ih]);
        sduthetap_layer[ifile][ih] = getStringFromTF1( fduthetap_layer[ifile][ih] );

        tanlacan[ih]->cd(2);
        duthetam_layer[ifile][ih]  -> SetName(Form("duthetam_layer_%i_%i",ifile,ih));
        gduthetam_layer[ifile][ih] =  getTGraphFromTH2(duthetam_layer[ifile][ih],thetabinsm,0,0,colors[ifile]);
        //gduthetam_layer[ifile][ih] -> SetMarkerColor(colors[ifile]);
        //gduthetam_layer[ifile][ih] -> SetLineColor(colors[ifile]);
        fduthetam_layer[ifile][ih] =  new TF1(Form("fduthetam_layer_%i_%i",ifile,ih),"pol1");
        fduthetam_layer[ifile][ih] -> SetLineColor(colors[ifile]);
        gduthetam_layer[ifile][ih] -> Fit(fduthetam_layer[ifile][ih]);
        sduthetam_layer[ifile][ih] = getStringFromTF1( fduthetam_layer[ifile][ih] );
        
        //du vs. dtantheta TProfile
        tanlacan[ih]->cd(3);
        duthetap_layer_tprof[ifile][ih]  -> SetMarkerColor(colors[ifile]);
        duthetap_layer_tprof[ifile][ih]  -> SetLineColor(colors[ifile]);
        duthetap_layer_tprof[ifile][ih]  -> SetMarkerStyle(8);
        duthetap_layer_tprof[ifile][ih]  -> SetMarkerSize( 0.5 );
        fduthetap_layer_tprof[ifile][ih] =  new TF1(Form("fduthetap_layer_tprof_%i_%i",ifile,ih),"pol1",thetapmin,thetapmax);
        fduthetap_layer_tprof[ifile][ih] -> SetLineColor(colors[ifile]);
        duthetap_layer_tprof[ifile][ih]  -> Fit(fduthetap_layer_tprof[ifile][ih],"R");
        duthetap_layer_tprof[ifile][ih]  -> GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
        duthetap_layer_tprof[ifile][ih]  -> GetYaxis()->SetTitle("<#Deltau> [#mum]");
        duthetap_layer_tprof[ifile][ih]  -> SetTitle(Form("%s (v+)",detnames[mysubdet].c_str()));
        duthetap_layer_tprof[ifile][ih]  -> GetXaxis()->SetRangeUser(thetapmin,thetapmax);
        sduthetap_layer_tprof[ifile][ih] = getStringFromTF1( fduthetap_layer_tprof[ifile][ih] );


        tanlacan[ih]->cd(4);
        duthetam_layer_tprof[ifile][ih]  -> SetMarkerColor(colors[ifile]);
        duthetam_layer_tprof[ifile][ih]  -> SetLineColor(colors[ifile]);
        duthetam_layer_tprof[ifile][ih]  -> SetMarkerStyle(8);
        duthetam_layer_tprof[ifile][ih]  -> SetMarkerSize( 0.5 );
        fduthetam_layer_tprof[ifile][ih] =  new TF1(Form("fduthetam_layer_tprof_%i_%i",ifile,ih),"pol1",thetammin,thetammax);
        fduthetam_layer_tprof[ifile][ih] -> SetLineColor(colors[ifile]);
        duthetam_layer_tprof[ifile][ih]  -> Fit(fduthetam_layer_tprof[ifile][ih],"R");
        duthetam_layer_tprof[ifile][ih]  -> GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
        duthetam_layer_tprof[ifile][ih]  -> GetYaxis()->SetTitle("<#Deltau> [#mum]");
        duthetam_layer_tprof[ifile][ih]  -> SetTitle(Form("%s (v-)",detnames[mysubdet].c_str()));
        duthetam_layer_tprof[ifile][ih]  -> GetXaxis()->SetRangeUser(thetammin,thetammax);
        sduthetam_layer_tprof[ifile][ih] = getStringFromTF1( fduthetam_layer_tprof[ifile][ih] );

      }

      if(drawlayers){
        
        tanlacan[ih]->cd(1);
          
        TMultiGraph *mgdtp=new TMultiGraph();
        for(unsigned int ifile = 0; ifile < nfiles ; ifile++){
          mgdtp->Add(gduthetap_layer[ifile][ih]);
        }
        mgdtp->SetTitle(Form("%s (v+)",detnames[mysubdet].c_str()));
        mgdtp->Draw("AP");
        mgdtp->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
        mgdtp->GetYaxis()->SetTitle("<#Deltau> [#mum]");
        mgdtp->SetTitle(Form("%s (v+)",detnames[mysubdet].c_str()));
        //mgdtp->GetXaxis()->SetLimits(-1.,1.);
        mgdtp->GetYaxis()->SetLimits(-20,20);
        
        for(unsigned int ifile = 0; ifile < nfiles ; ifile++){
          t->SetTextColor(colors[ifile]);
          t->DrawLatex(0.25,0.85-ifile*0.1,sduthetap_layer[ifile][ih].c_str());
        }    
          
        tanlacan[ih]->cd(2);
          
        TMultiGraph *mgdtm=new TMultiGraph();
        for(unsigned int ifile = 0; ifile < nfiles ; ifile++){
          mgdtm->Add(gduthetam_layer[ifile][ih]);
        }
        mgdtm->SetTitle(Form("%s (v-)",detnames[mysubdet].c_str()));
        mgdtm->Draw("AP");
        mgdtm->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
        mgdtm->GetYaxis()->SetTitle("<#Deltau> [#mum]");
        mgdtm->SetTitle(Form("%s (v-)",detnames[mysubdet].c_str()));
        //mgdtm->GetXaxis()->SetLimits(-1.,1.);
        mgdtm->GetYaxis()->SetLimits(-20,20);
          
        for(unsigned int ifile = 0; ifile < nfiles ; ifile++){
          t->SetTextColor(colors[ifile]);
          t->DrawLatex(0.25,0.85-ifile*0.1,sduthetam_layer[ifile][ih].c_str());
        }  
        
        tanlacan[ih]->cd(3);
        for(unsigned int ifile = 0; ifile < nfiles ; ifile++){
          if(ifile==0)     duthetap_layer_tprof[ifile][ih]->Draw();
          else             duthetap_layer_tprof[ifile][ih]->Draw("same");

          t->SetTextColor(colors[ifile]);
          t->DrawLatex(0.25,0.85-ifile*0.1,sduthetap_layer_tprof[ifile][ih].c_str());
        }
        
        tanlacan[ih]->cd(4);
        for(unsigned int ifile = 0; ifile < nfiles ; ifile++){
          if(ifile==0)     duthetam_layer_tprof[ifile][ih]->Draw();
          else             duthetam_layer_tprof[ifile][ih]->Draw("same");

          t->SetTextColor(colors[ifile]);
          t->DrawLatex(0.25,0.85-ifile*0.1,sduthetam_layer_tprof[ifile][ih].c_str());
        }
      }
      
      //calculate LA and delta W for each layer
      for(unsigned int ifile = 0; ifile < nfiles ; ifile++){

       
        bp    = fduthetap_layer[ifile][ih]->GetParameter(0);
        bm    = fduthetam_layer[ifile][ih]->GetParameter(0);
        
        mp    = fduthetap_layer[ifile][ih]->GetParameter(1);
        mm    = fduthetam_layer[ifile][ih]->GetParameter(1);
        
        berrp = fduthetap_layer[ifile][ih]->GetParError(0);
        berrm = fduthetam_layer[ifile][ih]->GetParError(0);
        
        merrp = fduthetap_layer[ifile][ih]->GetParError(1);
        merrm = fduthetam_layer[ifile][ih]->GetParError(1);
        
        b     = 0.5*(bp - bm);
        berr  = 0.5*(pow(berrp,2) + pow(berrm,2));
        
        m     = 0.5*(mp + mm);
        merr  = 0.5*sqrt(pow(merrp,2) + pow(merrm,2));
        
        
        tanla[ifile][ih]     -= b/(235. - m);
        tanlaerr[ifile][ih]  += pow(berr/(250.-m),2);
        tanlaerr[ifile][ih]   = sqrt(tanlaerr[ifile][ih]);
        
        deltaw[ifile][ih]     = m;
        deltawerr[ifile][ih]  = merr;

        deltawp[ifile][ih]    = mp;
        deltawperr[ifile][ih] = merrp;

        deltawm[ifile][ih]    = mm;
        deltawmerr[ifile][ih] = merrm;
 
      }
    }    

    //dtan(LA), dw for all layers and deco-peak offsets
    can[idx]=new TCanvas(Form("can_%i",idx),"dw_alllayers",1200,900);
    can[idx]->Divide(2,2);

    float layers[nlayers];
    float layerserr[nlayers];
  
    for(int ih = 0 ; ih < nlayers ; ih++){
      layers[ih]    = ih+1;
      layerserr[ih] = 0.;
    }
    
    TGraphErrors *gtanla[nfiles];
    TGraphErrors *gdw[nfiles];
    TGraphErrors *gdwp[nfiles];
    TGraphErrors *gdwm[nfiles];

    float tanlatemp[nlayers];
    float tanlaerrtemp[nlayers];
    float dwtemp[nlayers];
    float dwerrtemp[nlayers];
    float dwptemp[nlayers];
    float dwperrtemp[nlayers];
    float dwmtemp[nlayers];
    float dwmerrtemp[nlayers];
    
    for(unsigned int ifile = 0; ifile < nfiles ; ifile++){

      can[idx]->cd(1);

      for(int ilayer = 0 ; ilayer < nlayers ; ilayer++){
        tanlatemp[ilayer]     = tanla[ifile][ilayer];
        tanlaerrtemp[ilayer]  = tanlaerr[ifile][ilayer];
        dwtemp[ilayer]        = deltaw[ifile][ilayer];
        dwerrtemp[ilayer]     = deltawerr[ifile][ilayer];
      
      }
    
      gtanla[ifile]=new TGraphErrors(nlayers,layers,tanlatemp,layerserr,tanlaerrtemp);
      gtanla[ifile]->SetTitle(detnames[mysubdet].c_str());
      gtanla[ifile]->GetYaxis()->SetTitle("tan(#theta_{LA})");
      gtanla[ifile]->GetXaxis()->SetTitle("Layer");
      gtanla[ifile] -> SetMarkerColor(colors[ifile]);
      gtanla[ifile] -> SetLineColor(colors[ifile]);
      gtanla[ifile] -> SetMarkerStyle(8);
      gtanla[ifile] -> SetMarkerSize(0.5);
      gtanla[ifile] -> SetMinimum(0.06);
      gtanla[ifile] -> SetMaximum(0.1);
      
      if(ifile==0) gtanla[ifile]->Draw("AP");
      else         gtanla[ifile]->Draw("sameP");

      
      can[idx]->cd(2);

      gdw[ifile]=new TGraphErrors(nlayers,layers,dwtemp,layerserr,dwerrtemp);
      gdw[ifile]->SetTitle(detnames[mysubdet].c_str());
      gdw[ifile]->GetYaxis()->SetTitle("#DeltaW (#mum)");
      gdw[ifile]->GetXaxis()->SetTitle("Layer");
      gdw[ifile] -> SetMarkerColor(colors[ifile]);
      gdw[ifile] -> SetLineColor(colors[ifile]);
      gdw[ifile] -> SetMarkerStyle(8);
      gdw[ifile] -> SetMarkerSize(0.5);
      gdw[ifile] -> SetMinimum(-10);
      gdw[ifile] -> SetMaximum(40);
      
      if(ifile==0) gdw[ifile]->Draw("AP");
      else         gdw[ifile]->Draw("sameP");
    }


    if( filenames.size()>1 ){
      
      can[idx]->cd(3);
      TGraphErrors *gtanladiff = diffTGraph(gtanla[0],gtanla[1],Form("%s DECO - PEAK",detnames[mysubdet].c_str()),"Layer","tan(#theta_{LA})");
      gtanladiff->Draw("AP");
      
      can[idx]->cd(4);
      TGraphErrors *gdwdiff = diffTGraph(gdw[0],gdw[1],Form("%s DECO - PEAK",detnames[mysubdet].c_str()),"Layer","DeltaW (#mum)");
      gdwdiff->Draw("AP");

    }
    
    idx++;

    //dw for all layers split v+ vs. v-
    can[idx]=new TCanvas(Form("can_%i",idx),"dw_alllayers_split",1200,900);
    can[idx]->Divide(2,2);

    for(unsigned int ifile = 0; ifile < nfiles ; ifile++){

      for(int ilayer = 0 ; ilayer < nlayers ; ilayer++){ 
        dwptemp[ilayer]       = deltawp[ifile][ilayer];
        dwperrtemp[ilayer]    = deltawperr[ifile][ilayer];
        dwmtemp[ilayer]       = deltawm[ifile][ilayer];
        dwmerrtemp[ilayer]    = deltawmerr[ifile][ilayer];
      }

      can[idx]->cd(1);
    
      gdwp[ifile]=new TGraphErrors(nlayers,layers,dwptemp,layerserr,dwperrtemp);
      gdwp[ifile]->SetTitle(Form("%s (v+)",detnames[mysubdet].c_str()));
      gdwp[ifile]->GetYaxis()->SetTitle("#DeltaW (#mum)");
      gdwp[ifile]->GetXaxis()->SetTitle("Layer");
      gdwp[ifile] -> SetMarkerColor(colors[ifile]);
      gdwp[ifile] -> SetLineColor(colors[ifile]);
      gdwp[ifile] -> SetMarkerStyle(8);
      gdwp[ifile] -> SetMarkerSize(0.5);
      gdwp[ifile] -> SetMinimum(-50);
      gdwp[ifile] -> SetMaximum(30);
      
      if(ifile==0) gdwp[ifile]->Draw("AP");
      else         gdwp[ifile]->Draw("sameP");

      can[idx]->cd(2);
    
      gdwm[ifile]=new TGraphErrors(nlayers,layers,dwmtemp,layerserr,dwmerrtemp);
      gdwm[ifile]->SetTitle(Form("%s (v-)",detnames[mysubdet].c_str()));
      gdwm[ifile]->GetYaxis()->SetTitle("#DeltaW (#mum)");
      gdwm[ifile]->GetXaxis()->SetTitle("Layer");
      gdwm[ifile] -> SetMarkerColor(colors[ifile]);
      gdwm[ifile] -> SetLineColor(colors[ifile]);
      gdwm[ifile] -> SetMarkerStyle(8);
      gdwm[ifile] -> SetMarkerSize(0.5);
      gdwm[ifile] -> SetMinimum(-15);
      gdwm[ifile] -> SetMaximum(80);
      
      if(ifile==0) gdwm[ifile]->Draw("AP");
      else         gdwm[ifile]->Draw("sameP");

    }
    
    if( filenames.size()>1 ){
      
      can[idx]->cd(3);
      TGraphErrors *gdwpdiff = diffTGraph(gdwp[0],gdwp[1],Form("%s DECO - PEAK (v+)",detnames[mysubdet].c_str()),"Layer","DeltaW (#mum)");
      gdwpdiff->Draw("AP");
      
      can[idx]->cd(4);
      TGraphErrors *gdwmdiff = diffTGraph(gdwm[0],gdwm[1],Form("%s DECO - PEAK (v+)",detnames[mysubdet].c_str()),"Layer","DeltaW (#mum)");
      gdwmdiff->Draw("AP");
      
    }
    
    idx++;
    
  }












  /*
  if(draw("nstripsvstantrk",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
    
    plotHistsTH2(nstripstantrkpeak,
                 "TOB PEAK","tan(#theta_{trk})",
                 "nstrips",-2,2,0,20,0,can[idx],1);
    plotHistsTH2(nstripstantrkdeco,
                 "TOB DECO","tan(#theta_{trk})",
                 "nstrips",-2,2,0,20,0,can[idx],2);
    
    idx++;
  }

  if(draw("nstripsvstantrktgraph",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),800,600);
    can[idx]->cd();

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)",-1,1);
    //fnstripstantrkpeak->SetParameters(0.,1.,1.,1.,1.);
    //fnstripstantrkdeco->SetParameters(0.,1.,1.,1.,1.);

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+(x<[0])*((x-[0])*[2])+(x>0)*((x-[0])*[3])",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+(x<[0])*((x-[0])*[2])+(x>0)*((x-[0])*[3])",-1,1);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.,2.);

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+[2]*abs((x-[0]))",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+[2]*abs((x-[0]))",-1,1);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.,2.);

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+[2]*abs((x-[0]))",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+[2]*abs((x-[0]))",-1,1);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.,2.);
    //char* func="V-fit";

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","([1]+[2]*pow(x-[0],2))",-0.5,0.3);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","([1]+[2]*pow(x-[0],2))",-0.5,0.3);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.);

    TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak",
                                    "[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)+[5]*pow(x-[0],8)+[6]*pow(x-[0],10)+[7]*pow(x-[0],12)+[8]*pow(x-[0],14)",-1,1);
    TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco",
                                    "[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)+[5]*pow(x-[0],8)+[6]*pow(x-[0],10)+[7]*pow(x-[0],12)+[8]*pow(x-[0],14)",-1,1);
    fnstripstantrkpeak->SetParameters(0.,1.,1.,1.,1.,1.,1.,1.,1.);
    fnstripstantrkdeco->SetParameters(0.,1.,1.,1.,1.,1.,1.,1.,1.);
    char* func = "f(x) = #Sigma_{i=0}^{7} [ c_{i} (x-tan(#theta_{trk}^{MIN}))^{2i} ]";

    vector<float> bins;
    for(int ibin=0;ibin<201;ibin++) bins.push_back(ibin*0.01-1.);
    
    TGraphErrors *gnstripstantrkpeak = getTGraphFromTH2(nstripstantrkpeak,bins,0);
    fnstripstantrkpeak->SetLineColor(2);
    fnstripstantrkpeak->SetLineWidth(1);
    gnstripstantrkpeak->Fit(fnstripstantrkpeak,"R");
    gnstripstantrkpeak->SetMarkerColor(2);
    gnstripstantrkpeak->SetLineColor(2);
    gnstripstantrkpeak->SetMarkerStyle(8);
    gnstripstantrkpeak->SetMarkerSize(0.1);
    gnstripstantrkpeak->Draw("AP");
    gnstripstantrkpeak->SetTitle("TOB");
    gnstripstantrkpeak->GetXaxis()->SetTitle("<tan(#theta_{trk})>");
    gnstripstantrkpeak->GetYaxis()->SetTitle("<nstrips>");

    stringstream snstripstantrkpeak;
    snstripstantrkpeak<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkpeak->GetParameter(0),5)
                      <<" #pm "<<fround(fnstripstantrkpeak->GetParError(0),5)<<endl;

    TGraphErrors *gnstripstantrkdeco = getTGraphFromTH2(nstripstantrkdeco,bins,0);
    fnstripstantrkdeco->SetLineColor(4);
    fnstripstantrkdeco->SetLineWidth(1);
    gnstripstantrkdeco->Fit(fnstripstantrkdeco,"R");
    gnstripstantrkdeco->SetMarkerColor(4);
    gnstripstantrkdeco->SetLineColor(4);
    gnstripstantrkdeco->SetMarkerStyle(8);
    gnstripstantrkdeco->SetMarkerSize(0.1);
    gnstripstantrkdeco->Draw("sameP");

    stringstream snstripstantrkdeco;
    snstripstantrkdeco<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkdeco->GetParameter(0),5)
                      <<" #pm "<<fround(fnstripstantrkdeco->GetParError(0),5)<<endl;


    TLatex *t=new TLatex();
    t->SetNDC();
    t->SetTextSize(0.04);
    t->SetTextColor(2);
    t->DrawLatex(0.3,0.82,snstripstantrkpeak.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.3,0.75,snstripstantrkdeco.str().c_str());
    t->SetTextColor(1);
    t->DrawLatex(0.3,0.68,func);

    idx++;
  }

  float tanlapeak[8];
  float tanlaerrpeak[8];
  float tanladeco[8];
  float tanlaerrdeco[8];
  float dtpeak[8];
  float dterrpeak[8];
  float dtdeco[8];
  float dterrdeco[8];

  if(draw("nstripsvstantrktgraphdt",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,600);
    can[idx]->Divide(4,2);
    //can[idx]->cd();

    //TCanvas *c1=new TCanvas("c1","",800,600);
    //TCanvas *c2=new TCanvas("c2","",800,600);
    //c1->Divide(2,2);
    //c2->Divide(2,2);

    TF1* fnstripstantrkpeak_dt[8];
    TF1* fnstripstantrkdeco_dt[8];
    TGraphErrors *gnstripstantrkpeak_dt[8];
    TGraphErrors *gnstripstantrkdeco_dt[8];
    stringstream snstripstantrkpeak_dt[8];
    stringstream snstripstantrkdeco_dt[8];
    char* func = "V-fit";
    char* titles[]={
      "TOB (-20 ns < DT Time < -12 ns)",
      "TOB (-12 ns < DT Time < -8 ns)",
      "TOB (-8 ns < DT Time < -4 ns)",
      "TOB (-4 ns < DT Time < 0 ns)",
      "TOB (0 ns < DT Time < 4 ns)",
      "TOB (4 ns < DT Time < 8 ns)",
      "TOB (8 ns < DT Time < 12 ns)",
      "TOB (12 ns < DT Time < 20 ns)"
    };

    vector<float> bins2;
    bins2.clear();
    //for(int ibin2=0;ibin2<101;ibin2++) bins2.push_back(ibin2*0.02-1.);
    for(int ibin2=0;ibin2<201;ibin2++) bins2.push_back(ibin2*0.01-1.);

    for(int ih = 0 ; ih < 8 ; ih++ ){

      cout<<"DT Slice "<<ih<<endl;
      can[idx]->cd(ih+1);
      //if(ih<4) c1->cd(ih+1);
      //if(ih>3) c2->cd(ih-3);

      fnstripstantrkpeak_dt[ih] = new TF1(Form("fnstripstantrkpeak_dt_%i",ih),"[1]+[2]*abs((x-[0]))",-1,1);
      fnstripstantrkdeco_dt[ih] = new TF1(Form("fnstripstantrkdeco_dt_%i",ih),"[1]+[2]*abs((x-[0]))",-1,1);
      
      fnstripstantrkpeak_dt[ih]->SetParameters(-0.1,2.,2.,2.);
      fnstripstantrkdeco_dt[ih]->SetParameters(-0.1,2.,2.,2.);
      
      gnstripstantrkpeak_dt[ih] = getTGraphFromTH2(nstripstantrkpeak_dt[ih],bins2,0);
      fnstripstantrkpeak_dt[ih]->SetLineColor(2);
      fnstripstantrkpeak_dt[ih]->SetLineWidth(1);
      gnstripstantrkpeak_dt[ih]->Fit(fnstripstantrkpeak_dt[ih],"R");
      gnstripstantrkpeak_dt[ih]->SetMarkerColor(2);
      gnstripstantrkpeak_dt[ih]->SetLineColor(2);
      gnstripstantrkpeak_dt[ih]->SetMarkerStyle(8);
      gnstripstantrkpeak_dt[ih]->SetMarkerSize(0.1);
      gnstripstantrkpeak_dt[ih]->Draw("AP");
      gnstripstantrkpeak_dt[ih]->SetTitle(titles[ih]);
      gnstripstantrkpeak_dt[ih]->GetXaxis()->SetTitle("<tan(#theta_{trk})>");
      gnstripstantrkpeak_dt[ih]->GetYaxis()->SetTitle("<nstrips>");
      gnstripstantrkpeak_dt[ih]->SetMinimum(1);
      gnstripstantrkpeak_dt[ih]->SetMaximum(5);

      tanlapeak[ih]    = fnstripstantrkpeak_dt[ih]->GetParameter(0);
      tanlaerrpeak[ih] = fnstripstantrkpeak_dt[ih]->GetParError(0);
      dtpeak[ih]       = dttimepeak_dt[ih]->GetMean(1);
      dterrpeak[ih]    = dttimepeak_dt[ih]->GetRMS(1)/sqrt(dttimepeak_dt[ih]->GetEntries());

      snstripstantrkpeak_dt[ih]<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkpeak_dt[ih]->GetParameter(0),5)
                               <<" #pm "<<fround(fnstripstantrkpeak_dt[ih]->GetParError(0),5)<<endl;

      gnstripstantrkdeco_dt[ih] = getTGraphFromTH2(nstripstantrkdeco_dt[ih],bins2,0);
      fnstripstantrkdeco_dt[ih]->SetLineColor(4);
      fnstripstantrkdeco_dt[ih]->SetLineWidth(1);
      gnstripstantrkdeco_dt[ih]->Fit(fnstripstantrkdeco_dt[ih],"R");
      gnstripstantrkdeco_dt[ih]->SetMarkerColor(4);
      gnstripstantrkdeco_dt[ih]->SetLineColor(4);
      gnstripstantrkdeco_dt[ih]->SetMarkerStyle(8);
      gnstripstantrkdeco_dt[ih]->SetMarkerSize(0.1);
      gnstripstantrkdeco_dt[ih]->Draw("sameP");

      tanladeco[ih]    = fnstripstantrkdeco_dt[ih]->GetParameter(0);
      tanlaerrdeco[ih] = fnstripstantrkdeco_dt[ih]->GetParError(0);
      dtdeco[ih]       = dttimedeco_dt[ih]->GetMean(1);
      dterrdeco[ih]    = dttimedeco_dt[ih]->GetRMS(1)/sqrt(dttimedeco_dt[ih]->GetEntries());
      
      snstripstantrkdeco_dt[ih]<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkdeco_dt[ih]->GetParameter(0),5)
                               <<" #pm "<<fround(fnstripstantrkdeco_dt[ih]->GetParError(0),5)<<endl;
      
      
      TLatex *t=new TLatex();
      t->SetNDC();
      t->SetTextSize(0.04);
      t->SetTextColor(2);
      t->DrawLatex(0.3,0.82,snstripstantrkpeak_dt[ih].str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(0.3,0.75,snstripstantrkdeco_dt[ih].str().c_str());
      t->SetTextColor(1);
      t->DrawLatex(0.3,0.68,func);
    }
    
    idx++;
  }

  if(draw("tanladt",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),800,600);
    can[idx]->cd();

    TGraphErrors *glapeak=new TGraphErrors(8,dtpeak,tanlapeak,dterrpeak,tanlaerrpeak);
    TGraphErrors *gladeco=new TGraphErrors(8,dtdeco,tanladeco,dterrdeco,tanlaerrdeco);
    
    glapeak->SetMarkerColor(2);
    glapeak->SetLineColor(2);
    glapeak->SetMarkerStyle(8);
    glapeak->SetMarkerSize(1);
    glapeak->Draw("AP");
    glapeak->SetTitle("tan(#theta_{LA}) from Cluster Width");
    glapeak->GetYaxis()->SetTitle("tan(#theta_{LA})");
    glapeak->GetXaxis()->SetTitle("DT time (ns)");
    glapeak->SetMinimum(-0.11);
    gladeco->SetMaximum(-0.05);

    gladeco->SetMarkerColor(4);
    gladeco->SetLineColor(4);
    gladeco->SetMarkerStyle(8);
    gladeco->SetMarkerSize(1);
    gladeco->Draw("sameP");
   
    idx++;
  }
  */

  /*
  if(draw("duvsdtantheta_tgraph",cantitles)){
    bool addbpcorr = false; //add tgraph for BP-corrected deco
    float size = 0.1;

    vector<float> thetabins;
    //         thetabins.push_back(-1.);
    //         thetabins.push_back(-0.75);
    //         thetabins.push_back(-0.5);
    //         thetabins.push_back(-0.25);
    //         thetabins.push_back(0.);
    //         thetabins.push_back(0.25);
    //         thetabins.push_back(0.5);
    //         thetabins.push_back(0.75);
    //         thetabins.push_back(1.);

    thetabins.push_back(-0.9);
    thetabins.push_back(-0.7);
    thetabins.push_back(-0.5);
    thetabins.push_back(-0.3);
    thetabins.push_back(-0.1);
    thetabins.push_back(0.1);
    thetabins.push_back(0.3);
    thetabins.push_back(0.5);
    thetabins.push_back(0.7);
    thetabins.push_back(0.9);

    
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
  
    can[idx]->cd(1);
    
    duthetapeakp->SetName("duthetapeakp");
    TGraphErrors *gduthetapeakp = getTGraphFromTH2(duthetapeakp,thetabins,0);
    gduthetapeakp->SetMarkerColor(2);
    gduthetapeakp->SetLineColor(2);
    gduthetapeakp->SetMarkerStyle(8);
    gduthetapeakp->SetMarkerSize(size);
    TF1 *fduthetapeakp=new TF1("fduthetapeakp","pol1");
    fduthetapeakp->SetLineColor(2);
    gduthetapeakp->Fit(fduthetapeakp);
    stringstream sduthetapeakp;
    sduthetapeakp<<"m = "<<fround(fduthetapeakp->GetParameter(1),3)
                 <<" #pm " <<fround(fduthetapeakp->GetParError(1),3)
                 <<"   b = "<<fround(fduthetapeakp->GetParameter(0),3)
                 <<" #pm "<<fround(fduthetapeakp->GetParError(0),3)
                 <<endl;
    
    duthetadecop->SetName("duthetadecop");
    TGraphErrors *gduthetadecop = getTGraphFromTH2(duthetadecop,thetabins,0);
    gduthetadecop->SetMarkerColor(4);
    gduthetadecop->SetLineColor(4);
    gduthetadecop->SetMarkerStyle(8);
    gduthetadecop->SetMarkerSize(size);
    TF1 *fduthetadecop=new TF1("fduthetadecop","pol1");
    fduthetadecop->SetLineColor(4);
    gduthetadecop->Fit(fduthetadecop);
    stringstream sduthetadecop;
    sduthetadecop<<"m = "<<fround(fduthetadecop->GetParameter(1),3)
                 <<" #pm " <<fround(fduthetadecop->GetParError(1),3)
                 <<"   b = "<<fround(fduthetadecop->GetParameter(0),3)
                 <<" #pm "<<fround(fduthetadecop->GetParError(0),3)
                 <<endl;

    //add deco with BP correction------------------------------------------
    if(addbpcorr){
      TFile fdecobp("crabjobs/trial17/root/deco.root");
      fdecobp.cd();
      PeakDecoResiduals.cd();
      PeakDecoResiduals.cd();
      TH2F* duthetadecobp_vp = (TH2F*)duvsdtantheta_vp_TOB->Clone();
      TH2F* duthetadecobp_vm = (TH2F*)duvsdtantheta_vm_TOB->Clone();
      
      duthetadecobp_vp->SetName("duthetadecobp_vp");
      TGraphErrors *gduthetadecobp_vp = getTGraphFromTH2(duthetadecobp_vp,thetabins,0);
      gduthetadecobp_vp->SetMarkerColor(6);
      gduthetadecobp_vp->SetLineColor(6);
      gduthetadecobp_vp->SetMarkerStyle(8);
      gduthetadecobp_vp->SetMarkerSize(size);
      TF1 *fduthetadecobp_vp=new TF1("fduthetadecobp_vp","pol1");
      fduthetadecobp_vp->SetLineColor(6);
      gduthetadecobp_vp->Fit(fduthetadecobp_vp);
      stringstream sduthetadecobp_vp;
    
      sduthetadecobp_vp<<"m = "<<fround(fduthetadecobp_vp->GetParameter(1),3)
                       <<" #pm " <<fround(fduthetadecobp_vp->GetParError(1),3)
                       <<"   b = "<<fround(fduthetadecobp_vp->GetParameter(0),3)
                       <<" #pm "<<fround(fduthetadecobp_vp->GetParError(0),3)
                       <<endl;
    }

    TMultiGraph *mgdtp=new TMultiGraph();
    mgdtp->Add(gduthetapeakp);
    mgdtp->Add(gduthetadecop);
    if(addbpcorr) mgdtp->Add(gduthetadecobp_vp);
    mgdtp->SetTitle("TOB (v+)");
    mgdtp->Draw("AP");
    mgdtp->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
    mgdtp->GetYaxis()->SetTitle("<#Deltau> [#mum]");
    mgdtp->SetTitle("TOB");
    //mgdtp->GetXaxis()->SetLimits(thetabins.at(0),thetabins.at(thetabins.size())-1);
    //mgdtp->GetXaxis()->SetLimits(-0.25,0.25);
    //mgdtp->GetYaxis()->SetLimits(-30,10);
    mgdtp->SetMinimum(-20);
    mgdtp->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    //t->DrawLatex(0.5,0.8,"y = mx + b");
    t->SetTextColor(2);
    t->DrawLatex(0.25,0.85,sduthetapeakp.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.25,0.75,sduthetadecop.str().c_str());
    if(addbpcorr){
      t->SetTextColor(6);
      t->DrawLatex(0.25,0.65,sduthetadecobp_vp.str().c_str());
    }

    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);
    
    can[idx]->cd(2);
    duthetapeakm->SetName("duthetapeakm");
    TGraphErrors *gduthetapeakm = getTGraphFromTH2(duthetapeakm,thetabins,0);
    gduthetapeakm->SetMarkerColor(2);
    gduthetapeakm->SetLineColor(2);
    gduthetapeakm->SetMarkerStyle(8);
    gduthetapeakm->SetMarkerSize(size);
    TF1 *fduthetapeakm=new TF1("fduthetapeakm","pol1");
    fduthetapeakm->SetLineColor(2);
    gduthetapeakm->Fit(fduthetapeakm);
    stringstream sduthetapeakm;
    sduthetapeakm<<"m = "<<fround(fduthetapeakm->GetParameter(1),3)
                 <<" #pm " <<fround(fduthetapeakm->GetParError(1),3)
                 <<"   b = "<<fround(fduthetapeakm->GetParameter(0),3)
                 <<" #pm "<<fround(fduthetapeakm->GetParError(0),3)
                 <<endl;

    duthetadecom->SetName("duthetadecom");
    TGraphErrors *gduthetadecom = getTGraphFromTH2(duthetadecom,thetabins,0);
    gduthetadecom->SetMarkerColor(4);
    gduthetadecom->SetLineColor(4);
    gduthetadecom->SetMarkerStyle(8);
    gduthetadecom->SetMarkerSize(size);
    TF1 *fduthetadecom=new TF1("fduthetadecom","pol1");
    fduthetadecom->SetLineColor(4);
    gduthetadecom->Fit(fduthetadecom);
    stringstream sduthetadecom;
    
    sduthetadecom<<"m = "<<fround(fduthetadecom->GetParameter(1),3)
                 <<" #pm " <<fround(fduthetadecom->GetParError(1),3)
                 <<"   b = "<<fround(fduthetadecom->GetParameter(0),3)
                 <<" #pm "<<fround(fduthetadecom->GetParError(0),3)
                 <<endl;
    
    //add deco with BP correction------------------------------------------
    if(addbpcorr){
      duthetadecobp_vm->SetName("duthetadecobp_vm");
      TGraphErrors *gduthetadecobp_vm = getTGraphFromTH2(duthetadecobp_vm,thetabins,0);
      gduthetadecobp_vm->SetMarkerColor(6);
      gduthetadecobp_vm->SetLineColor(6);
      gduthetadecobp_vm->SetMarkerStyle(8);
      gduthetadecobp_vm->SetMarkerSize(size);
      TF1 *fduthetadecobp_vm=new TF1("fduthetadecobp_vm","pol1");
      fduthetadecobp_vm->SetLineColor(6);
      gduthetadecobp_vm->Fit(fduthetadecobp_vm);
      stringstream sduthetadecobp_vm;
      
      sduthetadecobp_vm<<"m = "<<fround(fduthetadecobp_vm->GetParameter(1),3)
                       <<" #pm " <<fround(fduthetadecobp_vm->GetParError(1),3)
                       <<"   b = "<<fround(fduthetadecobp_vm->GetParameter(0),3)
                       <<" #pm "<<fround(fduthetadecobp_vm->GetParError(0),3)
                       <<endl;
    }

    TMultiGraph *mgdtm=new TMultiGraph();
    mgdtm->Add(gduthetapeakm);
    mgdtm->Add(gduthetadecom);
    if(addbpcorr) mgdtm->Add(gduthetadecobp_vm);
    mgdtm->SetTitle("TOB (v-)");
    mgdtm->Draw("AP");
    mgdtm->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
    mgdtm->GetYaxis()->SetTitle("<#Deltau> [#mum]");
    mgdtm->SetTitle("TOB");
    //mgdtm->GetXaxis()->SetLimits(thetabins.at(0),thetabins.at(thetabins.size())-1);
    //mgdtm->GetXaxis()->SetLimits(-0.25,0.25);
    //mgdtm->GetYaxis()->SetLimits(-10,30);
    mgdtm->SetMinimum(-20);
    mgdtm->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    //t->DrawLatex(0.5,0.8,"y = mx + b");
    t->SetTextColor(2);
    t->DrawLatex(0.25,0.85,sduthetapeakm.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.25,0.75,sduthetadecom.str().c_str());
    if(addbpcorr){
      t->SetTextColor(6);
      t->DrawLatex(0.25,0.65,sduthetadecobp_vm.str().c_str());
    }


    float bpeakp = fduthetapeakp->GetParameter(0);
    float bpeakm = fduthetapeakm->GetParameter(0);
    float mpeakp = fduthetapeakp->GetParameter(1);
    float mpeakm = fduthetapeakm->GetParameter(1);
    float bdecop = fduthetadecop->GetParameter(0);
    float bdecom = fduthetadecom->GetParameter(0);
    float mdecop = fduthetadecop->GetParameter(1);
    float mdecom = fduthetadecom->GetParameter(1);

    float dw     = 0.5*( (mdecop - mpeakp) + (mdecom - mpeakm) );
    float db     = 0.5*( (bdecop - bpeakp) - (bdecom - bpeakm) );
    float dtanla = db/(235-dw);
    
    cout<<"dw "<<dw<<" dtanla "<<dtanla<<endl;

    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);
    idx++;
  }
*/


              
  if(printgif){
    for( int ican = 0 ; ican < idx ; ican++ ){
      can[ican]->Modified();
      can[ican]->Update(); 
      can[ican]->Print(Form( "plots/%s.gif",can[ican]->GetTitle() ));
    }
  }
              
  if(writeTFile){
    ofile->cd();
    ofile->Write();
    ofile->Close();
  }
}



void plotHists(TH1F** h, const unsigned int nhist, string title, string xtitle, float xmin, float xmax, int rebin, int fit){
   
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
    h[i]->SetTitle( title.c_str() );
    h[i]->GetXaxis()->SetTitle( xtitle.c_str() );
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

bool draw(string var, vector<string> cantitles){
  
  for(unsigned int i=0;i<cantitles.size();i++) {
    if(strcmp(var.c_str(),cantitles.at(i).c_str()) == 0) return true;
  }
  return false;
}

int getDTSlice(float t){

  int slice = -1;
  if(t > -20  && t < -12)  slice = 0;
  if(t > -12  && t < -8 )  slice = 1;
  if(t > -8   && t < -4 )  slice = 2;
  if(t > -4   && t < 0  )  slice = 3;
  if(t >  0   && t < 4  )  slice = 4;
  if(t >  4   && t < 8  )  slice = 5;
  if(t >  8   && t < 12 )  slice = 6;
  if(t >  12  && t < 20 )  slice = 7;

  return slice;
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

  //   r->Draw("sames");
//     gPad->Update(); 
//     TPaveStats* st2 = (TPaveStats*) r->GetListOfFunctions()->FindObject("stats");
//     st2->SetX1NDC(startingX);
//     st2->SetX2NDC(startingX+0.30);
//     st2->SetY1NDC(startingY);
//     st2->SetY2NDC(startingY+height);
//     st2->SetTextColor(4);
  }
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

TGraphErrors *getTGraphFromTH2(TH2F* h,vector<float> xbins, int method, int invert, int color){

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
      if(invert > 0) y[ibin] = -1 * y[ibin];
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
  g -> SetMarkerColor( color );
  g -> SetLineColor( color );
  g -> SetMarkerStyle( 8 );
  g -> SetMarkerSize( 0.5 );
  g->GetXaxis()->SetLimits(xbins.at(0),xbins.at(nbins));
  
  return g;
}


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


string getStringFromTF1( TF1* f){
  
  stringstream s;
  float dw  = f->GetParameter(1);
  float b   = f->GetParameter(0);
  float dwe = f->GetParError(1);
  float be  = f->GetParError(0);
  
  s <<"#DeltaW "         << fround(dw,2)
    <<" #pm "            << fround(dwe,2)
    <<" #Deltatan(LA) "  << fround(b/(235-dw),4)
    <<" #pm "            << fround(be/(235-dw),4) <<endl; 
  
  return s.str();
}
