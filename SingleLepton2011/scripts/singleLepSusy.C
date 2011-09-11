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
#include "Hootilities.h"
#include "singleLepSusy.h"
//#include "histtools.h"

using namespace std;

bool printgif_           = false;
bool alreadyInitialized_ = false;


TH1F* getHist( TChain* ch , char* var , TCut sel , char* histname , int nbins , float xmin , float xmax ){

  TH1F* h = new TH1F(histname,histname,nbins,xmin,xmax);
  h->Sumw2();
  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("TMath::Min(%s,%f) >> %s",var,xmax-0.001,histname),sel);
  delete ctemp;
  return h;
}

TH2F* getHist2D( TChain* ch , char* varx , char* vary , TCut sel , char* histname , 
		 int nbinsx , float xmin , float xmax , int nbinsy , float ymin , float ymax ){
  
  TH2F* h = new TH2F(histname,histname,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  h->Sumw2();

  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("TMath::Min(%s,%f):TMath::Min(%s,%f) >> %s",vary,ymax-0.001,varx,xmax-0.01,histname),sel);
  delete ctemp;
  return h;
}


float calculateHistError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

//--------------------------------------------------
// initialize data/MC samples
//--------------------------------------------------

void initialize(char* path){

  if( alreadyInitialized_ ){

    cout << "Resetting babies" << endl;

    data->Reset();
    ttall->Reset();
    ttl->Reset();
    ttll->Reset();
    ttltau->Reset();
    tttau->Reset();
    tttautau->Reset();
    ttotr->Reset();
    wjets->Reset();
    t2tt->Reset();

    mc.clear();
    mctex.clear();
    mclabels.clear();
  }

  else{

    data	= new TChain("t");
    ttall	= new TChain("t");
    ttl  	= new TChain("t");
    ttll	= new TChain("t");
    ttltau	= new TChain("t");
    tttau	= new TChain("t");
    tttautau	= new TChain("t");
    ttotr	= new TChain("t");
    wjets	= new TChain("t");
    t2tt	= new TChain("t");
  }

  cout << endl;
  cout << "Loading babies at       : " << path << endl;
  
  data->Add(Form("%s/data_smallTree.root",path));
  ttall->Add(Form("%s/ttall_smallTree.root",path));
  ttl->Add(Form("%s/ttl_smallTree.root",path));
  ttll->Add(Form("%s/ttll_smallTree.root",path));
  ttltau->Add(Form("%s/ttltau_smallTree.root",path));
  tttau->Add(Form("%s/tttau_smallTree.root",path));
  tttautau->Add(Form("%s/tttautau_smallTree.root",path));
  ttotr->Add(Form("%s/ttotr_smallTree.root",path));
  wjets->Add(Form("%s/wjets_smallTree.root",path));
  t2tt->Add(Form("%s/T2tt_smallTree.root",path));
  
  //mc.push_back(ttall);       mclabels.push_back("ttall");    
  mc.push_back(ttl);         mclabels.push_back("ttl");    
  mc.push_back(ttll);        mclabels.push_back("ttll");    
  mc.push_back(ttltau);      mclabels.push_back("ttltau");   
  mc.push_back(tttau);       mclabels.push_back("tttau");   
  mc.push_back(tttautau);    mclabels.push_back("tttautau");   
  mc.push_back(ttotr);       mclabels.push_back("ttotr");   
  mc.push_back(wjets);       mclabels.push_back("wjets");   
  //mc.push_back(t2tt);        mclabels.push_back("t2tt");   
  
  alreadyInitialized_ = true;
}

//------------------------------------------
// selection and weight to apply to babies
//------------------------------------------

TCut selection_TCut(){

  TCut njets1("njets >= 1");
  TCut njets2("njets >= 2");
  TCut njets3("njets >= 3");
  TCut njets4("njets >= 4");
  TCut met60("pfmet > 60");
  TCut met100("pfmet > 100");
  TCut ht300("ht > 300");
  TCut ht500("ht > 500");
  TCut btags0("nbtags==0");
  TCut btags2("nbtags>=2");

  TCut sel;
  sel    += njets3;
  sel    += met60;
  //sel    += met100;
  sel    += ht300;
  //sel    += ht500;
  //sel    += btags0;
  //sel    += btags2;

  cout << "Using selection         : " << sel.GetTitle() << endl;
 
  return sel;
}

TCut weight_TCut(){

  TCut weight("weight*ndavtxweight");

  cout << "Using weight            : " << weight.GetTitle() << endl;
  return weight; 
}


void plotHist( TH1F* h1 , TH1F* h2 , char* leg1 , char* leg2 , char* xtitle , bool residual){

  h1->Scale( 1. / h1->Integral() );
  h2->Scale( 1. / h2->Integral() );

  float max = h1->GetMaximum();
  if( h2->GetMaximum() > max ) max = h2->GetMaximum();
  h1->SetMaximum( 1.2 * max );

  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();

  if( residual ){
    fullpad = new TPad("fullpad","fullpad",0,0,1,1);
    fullpad->Draw();
    fullpad->cd();

    plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
    plotpad->Draw();
    plotpad->cd();
  }

  //gPad->SetGridx();
  //gPad->SetGridy();
  //gPad->SetLogy(1);
  
  h1->SetMarkerColor(4);
  h1->SetLineColor(4);
  h1->SetFillColor(4);
  h1->SetMarkerStyle(25);

  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetFillColor(2);
  h2->SetMarkerStyle(20);

  h1->GetXaxis()->SetTitle( xtitle );
  h1->DrawNormalized("E1");
  h2->DrawNormalized("sameE1");

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.DrawLine(100,0,100,1.2*max);

  TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
  leg->AddEntry(h1,leg1,"p");
  leg->AddEntry(h2,leg2,"p");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->Draw();

  int bin    = h1->FindBin(100);
  float eff1 = h1->Integral(bin,1000)/h1->Integral();
  float eff2 = h2->Integral(bin,1000)/h2->Integral();

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  text->DrawLatex(0.6,0.6,"njets #geq 3");
  text->SetTextColor(4);
  text->DrawLatex(0.6,0.5,Form("eff(M_{T}>100 GeV) = %.3f",eff1));	      
  text->SetTextColor(2);
  text->DrawLatex(0.6,0.4,Form("eff(M_{T}>100 GeV) = %.3f",eff2));

  if( residual ){
    fullpad->cd();

    respad = new TPad("respad","respad",0,0.8,1,1);
    respad->Draw();
    respad->cd();

    gPad->SetGridy();

    TH1F* ratio = (TH1F*) h2->Clone(Form("%s_ratio",h2->GetName()));
    ratio->Divide(h1);

    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetLabelSize(0.2);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetTitle("ratio    ");
    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetXaxis()->SetTitleSize(0);
    ratio->SetMarkerSize(0.7);
    ratio->Draw();

    TLine myline;
    myline.SetLineWidth(1);
    myline.DrawLine(h1->GetXaxis()->GetXmin(),1,h1->GetXaxis()->GetXmax(),1);

  }

}


//------------------------------------
// print yield table
//------------------------------------

void printYieldTable( char* path , bool latex = false ){

  gROOT->Reset();
  deleteHistos();

  initialize(path);
  initSymbols(latex);

  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  printYields( mc , mclabels , data , sel , weight , latex );
 
}

//--------------------------------------------------
// make data/MC plots
//--------------------------------------------------

void makePlots( char* path , bool printgif = false ){

  bool combine4 = false;

  initialize(path);

  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  vars.push_back("pfmet");     xt.push_back("pfmet (GeV)");      n.push_back(5); xi.push_back(0.); xf.push_back(250.);
  vars.push_back("ht");        xt.push_back("H_{T} (GeV)");      n.push_back(7); xi.push_back(0.); xf.push_back(700.);
  vars.push_back("dilmass");   xt.push_back("m(ll) (GeV)");      n.push_back(5); xi.push_back(0.); xf.push_back(100.);
  vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(5); xi.push_back(0.); xf.push_back(100.);
    
  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  bool residual = false;

  int canCounter = -1;
  bool log = false;

  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    if( ivar < 2 ) log = true;
    else           log = false;

    if( combine4 ){
      if( ivar % 4 == 0 ){
	canCounter++;
	can[canCounter] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),1400,1200);

	legpad[canCounter] = new TPad("legpad","legpad",12./14.,0,1,1);
	legpad[canCounter]->Draw();
	legpad[canCounter]->cd();

	TLegend *leg = getLegend( mc , mclabels , true , 0.2 , 0.3 , 0.8 , 0.7 );
	leg->SetTextSize(0.1);
	leg->SetBorderSize(1);
	leg->Draw();

	can[canCounter]->cd();

	plotpad[canCounter] = new TPad("plotpad","plotpad",0,0,12./14.,1);
	plotpad[canCounter]->Draw();
	plotpad[canCounter]->cd();

	plotpad[canCounter]->Divide(2,2);
	plotpad[canCounter]->cd(1);

      }else{
	plotpad[canCounter]->cd(1+ivar%4);
      }
    }else{
      can[ivar] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),600,600);
    }

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , !combine4 , log );

    if( printgif ) can[ivar]->Print(Form("../plots/%s.pdf",vars[ivar]));
  } 
}

void makeStandardPlots( char* path , bool sigregion = false ){

  bool residual = false;
  bool log      = false;

  cout << "Plot residual? " << residual << endl;
  cout << "Do log plot?   " << log      << endl;

  deleteHistos();

  initialize(path);

  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  if( sigregion ){
    TCut highmet = "pfmet>275 && htpf>300";
    TCut highht = "pfmet>200 && htpf>600";

    sel = sel + ( highmet || highht );
    cout << "Signal region: " << sel.GetTitle() << endl;
  }

  char* filename;
  if(  sigregion ) filename = "datamc_sig";
  else             filename = "datamc";

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  vars.push_back("lep1.pt()");  xt.push_back("lepton p_{T} (GeV)");	    n.push_back(20); xi.push_back(0.);   xf.push_back(200.);
  vars.push_back("lep1.eta()"); xt.push_back("lepton #eta")       ;	    n.push_back(20); xi.push_back(-3.);  xf.push_back(3.);
  vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	    n.push_back(20); xi.push_back(0.);   xf.push_back(400.);
  vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	            n.push_back(20); xi.push_back(-3);   xf.push_back( 3);
  vars.push_back("dphijm");     xt.push_back("#Delta#phi(max jet,pfmet)");  n.push_back(20); xi.push_back(0);    xf.push_back(3.2);
  vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");    	    n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");		    n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("y");          xt.push_back("y (GeV^{1/2})");  	    n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
  vars.push_back("ht");         xt.push_back("H_{T} (GeV)");		    n.push_back(20); xi.push_back(0.);   xf.push_back(1000.);
  vars.push_back("njets");      xt.push_back("njets");			    n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
  vars.push_back("nbtags");     xt.push_back("nbtags");			    n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
  vars.push_back("ndavtx");     xt.push_back("nDAVertices");		    n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
  vars.push_back("mt");         xt.push_back("M_{T}");			    n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("meff");       xt.push_back("meff");			    n.push_back(20); xi.push_back(0.);   xf.push_back(2000.);
   
  const unsigned int nvars = vars.size();
  
  //TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  //int canCounter = -1;

  TCanvas* canvas = new TCanvas("canvas","canvas",1100,750);
  gStyle->SetPaperSize(22,28);
  canvas->Print(Form("../plots/%s.ps[",filename));

  TLegend *leg = getLegend( mc , mclabels , true , 0.1 , 0.3 , 0.8 , 0.7 );
  leg->SetTextSize(0.1);
  leg->SetBorderSize(1);
  
  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    //can[ivar] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),1000,800);
    canvas->cd();

    legpad[ivar] = new TPad("legpad","legpad",0.8,0,1,1);
    legpad[ivar]->Draw();
    legpad[ivar]->cd();

    leg->Draw();

    //can[ivar]->cd();
    canvas->cd();

    plotpad[ivar] = new TPad("plotpad","plotpad",0,0,0.8,1);
    plotpad[ivar]->Draw();
    plotpad[ivar]->cd();
    
    plotpad[ivar]->Divide(2,2);

    plotpad[ivar]->cd(1);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==0") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "e"   );
    plotpad[ivar]->cd(2);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==1") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "m"   );
    //plotpad[ivar]->cd(3);
    //compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==2") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "em"  );
    plotpad[ivar]->cd(4);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel)              , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "all" );

    canvas->Print(Form("../plots/%s.ps",filename));
    canvas->Clear();

    //if( printgif ) can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
  } 

  canvas->Print(Form("../plots/%s.ps]",filename));
  canvas->Clear();
  
  gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.ps ../plots/%s.pdf",filename,filename));
}

void mt( char* path , bool latex = false ){

  gROOT->Reset();
  deleteHistos();

  initialize(path);
  initSymbols(latex);

  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  TCut njets3("njets >= 3");
  TCut ht1("ht<200");
  TCut ht2("ht>200 && ht<300");
  TCut ht3("ht>300");

  TH1F* w1 = getHist( wjets , "mt" , TCut(njets3) , "hw1" , 50 , 0 , 200 );
  TH1F* t1 = getHist( ttl   , "mt" , TCut(njets3) , "ht1" , 50 , 0 , 200 );

  TCanvas *c1 = new TCanvas();
  c1->cd();
  plotHist( t1 , w1 , "t#bar{t}#rightarrow l+jets" , "W+jets" , "M_{T} (GeV)" , true );  

}
