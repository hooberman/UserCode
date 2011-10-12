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
//#include "Hootilities.h"
#include "Hootilities.C"
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
    tt0l->Reset();
    tt1l->Reset();
    tt2l->Reset();
    ttl->Reset();
    ttll->Reset();
    ttltau->Reset();
    tttau->Reset();
    tttautau->Reset();
    ttotr->Reset();
    wjets->Reset();
    t2tt->Reset();
    t2ttA->Reset();
    t2ttB->Reset();
    t2ttC->Reset();
    qcd->Reset();
    tW->Reset();

    mc.clear();
    mctex.clear();
    mclabels.clear();
  }

  else{
    data	= new TChain("t");
    ttall	= new TChain("t");
    tt0l	= new TChain("t");
    tt1l	= new TChain("t");
    tt2l	= new TChain("t");
    ttl  	= new TChain("t");
    ttll	= new TChain("t");
    ttltau	= new TChain("t");
    tttau	= new TChain("t");
    tttautau	= new TChain("t");
    ttotr	= new TChain("t");
    wjets	= new TChain("t");
    t2tt	= new TChain("t");
    t2ttA	= new TChain("t");
    t2ttB	= new TChain("t");
    t2ttC	= new TChain("t");
    qcd	        = new TChain("t");
    tW	        = new TChain("t");
  }

  cout << endl;
  cout << "Loading babies at       : " << path << endl;
  
  data->Add(Form("%s/data_smallTree.root",path));
  ttall->Add(Form("%s/ttall_smallTree.root",path));
  tt0l->Add(Form("%s/tt0l_smallTree.root",path));
  tt1l->Add(Form("%s/tt1l_smallTree.root",path));
  tt2l->Add(Form("%s/tt2l_smallTree.root",path));
  ttl->Add(Form("%s/ttl_smallTree.root",path));
  ttll->Add(Form("%s/ttll_smallTree.root",path));
  ttltau->Add(Form("%s/ttltau_smallTree.root",path));
  tttau->Add(Form("%s/tttau_smallTree.root",path));
  tttautau->Add(Form("%s/tttautau_smallTree.root",path));
  ttotr->Add(Form("%s/ttotr_smallTree.root",path));
  wjets->Add(Form("%s/wjets_smallTree.root",path));
  qcd->Add(Form("%s/qcd_smallTree.root",path));
  tW->Add(Form("%s/tW_smallTree.root",path));
  //t2tt->Add(Form("%s/T2tt_few_smallTree.root",path));
  t2ttA->Add(Form("%s/T2tt_350_100_smallTree.root",path));
  t2ttB->Add(Form("%s/T2tt_450_100_smallTree.root",path));
  t2ttC->Add(Form("%s/T2tt_200_50_smallTree.root",path));
  
  //mc.push_back(ttall);       mclabels.push_back("ttall");    
  mc.push_back(tt1l);        mclabels.push_back("tt1l");    
  mc.push_back(tt2l);        mclabels.push_back("tt2l");    
  mc.push_back(tt0l);        mclabels.push_back("ttfake");    
  // mc.push_back(ttl);         mclabels.push_back("ttl");    
  // mc.push_back(ttll);        mclabels.push_back("ttll");    
  // mc.push_back(ttltau);      mclabels.push_back("ttltau");   
  // mc.push_back(tttau);       mclabels.push_back("tttau");   
  // mc.push_back(tttautau);    mclabels.push_back("tttautau");   
  // mc.push_back(ttotr);       mclabels.push_back("ttotr");   
  mc.push_back(wjets);       mclabels.push_back("wjets");   
  mc.push_back(qcd);         mclabels.push_back("qcd");   
  mc.push_back(tW);          mclabels.push_back("single top");   
  //mc.push_back(t2tt);        mclabels.push_back("T2tt");   
  mc.push_back(t2ttA);       mclabels.push_back("T2tt 350/100 X6");   
  mc.push_back(t2ttB);       mclabels.push_back("T2tt 450/100 X6");   
  //mc.push_back(t2ttC);       mclabels.push_back("T2tt 4.0");   
  
  alreadyInitialized_ = true;
}

//------------------------------------------
// selection and weight to apply to babies
//------------------------------------------

TCut selection_TCut(){

  TCut nlep1("ngoodlep == 1");
  TCut leppt("(leptype==0 && lep1.pt()>25)||(leptype==1 && lep1.pt()>20)");
  // TCut njets1("njets >= 1");
  // TCut njets2("njets >= 2");
  // TCut njets3("njets >= 3");
  // TCut njets4("njets >= 4");
  TCut njets1("ncalojets >= 1");
  TCut njets2("ncalojets >= 2");
  TCut njets3("ncalojets >= 3");
  TCut njets4("ncalojets >= 4");
  TCut met60("pfmet > 60");
  TCut met50("pfmet > 50");
  TCut met100("pfmet > 100");
  TCut metpresel("(leptype==0&&pfmet>30)||(leptype==1&&pfmet>20)");
  TCut ht300("ht > 300");
  TCut ht500("ht > 500");
  TCut dphi05("dphijm > 0.5");
  TCut btags0("nbtags==0");
  TCut btags2("nbctcm>=2");
  TCut btags1("nbctcm>=1");
  TCut mt100("mt>100");
  TCut mt150("mt>150");
  TCut mt200("mt>200");
  TCut trkreliso02("trkreliso10>0.2");
  TCut trkreliso01("trkreliso5>0.1");

  TCut sel;
  sel    += nlep1;
  sel    += leppt;
  sel    += njets3;
  sel    += met50;
  sel    += met100;
  //sel    += met60;
  //sel    += dphi05;
  //sel    += met100;
  //sel    += ht300;
  //sel    += ht500;
  //sel    += btags0;
  //sel    += btags1;
  sel    += btags2;
  //sel    += mt150;
  //sel    += mt150;
  //sel    += mt200;
  sel    += trkreliso01;
  //sel    += trkreliso02;

  cout << "Using selection         : " << sel.GetTitle() << endl;
 
  return sel;
}

TCut weight_TCut(){

  //TCut weight("weight*ndavtxweight");
  //TCut weight("weight * 0.204 * ndavtxweight");
  TCut weight("weight * 0.98 * ndavtxweight * 1.13");
  //TCut weight("weight * 5 * ndavtxweight");

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
  //text->DrawLatex(0.6,0.6,"njets #geq 3");
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

  bool combine6 = false;
  bool residual = false;
  bool log      = false;


  initialize(path);

  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  // vars.push_back("pfmet");        xt.push_back("pfmet (GeV)");      n.push_back(30);  xi.push_back(0.); xf.push_back(300.);
  // vars.push_back("htcalo");       xt.push_back("H_{T} (GeV)");      n.push_back(30);  xi.push_back(0.); xf.push_back(700.);
  // vars.push_back("ncalojets");    xt.push_back("H_{T} (GeV)");      n.push_back(10);  xi.push_back(0.); xf.push_back(10.);
  // vars.push_back("lep1.pt()");    xt.push_back("pfmet (GeV)");      n.push_back(5);   xi.push_back(0.); xf.push_back(250.);
  // vars.push_back("mt");           xt.push_back("H_{T} (GeV)");      n.push_back(7);   xi.push_back(0.); xf.push_back(700.);
  // vars.push_back("nbctcm");       xt.push_back("H_{T} (GeV)");      n.push_back(4);   xi.push_back(0.); xf.push_back(4.);

  // vars.push_back("pfmet");        xt.push_back("E_{T}^{miss} (GeV)");  	            n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  // vars.push_back("htcalo");       xt.push_back("H_{T} (GeV)");		            n.push_back(20); xi.push_back(0.);   xf.push_back(1000.);
  // vars.push_back("ncalojets");    xt.push_back("jet multiplicity");	            n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
  // vars.push_back("lep1.pt()");    xt.push_back("lepton p_{T} (GeV)");	            n.push_back(20); xi.push_back(0.);   xf.push_back(200.);
  // vars.push_back("mt");           xt.push_back("M_{T} (GeV)");		            n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  // vars.push_back("nbctcm");       xt.push_back("b-jet multiplicity");	            n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
  

  //vars.push_back("trkreliso5");    xt.push_back("track reliso");		            n.push_back(50); xi.push_back(0.);     xf.push_back(1.);
  //vars.push_back("mt");           xt.push_back("M_{T} (GeV)");		            n.push_back(16); xi.push_back(0.);     xf.push_back(400.);
  vars.push_back("mt");           xt.push_back("M_{T} (GeV)");		            n.push_back(11); xi.push_back(125.);   xf.push_back(400.);
  //vars.push_back("pfmet");        xt.push_back("E_{T}^{miss} (GeV)");  	            n.push_back(11); xi.push_back(70.);    xf.push_back(400.);

  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];


  int canCounter = -1;

  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     
    
    //if( ivar < 2 ) log = true;
    //else           log = false;

    if( combine6 ){
      if( ivar % 6 == 0 ){
	canCounter++;
	can[canCounter] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),2000,1200);
	
	legpad[canCounter] = new TPad("legpad","legpad",18./20.,0,1,1);
	legpad[canCounter]->Draw();
	legpad[canCounter]->cd();

	TLegend *leg = getLegend( mc , mclabels , true , 0.2 , 0.3 , 0.8 , 0.7 );
	leg->SetTextSize(0.1);
	leg->SetBorderSize(1);
	leg->Draw();

	can[canCounter]->cd();

	plotpad[canCounter] = new TPad("plotpad","plotpad",0,0,18./20.,1);
	plotpad[canCounter]->Draw();
	plotpad[canCounter]->cd();

	plotpad[canCounter]->Divide(3,2);
	plotpad[canCounter]->cd(1);

      }else{
	plotpad[canCounter]->cd(1+ivar%6);
      }
    }else{
      can[ivar] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),600,600);
    }

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , !combine6 , log );

    if( printgif && !combine6 ){
      can[ivar]->Print(Form("../plots/%s.pdf",vars[ivar]));
      can[ivar]->Print(Form("../plots/%s.eps",vars[ivar]));
      can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
    }
  } 

  if( printgif && combine6 ){
    can[0]->Print("../plots/makePlots.pdf");
    can[0]->Print("../plots/makePlots.eps");
    can[0]->Print("../plots/makePlots.png");
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

  vars.push_back("lep1.pt()");    xt.push_back("lepton p_{T} (GeV)");	            n.push_back(20); xi.push_back(0.);   xf.push_back(200.);
  vars.push_back("lep1.eta()");   xt.push_back("lepton #eta")       ;	            n.push_back(20); xi.push_back(-3.);  xf.push_back(3.);
  vars.push_back("jet.pt()");     xt.push_back("max jet p_{T} (GeV)");              n.push_back(20); xi.push_back(0.);   xf.push_back(400.);
  vars.push_back("jet.eta()");    xt.push_back("max jet #eta");	                    n.push_back(20); xi.push_back(-3);   xf.push_back( 3);
  vars.push_back("dphijm");       xt.push_back("#Delta#phi(max jet,pfmet)");        n.push_back(20); xi.push_back(0);    xf.push_back(3.2);
  vars.push_back("tcmet");        xt.push_back("tcmet (GeV)");    	            n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("pfmet");        xt.push_back("E_{T}^{miss} (GeV)");  	            n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("y");            xt.push_back("y (GeV^{1/2})");                    n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
  vars.push_back("htcalo");       xt.push_back("H_{T} (GeV)");		            n.push_back(20); xi.push_back(0.);   xf.push_back(1000.);
  vars.push_back("ncalojets");    xt.push_back("jet multiplicity");	            n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
  vars.push_back("nbctcm");       xt.push_back("b-jet multiplicity");	            n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
  vars.push_back("ndavtx");       xt.push_back("nDAVertices");		            n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
  vars.push_back("mt");           xt.push_back("M_{T} (GeV)");		            n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("mt");           xt.push_back("M_{T} (GeV)");		            n.push_back(20); xi.push_back(100.);   xf.push_back(300.);
  vars.push_back("mt");           xt.push_back("M_{T} (GeV)");		            n.push_back(15); xi.push_back(150.);   xf.push_back(300.);
  vars.push_back("meff");         xt.push_back("effective mass (GeV)");	            n.push_back(20); xi.push_back(0.);   xf.push_back(2000.);
  vars.push_back("trkreliso5");   xt.push_back("track rel iso (p_{T} > 5 GeV)");    n.push_back(20); xi.push_back(0.);   xf.push_back(1.);
  vars.push_back("trkreliso10");  xt.push_back("track rel iso (p_{T} > 10 GeV)");   n.push_back(20); xi.push_back(0.);   xf.push_back(1.);

  //vars.push_back("ht");         xt.push_back("H_{T} (GeV)");		              n.push_back(20); xi.push_back(0.);   xf.push_back(1000.);
  //vars.push_back("njets");      xt.push_back("jet multiplicity");	              n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
  //vars.push_back("npfjets40");  xt.push_back("jet multiplicity (p_{T} > 40 GeV)");  n.push_back(10); xi.push_back(0.);   xf.push_back(10.);

   
  const unsigned int nvars = vars.size();
  
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  TCanvas* canvas = new TCanvas("canvas","canvas",1100,750);
  gStyle->SetPaperSize(22,28);
  canvas->Print(Form("../plots/%s.ps[",filename));

  TLegend *leg = getLegend( mc , mclabels , true , 0.3 , 0.1 , 0.7 , 0.9 );
  leg->SetTextSize(0.1);
  leg->SetBorderSize(1);
  
  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    log = false;
    //if( strcmp(vars.at(ivar),"mt")==0) log = true;

    canvas->cd();

    plotpad[ivar] = new TPad("plotpad","plotpad",0,0,0.8,1);
    plotpad[ivar]->Draw();
    plotpad[ivar]->cd();
    
    plotpad[ivar]->Divide(2,2);

    plotpad[ivar]->cd(1);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==0") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "e"   );

    plotpad[ivar]->cd(2);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==1") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "m"   );

    plotpad[ivar]->cd(3);
    leg->SetTextSize(0.05);
    leg->Draw();

    plotpad[ivar]->cd(4);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel)              , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "all" );

    canvas->Print(Form("../plots/%s.ps",filename));
    canvas->Clear();

  } 

  canvas->Print(Form("../plots/%s.ps]",filename));
  canvas->Clear();
  
  gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.ps ../plots/%s.pdf",filename,filename));
}

void mt( char* path , bool printgif = false, bool latex = false ){

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
  TCut nleps1("nleps==1");
  TCut ntaus0("ntaus==0");
  TCut ntaus1("ntaus==1");
  TCut notfake("w1>0");

  TH1F* w1 = getHist( wjets , "mt" , TCut(sel+nleps1) , "hw1" , 50 , 0 , 200 );
  TH1F* t1 = getHist( ttall , "mt" , TCut(sel+nleps1) , "ht1" , 50 , 0 , 200 );

  //TH1F* w1 = getHist( wjets , "mt" , TCut(njets3) , "hw1" , 50 , 0 , 200 );
  //TH1F* t1 = getHist( ttl   , "mt" , TCut(njets3) , "ht1" , 50 , 0 , 200 );
  //TH1F* t1 = getHist( tttau   , "mt" , TCut(njets3) , "ht1" , 50 , 0 , 200 );

  TH1F* w2 = getHist( wjets , "mt" , TCut(sel+nleps1+ntaus0) , "hw2" , 50 , 0 , 200 );
  TH1F* t2 = getHist( ttall , "mt" , TCut(sel+nleps1+ntaus0) , "ht2" , 50 , 0 , 200 );

  TH1F* w3 = getHist( wjets , "mt" , TCut(sel+nleps1+ntaus1) , "hw3" , 20 , 0 , 200 );
  TH1F* t3 = getHist( ttall , "mt" , TCut(sel+nleps1+ntaus1) , "ht3" , 20 , 0 , 200 );
  

  TCanvas *c1 = new TCanvas();
  c1->cd();
  plotHist( t1 , w1 , "t#bar{t}#rightarrow l+jets" , "W+jets" , "M_{T} (GeV)" , true );  
  cout << "1-lepton        (tt,W) " << t1->GetEntries() << " " << w1->GetEntries() << endl; 

  TCanvas *c2 = new TCanvas();
  c2->cd();
  plotHist( t2 , w2 , "t#bar{t}#rightarrow l+jets" , "W+jets" , "M_{T} (GeV)" , true );  
  cout << "1-lepton, 0-tau (tt,W) " << t2->GetEntries() << " " << w2->GetEntries() << endl; 

  TCanvas *c3 = new TCanvas();
  c3->cd();
  plotHist( t3 , w3 , "t#bar{t}#rightarrow l+jets" , "W+jets" , "M_{T} (GeV)" , true );  
  cout << "1-lepton, 1-tau (tt,W) " << t3->GetEntries() << " " << w3->GetEntries() << endl; 
 
  if( printgif ){
    cout << "Printing plots" << endl;
    c1->Print("../plots/mt.png");
    c2->Print("../plots/mt_ntaus0.png");
    c3->Print("../plots/mt_ntaus1.png");
  }
}
