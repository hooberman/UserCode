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
#include "Hootilities.C"
#include "mmjj.h"
//#include "histtools.h"

using namespace std;

bool printgif_           = false;
bool alreadyInitialized_ = false;
const float lumi = 0.8;
 

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
    dy->Reset();

    mc.clear();
    mctex.clear();
    mclabels.clear();
  }

  else{

    data	= new TChain("t");
    ttall	= new TChain("t");
    dy  	= new TChain("t");
  }

  cout << endl;
  cout << "Loading babies at       : " << path << endl;
  

  data->Add(Form("%s/datamay10_smallTree.root",path));
  data->Add(Form("%s/dataPRv4_smallTree.root",path));
  data->Add(Form("%s/dataaug05_smallTree.root",path));
  data->Add(Form("%s/dataPRv6_smallTree.root",path));
  data->Add(Form("%s/data2011B_smallTree.root",path));

  ttall->Add(Form("%s/ttall_smallTree.root",path));
  dy->Add(Form("%s/Zjets_2jets_smallTree.root",path));
  

  /*
  data->Add("../output/V00-00-07/dataping_smallTree.root");
  ttall->Add("../output/V00-00-08/ttall_smallTree.root");
  dy->Add("../output/V00-00-08//Zjets_2jets_smallTree.root");
  */

  /*
  data->Add(Form("%s/data165_smallTree.root",path));
  data->Add(Form("%s/data166_smallTree.root",path));
  data->Add(Form("%s/data167_smallTree.root",path));
  data->Add(Form("%s/data168_smallTree.root",path));
  data->Add(Form("%s/dataPRv6_smallTree.root",path));

  ttall->Add(Form("%s/ttall_smallTree.root",path));
  dy->Add(Form("%s/DYtot_smallTree_tenPercent.root",path));
  */

  
  mc.push_back(ttall);       mclabels.push_back("ttall");    
  mc.push_back(dy);          mclabels.push_back("dy");    
  
  alreadyInitialized_ = true;
}

//------------------------------------------
// selection and weight to apply to babies
//------------------------------------------

TCut selection_TCut(){

  TCut njets2("njets>=2");
  TCut leadjet120("jet1.pt()>120");
  TCut mll50("dilmass>50");
  TCut zveto("dilmass<70 || dilmass>100");
  TCut zmass("dilmass>81&&dilmass<101");
  TCut pt2520("lep1.pt()>25 && lep2.pt()>20");
  TCut ngoodlep2("ngoodlep>=2");
  TCut pt3020("lep1.pt()>30 && lep2.pt()>20");
  TCut nbtags1("nbtags20>=1");
  //TCut nbtags1("nbtags33>=1");
  TCut iso1("iso1<0.15");
  TCut iso2("iso2<0.15");
  TCut id1("passid1==1");
  TCut id2("passid2==1");
  TCut mu2("ngoodmu>=2");
  TCut bump("mmjj>850 && mmjj<1150");
  TCut left("mmjj>650 && mmjj<850");
  TCut right("mmjj>1150");

  TCut sel;
  sel += pt3020;
  sel += njets2;
  //sel += zmass;
  sel += leadjet120;
  sel += mll50;
  sel += zveto;
  sel += ngoodlep2;
  // sel += "mmjj > 600 && mmjj < 1400";
  // sel += "mmjj > 850 && mmjj < 1150";
  sel += nbtags1;
  //sel += "dilmass>100";
  //sel += iso1;
  //sel += iso2;
  // sel += id1;
  // sel += id2;
  // sel += mu2;
  //sel += left;
  //sel += right;
  //sel += bump;
  //sel += "pfmet<20";


  cout << "Using selection         : " << sel.GetTitle() << endl;
 
  return sel;
}

TCut weight_TCut(){

  //TCut weight("weight*ndavtxweight");
  //TCut weight("weight * 0.204 * ndavtxweight");
  //TCut weight("weight * 0.8 * ndavtxweight");
  TCut weight("weight * ndavtxweight * 4.7");
  //TCut weight("weight * ndavtxweight * 2.0");
  //TCut weight("1");

  cout << "Using weight            : " << weight.GetTitle() << endl;
  return weight; 
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

  //vars.push_back("mmjj");             xt.push_back("m(#mu#mujj) (GeV)");           n.push_back(50); xi.push_back(350.); xf.push_back(1850.);
  //vars.push_back("mmjjdef");            xt.push_back("m(#mu#mujj) (GeV)");           n.push_back(25); xi.push_back(350.); xf.push_back(1850.);
  //vars.push_back("ndavtx");             xt.push_back("ndavtx");                      n.push_back(20); xi.push_back(0.); xf.push_back(20.);
  //vars.push_back("dilmass");          xt.push_back("m(#mu#mu) (GeV)");             n.push_back(50); xi.push_back(0.); xf.push_back(200.);
  vars.push_back("dilmass");          xt.push_back("m(#mu#mu) (GeV)");             n.push_back(30); xi.push_back(0.); xf.push_back(300.);
  //vars.push_back("lep1.pt()");        xt.push_back("primary muon p_{T} (GeV)");    n.push_back(10); xi.push_back(0.); xf.push_back(300.);
  //vars.push_back("lep2.pt()");        xt.push_back("secondary muon p_{T} (GeV)");  n.push_back(10); xi.push_back(0.); xf.push_back(200.);
  //vars.push_back("abs(jet1.eta())");    xt.push_back("1^{st} jet |#eta|");      n.push_back(6); xi.push_back(0); xf.push_back(3.);
  //vars.push_back("abs(jet2.eta())");    xt.push_back("2^{nd} jet |#eta|");      n.push_back(6); xi.push_back(0); xf.push_back(3.);

  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  bool residual = true;

  int canCounter = -1;
  bool log = true;

  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    TString tvar(vars[ivar]);
    tvar.ReplaceAll("()","");
    tvar.ReplaceAll(".","");
    const char* myvar = tvar;

    //if( ivar < 2 ) log = true;
    //else           log = false;

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
      can[ivar] = new TCanvas(Form("%s_can",myvar),Form("%s_can",myvar),600,600);
    }

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , !combine4 , log );

    //if( printgif ) can[ivar]->Print(Form("../plots/%s.pdf",vars[ivar]));
    if( printgif ){

      TString tvar(myvar);
      tvar.ReplaceAll("(","");
      tvar.ReplaceAll(")","");
      tvar.ReplaceAll(".","");
      const char* mtvar = tvar;

      can[ivar]->Print(Form("../plots/%s.eps",mtvar));
      can[ivar]->Print(Form("../plots/%s.gif",mtvar));
      can[ivar]->Print(Form("../plots/%s.png",mtvar));
      gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.eps ../plots/%s.pdf",mtvar,mtvar));
    }
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
  vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");        n.push_back(20); xi.push_back(0.);   xf.push_back(400.);
  vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	            n.push_back(20); xi.push_back(-3);   xf.push_back( 3);
  vars.push_back("dphijm");     xt.push_back("#Delta#phi(max jet,pfmet)");  n.push_back(20); xi.push_back(0);    xf.push_back(3.2);
  vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");    	    n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("pfmet");      xt.push_back("E_{T}^{miss} (GeV)");	    n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("y");          xt.push_back("y (GeV^{1/2})");  	    n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
  vars.push_back("ht");         xt.push_back("H_{T} (GeV)");		    n.push_back(20); xi.push_back(0.);   xf.push_back(1000.);
  vars.push_back("njets");      xt.push_back("jet multiplicity");	    n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
  vars.push_back("npfjets40");     xt.push_back("jet multiplicity (p_{T} > 40 GeV)");	    n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
  vars.push_back("nbtagstcm");  xt.push_back("b-jet multiplicity");	    n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
  vars.push_back("ndavtx");     xt.push_back("nDAVertices");		    n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
  vars.push_back("mt");         xt.push_back("M_{T} (GeV)");		    n.push_back(30); xi.push_back(0.);   xf.push_back(300.);
  vars.push_back("meff");       xt.push_back("effective mass (GeV)");	    n.push_back(20); xi.push_back(0.);   xf.push_back(2000.);
   
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
    if( strcmp(vars.at(ivar),"mt")==0) log = true;

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
