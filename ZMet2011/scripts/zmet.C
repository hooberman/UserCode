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
#include "TMath.h"
#include <sstream>
#include "Hootilities.h"
#include "zmet.h"
//#include "histtools.h"

using namespace std;

bool printgif_           = false;
bool alreadyInitialized_ = false;

//--------------------------------------------------
// initialize data/MC samples
//--------------------------------------------------

void initialize(char* path){

  if( alreadyInitialized_ ){

    cout << "Resetting babies" << endl;

    data->Reset();
    tt->Reset();
    zjets->Reset();

    mc.clear();
    mclabels.clear();
  }

  else{

    data	= new TChain("T1");
    tt	        = new TChain("T1");
    zjets	= new TChain("T1");
  }

  cout << endl;
  cout << "Loading babies at       : " << path << endl;
  
  data->Add(Form("%s/data_baby.root",path));
  tt->Add(Form("%s/ttbar_baby.root",path));
  zjets->Add(Form("%s/zjets_baby.root",path));

  mc.push_back(tt);     mclabels.push_back("tt");
  mc.push_back(zjets);  mclabels.push_back("zjets");

  alreadyInitialized_ = true;
}

//------------------------------------------
// selection and weight to apply to babies
//------------------------------------------

TCut selection_TCut(){

  TCut njets2("njets >= 2");
  TCut njets3("njets >= 3");
  TCut njets4("njets >= 4");
  TCut zpass("dilmass>81.0 && dilmass<101.0");
  TCut met60("pfmet > 60.0");
  TCut met120("pfmet > 120.0");
  TCut ht200("sumjetpt > 200.0");

  TCut sel = zpass + njets2 + ht200 + met60;

  cout << "Using selection         : " << sel.GetTitle() << endl;

  return sel;
}

TCut weight_TCut(){

  //TCut weight("weight");
  TCut weight("1");
  
  cout << "Using weight            : " << weight.GetTitle() << endl;
  return weight;
}


void plotHist( TH1F* h1 , TH1F* h2 , char* leg1 , char* leg2 , char* xtitle ){

  gPad->SetLogy(1);

  h1->GetXaxis()->SetTitle( xtitle );
  h1->DrawNormalized("hist");
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->DrawNormalized("sameE1");

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h1,leg1,"l");
  leg->AddEntry(h2,leg2,"p");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

}

//------------------------------------
// print yield table
//------------------------------------

void printYieldTable( char* path ){

  deleteHistos();

  initialize(path);
  
  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  printYields( mc , mclabels , data , sel , weight , false );
 
}


//--------------------------------------------------
// make data/MC plots
//--------------------------------------------------

void makePlots( char* path , bool printgif = false ){

  initialize(path);

  TCut sel    = selection_TCut();
  TCut weight = weight_TCut();

  //sel = sel + "y > 8.5 && ht > 300";

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  //vars.push_back("tcmet");     xt.push_back("tcmet (GeV)");      n.push_back(20); xi.push_back(0.); xf.push_back(200.);
  //vars.push_back("ht");        xt.push_back("H_{T} (GeV)");      n.push_back(20); xi.push_back(0.); xf.push_back(1000.);
  //vars.push_back("dilmass");   xt.push_back("M(ll) (GeV)");      n.push_back(30); xi.push_back(0.); xf.push_back(300.);
  //vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(30); xi.push_back(0.); xf.push_back(300.);
  //vars.push_back("njets");     xt.push_back("njets");            n.push_back(6);  xi.push_back(0);  xf.push_back(6);
  vars.push_back("dilep.mass()");     xt.push_back("dilmass (GeV)");    n.push_back(100);  xi.push_back(0);  xf.push_back(200);
  
  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
    
  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     
    can[ivar] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),600,600);

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] );

    if( printgif ) can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
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




