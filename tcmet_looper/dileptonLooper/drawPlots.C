#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include "TLegendEntry.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"

using namespace std;

//variables to plot
vector<char*>  vars;
vector<int>    rebin;
vector<char*>  xtitles; 
vector<char*>  mcprefix;
vector<TFile*> mcfiles;
TFile*         datafile;

vector<char*>  vars_nojets;
vector<int>    rebin_nojets;
vector<char*>  xtitles_nojets; 

int colors[]={2,5,2,3,6,7,8,9,10,11,12};

float dataintegral;
float mcintegral;
inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void initialize(){

  //variables to plot
  //vars.push_back("hpfmet");   xtitles.push_back("pfmet (GeV)"); rebin.push_back(4);
  //vars.push_back("htcmet");   xtitles.push_back("tcmet (GeV)"); rebin.push_back(4);
  vars.push_back("htcmetNew_calo");   xtitles.push_back("tcmet (GeV)"); rebin.push_back(4);
  vars.push_back("htcmetNew_pfc");    xtitles.push_back("tcmet (GeV)"); rebin.push_back(4);

  //vars_nojets.push_back("hdilMass");   xtitles_nojets.push_back("dilepton mass (GeV)");  rebin_nojets.push_back(1);
  //vars_nojets.push_back("hnjets");     xtitles_nojets.push_back("nJets");                rebin_nojets.push_back(1);
  //vars_nojets.push_back("hjetpt");     xtitles_nojets.push_back("jet p_{T} (GeV)");      rebin_nojets.push_back(5);
  
  //MC samples to include
  mcprefix.push_back("TTBar");
  mcprefix.push_back("ZJets");   

  for( unsigned i = 0 ; i < mcprefix.size() ; ++i ){
    cout << "Opening MC file  " << Form("output/%s_e_ttbarV1align_histos.root",mcprefix.at(i)) << endl;
    mcfiles.push_back( TFile::Open(Form("output/%s_e_ttbarV1align_histos.root",mcprefix.at(i)) ) );
  }

  char* datafilename = "output/data_e_ttbarV1align_histos.root";
  cout << "Opening data file " << datafilename << endl;
  datafile = TFile::Open( datafilename );
}

void formatHist(TH1F* h,char* var);
TLegend *getLegend();
THStack* getMCStack(char* varname, char* xtitle, char* histtype, char* jettype, int rebin = 1);
TH1F* getDataHist(char* varname, char* xtitle,char* leptype,char* jettype, int rebin = 1);
void drawOverlayPlot( TH1F* h , THStack* stack );

void drawPlots( bool printgif = false){
  
  initialize();

  const unsigned int nVars = vars.size();
  assert(vars.size() == xtitles.size());
  
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
 
  TCanvas *canArray[nVars][4];
  //TCanvas *canvas;

  TLegend *leg = getLegend();
  
  char* jetlabel[]={"_0j","_1j","_geq2j","_allj"};

  //looper over njets (for variables divided in njets bins)
  for( int iJ = 0 ; iJ < 4 ; iJ++ ){

    //loop over variables
    for(unsigned int iVar = 0 ; iVar < nVars ; iVar++ ){
      
      canArray[iVar][iJ] = new TCanvas(Form("%s_can%s",vars.at(iVar),jetlabel[iJ]),
                                       Form("%s_can%s",vars.at(iVar),jetlabel[iJ]),1200,600);
      canArray[iVar][iJ]->Divide(2,1);
            
      canArray[iVar][iJ]->cd(1);
      THStack* stackee = getMCStack(vars.at(iVar),xtitles.at(iVar),"ee",jetlabel[iJ],rebin.at(iVar));
      TH1F*    hee     = getDataHist(vars.at(iVar),xtitles.at(iVar),"ee",jetlabel[iJ],rebin.at(iVar));
      drawOverlayPlot( hee , stackee );
      formatHist( hee , vars.at(iVar) );
      gPad->SetLogy(1);
      leg->Draw();
      
      canArray[iVar][iJ]->cd(2);
      THStack* stackmm = getMCStack(vars.at(iVar),xtitles.at(iVar),"mm",jetlabel[iJ],rebin.at(iVar));
      TH1F*    hmm     = getDataHist(vars.at(iVar),xtitles.at(iVar),"mm",jetlabel[iJ],rebin.at(iVar));
      drawOverlayPlot( hmm , stackmm );
      formatHist( hmm , vars.at(iVar) );
      gPad->SetLogy(1);
      leg->Draw();
      
      if(printgif) canArray[iVar][iJ]->Print(Form("plots/%s%s.gif",vars.at(iVar),jetlabel[iJ]));
      
    }
  }


  TCanvas *canArray_nojets[nVars];

  //loop over variables
  for(unsigned int iVar = 0 ; iVar < vars_nojets.size() ; iVar++ ){
    
    canArray_nojets[iVar] = new TCanvas(Form("%s_can",vars_nojets.at(iVar)),
                                        Form("%s_can",vars_nojets.at(iVar)),1200,600);
    canArray_nojets[iVar]->Divide(2,1);
  
    canArray_nojets[iVar]->cd(1);
  
    THStack* stackee = getMCStack(vars_nojets.at(iVar),xtitles_nojets.at(iVar),"ee","",rebin_nojets.at(iVar));
    TH1F*    hee     = getDataHist(vars_nojets.at(iVar),xtitles_nojets.at(iVar),"ee","",rebin_nojets.at(iVar));
    drawOverlayPlot( hee , stackee );
    formatHist( hee , vars_nojets.at(iVar) );
    gPad->SetLogy(1);
    leg->Draw();
    
    canArray_nojets[iVar]->cd(2);
    THStack* stackmm = getMCStack(vars_nojets.at(iVar),xtitles_nojets.at(iVar),"mm","",rebin_nojets.at(iVar));
    TH1F*    hmm     = getDataHist(vars_nojets.at(iVar),xtitles_nojets.at(iVar),"mm","",rebin_nojets.at(iVar));
    drawOverlayPlot( hmm , stackmm );
    formatHist( hmm , vars_nojets.at(iVar) );
    gPad->SetLogy(1);
    leg->Draw();
    
    if(printgif) canArray_nojets[iVar]->Print(Form("plots/%s.gif",vars_nojets.at(iVar)));
    
  }
    
  







}

void drawOverlayPlot( TH1F* h , THStack* stack ){
  h->Draw("E1");
  stack->Draw("same");
  h->Draw("sameE1");
  h->Draw("sameaxis");
}


TH1F* getDataHist(char* varname, char* xtitle,char* leptype,char* jettype, int rebin){

  TH1F *h  = (TH1F*)datafile->Get(Form("%s_%s%s",varname,leptype,jettype));
  dataintegral = h->Integral( h->FindBin(30) , 100000 );

  if( rebin > 1 ) h->Rebin(rebin);

  h->SetMinimum( 0.001 );
  h->SetMarkerSize(0.5);
  h->GetXaxis()->SetTitle(xtitle);
  //h->SetTitle(leptitles[iLep].c_str());



  return h;
}

void formatHist(TH1F* h,char* var){
  TLine line;
  line.SetLineColor(4);
  line.SetLineWidth(2);
  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.04);

  if(strcmp(var,"hdilMass")==0){
    line.DrawLine( 76  , h->GetMinimum() , 76  , 2 * h->GetMaximum() );
    line.DrawLine( 106 , h->GetMinimum() , 106 , 2 * h->GetMaximum() );
  }
  if(strcmp(var,"htcmet")==0 || strcmp(var,"htcmetNew_calo")==0 || strcmp(var,"htcmetNew_pfc")==0){
    line.DrawLine( 30  , h->GetMinimum() , 30  , 2 * h->GetMaximum() );
    stringstream s1;
    s1 << "N(#slash{E_{T}} > 30 GeV): " << dataintegral ;
    stringstream s2;
    s2 << "predicted:       " << fround(mcintegral,2) ;
    t.DrawLatex( 0.45 , 0.8 ,  s1.str().c_str() );
    t.DrawLatex( 0.45 , 0.75 , s2.str().c_str() );
  }
  if(strcmp(var,"hpfmet")==0){
    line.DrawLine( 30  , h->GetMinimum() , 30  , 2 * h->GetMaximum() );
    stringstream s1;
    s1 << "N(#slash{E_{T}} > 30 GeV): " << dataintegral ;
    stringstream s2;
    s2 << "predicted:       " << fround(mcintegral,2) ;
    t.DrawLatex( 0.45 , 0.8 ,  s1.str().c_str() );
    t.DrawLatex( 0.45 , 0.75 , s2.str().c_str() );
  }
  if(strcmp(var,"htcmetNew_calo")==0){
    line.DrawLine( 30  , h->GetMinimum() , 30  , 2 * h->GetMaximum() );
  }
  if(strcmp(var,"htcmetNew_pfc")==0){
    line.DrawLine( 30  , h->GetMinimum() , 30  , 2 * h->GetMaximum() );
  }
}



TLegend *getLegend(){

  const unsigned int nMC = mcprefix.size();

  TLegend *leg =new TLegend(0.8,0.8,0.95,0.95);
  
  TH1F* hdummy[nMC];

  for(unsigned int i = 0 ; i < nMC ; i++){
    hdummy[i]=new TH1F(Form("hdummy_%i",i),"",1,0,1);
    hdummy[i]->SetLineColor(1);
    hdummy[i]->SetMarkerColor(colors[i]);
    hdummy[i]->SetFillColor(colors[i]);
    leg->AddEntry(hdummy[i],mcprefix.at(i),"f");
    //delete hdummy;
  }

  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  return leg;
}


THStack* getMCStack(char* varname, char* xtitle,char* leptype,char* jettype, int rebin){
  
  const unsigned int nMC = mcprefix.size();
  
  THStack *stack=new THStack("stack","");
  TH1F* h[nMC];

  TString varstring(varname);
  cout<<"Getting "<<Form("%s_%s%s",varname,leptype,jettype)<<endl;
 
  mcintegral = 0;

  for(unsigned int iMC=0;iMC<nMC;iMC++){
    
    h[iMC]  = (TH1F*)mcfiles.at(iMC)->Get(Form("%s_%s%s",varname,leptype,jettype));

    mcintegral += h[iMC]->Integral( h[iMC]->FindBin(30) , 1000000 );

    if( rebin > 1 ) h[iMC]->Rebin(rebin);
    h[iMC]->SetLineColor(1);
    h[iMC]->SetMarkerColor(colors[iMC]);
    h[iMC]->SetFillColor(colors[iMC]);

    //if(rebin(varname)>1) h[iMC]->Rebin(rebin(varname));
    stack->Add(h[iMC]);
  }
  
  cout << "integral > 30 " << mcintegral << endl;

  if(strcmp(leptype,"all")==0) stack->SetTitle(Form("%s (all leptons)",xtitle));
  if(strcmp(leptype,"ee")==0)  stack->SetTitle(Form("%s (ee)",xtitle));
  if(strcmp(leptype,"mm")==0)  stack->SetTitle(Form("%s (#mu#mu)",xtitle));
  if(strcmp(leptype,"em")==0)  stack->SetTitle(Form("%s (e#mu)",xtitle));

  return stack;
}

