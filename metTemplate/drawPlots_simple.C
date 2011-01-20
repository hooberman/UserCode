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

//char* iter = "e_ttbarV1align_aug30"; //1.9/pb
char* iter = "v5";

using namespace std;

struct histStruct{
  TH1F*    hist;
  THStack* stack;
};

//variables to plot
vector<char*>  vars;
vector<int>    rebin_;
vector<char*>  xtitles; 
vector<char*>  mcprefix;
vector<TFile*> mcfiles;
TFile*         datafile;
vector<char*>  mcfilenames;

int colors[]={2,5,3,6,7,8,9,10,11,12};

float dataintegral;
float mcintegral;
inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void initialize(){

  //variables to plot
  vars.push_back("hpfmet_all_allj");                     xtitles.push_back("pfmet (GeV)");                 rebin_.push_back(4);
  vars.push_back("hdilmass_all_allj");          xtitles.push_back("dilepton mass (GeV)");         rebin_.push_back(2);
  vars.push_back("hdilmass_ee_allj");          xtitles.push_back("dilepton mass (GeV)");         rebin_.push_back(2);
  vars.push_back("hdilmass_mm_allj");          xtitles.push_back("dilepton mass (GeV)");         rebin_.push_back(2);
  char* datafilename = Form("output/%s/babylooper_lepdata_skim_EGStitchedTemplate.root",iter);
  //mcprefix.push_back("W+jets");      mcfilenames.push_back(Form("output/%s/babylooper_WJets_PhotonJetTemplate.root",iter));
  //mcprefix.push_back("single top");  mcfilenames.push_back(Form("output/%s/babylooper_tW_PhotonJetTemplate.root",iter));
  //mcprefix.push_back("WW");          mcfilenames.push_back(Form("output/%s/babylooper_WW_PhotonJetTemplate.root",iter));
  //mcprefix.push_back("WW");          mcfilenames.push_back(Form("output/%s/babylooper_WW_PhotonJetTemplate.root",iter));
  //mcprefix.push_back("WZ");          mcfilenames.push_back(Form("output/%s/babylooper_WZ_PhotonJetTemplate.root",iter));
  //mcprefix.push_back("VV");          mcfilenames.push_back(Form("output/%s/babylooper_WW_PhotonJetTemplate.root",iter));
  mcprefix.push_back("t#bar{t}");    mcfilenames.push_back(Form("output/%s/babylooper_TTbar_PhotonJetTemplate.root",iter));
  mcprefix.push_back("Z+jets");      mcfilenames.push_back(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));
 

  for( unsigned i = 0 ; i < mcprefix.size() ; ++i ){
    cout << "Opening MC file  " << mcfilenames.at(i) << endl;
    mcfiles.push_back( TFile::Open( mcfilenames.at(i) ) );
  }

  cout << "Opening data file " << datafilename << endl;
  datafile = TFile::Open( datafilename );
}

void formatHist(TH1F* h,char* var);
TLegend *getLegend();
histStruct getMCStack(char* varname, char* xtitle, int rebin = 1, float ndata = -1);
TH1F* getDataHist(char* varname, char* xtitle, int rebin = 1);
void drawOverlayPlot( TH1F* h , THStack* stack );

void drawPlots( bool printgif = false){
  
  initialize();

  const unsigned int nVars = vars.size();
  assert(vars.size() == xtitles.size());
  
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
 
  TLegend *leg = getLegend();

  TCanvas *canArray[nVars];

  //loop over variables
  for(unsigned int iVar = 0 ; iVar < vars.size() ; iVar++ ){
    
    canArray[iVar] = new TCanvas(Form("%s_can",vars.at(iVar)), Form("%s_can",vars.at(iVar)),800,600);
    canArray[iVar]->cd();

    TH1F*      hall          = getDataHist( vars.at(iVar) , xtitles.at(iVar) , rebin_.at(iVar));  
    histStruct structall     = getMCStack(  vars.at(iVar) , xtitles.at(iVar) , rebin_.at(iVar));
    TH1F*      mcsumhistall  = structall.hist;
    THStack*   stackall      = structall.stack;

    drawOverlayPlot( hall , stackall );
    formatHist( hall , vars.at(iVar) );
    float KSall = hall->KolmogorovTest( mcsumhistall );
    stringstream sall;
    //sall << hall->GetTitle() << " (KS = " << fround(KSall,2) << ")";
    //sall << "Selected Dilepton + #geq2Jet  Events" << " (KS = " << fround(KSall,2) << ")";
    hall->SetTitle( "Selected Dilepton + #geq2Jet  Events" );
    gPad->SetLogy(1);
    leg->Draw();
    
    if(printgif) canArray[iVar]->Print(Form("plots/%s.gif",vars.at(iVar)));
    
  }
}

void drawOverlayPlot( TH1F* h , THStack* stack ){
  h->Draw("E1");
  stack->Draw("same");
  h->Draw("sameE1");
  h->Draw("sameaxis");
}


TH1F* getDataHist(char* varname, char* xtitle, int rebin){

  TString varstring(varname);

  TH1F *h  = (TH1F*)datafile->Get( varname );

  if( varstring.Contains("met") ){
    dataintegral = h->Integral( h->FindBin(30) , 100000 );
  }
  else if( varstring.Contains("dilmass") ){
    dataintegral = h->Integral( h->FindBin(76) , h->FindBin(106) - 1 );
  }

  if( rebin > 1 ) h->Rebin(rebin);

  h->SetMinimum( 0.1 );
  //h->SetMaximum( 1000 );
  h->SetMarkerSize(0.8);
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
  TString varstring(var);

  if( varstring.Contains("dilmass") ){
    line.DrawLine( 76  , h->GetMinimum() , 76  , 2 * h->GetMaximum() );
    line.DrawLine( 106 , h->GetMinimum() , 106 , 2 * h->GetMaximum() );

    stringstream s1;
    s1 << "N(Z): " << dataintegral ;
    stringstream s2;
    s2 << "predicted:    " << fround(mcintegral,2) ;
    cout << s1.str() << endl;
    cout << s2.str() << endl;
    //t.DrawLatex( 0.45 , 0.8 ,  s1.str().c_str() );
    //t.DrawLatex( 0.45 , 0.75 , s2.str().c_str() );
  }
  /*
  if( varstring.Contains("met") ){
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
    }*/
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


//THStack* getMCStack(char* varname, char* xtitle,char* leptype,char* jettype, int rebin){
histStruct getMCStack(char* varname, char* xtitle, int rebin, float ndata){

  const unsigned int nMC = mcprefix.size();
  
  THStack *stack=new THStack("stack","");
  TH1F* h[nMC];
  TH1F* sumhist = new TH1F();

  TString varstring(varname);
  //cout<<"Getting "<<Form("%s_%s%s",varname,leptype,jettype)<<endl;
 
  mcintegral = 0;
  
  for(unsigned int iMC=0;iMC<nMC;iMC++){
    
    h[iMC]  = (TH1F*)mcfiles.at(iMC)->Get(varname);

    if( varstring.Contains("met") ){
      mcintegral += h[iMC]->Integral( h[iMC]->FindBin(30) , 1000000 );
    }
    else if( varstring.Contains("dilmass") ){
      mcintegral += h[iMC]->Integral( h[iMC]->FindBin(76) , h[iMC]->FindBin(106) - 1 );
    }

    if( rebin > 1 ) h[iMC]->Rebin(rebin);
    h[iMC]->SetLineColor(1);
    h[iMC]->SetMarkerColor(colors[iMC]);
    h[iMC]->SetFillColor(colors[iMC]);
    
    if( iMC == 0 )    sumhist = (TH1F*) h[iMC]->Clone(Form("%s_MC_clone",h[iMC]->GetName()));
    else              sumhist->Add( h[iMC] );
    
  }

  for(unsigned int iMC=0;iMC<nMC;iMC++){
    if( ndata > 0 ){
      h[iMC]->Scale( ndata / sumhist->Integral() );
    }
    stack->Add(h[iMC]);
  }

  if( ndata > 0 ){
    mcintegral *= ndata/sumhist->Integral();
  }

  //cout << "integral > 30 " << mcintegral << endl;

//   if(strcmp(leptype,"all")==0) stack->SetTitle(Form("%s (all leptons)",xtitle));
//   if(strcmp(leptype,"ee")==0)  stack->SetTitle(Form("%s (ee)",xtitle));
//   if(strcmp(leptype,"mm")==0)  stack->SetTitle(Form("%s (#mu#mu)",xtitle));
//   if(strcmp(leptype,"em")==0)  stack->SetTitle(Form("%s (e#mu)",xtitle));

  histStruct myStruct;
  myStruct.hist  = sumhist;
  myStruct.stack = stack;

  //return stack;
  return myStruct;
}


























//   if(strcmp(var,"htcmet")==0 || strcmp(var,"htcmetNew_calo")==0 || strcmp(var,"htcmetNew_pfc")==0){
//     line.DrawLine( 30  , h->GetMinimum() , 30  , 2 * h->GetMaximum() );
//     stringstream s1;
//     s1 << "N(#slash{E_{T}} > 30 GeV): " << dataintegral ;
//     stringstream s2;
//     s2 << "predicted:       " << fround(mcintegral,2) ;
//     t.DrawLatex( 0.45 , 0.8 ,  s1.str().c_str() );
//     t.DrawLatex( 0.45 , 0.75 , s2.str().c_str() );
//   }
//   if(strcmp(var,"hpfmet")==0){
//     line.DrawLine( 30  , h->GetMinimum() , 30  , 2 * h->GetMaximum() );
//     stringstream s1;
//     s1 << "N(#slash{E_{T}} > 30 GeV): " << dataintegral ;
//     stringstream s2;
//     s2 << "predicted:       " << fround(mcintegral,2) ;
//     t.DrawLatex( 0.45 , 0.8 ,  s1.str().c_str() );
//     t.DrawLatex( 0.45 , 0.75 , s2.str().c_str() );
//   }
