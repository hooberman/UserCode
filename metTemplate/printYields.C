#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iomanip>

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

//-------------------------------------------------------
float scale               = 1;//30./11.06;
const int   nprec         = 2;            //ndigits precision
const int   width1        = 20;           //width 1
const int   width2        = 4;            //width 2
const bool  printdata     = true;         //print data yields?
const bool  showError     = false;         //show errors for MC?
const bool  noLumiScale   = false;        //for MC, show event counts
const bool  makeLatexPlot = false;         //plot in latex style
//const char* yield       = "yield";              //hist name
//const char* yield       = "yieldsig";              //hist name
//const char* yield       = "yield_g2j";              //hist name
//const char* yield       = "yield_g2j";      //hist name
const char* yield       = "yield_pfmet60";      //hist name
const char* iter        = "v5"; //iteration
//const char* iter        = "oct15th_v4_generalLepVeto"; //iteration
//-------------------------------------------------------

char* pm         = " +/- ";
char* delim      = "|";
char* delimstart = "|";
char* delimend   = "|";
char* ee         = "ee";
char* mm         = "mm";
char* em         = "em";

//-------------------------------------------------------

void print( TH1F* h , string label , bool error );
void printLine();

void printYields(){

  if( makeLatexPlot ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }

  vector<char*> mcsamples;
  //mcsamples.push_back("ttall"); 
  //mcsamples.push_back("DYee"); 
  //mcsamples.push_back("DYmm"); 
 
  //mcsamples.push_back("ttotr"); 
  mcsamples.push_back("ZJets"); 
  mcsamples.push_back("TTbar"); 
   mcsamples.push_back("WJets"); 
   mcsamples.push_back("WW"); 
   mcsamples.push_back("WZ"); 
   mcsamples.push_back("ZZ"); 
   mcsamples.push_back("tW");
  
  const unsigned int nmcsamples = mcsamples.size();

  vector<char*> susysamples;
  //susysamples.push_back("LM0");
  //susysamples.push_back("LM1");
  const unsigned int nsusysamples = susysamples.size();
 
  TH1F*  h      = new TH1F();  
  TH1F*  hmctot = new TH1F();
  TFile* mcfile[nmcsamples];

  printLine();
  
  //print header
  cout << delimstart << setw(width1) << "Sample"    << setw(width2)
       << delim      << setw(width1) << ee          << setw(width2)
       << delim      << setw(width1) << mm          << setw(width2)
       << delim      << setw(width1) << em          << setw(width2)
       << delim      << setw(width1) << "tot"       << setw(width2) 
       << delimend   << endl;
  
  printLine();

  //print SM MC samples
  for(unsigned int imcsample = 0 ; imcsample < nmcsamples ; imcsample++){

    mcfile[imcsample] = TFile::Open(Form("output/%s/babylooper_%s_PhotonJetTemplate.root",
                                         iter,mcsamples.at(imcsample)));

    h = (TH1F*) mcfile[imcsample]->Get(yield);
    if(h==0) continue;

    if( imcsample == 0 ) hmctot = (TH1F*) h->Clone();
    else                 hmctot->Add(h);
    
    print( h , mcsamples[imcsample] , showError );

  }

  printLine();
  
  //print sum of SM MC samples
  print( hmctot , "tot SM MC" , showError );

  if( printdata ){
    TFile *datafile = TFile::Open(Form("output/%s/babylooper_lepdata_skim_EGStitchedTemplate.root",iter));
    h = (TH1F*) datafile->Get(yield);
    if( h != 0 ){
      print( h , "data" , false );

      printLine();
   
    }
  }

  /*
  //print SUSY MC samples  
  if( nsusysamples > 0 ){
    
    for(unsigned int isusysample = 0 ; isusysample < nsusysamples ; isusysample++){
      
      h = (TH1F*) f->Get(Form("%s_%s",susysamples.at(isusysample),yield));
      if(h==0) continue;
      
      print( h , susysamples[isusysample] , showError );
      
    }
    
    cout << "------------------------------------------------" 
         << "------------------------------------------------" << endl;
  }
  */

}

void printLine(){

  if( makeLatexPlot ){
    cout << "\\hline" << endl;
  }
  else{
    cout << "------------------------------------------------" 
         << "------------------------------------------------" << endl;
  }
}


void print( TH1F* h , string label , bool error ){

  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  float lumiscale = 1;

  if( noLumiScale ) lumiscale = h->GetEntries() / h->Integral();

  if( label == "data" ){
    if( error ){
      see  << Form("%.0f",lumiscale*scale*h->GetBinContent(2)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(2));
      smm  << Form("%.0f",lumiscale*scale*h->GetBinContent(3)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(3));
      sem  << Form("%.0f",lumiscale*scale*h->GetBinContent(4)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(4));
      stot << Form("%.0f",lumiscale*scale*h->GetBinContent(1)) << pm << Form("%.0f",lumiscale*scale*h->GetBinError(1));
    }else{
      see  << Form("%.0f",lumiscale*scale*h->GetBinContent(2));
      smm  << Form("%.0f",lumiscale*scale*h->GetBinContent(3));
      sem  << Form("%.0f",lumiscale*scale*h->GetBinContent(4));
      stot << Form("%.0f",lumiscale*scale*h->GetBinContent(1));
    }
  }else{
    if( error ){
      see  << Form("%.2f",lumiscale*scale*h->GetBinContent(2)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(2));
      smm  << Form("%.2f",lumiscale*scale*h->GetBinContent(3)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(3));
      sem  << Form("%.2f",lumiscale*scale*h->GetBinContent(4)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(4));
      stot << Form("%.2f",lumiscale*scale*h->GetBinContent(1)) << pm << Form("%.2f",lumiscale*scale*h->GetBinError(1));
    }else{
      see  << Form("%.2f",lumiscale*scale*h->GetBinContent(2));
      smm  << Form("%.2f",lumiscale*scale*h->GetBinContent(3));
      sem  << Form("%.2f",lumiscale*scale*h->GetBinContent(4));
      stot << Form("%.2f",lumiscale*scale*h->GetBinContent(1));
    }
  }

  cout << delimstart << setw(width1) << label      << setw(width2)
       << delim      << setw(width1) << see.str()  << setw(width2)
       << delim      << setw(width1) << smm.str()  << setw(width2)
       << delim      << setw(width1) << sem.str()  << setw(width2)
       << delim      << setw(width1) << stot.str() << setw(width2)
       << delimend   << endl;
  
  
}

 
