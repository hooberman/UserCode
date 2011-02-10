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
#include "TCut.h"
#include "TRandom3.h"
#include <iomanip>

using namespace std;

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

//-------------------------------------------------------
float scale               = 1;//30./11.06;
const int   nprec         = 2;            //ndigits precision
const int   width1        = 15;           //width 1
const int   width2        = 4;            //width 2
const bool  printdata     = true;         //print data yields?
const bool  showError     = false;        //show errors for MC?
const bool  noLumiScale   = false;        //for MC, show event counts
const bool  makeLatexPlot = false;        //plot in latex style
const char* iter          = "v2";         //iteration

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

  enum selectionEnum { WW = 0, H130 = 1 , H160 = 2 , H200 = 3 };
  
  //selectionEnum selType = WW;
  //selectionEnum selType = H130;
  //selectionEnum selType = H160;
  selectionEnum selType = H200;


  //-------------------------------------------------------
  //selection
  //-------------------------------------------------------
  
  TCut met_projpt = "(event_type != 2 && met_projpt>35. ) || (event_type == 2 && met_projpt > 20.)";
  TCut pt2020     = "lephard_pt > 20 && lepsoft_pt > 20";
  TCut pt2010     = "lephard_pt > 20 && lepsoft_pt > 10";
  TCut jetveto    = "jets_num==0 && extralep_num==0 && lowptbtags_num==0 && softmu_num==0";
  TCut eetype     = "event_type==3";
  TCut mmtype     = "event_type==0";
  TCut emtype     = "event_type==2";
  TCut mll12      = "dil_mass > 12.";
  TCut h130sel    = "dil_dphi < 1.05 && dil_mass < 45. && lephard_pt > 25 && lepsoft_pt > 20";
  TCut h160sel    = "dil_dphi < 1.05 && dil_mass < 50. && lephard_pt > 30 && lepsoft_pt > 25";
  TCut h200sel    = "dil_dphi < 1.75 && dil_mass < 90. && lephard_pt > 40 && lepsoft_pt > 25";
  TCut weight     = "event_scale1fb * 0.0355";
  TCut wwsel      = pt2020 + met_projpt + jetveto + mll12;
  TCut sel;

  if( selType == WW ){
    cout << "Using WW selection" << endl;
    sel = wwsel; 
  }
  else if( selType == H130 ){
    cout << "Using Higgs 130 GeV selection" << endl;
    sel = wwsel + h130sel; 
  }
  else if( selType == H160 ){
    cout << "Using Higgs 160 GeV selection" << endl;
    sel = wwsel + h160sel; 
  }
  else if( selType == H200 ){
    cout << "Using Higgs 200 GeV selection" << endl;
    sel = wwsel + h200sel; 
  }
  else{
    cout << "Error, unrecognized selection " << selType << ", quitting" << endl;
    exit(0);
  }

//   TCut eeselweight  = (sel+eetype)*weight;
//   TCut mmselweight  = (sel+mmtype)*weight;
//   TCut emselweight  = (sel+emtype)*weight;
//   TCut allselweight =  sel        *weight;

  //-------------------------------------------------------

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
  mcsamples.push_back("WWTo2L2Nu");
  mcsamples.push_back("GluGluToWWTo4L");
  mcsamples.push_back("WZ");
  mcsamples.push_back("ZZ");
  mcsamples.push_back("TTJets");
  mcsamples.push_back("tW");
  mcsamples.push_back("WJetsToLNu");
  mcsamples.push_back("DY");
  if( selType == H130 ) mcsamples.push_back("Higgs130");
  if( selType == H160 ) mcsamples.push_back("Higgs160");
  if( selType == H200 ) mcsamples.push_back("Higgs200");
  
  const unsigned int nmcsamples = mcsamples.size();

  vector<char*> sigsamples;
  //sigsamples.push_back("130");
  const unsigned int nsigsamples = sigsamples.size();
 
  //TH1F*  h      = new TH1F();  
  TH1F*  hmctot = new TH1F();
 
  TChain* ch[nmcsamples];
  TH1F*   h[nmcsamples];

  //print SM MC samples
  for(unsigned int imc = 0 ; imc < nmcsamples ; imc++){

    ch[imc] = new TChain("Events");
    if( strcmp( mcsamples.at(imc) , "DY" ) == 0 ){
      ch[imc] -> Add( Form("babies/%s/DYToMuMuM20_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/DYToMuMuM10To20_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/DYToEEM20_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/DYToEEM10To20_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/DYToTauTauM20_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/DYToTauTauM10To20_PU_testFinal_baby.root",iter) );
    }
    else if( strcmp( mcsamples.at(imc) , "Higgs130" ) == 0 ){
      ch[imc] -> Add( Form("babies/%s/HToWWTo2L2NuM130_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/HToWWToLNuTauNuM130_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/HToWWTo2Tau2NuM130_PU_testFinal_baby.root",iter) );
    }
    else if( strcmp( mcsamples.at(imc) , "Higgs160" ) == 0 ){
      ch[imc] -> Add( Form("babies/%s/HToWWTo2L2NuM160_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/HToWWToLNuTauNuM160_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/HToWWTo2Tau2NuM160_PU_testFinal_baby.root",iter) );
    }
    else if( strcmp( mcsamples.at(imc) , "Higgs200" ) == 0 ){
      ch[imc] -> Add( Form("babies/%s/HToWWTo2L2NuM200_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/HToWWToLNuTauNuM200_PU_testFinal_baby.root",iter) );
      ch[imc] -> Add( Form("babies/%s/HToWWTo2Tau2NuM200_PU_testFinal_baby.root",iter) );
    }
    else{
      ch[imc] -> Add( Form("babies/%s/%s_PU_testFinal_baby.root",iter,mcsamples.at(imc)) );
    }

    h[imc] = new TH1F(Form("h_%i",imc),"",4,0,4);

    TH1F* hee  = new TH1F("hee", "",1,0,1);
    TH1F* hmm  = new TH1F("hmm", "",1,0,1);
    TH1F* hem  = new TH1F("hem", "",1,0,1);
    TH1F* hall = new TH1F("hall","",1,0,1);

    hee->Sumw2();
    hmm->Sumw2();
    hem->Sumw2();
    hall->Sumw2();
    
    ch[imc]->Draw("0.5>>hee", (sel+eetype)*weight);
    ch[imc]->Draw("0.5>>hmm", (sel+mmtype)*weight);
    ch[imc]->Draw("0.5>>hem", (sel+emtype)*weight);
    ch[imc]->Draw("0.5>>hall", sel    *weight);

    h[imc]->SetBinContent ( 2 , hee->GetBinContent(1)  );
    h[imc]->SetBinError   ( 2 , hee->GetBinError(1)    );
    h[imc]->SetBinContent ( 3 , hmm->GetBinContent(1)  );
    h[imc]->SetBinError   ( 3 , hmm->GetBinError(1)    );
    h[imc]->SetBinContent ( 4 , hem->GetBinContent(1)  );
    h[imc]->SetBinError   ( 4 , hem->GetBinError(1)    );
    h[imc]->SetBinContent ( 1 , hall->GetBinContent(1) );
    h[imc]->SetBinError   ( 1 , hall->GetBinError(1)   );

    if( !TString( mcsamples.at(imc) ).Contains("Higgs") ) {
      if( imc == 0 ) hmctot = (TH1F*) h[imc]->Clone();
      else           hmctot->Add(h[imc]);
    }

    delete hee;
    delete hmm;
    delete hem;
    delete hall;
  }


  printLine();
  
  //print header
  cout << delimstart << setw(width1) << "Sample"    << setw(width2)
       << delim      << setw(width1) << ee          << setw(width2)
       << delim      << setw(width1) << mm          << setw(width2)
       << delim      << setw(width1) << em          << setw(width2)
       << delim      << setw(width1) << "tot"       << setw(width2) 
       << delimend   << endl;
  
  printLine();

  for(unsigned int imc = 0 ; imc < nmcsamples ; imc++){
    if( TString( mcsamples.at(imc) ).Contains("Higgs") ) continue;
    print( h[imc] , mcsamples[imc] , showError );
  }

  printLine();
  
  //print sum of SM MC samples
  print( hmctot , "tot SM MC" , showError );

  printLine();

  for(unsigned int imc = 0 ; imc < nmcsamples ; imc++){
    if( !TString( mcsamples.at(imc) ).Contains("Higgs") ) continue;
    print( h[imc] , mcsamples[imc] , showError );
  }

  printLine();

//   if( printdata ){
//     TFile *datafile = TFile::Open(Form("output/%s/babylooper_lepdata_skim_EGStitchedTemplate.root",iter));
//     h = (TH1F*) datafile->Get(yield);
//     if( h != 0 ){
//       print( h , "data" , false );

//       printLine();
   
//     }
//   }

  /*
  //print SIG MC samples  
  if( nsigsamples > 0 ){
    
    for(unsigned int isigsample = 0 ; isigsample < nsigsamples ; isigsample++){
      
      h = (TH1F*) f->Get(Form("%s_%s",sigsamples.at(isigsample),yield));
      if(h==0) continue;
      
      print( h , sigsamples[isigsample] , showError );
      
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
    cout << "-------------------------------------------------------------" 
         << "------------------------------------------------------------" << endl;
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
