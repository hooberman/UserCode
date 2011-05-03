#include "Hootilities.h"
#include <algorithm>
#include <iostream>
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
#include <iomanip>


float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

void printLine( bool latex ){

  if( latex ){
    cout << "\\hline" << endl;
  }
  else{
    for( int i = 0 ; i < linelength ; ++i ) cout << "-";
    cout << endl;
  }
}

void printHeader(){

  cout << delimstart << setw(width1) << "Sample"    << setw(width2)
       << delim      << setw(width1) << ee          << setw(width2)
       << delim      << setw(width1) << mm          << setw(width2)
       << delim      << setw(width1) << em          << setw(width2)
       << delim      << setw(width1) << "tot"       << setw(width2) 
       << delimend   << endl;

}


void print( TH1F* h , string label ){

  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  if( label == "data" ){
    see  << Form( "%.0f" , h->GetBinContent(1) );
    smm  << Form( "%.0f" , h->GetBinContent(2) );
    sem  << Form( "%.0f" , h->GetBinContent(3) );
    stot << Form( "%.0f" , h->Integral()       );
  }else{
    see  << Form( "%.2f" , h->GetBinContent(1) ) << pm << Form( "%.2f" , h->GetBinError(1) );
    smm  << Form( "%.2f" , h->GetBinContent(2) ) << pm << Form( "%.2f" , h->GetBinError(2) );
    sem  << Form( "%.2f" , h->GetBinContent(3) ) << pm << Form( "%.2f" , h->GetBinError(3) );
    stot << Form( "%.2f" , h->Integral()       ) << pm << Form( "%.2f" , histError(h,1,4)  );
  }

  cout << delimstart << setw(width1) << label      << setw(width2)
       << delim      << setw(width1) << see.str()  << setw(width2)
       << delim      << setw(width1) << smm.str()  << setw(width2)
       << delim      << setw(width1) << sem.str()  << setw(width2)
       << delim      << setw(width1) << stot.str() << setw(width2)
       << delimend   << endl;
  
  
}


#include <TList.h>
#include <TIterator.h>

void deleteHistos() {
   // Delete all existing histograms in memory
   TObject* obj;
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();
   while ((obj=iter->Next())) {
     if (obj->IsA()->InheritsFrom(TH1::Class()) ||
         obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
   }
}

void initSymbols( bool latex ){

  //-------------------------------------------------------
  // table format
  //-------------------------------------------------------

  width1      = 15;
  width2      = 4;
  linelength  = (width1+width2)*5+1;

  //-------------------------------------------------------
  // symbols
  //-------------------------------------------------------
  
  if( latex ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }else{
    pm         = " +/- ";
    delim      = "|";
    delimstart = "|";
    delimend   = "|";
    ee         = "ee";
    mm         = "mm";
    em         = "em";
  }

}

void printYields( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , TCut sel , TCut weight , bool latex ){

  initSymbols( latex );

  TCanvas *ctemp = new TCanvas();

  printLine(latex);
  printHeader();
  printLine(latex);

  TH1F* hyield = new TH1F("hyield","yield",4,0,4);
  TH1F* hmctot = new TH1F("hmctot","hmctot",4,0,4);
  hyield->Sumw2();
  hmctot->Sumw2();

  TCut dil("(nels+nmus+ntaus)==2");
  TCut ll("(w1>0&&w2>0) && ntaus==0");
  TCut tau("(w1>0&&w2>0) && ntaus>0");
  TCut fake("!(w1>0&&w2>0)");
  TCut selclone = sel;

  //----------------------
  // print SM MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < chmc.size() ; imc++){

    if( TString(labels[imc]).Contains("LM") ) continue;

    sel = selclone;
    if     ( strcmp(labels[imc],"ttll")   == 0 ) sel = sel + ll;
    else if( strcmp(labels[imc],"tttau")  == 0 ) sel = sel + tau;
    else if( strcmp(labels[imc],"ttfake") == 0 ) sel = sel + fake;
    else if( strcmp(labels[imc],"ttdil")  == 0 ) sel = sel + dil;
    else if( strcmp(labels[imc],"ttotr")  == 0 ) sel = sel + !dil;

    chmc[imc]->Draw("leptype>>hyield",sel*weight);

    if( imc == 0 ) hmctot = (TH1F*) hyield->Clone();
    else           hmctot->Add(hyield);
    
    print( hyield , labels[imc] );

  }

  printLine(latex);

  //-------------------------------
  // print sum of SM MC samples
  //-------------------------------

  print( hmctot , "tot SM MC" );

  printLine(latex);
 
  chdata->Draw("leptype>>hyield",sel);

  print( hyield , "data" );
    
  printLine(latex);

  //----------------------
  // print SUSY MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < chmc.size() ; imc++){

    if( !TString(labels[imc]).Contains("LM") ) continue;

    chmc[imc]->Draw("leptype>>hyield",sel*weight);
    
    print( hyield , labels[imc] );

  }

  printLine(latex);


  
  delete ctemp;
}



void compareDataMC( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , char* var , 
		    TCut sel , TCut weight , int nbins ,  float xmin , float xmax ,  
		    char* xtitle , bool overlayData ){

  int colors[]={2,5,7,4,6,8,9};

  assert( chmc.size() == labels.size() );
  const unsigned int nmc = chmc.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist = new TH1F();
  TH1F*    mchist[nmc];
  TH1F*    datahist = new TH1F("datahist","datahist",nbins,xmin,xmax);

  TLegend *leg = new TLegend(0.8,0.65,0.96,0.94);

  if( overlayData ) leg->AddEntry(datahist,"data","p");
  
  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),nbins,xmin,xmax);

    chmc.at(imc)->Draw(Form("TMath::Min(%s,%f)>>mc_%i",var,xmax-0.01,imc),sel*weight);

    if( TString( labels.at(imc) ).Contains("LM") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    mcstack->Add( mchist[imc] );

    if( imc == 0 ) mctothist = (TH1F*) mchist[imc]->Clone();
    else           mctothist->Add(mchist[imc]);

    leg->AddEntry(mchist[imc],labels.at(imc),"f");

  }

  chdata->Draw(Form("TMath::Min(%s,%f)>>datahist",var,xmax-0.01),sel);

  if( overlayData ){
    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->Draw("E1");
    mcstack->Draw("same");
    datahist->Draw("sameE1");
    datahist->Draw("sameaxis");
  }
  else{
    mctothist->GetXaxis()->SetTitle(xtitle);
    mctothist->Draw();
    mcstack->Draw("same");
    mctothist->Draw("sameaxis");
  }

  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();


}
