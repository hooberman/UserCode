#ifndef OSSUSY
#define OSSUSY

#include "TChain.h"
#include "TCut.h"
#include <vector>

TH1F* getHist( TChain* ch , char* var , TCut sel , char* histname , int nbins , float xmin , float xmax );
TH2F* getHist2D( TChain* ch , char* varx , char* vary , TCut sel , char* histname , 
		 int nbinsx , float xmin , float xmax , int nbinsy , float ymin , float ymax );
float calculateHistError( TH1F* h , int minbin , int maxbin );

void  printABCDHeader();

//TH2F* doABCD( TChain *ch , char* label , TCut sel , TCut weight , TCut A , TCut B , TCut C , TCut D , TH1F* hyield = 0 );
TH2F* doABCD( TChain *ch , char* label , TCut sel , TCut weight , TCut A , TCut B , TCut C , TCut D );
void printABCDRow( char* sample, float A, float B, float C, float D, float dA, float dB, float dC, float dD );
void drawSquare(float x1, float y1, float x2, float y2, int color = 1);
void plotHist( TH1F* h1 , TH1F* h2 , char* leg1 , char* leg2 , char* xtitle , bool residual = false);
void plotVictoryHist( TH1F* hpred , TH1F* hobs , TH1F* hpred_data, TH1F* hobs_data,
		      string title, string xtitle , int rebin , float metcut);
void printVictoryRow( string sample , TH1F* hobs , TH1F* hpred , float cut , float mcratio = 1. );
pair<float,float> getMCratio( TH1F* hobs , TH1F* hpred , float cut );


void printVictoryHeader(){
  
  cout << delimstart << setw(width1) << "Sample"        << setw(width2)
       << delim      << setw(width1) << "Predicted"     << setw(width2)
       << delim      << setw(width1) << "Observed"      << setw(width2)
       << delim      << setw(width1) << "Obs/pred"      << setw(width2)
       << delimend   << endl;
}


TChain *data;
TChain *datasf;
TChain *datadf;
TChain *data2010;
TChain *ttall;
TChain *ttpowheg;
TChain *ttdil;
TChain *ttotr;
TChain *ttll;
TChain *tttau;
TChain *ttfake;
TChain *zjets;
TChain *dy;
TChain *dydata;
TChain *dytautau;
TChain *ww;
TChain *wz;
TChain *zz;
TChain *t;
TChain *wjets;
TChain *LM0;
TChain *LM1;
TChain *LM2;
TChain *LM3;
TChain *LM4;
TChain *LM5;
TChain *LM6;
TChain *LM7;
TChain *LM8;
TChain *LM9;
TChain *LM10;
TChain *LM11;
TChain *LM12;
TChain *LM13;
TChain *sm;
TChain *sm_LM1;
TChain *other;

vector<TChain*> mc;
vector<char*>   mclabels;
vector<char*>   mctex;


#endif
