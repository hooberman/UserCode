#ifndef ZMET
#define ZMET

#include "TChain.h"
#include "TCut.h"
#include <vector>

TH1F* getHist( TChain* ch , char* var , TCut sel , char* histname , int nbins , float xmin , float xmax );
TH2F* getHist2D( TChain* ch , char* varx , char* vary , TCut sel , char* histname , 
		 int nbinsx , float xmin , float xmax , int nbinsy , float ymin , float ymax );
float calculateHistError( TH1F* h , int minbin , int maxbin );

TChain *data;
TChain *tt;
TChain *zjets;

vector<TChain*> mc;
vector<char*>   mclabels;

#endif
