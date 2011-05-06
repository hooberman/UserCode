#ifndef HOOTILITIES
#define HOOTILITIES

#include "TChain.h"
#include "TCut.h"
#include <vector>

void compareDataMC( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , char* var , 
		    TCut sel , TCut weight , int nbins ,  float xmin , float xmax , char* xtitle , bool overlayData = true , bool drawLegend = true );

void printYields( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , TCut sel , TCut weight , bool latex = false );
void initSymbols(bool);
void  printLine(bool);
void deleteHistos();
TLegend *getLegend( vector<TChain*> chmc , vector<char*> labels , bool overlayData );

char* pm;         
char* delim;      
char* delimstart; 
char* delimend;   
char* ee;         
char* mm;         
char* em;         

int   width1;
int   width2;
int   linelength;

#endif
