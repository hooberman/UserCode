#ifndef __CINT__
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TProfile.h"
#include <sstream>
#include <vector>
#endif


void compareSumET( bool printgif = false){

  const unsigned int n = 5;
  TFile* f[n];
  TProfile* p[n];

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);

  for( unsigned int i = 0 ; i < n ; ++i ){

    float thresh = i * 0.4;

    f[i] = TFile::Open( Form( "root/ttbar_hb%f_histos.root" , thresh ) );
    p[i] = (TProfile*) f[i] -> Get("tmet");

    p[i]->SetMarkerSize(0.5);
    p[i]->SetMarkerColor( i + 1 );
    p[i]->SetLineColor( i + 1 );

    stringstream s;
    s << "E_{HB} > " << thresh << " GeV" << endl;

    leg->AddEntry( p[i] , s.str().c_str());

    if( i == 0 ) p[i]->Draw();
    else         p[i]->Draw("same");

  }

  leg->SetFillColor(0);
  leg->Draw();
}
