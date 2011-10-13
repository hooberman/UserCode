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
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <sstream>

using namespace std;

void makePlot( TChain *ch , TCut presel , TCut num , string title );

//------------------------------------------------------------------
// this function returns the minimum non-zero bin contents of a TH2
//------------------------------------------------------------------

float getMinimum( TH2F* h ){

  float min = 1e10;

  for(int xbin = 1 ; xbin <= h->GetXaxis()->GetNbins() ; xbin++ ){
    for(int ybin = 1 ; ybin <= h->GetYaxis()->GetNbins() ; ybin++ ){
      float val = h->GetBinContent(xbin,ybin);
      if( val < 1e-10 ) continue;
      if( val < min ) min = val;
    }
  }

  return min;
}


void makeEfficiencyPlots( bool print = false ){

  gStyle->SetPaintTextFormat(".2f");
  
  //--------------------------------
  // input file
  //--------------------------------

  TChain *ch = new TChain("t");
  ch->Add("../output/V00-01-01/T2tt_smallTree.root");

  //--------------------------------
  // selection
  //--------------------------------
  
  TCut nlep1("ngoodlep == 1");
  TCut leppt("(leptype==0 && lep1.pt()>25)||(leptype==1 && lep1.pt()>20)");
  TCut njets1("ncalojets >= 1");
  TCut njets2("ncalojets >= 2");
  TCut njets3("ncalojets >= 3");
  TCut njets4("ncalojets >= 4");
  TCut met60("pfmet > 60");
  TCut met50("pfmet > 50");
  TCut met100("pfmet > 100");
  TCut metpresel("(leptype==0&&pfmet>30)||(leptype==1&&pfmet>20)");
  TCut ht300("ht > 300");
  TCut ht500("ht > 500");
  TCut dphi05("dphijm > 0.5");
  TCut btags0("nbtags==0");
  TCut btags2("nbctcm>=2");
  TCut btags1("nbctcm>=1");
  TCut mt100("mt>100");
  TCut mt150("mt>150");
  TCut mt200("mt>200");
  TCut trkreliso02("trkreliso10>0.2");
  TCut trkreliso01("trkreliso5>0.1");

  TCut presel;
  presel    += nlep1;
  presel    += leppt;

  //-------------------------------------------
  // define cuts for efficiency measurement
  //-------------------------------------------

  vector<TCut>   cuts;
  vector<string> titles;
  vector<string> labels;

  cuts.push_back(njets3);            titles.push_back("njets #geq 3");                 labels.push_back("njets3");
  cuts.push_back(njets4);            titles.push_back("njets #geq 4");                 labels.push_back("njets4");
  cuts.push_back(met50);             titles.push_back("E_{T}^{miss} > 50 GeV");        labels.push_back("met50");
  cuts.push_back(trkreliso02);       titles.push_back("trkiso > 0.2 (p_{T} > 10 GeV)");  labels.push_back("trkiso02");
  cuts.push_back(trkreliso01);       titles.push_back("trkiso > 0.1 (p_{T} > 5 GeV)");   labels.push_back("trkiso01");
  cuts.push_back(btags1);            titles.push_back("#geq 1 b-tag");                 labels.push_back("nbtags1");
  cuts.push_back(btags2);            titles.push_back("#geq 2 b-tags");                labels.push_back("nbtags2");
  cuts.push_back(mt100);             titles.push_back("M_{T} > 100 GeV");              labels.push_back("mt100");
  cuts.push_back(mt150);             titles.push_back("M_{T} > 150 GeV");              labels.push_back("mt150");
  cuts.push_back(mt200);             titles.push_back("M_{T} > 200 GeV");              labels.push_back("mt200");

  const unsigned int nplots = cuts.size();

  //-------------------------------------------
  // draw plots
  //-------------------------------------------

  TCanvas *can[nplots];

  for( unsigned int i = 0 ; i < nplots ; i++){

    cout << "Plotting: " << cuts.at(i).GetTitle() << " " << titles.at(i) << endl;

    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),600,600);
    can[i]->cd();

    makePlot( ch , presel , cuts.at(i) , titles.at(i) );

    if( print ){
      can[i]->Print(Form("../plots/makeEfficiencyPlots_%s.eps",labels.at(i).c_str()));
      gROOT->ProcessLine(Form(".! ps2pdf ../plots/makeEfficiencyPlots_%s.eps ../plots/makeEfficiencyPlots_%s.pdf",labels.at(i).c_str(),labels.at(i).c_str()));
      can[i]->Print(Form("../plots/makeEfficiencyPlots_%s.png",labels.at(i).c_str()));
    }
  }

}


void makePlot( TChain *ch , TCut presel , TCut num , string title ){

  TH2F* hdenom = new TH2F("hdenom","hdenom",12,225,525,12,50,350);
  TH2F* hnum   = new TH2F("hnum"  ,"hnum"  ,12,225,525,12,50,350);

  ch->Draw("ml:mg>>hdenom",presel);
  ch->Draw("ml:mg>>hnum"  ,TCut(presel+num));

  hnum->Divide(hdenom);
  hnum->GetXaxis()->SetTitle("m( #tilde{t} ) (GeV)");
  hnum->GetYaxis()->SetTitle("m(#chi_{1}^{0}) (GeV)");
  hnum->SetMinimum( getMinimum(hnum) );
  hnum->DrawCopy("colz");
  hnum->DrawCopy("sametext");

  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.04);
  t.DrawLatex(.2,0.9,title.c_str());
  
  delete hdenom;
  delete hnum;

}
