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

using namespace std;

TGraph *observedLimit_OS2011(){

  const unsigned int no = 18;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 125;  yo[i] = 435;
  xo[++i] = 150;  yo[i] = 420;
  xo[++i] = 180;  yo[i] = 380;
  xo[++i] = 190;  yo[i] = 350;
  xo[++i] = 195;  yo[i] = 330;
  xo[++i] = 210;  yo[i] = 320;
  xo[++i] = 250;  yo[i] = 330;
  xo[++i] = 300;  yo[i] = 340;
  xo[++i] = 450;  yo[i] = 360;
  xo[++i] = 580;  yo[i] = 360;
  xo[++i] = 600;  yo[i] = 350;
  xo[++i] = 590;  yo[i] = 300;
  xo[++i] = 800;  yo[i] = 260;
  xo[++i] = 1000; yo[i] = 230;
  xo[++i] = 1400; yo[i] = 190;
  xo[++i] = 2000; yo[i] = 180;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph *expectedLimit_OS2011(){
  return observedLimit_OS2011();
}
  
TGraph *expectedLimitP1_OS2011(){

  const unsigned int no = 20;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 410;
  xo[++i] = 80;   yo[i] = 410;
  xo[++i] = 90;   yo[i] = 400;
  xo[++i] = 100;  yo[i] = 390;
  xo[++i] = 100;  yo[i] = 370;
  xo[++i] = 120;  yo[i] = 365;
  xo[++i] = 145;  yo[i] = 360;
  xo[++i] = 180;  yo[i] = 350;
  xo[++i] = 185;  yo[i] = 315;
  xo[++i] = 190;  yo[i] = 280;
  xo[++i] = 210;  yo[i] = 275;
  xo[++i] = 240;  yo[i] = 280;
  xo[++i] = 400;  yo[i] = 320;
  xo[++i] = 460;  yo[i] = 320;
  xo[++i] = 480;  yo[i] = 270;
  xo[++i] = 550;  yo[i] = 245;
  xo[++i] = 600;  yo[i] = 240;
  xo[++i] = 800;  yo[i] = 210;
  xo[++i] = 1200; yo[i] = 150;
  xo[++i] = 2000; yo[i] = 140;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}
  
TGraph *expectedLimitM1_OS2011(){

  const unsigned int no = 18;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 200;  yo[i] = 420;
  xo[++i] = 220;  yo[i] = 400;
  xo[++i] = 225;  yo[i] = 350;
  xo[++i] = 245;  yo[i] = 347;
  xo[++i] = 265;  yo[i] = 345;
  xo[++i] = 300;  yo[i] = 350;
  xo[++i] = 400;  yo[i] = 370;
  xo[++i] = 450;  yo[i] = 375;
  xo[++i] = 550;  yo[i] = 372;
  xo[++i] = 640;  yo[i] = 370;
  xo[++i] = 730;  yo[i] = 300;
  xo[++i] = 900;  yo[i] = 270;
  xo[++i] = 1000; yo[i] = 250;
  xo[++i] = 1200; yo[i] = 230;
  xo[++i] = 1400; yo[i] = 220;
  xo[++i] = 2000; yo[i] = 210;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph *observedLimitTheoryUp_OS2011(){
  return expectedLimitM1_OS2011();
}

TGraph *observedLimitTheoryDown_OS2011(){

  const unsigned int no = 17;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 430;
  xo[++i] = 100;  yo[i] = 420;
  xo[++i] = 110;  yo[i] = 410;
  xo[++i] = 120;  yo[i] = 390;
  xo[++i] = 180;  yo[i] = 360;
  xo[++i] = 190;  yo[i] = 310;
  xo[++i] = 230;  yo[i] = 305;
  xo[++i] = 320;  yo[i] = 310;
  xo[++i] = 400;  yo[i] = 320;
  xo[++i] = 500;  yo[i] = 330;
  xo[++i] = 550;  yo[i] = 270;
  xo[++i] = 600;  yo[i] = 265;
  xo[++i] = 800;  yo[i] = 240;
  xo[++i] = 1000; yo[i] = 220;
  xo[++i] = 1200; yo[i] = 180;
  xo[++i] = 1400; yo[i] = 170;
  xo[++i] = 2000; yo[i] = 170;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}


TGraph *expectedLimitTheoryUp_OS2011(){

  const unsigned int no = 18;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 180;  yo[i] = 420;
  xo[++i] = 200;  yo[i] = 400;
  xo[++i] = 205;  yo[i] = 350;
  xo[++i] = 225;  yo[i] = 347;
  xo[++i] = 245;  yo[i] = 345;
  xo[++i] = 300;  yo[i] = 350;
  xo[++i] = 400;  yo[i] = 370;
  xo[++i] = 450;  yo[i] = 375;
  xo[++i] = 550;  yo[i] = 372;
  xo[++i] = 640;  yo[i] = 350;
  xo[++i] = 670;  yo[i] = 310;
  xo[++i] = 800;  yo[i] = 270;
  xo[++i] = 1000; yo[i] = 240;
  xo[++i] = 1200; yo[i] = 220;
  xo[++i] = 1400; yo[i] = 200;
  xo[++i] = 2000; yo[i] = 200;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph* expectedLimitTheoryDown_OS2011(){
  return observedLimitTheoryDown_OS2011();
}


void CMSSM(){

  //----------------------------------------------------
  // get TH2's
  //----------------------------------------------------

  TFile *fnom = TFile::Open("/tas/benhoob/home/LandS/OS_LandS/fullShapeAnalysis/cards/V00-00-06/observed_limit.root");
  TFile *fup  = TFile::Open("/tas/benhoob/home/LandS/OS_LandS/fullShapeAnalysis/cards/V00-00-07/observed_limit.root");
  TFile *fdn  = TFile::Open("/tas/benhoob/home/LandS/OS_LandS/fullShapeAnalysis/cards/V00-00-08/observed_limit.root");

  TH2F* hobs   = (TH2F*) fnom->Get("hexcl");  
  TH2F* hexp   = (TH2F*) fnom->Get("hexp");  
  TH2F* hexpp1 = (TH2F*) fnom->Get("hexpp1");  
  TH2F* hexpm1 = (TH2F*) fnom->Get("hexpm1");  
  TH2F* hobsup = (TH2F*) fup->Get("hexcl");  
  TH2F* hobsdn = (TH2F*) fdn->Get("hexcl");  
  TH2F* hexpup = (TH2F*) fup->Get("hexp");  
  TH2F* hexpdn = (TH2F*) fdn->Get("hexp");  

  hobs->GetXaxis()->SetTitle("m_{0} [GeV]");
  hobs->GetYaxis()->SetTitle("m_{1/2} [GeV]");

  //----------------------------------------------------
  // no EWSB region
  //----------------------------------------------------

  double st_m0_tanBeta10[]  = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
  double st_m12_tanBeta10[] = {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0}; 
  TGraph* st_gr_tanBeta10 = new TGraph(25,st_m0_tanBeta10,st_m12_tanBeta10);
  
  st_gr_tanBeta10->SetFillColor(40);
  st_gr_tanBeta10->SetFillStyle(1001);

  TLatex *t = new TLatex();
  t->SetNDC();

  //----------------------------------------------------
  // draw plots
  //----------------------------------------------------
  

  //observed limit
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hobs->Draw("colz");
  TGraph* gobs = observedLimit_OS2011();
  gobs->Draw("c");
  gobs->Draw("p");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"observed limit");
  /*
  //expected limit
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexp->Draw("colz");
  TGraph* gexp = expectedLimit_OS2011();
  gexp->Draw("c");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"expected limit");
  
  //expected limit (+1)
  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  c3->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexpp1->Draw("colz");
  TGraph* gexpp1 = expectedLimitP1_OS2011();
  gexpp1->Draw("c");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"expected limit (-1#sigma)");

  //expected limit (-1)
  TCanvas *c4 = new TCanvas("c4","c4",1200,600);
  c4->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexpm1->Draw("colz");
  TGraph* gexpm1 = expectedLimitM1_OS2011();
  gexpm1->Draw("c");
  gexpm1->Draw("p");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"expected limit (+1#sigma)");

  //observed theory up
  TCanvas *c5 = new TCanvas("c5","c5",1200,600);
  c5->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hobsup->Draw("colz");
  TGraph* gobsup = observedLimitTheoryUp_OS2011();
  gobsup->Draw("c");
  gobsup->Draw("p");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"observed limit (theory UP)");

  //observed theory down
  TCanvas *c6 = new TCanvas("c6","c6",1200,600);
  c6->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hobsdn->Draw("colz");
  TGraph* gobsdn = observedLimitTheoryDown_OS2011();
  gobsdn->Draw("c");
  gobsdn->Draw("p");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"observed limit (theory DOWN)");

  //expected theory down
  TCanvas *c8 = new TCanvas("c8","c8",1200,600);
  c8->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexpup->Draw("colz");
  TGraph* gexpup = expectedLimitTheoryUp_OS2011();
  gexpup->Draw("c");
  gexpup->Draw("p");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"experved limit (theory UP)");

  //expected theory down
  TCanvas *c7 = new TCanvas("c7","c7",1200,600);
  c7->cd();
  gPad->SetRightMargin(0.1);
  gPad->SetGridx();
  gPad->SetGridy();
  hexpdn->Draw("colz");
  TGraph* gexpdn = expectedLimitTheoryDown_OS2011();
  gexpdn->Draw("c");
  gexpdn->Draw("p");
  st_gr_tanBeta10->Draw("fsame");
  t->DrawLatex(0.5,0.7,"experved limit (theory DOWN)");

  */
    


}
