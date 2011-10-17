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
#include "singleLepSusy.C"

using namespace std;

void saveCanvasPDF(TCanvas *can,char* path,char* title){

  can->cd();
  gStyle->SetPaperSize(16,8);
  can->Print(Form("%s/%s.ps",path,title));
  gROOT->ProcessLine(Form(".! ps2pdf %s/%s.ps %s/%s.pdf",path,title,path,title));
}

void makeAnalysisPlots( char* path , bool print = false ){

  //-----------------------------
  // load data/MC samples
  //-----------------------------

  initialize(path);

  //-----------------------------
  // selection
  //-----------------------------

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
  //TCut trkreliso02("trkreliso10>0.2");
  TCut trkveto_pt5_iso01("trkreliso5>0.1 || trkpt5<0");
  TCut trkveto_pt10_iso01("trkreliso10>0.1 || trkpt10<0");
  TCut weight("weight * 0.98 * ndavtxweight");

  TCut presel = nlep1 + leppt + njets3 + met50;

  // TCut sel;
  // sel    += nlep1;
  // sel    += leppt;
  // sel    += njets3;
  // sel    += met50;
  // //sel    += met100;
  // //sel    += met60;
  // //sel    += dphi05;
  // //sel    += met100;
  // //sel    += ht300;
  // //sel    += ht500;
  // //sel    += btags0;
  // //sel    += btags1;
  // sel    += btags2;
  // //sel    += mt150;
  // //sel    += mt150;
  // //sel    += mt200;
  // //sel    += trkreliso01;
  // //sel    += trkreliso02;

  //-----------------------------------------------------
  // plot MT passing vs. failing track veto cut
  //-----------------------------------------------------

  bool log      = true;
  bool combine6 = false;
  bool residual = true;

  TCut mtsel = presel+btags2;

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);

  c1->cd(1);
  compareDataMC( mc , mclabels , data , "mt" , TCut(mtsel+trkveto_pt5_iso01) , weight , 40 , 0 , 400 , "M_{T} (GeV)" , true , residual , !combine6 , log );

  c1->cd(2);
  compareDataMC( mc , mclabels , data , "mt" , TCut(mtsel+!trkveto_pt5_iso01) , weight , 40 , 0 , 400 , "M_{T} (GeV)" , true , residual , !combine6 , log );

  if( print ) saveCanvasPDF(c1,"../plots","mt_trkveto_pt5_reliso01");

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  compareDataMC( mc , mclabels , data , "mt" , TCut(mtsel+trkveto_pt10_iso01) , weight , 40 , 0 , 400 , "M_{T} (GeV)" , true , residual , !combine6 , log );

  c2->cd(2);
  compareDataMC( mc , mclabels , data , "mt" , TCut(mtsel+!trkveto_pt10_iso01) , weight , 40 , 0 , 400 , "M_{T} (GeV)" , true , residual , !combine6 , log );

  if( print ) saveCanvasPDF(c2,"../plots","mt_trkveto_pt10_reliso01");



}
