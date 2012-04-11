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
#include "TGraphAsymmErrors.h"
#include <sstream>
#include <iomanip>

using namespace std;

TGraphAsymmErrors* getGraph(TChain* ch, char* var, TCut denom, TCut num){

  TH1F* hall  = new TH1F("hall" ,"",20,0,1);
  TH1F* hpass = new TH1F("hpass","",20,0,1);

  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("min(%s,0.99)>>hall"  , var ) , denom       );
  ch->Draw(Form("min(%s,0.99)>>hpass" , var ) , denom + num );
  delete ctemp;

  cout << "Numerator   : " << TCut(num).GetTitle()       << endl;
  cout << "Entries     : " << hpass->GetEntries()        << endl;
  cout << "Denominator : " << TCut(denom+num).GetTitle() << endl;
  cout << "Entries     : " << hall->GetEntries()         << endl << endl;

  TGraphAsymmErrors* gr = new TGraphAsymmErrors();
  gr->BayesDivide(hpass,hall);

  delete hpass;
  delete hall;

  return gr;

}

void trigs(){

  TChain *chmu = new TChain("t");
  TChain *chel = new TChain("t");
  
  chmu->Add("../output/V00-00-02/muData.root");
  chel->Add("../output/V00-00-02/elData.root");

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  //------------------------------
  // single mu
  //------------------------------

  TH1F* hmuall  = new TH1F("hmuall" ,"",20,0,1);
  TH1F* hmupass = new TH1F("hmupass","",20,0,1);

  chmu->Draw("munoiso1_isopf>>hmuall" ,"mu24>0 && nmunoiso==1 && run>=190659");
  chmu->Draw("munoiso1_isopf>>hmupass","mu24>0 && nmunoiso==1 && munoiso1_mu24==1 && isomu24==1 && run>=190659");

  cout << "muall  " << hmuall->GetEntries() << endl;
  cout << "mupass " << hmupass->GetEntries() << endl;

  /*
  //------------------------------
  // single e
  //------------------------------

  TH1F* helall  = new TH1F("helall" ,"",20,0,1);
  TH1F* helpass = new TH1F("helpass","",20,0,1);

  //chel->Draw("elnoiso1_iso>>helall" ,"el8>0 && nelnoiso==1 && elnoiso1.pt()>30");
  //chel->Draw("elnoiso1_iso>>helpass","el8>0 && nelnoiso==1 && elnoiso1.pt()>30 && elnoiso1_wp80==1 && el27wp80==1");
  chel->Draw("elnoiso1_iso>>helall" ,"el8noisojet30>0 && nelnoiso==1 && elnoiso1.pt()>30");
  chel->Draw("elnoiso1_iso>>helpass","el8noisojet30>0 && nelnoiso==1 && elnoiso1.pt()>30 && elnoiso1_wp80==1 && el27wp80==1");

  cout << "elall  " << helall->GetEntries() << endl;
  cout << "elpass " << helpass->GetEntries() << endl;
  
  //------------------------------
  // e+jets
  //------------------------------

  TH1F* htopall  = new TH1F("htopall" ,"",20,0,1);
  TH1F* htoppass = new TH1F("htoppass","",20,0,1);

  chel->Draw("elnoiso1_iso>>htopall" ,"elnoisotrijet>0 && nelnoiso==1");
  chel->Draw("elnoiso1_iso>>htoppass","elnoisotrijet>0 && nelnoiso==1 && elnoiso1_top==1 && eltrijet==1");

  cout << "topall  " << htopall->GetEntries() << endl;
  cout << "toppass " << htoppass->GetEntries() << endl;

  delete ctemp;

  //------------------------------
  // draw stuff
  //------------------------------
  */
  TGraphAsymmErrors* grmu2 = getGraph(chmu,"munoiso1_isopf",TCut("mu24>0 && nmunoiso==1 && run>=190659"),TCut("munoiso1_mu24==1 && isomu24==1"));

  TCanvas *c1 = new TCanvas("c1","",1800,600);
  c1->Divide(3,1);

  c1->cd(1);
  TGraphAsymmErrors* grmu = new TGraphAsymmErrors();
  grmu->BayesDivide(hmupass,hmuall);
  grmu->SetMarkerSize(3);
  grmu->Draw("AP");



  grmu2->SetLineColor(2);
  grmu2->SetMarkerColor(2);
  grmu2->Draw("sameP");

  /*
  c1->cd(2);
  TGraphAsymmErrors* grel = new TGraphAsymmErrors();
  grel->BayesDivide(helpass,helall);
  grel->Draw("AP");

  c1->cd(3);
  TGraphAsymmErrors* grtop = new TGraphAsymmErrors();
  grtop->BayesDivide(htoppass,htopall);
  grtop->Draw("AP");
  */

  // TH1F* hratio = (TH1F*) helpass->Clone("hratio");
  // hratio->Divide(helall);
  // hratio->Draw("samehist");
}
