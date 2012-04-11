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

TGraphAsymmErrors* getGraph(TChain* ch, char* var, TCut denom, TCut num, int rebin = 1){

  TH1F* hall  = new TH1F("hall" ,"",20,0,1);
  TH1F* hpass = new TH1F("hpass","",20,0,1);

  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("max(%s,0.0001)>>hall"  , var ) , denom       );
  ch->Draw(Form("max(%s,0.0001)>>hpass" , var ) , denom + num );
  delete ctemp;

  cout << endl;
  cout << "Variable    : " << var                        << endl;
  cout << "Numerator   : " << TCut(num).GetTitle()       << endl;
  cout << "Entries     : " << hpass->GetEntries()        << endl;
  cout << "Denominator : " << TCut(denom+num).GetTitle() << endl;
  cout << "Entries     : " << hall->GetEntries()         << endl;
  cout << "Iso < 0.1   : " << hpass->Integral(1,2)/hall->Integral(1,2) << endl;
  cout << "Iso < 0.15  : " << hpass->Integral(1,3)/hall->Integral(1,3) << endl << endl;

  if( rebin > 1 ){
    hall->Rebin(rebin);
    hpass->Rebin(rebin);
  }

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
  // single e
  //------------------------------
  /*
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

  */

  //------------------------------
  // single muon
  //------------------------------

  /*
  TGraphAsymmErrors* grmu_iso    = getGraph(chmu,"munoiso1_iso"   , TCut("mu24>0 && nmunoiso==1 && run>=190659"),TCut("munoiso1_mu24==1 && isomu24==1"));
  TGraphAsymmErrors* grmu_isovtx = getGraph(chmu,"munoiso1_isovtx", TCut("mu24>0 && nmunoiso==1 && run>=190659"),TCut("munoiso1_mu24==1 && isomu24==1"));
  TGraphAsymmErrors* grmu_isopf  = getGraph(chmu,"munoiso1_isopf" , TCut("mu24>0 && nmunoiso==1 && run>=190659"),TCut("munoiso1_mu24==1 && isomu24==1"));
  TGraphAsymmErrors* grmu_isofj  = getGraph(chmu,"munoiso1_isofj" , TCut("mu24>0 && nmunoiso==1 && run>=190659"),TCut("munoiso1_mu24==1 && isomu24==1"));

  TCanvas *c1 = new TCanvas();

  c1->cd(1);

  grmu_iso->SetLineColor(1);
  grmu_iso->SetMarkerColor(1);
  grmu_isopf->SetLineColor(2);
  grmu_isopf->SetMarkerColor(2);
  grmu_isopf->SetMarkerStyle(25);
  grmu_isovtx->SetLineColor(4);
  grmu_isovtx->SetMarkerColor(4);
  grmu_isofj->SetLineColor(8);
  grmu_isofj->SetMarkerColor(8);
  grmu_isofj->SetMarkerStyle(26);

  grmu_iso->GetXaxis()->SetTitle("reliso");
  grmu_iso->Draw("AP");
  //grmu_isovtx->Draw("sameP");
  //grmu_isopf->Draw("sameP");
  //grmu_isofj->Draw("sameP");
  grmu_isopf->Draw("sameP");

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(grmu_iso,"det iso","lp");
  //leg->AddEntry(grmu_isovtx,"det iso vtx corr","lp");
  //leg->AddEntry(grmu_isofj,"det iso fast-jet","lp");
  leg->AddEntry(grmu_isopf,"pf iso","lp");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  c1->Print("../plots/mu.pdf");
  */

  //------------------------------
  // single electron
  //------------------------------

  /*
  TCut eldenom("el8noiso>0 && nelnoiso==1 && elnoiso1.pt()>30");
  TCut eldenom("el8noisojet30>0 && nelnoiso==1 && elnoiso1.pt()>30");
  TCut eldenom("el8noiso>0 && nelnoiso==1 && elnoiso1.pt()>30");
  TCut elnum("elnoiso1_wp80==1 && el27wp80==1");

  TGraphAsymmErrors* grel_iso     = getGraph(chel,"elnoiso1_iso"    , eldenom , elnum , 2 );
  TGraphAsymmErrors* grel_isofj   = getGraph(chel,"elnoiso1_isofj"  , eldenom , elnum , 2 );
  TGraphAsymmErrors* grel_isovtx  = getGraph(chel,"elnoiso1_isovtx" , eldenom , elnum , 2 );
  TGraphAsymmErrors* grel_isopf   = getGraph(chel,"elnoiso1_isopf"  , eldenom , elnum , 2 );

  grel_iso->SetLineColor(1);
  grel_iso->SetMarkerColor(1);
  grel_isopf->SetLineColor(2);
  grel_isopf->SetMarkerColor(2);
  grel_isopf->SetMarkerStyle(25);
  grel_isovtx->SetLineColor(4);
  grel_isovtx->SetMarkerColor(4);
  grel_isovtx->SetMarkerStyle(22);
  grel_isofj->SetLineColor(8);
  grel_isofj->SetMarkerColor(8);
  grel_isofj->SetMarkerStyle(26);

  TCanvas *c2 = new TCanvas();
  c2->cd();

  grel_iso->GetXaxis()->SetTitle("reliso");
  grel_iso->Draw("AP");
  grel_isovtx->Draw("sameP");
  //grel_isofj->Draw("sameP");
  //grel_isopf->Draw("sameP");

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry(grel_iso    , "det iso"          , "lp" );
  leg->AddEntry(grel_isovtx , "det iso vtx corr" , "lp" );
  //leg->AddEntry(grel_isofj  , "det iso fast-jet" , "lp" );
  //leg->AddEntry(grel_isopf  , "pf iso"           , "lp" );
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();
  */


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
