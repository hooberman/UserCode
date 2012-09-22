#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
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

using namespace std;

bool fit = false;

// void makePlot(TPad *pad , char* filename , TChain *ch, TCut sel, char* var, char* xtitle, int nbins, float xmin, float xmax, bool log, bool normalize );

/*
void makeAllPlots(TCanvas * can , char* filename , TChain *ch, TCut sel, char* var, char* xtitle, int nbins, float xmin, float xmax ){

  can->cd();
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.5,0.5,1.0);
  makePlot( pad1 , filename , ch , sel , "dilmass" , "M(ll) [GeV]" , 10 , 0 ,  100 ,  false , false );

  can->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.5,0.5,1.0,1.0);
  makePlot( pad2 , filename , ch , sel , "dilmass" , "M(ll) [GeV]" , 10 , 0 ,  100 ,  false , true );

  can->cd();
  TPad *pad3 = new TPad("pad3","pad3",0.0,0.0,0.5,0.5);
  makePlot( pad3 , filename , ch , sel , "dilmass" , "M(ll) [GeV]" , 10 , 0 ,  100 ,   true , false );

  can->cd();
  TPad *pad4 = new TPad("pad4","pad4",0.5,0.0,1.0,0.5);
  makePlot( pad4 , filename , ch , sel , "dilmass" , "M(ll) [GeV]" , 10 , 0 ,  100 ,   true ,  true );

  can->Print(Form("../plots/%s.ps",filename));  
}
*/


void makePlot(TCanvas *can , char* filename , TChain *ch, TCut sel, char* var, char* xtitle, int nbins, float xmin, float xmax ){

  can->cd();
  can->Divide(2,2);

  char* myvar = var;
  if( TString(var).Contains("mbb") ){
    myvar = "sqrt(  pow(bjet1->mass(),2) + pow(bjet2->mass(),2) + 2 * bjet1->E() * bjet2->E() - 2 * bjet1->Px() * bjet2->Px() - 2 * bjet1->Py() * bjet2->Py() - 2 * bjet1->Pz() * bjet2->Pz() )";
  }
  if( TString(var).Contains("st") ){
    myvar = "pfmet+ht40+lep1->pt()+lep2->pt()";
  }
  if( TString(var).Contains("lep1pt") ){
    myvar = "lep1->pt()";
  }
  if( TString(var).Contains("lep2pt") ){
    myvar = "lep2->pt()";
  }
  if( TString(var).Contains("lep1eta") ){
    myvar = "lep1->eta()";
  }
  if( TString(var).Contains("lep2eta") ){
    myvar = "lep2->eta()";
  }

  TCut ee("leptype==0 && ee==1");
  TCut mm("leptype==1 && (mm==1 || mmtk==1)");
  TCut em("leptype==2 && (em==1 || me==1)");
  TCut sf = ee||mm;

  // TPad* mainpad = new TPad(Form("%s_mainpad",var),Form("%s_mainpad",var),0.0,0.0,1.0,0.8);
  // mainpad->Draw();
  // mainpad->cd();
  // if( log ) mainpad->SetLogy();
  can->cd(1);

  TH1F* hsf = new TH1F(Form("%s_sf",var),Form("%s_sf",var),nbins,xmin,xmax);
  TH1F* hof = new TH1F(Form("%s_of",var),Form("%s_of",var),nbins,xmin,xmax);

  hsf->Sumw2();
  hof->Sumw2();
  
  ch->Draw(Form("min(%s,%f)>>%s_sf",myvar,xmax-0.001,var),sel+sf);
  ch->Draw(Form("min(%s,%f)>>%s_of",myvar,xmax-0.001,var),sel+em);

  hsf->SetLineColor(2);
  hsf->SetLineWidth(2);
  hsf->SetMarkerColor(2);

  hof->SetFillColor(7);
  hof->SetMarkerColor(4);

  hsf->GetXaxis()->SetTitle(xtitle);
  hsf->GetYaxis()->SetTitle("entries");

  TH1F* hratio = (TH1F*) hsf->Clone(Form("%s_ratio",var));

  float max = hsf->GetMaximum();
  if( hof->GetMaximum() > max ) max=hof->GetMaximum();
  hsf->SetMaximum(1.3*max);

  hsf->Draw("E1");
  hof->Draw("sameE2");
  hsf->Draw("sameE1");
  hsf->Draw("axissame");

  cout << "SF " << hsf->GetEntries() << endl;
  cout << "OF " << hof->GetEntries() << endl;

  can->cd(2);

  TH1F* hsf_norm = (TH1F*) hsf->Clone(Form("%s_norm",hsf->GetName()));
  TH1F* hof_norm = (TH1F*) hof->Clone(Form("%s_norm",hof->GetName()));

  hsf_norm->Scale(1.0/hsf_norm->Integral());
  hof_norm->Scale(1.0/hof_norm->Integral());

  float max_norm = hsf_norm->GetMaximum();
  if( hof_norm->GetMaximum() > max_norm ) max_norm=hof_norm->GetMaximum();
  hsf_norm->SetMaximum(1.3*max_norm);

  hsf_norm->GetYaxis()->SetTitle("arbitrary units");
  hsf_norm->Draw("E1");
  hof_norm->Draw("sameE2");
  hsf_norm->Draw("sameE1");
  hsf_norm->Draw("axissame");

  can->cd(4);

  // TPad* respad = new TPad(Form("%s_respad",var),Form("%s_respad",var),0.0,0.8,1.0,1.0);
  // respad->Draw();
  // respad->cd();
  // respad->SetGridy();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);

  //TH1F* hratio = (TH1F*) hsf->Clone(Form("%s_ratio",var));
  hratio->Divide(hof);
  // hratio->GetXaxis()->SetLabelSize(0.0);
  // hratio->GetYaxis()->SetLabelSize(0.2);
  // hratio->GetYaxis()->SetTitleSize(0.25);
  // hratio->GetYaxis()->SetTitleOffset(0.25);
  // hratio->GetYaxis()->SetNdivisions(5);
  hratio->GetYaxis()->SetTitle("SF / OF");
  // hratio->GetXaxis()->SetTitle("");
  // hratio->SetMinimum(0.0);
  // hratio->SetMaximum(2.0);
  //hratio->GetYaxis()->SetRangeUser(0.0,2.0);
  hratio->Draw();

  line.DrawLine(xmin,1.0,xmax,1.0);

  // TF1* fline = new TF1("fline","pol1",0,6);
  // fline->SetLineColor(4);
  // fline->SetLineWidth(2);
  // if( fit ) hratio->Fit(fline);

  can->cd(3);
  TH1F* hdiff = (TH1F*) hsf->Clone(Form("%s_ratio",var));
  hdiff->Add(hof,-1);
  hdiff->GetYaxis()->SetTitle("SF - OF");
  hdiff->Draw();

  line.DrawLine(xmin,0.0,xmax,0.0);

  can->Modified();
  can->Update();
  can->Print(Form("../plots/%s.ps",filename));  
  can->Clear();
}


void makeEdgePDF(){

  char* SR = "lowMET";
  //char* SR = "highMET";

  //char* mll = "lowMass";
  char* mll = "highMass";

  cout << "Signal region : " << SR  << endl;
  cout << "Dilepton mass : " << mll << endl;

  gStyle->SetErrorX(0.5);

  //------------------------------------------------------
  // add files
  //------------------------------------------------------

  TChain *data = new TChain("T1");
  data->Add("../output/V00-01-04/data_53X_baby_2jets_met100.root");
  data->Add("../output/V00-01-04/data_2012C_53X_baby_2jets_met100.root");

  TChain *tt = new TChain("T1");
  tt->Add("../output/V00-00-22/ttbar_baby.root");

  //------------------------------------------------------
  // selection
  //------------------------------------------------------

  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut pt2010("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut njets40_2("njets40>=2");
  TCut njets40_3("njets40>=3");
  TCut ht40_100("ht40>=100.0");
  TCut met100("pfmet>100.0");
  TCut met150("pfmet>150.0");
  TCut ee("leptype==0 && (ee==1 || isdata==0)");
  TCut mm("leptype==1 && (mm==1 || mmtk==1 || isdata==0)");
  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut runrange("isdata==0 || (run<197556 || run>198913)");
  TCut mll15to70("dilmass>15.0 && dilmass<70.0");
  TCut mll20to70("dilmass>20.0 && dilmass<70.0");
  TCut mll120("dilmass>120.0");
  TCut nb2("nbcsvm==2");

  TCut sel;

  //------------------------
  // high-MET SR (Aachen)
  //------------------------

  if( TString(SR).Contains("highMET") ){
    sel += runrange;
    sel += filters;
    sel += (ee||mm||em);
    sel += pt2010;
    sel += njets40_2;
    sel += ht40_100;
    sel += met150;

    if( TString(mll).Contains("lowMass") ) sel += mll15to70;
  }

  //------------------------
  // low-MET SR (ETH)
  //------------------------

  if( TString(SR).Contains("lowMET") ){
    sel += runrange;
    sel += filters;
    sel += (ee||mm||em);
    sel += pt2020;
    sel += njets40_3;
    sel += met100;

    if( TString(mll).Contains("lowMass") ) sel += mll20to70;

  }

  if( TString(mll).Contains("highMass") ) sel += mll120;


  cout << "Using selection " << sel.GetTitle() << endl;

  cout << "ee: " << data->GetEntries(sel+ee) << endl;
  cout << "mm: " << data->GetEntries(sel+mm) << endl;
  cout << "em: " << data->GetEntries(sel+em) << endl;

  bool log       = false;
  bool normalize = false;

  char* filename = Form("%s_%s",SR,mll);
  TCanvas* canvas = new TCanvas("canvas","canvas",1000,800);
  //gStyle->SetPaperSize(22,28);

  canvas->Print(Form("../plots/%s.ps[",filename));

  if( TString(mll).Contains("highMass") ){
    makePlot( canvas , filename , data , sel     , "dilmass" , "M(ll) [GeV]"  , 20 , 100 ,  300 );  filename = Form("%s_%s",SR,mll);
  }
  else{
    makePlot( canvas , filename , data , sel     , "dilmass" , "M(ll) [GeV]"  , 10 ,  0 ,  100 );  filename = Form("%s_%s",SR,mll);
  }
  /*
  makePlot( canvas , filename , data , sel  , "lep1pt"  , "1st lepton p_{T} [GeV]"  , 10 ,  0 ,  200 );    filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "lep2pt"  , "2st lepton p_{T} [GeV]"  , 10 ,  0 ,  200 );    filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "lep1eta" , "1st lepton #eta"         , 10 ,  -2.5 ,  2.5 ); filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "lep2eta" , "2st lepton #eta"         , 10 ,  -2.5 ,  2.5 ); filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "njets40" , "n_{jets}"           ,  8 ,  0 ,    8 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "ht40"    , "H_{T} [GeV]"        , 10 ,  0 , 1000 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "nbcsvm"  , "n b-tags"           ,  6 ,  0 ,    6 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "pfmet"   , "E_{T}^{miss} [GeV]" , 20 ,  0 ,  400 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "st"      , "S_{T} [GeV]"        , 15 ,  0 , 1500 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "mt2"     , "MT2 [GeV]"          , 15 ,  0 ,  150 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "mt2j"    , "MT2J [GeV]"         , 20 ,  0 ,  500 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "mlb1"    , "M(l1,any-b) [GeV]"  , 10 ,  0 ,  300 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "mlb2"    , "M(l2,any-b) [GeV]"  , 10 ,  0 ,  300 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel  , "mlbmin"  , "M(l,b)^{min} [GeV]" , 10 ,  0 ,  200 );  filename = Form("%s_%s",SR,mll);
  makePlot( canvas , filename , data , sel+nb2 , "mbb"  , "M(bb) [GeV]"        , 20 , 15 ,  415 );  filename = Form("%s_%s",SR,mll);
  */
  canvas->Print(Form("../plots/%s.ps]",filename));
  //gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.ps ../plots/%s.pdf",filename,filename));


}


