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

void makePlot(TChain *ch, TCut sel, char* var, char* xtitle, int nbins, float xmin, float xmax){

  char* myvar = var;
  if( TString(var).Contains("mbb") ){
    //myvar = "sqrt(  pow(bjet1.E()+bjet2.E(),2) - pow(bjet1.Px()+bjet2.Px(),2) - pow(bjet1.Py()+bjet2.Py(),2) - pow(bjet1.Pz()+bjet2.Pz(),2) )";
    myvar = "sqrt(  pow(bjet1->mass(),2) + pow(bjet2->mass(),2) + 2 * bjet1->E() * bjet2->E() - 2 * bjet1->Px() * bjet2->Px() - 2 * bjet1->Py() * bjet2->Py() - 2 * bjet1->Pz() * bjet2->Pz() )";
  }

  TCut ee("leptype==0 && ee==1");
  TCut mm("leptype==1 && mm==1");
  TCut em("leptype==2 && (em==1 || me==1)");
  TCut sf = ee||mm;

  TCanvas *can = new TCanvas(Form("%s_can",var),Form("%s_can",var),600,600);
  can->cd();

  TPad* mainpad = new TPad(Form("%s_mainpad",var),Form("%s_mainpad",var),0.0,0.0,1.0,0.8);
  mainpad->Draw();
  mainpad->cd();

  TH1F* hsf = new TH1F(Form("%s_sf",var),Form("%s_sf",var),nbins,xmin,xmax);
  TH1F* hof = new TH1F(Form("%s_of",var),Form("%s_of",var),nbins,xmin,xmax);

  hsf->Sumw2();
  hof->Sumw2();
  
  // ch->Draw(Form("min(%s,%f)>>%s_sf",myvar,xmax-0.001,var),sel+sf);
  // ch->Draw(Form("min(%s,%f)>>%s_of",myvar,xmax-0.001,var),sel+em);

  ch->Draw(Form("min(%s,%f)>>%s_sf",myvar,xmax-0.001,var),sel+ee);
  ch->Draw(Form("min(%s,%f)>>%s_of",myvar,xmax-0.001,var),sel+mm);

  ch->SetScanField(1000);
  //ch->Scan(Form("event:run:lumi:%s:leptype",myvar),sel,"col=20u");

  hsf->Scale(1.0/hsf->Integral());
  hof->Scale(1.0/hof->Integral());

  hsf->SetMaximum(0.0);
  float max = hsf->GetMaximum();
  if( hof->GetMaximum() > max ) max=hof->GetMaximum();
  hsf->SetMaximum(1.3*max);

  hsf->SetLineColor(2);
  hsf->SetLineWidth(2);
  hsf->SetMarkerColor(2);

  hof->SetFillColor(7);
  hof->SetMarkerColor(4);

  hsf->GetXaxis()->SetTitle(xtitle);
  // hsf->Draw("E1");
  // hof->Draw("samehist");
  // hsf->Draw("sameE1");

  hsf->GetYaxis()->SetTitle("A.U.");
  //hsf->SetMaximum(0.3);
  hsf->Draw("E1");
  hof->Draw("sameE2");
  hsf->Draw("sameE1");
  hsf->Draw("axissame");

  TLine line;
  //line.DrawLine(125,0,125,0.3);
    
  cout << "SF " << hsf->GetEntries() << endl;
  cout << "OF " << hof->GetEntries() << endl;

  can->cd();

  TPad* respad = new TPad(Form("%s_respad",var),Form("%s_respad",var),0.0,0.8,1.0,1.0);
  respad->Draw();
  respad->cd();
  respad->SetGridy();

  TH1F* hratio = (TH1F*) hsf->Clone(Form("%s_ratio",var));
  hratio->Divide(hof);
  hratio->GetXaxis()->SetLabelSize(0.0);
  hratio->GetYaxis()->SetLabelSize(0.2);
  hratio->GetYaxis()->SetTitleSize(0.25);
  hratio->GetYaxis()->SetTitleOffset(0.25);
  hratio->GetYaxis()->SetNdivisions(5);
  hratio->GetYaxis()->SetTitle("SF / OF");
  hratio->GetXaxis()->SetTitle("");
  hratio->SetMinimum(0.0);
  hratio->SetMaximum(2.0);
  hratio->GetYaxis()->SetRangeUser(0.6,1.4);
  hratio->Draw();

  TF1* fline = new TF1("fline","pol1",0,6);
  fline->SetLineColor(4);
  fline->SetLineWidth(2);
  hratio->Fit(fline);

  //can->Print(Form("../plots/edge_%s.pdf",var));
}


void makeEdgePlots(){

  gStyle->SetErrorX(0.5);

  TChain *data = new TChain("T1");
  //data->Add("../output/V00-00-24/dataskim2010_baby.root");
  //data->Add("../output/V00-00-24/ttbar_massiveb_baby.root");
  //data->Add("../output/V00-01-02/ttbar_53X_baby.root");
  //data->Add("../output/V00-01-02/data_53X_baby.root");
  data->Add("../output/V00-01-02/data_53X_baby.root");
  data->Add("../output/V00-01-02/data_2012C_53X_baby.root");

  TChain *tt = new TChain("T1");
  tt->Add("../output/V00-00-22/ttbar_baby.root");

  //TCut sel("njets40>=2 && ht40>100 && pfmet>150 && dilmass<70.0" );
  //TCut sel("nbcsvm>=1 && pfmet>50 && dilmass>12.0 && (dilmass<76.0||dilmass>106.0) && lep1.pt()>20.0 && lep2.pt()>20.0" );

  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut eta24("abs(lep1.eta())<2.4 && abs(lep2.eta())<2.4");
  TCut Zveto("dilmass<76.0 || dilmass>106.0");
  TCut mll120("dilmass>120.0");
  TCut nb1("nbcsvm>=1");
  TCut metcut("pfmet>100.0");
  TCut runrange("isdata==0 || (run<197556 || run>198913)");

  TCut sel;
  sel += pt2020;
  sel += eta24;
  sel += Zveto;
  sel += mll120;
  sel += nb1;
  sel += metcut;
  sel += runrange;

  cout << "Using selection " << sel.GetTitle() << endl;

  //TCut sel("njets40>=2 && ht40>100 && pfmet>150 && dilmass>81 && dilmass<101");
  TCut pflep("pflep1.pt()>20 && pflep2.pt()>10");
  TCut ee("leptype==0 && ee==1");
  TCut mm("leptype==1 && mm==1");
  TCut em("leptype==2 && (em==1 || me==1)");
  TCut eemc("leptype==0");
  TCut mmmc("leptype==1");
  TCut emmc("leptype==2");
  TCut nb2("nbcsvm==2");
  TCut weight("5.1 * vtxweight * trgeff * weight");

  cout << "ee: " << data->GetEntries(sel+ee) << endl;
  cout << "mm: " << data->GetEntries(sel+mm) << endl;
  cout << "em: " << data->GetEntries(sel+em) << endl;

  // makePlot( data , sel , "dilmass" , "M(ll) [GeV]" , 10 , 0 ,  100 );
  makePlot( data , sel     , "njets40" , "n_{jets}"           , 6 , 0 ,   6 );
  // makePlot( data , sel     , "ht40"    , "H_{T} [GeV]"        ,  5 , 0 , 1000 );
  // makePlot( data , sel     , "nbcsvm"  , "n b-tags"           ,  6 , 0 ,    6 );
  // makePlot( data , sel     , "mt2"     , "MT2 [GeV]"          , 15 , 0 ,  150 );
  // makePlot( data , sel     , "mt2j"    , "MT2J [GeV]"         , 20 , 0 ,  500 );
  // makePlot( data , sel     , "mlb1"    , "M(l1,any-b) [GeV]"  , 10 , 0 ,  300 );
  // makePlot( data , sel     , "mlb2"    , "M(l2,any-b) [GeV]"  , 10 , 0 ,  300 );
  // makePlot( data , sel     , "mlbmin"  , "M(l,b)^{min} [GeV]" , 10 , 0 ,  200 );
  // makePlot( data , sel+nb2 , "mbb"     , "M(bb) [GeV]"        , 20 , 15 ,  415 );


  /*
  // TH1F* hee = new TH1F("hee","",32,0,320);
  // TH1F* hmm = new TH1F("hmm","",32,0,320);
  // TH1F* hem = new TH1F("hem","",32,0,320);

  TH1F* hee = new TH1F("hee","",15,0,300);
  TH1F* hmm = new TH1F("hmm","",15,0,300);
  TH1F* hem = new TH1F("hem","",15,0,300);

  TH1F* heemc = new TH1F("heemc","",15,0,300);
  TH1F* hmmmc = new TH1F("hmmmc","",15,0,300);
  TH1F* hemmc = new TH1F("hemmc","",15,0,300);

  hee->Sumw2();
  hmm->Sumw2();
  hem->Sumw2();

  heemc->Sumw2();
  hmmmc->Sumw2();
  hemmc->Sumw2();

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  data->Draw("min(dilmass,319)>>hee",sel+ee);
  data->Draw("min(dilmass,319)>>hmm",sel+mm);
  data->Draw("min(dilmass,319)>>hem",sel+em);

  tt->Draw("min(dilmass,319)>>heemc",(sel+eemc)*weight);
  tt->Draw("min(dilmass,319)>>hmmmc",(sel+mmmc)*weight);
  tt->Draw("min(dilmass,319)>>hemmc",(sel+emmc)*weight);
  delete ctemp;

  hee->SetMaximum(100);
  hmm->SetMaximum(100);
  hem->SetMaximum(100);

  heemc->SetFillColor(5);
  hmmmc->SetFillColor(5);
  hemmc->SetFillColor(5);

  TCanvas *c1 = new TCanvas("c1","",1800,600);
  c1->Divide(3,1);

  c1->cd(1);
  hee->Draw();
  hee->GetXaxis()->SetTitle("M(ee) [GeV]");
  heemc->Draw("samehist");
  hee->Draw("sameE1");
  hee->Draw("axissame");
  

  c1->cd(2);
  hmm->Draw();
  hmm->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
  hmmmc->Draw("samehist");
  hmm->Draw("sameE1");
  hmm->Draw("axissame");

  c1->cd(3);
  hem->Draw();
  hem->GetXaxis()->SetTitle("M(e#mu) [GeV]");
  hemmc->Draw("samehist");
  hem->Draw("sameE1");
  hem->Draw("axissame");
  */

  /*
  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  
  TH1F* hdiff = (TH1F*) hee->Clone("hdiff");
  hdiff->Add(hmm);
  hdiff->Add(hem,-1);

  hdiff->GetXaxis()->SetTitle("dilepton mass [GeV]");
  hdiff->GetYaxis()->SetTitle("N(ee) + N(#mu#mu) - N(e#mu)");
  hdiff->Draw();

  TLine line;
  line.SetLineColor(2);
  line.DrawLine(0,0,150,0);
  */







}
