#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include "TLegendEntry.h"
#include "THStack.h"
#include "TStyle.h"
#include <sstream>
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TCut.h"
#include "TTree.h"

using namespace std;

//147217 *        81 *  60186645

void compareMet(){

  gStyle->SetOptStat(0);

  bool isData = true;

  char* file = "output/v8/ZJets_baby.root";
  if( isData ) file = "../../metTemplate/output/v8/lepdata_skim_baby.root";
  //if( isData ) file = "../output/V00-00-00/lepdata_skim_baby.root";
  //if( isData ) file = "output/uaf/lepdata_skim_baby.root";
  //if( isData ) file = "output/v8/lepdata_skim_nov4_baby_nov4json.root";
  //if( isData ) file = "output/uaf/lepdata_skim_baby.root";
  //if( isData ) file = "lepdata_skim_Nov4_baby.root";

  cout << "Adding " << file << endl;

  TFile *f = TFile::Open( file );

  TH1F* hee = new TH1F("hee","",75,0,150);
  TH1F* hmm = new TH1F("hmm","",75,0,150);
  TH1F* hem = new TH1F("hem","",75,0,150);

  TCut sel("njets>-1");
 
  TCut ee("leptype==0&&jetptll-ptll>-5&&jetptlt-ptlt>-5&&dilmass>81&&dilmass<101");
  //TCut ee("leptype==0&&dilmass>81&&dilmass<101");
  TCut mm("leptype==1&&nmatchedpfmuons==2&&dilmasspf>81&&dilmasspf<101");
  TCut em("leptype==2&&nmatchedpfmuons==1&&dilmass>81&&dilmass<101");

  //TCut ee("leptype==0&&jetptll-ptll>-5&&jetptlt-ptlt>-5");
  //TCut mm("leptype==1&&nmatchedpfmuons==2");
  //TCut em("leptype==2&&nmatchedpfmuons==1");

  //TCut weight("weight");
  TCut weight("1");

  TCanvas *c1 = new TCanvas();
  c1->cd();

  TPad *plotpad = new TPad("plotpad","plotpad",0.0,0.0,1.0,0.8);
  plotpad->Draw();
  plotpad->cd();

  TTree* T1 = (TTree*) f->Get("T1");

  T1->Draw("TMath::Min(pfmet,149.99)>>hee",(sel+ee)*weight);
  T1->Draw("TMath::Min(pfmet,149.99)>>hmm",(sel+mm)*weight);
  T1->Draw("TMath::Min(pfmet,149.99)>>hem",(sel+em)*weight);

  hee->Sumw2();
  hmm->Sumw2();
  hem->Sumw2();

  hmm->Draw();
  hmm->GetXaxis()->SetTitle("pfmet (GeV)");
  hee->SetLineColor(2);
  hee->SetMarkerColor(2);
  hem->SetLineColor(4);
  hem->SetMarkerColor(4);
  hmm->SetMarkerStyle(20);
  hem->SetMarkerStyle(24);
  hmm->SetMarkerSize(0.7);
  hee->Draw("samehist");
  hem->Draw("same");
  gPad->SetLogy(1);

  TLegend *leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->AddEntry(hee,"ee","l");
  leg->AddEntry(hmm,"#mu#mu");
  leg->AddEntry(hem,"e#mu");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

//   TLatex *t = new TLatex();
//   t->SetNDC();
//   if( isData ){
//     t->DrawLatex(0.4,0.6,  "CMS");
//     //t->DrawLatex(0.4,0.53, "Selected Z+#geq2jet Events (DATA)");
//     t->DrawLatex(0.4,0.53, "Selected Z Events (DATA)");
//     t->DrawLatex(0.4,0.46, "34.0 pb^{-1} at #sqrt{s} = 7 TeV");
//   }else{
//     t->DrawLatex(0.4,0.6, "CMS");
//     t->DrawLatex(0.4,0.53,"Selected Z+0jet Events (Z+jets MC)");
//   }

  c1->cd();

  TPad *pullpad = new TPad("pullpad","pullpad",0.0,0.8,1.0,1.0);
  pullpad->Draw();
  pullpad->cd();

  TH1F* hratio = (TH1F*) hmm->Clone();
  hratio->Divide(hee);
  hratio->Draw();
  gPad->SetGridy();
  hratio->SetMinimum(0.8);
  hratio->SetMaximum(1.6);
  hratio->GetYaxis()->SetLabelSize(0.15);
  hratio->GetYaxis()->SetNdivisions(7);
  hratio->GetYaxis()->SetTitle("#mu#mu / ee  ");
  hratio->GetXaxis()->SetTitle("");
  hratio->GetYaxis()->SetTitleSize(0.2);
  hratio->GetYaxis()->SetTitleOffset(0.2);

  float cut1 = 30;
  float cut2 = 60;
  float cut3 = 120;

  int bin1 = hee->FindBin(cut1);
  int bin2 = hee->FindBin(cut2);
  int bin3 = hee->FindBin(cut3);
  
  cout << "ee tot         " << hee->Integral() << endl;
  cout << "ee met>30  GeV " << hee->Integral(bin1,1000) << endl;
  cout << "ee met>60  GeV " << hee->Integral(bin2,1000) << endl;
  cout << "ee met>120 GeV " << hee->Integral(bin3,1000) << endl;

  cout << "mm tot         " << hmm->Integral() << endl;
  cout << "mm met>30  GeV " << hmm->Integral(bin1,1000) << endl;
  cout << "mm met>60  GeV " << hmm->Integral(bin2,1000) << endl;
  cout << "mm met>120 GeV " << hmm->Integral(bin3,1000) << endl;

  int width1 = 15;
  int width2 =  5;

  cout << "|" << setw(width1) << ""               << setw(width2)
       << "|" << setw(width1) << Form("N(met>%.0f GeV)",cut1)  << setw(width2)
       << "|" << setw(width1) << Form("N(met>%.0f GeV)",cut2)  << setw(width2)
       << "|" << setw(width1) << Form("N(met>%.0f GeV)",cut3)  << setw(width2) 
       << "|" << endl;

  cout << "|" << setw(width1) << "ee"                       << setw(width2)
       << "|" << setw(width1) << hee->Integral(bin1,1000)  << setw(width2)
       << "|" << setw(width1) << hee->Integral(bin2,1000)  << setw(width2)
       << "|" << setw(width1) << hee->Integral(bin3,1000) << setw(width2) 
       << "|" << endl;

  cout << "|" << setw(width1) << "mm"                       << setw(width2)
       << "|" << setw(width1) << hmm->Integral(bin1,1000)  << setw(width2)
       << "|" << setw(width1) << hmm->Integral(bin2,1000)  << setw(width2)
       << "|" << setw(width1) << hmm->Integral(bin3,1000) << setw(width2) 
       << "|" << endl;

  cout << "|" << setw(width1) << "em"                       << setw(width2)
       << "|" << setw(width1) << hem->Integral(bin1,1000)  << setw(width2)
       << "|" << setw(width1) << hem->Integral(bin2,1000)  << setw(width2)
       << "|" << setw(width1) << hem->Integral(bin3,1000) << setw(width2) 
       << "|" << endl;
}



