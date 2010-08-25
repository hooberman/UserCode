#include "TROOT.h"
#include "TH1.h"
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
#include "TLine.h"

using namespace std;



void drawPlots( bool printgif = false ){
  
  TFile *f = TFile::Open("output/data_e_ttbarV1align_histos.root");

  //electron/muon eta
  TH1F * heleta_metlt20 = (TH1F*) f->Get("heleta_metlt20");
  TH1F * heleta_metgt30 = (TH1F*) f->Get("heleta_metgt30");
  TH1F * hmueta_metlt20 = (TH1F*) f->Get("hmueta_metlt20");
  TH1F * hmueta_metgt30 = (TH1F*) f->Get("hmueta_metgt30");

  TH1F * hdphijetmet_metgt30_ee = (TH1F*) f->Get("hdphijetmet_metgt30_ee");
  TH1F * hdphijetmet_metlt20_ee = (TH1F*) f->Get("hdphijetmet_metlt20_ee");
  TH1F * hdphijetmet_metgt30_mm = (TH1F*) f->Get("hdphijetmet_metgt30_mm");
  TH1F * hdphijetmet_metlt20_mm = (TH1F*) f->Get("hdphijetmet_metlt20_mm");


  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);

  TLegend *leg = new TLegend(0.7,0.7,0.95,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  TLine line;
  line.SetLineColor(6);

  c1->cd(1);
  heleta_metlt20->Rebin(10);
  heleta_metlt20->SetTitle("");
  heleta_metgt30->Rebin(10);
  heleta_metlt20->SetLineColor(4);
  heleta_metlt20->Draw("hist");
  heleta_metgt30->SetLineColor(2);
  heleta_metgt30->SetMarkerColor(2);
  heleta_metlt20->SetMarkerSize(0);
  heleta_metgt30->Draw("sameE1");
  heleta_metlt20->GetXaxis()->SetRangeUser(-3,3);
  leg->AddEntry(heleta_metlt20,"tcmet < 20 GeV");
  leg->AddEntry(heleta_metgt30,"tcmet > 30 GeV");
  line.DrawLine(-1.479,0,-1.479,1.05*heleta_metlt20->GetMaximum());
  line.DrawLine(1.479, 0, 1.479,1.05*heleta_metlt20->GetMaximum());
  leg->Draw();

  c1->cd(2);
  hmueta_metlt20->Rebin(10);
  hmueta_metlt20->SetTitle("");
  hmueta_metgt30->Rebin(10);
  hmueta_metlt20->SetLineColor(4);
  hmueta_metlt20->Draw("hist");
  hmueta_metgt30->SetLineColor(2);
  hmueta_metgt30->SetMarkerColor(2);
  hmueta_metgt30->Draw("sameE1");
  hmueta_metlt20->GetXaxis()->SetRangeUser(-3,3);
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  hdphijetmet_metlt20_ee->Rebin(5);
  hdphijetmet_metlt20_ee->SetTitle("ee");
  hdphijetmet_metlt20_ee->GetXaxis()->SetTitle("#Delta#phi(jet,tcmet)");
  hdphijetmet_metgt30_ee->Rebin(5);
  hdphijetmet_metlt20_ee->SetLineColor(4);
  hdphijetmet_metlt20_ee->Draw("hist");
  hdphijetmet_metgt30_ee->SetLineColor(2);
  hdphijetmet_metgt30_ee->SetMarkerColor(2);
  hdphijetmet_metlt20_ee->SetMarkerSize(0);
  hdphijetmet_metgt30_ee->Draw("sameE1");
  hdphijetmet_metlt20_ee->GetXaxis()->SetRangeUser(-3,3);
  leg->Draw();

  c2->cd(2);
  hdphijetmet_metlt20_mm->Rebin(5);
  hdphijetmet_metlt20_mm->SetTitle("#mu#mu");
  hdphijetmet_metlt20_mm->GetXaxis()->SetTitle("#Delta#phi(jet,tcmet)");
  hdphijetmet_metgt30_mm->Rebin(5);
  hdphijetmet_metlt20_mm->SetLineColor(4);
  hdphijetmet_metlt20_mm->Draw("hist");
  hdphijetmet_metgt30_mm->SetLineColor(2);
  hdphijetmet_metgt30_mm->SetMarkerColor(2);
  hdphijetmet_metlt20_mm->SetMarkerSize(0);
  hdphijetmet_metgt30_mm->Draw("sameE1");
  hdphijetmet_metlt20_mm->GetXaxis()->SetRangeUser(-3,3);
  leg->Draw();

  if( printgif ){

    c1->Modified();
    c1->Update();
    c1->Print("plots/lepeta.gif");

    c2->Modified();
    c2->Update();
    c2->Print("plots/dphijetmet.gif");
  }

}
