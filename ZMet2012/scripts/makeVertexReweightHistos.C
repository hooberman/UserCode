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

void makeVertexReweightHistos(){

  char* Zfile        = (char*) "../output/V00-00-10/data_baby.root";
  char* photonfile   = (char*) "../photon_output/V00-00-07/DoubleElectron_templates.root";
  char* rootfilename = (char*) "vtxreweight_24fb.root";

  // Z
  TChain *chZ = new TChain("T1");
  chZ->Add(Zfile);

  TH1F* hZ   = new TH1F("hZ","",50,0,50);
  chZ->Draw("nvtx>>hZ","dilmass>81 && dilmass<101 && njets>=2");

  hZ->SetLineColor(2);
  hZ->SetLineWidth(3);
  hZ->Sumw2();

  // photons
  TFile *f = TFile::Open(photonfile);

  TH1F* h20  = (TH1F*) f->Get("hnvtxPt20");
  TH1F* h30  = (TH1F*) f->Get("hnvtxPt30");
  TH1F* h50  = (TH1F*) f->Get("hnvtxPt50");
  TH1F* h70  = (TH1F*) f->Get("hnvtxPt70");
  TH1F* h90  = (TH1F*) f->Get("hnvtxPt90");
  TH1F* hall = (TH1F*) f->Get("hnvtxAll");

  h20->Rebin(5);
  h30->Rebin(5);
  h50->Rebin(5);
  h70->Rebin(5);
  h90->Rebin(5);
  hall->Rebin(5);
  hZ->Rebin(5);

  h20-> Scale(1.0 / h20->Integral());
  h30-> Scale(1.0 / h30->Integral());
  h50-> Scale(1.0 / h50->Integral());
  h70-> Scale(1.0 / h70->Integral());
  h90-> Scale(1.0 / h90->Integral());
  hall->Scale(1.0 / hall->Integral());
  hZ->  Scale(1.0 / hZ->Integral());

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.05);

  h20->SetLineColor(2);
  h30->SetLineColor(3);
  h50->SetLineColor(4);
  h70->SetLineColor(5);
  h90->SetLineColor(6);

  hall->SetLineWidth(3);
  
  hall->GetXaxis()->SetRangeUser(0,200);
  hall->GetXaxis()->SetTitle("N_{VTX}");
  hall->Draw("hist");
  h20->Draw("samehist");
  h30->Draw("samehist");
  h50->Draw("samehist");
  h70->Draw("samehist");
  h90->Draw("samehist");
  hall->Draw("samehist");
  hZ->Draw("samehist");


  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(h20,"HLT20","l");
  leg->AddEntry(h30,"HLT30","l");
  leg->AddEntry(h50,"HLT50","l");
  leg->AddEntry(h70,"HLT70","l");
  leg->AddEntry(h90,"HLT90","l");
  leg->AddEntry(hall,"all photons","l");
  leg->AddEntry(hZ  ,"Z","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TH1F* hratio20 = (TH1F*) hZ->Clone("hratio20");
  TH1F* hratio30 = (TH1F*) hZ->Clone("hratio30");
  TH1F* hratio50 = (TH1F*) hZ->Clone("hratio50");
  TH1F* hratio70 = (TH1F*) hZ->Clone("hratio70");
  TH1F* hratio90 = (TH1F*) hZ->Clone("hratio90");

  hratio20->Divide(h20);
  hratio30->Divide(h30);
  hratio50->Divide(h50);
  hratio70->Divide(h70);
  hratio90->Divide(h90);

  TCanvas *c2 = new TCanvas();
  c2->cd();

  hratio20->SetLineColor(2);
  hratio30->SetLineColor(3);
  hratio50->SetLineColor(4);
  hratio70->SetLineColor(5);
  hratio90->SetLineColor(6);

  hratio20->Draw("hist");
  hratio30->Draw("samehist");
  hratio50->Draw("samehist");
  hratio70->Draw("samehist");
  hratio90->Draw("samehist");

  TFile *fratio = TFile::Open(rootfilename,"RECREATE");
  fratio->cd();
  hratio20->Write();
  hratio30->Write();
  hratio50->Write();
  hratio70->Write();
  hratio90->Write();
  fratio->Close();
}
