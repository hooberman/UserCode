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

void makeVertexReweightHistos( bool doPhoton = false ){

  if( doPhoton )  cout << "Using triggers from Photon dataset"         << endl;
  else            cout << "Using triggers from DoubleElectron dataset" << endl;

  char* Zfile1        = (char*) "../output/V00-01-04/data_53X_baby.root";
  char* Zfile2        = (char*) "../output/V00-01-04/data_2012C_53X_baby.root";

  //char* photonfile;
  char* photonbaby;
  char* rootfilename;

  if( doPhoton ){
    //photonfile   = (char*) "../photon_output/V00-00-12/Photon_templates.root";
    photonbaby   = (char*) "../photon_output/V00-00-12/Photon_baby.root";
    rootfilename = (char*) "vtxreweight_Photon_5p1fb.root";
  }

  else{
    //photonfile   = (char*) "../photon_output/V00-01-00/DoubleElectron_templates.root";
    photonbaby   = (char*) "../photon_output/V00-01-00/DoubleElectron_baby_2jets.root";
    rootfilename = (char*) "vtxreweight_DoubleElectron_9p2fb.root";
  }

  cout << "Z files:"      << endl;
  cout << Zfile1          << endl;
  cout << Zfile2          << endl;
  cout << "photon files:" << endl;
  cout << photonbaby      << endl;
  cout << "output file"   << endl;
  cout << rootfilename    << endl;

  //---------------------------
  // Z
  //---------------------------

  TChain *chZ = new TChain("T1");
  chZ->Add(Zfile1);
  chZ->Add(Zfile2);

  TCut Zselection("dilmass>81 && dilmass<101 && njets>=2 && (run<197556 || run>198913)");
  cout << "Using Z selection " << Zselection.GetTitle() << endl;

  TH1F* hZ   = new TH1F("hZ","",50,0,50);
  chZ->Draw("nvtx>>hZ",Zselection);

  //---------------------------
  // photon
  //---------------------------

  TChain* chPhoton = new TChain("T1");
  chPhoton->Add(photonbaby);

  TCut njets("njets>=2");
  TCut etg("etg>=20");
  TCut etag("abs(etag)<2");
  TCut hoe("hoe<0.1");
  TCut pix("photon_pixelseed==0");
  TCut emfrac("jetneutralemfrac>0.7");
  TCut pfjet("jetpt - etg >= -5");
  TCut cjet ("calojetpt - etg >= -5");
  TCut elveto("elveto==0");
  TCut lepveto("maxleppt<20");
  TCut dphi("acos(cos( phig-pfmetphi ) ) > 0.14");
  TCut alltrig;
  TCut runrange("(run<197556 || run>198913)");
  
  if( doPhoton) alltrig = TCut("hlt20>0 || hlt30>0 || hlt50>0 || hlt75>0 || hlt90>0");
  else          alltrig = TCut("hgg22>0 || hgg36>0 || hgg50>0 || hgg75>0 || hgg90>0");

  TCut photonSelection;
  photonSelection += njets;
  photonSelection += etg;
  photonSelection += etag;
  photonSelection += hoe;
  photonSelection += pix;
  photonSelection += emfrac;
  photonSelection += pfjet;
  photonSelection += cjet;
  photonSelection += elveto;
  photonSelection += lepveto;
  photonSelection += dphi;
  photonSelection += alltrig;
  photonSelection += runrange;

  cout << "Using photon selection" << photonSelection.GetTitle() << endl;

  TCut hlt90    , hlt75    , hlt50    , hlt30    , hlt20;
  TCut weight90 , weight75 , weight50 , weight30 , weight20;

  if( doPhoton ){
    hlt90 = TCut("hlt90>0");
    hlt75 = TCut("hlt75>0 && hlt90<1");
    hlt50 = TCut("hlt50>0 && hlt75<1 && hlt90<1");
    hlt30 = TCut("hlt30>0 && hlt50<1 && hlt75<1 && hlt90<1");
    hlt20 = TCut("hlt20>0 && hlt30<1 && hlt50<1 && hlt75<1 && hlt90<1");

    weight90 = TCut("hlt90");
    weight75 = TCut("hlt75");
    weight50 = TCut("hlt50");
    weight30 = TCut("hlt30");
    weight20 = TCut("hlt20");
  }

  else{
    hlt90 = TCut("hgg90>0");
    hlt75 = TCut("hgg75>0 && hgg90<1");
    hlt50 = TCut("hgg50>0 && hgg75<1 && hgg90<1");
    hlt30 = TCut("hgg36>0 && hgg50<1 && hgg75<1 && hgg90<1");
    hlt20 = TCut("hgg22>0 && hgg36<1 && hgg50<1 && hgg75<1 && hgg90<1");

    weight90 = TCut("hgg90");
    weight75 = TCut("hgg75");
    weight50 = TCut("hgg50");
    weight30 = TCut("hgg36");
    weight20 = TCut("hgg22");
  }

  hZ->SetLineColor(2);
  hZ->SetLineWidth(3);
  hZ->Sumw2();

  // photons
  // TFile *f = TFile::Open(photonfile);

  // TH1F* h20  = (TH1F*) f->Get("hnvtxPt20");
  // TH1F* h30  = (TH1F*) f->Get("hnvtxPt30");
  // TH1F* h50  = (TH1F*) f->Get("hnvtxPt50");
  // TH1F* h70  = (TH1F*) f->Get("hnvtxPt70");
  // TH1F* h90  = (TH1F*) f->Get("hnvtxPt90");
  // TH1F* hall = (TH1F*) f->Get("hnvtxAll");

  TH1F* h20  = new TH1F("h20" ,"h20",50,0,50);
  TH1F* h30  = new TH1F("h30" ,"h30",50,0,50);
  TH1F* h50  = new TH1F("h50" ,"h50",50,0,50);
  TH1F* h70  = new TH1F("h70" ,"h70",50,0,50);
  TH1F* h90  = new TH1F("h90" ,"h90",50,0,50);
  TH1F* hall = new TH1F("hall","hall",50,0,50);

  chPhoton->Draw("nvtx>>h20", (photonSelection + hlt20) * weight20);
  chPhoton->Draw("nvtx>>h30", (photonSelection + hlt30) * weight30);
  chPhoton->Draw("nvtx>>h50", (photonSelection + hlt50) * weight50);
  chPhoton->Draw("nvtx>>h70", (photonSelection + hlt75) * weight75);
  chPhoton->Draw("nvtx>>h90", (photonSelection + hlt90) * weight90);

  hall->Add(h20);
  hall->Add(h30);
  hall->Add(h50);
  hall->Add(h70);
  hall->Add(h90);

  cout << "pt20  " << h20->GetEntries()  << " " << h20->Integral() << endl;
  cout << "pt30  " << h30->GetEntries()  << " " << h30->Integral() << endl;
  cout << "pt50  " << h50->GetEntries()  << " " << h50->Integral() << endl;
  cout << "pt70  " << h70->GetEntries()  << " " << h70->Integral() << endl;
  cout << "pt90  " << h90->GetEntries()  << " " << h90->Integral() << endl;
  cout << "ptall " << hall->GetEntries() << " " << hall->Integral() << endl;

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
