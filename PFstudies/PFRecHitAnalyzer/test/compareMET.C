#ifndef __CINT__
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TProfile.h"
#include <sstream>
#include <vector>
#endif

void plotHists(TH1* h1, TH1* h2, TH1* h3, TH1* h4, char* leg1, char* leg2, char* leg3, char* leg4, char* title, char* xtitle, int print = 0, int rebin = 1);
void plotHists(TH1* h1, TH1* h2, TH1* h3,          char* leg1, char* leg2, char* leg3,             char* title, char* xtitle, int print = 0, int rebin = 1);
void plotHists(TH1* h1, TH1* h2,                   char* leg1, char* leg2,                         char* title, char* xtitle, int print = 0, int rebin = 1);


int colors[] = { 2,1,4};
float metval = 40;

void compareMET( bool printgif = false){

  //TFile *file   = TFile::Open( "root/qcd_histos.root" );
  //char* filename = "qcd";
  //char* title    = "QCD";

  TFile *file   = TFile::Open( "root/ttbar_matchsumet_histos.root" );
  char* filename = "ttbar";
  char* title    = "t#bar{t}";


  
  char* leg1 = "calo  MET";
  char* leg2 = "PFR MET (nothresh)";
  char* leg3 = "PFR MET";
  
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();

  TH1F *h_calomet           = (TH1F*) file->Get("hcalo_met");
  TH1F *h_pfrmet_nothresh   = (TH1F*) file->Get("hpfr_met_nothresh");
  TH1F *h_pfrmet            = (TH1F*) file->Get("hpfr_met");
  TH1F *h_pfcmet            = (TH1F*) file->Get("hpfc_met");
  
  plotHists(h_calomet, h_pfrmet_nothresh, h_pfrmet, h_pfcmet,
            leg1,leg2,leg3, "PFC MET", title,"MET (GeV)",0,2);
  
  if( printgif) c1->Print(Form( "plots/met_%s.gif" , filename ));

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->cd();

  TH1F *h_calosumet           = (TH1F*) file->Get("hcalo_sumet");
  TH1F *h_pfrsumet_nothresh   = (TH1F*) file->Get("hpfr_sumet_nothresh");
  TH1F *h_pfrsumet            = (TH1F*) file->Get("hpfr_sumet");
  TH1F *h_pfcsumet            = (TH1F*) file->Get("hpfc_sumet");
  
  plotHists(h_calosumet, h_pfrsumet_nothresh, h_pfrsumet, h_pfcsumet, 
            leg1,leg2,leg3, "PFC MET", title,  "sumET (GeV)",0,2);
  
  if( printgif) c2->Print(Form( "plots/sumet_%s.gif" , filename ));


  TCanvas *c3 = new TCanvas("c3a","",1200,800);
  c3->Divide(3,2);

  //EB
  c3->cd(1);

  TH1F *h_caloebsumet           = (TH1F*) file->Get("hcalo_ebsumet");
  TH1F *h_pfrebsumet            = (TH1F*) file->Get("hpfr_ebsumet");
  TH1F *h_pfcebsumet            = (TH1F*) file->Get("hpfc_ebsumet");
  
  plotHists(h_caloebsumet, h_pfrebsumet,h_pfcebsumet, "calo","pfrechits", "pfclusters", title,"EB sumET (GeV)",0,2);
  
  //EE
  c3->cd(4);

  TH1F *h_caloeesumet           = (TH1F*) file->Get("hcalo_eesumet");
  TH1F *h_pfreesumet            = (TH1F*) file->Get("hpfr_eesumet");
  TH1F *h_pfceesumet            = (TH1F*) file->Get("hpfc_eesumet");

  plotHists(h_caloeesumet, h_pfreesumet, h_pfceesumet, "calo","pfrechits","pfclusters", title,"EE sumET (GeV)",0,2);
  
  //HB
  c3->cd(2);

  TH1F *h_calohbsumet           = (TH1F*) file->Get("hcalo_hbsumet");
  TH1F *h_pfrhbsumet            = (TH1F*) file->Get("hpfr_hbsumet");
  TH1F *h_pfchbsumet            = (TH1F*) file->Get("hpfc_hbsumet");

  plotHists(h_calohbsumet, h_pfrhbsumet, h_pfchbsumet, "calo","pfrechits","pfclusters",  title,"HB sumET (GeV)",0,2);
  
  //HE
  c3->cd(5);

  TH1F *h_calohesumet           = (TH1F*) file->Get("hcalo_hesumet");
  TH1F *h_pfrhesumet            = (TH1F*) file->Get("hpfr_hesumet");
  TH1F *h_pfchesumet            = (TH1F*) file->Get("hpfc_hesumet");

  plotHists(h_calohesumet, h_pfrhesumet, h_pfchesumet, "calo","pfrechits","pfclusters", title,"HE sumET (GeV)",0,2);
  
  //HFH
  c3->cd(3);

  TH1F *h_calohfhsumet           = (TH1F*) file->Get("hcalo_hfhsumet");
  TH1F *h_pfrhfhsumet            = (TH1F*) file->Get("hpfr_hfhsumet");
  TH1F *h_pfchfhsumet            = (TH1F*) file->Get("hpfc_hfhsumet");

  plotHists(h_calohfhsumet, h_pfrhfhsumet, h_pfchfhsumet, "calo","pfrechits","pfclusters", title,"HFH sumET (GeV)",0,2);
  
  //HFE
  c3->cd(6);

  TH1F *h_calohfesumet           = (TH1F*) file->Get("hcalo_hfesumet");
  TH1F *h_pfrhfesumet            = (TH1F*) file->Get("hpfr_hfesumet");
  TH1F *h_pfchfesumet            = (TH1F*) file->Get("hpfc_hfesumet");

  plotHists(h_calohfesumet, h_pfrhfesumet, h_pfchfesumet, "calo","pfrechits","pfclusters", title,"HFE sumET (GeV)",0,2);
  
  if( printgif) c3->Print(Form( "plots/sumet_subdets_%s.gif" , filename ));

  TH1F *hebe  = (TH1F*) file->Get("hebe");
  TH1F *heee  = (TH1F*) file->Get("heee");
  TH1F *hhbe  = (TH1F*) file->Get("hhbe");
  TH1F *hhee  = (TH1F*) file->Get("hhee");
  TH1F *hhfhe = (TH1F*) file->Get("hhfhe");
  TH1F *hhfee = (TH1F*) file->Get("hhfee");


  TCanvas *c4=new TCanvas("c4","",1200,800);
  c4->Divide(3,2);
  c4->cd(1);
  gPad->SetLogy(1);
  hebe->Draw();
  c4->cd(4);
  gPad->SetLogy(1);
  heee->Draw();
  c4->cd(2);
  gPad->SetLogy(1);
  hhbe->Draw();
  c4->cd(5);
  gPad->SetLogy(1);
  hhee->Draw();
  c4->cd(3);
  gPad->SetLogy(1);
  hhfhe->Draw();
  c4->cd(6);
  gPad->SetLogy(1);
  hhfee->Draw();

  if( printgif) c4->Print(Form( "plots/pfrechit_et_%s.gif" , filename ));

  TCanvas *c5=new TCanvas("c5","",1200,800);
  c5->cd();
  TProfile *tmet = (TProfile*) file->Get("tmet");
  tmet->Draw();

  if( printgif) c5->Print(Form( "plots/tmet_%s.gif" , filename ));

  //==============met in subdetectors==========
 
  TCanvas *c6 = new TCanvas("c6","",1200,800);
  c6->Divide(3,2);

  //EB
  c6->cd(1);

  //TH1F *h_caloebmet           = (TH1F*) file->Get("hcalo_ebmet");
  TH1F *h_pfrebmet            = (TH1F*) file->Get("hpfr_ebmet");
  TH1F *h_pfcebmet            = (TH1F*) file->Get("hpfc_ebmet");
  
  plotHists(h_pfrebmet,h_pfcebmet, "pfrechits", "pfclusters", title,"EB met (GeV)",0,2);
  
  //EE
  c6->cd(4);

  //TH1F *h_caloeemet           = (TH1F*) file->Get("hcalo_eemet");
  TH1F *h_pfreemet            = (TH1F*) file->Get("hpfr_eemet");
  TH1F *h_pfceemet            = (TH1F*) file->Get("hpfc_eemet");

  plotHists(h_pfreemet, h_pfceemet, "pfrechits","pfclusters", title,"EE met (GeV)",0,2);
  
  //HB
  c6->cd(2);

  //TH1F *h_calohbmet           = (TH1F*) file->Get("hcalo_hbmet");
  TH1F *h_pfrhbmet            = (TH1F*) file->Get("hpfr_hbmet");
  TH1F *h_pfchbmet            = (TH1F*) file->Get("hpfc_hbmet");

  plotHists(h_pfrhbmet, h_pfchbmet, "pfrechits","pfclusters",  title,"HB met (GeV)",0,2);
  
  //HE
  c6->cd(5);

  //TH1F *h_calohemet           = (TH1F*) file->Get("hcalo_hemet");
  TH1F *h_pfrhemet            = (TH1F*) file->Get("hpfr_hemet");
  TH1F *h_pfchemet            = (TH1F*) file->Get("hpfc_hemet");

  plotHists(h_pfrhemet, h_pfchemet, "pfrechits","pfclusters", title,"HE met (GeV)",0,2);
  
  //HFH
  c6->cd(3);

  //TH1F *h_calohfhmet           = (TH1F*) file->Get("hcalo_hfhmet");
  TH1F *h_pfrhfhmet            = (TH1F*) file->Get("hpfr_hfhmet");
  TH1F *h_pfchfhmet            = (TH1F*) file->Get("hpfc_hfhmet");

  plotHists(h_pfrhfhmet, h_pfchfhmet, "pfrechits","pfclusters", title,"HFH met (GeV)",0,2);
  
  //HFE
  c6->cd(6);

  //TH1F *h_calohfemet           = (TH1F*) file->Get("hcalo_hfemet");
  TH1F *h_pfrhfemet            = (TH1F*) file->Get("hpfr_hfemet");
  TH1F *h_pfchfemet            = (TH1F*) file->Get("hpfc_hfemet");

  plotHists(h_pfrhfemet, h_pfchfemet, "pfrechits","pfclusters", title,"HFE met (GeV)",0,2);
 





}


void plotHists(TH1* h1, TH1* h2, TH1* h3, TH1* h4, 
               char* leg1, char* leg2, char* leg3, char* leg4, 
               char* title, char* xtitle, int print, int rebin){

  if(rebin > 1){
    h1->Rebin(rebin);    
    h2->Rebin(rebin);    
    h3->Rebin(rebin);
    h4->Rebin(rebin);
  }
	   
  gPad->SetLogy(1);
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(xtitle);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(4);
  h4->SetLineColor(7);
  h1->SetMarkerColor(1);
  h2->SetMarkerColor(2);
  h3->SetMarkerColor(4);
  h4->SetMarkerColor(7);
  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h3->SetMarkerSize(0);
  h4->SetMarkerSize(0);
  h1->Draw();
//   h2->Draw("samep");
//   h3->Draw("samep");
//   h4->Draw("samep");
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");

  float max=h1->GetMaximum();
  if(h2->GetMaximum()>max)max=h2->GetMaximum();
  if(h3->GetMaximum()>max)max=h3->GetMaximum();
  if(h4->GetMaximum()>max)max=h4->GetMaximum();
  h1->SetMaximum(1.5*max);

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.85);
  leg->AddEntry(h1,leg1);
  leg->AddEntry(h2,leg2);
  leg->AddEntry(h3,leg3);
  leg->AddEntry(h4,leg4);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
  
  if( print > 0 ){
    stringstream s1;
    stringstream s2;
    stringstream s3;
    stringstream s4;
    
    int minbin = h1->FindBin( metval );
    int maxbin = h1->GetXaxis()->GetNbins();

    s1<<"#epsilon(met>" << metval << ") = " << h1->Integral(minbin,maxbin) / h1->Integral() << endl;
    s2<<"#epsilon(met>" << metval << ") = " << h2->Integral(minbin,maxbin) / h2->Integral() << endl;
    s3<<"#epsilon(met>" << metval << ") = " << h3->Integral(minbin,maxbin) / h3->Integral() << endl;
    s4<<"#epsilon(met>" << metval << ") = " << h4->Integral(minbin,maxbin) / h4->Integral() << endl;

    TLatex *l1=new TLatex();
    TLatex *l2=new TLatex();
    TLatex *l3=new TLatex();
    TLatex *l4=new TLatex();

    l1->SetTextSize(0.03);
    l2->SetTextSize(0.03);
    l3->SetTextSize(0.03);
    l4->SetTextSize(0.03);

    l1->SetNDC();
    l2->SetNDC();
    l3->SetNDC();
    l4->SetNDC();

    l1->SetTextColor(1);
    l2->SetTextColor(2);
    l3->SetTextColor(4);
    l4->SetTextColor(7);

    l1->DrawLatex(0.15,0.5, s1.str().c_str());
    l2->DrawLatex(0.15,0.45,s2.str().c_str());
    l3->DrawLatex(0.15,0.4, s3.str().c_str());
    l4->DrawLatex(0.15,0.35,s4.str().c_str());

  }
}

void plotHists(TH1* h1, TH1* h2, TH1* h3, 
               char* leg1, char* leg2, char* leg3, 
               char* title, char* xtitle, int print, int rebin){

  if(rebin > 1){
    h1->Rebin(rebin);
    h2->Rebin(rebin);
    h3->Rebin(rebin);
  }
	 
  gPad->SetLogy(1);
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(xtitle);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(4);
  h1->SetMarkerColor(1);
  h2->SetMarkerColor(2);
  h3->SetMarkerColor(4);
  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0.5);
  h3->SetMarkerSize(0.5);
  h1->Draw();
  h2->Draw("samep");
  h3->Draw("samep");
  
  float max=h1->GetMaximum();
  if(h2->GetMaximum()>max)max=h2->GetMaximum();
  if(h3->GetMaximum()>max)max=h3->GetMaximum();
  h1->SetMaximum(1.5*max);

  TLegend *leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->AddEntry(h1,leg1);
  leg->AddEntry(h2,leg2);
  leg->AddEntry(h3,leg3);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
  
  if( print > 0 ){
    stringstream s1;
    stringstream s2;
    stringstream s3;
    
    int minbin = h1->FindBin( metval );
    int maxbin = h1->GetXaxis()->GetNbins();

    if( print == 1 ){
      //s1<<"#epsilon(met>" << metval << ") = " << h1->Integral(minbin,maxbin) / h1->Integral() << endl;
      //s2<<"#epsilon(met>" << metval << ") = " << h2->Integral(minbin,maxbin) / h2->Integral() << endl;
      //s3<<"#epsilon(met>" << metval << ") = " << h3->Integral(minbin,maxbin) / h3->Integral() << endl;
      s1<<"N(met>" << metval << ") = " << h1->Integral(minbin,maxbin) << endl;
      s2<<"N(met>" << metval << ") = " << h2->Integral(minbin,maxbin) << endl;
      s3<<"N(met>" << metval << ") = " << h3->Integral(minbin,maxbin) << endl;
    }else if ( print == 2 ){
      s1<<"RMS = " << h1->GetRMS(1) << " GeV" << endl;
      s2<<"RMS = " << h2->GetRMS(1) << " GeV" << endl;
      s3<<"RMS = " << h3->GetRMS(1) << " GeV" << endl;
    }

    TLatex *l1=new TLatex();
    TLatex *l2=new TLatex();
    TLatex *l3=new TLatex();

    l1->SetTextSize(0.03);
    l2->SetTextSize(0.03);
    l3->SetTextSize(0.03);

    l1->SetNDC();
    l2->SetNDC();
    l3->SetNDC();

    l1->SetTextColor(1);
    l2->SetTextColor(2);
    l3->SetTextColor(4);

    l1->DrawLatex(0.55,0.5, s1.str().c_str());
    l2->DrawLatex(0.55,0.45,s2.str().c_str());
    l3->DrawLatex(0.55,0.4, s3.str().c_str());
    
  }

}

void plotHists(TH1* h1, TH1* h2, char* leg1, char* leg2, char* title, char* xtitle, int print, int rebin){

  if(rebin > 1){
    h1->Rebin(rebin);
    h2->Rebin(rebin);
  }
	 
  gPad->SetLogy(1);
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(xtitle);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h1->SetMarkerColor(1);
  h2->SetMarkerColor(2);
  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0.5);
  h1->Draw();
  h2->Draw("samep");

  float max=h1->GetMaximum();
  if(h2->GetMaximum()>max)max=h2->GetMaximum();
  h1->SetMaximum(1.5*max);

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.85);
  leg->AddEntry(h1,leg1);
  leg->AddEntry(h2,leg2);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
  
  if( print > 0 ){
    stringstream s1;
    stringstream s2;
    
    int minbin = h1->FindBin( metval );
    int maxbin = h1->GetXaxis()->GetNbins();

    if( print == 1 ){
      //s1<<"#epsilon(met>" << metval << ") = " << h1->Integral(minbin,maxbin) / h1->Integral() << endl;
      //s2<<"#epsilon(met>" << metval << ") = " << h2->Integral(minbin,maxbin) / h2->Integral() << endl;
      s1<<"N(met>" << metval << ") = " << h1->Integral(minbin,maxbin) << endl;
      s2<<"N(met>" << metval << ") = " << h2->Integral(minbin,maxbin) << endl;

    }else if ( print == 2 ){
      s1<<"RMS = " << h1->GetRMS(1) << " GeV" << endl;
      s2<<"RMS = " << h2->GetRMS(1) << " GeV" << endl;
    }

    TLatex *l1=new TLatex();
    TLatex *l2=new TLatex();
    
    l1->SetTextSize(0.03);
    l2->SetTextSize(0.03);
    
    l1->SetNDC();
    l2->SetNDC();
    
    l1->SetTextColor(1);
    l2->SetTextColor(2);
    
    l1->DrawLatex(0.55,0.5, s1.str().c_str());
    l2->DrawLatex(0.55,0.45,s2.str().c_str());
        
  }

}

