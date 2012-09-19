#include "Utils/SMS_utils.C"
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
#include <sstream>
#include <iomanip>

using namespace std;

void combinePlots(char* sample , bool print = false){
 
  //TFile *file = TFile::Open("T5zz_histos.root");
  TFile *file = TFile::Open(Form("%s_histos.root",sample));
  
  cout << "Opening " << Form("%s_histos.root",sample) << endl;

  const unsigned int ncuts = 3;

  TH2F* hxsec[ncuts];
  TH2F* hxsec_exp[ncuts];
  TH2F* hxsec_best;
  TH2F* hbest;
  TH2F* hxsec_best_exp;

  for( unsigned int i = 0 ; i < ncuts ; ++i ){
    hxsec[i]     = (TH2F*) file->Get(Form("hxsec_%i",i));
    hxsec_exp[i] = (TH2F*) file->Get(Form("hxsec_exp_%i",i));

    if( i == 0 ){
      hxsec_best = (TH2F*) hxsec[i]->Clone("hxsec_best");
      hxsec_best->Reset();
      hbest = (TH2F*) hxsec[i]->Clone("hbest");
      hbest->Reset();
      hxsec_best_exp = (TH2F*) hxsec[i]->Clone("hxsec_best_exp");
      hxsec_best_exp->Reset();
    }


  }

  for( int xbin = 1 ; xbin <= hxsec_best->GetXaxis()->GetNbins() ; ++xbin ){
    for( int ybin = 1 ; ybin <= hxsec_best->GetYaxis()->GetNbins() ; ++ybin ){

      cout << "xbin ybin " << xbin << " " << ybin << endl;
     
      if( hxsec_exp[0]->GetBinContent(xbin,ybin) < 1e-10 ) continue;
      
      int   best_ul    = 0;
      float min_exp_ul = hxsec_exp[0]->GetBinContent(xbin,ybin);

      cout << "exp0 " << hxsec_exp[0]->GetBinContent(xbin,ybin) << endl;

      for( unsigned int i = 1 ; i < ncuts ; ++i ){

	cout << "exp" << i << " " << hxsec_exp[i]->GetBinContent(xbin,ybin) << endl;

	if( hxsec_exp[i]->GetBinContent(xbin,ybin) < min_exp_ul ){
	  min_exp_ul = hxsec_exp[i]->GetBinContent(xbin,ybin);
	  best_ul    = i;
	}

      }

      hxsec_best->SetBinContent(xbin,ybin,hxsec[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_exp->SetBinContent(xbin,ybin,hxsec_exp[best_ul]->GetBinContent(xbin,ybin));
      hbest->SetBinContent(xbin,ybin,best_ul+1);

      cout << "best ul " << best_ul << " " << min_exp_ul << endl;
      cout << "obs limit " << hxsec[best_ul]->GetBinContent(xbin,ybin) << endl << endl << endl;
      
 
    }
  }

  TCanvas *can = new TCanvas("can","",1500,500);
  can->cd();
  can->Divide(3,1);

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // efficiency
  //-------------------------------
  
  can->cd(1);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  TH2F* heff = (TH2F*) file->Get("heff_0");
  heff->GetXaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetLabelSize(0.035);
  heff->GetZaxis()->SetLabelSize(0.035);
  heff->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  heff->GetXaxis()->SetTitle("gluino mass (GeV)");
  heff->GetZaxis()->SetTitle("efficiency");
  heff->Draw("colz");
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
  t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
  t->DrawLatex(0.2,0.71,"E_{T}^{miss} > 100 GeV");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

  //-------------------------------
  // cross section limit
  //-------------------------------

  can->cd(2);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  hxsec_best->GetXaxis()->SetLabelSize(0.035);
  hxsec_best->GetYaxis()->SetLabelSize(0.035);
  hxsec_best->GetZaxis()->SetLabelSize(0.035);
  hxsec_best->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  hxsec_best->GetXaxis()->SetTitle("gluino mass (GeV)");
  hxsec_best->GetZaxis()->SetTitle("#sigma upper limit");
  hxsec_best->Draw("colz");
  hxsec_best->SetMinimum(0.01);
  hxsec_best->SetMaximum(10);
  
  TGraph* gr_excl      = getRefXsecGraph(hxsec_best, "T5zz", 1.0);
  TGraph* gr_excl_down = getRefXsecGraph(hxsec_best, "T5zz", 1./3.);
  TGraph* gr_excl_up   = getRefXsecGraph(hxsec_best, "T5zz", 3.);
  
  gr_excl->SetLineWidth(1);
  gr_excl_up->SetLineWidth(1);
  gr_excl_down->SetLineWidth(1);
  gr_excl_up->SetLineStyle(2);
  gr_excl_down->SetLineStyle(3);
  gr_excl->Draw("same");
  gr_excl_up->Draw("same");
  gr_excl_down->Draw("same");
  
  TLegend *leg = new TLegend(0.2,0.53,0.53,0.67);
  leg->AddEntry(gr_excl,     "#sigma^{NLO-QCD}","l");
  leg->AddEntry(gr_excl_up,  "3 #times #sigma^{NLO-QCD}","l");
  leg->AddEntry(gr_excl_down,"1/3 #times #sigma^{NLO-QCD}","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
  t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
  t->DrawLatex(0.2,0.71,"best limit");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");


  //-------------------------------
  // best signal region
  //-------------------------------
  
  can->cd(3);
  hbest->Draw("colz");
  hbest->Draw("sametext");
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  hbest->SetMaximum(4);
  hbest->GetXaxis()->SetLabelSize(0.035);
  hbest->GetYaxis()->SetLabelSize(0.035);
  hbest->GetZaxis()->SetLabelSize(0.035);
  hbest->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  hbest->GetXaxis()->SetTitle("gluino mass (GeV)");
  hbest->GetZaxis()->SetTitle("best signal region");

  t->SetTextSize(0.04);
  t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
  t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary      #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");

  t->DrawLatex(0.2,0.60,"1 = E_{T}^{miss} > 100 GeV");
  t->DrawLatex(0.2,0.55,"2 = E_{T}^{miss} > 200 GeV");
  t->DrawLatex(0.2,0.50,"3 = E_{T}^{miss} > 300 GeV");

  if( print ){
    can->Print("../../plots/SMS.eps");
    can->Print("../../plots/SMS.pdf");
    can->Print("../../plots/SMS.png");
  }

  //-------------------------------
  // excluded region
  //-------------------------------

  TFile *xsecfile = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("gluino");

  TH2F* hexcluded   = new TH2F("hexcluded","hexcluded", 48,0,1200,48,0,1200);
  TH2F* hexcluded13 = new TH2F("hexcluded13","hexcluded13", 48,0,1200,48,0,1200);
  TH2F* hexcluded3  = new TH2F("hexcluded3","hexcluded3", 48,0,1200,48,0,1200);
  
  for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

      float xsecul = hxsec_best->GetBinContent(ibin,jbin);

      if( xsecul < 1.e-10 ) continue;

      float mg = hexcluded->GetXaxis()->GetBinCenter(ibin)-12.5;
      float ml = hexcluded->GetYaxis()->GetBinCenter(jbin)-12.5;

      int   bin  = refxsec->FindBin(mg);
      float xsec = refxsec->GetBinContent(bin);

      hexcluded->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul )   hexcluded->SetBinContent(ibin,jbin,1);

      hexcluded3->SetBinContent(ibin,jbin,0);
      if( 3 * xsec > xsecul )   hexcluded3->SetBinContent(ibin,jbin,1);

      hexcluded13->SetBinContent(ibin,jbin,0);
      if( (1./3.) * xsec > xsecul )   hexcluded13->SetBinContent(ibin,jbin,1);

      //cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
    }
  }

  TCanvas *c2 = new TCanvas("c2","c2",1500,500);
  c2->Divide(3,1);

  t->SetTextSize(0.07);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded13->Draw("colz");
  hexcluded13->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded13->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded13->Draw("colz");
  gr_excl_down->Draw();
  t->DrawLatex(0.3,0.8,"1/3 #times #sigma^{NLO-QCD}");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded->Draw("colz");
  hexcluded->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded->Draw("colz");
  gr_excl->Draw();
  t->DrawLatex(0.3,0.8,"#sigma^{NLO-QCD}");

  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded3->Draw("colz");
  hexcluded3->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded3->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded3->Draw("colz");
  gr_excl_up->Draw();
  t->DrawLatex(0.3,0.8,"3 #times #sigma^{NLO-QCD}");


  TFile* fout = TFile::Open(Form("%s_combinePlots.root",sample),"RECREATE");
  fout->cd();
  hxsec_best->Write();
  hxsec_best_exp->Write();
  gr_excl->Write();
  fout->Close();


  // TCanvas *can = new TCanvas("can","can",1800,1200);
  // can->cd();
  // can->Divide(3,2);

  // can->cd(i+1);
  // gPad->SetRightMargin(0.2);
  // hxsec[i]->Draw("colz");
  // can->cd(i+4);
  // gPad->SetRightMargin(0.2);
  // hxsec_exp[i]->SetMaximum(10);
  // hxsec_exp[i]->Draw("colz");

}
