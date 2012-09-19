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

float getObservedLimit( int metcut , float seff , bool do3jets );
float getExpectedLimit( int metcut , float seff , bool do3jets );

using namespace std;

void SMS(char* sample , bool print = false){

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------
  
  const float denom    = 20000;
  const float lumi     = 4980;
  const char* filename = Form("../../output/V00-02-04/%s_baby.root",sample);
  bool do3jets = false;

  cout << "Using file        " << filename << endl;
  cout << "Using denominator " << denom    << " events" << endl;
  cout << "Using lumi        " << lumi     << " pb-1" << endl;
  cout << "Doing 3 jets?     " << do3jets  << endl;

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("T1");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("gluino");

  //--------------------------------------------------
  // preselection
  //--------------------------------------------------

  TCut zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets >= 2");
  TCut njets3("njets >= 3");
  TCut sf("leptype==0||leptype==1");
  TCut met30 ("pfmet>30");
  TCut met60 ("pfmet>60");
  TCut met100("pfmet>100");
  TCut met200("pfmet>200");
  TCut met300("pfmet>300");

  TCut presel  = zmass + njets2 + sf;
  if( do3jets ) presel  = zmass + njets3 + sf;

  cout << "Using selection   " << presel.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<int>     cuts;

  sigcuts.push_back(TCut(presel+met100));     signames.push_back("E_{T}^{miss} > 100 GeV");     labels.push_back("met100"); cuts.push_back(100);
  sigcuts.push_back(TCut(presel+met200));     signames.push_back("E_{T}^{miss} > 200 GeV");     labels.push_back("met200"); cuts.push_back(200);
  sigcuts.push_back(TCut(presel+met300));     signames.push_back("E_{T}^{miss} > 300 GeV");     labels.push_back("met300"); cuts.push_back(300);

  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* heff[nsig];
  TH2F* heffup[nsig];
  TH2F* heffdn[nsig];
  TH2F* hxsec[nsig];
  TH2F* hxsec_exp[nsig];
  TH2F* hexcl[nsig];
  TH2F* hjes[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    TString jesup(sigcuts.at(i));
    jesup.ReplaceAll("njets","njetsup");
    jesup.ReplaceAll("pfmet","pfmetup");

    TString jesdn(sigcuts.at(i));
    jesdn.ReplaceAll("njets","njetsdn");
    jesdn.ReplaceAll("pfmet","pfmetdn");

    TCut jesupcut(jesup);
    TCut jesdncut(jesdn);

    cout << endl << endl;
    cout << "Signal region : " << labels.at(i)  << endl;
    cout << "Selection     : " << sigcuts.at(i) << endl;
    cout << "Selection up  : " << jesupcut      << endl;
    cout << "Selection dn  : " << jesdncut      << endl;

    heff[i]      = new TH2F(Form("heff_%i",i)        , Form("heff_%i",i)       , 48,0,1200,48,0,1200);
    heffup[i]    = new TH2F(Form("heffup_%i",i)      , Form("heffup_%i",i)     , 48,0,1200,48,0,1200);
    heffdn[i]    = new TH2F(Form("heffdn_%i",i)      , Form("heffdn_%i",i)     , 48,0,1200,48,0,1200);
    hxsec[i]     = new TH2F(Form("hxsec_%i",i)       , Form("hxsec_%i",i)      , 48,0,1200,48,0,1200);
    hxsec_exp[i] = new TH2F(Form("hxsec_exp_%i",i)   , Form("hxsec_exp_%i",i)  , 48,0,1200,48,0,1200);
    hexcl[i]     = new TH2F(Form("hexcl_%i",i)       , Form("hexcl_%i",i)      , 48,0,1200,48,0,1200);
    hjes[i]      = new TH2F(Form("hjes_%i",i)        , Form("hjes_%i",i)       , 48,0,1200,48,0,1200);

    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i));
    heff[i]->Scale(0.95/denom);

    ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut);
    heffup[i]->Scale(0.95/denom);

    ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut);
    heffdn[i]->Scale(0.95/denom);

    for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin)-12.5;
	float ml = heff[i]->GetYaxis()->GetBinCenter(jbin)-12.5;

	float eff    = heff[i]->GetBinContent(ibin,jbin);
	float effup  = heffup[i]->GetBinContent(ibin,jbin);
	float effdn  = heffdn[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	float dup    = effup/eff-1;
	float ddn    = 1-effdn/eff;
	float djes   = 0.5 * (dup+ddn);
	hjes[i]->SetBinContent(ibin,jbin,djes);

	float toterr = sqrt( 0.022*0.022 + 0.05*0.05 + djes*djes );

	float this_ul = getObservedLimit( cuts.at(i) , toterr , do3jets );
	float xsecul  = this_ul / ( lumi * eff * 0.19 );

	float this_ul_exp = getExpectedLimit( cuts.at(i) , toterr , do3jets );
	float xsecul_exp  = this_ul_exp / ( lumi * eff * 0.19 );

	if( eff > 0 ){
	  hxsec[i]->SetBinContent(ibin,jbin, xsecul );
	  hxsec_exp[i]->SetBinContent(ibin,jbin, xsecul_exp );
	}

	int   bin = refxsec->FindBin(mg);
	float xsec = refxsec->GetBinContent(bin);

	hexcl[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul )   hexcl[i]->SetBinContent(ibin,jbin,1);
	//cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
      }
    }
  }

  delete ctemp;

  cout << endl << endl;

  //--------------------------------------------------
  // make pretty pictures
  //--------------------------------------------------
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);

  TCanvas* can[nsig];

  for( unsigned int i = 0 ; i < nsig ; ++i ){
  
    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    //can[i]->Divide(2,1);
    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1800,600);
    //can[i]->Divide(3,1);
    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,1200);
    can[i]->Divide(2,2);

    //-------------------------------
    // efficiency
    //-------------------------------

    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    heff[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    heff[i]->GetZaxis()->SetTitle("efficiency");
    heff[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    //-------------------------------
    // cross section
    //-------------------------------
  
    can[i]->cd(2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    gPad->SetLogz();
    hxsec[i]->GetXaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hxsec[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hxsec[i]->GetZaxis()->SetTitle("#sigma upper limit");
    hxsec[i]->Draw("colz");
    hxsec[i]->SetMinimum(0.01);
    hxsec[i]->SetMaximum(10);

    TGraph* gr_excl      = getRefXsecGraph(hxsec[i], "T5zz", 1.0);
    TGraph* gr_excl_down = getRefXsecGraph(hxsec[i], "T5zz", 1./3.);
    TGraph* gr_excl_up   = getRefXsecGraph(hxsec[i], "T5zz", 3.);

    gr_excl->SetLineWidth(2);
    gr_excl_up->SetLineWidth(2);
    gr_excl_down->SetLineWidth(2);
    gr_excl_up->SetLineStyle(2);
    gr_excl_down->SetLineStyle(3);
    gr_excl->Draw("same");
    gr_excl_up->Draw("same");
    gr_excl_down->Draw("same");

    TLegend *leg = new TLegend(0.2,0.53,0.53,0.67);
    leg->AddEntry(gr_excl,"#sigma^{prod} = #sigma^{NLO-QCD}","l");
    leg->AddEntry(gr_excl_up,"#sigma^{prod} = 3 #times #sigma^{NLO-QCD}","l");
    leg->AddEntry(gr_excl_down,"#sigma^{prod} = 1/3 #times #sigma^{NLO-QCD}","l");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    //-------------------------------
    // excluded points
    //-------------------------------

    can[i]->cd(3);    
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hexcl[i]->GetZaxis()->SetTitle("excluded points");
    hexcl[i]->Draw("colz");
    //gr_excl->Draw("same");

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    //-------------------------------
    // JES uncertainty
    //-------------------------------

    can[i]->cd(4);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hjes[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hjes[i]->GetXaxis()->SetTitle("gluino mass (GeV)");
    hjes[i]->GetZaxis()->SetTitle("JES uncertainty");
    hjes[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    if( print ){
      can[i]->Print(Form("../../plots/%s.pdf",labels.at(i).c_str()));
      //can[i]->Print(Form("../plots/%s.eps",labels.at(i).c_str()));
      //gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.eps  ../plots/%s.pdf",labels.at(i).c_str(),labels.at(i).c_str()));
    }

    int bin = heff[i]->FindBin(600,200);
    cout << "efficiency (600,200) " << heff[i]->GetBinContent(bin) << endl;
    bin = heff[i]->FindBin(725,425);
    cout << "efficiency (600,200) " << heff[i]->GetBinContent(bin) << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin) << endl;
    cout << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    cout << "JES                  " << hjes[i]->GetBinContent(bin) << endl;
    cout << "tot err              " << sqrt(pow(hjes[i]->GetBinContent(bin),2)+0.06*0.06+0.05*0.05) << endl;
    cout << endl << endl;
  }
  
  TFile *outfile = TFile::Open(Form("%s_histos.root",sample),"RECREATE");
  outfile->cd();
  for( unsigned int i = 0 ; i < nsig ; ++i ){
    hxsec[i]->Write();
    heff[i]->Write();
    hxsec_exp[i]->Write();
  }
  outfile->Close();

}



float getObservedLimit( int metcut , float seff , bool do3jets ){

  float ul = 999;

  if( do3jets ){
    cout << "ERROR!!!! DOING 3 JET!!!" << endl;
    exit(0);

    if( metcut == 100 ){
      if(seff >= 0.0 && seff < 0.1) ul = 42.2;
      if(seff >= 0.1 && seff < 0.2) ul = 45.3;
      if(seff >= 0.2 && seff < 0.3) ul = 48.6;
      if(seff >= 0.3 && seff < 0.4) ul = 53.1;
      if(seff >= 0.4 && seff < 0.5) ul = 57.4;
      if(seff >= 0.5 && seff < 0.6) ul = 63.0;
      if(seff >= 0.6 && seff < 0.7) ul = 68.1;
      if(seff >= 0.7 && seff < 0.8) ul = 73.8;
      if(seff >= 0.8 && seff < 0.9) ul = 80.7;
    }
    else if( metcut == 200 ){
      if(seff >= 0.0 && seff < 0.1) ul = 6.7;
      if(seff >= 0.1 && seff < 0.2) ul = 7.3;
      if(seff >= 0.2 && seff < 0.3) ul = 7.5;
      if(seff >= 0.3 && seff < 0.4) ul = 7.5;
      if(seff >= 0.4 && seff < 0.5) ul = 7.8;
      if(seff >= 0.5 && seff < 0.6) ul = 8.6;
      if(seff >= 0.6 && seff < 0.7) ul = 9.0;
      if(seff >= 0.7 && seff < 0.8) ul = 9.9;
      if(seff >= 0.8 && seff < 0.9) ul = 10.6;
    }
    else if( metcut == 300 ){
      if(seff >= 0.0 && seff < 0.1) ul = 2.8;
      if(seff >= 0.1 && seff < 0.2) ul = 2.7;
      if(seff >= 0.2 && seff < 0.3) ul = 2.8;
      if(seff >= 0.3 && seff < 0.4) ul = 3.0;
      if(seff >= 0.4 && seff < 0.5) ul = 3.2;
      if(seff >= 0.5 && seff < 0.6) ul = 3.4;
      if(seff >= 0.6 && seff < 0.7) ul = 3.7;
      if(seff >= 0.7 && seff < 0.8) ul = 3.8;
      if(seff >= 0.8 && seff < 0.9) ul = 4.2;
    }  
    else{
      cout << "ERROR! unrecognized met cut " << metcut << ", quitting" << endl;
      exit(0);
    }
  }
  else{
    if( metcut == 100 ){
      if(seff >= 0.0 && seff < 0.1) ul = 57.6;
      if(seff >= 0.1 && seff < 0.2) ul = 59.9;
      if(seff >= 0.2 && seff < 0.3) ul = 64.5;
      if(seff >= 0.3 && seff < 0.4) ul = 69.0;
      if(seff >= 0.4 && seff < 0.5) ul = 74.5;
      if(seff >= 0.5 && seff < 0.6) ul = 79.9;
      if(seff >= 0.6 && seff < 0.7) ul = 86.1;
      if(seff >= 0.7 && seff < 0.8) ul = 92.1;
      if(seff >= 0.8 && seff < 0.9) ul = 101.2;
    }
    else if( metcut == 200 ){
      if(seff >= 0.0 && seff < 0.1) ul = 8.3;
      if(seff >= 0.1 && seff < 0.2) ul = 8.9;
      if(seff >= 0.2 && seff < 0.3) ul = 8.8;
      if(seff >= 0.3 && seff < 0.4) ul = 9.5;
      if(seff >= 0.4 && seff < 0.5) ul = 10.1;
      if(seff >= 0.5 && seff < 0.6) ul = 10.7;
      if(seff >= 0.6 && seff < 0.7) ul = 11.2;
      if(seff >= 0.7 && seff < 0.8) ul = 11.6;
      if(seff >= 0.8 && seff < 0.9) ul = 12.3;
    }
    else if( metcut == 300 ){
      if(seff >= 0.0 && seff < 0.1) ul = 2.8;
      if(seff >= 0.1 && seff < 0.2) ul = 2.8;
      if(seff >= 0.2 && seff < 0.3) ul = 2.8;
      if(seff >= 0.3 && seff < 0.4) ul = 3.0;
      if(seff >= 0.4 && seff < 0.5) ul = 3.1;
      if(seff >= 0.5 && seff < 0.6) ul = 3.3;
      if(seff >= 0.6 && seff < 0.7) ul = 3.5;
      if(seff >= 0.7 && seff < 0.8) ul = 3.6;
      if(seff >= 0.8 && seff < 0.9) ul = 3.6;
    }  
  }

  if( ul > 998 ){
    cout << "Error ul " << ul << " metcut " << metcut << " SEFF " << seff << endl;
  }

  return ul;
}



float getExpectedLimit( int metcut , float seff , bool do3jets ){

  float ul = 999;

  if( do3jets ){
    if( metcut == 100 ){
      cout << "ERROR DOING 3 JETS!!!" << endl;
      if(seff >= 0.0 && seff < 0.1) ul = 36.0;
      if(seff >= 0.1 && seff < 0.2) ul = 38.4;
      if(seff >= 0.2 && seff < 0.3) ul = 40.6;
      if(seff >= 0.3 && seff < 0.4) ul = 44.6;
      if(seff >= 0.4 && seff < 0.5) ul = 48.4;
      if(seff >= 0.5 && seff < 0.6) ul = 51.2;
      if(seff >= 0.6 && seff < 0.7) ul = 57.3;
      if(seff >= 0.7 && seff < 0.8) ul = 61.0;
      if(seff >= 0.8 && seff < 0.9) ul = 66.7;
    }
    else if( metcut == 200 ){
      if(seff >= 0.0 && seff < 0.1) ul = 8.3;
      if(seff >= 0.1 && seff < 0.2) ul = 9.3;
      if(seff >= 0.2 && seff < 0.3) ul = 9.2;
      if(seff >= 0.3 && seff < 0.4) ul = 10.1;
      if(seff >= 0.4 && seff < 0.5) ul = 12.9;
      if(seff >= 0.5 && seff < 0.6) ul = 12.2;
      if(seff >= 0.6 && seff < 0.7) ul = 13.9;
      if(seff >= 0.7 && seff < 0.8) ul = 14.9;
      if(seff >= 0.8 && seff < 0.9) ul = 18.0;
    }
    else if( metcut == 300 ){
      if(seff >= 0.0 && seff < 0.1) ul = 4.0;
      if(seff >= 0.1 && seff < 0.2) ul = 4.2;
      if(seff >= 0.2 && seff < 0.3) ul = 4.3;
      if(seff >= 0.3 && seff < 0.4) ul = 4.5;
      if(seff >= 0.4 && seff < 0.5) ul = 4.8;
      if(seff >= 0.5 && seff < 0.6) ul = 5.1;
      if(seff >= 0.6 && seff < 0.7) ul = 5.3;
      if(seff >= 0.7 && seff < 0.8) ul = 5.6;
      if(seff >= 0.8 && seff < 0.9) ul = 5.9;
    }  
    else{
      cout << "ERROR! unrecognized met cut " << metcut << ", quitting" << endl;
      exit(0);
    }
  }
  else{
    if( metcut == 100 ){
      if(seff >= 0.0 && seff < 0.1) ul = 61.3;
      if(seff >= 0.1 && seff < 0.2) ul = 63.8;
      if(seff >= 0.2 && seff < 0.3) ul = 68.1;
      if(seff >= 0.3 && seff < 0.4) ul = 73.6;
      if(seff >= 0.4 && seff < 0.5) ul = 79.5;
      if(seff >= 0.5 && seff < 0.6) ul = 86.3;
      if(seff >= 0.6 && seff < 0.7) ul = 93.3;
      if(seff >= 0.7 && seff < 0.8) ul = 100.3;
      if(seff >= 0.8 && seff < 0.9) ul = 109.5;
    }
    else if( metcut == 200 ){
      if(seff >= 0.0 && seff < 0.1) ul = 11.4;
      if(seff >= 0.1 && seff < 0.2) ul = 11.8;
      if(seff >= 0.2 && seff < 0.3) ul = 12.1;
      if(seff >= 0.3 && seff < 0.4) ul = 13.5;
      if(seff >= 0.4 && seff < 0.5) ul = 13.7;
      if(seff >= 0.5 && seff < 0.6) ul = 14.8;
      if(seff >= 0.6 && seff < 0.7) ul = 16.1;
      if(seff >= 0.7 && seff < 0.8) ul = 17.2;
      if(seff >= 0.8 && seff < 0.9) ul = 18.2;
    }
    else if( metcut == 300 ){
      if(seff >= 0.0 && seff < 0.1) ul = 4.8;
      if(seff >= 0.1 && seff < 0.2) ul = 4.9;
      if(seff >= 0.2 && seff < 0.3) ul = 5.1;
      if(seff >= 0.3 && seff < 0.4) ul = 5.3;
      if(seff >= 0.4 && seff < 0.5) ul = 5.7;
      if(seff >= 0.5 && seff < 0.6) ul = 5.9;
      if(seff >= 0.6 && seff < 0.7) ul = 6.3;
      if(seff >= 0.7 && seff < 0.8) ul = 6.8;
      if(seff >= 0.8 && seff < 0.9) ul = 7.0;
    }  
    else{
      cout << "ERROR! unrecognized met cut " << metcut << ", quitting" << endl;
      exit(0);
    }
  }

  if( ul > 998 ){
    cout << "Error ul " << ul << " metcut " << metcut << " SEFF " << seff << endl;
  }

  return ul;


}




