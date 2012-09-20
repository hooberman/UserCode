//#include "Utils/SMS_utils.C"
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

float getObservedUpperLimit( int metcut );
float getExpectedUpperLimit( int metcut );

using namespace std;

void SMS(char* sample = "wzsms" , bool print = false){

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------
  
  const float denom    = 100000;
  const float lumi     =   9200;
  const char* filename = Form("../../output/V00-01-05/%s_baby_oldIso.root",sample);

  cout << "Using file        " << filename << endl;
  cout << "Using denominator " << denom    << " events" << endl;
  cout << "Using lumi        " << lumi     << " pb-1" << endl;

  //--------------------------------------------------
  // set text stuff
  //--------------------------------------------------

  char* xtitle  = "m_{#chi_{2}^{0}} = m_{#chi_{1}^{#pm}} [GeV]";
  char* ytitle  = "m_{#chi_{1}^{0}} [GeV]";
  char* label   = "CMS Preliminary #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 9.2 fb^{-1}";
  char* process = "pp #rightarrow#chi_{2}^{0} #chi_{1}^{#pm} #rightarrow WZ+E_{T}^{miss}";

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("T1");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("C1N2_8TeV_finer.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("C1N2_8TeV_NLO");

  //--------------------------------------------------
  // preselection
  //--------------------------------------------------

  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut Zmass("dilmass>81 && dilmass<101");
  TCut sf = "leptype<2";
  TCut njets2("njets>=2");
  TCut bveto("nbcsvm==0");
  TCut nlep2("nlep==2");
  TCut mjj("mjj>70.0 && mjj<110.0");

  TCut presel;
  presel += pt2020;
  presel += Zmass;
  presel += sf;
  presel += njets2;
  presel += bveto;
  presel += mjj;
  presel += nlep2;

  TCut weight("vtxweight * trgeff");

  cout << "Using selection   " << presel.GetTitle() << endl;
  cout << "Using weight      " << weight.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<int>     cuts;

  TCut met60 ("pfmet > 60.0");
  TCut met80 ("pfmet > 80.0");
  TCut met100("pfmet > 100.0");
  TCut met120("pfmet > 120.0");
  TCut met140("pfmet > 140.0");
  TCut met160("pfmet > 160.0");
  TCut met180("pfmet > 180.0");
  TCut met200("pfmet > 200.0");

  // sigcuts.push_back(TCut(presel+met60));      signames.push_back("E_{T}^{miss} > 60 GeV");      labels.push_back("met60");  cuts.push_back(60);
  // sigcuts.push_back(TCut(presel+met80));      signames.push_back("E_{T}^{miss} > 80 GeV");      labels.push_back("met80");  cuts.push_back(80);
  sigcuts.push_back(TCut(presel+met100));     signames.push_back("E_{T}^{miss} > 100 GeV");     labels.push_back("met100"); cuts.push_back(100);
  // sigcuts.push_back(TCut(presel+met120));     signames.push_back("E_{T}^{miss} > 120 GeV");     labels.push_back("met120"); cuts.push_back(120);
  // sigcuts.push_back(TCut(presel+met140));     signames.push_back("E_{T}^{miss} > 140 GeV");     labels.push_back("met140"); cuts.push_back(140);
  // sigcuts.push_back(TCut(presel+met160));     signames.push_back("E_{T}^{miss} > 160 GeV");     labels.push_back("met160"); cuts.push_back(160);
  // sigcuts.push_back(TCut(presel+met180));     signames.push_back("E_{T}^{miss} > 180 GeV");     labels.push_back("met180"); cuts.push_back(180);
  // sigcuts.push_back(TCut(presel+met200));     signames.push_back("E_{T}^{miss} > 200 GeV");     labels.push_back("met200"); cuts.push_back(200);

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

    int   nx   =  26;
    float xmin =  95;
    float xmax = 355;

    int   ny   =  36;
    float ymin =  -5;
    float ymax = 355;
    
    heff[i]      = new TH2F(Form("heff_%i",i)        , Form("heff_%i",i)       ,  nx , xmin , xmax , ny , ymin , ymax );
    heffup[i]    = new TH2F(Form("heffup_%i",i)      , Form("heffup_%i",i)     ,  nx , xmin , xmax , ny , ymin , ymax );
    heffdn[i]    = new TH2F(Form("heffdn_%i",i)      , Form("heffdn_%i",i)     ,  nx , xmin , xmax , ny , ymin , ymax );
    hxsec[i]     = new TH2F(Form("hxsec_%i",i)       , Form("hxsec_%i",i)      ,  nx , xmin , xmax , ny , ymin , ymax );
    hxsec_exp[i] = new TH2F(Form("hxsec_exp_%i",i)   , Form("hxsec_exp_%i",i)  ,  nx , xmin , xmax , ny , ymin , ymax );
    hexcl[i]     = new TH2F(Form("hexcl_%i",i)       , Form("hexcl_%i",i)      ,  nx , xmin , xmax , ny , ymin , ymax );
    hjes[i]      = new TH2F(Form("hjes_%i",i)        , Form("hjes_%i",i)       ,  nx , xmin , xmax , ny , ymin , ymax );

    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i)*weight);
    heff[i]->Scale(1.0/denom);

    // ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut*weight);
    // heffup[i]->Scale(1.0/denom);

    // ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut*weight);
    // heffdn[i]->Scale(1.0/denom);

    for( unsigned int ibin = 1 ; ibin <= nx ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= ny ; jbin++ ){

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin);
	float ml = heff[i]->GetYaxis()->GetBinCenter(jbin);

	float eff    = heff[i]->GetBinContent(ibin,jbin);
	// float effup  = heffup[i]->GetBinContent(ibin,jbin);
	// float effdn  = heffdn[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	// float dup    = effup/eff-1;
	// float ddn    = 1-effdn/eff;
	// float djes   = 0.5 * (dup+ddn);
	// hjes[i]->SetBinContent(ibin,jbin,djes);

	// float toterr = sqrt( 0.04*0.04 + 0.05*0.05 + 0.05*0.05 + djes*djes );

	//float this_ul = getObservedLimit( cuts.at(i) , toterr );
	float this_ul = getObservedUpperLimit( cuts.at(i) );
	float xsecul  = this_ul / ( lumi * eff );

	//float this_ul_exp = getExpectedLimit( cuts.at(i) , toterr );
	float this_ul_exp = getExpectedUpperLimit( cuts.at(i) );
	float xsecul_exp  = this_ul_exp / ( lumi * eff );

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
    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1000,1000);
    can[i]->Divide(2,2);

    //-------------------------------
    // efficiency
    //-------------------------------

    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->Scale(1000);
    heff[i]->GetYaxis()->SetTitle(ytitle);
    heff[i]->GetXaxis()->SetTitle(xtitle);
    heff[i]->GetZaxis()->SetTitle("efficiency (10^{-3})");
    heff[i]->Draw("colz");

    
    t->DrawLatex(0.2 ,0.83,process);
    t->DrawLatex(0.2 ,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,label);

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
    hxsec[i]->SetMaximum(100);

    // TGraph* gr_excl      = getRefXsecGraph(hxsec[i], "T5zz", 1.0);
    // TGraph* gr_excl_down = getRefXsecGraph(hxsec[i], "T5zz", 1./3.);
    // TGraph* gr_excl_up   = getRefXsecGraph(hxsec[i], "T5zz", 3.);

    // gr_excl->SetLineWidth(2);
    // gr_excl_up->SetLineWidth(2);
    // gr_excl_down->SetLineWidth(2);
    // gr_excl_up->SetLineStyle(2);
    // gr_excl_down->SetLineStyle(3);
    // gr_excl->Draw("same");
    // gr_excl_up->Draw("same");
    // gr_excl_down->Draw("same");

    // TLegend *leg = new TLegend(0.2,0.53,0.53,0.67);
    // leg->AddEntry(gr_excl,"#sigma^{prod} = #sigma^{NLO-QCD}","l");
    // leg->AddEntry(gr_excl_up,"#sigma^{prod} = 3 #times #sigma^{NLO-QCD}","l");
    // leg->AddEntry(gr_excl_down,"#sigma^{prod} = 1/3 #times #sigma^{NLO-QCD}","l");
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->Draw();

    // t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    // t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    // t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    // t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    t->DrawLatex(0.2 ,0.83,process);
    t->DrawLatex(0.2 ,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,label);

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

    // t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    // t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    // t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    // t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    t->DrawLatex(0.2 ,0.83,process);
    t->DrawLatex(0.2 ,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,label);

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

    t->DrawLatex(0.2 ,0.83,process);
    t->DrawLatex(0.2 ,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,label);

    // t->DrawLatex(0.2,0.83,"pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow 2j+#chi_{2}^{0}, #chi_{2}^{0} #rightarrow Z #chi_{1}^{0}");
    // t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    // t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    // t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.7 fb^{-1}");

    if( print ){
      can[i]->Print(Form("../../plots/%s.pdf",labels.at(i).c_str()));
      //can[i]->Print(Form("../plots/%s.eps",labels.at(i).c_str()));
      //gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.eps  ../plots/%s.pdf",labels.at(i).c_str(),labels.at(i).c_str()));
    }

    int bin = heff[i]->FindBin(250,50);
    cout << "efficiency (250,50)  " << heff[i]->GetBinContent(bin) << endl;
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


float getObservedUpperLimit( int metcut ){
  float ul = 9999.;
  if(metcut == 60)  ul = 161.1;
  if(metcut == 80)  ul = 26.6;
  if(metcut == 100) ul = 12.6;
  if(metcut == 120) ul = 9.1;
  if(metcut == 140) ul = 7.2;
  if(metcut == 160) ul = 5.8;
  if(metcut == 180) ul = 4.4;
  if(metcut == 200) ul = 3.6;
  return ul;
}


float getExpectedUpperLimit( int metcut ){
  float ul = 9999.;
  if(metcut == 60)  ul = 160.5;
  if(metcut == 80)  ul = 26.4;
  if(metcut == 100) ul = 12.9;
  if(metcut == 120) ul = 9.4;
  if(metcut == 140) ul = 7.2;
  if(metcut == 160) ul = 5.9;
  if(metcut == 180) ul = 4.8;
  if(metcut == 200) ul = 3.9;
  return ul;
}



