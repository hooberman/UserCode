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

//-------------------------------------------
// THING TO ADD
//-------------------------------------------
/*

2) b-tagging uncertanties

*/


float getObservedLimit( float seff , string SR );
float getExpectedLimit( float seff , string SR );

using namespace std;

void SMS(char* sample , int x = 1, bool print = false){

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------
  
  const float denom    = 50000;
  const float lumi     = 4980;
  const char* filename = Form(Form("/tas/vimartin/SingleLepton2011/output/V00-04-08/%s_smallTree.root",sample));
  const float btagerr  = 0.04;

  cout << "Using file        " << filename << endl;
  cout << "Using denominator " << denom    << " events" << endl;
  cout << "Using lumi        " << lumi     << " pb-1" << endl;

  char* label = (char*)"";
  if     ( TString(sample).Contains("T2tt") ) label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow t+#chi_{1}^{0}";
  else if( TString(sample).Contains("T2bw") ){
    if( x == 25 ){
      cout << "Doing x=0.25 slice" << endl;
      label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#chi_{1}^{#pm}, x=0.25";
    }
    else if( x == 50 ){
      cout << "Doing x=0.50 slice" << endl;
      label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#chi_{1}^{#pm}, x=0.50";
    }
    else if( x == 75 ){
      cout << "Doing x=0.75 slice" << endl;
      label = (char*)"pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow b+#chi_{1}^{#pm}, x=0.75";
    }
    else{
      cout << "ERROR! unrecognized x value " << x << endl;
      exit(0);
    }
  }
  else{
    cout << "ERROR! unrecognized sample " << sample << endl;
    exit(0);
  }

  //if( TString(sample).Contains("T2tt") ) label = "";

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("t");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("stop_xsec.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("h_stop_xsec");

  //--------------------------------------------------
  // preselection
  //--------------------------------------------------

  TCut rho("rhovor>0 && rhovor<40");
  TCut goodlep("ngoodlep > 0 && leptype==1 && lep1->Pt()>30 && abs(lep1->Eta())<2.1");
  TCut njets4("npfjets30 >= 4");
  TCut btag1("nbtagsssvcorr>=1");
  TCut isotrk("pfcandpt10 > 9998. || pfcandiso10 > 0.1");
  TCut SRA("t1metphicorr > 100 && t1metphicorrmt > 150");
  TCut SRB("t1metphicorr > 150 && t1metphicorrmt > 120");

  TCut weight("ndavtxweight * mutrigweight");
  //TCut weight("1");

  TCut presel  = rho + goodlep + njets4 + btag1 + isotrk;

  if( x == 25 ) presel += TCut("x==0.25");
  if( x == 50 ) presel += TCut("x==0.50");
  if( x == 75 ) presel += TCut("x==0.75");

  cout << "Using pre-selection   " << presel.GetTitle() << endl;
  cout << "Using weight          " << weight.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<int>     cuts;

  sigcuts.push_back(TCut(presel+SRA));  signames.push_back("SRA");  labels.push_back("SRA");  cuts.push_back(1);
  //sigcuts.push_back(TCut(presel+SRB));  signames.push_back("SRB");  labels.push_back("SRB");  cuts.push_back(1);

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
    jesup.ReplaceAll("npfjets30"       , "njetsUp");
    jesup.ReplaceAll("t1metphicorr "   , "t1metphicorrup ");
    jesup.ReplaceAll("t1metphicorrmt " , "t1metphicorrmtup ");

    TString jesdn(sigcuts.at(i));
    jesdn.ReplaceAll("npfjets30"       , "njetsDown");
    jesdn.ReplaceAll("t1metphicorr "   , "t1metphicorrdn ");
    jesdn.ReplaceAll("t1metphicorrmt " , "t1metphicorrmtdn ");

    TCut jesupcut(jesup);
    TCut jesdncut(jesdn);

    cout << endl << endl;
    cout << "Signal region : " << labels.at(i)  << endl;
    cout << "Selection     : " << sigcuts.at(i) << endl;
    cout << "Selection up  : " << jesupcut      << endl;
    cout << "Selection dn  : " << jesdncut      << endl;

    int   nbinsx  = 17;
    float xmin    = 187.5;
    float xmax    = 612.5;
    int   nbinsy  = 13;
    float ymin    = -12.5;
    float ymax    = 312.5;

    if( TString(sample).Contains("T2bw") ){
      nbinsx =   9;
      xmin   = 175;
      xmax   = 625;
      nbinsy =   7;
      ymin   = -25;
      ymax   = 352;
    }

    heff[i]      = new TH2F(Form("heff_%i",i)        , Form("heff_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffup[i]    = new TH2F(Form("heffup_%i",i)      , Form("heffup_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffdn[i]    = new TH2F(Form("heffdn_%i",i)      , Form("heffdn_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec[i]     = new TH2F(Form("hxsec_%i",i)       , Form("hxsec_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_exp[i] = new TH2F(Form("hxsec_exp_%i",i)   , Form("hxsec_exp_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl[i]     = new TH2F(Form("hexcl_%i",i)       , Form("hexcl_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hjes[i]      = new TH2F(Form("hjes_%i",i)        , Form("hjes_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    ch->Draw(Form("ml:mg>>heff_%i",i),sigcuts.at(i)*weight);
    heff[i]->Scale(1.0/denom);

    ch->Draw(Form("ml:mg>>heffup_%i",i),jesupcut*weight);
    heffup[i]->Scale(1.0/denom);

    ch->Draw(Form("ml:mg>>heffdn_%i",i),jesdncut*weight);
    heffdn[i]->Scale(1.0/denom);

    for( unsigned int ibin = 1 ; ibin <= 25 ; ibin++ ){
      for( unsigned int jbin = 1 ; jbin <= 25 ; jbin++ ){

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin);
	//float ml = heff[i]->GetYaxis()->GetBinCenter(jbin);

	float eff    = heff[i]->GetBinContent(ibin,jbin);
	float effup  = heffup[i]->GetBinContent(ibin,jbin);
	float effdn  = heffdn[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	float dup    = effup/eff-1;
	float ddn    = 1-effdn/eff;
	float djes   = 0.5 * (dup+ddn);

	//djes = 0.10;

	hjes[i]->SetBinContent(ibin,jbin,djes);

	// lumi, lepton selection, trigger, b-tagging, JES

	float toterr  = sqrt( 0.022*0.022 + 0.02*0.02 + 0.03*0.03 + btagerr*btagerr + djes*djes );

	float this_ul = getObservedLimit( toterr , labels.at(i) );
	float xsecul  = this_ul / ( lumi * eff );

	float this_ul_exp = getExpectedLimit( toterr , labels.at(i) );
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
    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,1200);
    can[i]->Divide(2,2);

    //-------------------------------
    // efficiency
    //-------------------------------

    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    heff[i]->Scale(100);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    heff[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    heff[i]->GetZaxis()->SetTitle("efficiency (%)");
    heff[i]->Draw("colz");
    heff[i]->Draw("sametext");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");

    //-------------------------------
    // cross section
    //-------------------------------
  
    can[i]->cd(2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    gPad->SetLogz();
  
    hxsec_exp[i]->GetXaxis()->SetLabelSize(0.035);
    hxsec_exp[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hxsec_exp[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hxsec_exp[i]->GetZaxis()->SetTitle("#sigma upper limit");
    hxsec_exp[i]->Draw("colz");
    gStyle->SetPaintTextFormat(".2f");
    hxsec_exp[i]->Draw("sametext");
    hxsec_exp[i]->SetMinimum(0.01);
    hxsec_exp[i]->SetMaximum(100);

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

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");

    //-------------------------------
    // excluded points
    //-------------------------------

    can[i]->cd(3);    
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hexcl[i]->GetZaxis()->SetTitle("excluded points");
    hexcl[i]->Draw("colz");
    //gr_excl->Draw("same");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");

    //-------------------------------
    // JES uncertainty
    //-------------------------------

    can[i]->cd(4);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hjes[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hjes[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hjes[i]->GetZaxis()->SetTitle("JES uncertainty");
    hjes[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");

    if( print ){
      can[i]->Print(Form("../../plots/%s_%s.pdf",sample,labels.at(i).c_str()));
      //can[i]->Print(Form("../plots/%s.eps",labels.at(i).c_str()));
      //gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.eps  ../plots/%s.pdf",labels.at(i).c_str(),labels.at(i).c_str()));
    }

    int bin = heff[i]->FindBin(300,50);

    float toterr = sqrt(pow(hjes[i]->GetBinContent(bin),2)+0.022*0.022 + 0.02*0.02 + 0.03*0.03 + btagerr*btagerr);
    cout << "efficiency (300,50)  " << heff[i]->GetBinContent(bin) << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin) << endl;
    cout << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    cout << "JES                  " << hjes[i]->GetBinContent(bin) << endl;
    cout << "tot err              " << toterr << endl;
    cout << "obs limit            " << getObservedLimit(toterr,labels.at(i)) << endl;
    cout << "exp limit            " << getExpectedLimit(toterr,labels.at(i)) << endl;
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


//-----------------------------------------
// the following UL's correspdond to:
// SRA: observed 116, predicted 112.7 8.3% uncertainty
// SRB: observed  79, predicted 58.75 9.1% uncertainty
// predicted bkg  = 116, 10% uncertainty
//-----------------------------------------

float getObservedLimit( float seff , string SR ){
  float ul = 9999.;

  if( TString(SR).Contains("SRA") ){
    if(seff >= 0.00 && seff < 0.05) ul = 30.3;
    if(seff >= 0.05 && seff < 0.10) ul = 30.7;
    if(seff >= 0.10 && seff < 0.15) ul = 31.6;
    if(seff >= 0.15 && seff < 0.20) ul = 32.3;
    if(seff >= 0.20 && seff < 0.25) ul = 36.0;
  }

  else if( TString(SR).Contains("SRB") ){
    if(seff >= 0.00 && seff < 0.05) ul = 38.3;
    if(seff >= 0.05 && seff < 0.10) ul = 39.3;
    if(seff >= 0.10 && seff < 0.15) ul = 40.4;
    if(seff >= 0.15 && seff < 0.20) ul = 42.1;
    if(seff >= 0.20 && seff < 0.25) ul = 47.7;
  }

  else{
    cout << "ERROR! unrecognized SR name " << SR << endl;
    exit(0);
  }

  if( ul > 1000 ){
    cout << "WARNING: signal efficiency uncertainty " << seff << ", returning ul 9999" << endl;
  }

  return ul;
}


float getExpectedLimit( float seff , string SR ){
  float ul = 9999.;

  if( TString(SR).Contains("SRA") ){
    if(seff >= 0.00 && seff < 0.05) ul = 27.9;
    if(seff >= 0.05 && seff < 0.10) ul = 28.2;
    if(seff >= 0.10 && seff < 0.15) ul = 28.8;
    if(seff >= 0.15 && seff < 0.20) ul = 29.5;
    if(seff >= 0.20 && seff < 0.25) ul = 32.4;
  }

  else if( TString(SR).Contains("SRB") ){
    if(seff >= 0.00 && seff < 0.05) ul = 20.1;
    if(seff >= 0.05 && seff < 0.10) ul = 22.2;
    if(seff >= 0.10 && seff < 0.15) ul = 20.8;
    if(seff >= 0.15 && seff < 0.20) ul = 21.4;
    if(seff >= 0.20 && seff < 0.25) ul = 24.0;
  }

  else{
    cout << "ERROR! unrecognized SR name " << SR << endl;
    exit(0);
  }

  if( ul > 1000 ){
    cout << "WARNING: signal efficiency uncertainty " << seff << ", returning ul 9999" << endl;
  }

  return ul;
}


// float getExpectedP1UpperLimit( float seff ){
//   float ul = 9999.;
//   if(seff >= 0.00 && seff < 0.05) ul = 41.5;
//   if(seff >= 0.05 && seff < 0.10) ul = 42.7;
//   if(seff >= 0.10 && seff < 0.15) ul = 44.3;
//   if(seff >= 0.15 && seff < 0.20) ul = 46.2;
//   if(seff >= 0.20 && seff < 0.25) ul = 54.3;
//   if(seff >= 0.25 && seff < 0.30) ul = 51.3;
//   return ul;
// }


// float getExpectedM1UpperLimit( float seff ){
//   float ul = 9999.;
//   if(seff >= 0.00 && seff < 0.05) ul = 21.6;
//   if(seff >= 0.05 && seff < 0.10) ul = 21.9;
//   if(seff >= 0.10 && seff < 0.15) ul = 22.4;
//   if(seff >= 0.15 && seff < 0.20) ul = 22.7;
//   if(seff >= 0.20 && seff < 0.25) ul = 24.1;
//   if(seff >= 0.25 && seff < 0.30) ul = 23.6;
//   return ul;
// }
