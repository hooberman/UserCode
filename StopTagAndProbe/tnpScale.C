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
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;

int iplot = 0;

void printline(TH2F* h2)
{
  for (unsigned int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {

    Float_t min = h2->GetXaxis()->GetBinLowEdge(x);
    Float_t max = min + h2->GetXaxis()->GetBinWidth(x);
    printf("  %4.1f - %4.1f  & ", min, max);

    for (unsigned int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
      {
	Float_t eff = h2->GetBinContent(x, y);
	Float_t err = h2->GetBinError(x, y);
	if (y == h2->GetYaxis()->GetNbins())
	  printf("\t%4.4f $\\pm$ %4.4f \\\\", eff, err);
	else
	  printf("\t%4.4f $\\pm$ %4.4f & ", eff, err);
      }
    printf("\n");
  }
}

void printHisto( TChain *data , TChain *mc , TCut num , TCut denom , char* var , int nbins , float xmin , float xmax ){

  TH1F* hpass_data = new TH1F(Form("hpass_data_%i",iplot),Form("hpass_data_%i",iplot),nbins,xmin,xmax);
  TH1F* hall_data  = new TH1F(Form("hall_data_%i" ,iplot),Form("hall_data_%i" ,iplot),nbins,xmin,xmax);
  TH1F* hpass_mc   = new TH1F(Form("hpass_mc_%i"  ,iplot),Form("hpass_mc_%i"  ,iplot),nbins,xmin,xmax);
  TH1F* hall_mc    = new TH1F(Form("hall_mc_%i"   ,iplot),Form("hall_mc_%i"   ,iplot),nbins,xmin,xmax);

  data->Draw(Form("min(%s,%f)>>hpass_data_%i"  , var,xmax-0.0001,iplot),denom+num);
  data->Draw(Form("min(%s,%f)>>hall_data_%i"   , var,xmax-0.0001,iplot),denom);
  mc->Draw  (Form("min(%s,%f)>>hpass_mc_%i"    , var,xmax-0.0001,iplot),denom+num);
  mc->Draw  (Form("min(%s,%f)>>hall_mc_%i"     , var,xmax-0.0001,iplot),denom);

  TGraphAsymmErrors *grdata = new TGraphAsymmErrors();
  grdata->BayesDivide(hpass_data,hall_data);


  TCanvas *can = new TCanvas(Form("can_%i",iplot),Form("can_%i",iplot),600,600);
  can->cd();
  grdata->Draw("AP");

  iplot ++;
}


void tnpScale() {

  //----------------------------------------
  // Files
  //----------------------------------------

  char* version = "V00-00-00";

  TChain *chmc = new TChain("leptons");
  chmc->Add(Form("smurf/%s/dymm_test.root",version));

  TChain *chdata = new TChain("leptons");
  //chdata->Add(Form("smurf/%s/data_May10.root"      , version));
  //chdata->Add(Form("smurf/%s/data_PRv4.root"       , version));
  chdata->Add(Form("smurf/%s/data_Aug05.root"      , version));
  //chdata->Add(Form("smurf/%s/data_PRv6.root"       , version));
  //chdata->Add(Form("smurf/%s/data_2011B-V33.root"  , version));
  //chdata->Add(Form("smurf/%s/data_2011B-V34.root"  , version));
  //  chdata->Add("smurf/data_test.root");

  //----------------------------------------
  // bins 
  //----------------------------------------

  float ptbin[] = {10., 15., 20., 30., 40., 50., 7000.};
  float etabin[] = {0, 0.8, 1.479, 2.0, 2.5};
  int nptbin=6;
  int netabin=4;

  //
  // histogram
  //
  //deno
  TH2F *hmcid_deno 	= new TH2F("hmcid_deno", "hmcid_deno", nptbin, ptbin, netabin, etabin);
  TH2F *hmciso_deno 	= new TH2F("hmciso_deno", "hmciso_deno", nptbin, ptbin, netabin, etabin);
  TH2F *hdataid_deno 	= new TH2F("hdataid_deno", "hdataid_deno", nptbin, ptbin, netabin, etabin);
  TH2F *hdataiso_deno	= new TH2F("hdataiso_deno", "hdataiso_deno", nptbin, ptbin, netabin, etabin);
  hmcid_deno->Sumw2();
  hmciso_deno->Sumw2();
  hdataid_deno->Sumw2();
  hdataiso_deno->Sumw2();
  //num
  TH2F *hmcid_num 	= new TH2F("hmcid_num", "hmcid_num", nptbin, ptbin, netabin, etabin);
  TH2F *hmciso_num 	= new TH2F("hmciso_num", "hmciso_num", nptbin, ptbin, netabin, etabin);
  TH2F *hdataid_num 	= new TH2F("hdataid_num", "hdataid_num", nptbin, ptbin, netabin, etabin);
  TH2F *hdataiso_num 	= new TH2F("hdataiso_num", "hdataiso_num", nptbin, ptbin, netabin, etabin);
  hmcid_num->Sumw2();
  hmciso_num->Sumw2();
  hdataid_num->Sumw2();
  hdataiso_num->Sumw2();
  // eff
  TH2F *hmcid 	= new TH2F("hmcid", "hmcid", nptbin, ptbin, netabin, etabin);
  TH2F *hmciso 	= new TH2F("hmciso", "hmciso", nptbin, ptbin, netabin, etabin);
  TH2F *hdataid 	= new TH2F("hdataid", "hdataid", nptbin, ptbin, netabin, etabin);
  TH2F *hdataiso 	= new TH2F("hdataiso", "hdataiso", nptbin, ptbin, netabin, etabin);
  hmcid->Sumw2();
  hmciso->Sumw2();
  hdataid->Sumw2();
  hdataiso->Sumw2();
  // SF
  TH2F *hsfid 	= new TH2F("hsfid", "hsfid", nptbin, ptbin, netabin, etabin);
  TH2F *hsfiso 	= new TH2F("hsfiso", "hsfiso", nptbin, ptbin, netabin, etabin);
  hsfid->Sumw2();
  hsfiso->Sumw2();

  //------------------------------------------
  // cuts
  //------------------------------------------

  //TString tnpcut 	= "abs(tagAndProbeMass-91)<15&&(eventSelection&1)==1&&HLT_Ele27_WP80_tag>0&&qProbe*qTag<0"; 
  //TString tnpcut 	= "abs(tagAndProbeMass-91)<15&&(eventSelection&2)==2&&qProbe*qTag<0"; 

  TString tnpcut 	= "abs(tagAndProbeMass-91)<15&&(eventSelection&2)==2&&HLT_IsoMu30_eta2p1_tag>0&&qProbe*qTag<0"; 
  TString mufo 	        = "(leptonSelection&32768)==32768";    // mu fo
  TString muid 	        = "(leptonSelection&65536)==65536";    // mu id 
  TString muiso 	= "(leptonSelection&131072)==131072";  // mu iso 
  TString elfo        	= "(leptonSelection&4)==4";            // ele fo 
  TString elid  	= "(leptonSelection&8)==8";            // ele id 
  TString eliso 	= "(leptonSelection&16)==16";          // ele iso

  //printHisto( chdata , chmc , muiso , TCut

  cout << "Total MC yields 	: " << chmc->GetEntries(tnpcut) << endl;
  cout << "Total DATA yields 	: " << chdata->GetEntries(tnpcut) << endl;

  //
  // Fill histograms
  //
  // chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_deno", 	tnpcut+"&&"+eliso,				"goff");
  // chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_num", 		tnpcut+"&&"+eliso+"&&"+elid,	"goff");
  // chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_deno", 	tnpcut+"&&"+elid,				"goff");
  // chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_num", 	tnpcut+"&&"+elid+"&&"+eliso,	"goff");
  // chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_deno", 	tnpcut+"&&"+eliso,				"goff");
  // chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_num", 		tnpcut+"&&"+eliso+"&&"+elid,	"goff");
  // chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_deno", 	tnpcut+"&&"+elid,				"goff");
  // chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_num", 	tnpcut+"&&"+elid+"&&"+eliso,	"goff");

  chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_deno", 	tnpcut+"&&"+muiso,	      	"goff");
  chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_num", 		tnpcut+"&&"+muiso+"&&"+muid,	"goff");
  chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_deno", 	tnpcut+"&&"+muid,	      	"goff");
  chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_num", 	tnpcut+"&&"+muid+"&&"+muiso,	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_deno", 	tnpcut+"&&"+muiso,		"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_num", 	tnpcut+"&&"+muiso+"&&"+muid,	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_deno", 	tnpcut+"&&"+muid,	       	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_num", 	tnpcut+"&&"+muid+"&&"+muiso,	"goff");

  // get efficiencies 
  hmcid->Divide(hmcid_num,hmcid_deno,1,1,"B");
  hmciso->Divide(hmciso_num,hmciso_deno,1,1,"B");
  hdataid->Divide(hdataid_num,hdataid_deno,1,1,"B");
  hdataiso->Divide(hdataiso_num,hdataiso_deno,1,1,"B");
	
  // get scale factors
  hsfid->Divide(hdataid, hmcid, 1, 1);
  hsfiso->Divide(hdataiso, hmciso, 1, 1);

  // Draw histograms	
  //hmcid->Draw("text");
	
  // print table
  cout << " ------ MC ID ----- " << endl;
  printline(hmcid);
  cout << " ------ MC ISO ----- " << endl;
  printline(hmciso);
  cout << " ------ DATA ID ----- " << endl;
  printline(hdataid);
  cout << " ------ DATA ISO ----- " << endl;
  printline(hdataiso);
  cout << " ------ Scale Factor ID ----- " << endl;
  printline(hsfid);
  cout << " ------ Scale Factor ISO ----- " << endl;
  printline(hsfiso);
	
}
