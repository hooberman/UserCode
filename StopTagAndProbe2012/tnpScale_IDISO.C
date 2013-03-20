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

//------------------------------
// print header for tables
//------------------------------

void printHeader( int leptype , char* type ){

  cout << endl << endl;
  cout << "\\hline" << endl << "\\hline" << endl;
  if( leptype == 0 ){
    cout << type << " & & \\\\" << endl;
    cout <<  "\\pt\\ range [GeV] & $|\\eta|<0.8$ & $0.8<|\\eta|<1.4442$ \\\\" << endl;
  }
  if( leptype == 1 ){
    cout << type << "& & & \\\\" << endl;
    cout <<  "\\pt\\ range [GeV] & $|\\eta|<0.8$ & $0.8<|\\eta|<1.5$ & $1.5<|\\eta|<2.1$ \\\\" << endl;
  }
  cout << "\\hline" << endl;
}

//------------------------------
// print line in tables
//------------------------------

void printline(TH2F* h2)
{
  for (int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {

    Float_t min = h2->GetXaxis()->GetBinLowEdge(x);
    Float_t max = min + h2->GetXaxis()->GetBinWidth(x);
    printf("  %4.0f - %4.0f  & ", min, max);

    for (int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
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

//------------------------------
// make data/MC comparison plot
//------------------------------

void printHisto( TCanvas *can , TChain *data , TChain *mc , TCut num , TCut denom , char* var , int nbins , float xmin , float xmax , char* xtitle , char* ytitle){

  can->cd();

  TPad *plotpad = new TPad("plotpad","plotpad",0.0,0.0,1.0,0.8);
  plotpad->Draw();
  plotpad->cd();

  // float ptbin[] = {30., 40., 60. , 80. , 100. , 120. , 150. , 200 , 300 };
  // int   nptbin  = 8;

  float ptbin[] = { 20., 30. , 40. , 50. , 60. , 80.0 , 100.0 , 150.0 , 200.0 , 300.0 , 350.0 };
  int   nptbin  = 10;

  // float ptbin[] = { 20., 30. , 40. , 50. , 60. , 80.0 , 100.0 , 150.0 , 200.0 , 300.0 , 500.0 , 1000.0 , 1100.0 };
  // int   nptbin  = 12;

  xmax=ptbin[nptbin];

  cout << "Using xmax " << xmax << endl;

  //TH1F* hpass   = new TH1F(Form("hpass_%i",iplot),Form("hpass_%i",iplot),nptbin,ptbin);
  //TH1F* hall    = new TH1F(Form("hall_%i" ,iplot),Form("hall_%i" ,iplot),nptbin,ptbin);

  // TH1F* hpass_data = new TH1F(Form("hpass_data_%i",iplot),Form("hpass_data_%i",iplot),nbins,xmin,xmax);
  // TH1F* hall_data  = new TH1F(Form("hall_data_%i" ,iplot),Form("hall_data_%i" ,iplot),nbins,xmin,xmax);
  // TH1F* hpass_mc   = new TH1F(Form("hpass_mc_%i"  ,iplot),Form("hpass_mc_%i"  ,iplot),nbins,xmin,xmax);
  // TH1F* hall_mc    = new TH1F(Form("hall_mc_%i"   ,iplot),Form("hall_mc_%i"   ,iplot),nbins,xmin,xmax);

  TH1F* hpass_data   = new TH1F(Form("hpass_data_%i"  ,iplot),Form("hpass_data_%i"  ,iplot),nptbin,ptbin);
  TH1F* hall_data    = new TH1F(Form("hall_data_%i"   ,iplot),Form("hall_data_%i"   ,iplot),nptbin,ptbin);
  TH1F* hratio_data  = new TH1F(Form("hratio_data_%i" ,iplot),Form("hratio_data_%i" ,iplot),nptbin,ptbin);
  TH1F* hpass_mc     = new TH1F(Form("hpass_mc_%i"    ,iplot),Form("hpass_mc_%i"    ,iplot),nptbin,ptbin);
  TH1F* hall_mc      = new TH1F(Form("hall_mc_%i"     ,iplot),Form("hall_mc_%i"     ,iplot),nptbin,ptbin);
  TH1F* hratio_mc    = new TH1F(Form("hratio_mc_%i"   ,iplot),Form("hratio_mc_%i"   ,iplot),nptbin,ptbin);
  TH1F* hsf          = new TH1F(Form("hsf_%i"         ,iplot),Form("hsf_%i"         ,iplot),nptbin,ptbin);

  hpass_data->Sumw2();
  hall_data->Sumw2();
  hratio_data->Sumw2();
  hpass_mc->Sumw2();
  hall_mc->Sumw2();
  hratio_mc->Sumw2();
  hsf->Sumw2();

  //TCanvas *can = new TCanvas(Form("can_%i",iplot),Form("can_%i",iplot),600,600);
  //can->cd();
  
  data->Draw(Form("min(%s,%f)>>hpass_data_%i"  , var,xmax-0.0001,iplot),denom+num);
  data->Draw(Form("min(%s,%f)>>hall_data_%i"   , var,xmax-0.0001,iplot),denom);
  mc->Draw  (Form("min(%s,%f)>>hpass_mc_%i"    , var,xmax-0.0001,iplot),denom+num);
  mc->Draw  (Form("min(%s,%f)>>hall_mc_%i"     , var,xmax-0.0001,iplot),denom);

  TGraphAsymmErrors *grdata = new TGraphAsymmErrors();
  grdata->BayesDivide(hpass_data,hall_data);

  TGraphAsymmErrors *grmc = new TGraphAsymmErrors();
  grmc->BayesDivide(hpass_mc,hall_mc);

  cout << "data all  " << hall_data->GetBinContent(2) << endl;
  cout << "data pass " << hpass_data->GetBinContent(2) << endl;
  cout << "data eff  " << hpass_data->GetBinContent(2) / hall_data->GetBinContent(2) << endl;

  Double_t x;
  Double_t y;
  grdata->GetPoint(1,x,y);
  cout << "data eff2 " << y << endl;

  cout << "MC all  " << hall_mc->GetBinContent(2) << endl;
  cout << "MC pass " << hpass_mc->GetBinContent(2) << endl;
  cout << "MC eff  " << hpass_mc->GetBinContent(2) / hall_mc->GetBinContent(2) << endl;

  Double_t xmc;
  Double_t ymc;
  grmc->GetPoint(1,xmc,ymc);
  cout << "MC eff2 " << ymc << endl;

  gPad->SetGridx();
  gPad->SetGridy();

  grdata->SetMarkerColor(2);
  grdata->SetLineColor(2);
  grmc->SetMarkerColor(4);
  grmc->SetLineColor(4);
  grmc->SetMarkerStyle(25);

  grdata->GetXaxis()->SetRangeUser(20,350);
  if( TString(ytitle).Contains("iso") ) grdata->GetYaxis()->SetRangeUser(0.8,1.0);
  if( TString(ytitle).Contains("ID")  ) grdata->GetYaxis()->SetRangeUser(0.5,1.0);
  grdata->GetXaxis()->SetTitle(xtitle);
  grdata->GetYaxis()->SetTitle(ytitle);
  grdata->Draw("AP");
  grmc->Draw("sameP");

  TLegend *leg = new TLegend(0.6,0.15,0.75,0.3);
  leg->AddEntry(grdata ,"data","lp");
  leg->AddEntry(grmc   ,"mc"  ,"lp");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();

  if( TString(denom.GetTitle()).Contains("njets>=0") ) t->DrawLatex(0.2,0.2,"n_{jets} #geq 0");
  if( TString(denom.GetTitle()).Contains("njets>=1") ) t->DrawLatex(0.2,0.2,"n_{jets} #geq 1");
  if( TString(denom.GetTitle()).Contains("njets>=2") ) t->DrawLatex(0.2,0.2,"n_{jets} #geq 2");
  if( TString(denom.GetTitle()).Contains("njets>=3") ) t->DrawLatex(0.2,0.2,"n_{jets} #geq 3");
  if( TString(denom.GetTitle()).Contains("njets>=4") ) t->DrawLatex(0.2,0.2,"n_{jets} #geq 4");

  can->cd();

  TPad *respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
  respad->Draw();
  respad->cd();
  respad->SetGridy();

  // TGraphAsymmErrors* gr_ratio = (TGraphAsymmErrors*) grdata->Clone("gr_ratio");
  // gr_ratio->Divide(grmc);
  // gr_ratio->Draw();

  hratio_data->Divide(hpass_data,hall_data,1,1,"B");
  hratio_mc  ->Divide(hpass_mc  ,hall_mc,1,1,"B");
  hsf        ->Divide(hratio_data,hratio_mc,1,1);

  if( TString(ytitle).Contains("iso") )   hsf->GetYaxis()->SetRangeUser(0.9,1.1);
  if( TString(ytitle).Contains("ID") )    hsf->GetYaxis()->SetRangeUser(0.6,1.1);

  hsf->GetXaxis()->SetRangeUser(20,350);
  hsf->GetYaxis()->SetNdivisions(3);
  hsf->GetYaxis()->SetLabelSize(0.2);
  hsf->GetXaxis()->SetLabelSize(0.0);
    
  hsf->Draw("E1");

  iplot ++;
}

//------------------------------
// main function
//------------------------------

void tnpScale_IDISO( int leptype = 1, bool printplot = false ) {

  cout << endl;
  cout << "-------------------" << endl;
  if     ( leptype == 0 ) cout << "Doing electrons" << endl;
  else if( leptype == 1 ) cout << "Doing muons"      << endl;
  else{
    cout << "ERROR! unrecognized leptype " << leptype << endl;
    exit(0);
  }
  cout << "-------------------" << endl;

  //----------------------------------------
  // Files
  //----------------------------------------

  TChain *chmc   = new TChain("leptons");
  TChain *chdata = new TChain("leptons");

  char* version = (char*) "V00-00-06";
  char* suffix = "";
  //char* suffix = "_2jets";
  //char* suffix = "_probept100";


  chmc->  Add(Form("smurf/ZJetsFull_%s/merged%s.root",version,suffix));

  if( leptype == 1 ){
    chdata->Add(Form("smurf/SingleMu2012AFull_%s/merged_json%s.root",version,suffix));
    chdata->Add(Form("smurf/SingleMu2012BFull_%s/merged_json%s.root",version,suffix));
    chdata->Add(Form("smurf/SingleMu2012CFull_%s/merged_json%s.root",version,suffix));
    chdata->Add(Form("smurf/SingleMu2012DFull_%s/merged_json%s.root",version,suffix));
  }
  else{
    chdata->Add(Form("smurf/SingleEl2012AFull_%s/merged_json%s.root",version,suffix));
    chdata->Add(Form("smurf/SingleEl2012BFull_%s/merged_json%s.root",version,suffix));
    chdata->Add(Form("smurf/SingleEl2012CFull_%s/merged_json%s.root",version,suffix));
    chdata->Add(Form("smurf/SingleEl2012DFull_%s/merged_json%s.root",version,suffix));
  }

  //----------------------------------------
  // bins 
  //----------------------------------------

  // float ptbin[]  = {10., 15., 20., 30., 40., 50., 7000.};
  // float ptbin[]  = { 30. , 40. , 50. , 60. , 80.0 , 100.0 , 120.0 , 150.0 , 7000.};
  // float etabin[] = {0, 0.8, 1.5, 2.1};
  // int nptbin=8;
  // int netabin=3;

  // float ptbin[]  = { 20., 30. , 40. , 50. , 60. , 80.0 , 100.0 , 150.0 , 200.0 , 300.0 , 500.0 , 1000.0 , 10000000.0};
  // int   nptbin   = 12;

  // float etabin[] = {0,2.1};
  // int   netabin  = 1;

  float ptbin[] = { 20., 30. , 40. , 50. , 60. , 80.0 , 100.0 , 150.0 , 200.0 , 300.0, 10000.0};
  int   nptbin  = 10;

  float etabin[4];
  int   netabin = 0;

  if( leptype == 1 ){
    cout << "DOING MUON ETA BINS" << endl;
    netabin=3;
    etabin[0] = 0.0;
    etabin[1] = 0.8;
    etabin[2] = 1.5;
    etabin[3] = 2.1;
  }

  if( leptype == 0 ){
    cout << "DOING ELECTRON ETA BINS" << endl;
    netabin=2;
    etabin[0] = 0.0;
    etabin[1] = 0.8;
    etabin[2] = 1.4442;

    // netabin=1;
    // etabin[0] = 0.0;
    // etabin[1] = 1.4442;
  }


  //deno
  TH2F *hmcid_deno 	= new TH2F("hmcid_deno"   , "hmcid_deno"   , nptbin, ptbin, netabin, etabin);
  TH2F *hmciso_deno 	= new TH2F("hmciso_deno"  , "hmciso_deno"  , nptbin, ptbin, netabin, etabin);
  TH2F *hdataid_deno 	= new TH2F("hdataid_deno" , "hdataid_deno" , nptbin, ptbin, netabin, etabin);
  TH2F *hdataiso_deno	= new TH2F("hdataiso_deno", "hdataiso_deno", nptbin, ptbin, netabin, etabin);
  hmcid_deno->Sumw2();
  hmciso_deno->Sumw2();
  hdataid_deno->Sumw2();
  hdataiso_deno->Sumw2();

  //num
  TH2F *hmcid_num 	= new TH2F("hmcid_num"    , "hmcid_num"    , nptbin, ptbin, netabin, etabin);
  TH2F *hmciso_num 	= new TH2F("hmciso_num"   , "hmciso_num"   , nptbin, ptbin, netabin, etabin);
  TH2F *hdataid_num 	= new TH2F("hdataid_num"  , "hdataid_num"  , nptbin, ptbin, netabin, etabin);
  TH2F *hdataiso_num 	= new TH2F("hdataiso_num" , "hdataiso_num" , nptbin, ptbin, netabin, etabin);
  hmcid_num->Sumw2();
  hmciso_num->Sumw2();
  hdataid_num->Sumw2();
  hdataiso_num->Sumw2();

  // eff
  TH2F *hmcid 	        = new TH2F("hmcid"        , "hmcid"        , nptbin, ptbin, netabin, etabin);
  TH2F *hmciso 	        = new TH2F("hmciso"       , "hmciso"       , nptbin, ptbin, netabin, etabin);
  TH2F *hdataid 	= new TH2F("hdataid"      , "hdataid"      , nptbin, ptbin, netabin, etabin);
  TH2F *hdataiso 	= new TH2F("hdataiso"     , "hdataiso"     , nptbin, ptbin, netabin, etabin);
  hmcid->Sumw2();
  hmciso->Sumw2();
  hdataid->Sumw2();
  hdataiso->Sumw2();

  // SF
  TH2F *hsfid 	= new TH2F("hsfid"  , "hsfid" , nptbin, ptbin, netabin, etabin);
  TH2F *hsfiso 	= new TH2F("hsfiso" , "hsfiso", nptbin, ptbin, netabin, etabin);
  hsfid->Sumw2();
  hsfiso->Sumw2();

  // TCuts
  TCut muid ("(leptonSelection&65536)==65536");     // mu id 
  TCut muiso("(leptonSelection&131072)==131072");   // mu iso 
  TCut elid ("(leptonSelection&8)==8");             // ele id 
  TCut eliso("(leptonSelection&16)==16");           // ele iso

  TCut zmass("abs(tagAndProbeMass-91)<15");
  TCut tightzmass("abs(tagAndProbeMass-91)<5");
  TCut os("qProbe*qTag<0");

  TCut mutnp("(eventSelection&2)==2");
  TCut mutnptrig("HLT_IsoMu24_tag > 0");

  TCut eltnp("(eventSelection&1)==1");
  TCut eltnptrig("HLT_Ele27_WP80_tag > 0");

  TCut tag_eta21("abs(tag->eta())<2.1");
  TCut tag_eta14("abs(tag->eta())<1.4442");
  TCut tag_pt30("tag->pt()>30.0");
  TCut probe_eta21("abs(probe->eta())<2.1");
  TCut probe_eta14("abs(probe->eta())<1.4442");

  TCut mutrk("mutrk==1");

  TCut met30("met<30");
  TCut nbl0("nbl==0");

  TCut njets0("njets>=0");
  TCut njets1("njets>=1");
  TCut njets2("njets>=2");
  TCut njets3("njets>=3");
  TCut njets4("njets>=4");

  TCut mud0("mud0 < 0.02");
  TCut mudz("mudz < 0.5");

  // TCut mud0("mud0 < 0.2");
  // TCut mudz("mudz < 1.0");

  //TCut tnpcut   = "abs(tagAndProbeMass-91)<15 && (eventSelection&2)==2 && HLT_IsoMu30_eta2p1_tag>0 && qProbe*qTag<0 && abs(tag->eta())<2.1 && tag->pt()>30.0"; 

  TCut tnpcut;
  tnpcut += zmass;
  //tnpcut += tightzmass;
  tnpcut += os;
  tnpcut += tag_pt30;

  tnpcut += met30;
  tnpcut += nbl0;

  TCut  lepid;
  TCut  lepiso;
  char* lepchar = "";

  if( leptype == 0 ){
    tnpcut += tag_eta14;
    tnpcut += probe_eta14;
    tnpcut += eltnp;
    tnpcut += eltnptrig;
    lepid   = TCut(elid);
    lepiso  = TCut(eliso);
    lepchar = "el";
  }
  else if( leptype == 1 ){
    tnpcut += tag_eta21;
    tnpcut += probe_eta21;
    //tnpcut += mutrk;
    //tnpcut += mud0;
    //tnpcut += mudz;
    tnpcut += mutnp;
    tnpcut += mutnptrig;
    lepid   = TCut(muid);
    lepiso  = TCut(muiso);
    lepchar = "mu";
  }
  
  //tnpcut += njets2;
  cout << "Selection  : " << tnpcut.GetTitle()          << endl;
  cout << "Ndata      : " << chdata->GetEntries(tnpcut) << endl;
  cout << "NMC        : " << chmc->GetEntries(tnpcut)   << endl;
  cout << "ID cut     : " << lepid.GetTitle()           << endl;
  cout << "iso cut    : " << lepiso.GetTitle()          << endl;

  chmc->  Draw("abs(probe->eta()):probe->pt()>>hmcid_deno", 	tnpcut+lepiso,	      	"goff");
  chmc->  Draw("abs(probe->eta()):probe->pt()>>hmcid_num", 	tnpcut+lepiso+lepid,	"goff");
  chmc->  Draw("abs(probe->eta()):probe->pt()>>hmciso_deno", 	tnpcut+lepid,	      	"goff");
  chmc->  Draw("abs(probe->eta()):probe->pt()>>hmciso_num", 	tnpcut+lepid+lepiso,	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_deno", 	tnpcut+lepiso,		"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_num", 	tnpcut+lepiso+lepid,	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_deno", 	tnpcut+lepid,	       	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_num", 	tnpcut+lepid+lepiso,	"goff");

  // get efficiencies 
  hmcid->Divide(hmcid_num,hmcid_deno,1,1,"B");
  hmciso->Divide(hmciso_num,hmciso_deno,1,1,"B");
  hdataid->Divide(hdataid_num,hdataid_deno,1,1,"B");
  hdataiso->Divide(hdataiso_num,hdataiso_deno,1,1,"B");

  // hmcid->Divide(hmcid_num,hmcid_deno,1,1);
  // hmciso->Divide(hmciso_num,hmciso_deno,1,1);
  // hdataid->Divide(hdataid_num,hdataid_deno,1,1);
  // hdataiso->Divide(hdataiso_num,hdataiso_deno,1,1);
	
  // get scale factors
  hsfid->Divide(hdataid, hmcid, 1, 1);
  hsfiso->Divide(hdataiso, hmciso, 1, 1);

  // Draw histograms	
  //hmcid->Draw("text");


  printHeader( leptype , "MC ID" );
  printline(hmcid);

  printHeader( leptype , "MC ISO" );
  printline(hmciso);

  printHeader( leptype , "DATA ID" );
  printline(hdataid);

  printHeader( leptype , "DATA ISO" );
  printline(hdataiso);

  printHeader( leptype , "Scale Factor ID" );
  printline(hsfid);

  printHeader( leptype , "Scale Factor ISO" );
  printline(hsfiso);

  cout << "\\hline" << endl << "\\hline" << endl;

  TCanvas *c_iso[10];
  TCanvas *c_id[10];

  for( int i = 0 ; i < 1 ; i++ ){

    TCut mysel;
    if     ( i==0 ) mysel = TCut(tnpcut+njets0);
    else if( i==1 ) mysel = TCut(tnpcut+njets1);
    else if( i==2 ) mysel = TCut(tnpcut+njets2);
    else if( i==3 ) mysel = TCut(tnpcut+njets3);
    else if( i==4 ) mysel = TCut(tnpcut+njets4);

    c_iso[i] = new TCanvas(Form("c_iso_%i",i),Form("c_iso_%i",i),600,600);
    c_iso[i]->cd();
    printHisto( c_iso[i] , chdata , chmc , TCut(lepiso) , TCut(mysel+lepid) , "probe.pt()" , 10 , 0.0 , 340.0 , "lepton p_{T} [GeV]" , "iso efficiency" );
    if( printplot ) c_iso[i]->Print(Form("plots/%s_iso_njets%i.pdf",lepchar,i));

    c_id[i] = new TCanvas(Form("c_id_%i",i),Form("c_id_%i",i),600,600);
    c_id[i]->cd();
    printHisto( c_id[i] , chdata , chmc , TCut(lepid) , TCut(mysel+lepiso) , "probe.pt()" , 10 , 0.0 , 340.0 , "lepton p_{T} [GeV]" , "ID efficiency" );
    if( printplot ) c_id[i]->Print(Form("plots/%s_id_njets%i.pdf",lepchar,i));

  }

  /*

  //---------------------------
  // tag cuts
  //---------------------------

  TCut zmass("abs(tagAndProbeMass-91)<15");
  TCut eltnp("(eventSelection&1)==1");
  TCut mutnp("(eventSelection&2)==2");
  TCut os("qProbe*qTag<0");
  TCut tag_eta21("abs(tag->eta())<2.1");
  TCut tag_eta25("abs(tag->eta())<2.5");
  TCut njets1("njets>=1");
  TCut njets2("njets>=2");
  TCut njets3("njets>=3");
  TCut njets4("njets>=4");
  TCut tag_pt30("tag->pt()>30.0");
  TCut met30("met<30");
  TCut met20("met<20");
  TCut nbm0("nbm==0");
  TCut nbl0("nbl==0");
  TCut mt30("mt<30");
  TCut eltnptrig("HLT_TNP_tag > 0 || HLT_TNPel_tag > 0");
  TCut mutnptrig("HLT_IsoMu30_eta2p1_tag > 0");

  //---------------------------
  // tag cuts
  //---------------------------

  TCut mufo 	= "(leptonSelection&32768)==32768";    // mu fo
  TCut elfo     = "(leptonSelection&4)==4";            // ele fo 
  TCut elid  	= "(leptonSelection&8)==8";            // ele id 
  TCut eliso 	= "(leptonSelection&16)==16";          // ele iso
  TCut probept  = "probe->pt()>30";                    // probe pt
  TCut drprobe  = "drprobe<0.05";                      // dR(probe,pfcandidate)

  TCut eltnpcut;
  eltnpcut += zmass;
  eltnpcut += os;
  eltnpcut += eltnp;
  eltnpcut += tag_eta25;
  //eltnpcut += njets2;
  eltnpcut += tag_pt30;
  eltnpcut += eltnptrig;
  eltnpcut += met30;
  // eltnpcut += mt30;
  eltnpcut += nbl0;
  
  eltnpcut += elid;
  eltnpcut += probept;
  eltnpcut += drprobe;

  TCut mutnpcut;
  mutnpcut += zmass;
  mutnpcut += os;
  mutnpcut += mutnp;
  mutnpcut += tag_eta21;
  //mutnpcut += njets2;
  mutnpcut += tag_pt30;
  mutnpcut += mutnptrig;
  mutnpcut += met30;
  // mutnpcut += mt30;
  mutnpcut += nbl0;

  mutnpcut += muid;
  mutnpcut += probept;
  mutnpcut += drprobe;


  //eltnpcut += njets2;
  //eltnpcut += njets3;
  //eltnpcut += nbm0;
  //eltnpcut += mt30;
  //eltnpcut += met20;

  //TCut eltnpcut 	 = "abs(tagAndProbeMass-91)<15 && (eventSelection&1)==1 && qProbe*qTag<0 && abs(tag->eta())<2.5 && njets>=4 && tag->pt()>30.0 && met<30.0 && nbm==0 && mt<30"; 
  //TCut mutnpcut 	 = "abs(tagAndProbeMass-91)<15 && (eventSelection&2)==2 && HLT_IsoMu30_eta2p1_tag>0 && qProbe*qTag<0 && abs(tag->eta())<2.1 && njets>=4 && tag->pt()>30.0"; 

  TCut vtxweight = "vtxweight";

  cout << "Electrons:" << endl;
  cout << "Total MC yields 	: " << chmc->GetEntries(eltnpcut) << endl;
  cout << "Total DATA yields 	: " << chdata->GetEntries(eltnpcut) << endl;

  cout << "Muons:" << endl;
  cout << "Total MC yields 	: " << chmc->GetEntries(mutnpcut) << endl;
  cout << "Total DATA yields 	: " << chdata->GetEntries(mutnpcut) << endl;


  //TCut njets    = "njets>=2";
  TCut tkisoold = "tkisoold/probe->pt()>0.1";
  TCut tkisonew = "tkisonew/probe->pt()>0.1";

  //-----------------------------------------
  // check nvtx data vs. MC
  //-----------------------------------------

  TH1F *hnvtx_mc   = new TH1F("hnvtx_mc"  ,"",30,0,30);
  TH1F *hnvtx_data = new TH1F("hnvtx_data","",30,0,30);

  hnvtx_mc->Sumw2();
  hnvtx_data->Sumw2();

  chdata->Draw("nvtx>>hnvtx_data",(eltnpcut||mutnpcut));
  chmc->Draw("nvtx>>hnvtx_mc",(eltnpcut||mutnpcut)*vtxweight);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  hnvtx_mc->SetLineColor(2);
  hnvtx_mc->SetMarkerColor(2);
  hnvtx_data->SetLineColor(4);
  hnvtx_data->SetMarkerColor(4);
  
  hnvtx_data->GetXaxis()->SetTitle("N_{VTX}");
  hnvtx_data->DrawNormalized();
  hnvtx_mc->DrawNormalized("same");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(hnvtx_data,"data","lp");
  leg->AddEntry(hnvtx_mc,"MC","lp");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  if( printplot ) c1->Print("plots/nvtx.pdf");



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

  */
	
}
