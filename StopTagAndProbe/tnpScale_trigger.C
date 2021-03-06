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
  for (int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {

    Float_t min = h2->GetXaxis()->GetBinLowEdge(x);
    Float_t max = min + h2->GetXaxis()->GetBinWidth(x);
    printf("  %4.1f - %4.1f  & ", min, max);

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

float getBinomialError( float num , float den ){

  TGraphAsymmErrors* grtemp = new TGraphAsymmErrors();

  TH1F* hnum = new TH1F("hnum","",1,0,1);
  TH1F* hden = new TH1F("hden","",1,0,1);

  hnum->SetBinContent(1,num);
  hden->SetBinContent(1,den);

  grtemp->BayesDivide(hnum,hden);

  float err = 0.5 * ( grtemp->GetErrorYlow(0) + grtemp->GetErrorYhigh(0) );

  delete hnum;
  delete hden;

  return err;
}


void plotDistribution( TChain* data , TChain *mc , TCut sel , TCut vtxweight , char* var , int nbins , float xmin , float xmax , char* xtitle , char* plottitle = "" , bool printplot = false , bool residual = false , bool log = false ){

  //--------------------------------------
  // define histograms and TGraphs
  //--------------------------------------

  TH1F* hdata      = new TH1F(Form("hdata_%i"     , iplot),Form("hdata_%i"    , iplot),nbins,xmin,xmax);
  TH1F* hmc        = new TH1F(Form("hmc_%i"       , iplot),Form("hmc_%i"      , iplot),nbins,xmin,xmax);
  TH1F* hmc_novtx  = new TH1F(Form("hmc_novtx_%i" , iplot),Form("hmc_novtx%i" , iplot),nbins,xmin,xmax);

  hdata->Sumw2();
  hmc->Sumw2();

  TGraphAsymmErrors* grdata = new TGraphAsymmErrors();
  TGraphAsymmErrors* grmc   = new TGraphAsymmErrors();

  TH1F* hdata_denom = new TH1F(Form("hdata_denom_%i",iplot),"",nbins,xmin,xmax);
  TH1F* hmc_denom   = new TH1F(Form("hmc_denom_%i"  ,iplot),"",nbins,xmin,xmax);

  //--------------------------------------
  // set up canvas and pads
  //--------------------------------------

  TCanvas *can = new TCanvas(Form("can_%i",iplot),Form("can_%i",iplot),600,600);
  can->cd();
  if( log ) gPad->SetLogy();

  TPad *mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);

  if( residual ){
    mainpad->Draw();
    mainpad->cd();
    if( log ) mainpad->SetLogy();
  }

  //--------------------------------------
  // fill histos and TGraphs
  //--------------------------------------

  data->Draw(Form("min(%s,%f)>>hdata_%i"     , var,xmax-0.0001,iplot),sel);
  mc  ->Draw(Form("min(%s,%f)>>hmc_%i"       , var,xmax-0.0001,iplot),sel*vtxweight);
  mc  ->Draw(Form("min(%s,%f)>>hmc_novtx_%i" , var,xmax-0.0001,iplot),sel);

  for( int ibin = 1 ; ibin <= nbins ; ibin++ ){
    hdata_denom->SetBinContent(ibin,hdata->Integral());
    hmc_denom->SetBinContent(ibin,hmc->Integral());
  }

  grdata->BayesDivide(hdata,hdata_denom);
  grmc->BayesDivide(hmc_novtx,hmc_denom);

  //--------------------------------------
  // get efficiencies and errors
  //--------------------------------------

  /*
  float ndata1     = (float) hdata->GetBinContent(1);
  float ndata      = (float) hdata->Integral();
  float effdata    = 1-ndata1 / ndata;

  // TGraphAsymmErrors* grdata_temp = new TGraphAsymmErrors();
  // TH1F* hdata_num_temp = new TH1F(Form("hdata_num_temp_%i",iplot),"",1,0,1);
  // TH1F* hdata_den_temp = new TH1F(Form("hdata_den_temp_%i",iplot),"",1,0,1);
  // hdata_num_temp->SetBinContent(1,ndata-ndata1);
  // hdata_den_temp->SetBinContent(1,ndata);
  // grdata_temp->BayesDivide(hdata_num_temp,hdata_den_temp);

  //float effdataerr = sqrt(ndata1) / ndata;
  float effdataerr = 0.5 * ( grdata->GetErrorYlow(0) + grdata->GetErrorYhigh(0) );
  //float effdataerr = 0.5 * ( grdata_temp->GetErrorYlow(0) + grdata_temp->GetErrorYhigh(0) );

  float nmc1       = (float) hmc->GetBinContent(1);
  float nmc        = (float) hmc->Integral();
  float effmc      = 1-nmc1 / nmc;
  //float effmcerr   = hmc->GetBinError(1) / nmc;
  float effmcerr   = 0.5 * ( grmc->GetErrorYlow(0) + grmc->GetErrorYhigh(0) );


  float datatot = hdata->Integral();
  float mctot   = hmc->Integral();
  
  cout << endl;
  cout << plottitle << endl;

  cout << "Data eff  " << Form("%.2f +/- %.3f",effdata,effdataerr) << endl;
  cout << "MC   eff  " << Form("%.2f +/- %.3f",effmc  ,effmcerr)   << endl;
  cout << "Data/MC   " << Form("%.2f +/- %.2f",ratio  ,ratioerr)   << endl;
  */

  float ndata    = hdata->Integral();
  float ndata1   = hdata->Integral(2,20);
  float ndata2   = hdata->Integral(3,20);
  float ndata3   = hdata->Integral(4,20);
  float ndata4   = hdata->Integral(5,20);
  float ndata5   = hdata->Integral(6,20);

  float nmc      = hmc->Integral();
  float nmc1     = hmc->Integral(2,20);
  float nmc2     = hmc->Integral(3,20);
  float nmc3     = hmc->Integral(4,20);
  float nmc4     = hmc->Integral(5,20);
  float nmc5     = hmc->Integral(6,20);

  float effdata1 = ndata1/ndata;
  float effdata2 = ndata2/ndata;
  float effdata3 = ndata3/ndata;
  float effdata4 = ndata4/ndata;
  float effdata5 = ndata5/ndata;

  float effmc1   = nmc1/nmc;
  float effmc2   = nmc2/nmc;
  float effmc3   = nmc3/nmc;
  float effmc4   = nmc4/nmc;
  float effmc5   = nmc5/nmc;

  float effdata1err = getBinomialError(ndata1,ndata);
  float effdata2err = getBinomialError(ndata2,ndata);
  float effdata3err = getBinomialError(ndata3,ndata);
  float effdata4err = getBinomialError(ndata4,ndata);
  float effdata5err = getBinomialError(ndata5,ndata);

  float effmc1err   = getBinomialError(nmc1,nmc);
  float effmc2err   = getBinomialError(nmc2,nmc);
  float effmc3err   = getBinomialError(nmc3,nmc);
  float effmc4err   = getBinomialError(nmc4,nmc);
  float effmc5err   = getBinomialError(nmc5,nmc);

  float ratio1      = effdata1/effmc1;
  float ratio2      = effdata2/effmc2;
  float ratio3      = effdata3/effmc3;
  float ratio4      = effdata4/effmc4;
  float ratio5      = effdata5/effmc5;

  float ratio1err   = ratio1 * sqrt(pow(effdata1err/effdata1,2)+pow(effmc1err/effmc1,2));
  float ratio2err   = ratio2 * sqrt(pow(effdata2err/effdata2,2)+pow(effmc2err/effmc2,2));
  float ratio3err   = ratio3 * sqrt(pow(effdata3err/effdata3,2)+pow(effmc3err/effmc3,2));
  float ratio4err   = ratio4 * sqrt(pow(effdata4err/effdata4,2)+pow(effmc4err/effmc4,2));
  float ratio5err   = ratio5 * sqrt(pow(effdata5err/effdata5,2)+pow(effmc5err/effmc5,2));

  cout << endl << endl << plottitle << endl;

  int left = 20;


  // char* delimstart = "|";
  // char* delim      = "|";
  // char* delimend   = "|";
  // char* pm         = "+/-";

  char* delimstart = "";
  char* delim      = "&";
  char* delimend   = "\\\\";
  char* pm         = "$\\pm$";

  cout << delimstart << setw(10) << "" << setw(4)
       << delim << setw(left) << "$>$ 1 GeV" << setw(4)
       << delim << setw(left) << "$>$ 2 GeV" << setw(4)
       << delim << setw(left) << "$>$ 3 GeV" << setw(4) 
       << delim << setw(left) << "$>$ 4 GeV" << setw(4)
       << delim << setw(left) << "$>$ 5 GeV" << setw(4) 
       << delimend << endl;

  cout << delimstart << setw(10) << "data" << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effdata1,pm,effdata1err) << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effdata2,pm,effdata2err) << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effdata3,pm,effdata3err) << setw(4) 
       << delim << setw(left) << Form("%.3f %s %.4f",effdata4,pm,effdata4err) << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effdata5,pm,effdata5err) << setw(4) 
       << delimend << endl;

  cout << delimstart << setw(10) << "mc" << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effmc1,pm,effmc1err) << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effmc2,pm,effmc2err) << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effmc3,pm,effmc3err) << setw(4) 
       << delim << setw(left) << Form("%.3f %s %.4f",effmc4,pm,effmc4err) << setw(4)
       << delim << setw(left) << Form("%.3f %s %.4f",effmc5,pm,effmc5err) << setw(4) 
       << delimend << endl;

  cout << delimstart << setw(10) << "data/mc" << setw(4)
       << delim << setw(left) << Form("%.2f %s %.2f",ratio1,pm,ratio1err) << setw(4)
       << delim << setw(left) << Form("%.2f %s %.2f",ratio2,pm,ratio2err) << setw(4)
       << delim << setw(left) << Form("%.2f %s %.2f",ratio3,pm,ratio3err) << setw(4) 
       << delim << setw(left) << Form("%.2f %s %.2f",ratio4,pm,ratio4err) << setw(4)
       << delim << setw(left) << Form("%.2f %s %.2f",ratio5,pm,ratio5err) << setw(4) 
       << delimend << endl;

  //--------------------------------------
  // draw stuff
  //--------------------------------------

  hdata->Scale(1.0/hdata->Integral());
  hmc->Scale(1.0/hmc->Integral());

  if( log ) hmc->GetYaxis()->SetRangeUser(0.0001,5);  
  else      hmc->GetYaxis()->SetRangeUser(0.0,1);  

  hmc->GetXaxis()->SetTitle(xtitle);
  hmc->SetLineColor(2);
  hmc->SetMarkerColor(2);
  hmc->DrawNormalized("hist");
  hmc->DrawNormalized("sameE1");
  hdata->SetLineColor(4);
  hdata->SetMarkerColor(4);
  hdata->Draw("sameE1");

  grdata->SetLineColor(6);
  grmc->SetLineColor(7);
  //grdata->Draw("sameP");
  //grmc->Draw("sameP");

  TLegend *leg = new TLegend(0.6,0.7,0.8,0.9);
  leg->AddEntry(hdata , "data" , "lp");
  leg->AddEntry(hmc   , "MC"   , "lp");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);			       
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();

  if( TString(plottitle).Contains("el") ) t->DrawLatex(0.6,0.6,"electrons");
  if( TString(plottitle).Contains("mu") ) t->DrawLatex(0.6,0.6,"muons");

  if( TString(plottitle).Contains("0j") ) t->DrawLatex(0.6,0.5,"n_{jets} #geq 0");
  if( TString(plottitle).Contains("1j") ) t->DrawLatex(0.6,0.5,"n_{jets} #geq 1");
  if( TString(plottitle).Contains("2j") ) t->DrawLatex(0.6,0.5,"n_{jets} #geq 2");
  if( TString(plottitle).Contains("3j") ) t->DrawLatex(0.6,0.5,"n_{jets} #geq 3");
  if( TString(plottitle).Contains("4j") ) t->DrawLatex(0.6,0.5,"n_{jets} #geq 4");

  //--------------------------------------
  // draw residual plots
  //--------------------------------------

  if( residual ){
    can->cd();
  
    TPad *respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
    respad->Draw();
    respad->cd();
    respad->SetGridy();

    TH1F* hratio = (TH1F*) hdata->Clone(Form("hratio_%i",iplot));
    hratio->Divide(hmc);

    hratio->SetMarkerColor(1);
    hratio->SetLineColor(1);
    hratio->Draw();
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetLabelSize(0.2);
    hratio->GetXaxis()->SetLabelSize(0.0);
  
    TLine line;
    line.DrawLine(xmin,1.0,xmax,1.0);
  }
  
  //data->Scan("run:lumi:event:probe->pt():probe->eta():tkisonew:met:mt:njets:nbl:nbm",sel+"tkisonew>20");
  //data->Scan("run:lumi:event:probe->pt():probe->eta():tkisonew:met:mt:njets:nbl:nbm",sel);

  if( printplot ) can->Print(Form("plots/%s.pdf",plottitle));

  iplot++;

  // TCanvas *c2 = new TCanvas();
  // c2->cd();
  // grdata->Draw("AP");

}

TGraphAsymmErrors* getEfficiencyGraph( TChain* ch , TCut num , TCut denom , char* var , int nbins , float xmin , float xmax , char* xtitle , char* ytitle){


  //float ptbin[] = {10., 15., 20., 30., 40., 50., 7000.};
  float ptbin[] = { 20.0 , 25.0 , 30.0 , 35.0 , 40.0 , 45.0 , 50.0 , 55.0 , 60.0 , 70.0 , 80.0 , 90.0 , 100.0 , 150.0 , 200.0 , 250.0 , 300.0};
  int   nptbin  = 16;

  TH1F* hpass   = new TH1F(Form("hpass_%i",iplot),Form("hpass_%i",iplot),nptbin,ptbin);
  TH1F* hall    = new TH1F(Form("hall_%i" ,iplot),Form("hall_%i" ,iplot),nptbin,ptbin);

  //TH1F* hpass   = new TH1F(Form("hpass_%i",iplot),Form("hpass_%i",iplot),nbins,xmin,xmax);
  //TH1F* hall    = new TH1F(Form("hall_%i" ,iplot),Form("hall_%i" ,iplot),nbins,xmin,xmax);

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();  
  ch->Draw(Form("min(%s,%f)>>hpass_%i"  , var,xmax-0.0001,iplot),denom+num);
  ch->Draw(Form("min(%s,%f)>>hall_%i"   , var,xmax-0.0001,iplot),denom);
  delete ctemp;

  /*
  cout << "PASS " << hpass->Integral() << endl;
  cout << "ALL  " << hall->Integral() << endl;
  cout << "EFF  " << (float) hpass->Integral() / (float) hall->Integral() << endl;

  TCanvas *c1 = new TCanvas();
  c1->cd();
  hall->Draw();
  hpass->SetLineColor(2);
  hpass->Draw("samehist");
  */

  TGraphAsymmErrors *gr = new TGraphAsymmErrors();
  gr->BayesDivide(hpass,hall);

  gr->GetXaxis()->SetTitle(xtitle);
  gr->GetYaxis()->SetTitle(ytitle);

  iplot++;

  return gr;
}

void printHisto( TChain *data , TChain *mc , TCut num , TCut denom , char* var , int nbins , float xmin , float xmax , char* xtitle , char* ytitle){

  TH1F* hpass_data = new TH1F(Form("hpass_data_%i",iplot),Form("hpass_data_%i",iplot),nbins,xmin,xmax);
  TH1F* hall_data  = new TH1F(Form("hall_data_%i" ,iplot),Form("hall_data_%i" ,iplot),nbins,xmin,xmax);
  TH1F* hpass_mc   = new TH1F(Form("hpass_mc_%i"  ,iplot),Form("hpass_mc_%i"  ,iplot),nbins,xmin,xmax);
  TH1F* hall_mc    = new TH1F(Form("hall_mc_%i"   ,iplot),Form("hall_mc_%i"   ,iplot),nbins,xmin,xmax);

  TCanvas *can = new TCanvas(Form("can_%i",iplot),Form("can_%i",iplot),600,600);
  can->cd();
  
  data->Draw(Form("min(%s,%f)>>hpass_data_%i"  , var,xmax-0.0001,iplot),denom+num);
  data->Draw(Form("min(%s,%f)>>hall_data_%i"   , var,xmax-0.0001,iplot),denom);
  mc->Draw  (Form("min(%s,%f)>>hpass_mc_%i"    , var,xmax-0.0001,iplot),denom+num);
  mc->Draw  (Form("min(%s,%f)>>hall_mc_%i"     , var,xmax-0.0001,iplot),denom);

  TGraphAsymmErrors *grdata = new TGraphAsymmErrors();
  grdata->BayesDivide(hpass_data,hall_data);

  TGraphAsymmErrors *grmc = new TGraphAsymmErrors();
  grmc->BayesDivide(hpass_mc,hall_mc);

  cout << "data all  " << hall_data->GetBinContent(8) << endl;
  cout << "data pass " << hpass_data->GetBinContent(8) << endl;
  cout << "data eff  " << hpass_data->GetBinContent(8) / hall_data->GetBinContent(8) << endl;

  Double_t x;
  Double_t y;
  grdata->GetPoint(7,x,y);
  cout << "data eff2 " << y << endl;

  gPad->SetGridx();
  gPad->SetGridy();

  grdata->SetMarkerColor(2);
  grdata->SetLineColor(2);
  grmc->SetMarkerColor(4);
  grmc->SetLineColor(4);
  grmc->SetMarkerStyle(25);

  grdata->GetXaxis()->SetTitle(xtitle);
  grdata->GetYaxis()->SetTitle(ytitle);
  grdata->Draw("AP");
  grmc->Draw("sameP");

  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);
  leg->AddEntry(grdata ,"data","lp");
  leg->AddEntry(grmc   ,"mc"  ,"lp");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  iplot ++;
}


void tnpScale_trigger( bool printplot = false ) {

  //----------------------------------------
  // Files
  //----------------------------------------

  char* version = (char*) "V00-00-07";

  TChain *chmc   = new TChain("leptons");
  TChain *chdata = new TChain("leptons");

  //char* suffix = "";
  char* suffix = "_2jets";

  //chmc->Add(Form("smurf/%s/dymm_test%s.root" , version , suffix));
  //chmc->Add(Form("smurf/%s/dymm_test%s_INCOMPLETE.root" , version , suffix));
  //chmc->Add(Form("smurf/%s/dymm_testskim%s.root" , version , suffix));
  //chmc->Add(Form("smurf/%s/dymm_test%s.root" , version , suffix));
  
  // chdata->Add(Form("smurf/%s/data_DoubleElectron_May10%s.root"    , version , suffix));
  // chdata->Add(Form("smurf/%s/data_DoubleElectron_PRv4%s.root"     , version , suffix));
  // chdata->Add(Form("smurf/%s/data_DoubleElectron_PRv6%s.root"     , version , suffix));
  // chdata->Add(Form("smurf/%s/data_DoubleElectron_Aug05%s.root"    , version , suffix));
  // chdata->Add(Form("smurf/%s/data_DoubleElectron_B30%s.root"      , version , suffix));
  // chdata->Add(Form("smurf/%s/data_DoubleElectron_B34%s.root"      , version , suffix));

  chdata->Add(Form("smurf/%s/data_SingleMu_May10%s.root"          , version , suffix));
  chdata->Add(Form("smurf/%s/data_SingleMu_PRv4%s.root"           , version , suffix));
  chdata->Add(Form("smurf/%s/data_SingleMu_Aug05%s.root"          , version , suffix));
  chdata->Add(Form("smurf/%s/data_SingleMu_PRv6%s.root"           , version , suffix));
  chdata->Add(Form("smurf/%s/data_SingleMu_B30%s.root"            , version , suffix));
  chdata->Add(Form("smurf/%s/data_SingleMu_B34%s.root"            , version , suffix));

  //----------------------------------------
  // bins 
  //----------------------------------------

  //float ptbin[] = {10., 15., 20., 30., 40., 50., 7000.};
  float ptbin[] = {30., 40., 10000.};
  float etabin[] = {0, 0.8, 1.5, 2.1};
  int nptbin=2;
  int netabin=3;

  TH2F *hdataid_deno 	= new TH2F("hdataid_deno", "hdataid_deno", nptbin, ptbin, netabin, etabin);
  TH2F *hdataid_num 	= new TH2F("hdataid_num" , "hdataid_num" , nptbin, ptbin, netabin, etabin);
  TH2F *hdataid_eff 	= new TH2F("hdataid_eff" , "hdataid_eff" , nptbin, ptbin, netabin, etabin);


  hdataid_deno->Sumw2();
  hdataid_num->Sumw2();
  hdataid_eff->Sumw2();

  //
  // histogram
  //
  //deno
  // TH2F *hmcid_deno 	= new TH2F("hmcid_deno", "hmcid_deno", nptbin, ptbin, netabin, etabin);
  // TH2F *hmciso_deno 	= new TH2F("hmciso_deno", "hmciso_deno", nptbin, ptbin, netabin, etabin);
  // TH2F *hdataid_deno 	= new TH2F("hdataid_deno", "hdataid_deno", nptbin, ptbin, netabin, etabin);
  // TH2F *hdataiso_deno	= new TH2F("hdataiso_deno", "hdataiso_deno", nptbin, ptbin, netabin, etabin);
  // hmcid_deno->Sumw2();
  // hmciso_deno->Sumw2();
  // hdataid_deno->Sumw2();
  // hdataiso_deno->Sumw2();
  // //num
  // TH2F *hmcid_num 	= new TH2F("hmcid_num", "hmcid_num", nptbin, ptbin, netabin, etabin);
  // TH2F *hmciso_num 	= new TH2F("hmciso_num", "hmciso_num", nptbin, ptbin, netabin, etabin);
  // TH2F *hdataid_num 	= new TH2F("hdataid_num", "hdataid_num", nptbin, ptbin, netabin, etabin);
  // TH2F *hdataiso_num 	= new TH2F("hdataiso_num", "hdataiso_num", nptbin, ptbin, netabin, etabin);
  // hmcid_num->Sumw2();
  // hmciso_num->Sumw2();
  // hdataid_num->Sumw2();
  // hdataiso_num->Sumw2();
  // // eff
  // TH2F *hmcid 	= new TH2F("hmcid", "hmcid", nptbin, ptbin, netabin, etabin);
  // TH2F *hmciso 	= new TH2F("hmciso", "hmciso", nptbin, ptbin, netabin, etabin);
  // TH2F *hdataid 	= new TH2F("hdataid", "hdataid", nptbin, ptbin, netabin, etabin);
  // TH2F *hdataiso 	= new TH2F("hdataiso", "hdataiso", nptbin, ptbin, netabin, etabin);
  // hmcid->Sumw2();
  // hmciso->Sumw2();
  // hdataid->Sumw2();
  // hdataiso->Sumw2();
  // // SF
  // TH2F *hsfid 	= new TH2F("hsfid", "hsfid", nptbin, ptbin, netabin, etabin);
  // TH2F *hsfiso 	= new TH2F("hsfiso", "hsfiso", nptbin, ptbin, netabin, etabin);
  // hsfid->Sumw2();
  // hsfiso->Sumw2();


  //---------------------------
  // tag cuts
  //---------------------------

  TCut zmass("abs(tagAndProbeMass-91)<15");
  TCut eltnp("(eventSelection&1)==1");
  TCut mutnp("(eventSelection&2)==2");
  TCut os("qProbe*qTag<0");
  TCut tag_eta21("abs(tag->eta())<2.1");
  TCut probe_eta21("abs(probe->eta())<2.1");
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
  TCut tag_trig("HLT_IsoMu30_eta2p1_tag > 0");
  TCut probe_trig("HLT_IsoMu30_eta2p1_probe > 0");

  //---------------------------
  // tag cuts
  //---------------------------

  TCut mufo 	   = "(leptonSelection&32768)==32768";    // mu fo
  TCut muid 	   = "(leptonSelection&65536)==65536";    // mu id 
  TCut muiso 	   = "(leptonSelection&131072)==131072";  // mu iso 
  TCut elfo        = "(leptonSelection&4)==4";            // ele fo 
  TCut elid  	   = "(leptonSelection&8)==8";            // ele id 
  TCut eliso 	   = "(leptonSelection&16)==16";          // ele iso
  TCut probept     = "probe->pt()>30";                    // probe pt
  TCut probe_pt20  = "probe->pt()>20";                    // probe pt
  TCut drprobe     = "drprobe<0.05";                      // dR(probe,pfcandidate)

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
  //mutnpcut += njets2;

  mutnpcut += tag_trig;
  mutnpcut += tag_pt30;
  mutnpcut += tag_eta21;

  mutnpcut += probe_pt20;
  mutnpcut += probe_eta21;

  mutnpcut += met30;
  mutnpcut += nbl0;

  mutnpcut += muid;
  mutnpcut += muiso;

  //mutnpcut += probept;
  //mutnpcut += drprobe;

  // mutnpcut += mt30;

  //
  //eltnpcut += njets3;
  //eltnpcut += nbm0;
  //eltnpcut += mt30;
  //eltnpcut += met20;

  //TCut eltnpcut 	 = "abs(tagAndProbeMass-91)<15 && (eventSelection&1)==1 && qProbe*qTag<0 && abs(tag->eta())<2.5 && njets>=4 && tag->pt()>30.0 && met<30.0 && nbm==0 && mt<30"; 
  //TCut mutnpcut 	 = "abs(tagAndProbeMass-91)<15 && (eventSelection&2)==2 && HLT_IsoMu30_eta2p1_tag>0 && qProbe*qTag<0 && abs(tag->eta())<2.1 && njets>=4 && tag->pt()>30.0"; 

  //TCut vtxweight = "vtxweight";

  // cout << "Electrons:" << endl;
  // cout << "Total MC yields 	: " << chmc->GetEntries(eltnpcut) << endl;
  // cout << "Total DATA yields 	: " << chdata->GetEntries(eltnpcut) << endl;

  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_deno" , TCut(mutnpcut)            );
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_num"  , TCut(mutnpcut+probe_trig) );	

  hdataid_eff->Divide(hdataid_num,hdataid_deno,1,1);

  cout << endl << "Denominator" << endl;
  printline(hdataid_deno);
  cout << endl << "Numerator" << endl;
  printline(hdataid_num);
  cout << endl << "Efficiency" << endl;
  printline(hdataid_eff);
  
  cout << "Muons:" << endl;
  cout << "Selection            : " << mutnpcut.GetTitle() << endl;
  //cout << "Total MC yields 	: " << chmc->GetEntries(mutnpcut) << endl;
  cout << "Total DATA yields 	: " << chdata->GetEntries(mutnpcut) << endl;

  TCut etabin1("abs(probe->eta()) < 0.8");
  TCut etabin2("abs(probe->eta()) > 0.8 && abs(probe->eta())<1.5");
  TCut etabin3("abs(probe->eta()) > 1.5 && abs(probe->eta())<2.1");
  TCut etabin4("abs(probe->eta()) > 2.1 && abs(probe->eta())<2.4");

  TGraphAsymmErrors* gr   = getEfficiencyGraph( chdata , probe_trig , mutnpcut               , "probe->pt()" , 28,20,300,"probe p_{T} [GeV]","trigger efficiency");
  TGraphAsymmErrors* gr1  = getEfficiencyGraph( chdata , probe_trig , TCut(mutnpcut+etabin1) , "probe->pt()" , 28,20,300,"probe p_{T} [GeV]","trigger efficiency");
  TGraphAsymmErrors* gr2  = getEfficiencyGraph( chdata , probe_trig , TCut(mutnpcut+etabin2) , "probe->pt()" , 28,20,300,"probe p_{T} [GeV]","trigger efficiency");
  TGraphAsymmErrors* gr3  = getEfficiencyGraph( chdata , probe_trig , TCut(mutnpcut+etabin3) , "probe->pt()" , 28,20,300,"probe p_{T} [GeV]","trigger efficiency");
  TGraphAsymmErrors* gr4  = getEfficiencyGraph( chdata , probe_trig , TCut(mutnpcut+etabin4) , "probe->pt()" , 28,20,300,"probe p_{T} [GeV]","trigger efficiency");

  TCanvas *c1 = new TCanvas();
  c1->cd();

  //gr->Draw("AP");

  gr1->SetLineColor(1);
  gr2->SetLineColor(2);
  gr3->SetLineColor(4);
  gr4->SetLineColor(6);

  gr1->SetMarkerColor(1);
  gr2->SetMarkerColor(2);
  gr3->SetMarkerColor(4);
  gr4->SetMarkerColor(6);

  gr1->Draw("AP");
  gr2->Draw("sameP");
  gr3->Draw("sameP");
  //gr4->Draw("sameP");

  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);
  leg->AddEntry(gr1,"|#eta| < 0.8","lp");
  leg->AddEntry(gr2,"|#eta| 0.8 - 1.5","lp");
  leg->AddEntry(gr3,"|#eta| 1.5 - 2.1","lp");
  //leg->AddEntry(gr4,"|#eta| > 2.1");

  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  /*
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
  */

  //--------------------------
  // absolute track isolation
  //--------------------------
  /*  
  bool residual = true;
  bool log      = true;
  
  char* var = "tkisonewnoveto";
  
  plotDistribution( chdata , chmc , TCut(eltnpcut)        , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_0j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut)        , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_0j" , printplot , residual , log );
  
  plotDistribution( chdata , chmc , TCut(eltnpcut+njets1) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_1j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut+njets1) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_1j" , printplot , residual , log );
  
  plotDistribution( chdata , chmc , TCut(eltnpcut+njets2) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_2j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut+njets2) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_2j" , printplot , residual , log );

  plotDistribution( chdata , chmc , TCut(eltnpcut+njets3) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_3j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut+njets3) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_3j" , printplot , residual , log );

  plotDistribution( chdata , chmc , TCut(eltnpcut+njets4) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_4j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut+njets4) , vtxweight , var , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_4j" , printplot , residual , log );
  */
  //--------------------------
  // relative track isolation
  //--------------------------

  // plotDistribution( chdata , chmc , TCut(eltnpcut) , vtxweight , "tkisonew/probe->pt()" , 10 , 0 , 1 , "rel tkiso [GeV]"    , "el_tkrliso_2j" , printplot );
  // plotDistribution( chdata , chmc , TCut(mutnpcut) , vtxweight , "tkisonew/probe->pt()" , 10 , 0 , 1 , "rel tkiso [GeV]"    , "mu_tkrliso_2j" , printplot );

  // plotDistribution( chdata , chmc , TCut(eltnpcut+njets3) , vtxweight , "tkisonew/probe->pt()" , 10 , 0 , 1 , "rel tkiso [GeV]"    , "el_tkrliso_3j" , printplot );
  // plotDistribution( chdata , chmc , TCut(mutnpcut+njets3) , vtxweight , "tkisonew/probe->pt()" , 10 , 0 , 1 , "rel tkiso [GeV]"    , "mu_tkrliso_3j" , printplot );

  // plotDistribution( chdata , chmc , TCut(eltnpcut+njets4) , vtxweight , "tkisonew/probe->pt()" , 10 , 0 , 1 , "rel tkiso [GeV]"    , "el_tkrliso_4j" , printplot );
  // plotDistribution( chdata , chmc , TCut(mutnpcut+njets4) , vtxweight , "tkisonew/probe->pt()" , 10 , 0 , 1 , "rel tkiso [GeV]"    , "mu_tkrliso_4j" , printplot );




  // plotDistribution( chdata , chmc , TCut(eltnpcut) , vtxweight , "tkisoold" , 10 , 0 , 10 , "trkiso (old) [GeV]" , "el_tkiso_old" , printplot );
  // plotDistribution( chdata , chmc , TCut(mutnpcut) , vtxweight , "tkisoold" , 10 , 0 , 10 , "trkiso (old) [GeV]" , "mu_tkiso_old" , printplot );

  //printHisto( chdata , chmc , TCut(muiso) , TCut(tnpcut+muid) , "probe.pt()" , 10 , 0.0 , 100.0 , "lepton p_{T} [GeV]" , "iso efficiency" );
  //printHisto( chdata , chmc , tkisoold , TCut(tnpcut+muid+njets) , "probe.pt()" , 10 , 0.0 , 100.0 , "lepton p_{T} [GeV]" , "tkiso(old) efficiency" );
  //printHisto( chdata , chmc , tkisonew , TCut(tnpcut+muid+njets) , "probe.pt()" , 10 , 0.0 , 100.0 , "lepton p_{T} [GeV]" , "tkiso(new) efficiency" );


  /*


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

  TCut tnpcut   	 = "abs(tagAndProbeMass-91)<15 && (eventSelection&2)==2 && HLT_IsoMu30_eta2p1_tag>0 && qProbe*qTag<0 && abs(tag->eta())<2.1 && njets>=4 && tag->pt()>30.0"; 

  chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_deno", 	tnpcut+muiso,	      	"goff");
  chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_num", 	tnpcut+muiso+muid,	"goff");
  chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_deno", 	tnpcut+muid,	      	"goff");
  chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_num", 	tnpcut+muid+muiso,	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_deno", 	tnpcut+muiso,		"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_num", 	tnpcut+muiso+muid,	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_deno", 	tnpcut+muid,	       	"goff");
  chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_num", 	tnpcut+muid+muiso,	"goff");

  // get efficiencies 
  // hmcid->Divide(hmcid_num,hmcid_deno,1,1,"B");
  // hmciso->Divide(hmciso_num,hmciso_deno,1,1,"B");
  // hdataid->Divide(hdataid_num,hdataid_deno,1,1,"B");
  // hdataiso->Divide(hdataiso_num,hdataiso_deno,1,1,"B");

  hmcid->Divide(hmcid_num,hmcid_deno,1,1);
  hmciso->Divide(hmciso_num,hmciso_deno,1,1);
  hdataid->Divide(hdataid_num,hdataid_deno,1,1);
  hdataiso->Divide(hdataiso_num,hdataiso_deno,1,1);
	
  // get scale factors
  hsfid->Divide(hdataid, hmcid, 1, 1);
  hsfiso->Divide(hdataiso, hmciso, 1, 1);

  // Draw histograms	
  //hmcid->Draw("text");
  */

  /*
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
  */
	
}
