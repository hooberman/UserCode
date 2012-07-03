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

  float ratio    = effdata/effmc;
  float ratioerr = ratio * sqrt(pow(effdataerr/effdata,2)+pow(effmcerr/effmc,2));

  float datatot = hdata->Integral();
  float mctot   = hmc->Integral();
  
  cout << endl;
  cout << plottitle << endl;
  cout << "Data eff  " << Form("%.2f +/- %.3f",effdata,effdataerr) << endl;
  cout << "MC   eff  " << Form("%.2f +/- %.3f",effmc  ,effmcerr)   << endl;
  cout << "Data/MC   " << Form("%.2f +/- %.2f",ratio  ,ratioerr)   << endl;

  //--------------------------------------
  // draw stuff
  //--------------------------------------

  hdata->Scale(1.0/hdata->Integral());
  hmc->Scale(1.0/hmc->Integral());

  if( log ) hmc->GetYaxis()->SetRangeUser(0.001,3);  
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

    TH1F* hratio = (TH1F*) hdata->Clone(Form("hratio_%i",iplot));
    hratio->Divide(hmc);

    hratio->SetMarkerColor(1);
    hratio->SetLineColor(1);
    hratio->Draw();
    hratio->GetYaxis()->SetRangeUser(0,2);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetLabelSize(0.2);
    hratio->GetXaxis()->SetLabelSize(0.0);
  
    TLine line;
    line.DrawLine(xmin,1.0,xmax,1.0);
  }
  
  //data->Scan("run:lumi:event:probe->pt():probe->eta():tkisonew:met:mt:njets:nbl:nbm",sel+"tkisonew>20");

  if( printplot ) can->Print(Form("plots/%s.pdf",plottitle));

  iplot++;

  // TCanvas *c2 = new TCanvas();
  // c2->cd();
  // grdata->Draw("AP");

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


void tnpScale( bool printplot = false ) {

  //----------------------------------------
  // Files
  //----------------------------------------

  char* version = (char*) "V00-00-05";

  TChain *chmc   = new TChain("leptons");
  TChain *chdata = new TChain("leptons");

  //char* suffix = "";
  char* suffix = "_2jets";

  chmc->Add(Form("smurf/%s/dymm_testskim%s.root" , version , suffix));
  
  chdata->Add(Form("smurf/%s/data_DoubleElectron_May10%s.root"    , version , suffix));
  chdata->Add(Form("smurf/%s/data_DoubleElectron_PRv4%s.root"     , version , suffix));
  chdata->Add(Form("smurf/%s/data_DoubleElectron_PRv6%s.root"     , version , suffix));
  chdata->Add(Form("smurf/%s/data_DoubleElectron_Aug05%s.root"    , version , suffix));
  chdata->Add(Form("smurf/%s/data_DoubleElectron_B30%s.root"      , version , suffix));
  chdata->Add(Form("smurf/%s/data_DoubleElectron_B34%s.root"      , version , suffix));

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
  float ptbin[] = {10., 15., 20., 30., 40., 50., 60. , 70. , 80.0 , 100.0 , 7000.};
  float etabin[] = {0, 0.8, 1.479, 2.0, 2.5};
  int nptbin=10;
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


  //---------------------------
  // tag cuts
  //---------------------------

  TCut zmass("abs(tagAndProbeMass-91)<15");
  TCut eltnp("(eventSelection&1)==1");
  TCut mutnp("(eventSelection&2)==2");
  TCut os("qProbe*qTag<0");
  TCut tag_eta21("abs(tag->eta())<2.1");
  TCut tag_eta25("abs(tag->eta())<2.5");
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
  TCut muid 	= "(leptonSelection&65536)==65536";    // mu id 
  TCut muiso 	= "(leptonSelection&131072)==131072";  // mu iso 
  TCut elfo     = "(leptonSelection&4)==4";            // ele fo 
  TCut elid  	= "(leptonSelection&8)==8";            // ele id 
  TCut eliso 	= "(leptonSelection&16)==16";          // ele iso
  TCut probept  = "probe->pt()>30";                    // probe pt

  TCut eltnpcut;
  eltnpcut += zmass;
  eltnpcut += os;
  eltnpcut += eltnp;
  eltnpcut += tag_eta25;
  eltnpcut += njets2;
  eltnpcut += tag_pt30;
  eltnpcut += eltnptrig;
  eltnpcut += met30;
  eltnpcut += nbl0;

  eltnpcut += elid;
  eltnpcut += probept;

  TCut mutnpcut;
  mutnpcut += zmass;
  mutnpcut += os;
  mutnpcut += mutnp;
  mutnpcut += tag_eta21;
  mutnpcut += njets2;
  mutnpcut += tag_pt30;
  mutnpcut += mutnptrig;
  //mutnpcut += met30;
  //mutnpcut += nbl0;

  mutnpcut += muid;
  mutnpcut += probept;

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

  //--------------------------
  // absolute track isolation
  //--------------------------
  
  bool residual = true;
  bool log      = true;

  plotDistribution( chdata , chmc , TCut(eltnpcut)        , vtxweight , "tkisonew" , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_2j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut)        , vtxweight , "tkisonew" , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_2j" , printplot , residual , log );

  plotDistribution( chdata , chmc , TCut(eltnpcut+njets3) , vtxweight , "tkisonew" , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_3j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut+njets3) , vtxweight , "tkisonew" , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_3j" , printplot , residual , log );

  plotDistribution( chdata , chmc , TCut(eltnpcut+njets4) , vtxweight , "tkisonew" , 10 , 0 , 10 , "abs tkiso [GeV]"    , "el_tkiso_4j" , printplot , residual , log );
  plotDistribution( chdata , chmc , TCut(mutnpcut+njets4) , vtxweight , "tkisonew" , 10 , 0 , 10 , "abs tkiso [GeV]"    , "mu_tkiso_4j" , printplot , residual , log );

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
