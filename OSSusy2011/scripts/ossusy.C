#include <algorithm>
#include <iostream>
#include <iomanip>
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
#include "Hootilities.h"
#include "ossusy.h"
//#include "histtools.h"

using namespace std;

bool printgif_           = false;
bool alreadyInitialized_ = false;

//--------------------------------------------------
// initialize data/MC samples
//--------------------------------------------------

void initialize(char* path){

  if( alreadyInitialized_ ){

    cout << "Resetting babies" << endl;

    data->Reset();
    ttall->Reset();
    ttpowheg->Reset();
    ttdil->Reset();
    ttotr->Reset();
    ttll->Reset();
    tttau->Reset();
    ttfake->Reset();
    zjets->Reset();
    dy->Reset();
    ww->Reset();
    wz->Reset();
    zz->Reset();
    t->Reset();
    wjets->Reset();
    LM0->Reset();
    LM1->Reset();
    sm->Reset();
    sm_LM1->Reset();

    mc.clear();
    mctex.clear();
    mclabels.clear();
  }

  else{

    data	= new TChain("t");
    ttall	= new TChain("t");
    ttpowheg	= new TChain("t");
    ttdil	= new TChain("t");
    ttotr	= new TChain("t");
    ttll	= new TChain("t");
    tttau	= new TChain("t");
    ttfake	= new TChain("t");
    zjets      	= new TChain("t");
    dy		= new TChain("t");
    ww		= new TChain("t");
    wz		= new TChain("t");
    zz		= new TChain("t");
    t		= new TChain("t");
    wjets	= new TChain("t");
    LM0 	= new TChain("t");
    LM1 	= new TChain("t");
    sm          = new TChain("t");
    sm_LM1      = new TChain("t");
  }

  cout << endl;
  cout << "Loading babies at       : " << path << endl;
  
  data->Add(Form("%s/data_smallTree.root",path));
  //data->Add(Form("%s/data_smallTree_allpr_nojson.root",path));
  ttall->Add(Form("%s/ttall_smallTree.root",path));
  ttpowheg->Add(Form("%s/ttpowheg_smallTree.root",path));
  ttdil->Add(Form("%s/ttdil_smallTree.root",path));
  ttotr->Add(Form("%s/ttotr_smallTree.root",path));
  ttdil->Add(Form("%s/ttall_smallTree.root",path));
  ttotr->Add(Form("%s/ttall_smallTree.root",path));
  ttll->Add(Form("%s/ttall_smallTree.root",path));
  tttau->Add(Form("%s/ttall_smallTree.root",path));
  ttfake->Add(Form("%s/ttall_smallTree.root",path));
  zjets->Add(Form("%s/Zjets_smallTree.root",path));
  dy->Add(Form("%s/DYtot_smallTree.root",path));
  ww->Add(Form("%s/ww_smallTree.root",path));
  wz->Add(Form("%s/wz_smallTree.root",path));
  zz->Add(Form("%s/zz_smallTree.root",path));
  t->Add(Form("%s/tW_smallTree.root",path));
  wjets->Add(Form("%s/wjetsMG_smallTree.root",path));
  LM0->Add(Form("%s/LM0_smallTree.root",path));
  LM1->Add(Form("%s/LM1_smallTree.root",path));

  mc.push_back(ttall);     mclabels.push_back("ttall");
  //mc.push_back(ttdil);     mclabels.push_back("ttdil");
  //mc.push_back(ttotr);     mclabels.push_back("ttotr");
  //mc.push_back(ttpowheg);  mclabels.push_back("ttpowheg");
  //mc.push_back(ttotr);    mclabels.push_back("ttotr");
  //mc.push_back(ttll);     mclabels.push_back("ttll");    mctex.push_back("$t\\bar{b}\\rightarrow\\ell^+\\ell^-$");
  //mc.push_back(tttau);    mclabels.push_back("tttau");   mctex.push_back("$t\\bar{b}\\rightarrow\\ell^{\\pm}\\tau^{\\mp}$");
  //mc.push_back(ttfake);   mclabels.push_back("ttfake");  mctex.push_back("$t\\bar{b}\\rightarrow$fake");
  mc.push_back(wjets);    mclabels.push_back("wjets");   mctex.push_back("$W^{\\pm}$+jets");
  // mc.push_back(zjets);    mclabels.push_back("zjets"); mctex.push_back("$Z^0$+jets");
  mc.push_back(dy);       mclabels.push_back("DY");      mctex.push_back("DY");
  mc.push_back(ww);       mclabels.push_back("WW");      mctex.push_back("W^+W^-");
  mc.push_back(wz);       mclabels.push_back("WZ");      mctex.push_back("W^{\\pm}Z^0");
  mc.push_back(zz);       mclabels.push_back("ZZ");      mctex.push_back("Z^0Z^0");
  mc.push_back(t);        mclabels.push_back("t");       mctex.push_back("single top");
  //mc.push_back(LM0);      mclabels.push_back("LM0");     mctex.push_back("LM0");
  //mc.push_back(LM1);      mclabels.push_back("LM1");     mctex.push_back("LM1");
  
  sm->Add(Form("%s/ttall_smallTree.root",path));
  //sm->Add(Form("%s/ttdil_smallTree.root",path));
  //sm->Add(Form("%s/ttotr_smallTree.root",path));
  sm->Add(Form("%s/DYtot_smallTree.root",path));
  sm->Add(Form("%s/ww_smallTree.root",path));
  sm->Add(Form("%s/wz_smallTree.root",path));
  sm->Add(Form("%s/zz_smallTree.root",path));
  sm->Add(Form("%s/tW_smallTree.root",path));
  sm->Add(Form("%s/wjets_smallTree.root",path));

  sm_LM1->Add(Form("%s/ttall_smallTree.root",path));
  //sm_LM1->Add(Form("%s/ttdil_smallTree.root",path));
  //sm_LM1->Add(Form("%s/ttotr_smallTree.root",path));
  sm_LM1->Add(Form("%s/DYtot_smallTree.root",path));
  sm_LM1->Add(Form("%s/ww_smallTree.root",path));
  sm_LM1->Add(Form("%s/wz_smallTree.root",path));
  sm_LM1->Add(Form("%s/zz_smallTree.root",path));
  sm_LM1->Add(Form("%s/tW_smallTree.root",path));
  sm_LM1->Add(Form("%s/wjets_smallTree.root",path));
  sm_LM1->Add(Form("%s/LM1_smallTree.root",path));  

  alreadyInitialized_ = true;
}

//------------------------------------------
// selection and weight to apply to babies
//------------------------------------------

TCut selection_TCut( bool highpt ){

  TCut njets2("npfjets >= 2");
  TCut njets3("npfjets >= 3");
  TCut njets4("npfjets >= 4");
  TCut zveto("passz == 0");
  TCut zpass("(leptype==0||leptype==1) && dilep.mass()>76 && dilep.mass()<106");
  TCut met50("pfmet > 50");
  TCut ht100("htpf > 100");
  TCut ht200("htpf > 200");
  TCut pt2010("lep1.pt()>20 && lep2.pt()>10");
  TCut pt1010("lep1.pt()>10 && lep2.pt()>10");
  TCut pt105 ("lep1.pt()>10 && lep2.pt()>5" );
  TCut pt55  ("lep1.pt()>5  && lep2.pt()>5" );
  TCut relnt015 ("isont1<0.15&&isont2<0.15");
  TCut isolep1_t10_015("(isont1*lep1.pt())/max(lep1.pt(),20)<0.15");
  TCut isolep2_t10_015("(isont2*lep2.pt())/max(lep2.pt(),20)<0.15");
  TCut iso_t10_015 = isolep1_t10_015 + isolep2_t10_015;
  TCut ll ("w1>0 && w2>0");
  TCut lep2tight("lep2.pt()>10 || isont2<0.15");
  TCut goodrun("json==1");

  TCut highptsel = zveto + njets2 + met50 + ht100 + pt2010 ;
  TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + !pt2010 + lep2tight;

  //TCut highptsel = zveto + "npfjetspv>1" + met50 + "htpfpv>100" + pt2010 ;
  //TCut highptsel = zveto + pt2010 + "dilmass>50";
  //TCut highptsel = "dilmass>81&&dilmass<101";
  //TCut highptsel = zveto + "njetsuncor>1" + met50 + "htuncor>100." + pt2010 ;
  //TCut highptsel = zveto + "njetsoffset>1" + met50 + "htoffset>100." + pt2010 ;
  //TCut highptsel = zveto + njets2 + met50 + ht100 + pt2010 + "ndavtx>4";
  //TCut highptsel = zveto + met50 + pt2010 ;
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + lep2tight;
  //TCut highptsel = "";
  //TCut lowptsel  = ht200 + "dilep.mass()>50";
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + !pt1010 + ll;
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt1010 + !pt2010;
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + !pt1010 + relnt015;
  //TCut highptsel = "leptype==0 && pfmet>50";
  //TCut highptsel = pt2010 + met50 + njets2 + ht200 + zveto;
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + !pt1010 + lep2tight;
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt1010 + !pt2010;
  //TCut highptsel  = met50+njets2+ht100+pt2010+zveto;

  if( highpt ){
    cout << "Using high pt selection : " << highptsel.GetTitle() << endl;
    return highptsel;
  }
  else{
    cout << "Using low pt selection  : " << lowptsel.GetTitle() << endl;
    return lowptsel;
  }

  return 0;
}

TCut weight_TCut(){

  TCut weight("weight * ndavtxweight");
  //TCut weight("weight*(186./43.)*1.17*ndavtxweight");
  //TCut weight("weight");
  //TCut weight("1");

  cout << "Using weight            : " << weight.GetTitle() << endl;
  return weight;
}

void lowptiso( char* path ){

  deleteHistos();
  
  initialize(path);

  TCut sel    = selection_TCut(false);
  TCut weight = weight_TCut();

  TCut relnt015 ("isont1<0.15&&isont2<0.15");  
  TCut isolep1_t10_015("(isont1*lep1.pt())/max(lep1.pt(),20)<0.15");
  TCut isolep2_t10_015("(isont2*lep2.pt())/max(lep2.pt(),20)<0.15");
  TCut ll ("(w1>0&&w2>0) && (w1<3&&w2<3)");
  TCut fake("!(w1>0&&w2>0)");

  TH1F* hll   = new TH1F("hll",  "",1,0,1);  hll->Sumw2();
  TH1F* hfake = new TH1F("hfake","",1,0,1);  hfake->Sumw2();

  const unsigned int n = 6;

  float signt[n];
  float bkgnt[n];
  float sigerrnt[n];
  float bkgerrnt[n];

  float sigt10[n];
  float bkgt10[n];
  float sigerrt10[n];
  float bkgerrt10[n];

  float sig[n];
  float bkg[n];
  float sigerr[n];
  float bkgerr[n];

  TCut selclone = sel;

  TChain* tt = ttall;

  for( unsigned int i = 0 ; i < n ; ++i ){
    
    float cut = 0.05 + i * 0.02;

    //-------------------------------
    // non-truncated
    //-------------------------------

    sel = selclone + TCut(Form("isont1<%.2f && isont2<%.2f",cut,cut));

    tt->Draw("0.5>>hll",  (sel+ll  )*weight);
    tt->Draw("0.5>>hfake",(sel+fake)*weight);
    
    signt[i]    = hll->GetBinContent(1);
    bkgnt[i]    = hfake->GetBinContent(1);
    sigerrnt[i] = hll->GetBinError(1);
    bkgerrnt[i] = hfake->GetBinError(1);

    //-------------------------------
    // truncated
    //-------------------------------

    sel = selclone + TCut(Form("iso1<%.2f && iso2<%.2f",cut,cut));

    tt->Draw("0.5>>hll",  (sel+ll  )*weight);
    tt->Draw("0.5>>hfake",(sel+fake)*weight);
    
    sig[i]    = hll->GetBinContent(1);
    bkg[i]    = hfake->GetBinContent(1);
    sigerr[i] = hll->GetBinError(1);
    bkgerr[i] = hfake->GetBinError(1);

    //-------------------------------
    // truncated
    //-------------------------------

    sel = selclone + TCut(Form("(isont1*lep1.pt())/max(lep1.pt(),10)<%.2f && (isont2*lep2.pt())/max(lep2.pt(),10)<%.2f",cut,cut));

    tt->Draw("0.5>>hll",  (sel+ll  )*weight);
    tt->Draw("0.5>>hfake",(sel+fake)*weight);
    
    sigt10[i]    = hll->GetBinContent(1);
    bkgt10[i]    = hfake->GetBinContent(1);
    sigerrt10[i] = hll->GetBinError(1);
    bkgerrt10[i] = hfake->GetBinError(1);


    //cout << "cut    : " << sel.GetTitle() << endl;
    //cout << "ll     : " << hll->GetBinContent(1)   << " +/- " << hll->GetBinError(1) << endl;
    //cout << "fake   : " << hfake->GetBinContent(1) << " +/- " << hfake->GetBinError(1) << endl;


    
  }

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  TGraphErrors* gr    = new TGraphErrors(n,sig,bkg,sigerr,bkgerr);
  TGraphErrors* grnt  = new TGraphErrors(n,signt,bkgnt,sigerrnt,bkgerrnt);
  TGraphErrors* grt10 = new TGraphErrors(n,sigt10,bkgt10,sigerrt10,bkgerrt10);

  grnt->SetLineColor(2);
  grnt->SetMarkerColor(2);
  grnt->SetMarkerStyle(25);

  grt10->SetLineColor(4);
  grt10->SetMarkerColor(4);
  grt10->SetMarkerStyle(23);

  gr->GetXaxis()->SetLimits(0,30);
  gr->GetYaxis()->SetLimits(0,30);
  gr->GetYaxis()->SetRangeUser(0,30);

  gr->Draw("AP");
  grnt->Draw("sameP");
  //grt10->Draw("sameP");

  gr->GetXaxis()->SetTitle("ttll yield/fb^{-1}");
  gr->GetYaxis()->SetTitle("ttfake yield/fb^{-1}");

  TLegend *leg = new TLegend(0.2,0.6,0.6,0.8);
  leg->AddEntry(gr,"denom = max(pt,20)","p");
  leg->AddEntry(grnt,"denom = pt","p");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();


}

void plotHist( TH1F* h1 , TH1F* h2 , char* leg1 , char* leg2 , char* xtitle ){

  gPad->SetLogy(1);

  h1->GetXaxis()->SetTitle( xtitle );
  h1->DrawNormalized("hist");
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->DrawNormalized("sameE1");

  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h1,leg1,"l");
  leg->AddEntry(h2,leg2,"p");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

}

void sigregion( char* path , bool printgif = false ){

  deleteHistos();
  
  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;
  
  initialize(path);
  
  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();
  
  float pfmetcut = 300;
  float htcut    = 300;

  TH2F* hyield = new TH2F("hyield","",14,0,1400,10,0,500);
  TH1F* h = new TH1F("h","",1,0,1);

  for( unsigned int imet = 0 ; imet < 10 ; imet++ ){
    for( unsigned int iht = 0 ; iht < 14 ; iht++ ){
      
      pfmetcut = imet  * 50;
      htcut    = iht   * 100;

      //ttall->Draw("0.5>>h",(sel+Form("pfmet>%.0f && htpf>%.0f",pfmetcut,htcut))*weight*"1000./90.");
      ttpowheg->Draw("0.5>>h",(sel+Form("pfmet>%.0f && htpf>%.0f",pfmetcut,htcut))*weight*"1000./90.");

      float yield = h->Integral();
      
      //int metbin = hyield->FindBin(pfmetcut);
      //int htbin  = hyield->FindBin(htcut);

      cout << endl;
      cout << "met>" << pfmetcut << " ht>" << htcut << " yield " << yield << endl;
      

      hyield->SetBinContent( iht+1 , imet+1 , yield );
    }
  }

  TCanvas *c1 = new TCanvas();
  c1->cd();

  gStyle->SetPaintTextFormat(".0f");

  hyield->GetXaxis()->SetTitle("H_{T} cut (GeV)");
  hyield->GetYaxis()->SetTitle("MET cut (GeV)");
  hyield->SetMaximum(2);
  hyield->Draw("colz");
  hyield->Draw("sametext");

}

//------------------------------------
// compare ttbar samples
//------------------------------------

void compareTT( char* path , bool printgif = false ){

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  
  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  TH1F* njets_mg  = getHist( ttall    , "npfjets" , TCut(sel*weight) , "njets_mg"  , 10 , 0 , 10 );
  TH1F* njets_pow = getHist( ttpowheg , "npfjets" , TCut(sel*weight) , "njets_pow" , 10 , 0 , 10 );

  TH1F* ht_mg     = getHist( ttall    , "htpf"    , TCut(sel*weight) , "ht_mg"     , 10 , 0 , 1000 );
  TH1F* ht_pow    = getHist( ttpowheg , "htpf"    , TCut(sel*weight) , "ht_pow"    , 10 , 0 , 1000 );

  TH1F* met_mg    = getHist( ttall    , "pfmet"   , TCut(sel*weight) , "met_mg"     , 10 , 0 , 500 );
  TH1F* met_pow   = getHist( ttpowheg , "pfmet"   , TCut(sel*weight) , "met_pow"    , 10 , 0 , 500 );

  TH1F* y_mg      = getHist( ttall    , "y"       , TCut(sel*weight) , "y_mg"       , 10 , 0 , 30 );
  TH1F* y_pow     = getHist( ttpowheg , "y"       , TCut(sel*weight) , "y_pow"      , 10 , 0 , 30 );

  TCanvas *c1 = new TCanvas("c1","",1200,1200);
  c1->Divide(2,2);

  c1->cd(1);
  plotHist( njets_pow , njets_mg , "powheg" , "madgraph" , "npfjets" );
  c1->cd(2);
  plotHist( ht_pow , ht_mg , "powheg" , "madgraph" , "H_{T} (GeV)" );
  c1->cd(3);
  plotHist( met_pow , met_mg , "powheg" , "madgraph" , "pfmet (GeV)" );
  c1->cd(4);
  plotHist( y_pow , y_mg , "powheg" , "madgraph" , "y (GeV^{1/2})" );

  if( printgif ) c1->Print("../plots/tt_comparison.gif");
}


//------------------------------------
// print yield table
//------------------------------------

void printYieldTable( char* path , bool latex = false ){

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  
  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  
  printYields( mc , mclabels , data , sel , weight , latex );
 
}


void njets( char* path , bool printplot = false ){

  deleteHistos();
  initialize(path);
  
  //char* var = "njetsuncor";
  //char* var = "npfjets";
  char* var = "htuncor";

  char* sample = "data";
  TChain *ch = data;

  TH1F* data_v1  = getHist( ch , var , TCut("dilmass>81&&dilmass<101&&(leptype==0||leptype==1)&&ndavtx<4")           , "data_v1"  , 10 , 0 , 500 );
  TH1F* data_v2  = getHist( ch , var , TCut("dilmass>81&&dilmass<101&&(leptype==0||leptype==1)&&ndavtx>3&&ndavtx<8") , "data_v2"  , 10 , 0 , 500 );
  TH1F* data_v3  = getHist( ch , var , TCut("dilmass>81&&dilmass<101&&(leptype==0||leptype==1)&&ndavtx>7")           , "data_v3"  , 10 , 0 , 500 );

  data_v2->SetLineColor(2);
  data_v2->SetMarkerColor(2);
  data_v2->SetMarkerStyle(23);

  data_v3->SetLineColor(4);
  data_v3->SetMarkerColor(4);
  data_v3->SetMarkerStyle(25);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  gPad->SetLogy();
  data_v1->GetXaxis()->SetTitle(var);
  data_v1->DrawNormalized();
  data_v2->DrawNormalized("same");
  data_v3->DrawNormalized("same");

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(data_v1,"nvtx #leq 3");
  leg->AddEntry(data_v2,"4 < nvtx #leq 7");
  leg->AddEntry(data_v3,"nvtx > 7");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  if( printplot ) c1->Print(Form("../plots/%s_%s.gif",sample,var));
}



//------------------------------------
// do ABCD method
//------------------------------------

void plotABCD( char* path , bool printplot = false , bool latex = false ){

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);


  initSymbols(latex);
  
  width1     = 17;
  width2     = 3;
  linelength = (width1+width2)*7+1;

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  //----------------------
  // 2010 signal region
  //----------------------
  
  float x1=125;
  if( !highpt ) x1 = 200;
  float x2=300;
  float x3=300;
  float x4=1500;
  
  float y1 = 4.5;
  float y2 = 8.5;
  float y3 = 8.5;
  float y4 = 30;
  
  // float x1=125;
  // float x2=300;
  // float x3=300;
  // float x4=1500;
  
  // float y1 = 4.5;
  // float y2 = 13;
  // float y3 = 13;
  // float y4 = 30;
  
  TCut A(Form("htpf>%.0f && htpf<%.0f && y>%.1f",x1,x2,y3));
  TCut B(Form("htpf>%.0f && htpf<%.0f && y>%.1f && y<%.1f",x1,x2,y1,y2));
  TCut C(Form("htpf>%.0f && y>%.1f && y<%.1f",x3,y1,y2));
  TCut D(Form("htpf>%.0f && y>%.1f",x3,y3));

  cout << endl;
  cout << "Region A : " << A.GetTitle() << endl;
  cout << "Region B : " << B.GetTitle() << endl;
  cout << "Region C : " << C.GetTitle() << endl;
  cout << "Region D : " << D.GetTitle() << endl;
  cout << endl;

  const unsigned int nmc = mc.size();
  TH2F* hmc[nmc];

  printLine(latex);
  printABCDHeader();
  printLine(latex);

  TChain* mctot = new TChain("t");

  for( unsigned int i = 0 ; i < mc.size() ; ++i ){
    hmc[i] = doABCD( mc[i] , mclabels[i] , sel , weight , A , B , C , D );
    mctot->Add(mc[i]);
  }

  printLine(latex);

  TH2F* hmctot = doABCD( mctot , "Total SM MC" , sel , weight , A , B , C , D );

  printLine(latex);

  TH2F* hdata = doABCD( data , "data" , sel , TCut("1") , A , B , C , D );

  printLine(latex);

  //-------------------------
  // draw TH2 data/MC
  //-------------------------

  TCanvas *abcd_can = new TCanvas();
  abcd_can->cd();

  hmctot->Draw("box");
  hmctot->SetLineColor(4);
  hmctot->SetMarkerColor(0);
  hmctot->SetFillColor(0);
  hmctot->GetYaxis()->SetTitle("y [ #sqrt{GeV} ]");
  hmctot->GetYaxis()->SetTitleOffset(1);
  hmctot->GetXaxis()->SetTitle("H_{T} [ GeV ]");
  hdata->SetMarkerColor(2);
  hdata->SetMarkerSize(0.8);
  hdata->Draw("same");


  //-------------------------
  // draw legend
  //-------------------------
  
  TLegend *legall = new TLegend(0.7,0.7,0.9,0.9);     
  TH1F* hdummymc = (TH1F*) hmctot->Clone();
  hdummymc->SetMarkerStyle(25);
  hdummymc->SetMarkerSize(1);
  hdummymc->SetMarkerColor(4);
  legall->AddEntry(hdummymc,"SM MC","p");
  legall->AddEntry(hdata,"Data","p");
  legall->SetFillColor(0);
  legall->SetBorderSize(1);
  legall->Draw();


  //--------------------------
  // draw ABCD regions
  //--------------------------

  TLatex *text=new TLatex();
  text->SetTextSize(0.05);
  text->SetTextColor(1);
    
  drawSquare(x1,y1,x2,y2);
  drawSquare(x3,y1,x4,y2);
  drawSquare(x1,y3,x2,y4);
  drawSquare(x3,y3,x4,y4);
  
  text->DrawLatex(x1+50,y1+1.5,"B");
  text->DrawLatex(x1+50,y3+10,"A");
  text->DrawLatex(x3+50,y1+1.5,"C");
  text->DrawLatex(x3+50,y3+10,"D");

  
  text->SetNDC();
  text->SetTextSize(0.037);
  text->DrawLatex(0.35,0.85,"CMS");
  text->DrawLatex(0.35,0.80,"153 pb^{-1} at #sqrt{s} = 7 TeV");
  text->DrawLatex(0.35,0.75,"Events with ee/#mu#mu/e#mu");

  if( printplot ){
    abcd_can->Print("../plots/abcd.png");
    abcd_can->Print("../plots/abcd.pdf");
  }

}


pair<float,float> getObservedOverPredicted( TChain *ch , TCut sel , TCut weight , 
					    TCut A , TCut B , TCut C , TCut D ){

  TH1F* hA = new TH1F("hA","",1,0,1); hA->Sumw2();
  TH1F* hB = new TH1F("hB","",1,0,1); hB->Sumw2();
  TH1F* hC = new TH1F("hC","",1,0,1); hC->Sumw2();
  TH1F* hD = new TH1F("hD","",1,0,1); hD->Sumw2();

  TCanvas *abcd_temp = new TCanvas();
  ch->Draw("0.5>>hA",(sel+A)*weight);
  ch->Draw("0.5>>hB",(sel+B)*weight);
  ch->Draw("0.5>>hC",(sel+C)*weight);
  ch->Draw("0.5>>hD",(sel+D)*weight);
  delete abcd_temp;

  float nA = hA->GetBinContent(1);
  float nB = hB->GetBinContent(1);
  float nC = hC->GetBinContent(1);
  float nD = hD->GetBinContent(1);

  float eA = hA->GetBinError(1);
  float eB = hB->GetBinError(1);
  float eC = hC->GetBinError(1);
  float eD = hD->GetBinError(1);

  delete hA;
  delete hB;
  delete hC;
  delete hD;

  cout << "A: " << Form("%.1f +- %.1f",nA,eA) << endl;
  cout << "B: " << Form("%.1f +- %.1f",nB,eB) << endl;
  cout << "C: " << Form("%.1f +- %.1f",nC,eC) << endl;
  cout << "D: " << Form("%.1f +- %.1f",nD,eD) << endl;


  if( nA > 0 && nB > 0 && nC > 0 && nD > 0 ){

    float op    = nD / ( ( nA * nC ) / nB );

    float operr = op * sqrt( pow(eA/nA,2) + pow(eB/nB,2) + pow(eC/nC,2) + pow(eD/nD,2) );

    return make_pair( op , operr );

  }

  else{
    return make_pair( -99 , -99 );
  }

  return make_pair(-9999,-9999);
}

void abcdClosure( char* path ){

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();
  
  float x1=125;
  float x2=300;
  float x3=300;
  float x4=1500;
  
  float y1 = 4.5;
  float y2 = 8.5;
  float y3 = 8.5;
  float y4 = 30;

  TCut A("");
  TCut B("");
  TCut C("");
  TCut D("");
    
  const unsigned int n = 10;

  float x[n];
  float xerr[n];
  float mg[n];
  float pow[n];
  float mgerr[n];
  float powerr[n];

  bool y_vary = false;

  for( int i = 0 ; i < n ; ++i ){

    if( y_vary ){
      float ycut = 6.5 + i;

      y2 = ycut;
      y3 = ycut;

      x[i]    = ycut;
      xerr[i] = 0.0;
    }

    else{
      float htcut = 200 + 50 * i;
      
      x2 = htcut;
      x3 = htcut;
      
      x[i]    = htcut;
      xerr[i] = 0.0;
    }


    A = TCut(Form("ht>%.0f && ht<%.0f && y>%.1f",x1,x2,y3));
    B = TCut(Form("ht>%.0f && ht<%.0f && y>%.1f && y<%.1f",x1,x2,y1,y2));
    C = TCut(Form("ht>%.0f && y>%.1f && y<%.1f",x3,y1,y2));
    D = TCut(Form("ht>%.0f && y>%.1f",x3,y3));

    //A = TCut(Form("ht>%.0f && ht<%.0f && pfmet/pow(htpf,0.45)>%.1f",x1,x2,y3));
    //B = TCut(Form("ht>%.0f && ht<%.0f && pfmet/pow(htpf,0.45)>%.1f && pfmet/pow(htpf,0.45)<%.1f",x1,x2,y1,y2));
    //C = TCut(Form("ht>%.0f && pfmet/pow(htpf,0.45)>%.1f && pfmet/pow(htpf,0.45)<%.1f",x3,y1,y2));
    //D = TCut(Form("ht>%.0f && pfmet/pow(htpf,0.45)>%.1f",x3,y3));

    cout << endl;
    cout << "Region A : " << A.GetTitle() << endl;
    cout << "Region B : " << B.GetTitle() << endl;
    cout << "Region C : " << C.GetTitle() << endl;
    cout << "Region D : " << D.GetTitle() << endl;
    cout << endl;

    pair<float,float> op_mg_pair  = getObservedOverPredicted( ttall , sel , weight , A , B , C , D );
    float opmg                    = op_mg_pair.first;
    float opmgerr                 = op_mg_pair.second;

    pair<float,float> op_pow_pair = getObservedOverPredicted( ttpowheg , sel , weight , A , B , C , D );
    float oppow                   = op_pow_pair.first;
    float oppowerr                = op_pow_pair.second;
    
    cout << "O/P (MG)  : " << opmg  << " +/- " << opmgerr  << endl;
    cout << "O/P (pow) : " << oppow << " +/- " << oppowerr << endl;

    mg[i]       = opmg;
    mgerr[i]    = opmgerr;
    pow[i]      = oppow;
    powerr[i]   = oppowerr;
    
  }

  TGraphErrors *gr_mg  = new TGraphErrors(n,x,mg ,xerr,mgerr);
  TGraphErrors *gr_pow = new TGraphErrors(n,x,pow,xerr,powerr);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  gPad->SetGridx(1);
  gPad->SetGridy(1);

  gr_mg->SetLineColor(2);
  gr_mg->SetMarkerColor(2);
  gr_mg->SetMarkerStyle(20);

  gr_pow->SetLineColor(4);
  gr_pow->SetMarkerColor(4);
  gr_pow->SetMarkerStyle(25);

  gr_mg->GetXaxis()->SetTitle("y cut (GeV^{1/2})");
  gr_mg->GetYaxis()->SetTitle("observed / predicted");
  gr_mg->Draw("AP");
  gr_pow->Draw("sameP");

  TLegend *leg = new TLegend(0.3,0.5,0.5,0.7);
  leg->AddEntry(gr_mg,"madgraph","p");
  leg->AddEntry(gr_pow,"powheg","p");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();


}

TH2F* doABCD( TChain *ch , char* label , TCut sel , TCut weight , TCut A , TCut B , TCut C , TCut D ){

  TCut dil("(nels+nmus+ntaus)==2");
  TCut ll ("(w1>0&&w2>0) && (w1<3&&w2<3)");
  TCut tau("(w1>0&&w2>0) && (w1>2||w2>2)");
  TCut fake("!(w1>0&&w2>0)");

  if     ( strcmp(label,"ttll")   == 0 ) sel = sel + ll;
  else if( strcmp(label,"tttau")  == 0 ) sel = sel + tau;
  else if( strcmp(label,"ttfake") == 0 ) sel = sel + fake;
  else if( strcmp(label,"ttdil")  == 0 ) sel = sel + dil;
  else if( strcmp(label,"ttotr")  == 0 ) sel = sel + !dil;

  TH1F* hA = new TH1F("hA","",1,0,1); hA->Sumw2();
  TH1F* hB = new TH1F("hB","",1,0,1); hB->Sumw2();
  TH1F* hC = new TH1F("hC","",1,0,1); hC->Sumw2();
  TH1F* hD = new TH1F("hD","",1,0,1); hD->Sumw2();
  TH2F* h  = new TH2F(Form("%s_abcd",label),Form("%s_abcd",label),60,0,1500,60,0,30); h->Sumw2();

  TCanvas *abcd_temp = new TCanvas();
  ch->Draw("0.5>>hA",(sel+A)*weight);
  ch->Draw("0.5>>hB",(sel+B)*weight);
  ch->Draw("0.5>>hC",(sel+C)*weight);
  ch->Draw("0.5>>hD",(sel+D)*weight);
  ch->Draw(Form("y:ht>>%s_abcd",label),sel*weight);
  delete abcd_temp;

  float nA = hA->GetBinContent(1);
  float nB = hB->GetBinContent(1);
  float nC = hC->GetBinContent(1);
  float nD = hD->GetBinContent(1);

  float eA = hA->GetBinError(1);
  float eB = hB->GetBinError(1);
  float eC = hC->GetBinError(1);
  float eD = hD->GetBinError(1);

  printABCDRow( label , nA , nB , nC , nD , eA , eB , eC , eD );

  delete hA;
  delete hB;
  delete hC;
  delete hD;

  return h;

}



void printABCDRow( char* sample, float A, float B, float C, float D, float dA, float dB, float dC, float dD ){
  
  float pred    = B > 0 ? A * C / B : 0.;
  float prederr = (A > 0 && B > 0 && D > 0 ) ? pred * sqrt( pow( dA/A , 2 ) + pow( dB/B , 2) + pow( dC/C , 2) ) : 0.;
  

  float ratio    = pred > 0 ? D / pred : 0;
  float ratioerr = (A>0 && B>0 && C>0 && D>0) ? ratio * sqrt( pow(dA/A,2) + pow(dB/B,2) + pow(dC/C,2) + pow(dD/D,2) ) : 0;

  stringstream sA;
  stringstream sB;
  stringstream sC;
  stringstream sD;
  stringstream sPred;
  stringstream sRatio;
    
  if( sample == "data"){
    sA     << A;
    sB     << B;
    sC     << C;
    sD     << D;
    sPred  << Form("%.1f %s %.1f",  pred , pm ,  prederr );
    sRatio << Form("%.2f %s %.2f", ratio , pm , ratioerr );
  }
  else{
    sA     << Form("%.1f %s %.1f",  A    , pm , dA       );
    sB     << Form("%.1f %s %.1f",  B    , pm , dB       );
    sC     << Form("%.1f %s %.1f",  C    , pm , dC       );
    sD     << Form("%.2f %s %.2f",  D    , pm , dD       );
    sPred  << Form("%.2f %s %.2f",  pred , pm , prederr  );
    sRatio << Form("%.2f %s %.2f", ratio , pm , ratioerr );
  }
  
  cout  << delimstart << setw(width1) << sample        << setw(width2)
        << delim      << setw(width1) << sA.str()      << setw(width2)
        << delim      << setw(width1) << sB.str()      << setw(width2)
        << delim      << setw(width1) << sC.str()      << setw(width2)
        << delim      << setw(width1) << sD.str()      << setw(width2)
        << delim      << setw(width1) << sPred.str()   << setw(width2) 
        << delim      << setw(width1) << sRatio.str()  << setw(width2)  << delimend << endl;
}


void printABCDHeader(){

  cout  << delimstart << setw(width1) << "sample"   << setw(width2) 
        << delim      << setw(width1) << "A"        << setw(width2) 
        << delim      << setw(width1) << "B"        << setw(width2) 
        << delim      << setw(width1) << "C"        << setw(width2) 
        << delim      << setw(width1) << "D"        << setw(width2) 
        << delim      << setw(width1) << "pred"     << setw(width2) 
        << delim      << setw(width1) << "obs/pred" << setw(width2) 
        << delimend << endl;

}

void drawSquare(float x1, float y1, float x2, float y2, int color){

  TLine *line = new TLine();
  line->SetLineColor(color);
  line->SetLineWidth(2);

  line->DrawLine(x1,y1,x2,y1);
  line->DrawLine(x2,y1,x2,y2);
  line->DrawLine(x2,y2,x1,y2);
  line->DrawLine(x1,y2,x1,y1);

  delete line;
  
}


//--------------------------------------------------
// make data/MC plots
//--------------------------------------------------

void leptonpt( char* path , bool printgif = false ){

  deleteHistos();
  
  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;
  
  initialize(path);


  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  TCanvas *lepcan = new TCanvas("lepcan","lepcan",1200,1200);
  lepcan->Divide(2,2);

  lepcan->cd(1);
  compareDataMC( mc , mclabels , data , "lep1.pt()" , TCut(sel+"abs(id1)==11") , weight , 5 , 10 , 110 , "leading electron p_{T} (GeV)"  );

  lepcan->cd(2);
  compareDataMC( mc , mclabels , data , "lep1.pt()" , TCut(sel+"abs(id1)==13") , weight , 5 , 10 , 110 , "leading muon p_{T} (GeV)"      );

  lepcan->cd(3);
  compareDataMC( mc , mclabels , data , "lep2.pt()" , TCut(sel+"abs(id2)==11") , weight , 4 , 0 , 20 , "trailing electron p_{T} (GeV)" );

  lepcan->cd(4);
  compareDataMC( mc , mclabels , data , "lep2.pt()" , TCut(sel+"abs(id2)==13") , weight , 4 , 0 , 20 , "trailing muon p_{T} (GeV)"     );

  if( printgif ) lepcan->Print("../plots/leppt.png");

}

//--------------------------------------------------
// make data/MC plots
//--------------------------------------------------

void makePlots( char* path , bool printgif = false ){

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  //sel = sel + "y > 8.5 && ht > 300";

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  if( highpt ){
    //vars.push_back("pfmet");     xt.push_back("pfmet (GeV)");      n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    //vars.push_back("y");         xt.push_back("y #equiv MET  /  #sqrt{H_{T}} (GeV^{1/2})");    n.push_back(10); xi.push_back(0.); xf.push_back(20.);
    //vars.push_back("htpf");      xt.push_back("H_{T} (GeV)");      n.push_back(10); xi.push_back(0.); xf.push_back(1000.);
    //vars.push_back("dilmass");   xt.push_back("M(ll) (GeV)");      n.push_back(60); xi.push_back(0.); xf.push_back(300.);
    //vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(10); xi.push_back(0.); xf.push_back(300.);
    //vars.push_back("npfjets");   xt.push_back("njets");            n.push_back(10); xi.push_back(0.); xf.push_back(10.);
    vars.push_back("ndavtx");      xt.push_back("nDAVertices");      n.push_back(20); xi.push_back(0.); xf.push_back(20.);
    //vars.push_back("htoffset");  xt.push_back("L1Offset-H_{T} (GeV)");    n.push_back(10); xi.push_back(0.); xf.push_back(1000.);
    //vars.push_back("htuncor");   xt.push_back("uncorrected H_{T} (GeV)");    n.push_back(10); xi.push_back(0.); xf.push_back(1000.);
  }else{
    vars.push_back("pfmet");     xt.push_back("pfmet (GeV)");      n.push_back(5); xi.push_back(0.); xf.push_back(250.);
    vars.push_back("htpf");      xt.push_back("H_{T} (GeV)");      n.push_back(7); xi.push_back(0.); xf.push_back(700.);
    vars.push_back("dilmass");   xt.push_back("M(ll) (GeV)");      n.push_back(5); xi.push_back(0.); xf.push_back(100.);
    vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(5); xi.push_back(0.); xf.push_back(100.);
  } 

  //vars.push_back("njets");     xt.push_back("njets");            n.push_back(6);  xi.push_back(0);  xf.push_back(6);
  //vars.push_back("dilep.mass()");     xt.push_back("dilmass (GeV)");    n.push_back(100);  xi.push_back(0);  xf.push_back(200);
  
  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  bool combine4 = false;
  int canCounter = -1;
  
  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    if( combine4 ){
      if( ivar % 4 == 0 ){
	canCounter++;
	can[canCounter] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),1400,1200);

	legpad[canCounter] = new TPad("legpad","legpad",12./14.,0,1,1);
	legpad[canCounter]->Draw();
	legpad[canCounter]->cd();

	TLegend *leg = getLegend( mc , mclabels , true , 0.2 , 0.3 , 0.8 , 0.7 );
	leg->SetTextSize(0.1);
	leg->SetBorderSize(1);
	leg->Draw();

	can[canCounter]->cd();

	plotpad[canCounter] = new TPad("plotpad","plotpad",0,0,12./14.,1);
	plotpad[canCounter]->Draw();
	plotpad[canCounter]->cd();

	plotpad[canCounter]->Divide(2,2);
	plotpad[canCounter]->cd(1);

      }else{
	plotpad[canCounter]->cd(1+ivar%4);
      }
    }else{
      can[ivar] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),600,600);
    }

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , !combine4 );

    if( printgif ) can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
  } 
}



void makeStandardPlots( char* path , bool sigregion = false , bool printgif = false ){

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  if( sigregion )  sel = sel + "y > 8.5 && htpf > 300";

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  if( highpt ){
    vars.push_back("lep1.pt()");  xt.push_back("max lepton p_{T} (GeV)");	n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("lep2.pt()");  xt.push_back("min lepton p_{T} (GeV)");	n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("lep1.eta()"); xt.push_back("max lepton #eta")       ;	n.push_back(10); xi.push_back(-5.); xf.push_back(5.);
    // vars.push_back("lep2.eta()"); xt.push_back("min lepton #eta");	        n.push_back(10); xi.push_back(-5.); xf.push_back(5.);
    // vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	        n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	                n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("dphijm");     xt.push_back("#phi(max jet,pfmet)");	        n.push_back(10); xi.push_back(0.); xf.push_back(3.2);
    // vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");			n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");			n.push_back(10); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("y");          xt.push_back("y (GeV^{1/2})");		n.push_back(10); xi.push_back(0.); xf.push_back(20.);
    // vars.push_back("htpf");       xt.push_back("H_{T} (GeV)");			n.push_back(10); xi.push_back(0.); xf.push_back(1000.);
    // vars.push_back("dilmass");    xt.push_back("M(ll) (GeV)");			n.push_back(10); xi.push_back(0.); xf.push_back(300.);
    // vars.push_back("dilpt");      xt.push_back("p_{T}(ll) (GeV)");		n.push_back(10); xi.push_back(0.); xf.push_back(300.);
    // vars.push_back("npfjets");    xt.push_back("npfjets");			n.push_back(10); xi.push_back(0.); xf.push_back(10.);
    // vars.push_back("nbtags");     xt.push_back("nbtags");			n.push_back(10); xi.push_back(0.); xf.push_back(10.);
    // vars.push_back("ndavtx");     xt.push_back("nDAVertices");			n.push_back(20); xi.push_back(0.); xf.push_back(20.);
    // vars.push_back("mt2jcore");   xt.push_back("MT2J");				n.push_back(10); xi.push_back(0.); xf.push_back(500.);
    // vars.push_back("mt2core");    xt.push_back("MT2");				n.push_back(10); xi.push_back(0.); xf.push_back(500.);
    // vars.push_back("meff");       xt.push_back("meff");				n.push_back(10); xi.push_back(0.); xf.push_back(500.);
    
  }else{

  } 
  
  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  int canCounter = -1;
  
  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    canCounter++;
    can[canCounter] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),1400,1200);

    legpad[canCounter] = new TPad("legpad","legpad",12./14.,0,1,1);
    legpad[canCounter]->Draw();
    legpad[canCounter]->cd();

    TLegend *leg = getLegend( mc , mclabels , true , 0.2 , 0.3 , 0.8 , 0.7 );
    leg->SetTextSize(0.1);
    leg->SetBorderSize(1);
    leg->Draw();

    can[canCounter]->cd();

    plotpad[canCounter] = new TPad("plotpad","plotpad",0,0,12./14.,1);
    plotpad[canCounter]->Draw();
    plotpad[canCounter]->cd();
    
    plotpad[canCounter]->Divide(2,2);

    plotpad[canCounter]->cd(1);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==0") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , false );
    plotpad[canCounter]->cd(2);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==1") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , false );
    plotpad[canCounter]->cd(3);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==2") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , false );
    plotpad[canCounter]->cd(4);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel)              , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , false );

    if( printgif ) can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
  } 
}


//--------------------------------------------------
// make data/MC plots
//--------------------------------------------------

void ofsubtraction( char* path , bool printgif = false ){

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut njets2("npfjets >= 2");
  TCut met100("pfmet > 100");
  TCut ht100("htpf > 100");
  TCut pt2010("lep1.pt()>20 && lep2.pt()>10");
  TCut eetype("leptype==0");
  TCut mmtype("leptype==1");
  TCut emtype("leptype==2");
  TCut sf=eetype||mmtype;
  TCut of=emtype;
  TCut weight("weight * ndavtxweight");

  TCut ofsel = njets2 + met100 + ht100 + pt2010;

  cout << "Selection: " << ofsel.GetTitle() << endl;

  TCanvas *can_of = new TCanvas("can_of","can_of",1200,600);
  can_of->Divide(2,1);
  
  can_of->cd(1);
  compareDataMC( mc , mclabels , data , "dilmass" , TCut(ofsel+sf) , weight , 10 , 0 , 200 , "M(ll) (GeV)" , true );

  can_of->cd(2);
  compareDataMC( mc , mclabels , data , "dilmass" , TCut(ofsel+of) , weight , 10 , 0 , 200 , "M(ll) (GeV)" , true );

  TH1F* data_sf  = getHist( data   , "dilmass" , TCut((ofsel+sf)) , "data_sf"  , 10 , 0 , 200 );
  TH1F* data_of  = getHist( data   , "dilmass" , TCut((ofsel+of)) , "data_of"  , 10 , 0 , 200 );

  TH1F* sm_sf    = getHist( sm     , "dilmass" , TCut((ofsel+sf)*weight) , "sm_sf"    , 10 , 0 , 200 );
  TH1F* sm_of    = getHist( sm     , "dilmass" , TCut((ofsel+of)*weight) , "sm_of"    , 10 , 0 , 200 );
  TH1F* smLM1_sf = getHist( sm_LM1 , "dilmass" , TCut((ofsel+sf)*weight) , "smLM1_sf" , 10 , 0 , 200 );
  TH1F* smLM1_of = getHist( sm_LM1 , "dilmass" , TCut((ofsel+of)*weight) , "smLM1_of" , 10 , 0 , 200 );

  TH1F* data_fs  = (TH1F*) data_sf ->Clone("data_fs");
  TH1F* sm_fs    = (TH1F*) sm_sf   ->Clone("sm_fs");
  TH1F* smLM1_fs = (TH1F*) smLM1_sf->Clone("smLM1_fs");

  data_fs   ->Add(data_of,-1);
  sm_fs   ->Add(sm_of,-1);
  smLM1_fs->Add(smLM1_of,-1);

  TCanvas *can_of2 = new TCanvas("can_of2","can_of2",600,600);
  can_of2->cd();
  gPad->SetGridy();

  sm_fs->SetLineColor(4);
  sm_fs->SetLineWidth(2);
  sm_fs->SetMarkerColor(4);
  //sm_fs->SetFillColor(4);
  sm_fs->SetMarkerStyle(20);
  smLM1_fs->SetLineColor(2);
  smLM1_fs->SetLineWidth(2);
  smLM1_fs->SetMarkerColor(2);
  //smLM1_fs->SetFillColor(2);
  smLM1_fs->SetMarkerStyle(25);
  smLM1_fs->GetXaxis()->SetTitle("M(ll) (GeV)");
  smLM1_fs->GetYaxis()->SetTitle("SF - OF Yield");

  smLM1_fs->SetMinimum(-12);
  smLM1_fs->SetMaximum(12);
  smLM1_fs->Draw("hist");
  sm_fs->Draw("samehist");
  data_fs->SetMarkerSize(1.5);
  data_fs->SetLineWidth(2);
  data_fs->Draw("same");

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(data_fs,"data");
  leg->AddEntry(sm_fs,"SM MC","l");
  leg->AddEntry(smLM1_fs,"SM+LM1 MC","l");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
}


TH1F* getHist( TChain* ch , char* var , TCut sel , char* histname , int nbins , float xmin , float xmax ){

  TH1F* h = new TH1F(histname,histname,nbins,xmin,xmax);
  h->Sumw2();
  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("TMath::Min(%s,%f) >> %s",var,xmax-0.001,histname),sel);
  delete ctemp;
  return h;
}

TH2F* getHist2D( TChain* ch , char* varx , char* vary , TCut sel , char* histname , 
		 int nbinsx , float xmin , float xmax , int nbinsy , float ymin , float ymax ){
  
  TH2F* h = new TH2F(histname,histname,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  h->Sumw2();

  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("TMath::Min(%s,%f):TMath::Min(%s,%f) >> %s",vary,ymax-0.001,varx,xmax-0.01,histname),sel);
  delete ctemp;
  return h;
}


float calculateHistError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}




