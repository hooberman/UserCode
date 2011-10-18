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
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <sstream>
#include "Hootilities.C"
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
    data2010->Reset();
    datasf->Reset();
    datadf->Reset();
    ttall->Reset();
    ttpowheg->Reset();
    ttdil->Reset();
    ttotr->Reset();
    ttll->Reset();
    tttau->Reset();
    ttfake->Reset();
    zjets->Reset();
    dy->Reset();
    dydata->Reset();
    dytautau->Reset();
    ww->Reset();
    wz->Reset();
    zz->Reset();
    vv->Reset();
    t->Reset();
    wjets->Reset();
    LM0->Reset();
    LM1->Reset();
    LM2->Reset();
    LM3->Reset();
    LM4->Reset();
    LM5->Reset();
    LM6->Reset();
    LM7->Reset();
    LM8->Reset();
    LM9->Reset();
    LM10->Reset();
    LM11->Reset();
    LM12->Reset();
    LM13->Reset();
    T2tt->Reset();
    sm->Reset();
    sm_LM1->Reset();
    other->Reset();

    mc.clear();
    mctex.clear();
    mclabels.clear();
  }

  else{

    data	= new TChain("t");
    data2010	= new TChain("t");
    datasf	= new TChain("t");
    datadf	= new TChain("t");
    ttall	= new TChain("t");
    ttpowheg	= new TChain("t");
    ttdil	= new TChain("t");
    ttotr	= new TChain("t");
    ttll	= new TChain("t");
    tttau	= new TChain("t");
    ttfake	= new TChain("t");
    zjets      	= new TChain("t");
    dy		= new TChain("t");
    dydata	= new TChain("t");
    dytautau	= new TChain("t");
    ww		= new TChain("t");
    wz		= new TChain("t");
    zz		= new TChain("t");
    vv		= new TChain("t");
    t		= new TChain("t");
    wjets	= new TChain("t");
    LM0 	= new TChain("t");
    LM1 	= new TChain("t");
    LM2 	= new TChain("t");
    LM3 	= new TChain("t");
    LM4 	= new TChain("t");
    LM5 	= new TChain("t");
    LM6 	= new TChain("t");
    LM7 	= new TChain("t");
    LM8 	= new TChain("t");
    LM9 	= new TChain("t");
    LM10 	= new TChain("t");
    LM11 	= new TChain("t");
    LM12 	= new TChain("t");
    LM13 	= new TChain("t");
    T2tt 	= new TChain("t");
    sm          = new TChain("t");
    sm_LM1      = new TChain("t");
    other       = new TChain("t");
  }

  cout << endl;
  cout << "Loading babies at       : " << path << endl;
  
  data->Add(Form("%s/data_smallTree.root",path));
  //data->Add(Form("%s/data_smallTree_hbhe.root",path));
  //data->Add(Form("%s/data_smallTree_allpr_nojson.root",path));
  dydata->Add(Form("%s/data_smallTree.root",path));
  datasf->Add(Form("%s/data_smallTree_singleFake.root",path));
  datadf->Add(Form("%s/data_smallTree_doubleFake.root",path));
  data2010->Add("/tas03/home/benhoob/OSSusy/output_38X/nov5th_v6_skim/dataskim_smallTree.root");
  ttall->Add(Form("%s/ttall_smallTree.root",path));
  //cout << "ADDING TTBAR 10% SAMPLE" << endl;
  //ttall->Add(Form("%s/ttall_smallTree_tenPercent.root",path));
  ttpowheg->Add(Form("%s/ttpowheg_smallTree.root",path));
  ttdil->Add(Form("%s/ttll_smallTree.root",path));
  ttdil->Add(Form("%s/tttau_smallTree.root",path));
  ttotr->Add(Form("%s/ttotr_smallTree.root",path));
  ttotr->Add(Form("%s/ttall_smallTree.root",path));
  ttll->Add(Form("%s/ttll_smallTree.root",path));
  tttau->Add(Form("%s/tttau_smallTree.root",path));
  ttfake->Add(Form("%s/ttfake_smallTree.root",path));
  zjets->Add(Form("%s/Zjets_smallTree.root",path));
  dy->Add(Form("%s/DYtot_smallTree.root",path));
  //dy->Add(Form("%s/DYtot_njets2_smallTree.root",path));
  dytautau->Add(Form("%s/DYtot_smallTree.root",path));
  ww->Add(Form("%s/ww_smallTree.root",path));
  wz->Add(Form("%s/wz_smallTree.root",path));
  zz->Add(Form("%s/zz_smallTree.root",path));
  vv->Add(Form("%s/ww_smallTree.root",path));
  vv->Add(Form("%s/wz_smallTree.root",path));
  vv->Add(Form("%s/zz_smallTree.root",path));
  t->Add(Form("%s/tW_smallTree.root",path));
  wjets->Add(Form("%s/wjetsMG_smallTree.root",path));
  //wjets->Add(Form("%s/wjets_smallTree.root",path));
  LM0->Add(Form("%s/LM0_smallTree.root",path));
  LM1->Add(Form("%s/LM1_smallTree.root",path));
  LM2->Add(Form("%s/LM2_smallTree.root",path));
  LM3->Add(Form("%s/LM3_smallTree.root",path));
  LM4->Add(Form("%s/LM4_smallTree.root",path));
  LM5->Add(Form("%s/LM5_smallTree.root",path));
  LM6->Add(Form("%s/LM6_smallTree.root",path));
  LM7->Add(Form("%s/LM7_smallTree.root",path));
  LM8->Add(Form("%s/LM8_smallTree.root",path));
  LM9->Add(Form("%s/LM9_smallTree.root",path));
  //LM10->Add(Form("%s/LM10_smallTree.root",path));
  LM11->Add(Form("%s/LM11_smallTree.root",path));
  LM12->Add(Form("%s/LM12_smallTree.root",path));
  LM13->Add(Form("%s/LM13_smallTree.root",path));
  T2tt->Add(Form("%s/T2tt_250_50_smallTree.root",path));
  other->Add(Form("%s/ww_smallTree.root",path));
  other->Add(Form("%s/wz_smallTree.root",path));
  other->Add(Form("%s/zz_smallTree.root",path));
  other->Add(Form("%s/tW_smallTree.root",path));

  //------------------------------
  // reverse order: for plotting
  //------------------------------
  

  //mc.push_back(t);        mclabels.push_back("t");           mctex.push_back("single top");
  //mc.push_back(zz);       mclabels.push_back("ZZ");          mctex.push_back("Z^0Z^0");
  //mc.push_back(wz);       mclabels.push_back("WZ");          mctex.push_back("W^{\\pm}Z^0");
  //mc.push_back(ww);       mclabels.push_back("WW");          mctex.push_back("W^+W^-");
  //mc.push_back(vv);       mclabels.push_back("VV");          mctex.push_back("VV");
  //mc.push_back(dy);       mclabels.push_back("DY");          mctex.push_back("DY");
  //mc.push_back(ttall);    mclabels.push_back("ttall");       mctex.push_back("tt");
  /*
  mc.push_back(ttll);     mclabels.push_back("ttll");    mctex.push_back("$t\\bar{b}\\rightarrow\\ell^+\\ell^-$");
  //mc.push_back(tttau);    mclabels.push_back("tttau");   mctex.push_back("$t\\bar{b}\\rightarrow\\ell^{\\pm}\\tau^{\\mp}$");
  mc.push_back(ttfake);   mclabels.push_back("ttfake");  mctex.push_back("$t\\bar{b}\\rightarrow$fake");
  mc.push_back(dy);       mclabels.push_back("DY");          mctex.push_back("DY");
  mc.push_back(ww);       mclabels.push_back("WW");          mctex.push_back("W^+W^-");
  mc.push_back(wz);       mclabels.push_back("WZ");          mctex.push_back("W^{\\pm}Z^0");
  mc.push_back(zz);       mclabels.push_back("ZZ");          mctex.push_back("Z^0Z^0");
  mc.push_back(t);        mclabels.push_back("t");           mctex.push_back("single top");
  mc.push_back(wjets);    mclabels.push_back("wjets");       mctex.push_back("$W^{\\pm}$+jets");

  mc.push_back(LM1);      mclabels.push_back("LM1");         mctex.push_back("LM1");
  mc.push_back(LM3);      mclabels.push_back("LM3");         mctex.push_back("LM3");
  mc.push_back(LM6);      mclabels.push_back("LM6");         mctex.push_back("LM6");
  */

  //mc.push_back(ttall);     mclabels.push_back("ttall");
  //mc.push_back(ttdil);     mclabels.push_back("ttdil");
  //mc.push_back(ttotr);     mclabels.push_back("ttotr");
  //mc.push_back(ttpowheg);  mclabels.push_back("ttpowheg");
  // //mc.push_back(ttotr);    mclabels.push_back("ttotr");  

  //mc.push_back(ttall);    mclabels.push_back("ttall");       mctex.push_back("tt");
  //mc.push_back(dy);       mclabels.push_back("DY");          mctex.push_back("DY");

  // mc.push_back(ttll);     mclabels.push_back("ttll");    mctex.push_back("$t\\bar{b}\\rightarrow\\ell^+\\ell^-$");
  // mc.push_back(tttau);    mclabels.push_back("tttau");   mctex.push_back("$t\\bar{b}\\rightarrow\\ell^{\\pm}\\tau^{\\mp}$");
  // mc.push_back(ttfake);   mclabels.push_back("ttfake");  mctex.push_back("$t\\bar{b}\\rightarrow$fake");

  //mc.push_back(dydata);   mclabels.push_back("DYdata");      mctex.push_back("DYdata"); 


  //mc.push_back(ttall);    mclabels.push_back("ttall");       mctex.push_back("tt");
  mc.push_back(ttll);     mclabels.push_back("ttll");    mctex.push_back("$t\\bar{b}\\rightarrow\\ell^+\\ell^-$");
  mc.push_back(tttau);    mclabels.push_back("tttau");   mctex.push_back("$t\\bar{b}\\rightarrow\\ell^{\\pm}\\tau^{\\mp}$");
  mc.push_back(ttfake);   mclabels.push_back("ttfake");  mctex.push_back("$t\\bar{b}\\rightarrow$fake");
  //mc.push_back(wjets);    mclabels.push_back("wjets");       mctex.push_back("$W^{\\pm}$+jets");
  //mc.push_back(datasf);   mclabels.push_back("single fakes");      mctex.push_back("single fakes");
  //mc.push_back(datadf);   mclabels.push_back("double fakes");      mctex.push_back("double fakes");
  //mc.push_back(zjets);    mclabels.push_back("zjets");     mctex.push_back("$Z^0$+jets");
  //mc.push_back(dydata);   mclabels.push_back("DYdata");      mctex.push_back("DYdata"); 
  mc.push_back(dy);       mclabels.push_back("DY");          mctex.push_back("DY");
  // mc.push_back(dytautau); mclabels.push_back("DYtautau");    mctex.push_back("DYtautau");
  //mc.push_back(ww);       mclabels.push_back("WW");          mctex.push_back("W^+W^-");
  //mc.push_back(wz);       mclabels.push_back("WZ");          mctex.push_back("W^{\\pm}Z^0");
  //mc.push_back(zz);       mclabels.push_back("ZZ");          mctex.push_back("Z^0Z^0");
  mc.push_back(vv);       mclabels.push_back("VV");          mctex.push_back("VV");
  mc.push_back(t);        mclabels.push_back("t");           mctex.push_back("single top");
  mc.push_back(wjets);    mclabels.push_back("wjets");       mctex.push_back("$W^{\\pm}$+jets");
  mc.push_back(T2tt);     mclabels.push_back("T2tt 250/50"); mctex.push_back("T2tt 250/50");


  //mc.push_back(other);    mclabels.push_back("other");   mctex.push_back("WW/WZ/ZZ/t");
  //mc.push_back(LM0);      mclabels.push_back("LM0");     mctex.push_back("LM0");
  //mc.push_back(LM1);      mclabels.push_back("LM1");     mctex.push_back("LM1");
  //mc.push_back(LM2);      mclabels.push_back("LM2");     mctex.push_back("LM2");
  //mc.push_back(LM3);      mclabels.push_back("LM3");     mctex.push_back("LM3");
  //mc.push_back(LM4);      mclabels.push_back("LM4");     mctex.push_back("LM4");
  //mc.push_back(LM5);      mclabels.push_back("LM5");     mctex.push_back("LM5");
  //mc.push_back(LM6);      mclabels.push_back("LM6");     mctex.push_back("LM6");
  //mc.push_back(LM7);      mclabels.push_back("LM7");     mctex.push_back("LM7");
  //mc.push_back(LM8);      mclabels.push_back("LM8");     mctex.push_back("LM8");
  //mc.push_back(LM9);      mclabels.push_back("LM9");     mctex.push_back("LM9");
  //mc.push_back(LM10);     mclabels.push_back("LM10");    mctex.push_back("LM10");
  //mc.push_back(LM11);     mclabels.push_back("LM11");    mctex.push_back("LM11");
  //mc.push_back(LM12);     mclabels.push_back("LM12");    mctex.push_back("LM12");
  //mc.push_back(LM13);     mclabels.push_back("LM13");    mctex.push_back("LM13");
  

  int nmc = mc.size();
  for( int imc = 0 ; imc < nmc ; ++imc ) sm->Add(mc[imc]);

  // sm->Add(Form("%s/ttall_smallTree.root",path));
  // sm->Add(Form("%s/wjetsMG_smallTree.root",path));
  // //sm->Add(Form("%s/ttdil_smallTree.root",path));
  // //sm->Add(Form("%s/ttotr_smallTree.root",path));
  // sm->Add(Form("%s/DYtot_smallTree.root",path));
  // sm->Add(Form("%s/ww_smallTree.root",path));
  // sm->Add(Form("%s/wz_smallTree.root",path));
  // sm->Add(Form("%s/zz_smallTree.root",path));
  // sm->Add(Form("%s/tW_smallTree.root",path));


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

  TCut njets1("npfjets >= 1");
  TCut njets2("npfjets >= 2");
  TCut njets3("npfjets >= 3");
  TCut njets4("npfjets >= 4");
  TCut zveto("passz == 0");
  TCut zpass("dilep.mass()>76 && dilep.mass()<106");
  TCut met100("pfmet > 100");
  TCut met50("pfmet > 50");
  TCut met30("pfmet > 30");
  TCut ht100("htpf > 100");
  TCut ht200("htpf > 200");
  TCut ht250("htpf > 250");
  TCut ht300("htpf > 300");
  TCut pt2020("lep1.pt()>20 && lep2.pt()>20");
  TCut pt2010("lep1.pt()>20 && lep2.pt()>10");
  TCut pt1010("lep1.pt()>10 && lep2.pt()>10");
  TCut pt2015("lep1.pt()>20 && lep2.pt()>15");
  TCut pt105 ("lep1.pt()>10 && lep2.pt()>5" );
  TCut pt55  ("lep1.pt()>5  && lep2.pt()>5" );
  TCut relnt015 ("isont1<0.15&&isont2<0.15");
  TCut isolep1_t10_015("(isont1*lep1.pt())/max(lep1.pt(),20)<0.15");
  TCut isolep2_t10_015("(isont2*lep2.pt())/max(lep2.pt(),20)<0.15");
  TCut iso_t10_015 = isolep1_t10_015 + isolep2_t10_015;
  TCut ll ("w1>0 && w2>0");
  TCut lep2tight("lep2.pt()>10 || isont2<0.15");
  TCut goodrun("json==1");
  TCut hbhe("hbhe == 1");
  TCut nbtags1("nbtags > 0");
  TCut nbtags2("nbtags > 1");
  TCut trigee("leptype==0");
  TCut trigmm("leptype==1 * 0.90");
  TCut trigem("leptype==2 * 0.95");
  TCut trig=trigee||trigmm||trigem;

  TCut highptsel = zveto + njets2 + met50 + ht100 + pt2010;
  TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + !pt2010;

  TCut highmet ("pfmet>275 && htpf>300");
  TCut highht  ("pfmet>200 && htpf>600");
  TCut sig1    ("pfmet>275 && htpf>300 && htpf<600");
  TCut sig2    ("pfmet>275 && htpf>600");
  TCut sig3    ("pfmet>200 && pfmet<275 && htpf>600");
  
  // highptsel += highmet;
  // highptsel += highht;
  // highptsel += sig1;
  // highptsel += sig2;
  // highptsel += sig3;



  //----------------------------------------
  // selection for mlljj discrepancy
  //----------------------------------------

  // highptsel += "pfmet>275 && ht>300";
  // TCut eetype("leptype==0");
  // TCut mmtype("leptype==1");
  // TCut emtype("leptype==2");
  // TCut mynjets2("njets>=2");
  // TCut leadjet120("jet.pt()>120");
  // TCut mll50("dilmass>50");
  // TCut myzveto("dilmass<70 || dilmass>100 || leptype==2");
  // TCut pt2520("lep1.pt()>25 && lep2.pt()>20");
  // TCut mynbtags1("nbtagstcl>=1");

  // TCut sel;
  // sel += pt2520;
  // sel += mmtype;
  // sel += mynjets2;
  // sel += leadjet120;
  // sel += mll50;
  // sel += myzveto;
  // //sel += mynbtags1;

  // highptsel = sel;

  //highptsel = zveto + njets2 + met50 + ht200 + pt2010;
  //lowptsel  = zveto + njets2 + met50 + ht200 + pt2010;
  //highptsel = met50 + njets2 + ht100 + zveto;

  //lowptsel = ht200 + zpass;

  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + pt2010;

  //lowptsel = lowptsel + "dilmass<76||dilmass>106";
  //highptsel = highptsel + "pfmet>250 && htpf>250";

  //TCut highptsel = "dilmass>76&&dilmass<106";
  //TCut sig = "y >  8.5 && htpf>300";
  //TCut sig = "y > 13.0 && htpf>300";
  //TCut sig = "y >  8.5 && htpf>300";
  //TCut sig = "pfmet>275 && htpf>300";
  //TCut sig = "pfmet>200 && htpf>600";
  
  //TCut sr2010  = "y>8.5 && htpf>300";
  //lowptsel  = lowptsel  + sig;


  //TCut highptsel  = zveto + njets2 + met50 + ht200 + pt2010;
  //TCut lowptsel   = zveto + njets2 + met50 + ht200 + pt2010;

  //TCut highptsel = zveto + njets2 + met50 + ht100 + pt2010 + "pfmet>250&&htpf>300";
  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt105 + lep2tight;
  //TCut highptsel = met50 + njets2 + ht100 + zveto;

  //TCut lowptsel  = zveto + njets2 + met50 + ht200 + pt1010 + !pt2010 + lep2tight;
  //TCut sigyht("y>14 && htpf>300");
  //TCut sigmetht("pfmet>290 && ht>300");
  //highptsel = highptsel + sigyht;
  //highptsel = highptsel + sigmetht;
  //TCut highptsel = zveto + njets2 + pt2020 +met30;

  //TCut highptsel = zpass + pt2020;
  //TCut highptsel = zveto + njets2 + met30 + pt2020 + "nbtags>0";
  //TCut highptsel = !zveto + njets2  + pt2020 + "leptype==1";
  //TCut highptsel = zveto + pt2020 + njets2 + met30;
  //TCut highptsel = zveto + njets2 + pt2020;
  //TCut highptsel = pt2020 + "dilmass>76 && dilmass<106";
  //TCut highptsel = "run <= 161312";
  //TCut highptsel = zveto + pt2020;
  //TCut highptsel = pt2020 + zpass;
  //TCut highptsel = zveto + njets2 + "pfmet>75" + ht100 + pt2010;
  //TCut highptsel = zveto + pt2015;
  //TCut highptsel = zveto + met50+ pt2010 + njets2 + ht100 + "nbtags>0";
  //TCut highptsel = pt2015;// + "dilmass<76 || dilmass>106";;
  //TCut highptsel = zveto + pt2010 + "pfmet>100";
  //TCut highptsel = zpass;
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


  //TCut weight("weight * ndavtxweight * (1000./191.)");
  //TCut weight("weight * trgeff * 0.976 * 1.13");
  TCut weight("weight * 0.976 * trgeff * ndavtxweight * 1.13");
  //TCut weight("weight * 2.05");
  //TCut weight("weight * trgeff * 0.976");
  //TCut weight("weight * ndavtxweight * trgeff * 0.349");
  //TCut weight("weight * ndavtxweight * trgeff * 0.715");
  //TCut weight("weight  * trgeff * 0.349");
  //TCut weight("weight * ndavtxweight * trgeff");
  //TCut weight("weight * ndavtxweight");
  //TCut weight("weight * ndavtxweight * (90./191.)");
  //TCut weight("weight * (90./191.)");
  //TCut weight("weight");
  //TCut weight("1");
  //weight = weight * trig;

  cout << "Using weight            : " << weight.GetTitle() << endl;
  return weight;
}


void HBHE( char* path ){

  deleteHistos();
  
  initialize(path);

  TCut sel("dilmass>76&&dilmass<106&&lep1.pt()>20&&lep2.pt()>20&&leptype<2");
  TCut hbhe("hbhe==1");

  TH1F* met_pass   = getHist( data  , "pfmet"  , TCut(sel+hbhe)   , "met_pass"   , 50 , 0 , 100 );
  TH1F* met_fail   = getHist( data  , "pfmet"  , TCut(sel+!hbhe)  , "met_fail"   , 50 , 0 , 100 );

  cout << "FAIL " << met_fail->GetEntries() << endl;
  cout << "PASS " << met_pass->GetEntries() << endl;

  met_pass->Scale(1./met_pass->Integral());
  met_fail->Scale(1./met_fail->Integral());

  TCanvas *can = new TCanvas();
  can->cd();

  met_pass->GetXaxis()->SetTitle("pfmet (GeV)");
  met_pass->GetYaxis()->SetTitle("A.U.");
  met_pass->Draw("hist");
  met_fail->SetMarkerColor(2);
  met_fail->SetLineColor(2);
  met_fail->Draw("sameE1");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(met_pass,"pass HBHE","l");
  leg->AddEntry(met_fail,"fail HBHE","p");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  
  
}

void scaleFactors( char* path ){

  deleteHistos();
  
  initialize(path);

  TCut eetype("leptype==0");
  TCut mmtype("leptype==1");
  TCut zmass("dilmass>76&&dilmass<106&&lep1.pt()>20&&lep2.pt()>10");
  //TCut weight("weight*ndavtxweight");
  //TCut eeweight("weight*ndavtxweight*1.0");
  //TCut mmweight("weight*ndavtxweight*0.9");
  //TCut weight("weight");
  //TCut eeweight("weight");
  //TCut mmweight("weight");
  TCut weight  ("weight*ndavtxweight*trgeff*0.976");
  TCut eeweight("weight*ndavtxweight*trgeff*0.976");
  TCut mmweight("weight*ndavtxweight*trgeff*0.976");

  TH1F* hee = new TH1F("hee","",1,0,1);
  TH1F* hmm = new TH1F("hmm","",1,0,1);
  hee->Sumw2();
  hmm->Sumw2();

  TCanvas *ctemp = new TCanvas();
  sm->Draw("0.5>>hee",(eetype+zmass)*eeweight);
  sm->Draw("0.5>>hmm",(mmtype+zmass)*mmweight);
  
  float neeMC = hee->Integral();
  float nmmMC = hmm->Integral();

  data->Draw("0.5>>hee",eetype+zmass);
  data->Draw("0.5>>hmm",mmtype+zmass);
  delete ctemp;

  float needata = hee->Integral();
  float nmmdata = hmm->Integral();
  
  cout << "ee data    : " << needata       << endl;
  cout << "ee MC      : " << neeMC         << endl;
  cout << "mm data    : " << nmmdata       << endl;
  cout << "mm MC      : " << nmmMC         << endl;
  cout << "ee data/MC : " << needata/neeMC << endl;
  cout << "mm data/MC : " << nmmdata/nmmMC << endl;

  TLatex *text = new TLatex();
  text->SetTextSize(0.04);
  //text->SetTextColor(4);
  text->SetNDC();

  TCanvas *zcan = new TCanvas("zcan","zcan",1200,600);
  zcan->Divide(2,1);

  zcan->cd(1);
  compareDataMC( mc , mclabels , data , "dilmass" , eetype , weight , 100 , 0 , 200 , "M(e^{+}e^{-}) (GeV)"     , true, false, true, true, "ee" );
  //text->DrawLatex(0.2,0.75,Form("N_{DATA} = %.0f",needata));
  //text->DrawLatex(0.2,0.70,Form("N_{MC}   = %.0f",neeMC));
  //text->DrawLatex(0.2,0.65,Form("N_{DATA}/N_{MC} = %.2f",needata/neeMC));
  

  zcan->cd(2);
  compareDataMC( mc , mclabels , data , "dilmass" , mmtype , weight , 100 , 0 , 200 , "M(#mu^{+}#mu^{-}) (GeV)" , true, false, true, true, "mm" );
  //text->DrawLatex(0.2,0.75,Form("N_{DATA} = %.0f",nmmdata));
  //text->DrawLatex(0.2,0.70,Form("N_{MC}   = %.0f",nmmMC));
  //text->DrawLatex(0.2,0.65,Form("N_{DATA}/N_{MC} = %.2f",nmmdata/nmmMC));

}

void OF_LM( char* path , bool latex = false ){

  deleteHistos();
  
  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  TH1F* htemp = new TH1F("htemp","",1,0,1);
  htemp->Sumw2();

  TCut sig;
  TCut sftype("leptype==0||leptype==1");
  TCut dftype("leptype==2");

  //TChain *ch = LM1;
  //TChain *ch = LM3;
  TChain *ch = LM6;

  for( int i = 0 ; i < 2 ; ++i ){

    if( i == 0 ) sig = TCut("pfmet>275 && htpf>300");
    if( i == 1 ) sig = TCut("pfmet>200 && htpf>600");


    ch->Draw("0.5>>htemp",(sel+sig+sftype)*weight);
    float sf    = htemp->GetBinContent(1);
    float sferr = htemp->GetBinError(1);
    
    ch->Draw("0.5>>htemp",(sel+sig+dftype)*weight);
    float df    = htemp->GetBinContent(1);
    float dferr = htemp->GetBinError(1);
    
    float diff    = sf - df;
    float differr = sqrt(sferr*sferr+dferr*dferr);
    
    cout << endl;
    cout << "sig  " << sig.GetTitle() << endl;
    cout << "SF   " << Form("%.2f%s%.2f",sf,pm,sferr) << endl;
    cout << "DF   " << Form("%.2f%s%.2f",df,pm,dferr) << endl;
    cout << "diff " << Form("%.2f%s%.2f",diff,pm,differr) << endl;
    cout << endl;
  }

}



void plotMetHt( char* path , bool print = false ){

  float htmax = 2000.;

  deleteHistos();
  
  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  cout << "Trig eff correction     : 100% (ee), 90% (mm), 95% (em)" << endl;
  
  TCut sel    = selection_TCut(highpt);

  //TH2F* h = new TH2F("h","",10000,0,htmax,10000,0,500);
  TH2F* h  = new TH2F("h","",100,0,htmax,100,0,500);
  TH2F* h2 = new TH2F("h2","",100,0,htmax,100,0,30);

  TCanvas *ctemp = new TCanvas();
  data->Draw("min(pfmet,499.99):min(htpf,1499.99)>>h",sel);
  //ttall->Draw("min(pfmet,499.99):min(htpf,1499.99)>>h",sel);
  //ttall->Draw("min(y,29.99):min(htpf,1499.99)>>h2",sel);
  delete ctemp;

  TCanvas *c1 = new TCanvas();
  c1->cd();

  h->GetXaxis()->SetTitle("H_{T} [GeV]");
  h->GetYaxis()->SetTitle("E_{T}^{miss} [GeV]");
  h->GetXaxis()->SetNdivisions(5);

  h->Draw();

  drawSquare(600,200,htmax,500,2);
  drawSquare(300,275,htmax,500,4);

  TBox box;
  box.SetLineColor(2);
  box.SetLineWidth(3);
  box.SetFillColor(2);
  box.SetFillStyle(3004);
  box.DrawBox(600,200,htmax,500);
  box.SetLineColor(4);
  box.SetLineWidth(3);
  box.SetFillColor(4);
  box.SetFillStyle(3003);
  box.DrawBox(300,275,htmax,500);



  TLatex text;
  text.SetTextSize(0.03);
  text.SetNDC();
  
  text.DrawLatex(0.6,0.27,"CMS Preliminary");
  text.DrawLatex(0.6,0.23,"#sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 0.98 fb^{-1}");
  text.DrawLatex(0.6,0.19,"Events with ee/#mu#mu/e#mu");

  /*
  TCanvas *c2 = new TCanvas();
  c2->cd();
  h2->Draw();

  h2->GetXaxis()->SetTitle("H_{T} [GeV]");
  h2->GetYaxis()->SetTitle("y [GeV^{1/2}]");

  TF1* fmet1 = new TF1("fmet1","[0]/sqrt(x)",300,htmax);
  fmet1->SetParameter(0,275);
  fmet1->SetLineColor(4);
  fmet1->SetLineWidth(4);
  fmet1->Draw("same");

  TLine line1;
  line1.SetLineColor(4);
  line1.SetLineWidth(4);
  line1.DrawLine(300,275/sqrt(300),300,30);    
  line1.DrawLine(300,30,htmax,30);    
  line1.DrawLine(htmax,275/sqrt(htmax),htmax,30);    

  TF1* fmet2 = new TF1("fmet2","[0]/sqrt(x)",600,htmax);
  fmet2->SetParameter(0,200);
  fmet2->SetLineColor(2);
  fmet2->SetLineWidth(4);
  fmet2->Draw("same");

  TLine line2;
  line2.SetLineColor(2);
  line2.SetLineWidth(4);
  line2.DrawLine(600,200/sqrt(600),600,30);    
  line2.DrawLine(600,30,htmax,30);    
  line2.DrawLine(htmax,200/sqrt(htmax),htmax,30);    
  */


  //if( print ) c1->Print("../plots/met_ht_349pb.pdf");
  if( print ){
    c1->Print("../plots/met_ht_098fb.pdf");
    c1->Print("../plots/met_ht_098fb.eps");
    c1->Print("../plots/met_ht_098fb.C");
  }
}
void Routin( char* path ){

  deleteHistos();
  
  initialize(path);

  TH1F* htemp = new TH1F("htemp","htemp",1,0,1);
  htemp->Sumw2();

  TCut zpass("dilmass>76 && dilmass<106");
  TCut presel("npfjets>1 && htpf>100");
  TCut eetype("nels==2");
  TCut mmtype("nmus==2");
  TCut weight("weight*ndavtxweight");

  const unsigned int npoints = 13;

  float R_ee[npoints];
  float R_ee_err[npoints];
  float R_mm[npoints];
  float R_mm_err[npoints];
  float met[npoints];
  float meterr[npoints];

  //------------------------------------
  // ee
  //------------------------------------

  for( unsigned int i = 0 ; i < npoints ; ++i ){

    met[i]    = i*5;
    meterr[i] = 0.;

    TCut metcut(Form("pfmet>%i && pfmet<%i",i*5,(i+1)*5));

    dy->Draw("0.5>>htemp",(presel+eetype+zpass+metcut)*weight);
    float ee_in      = htemp->GetBinContent(1);
    float ee_in_err  = htemp->GetBinError(1);
    
    dy->Draw("0.5>>htemp",(presel+eetype+!zpass+metcut)*weight);
    float ee_out     = htemp->GetBinContent(1);
    float ee_out_err = htemp->GetBinError(1);
    
    float ree     = ee_out / ee_in;
    float ree_err = ree * sqrt( pow( ee_in_err/ee_in , 2 ) + pow( ee_out_err/ee_out , 2 ) );
    
    cout << "ee " << metcut.GetTitle() << endl;
    cout << "in      " << ee_in  << " +/- " << ee_in_err  << endl;
    cout << "out     " << ee_out << " +/- " << ee_out_err << endl;
    cout << "Rout/in " << ree    << " +/- " << ree_err    << endl;
    
    R_ee[i]     = ree;
    R_ee_err[i] = ree_err;
  }

  //------------------------------------
  // mm
  //------------------------------------

  for( unsigned int i = 0 ; i < npoints ; ++i ){

    TCut metcut(Form("pfmet>%i && pfmet<%i",i*5,(i+1)*5));

    dy->Draw("0.5>>htemp",(presel+mmtype+zpass+metcut)*weight);
    float mm_in      = htemp->GetBinContent(1);
    float mm_in_err  = htemp->GetBinError(1);
    
    dy->Draw("0.5>>htemp",(presel+mmtype+!zpass+metcut)*weight);
    float mm_out     = htemp->GetBinContent(1);
    float mm_out_err = htemp->GetBinError(1);
    
    float rmm     = mm_out / mm_in;
    float rmm_err = rmm * sqrt( pow( mm_in_err/mm_in , 2 ) + pow( mm_out_err/mm_out , 2 ) );
    
    cout << "mm " << metcut.GetTitle() << endl;
    cout << "in      " << mm_in  << " +/- " << mm_in_err  << endl;
    cout << "out     " << mm_out << " +/- " << mm_out_err << endl;
    cout << "Rout/in " << rmm    << " +/- " << rmm_err    << endl;
    
    R_mm[i]     = rmm;
    R_mm_err[i] = rmm_err;
  }


  TCanvas *rcan = new TCanvas();
  rcan->cd();
  gPad->SetGridy();

  TGraphErrors *gree = new TGraphErrors(npoints,met,R_ee,meterr,R_ee_err);
  TGraphErrors *grmm = new TGraphErrors(npoints,met,R_mm,meterr,R_mm_err);

  gree->SetMarkerColor(2);
  gree->SetLineColor(2);
  gree->SetMarkerStyle(21);

  grmm->SetMarkerColor(4);
  grmm->SetLineColor(4);
  grmm->SetMarkerStyle(25);

  gree->Draw("AP");
  grmm->Draw("sameP");

  gree->GetXaxis()->SetTitle("MET interval (GeV)");
  gree->GetYaxis()->SetTitle("R_{out/in}");

  TLegend *leg = new TLegend(0.25,0.65,0.45,0.85);
  leg->AddEntry(gree,"ee","p");
  leg->AddEntry(grmm,"#mu#mu","p");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  TLine line;
  line.SetLineColor(8);
  line.SetLineWidth(2);
  line.DrawLine(0,0.13,60,0.13);
  line.SetLineStyle(2);
  line.DrawLine(0,0.20,60,0.20);
  line.DrawLine(0,0.06,60,0.06);
  

  // TH1F* hee = new TH1F("hee","hee",150,0,300);
  // TH1F* hmm = new TH1F("hmm","hmm",150,0,300);
  // hee->Sumw2();
  // hmm->Sumw2();

  // TCanvas *can = new TCanvas("can","",1200,600);
  // can->Divide(2,1);
  // can->cd(1);
  // dy->Draw("dilmass>>hee",(presel+ee)*weight);
  // hee->Draw();
  // hee->GetXaxis()->SetTitle("M(e^{+}e^{-}) (GeV)");
  // can->cd(1);
  // dy->Draw("dilmass>>hmm",(presel+mm)*weight);
  // hmm->Draw();
  // hmm->GetXaxis()->SetTitle("M(#mu^{+}#mu^{-}) (GeV)");


  /*
  float ndata_ee_in = data->GetEntries("pass&&passz&&leptype==0");
  float ndata_mm_in = data->GetEntries("pass&&passz&&leptype==1");
  float ndata_em_in = data->GetEntries("pass&&passz&&leptype==2");

  float ndata_ee_tot = data->GetEntries("passz&&leptype==0");
  float ndata_mm_tot = data->GetEntries("passz&&leptype==1");

  
  float k = sqrt( ndata_mm_tot / ndata_ee_tot );

  cout << "k       " << k << endl;
  cout << "nee     " << ndata_ee_in << endl;
  cout << "nmm     " << ndata_mm_in << endl;
  cout << "nem     " << ndata_em_in << endl;


  float neepred    = ree * ( ndata_ee_in - (1/k) * ndata_em_in );
  float neeprederr = ree_err * ( ndata_ee_in - (1/k) * ndata_em_in );
  float nmmpred    = rmm * ( ndata_mm_in - k     * ndata_em_in );
  float nmmprederr = rmm_err * ( ndata_mm_in - k     * ndata_em_in );

  cout << "nee pred " << neepred << " +/- " << neeprederr << endl;
  cout << "nmm pred " << nmmpred << " +/- " << nmmprederr << endl;
  */

}



void jetzbalance( char* path ){

  deleteHistos();
  
  initialize(path);

  TCut sel("dilmass>81 && dilmass<101 && npfjets==1 && njets15==1 && (leptype==0||leptype==1) && dilpt>50");
  TCut weight("weight*ndavtxweight");

  TH1F* data_F23   = getHist( data  , "dilpt-ptjetF23"  , TCut(sel)        , "data_F23"   , 25 , -50 , 50 );
  TH1F* zjets_F23  = getHist( zjets , "dilpt-ptjetF23"  , TCut(sel*weight) , "zjets_F23"  , 25 , -50 , 50 );

  TH1F* data_O23   = getHist( data  , "dilpt-ptjetO23"  , TCut(sel)        , "data_O23"   , 25 , -50 , 50 );
  TH1F* zjets_O23  = getHist( zjets , "dilpt-ptjetO23"  , TCut(sel*weight) , "zjets_O23"  , 25 , -50 , 50 );

  TH1F* data_23    = getHist( data  , "dilpt-ptjet23"   , TCut(sel)        , "data_23"    , 25 , -50 , 50 );
  TH1F* zjets_23   = getHist( zjets , "dilpt-ptjet23"   , TCut(sel*weight) , "zjets_23"   , 25 , -50 , 50 );

  TH1F* data_raw   = getHist( data  , "dilpt-ptjetraw"  , TCut(sel)        , "data_raw"   , 25 , -50 , 50 );
  TH1F* zjets_raw  = getHist( zjets , "dilpt-ptjetraw"  , TCut(sel*weight) , "zjets_raw"  , 25 , -50 , 50 );

  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->cd();
  plotHist( zjets_F23 , data_F23 , "Z+jets MC" , "data" , "p_{T}(Z)-p_{T}(L1FastL2L3-jet)" );

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd();
  plotHist( zjets_O23 , data_O23 , "Z+jets MC" , "data" , "p_{T}(Z)-p_{T}(L1OffsetL2L3-jet)" );
  
  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->cd();
  plotHist( zjets_23  , data_23  , "Z+jets MC" , "data" , "p_{T}(Z)-p_{T}(L2L3-jet)" );
  
  TCanvas *c4 = new TCanvas("c4","",600,600);
  c4->cd();
  plotHist( zjets_raw , data_raw , "Z+jets MC" , "data" , "p_{T}(Z)-p_{T}(raw-jet)" );
  



}

void OOT(char* path ){

  deleteHistos();
  
  initialize(path);

  TCut sel("lep1.pt()>20 && lep2.pt()>20 && leptype==2");
  TCut v1("run<=161312");
  TCut v2("run>161312");
  
  TH1F* pfmet_v1  = getHist( data , "pfmet" , sel+v1 , "pfmet_v1"  , 20 , 0 , 500 );
  TH1F* pfmet_v2  = getHist( data , "pfmet" , sel+v2 , "pfmet_v2"  , 20 , 0 , 500 );

  pfmet_v1->Scale(1./pfmet_v1->Integral());
  pfmet_v2->Scale(1./pfmet_v2->Integral());

  TCanvas *c1 = new TCanvas();
  c1->cd();
  plotHist( pfmet_v1 , pfmet_v2 , "run #leq 161312" , "run > 161312" , "pfmet (GeV)" , true);
  
}


void compare2010(char* path ){

  deleteHistos();
  
  initialize(path);

  TCut sel2011("lep1.pt()>20 && lep2.pt()>20 && leptype==2");
  TCut sel2010("ptl1>20 && ptl2>20 && leptype==2");

  float n2010 = data2010->GetEntries(sel2010);
  float n2011 = data->GetEntries(sel2011);

  float scale = (191./34.) * 0.95;
  
  float ratio = n2011/(n2010*scale);
  float error = ratio * sqrt( 1./n2011 + 1./n2010 );
  
  cout << "n2010          " << n2010 << endl;
  cout << "n2010 (scaled) " << n2010 * scale << endl;
  cout << "n2011          " << n2011 << endl;
  cout << "2011/2010      " << Form("%.2f +/- %.2f",ratio,error) << endl;

  TH1F* pfmet_2010  = getHist( data2010 , "pfmet" , sel2010 , "pfmet_2010"  , 10 , 0 , 300 );
  TH1F* pfmet_2011  = getHist( data     , "pfmet" , sel2011 , "pfmet_2011"  , 10 , 0 , 300 );

  pfmet_2010->Scale(1./pfmet_2010->Integral());
  pfmet_2011->Scale(1./pfmet_2011->Integral());

  TCanvas *c1 = new TCanvas();
  c1->cd();
  plotHist( pfmet_2010 , pfmet_2011 , "2010 data" , "2011 data" , "pfmet (GeV)" , true);

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

void plotHist( TH1F* h1 , TH1F* h2 , char* leg1 , char* leg2 , char* xtitle , bool residual){


  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();

  if( residual ){
    fullpad = new TPad("fullpad","fullpad",0,0,1,1);
    fullpad->Draw();
    fullpad->cd();

    plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
    plotpad->Draw();
    plotpad->cd();
  }

  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetLogy(1);
  
  //--------------------------------
  // TF1* f1 = new TF1("f1","gaus",-20,20);
  // TF1* f2 = new TF1("f2","gaus",-20,20);

  // h1->Fit(f1,"R");
  // h2->Fit(f2,"R");

  //--------------------------------
  //h1->SetMarkerColor(4);
  //h1->SetLineColor(4);
  //h1->SetFillColor(4);
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  //h2->SetFillColor(2);
  h2->SetMarkerStyle(21);
  h1->GetXaxis()->SetTitle( xtitle );
  //h1->DrawNormalized("hist");
  h1->DrawNormalized("E1");
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->DrawNormalized("samehist");
  h2->DrawNormalized("sameE1");

  h1->SetMinimum(0);
  h1->SetMaximum(0.6);

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.DrawLine(30,0,30,0.6);

  TLegend *leg = new TLegend(0.75,0.75,0.95,0.95);
  leg->AddEntry(h1,leg1,"p");
  leg->AddEntry(h2,leg2,"l");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  // TLatex *text = new TLatex();
  // text->SetNDC();
  // text->SetTextSize(0.04);
  // text->DrawLatex(0.2,0.8,Form("<x> = %.1f GeV",h1->GetMean(1)));
  // text->SetTextColor(2);
  // text->DrawLatex(0.2,0.7,Form("<x> = %.1f GeV",h2->GetMean(1)));

  int bin50 = h1->FindBin(30);

  float eff1 = h1->Integral(bin50,1000)/h1->Integral();
  float eff2 = h2->Integral(bin50,1000)/h2->Integral();

  // float num1 =  h1->Integral(bin50,100);
  // float num2 =  h2->Integral(bin50,100);
  // float den1 =  h1->Integral();
  // float den2 =  h2->Integral();


  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  text->DrawLatex(0.5,0.6,Form("eff(MET>30 GeV) = %.2f",eff1));	      
  text->SetTextColor(2);
  text->DrawLatex(0.5,0.5,Form("eff(MET>30 GeV) = %.2f",eff2));



  if( residual ){
    fullpad->cd();

    respad = new TPad("respad","respad",0,0.8,1,1);
    respad->Draw();
    respad->cd();

    gPad->SetGridy();

    TH1F* ratio = (TH1F*) h2->Clone(Form("%s_ratio",h2->GetName()));
    ratio->Divide(h1);

    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetLabelSize(0.2);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetTitle("data/MC  ");
    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetXaxis()->SetTitleSize(0);
    ratio->SetMarkerSize(0.7);
    ratio->Draw();

    TLine myline;
    myline.SetLineWidth(1);
    myline.DrawLine(h1->GetXaxis()->GetXmin(),1,h1->GetXaxis()->GetXmax(),1);

  }

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

  if( printgif ) c1->Print("../plots/sigregion.gif");
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


void signalTable( char* path , bool latex = false ){

  deleteHistos();
  initialize(path);
  initSymbols(latex);

  TCut njets2("npfjets >= 2");
  TCut zveto("passz == 0");
  TCut met50("pfmet > 50");
  TCut ht100("htpf > 100");
  TCut pt2010("lep1.pt()>20 && lep2.pt()>10");
  TCut sig("y>8.5 && htpf > 300");
  //TCut sig("pfmet>250. && htpf > 300");

  TCut sigsel = zveto + njets2 + met50 + ht100 + pt2010 + sig;

  cout << "Using selection: " << sigsel.GetTitle() << endl;
  char* colformat = "precision=3 col=7.6::10.10:";

  data->Scan("run:lumi:event:lep1.pt():lep2.pt():leptype:pfmet:npfjets:htpf:dilmass:y:nbtags",sigsel,colformat);
}


void doOFCounting( char* path , bool latex = false ){

  gROOT->Reset();
  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  //------------------------------------------

  //TChain* ch = ttall;
  //TChain* ch = ttpowheg;
  //float R   = 1.14;
  //float Rup = 1.19;
  //float Rdn = 1.09;

  TChain* ch = data;
  float R   = 1.13;
  float Rup = 1.18;
  float Rdn = 1.08;
  weight = TCut("");

  sel = sel + "dilmass<76 || dilmass>106";
  //sel = sel + "y > 8.5 && htpf > 300.";
  //sel = sel + "y > 13 && htpf > 300.";
  //sel = sel + "y > 8.5 && htpf > 600.";
  sel = sel + "pfmet > 200 && htpf > 600.";
  //sel = sel + "pfmet > 275 && htpf > 300.";

  //------------------------------------------

  TH1F* htemp = new TH1F("htemp","",1,0,1);
  htemp->Sumw2();

  TCut eetype("leptype==0");
  TCut mmtype("leptype==1");
  TCut emtype("leptype==2");

  ch->Draw("0.5>>htemp",(sel+eetype)*weight);
  float nee    = htemp->GetBinContent(1);
  float neeerr = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+mmtype)*weight);
  float nmm    = htemp->GetBinContent(1);
  float nmmerr = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+emtype)*weight);
  float nem    = htemp->GetBinContent(1);
  float nemerr = htemp->GetBinError(1);

  cout << Form("ee  : %.2f %s %.2f",nee,pm,neeerr) << endl;
  cout << Form("mm  : %.2f %s %.2f",nmm,pm,nmmerr) << endl;
  cout << Form("em  : %.2f %s %.2f",nem,pm,nemerr) << endl;

  float delta     = R * nee + (1./R) * nmm - nem;
  float deltastat = sqrt( pow(R*neeerr,2) + pow((1./R)*nmmerr,2)+ pow(nemerr,2));
  float deltaup   = Rup * nee + (1./Rup) * nmm - nem;
  float deltadn   = Rdn * nee + (1./Rdn) * nmm - nem;
  float d1        = fabs( deltaup - delta );
  float d2        = fabs( deltadn - delta );
  float deltasyst = 0.5 * ( d1 + d2 );

  cout << Form("%.0f & %.0f & %.0f &  %.1f%s%.1f (stat)%s%.1f (syst) \\\\",nee,nmm,nem,delta,pm,deltastat,pm,deltasyst) << endl;


}


//------------------------------------
// print yield table
//------------------------------------

void jetMet( char* path , bool latex = false ){

  gROOT->Reset();
  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  //TCut nom = zveto + njets2 + met50 + ht100 + pt2010;
  

 
}


//------------------------------------
// print yield table
//------------------------------------

void printYieldTable( char* path , bool latex = false ){

  gROOT->Reset();
  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);

  cout << "Trig eff correction     : 100% (ee), 90% (mm), 95% (em)" << endl;
  
  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  printYields( mc , mclabels , data , sel , weight , latex );
 
}


void getK( char* path , bool latex = false ){

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  TCut jet_control("htpf>125&&htpf<300");
  TCut jet_ht300  ("htpf>300");
  TCut jet_ht600  ("htpf>600");
  
  TCut lt50("dilpt<50");
  TCut gt50("dilpt>50");

  //----------------------------
  // control region
  //----------------------------

  float ndata_control_lo = data->GetEntries(sel+jet_control+lt50);
  float ndata_control_hi = data->GetEntries(sel+jet_control+gt50);
  float ndata_control    = data->GetEntries(sel+jet_control);

  float kdata_control     = 1 + ndata_control_lo / ndata_control_hi;
  float kdata_control_err = kdata_control * sqrt( 1./ndata_control_lo + 1./ndata_control_hi );

  cout << "k(data) control : " << ndata_control << " / " << ndata_control_hi << " "
       << Form("%.2f%s%.2f",kdata_control,pm,kdata_control_err) << endl;

  //----------------------------
  // HT > 300 GeV region
  //----------------------------

  float ndata_ht300_lo = data->GetEntries(sel+jet_ht300+lt50);
  float ndata_ht300_hi = data->GetEntries(sel+jet_ht300+gt50);
  float ndata_ht300    = data->GetEntries(sel+jet_ht300);

  float kdata_ht300     = 1 + ndata_ht300_lo / ndata_ht300_hi;
  float kdata_ht300_err = kdata_ht300 * sqrt( 1./ndata_ht300_lo + 1./ndata_ht300_hi );

  cout << "k(data) ht300 : " << ndata_ht300 << " / " << ndata_ht300_hi << " "
       << Form("%.2f%s%.2f",kdata_ht300,pm,kdata_ht300_err) << endl;

  //----------------------------
  // HT > 600 GeV region
  //----------------------------

  float ndata_ht600_lo = data->GetEntries(sel+jet_ht600+lt50);
  float ndata_ht600_hi = data->GetEntries(sel+jet_ht600+gt50);
  float ndata_ht600    = data->GetEntries(sel+jet_ht600);

  float kdata_ht600     = 1 + ndata_ht600_lo / ndata_ht600_hi;
  float kdata_ht600_err = kdata_ht600 * sqrt( 1./ndata_ht600_lo + 1./ndata_ht600_hi );

  cout << "k(data) ht600 : " << ndata_ht600 << " / " << ndata_ht600_hi << " "
       << Form("%.2f%s%.2f",kdata_ht600,pm,kdata_ht600_err) << endl;

}


//--------------------------------------------------------
// estimate DY contamination in pt(ll) control region
//--------------------------------------------------------

pair<float,float>  doDYestimate( TChain *ch , TCut sel , float metcutval ){

  TCut dilptcut(Form("dilpt>%.1f",metcutval));
  sel = sel + dilptcut;
  cout << "DY estimate: selection " << sel.GetTitle() << endl;

  //-------------------
  // set parameters
  //-------------------

  bool  verbose = true;
  float R       = 0.13;
  float R_err   = 0.07;
  float k       = 1.13;

  if( verbose ){
    cout << "R : " << R << " +/- " << R_err << endl;
    cout << "k : " << k << endl;
  }

  //-------------------
  // invert Z-veto
  //-------------------

  TString selstring(sel.GetTitle());
  selstring.ReplaceAll("passz == 0","dilmass>76&&dilmass<106");
  TCut newsel = TCut(selstring);

  if( verbose ){
    cout << "Pre  : " << sel.GetTitle() << endl;
    cout << "Post : " << newsel.GetTitle() << endl;
  }

  //-------------------
  // get data yields
  //-------------------

  float ndata_ee_in = ch->GetEntries(newsel+"leptype==0");
  float ndata_mm_in = ch->GetEntries(newsel+"leptype==1");
  float ndata_em_in = ch->GetEntries(newsel+"leptype==2");

  if( verbose ){
    cout << "nee     " << ndata_ee_in << endl;
    cout << "nmm     " << ndata_mm_in << endl;
    cout << "nem     " << ndata_em_in << endl;
  }

  //-------------------
  // do DY estimate
  //-------------------

  float neepred     = R     * ( ndata_ee_in - (0.5/k)  * ndata_em_in );
  float nmmpred     = R     * ( ndata_mm_in - k/2.     * ndata_em_in );

  float neeprederr  = 0.;
  float nmmprederr  = 0.;

  neeprederr  += pow( R_err * ( ndata_ee_in - (0.5/k) * ndata_em_in ) , 2 );
  nmmprederr  += pow( R_err * ( ndata_mm_in - k/2.    * ndata_em_in ) , 2 );

  neeprederr  += pow( R * sqrt( ndata_ee_in + pow(0.5/k,2)  * ndata_em_in ) , 2 );
  nmmprederr  += pow( R * sqrt( ndata_mm_in + pow(k/2.,2)   * ndata_em_in ) , 2 );

  neeprederr = sqrt(neeprederr);
  nmmprederr = sqrt(nmmprederr);

  float ntotpred    = neepred + nmmpred;
  float ntotprederr = sqrt( pow(neeprederr,2) + pow(nmmprederr,2) );

  if( verbose ){
    cout << "nee   pred " << Form("%.2f%s%.2f",neepred,pm,neeprederr) << endl;
    cout << "nmm   pred " << Form("%.2f%s%.2f",nmmpred,pm,nmmprederr) << endl;
    cout << "ntot  pred " << Form("%.2f%s%.2f",ntotpred,pm,ntotprederr) << endl;
  }

  return make_pair( ntotpred , ntotprederr );

}

//------------------------------------
// victory
//------------------------------------

void victory( char* path , bool latex = false , string sigregion = "highmet" , bool print = false ){

  //-----------------
  // config params
  //-----------------

  const int nbins  = 16;
  const float xmin = 0;
  const float xmax = 400;

  TCut jetcut;
  float metcutval;
  float kc;
  float kcerr;
  float kdata;
  float kdataerr;
  char* htrange;
  float ymax;
  float ymin;
  int   rebin = 1;

  //string sigregion = "control";
  //string sigregion = "2010";
  //string sigregion = "highmet";
  //string sigregion = "highht";

  if( sigregion == "control" ){
    jetcut    = TCut("htpf>125 && htpf<300");
    metcutval = 200.;
    
    kc        = 1.3;
    kcerr     = 0.2;

    kdata     = 1.65;
    kdataerr  = 0.08;

    htrange   = "125 < H_{T} < 300 GeV";
    ymin      = 0.01;
    ymax      = 50000;
  }
  
  else if( sigregion == "highmet" ){
    jetcut = TCut("htpf>300");  
    metcutval = 275.0;

    kc        = 1.5;
    kcerr     = 0.5;

    //kdata     = 1.48;
    //kdataerr  = 0.17;
    kdata     = 1.43;
    kdataerr  = 0.18;

    htrange   = "H_{T} > 300 GeV";
    ymin      = 0.01;
    ymax      = 5000;
  }

  else if( sigregion == "highht" ){
    jetcut    = TCut("htpf>600");  
    metcutval = 200.0;

    //kc        = 1.3;
    kc        = 1.4;
    kcerr     = 0.4;

    //kdata     = 1.32;
    //kdataerr  = 0.20;
    kdata     = 1.42;
    kdataerr  = 0.58;
 
    htrange   = "H_{T} > 600 GeV";
    ymin      = 0.1;
    ymax      = 100;

    rebin     = 2;
  }

  //TCut metcut(Form("pfmet > %.0f",metcutval));

  const bool kFromMC = true;

  //---------------------------------
  // setup
  //---------------------------------

  gROOT->Reset();
  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);
  linelength = (width1+width2)*4+1;

  const unsigned int nsamples = mc.size();

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();
  //weight = TCut("1");

  cout << "jet cut                 : " << jetcut.GetTitle() << endl;
  //cout << "met cut                 : " << metcut.GetTitle() << endl;

  //-------------------
  // extract K
  //-------------------

  TH1F *hall  = new TH1F("hall" ,"hall" ,1,0,1);
  TH1F *hfail = new TH1F("hfail","hfail",1,0,1);
  TH1F *hpass = new TH1F("hpass","hpass",1,0,1);

  hall->Sumw2();
  hfail->Sumw2();
  hpass->Sumw2();

  TCut dilpt50("dilpt>50.");

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  if( kFromMC ){
    TChain *ch_mctot = new TChain();

    for( unsigned int i = 0 ; i < nsamples ; ++i ){ 
      ch_mctot->Add(mc.at(i));
    }

    ch_mctot->Draw("0.5>>hall"  , (sel+jetcut        ) *weight);
    ch_mctot->Draw("0.5>>hpass" , (sel+jetcut+dilpt50) *weight);
    ch_mctot->Draw("0.5>>hfail" , (sel+jetcut+!dilpt50)*weight);
  }
  else{
    data->Draw("0.5>>hall"  , (sel+jetcut        ));
    data->Draw("0.5>>hpass" , (sel+jetcut+dilpt50));
    data->Draw("0.5>>hfail" , (sel+jetcut+!dilpt50));
  }

  float all     = hall->GetBinContent(1);
  float pass    = hpass->GetBinContent(1);
  //float allerr  = hall->GetBinError(1);
  float passerr = hpass->GetBinError(1);
  float fail    = hfail->GetBinContent(1);
  float failerr = hfail->GetBinError(1);

  float k       = all/pass;
  float kerr    = k * sqrt( pow(failerr/fail,2) + pow(passerr/pass,2) );

  if( kFromMC )
    cout << "K (from MC)             : " << Form("%.2f / %.2f = %.2f +/- %.2f",all,pass,k,kerr) << endl;
  else
    cout << "K (from data)           : " << Form("%.2f / %.2f = %.2f +/- %.2f",all,pass,k,kerr) << endl;

  //------------------------
  // make MC histos
  //------------------------

  TH1F* hmet[nsamples];
  TH1F* hdilpt[nsamples];

  TH1F* hmet_mctot     = new TH1F("hmet_mctot"    , "hmet_mctot"   , nbins,xmin,xmax);     
  TH1F* hdilpt_mctot   = new TH1F("hdilpt_mctot"  , "hdilpt_mctot" , nbins,xmin,xmax);     

  hmet_mctot->Sumw2();    
  hdilpt_mctot->Sumw2();    

  TH1F* hmet_temp     = new TH1F("hmet_temp"    , "" , nbins,xmin,xmax);     
  TH1F* hdilpt_temp   = new TH1F("hdilpt_temp"  , "" , nbins,xmin,xmax);     

  hmet_temp->Sumw2();
  hdilpt_temp->Sumw2();

  for( unsigned int i = 0 ; i < nsamples ; ++i ){ 

    mc.at(i)->Draw(Form("TMath::Min(pfmet,%f)>>hmet_temp"    , xmax-0.01 ) , (sel+jetcut)*weight );
    mc.at(i)->Draw(Form("TMath::Min(dilpt,%f)>>hdilpt_temp"  , xmax-0.01 ) , (sel+jetcut)*weight );

    hmet[i]   = (TH1F*) hmet_temp->Clone(   Form("%s_hmety"    , mclabels.at(i) ) );
    hdilpt[i] = (TH1F*) hdilpt_temp->Clone( Form("%s_hdilpty"  , mclabels.at(i) ) );

    //hdilpt[i]->Sumw2();
    //hmet[i]->Sumw2();

    hdilpt[i]->Scale( k );
    hmet_temp->Reset();
    hdilpt_temp->Reset();
    
    hmet_mctot->Add(hmet[i]);
    hdilpt_mctot->Add(hdilpt[i]);
  }

  //-------------------
  // make data histos
  //-------------------
  
  TH1F* hmet_data          =  new TH1F("hmet_data",   "", nbins,xmin,xmax);
  TH1F* hdilpt_data        =  new TH1F("hdilpt_data", "", nbins,xmin,xmax); 
  
  hdilpt_data->Sumw2();
  hmet_data->Sumw2();
  
  data->Draw(Form("TMath::Min(pfmet,%f)>>hmet_data"   , xmax-0.01 ) ,  sel+jetcut );
  data->Draw(Form("TMath::Min(dilpt,%f)>>hdilpt_data" , xmax-0.01 ) ,  sel+jetcut );
  
  //cout << "Using K from MC" << endl;
  //float kdata = k;
  int ndprime = (int) hdilpt_data->Integral(hdilpt_data->FindBin(metcutval),10000);
  hdilpt_data->Scale(kdata);

  //--------------------
  // draw plot
  //--------------------
  
  TCanvas *cvic=new TCanvas("cvic","",800,600);
  cvic->cd();  
  
  TPad *plotpad = new TPad("plotpad","",0,0,1,0.9);
  plotpad->Draw();
  plotpad->cd();


  plotVictoryHist( hdilpt_mctot , hmet_mctot , hdilpt_data , hmet_data , 
		   "TITLE" , "E_{T}^{miss} (GeV)" , rebin , metcutval , htrange , ymin, ymax);

  //if( printgif ) c1->Print("plots/victory.gif");

  
  cvic->cd();
  TPad *titlepad = new TPad("titlepad","",0,0.9,1,1.0);
  titlepad->Draw(); 
  titlepad->cd();
  
  TLatex *title = new TLatex();
  title->SetTextSize(0.5);
  title->SetNDC();
  title->DrawLatex(0.18,0.15,"CMS Preliminary");
  title->DrawLatex(0.58,0.15,"#sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 0.98 fb^{-1}");
  //title->DrawLatex(0.58,0.15,"L_{int} = 976 pb^{-1},  #sqrt{s} = 7 TeV");

  //-------------------------------
  // get KC
  //-------------------------------

  pair<float,float> mcratio_pair = getMCratio(  hmet_mctot , hdilpt_mctot , metcutval );  
  float mcratio    = mcratio_pair.first;
  float mcratioerr = mcratio_pair.second;

  cout << "KC                      : " << Form("%.2f +/- %.2f",mcratio,mcratioerr) << endl;

  //-------------------------------
  // print table
  //-------------------------------

  
  cout << endl << endl;

  printLine(latex);
  printVictoryHeader();
  printLine(latex);

  for( unsigned int i = 0 ; i < nsamples ; ++i ){
    printVictoryRow( mclabels.at(i) , hmet[i] , hdilpt[i] , metcutval );
  }

  printLine(latex);
  printVictoryRow( "total MC" , hmet_mctot , hdilpt_mctot , metcutval );
  printLine(latex);

  //------------------------------
  // scale data by MC obs/pred
  //------------------------------

  printVictoryRow( "data" , hmet_data , hdilpt_data , metcutval , mcratio );
  printLine(latex);
  


  if( print ){
    
    cvic->Print(Form("../plots/victory_098fb_%s.pdf" , sigregion.c_str()));
    cvic->Print(Form("../plots/victory_098fb_%s.eps" , sigregion.c_str()));
    cvic->Print(Form("../plots/victory_098fb_%s.C"   , sigregion.c_str()));

  }


  //-------------------
  // do DY estimate
  //-------------------

  pair<float,float> dy_pair = doDYestimate( data , TCut(sel+jetcut) , metcutval );
  float dypred = dy_pair.first;
  float dyerr  = dy_pair.second;

  float pred        = kdata * kc * ( ndprime - dypred );
  float predstaterr = kdata * kc * sqrt(ndprime);
  float predsysterr = 0;

  predsysterr += pow( kdata    * kc    * dyerr              , 2 );
  predsysterr += pow( kdata    * kcerr * (ndprime - dypred) , 2 );
  predsysterr += pow( kdataerr * kc    * (ndprime - dypred) , 2 );
  predsysterr = sqrt(predsysterr);

  cout << endl << endl;
  cout << "N(D') " << ndprime << endl;
  cout << "N(DY) " << Form("%.2f%s%.2f",dypred,pm,dyerr)       << endl;
  cout << "K     " << Form("%.2f%s%.2f",kdata,pm,kdataerr)     << endl;
  cout << "KC    " << Form("%.2f%s%.2f",kc,pm,kcerr)           << endl;
  cout << "NP    " << Form("%.1f%s%.1f (stat)%s%.1f (syst)",pred,pm,predstaterr,pm,predsysterr) << endl;

  cout << endl;
  cout << ndprime << " & " 
       << Form("%.1f%s%.1f",dypred,pm,dyerr)   << " & " 
       << Form("%.1f%s%.1f",kdata,pm,kdataerr) << " & "
       << Form("%.1f%s%.1f",kc,pm,kcerr)       << " & "
       << Form("%.1f%s%.1f (stat)%s%.1f (syst)",pred,pm,predstaterr,pm,predsysterr) << " \\\\" << endl;


}


TH1F* getRandHist(TH1F* hin, string sample, int iter){

  TH1F* hout = (TH1F*) hin->Clone(Form("%s_clone_%i",hin->GetName(),iter));

  TRandom r;
  r.SetSeed(iter+1);

  for( int ibin = 1 ; ibin < hout->GetNbinsX() ; ibin++ ){
    float mean = hout->GetBinContent(ibin);
    float err  = hout->GetBinError(ibin);
    float rand;
    if( sample == "data" ) rand = r.Poisson(mean);
    else                   rand = r.Gaus(mean,err);
    if( rand < 0 ) rand = 0;
    hout->SetBinContent(ibin,rand);

    // cout << endl ;
    // cout << "bin " << ibin << endl;
    // cout << "mean " << mean << endl;
    // cout << "err  " << err  << endl;
    // cout << "rand " << rand << endl;
  }
  
  return hout;
  
}

void ABCDprime_closure(){

  pair<float,float> result;

  const unsigned int n = 9;
  //const unsigned int n = 15;

  float op[n];
  float operr[n];
  float met[n];
  float meterr[n];
  float ht[n];
  float hterr[n];
  
  for( unsigned int i = 0 ; i < n ; ++i ){

    ht[i]     = 300;
    hterr[i]  = 0.;
    met[i]    = 200. + 25. * i;
    meterr[i] = 0.;

    //ht[i]     = 300 + 50 * i;
    //hterr[i]  = 0.;
    //met[i]    = 200;
    //meterr[i] = 0;

    result = ABCDprime("../output/V00-01-03/highpt","ttpowheg",false,met[i],ht[i]);
    op[i]    = result.first;
    operr[i] = result.second;
  
    cout << "O/P " << Form("%.2f%s%.2f",op[i],pm,operr[i]) << endl;
 
  }

  TCanvas *can = new TCanvas();
  can->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  TGraphErrors *gr = new TGraphErrors(n,met,op,meterr,operr);
  gr->GetXaxis()->SetTitle("MET cut (GeV)");
  //TGraphErrors *gr = new TGraphErrors(n,ht,op,hterr,operr);
  //gr->GetXaxis()->SetTitle("H_{T} cut (GeV)");
  gr->GetYaxis()->SetTitle("obs / pred");

  gr->Draw("AP");

}

pair<float,float> ABCDprime( char* path , string sample , bool latex , float metcutval , float htcutval , float metcutvalmax ,  float htcutvalmax , bool print  ){

  //--------------------------------------------------
  // choose signal region plane: y_ht OR met_ht
  //--------------------------------------------------
 
  //string sigregion = "y_ht";
  string sigregion = "met_ht";

  //--------------------------------------------------
  // parameters for stat error estimation
  //--------------------------------------------------
 
  bool doStatError = true;
  int  nstats      = 20;

  //--------------------------------------------------
  // define the HT and MET variables
  //--------------------------------------------------

  char* htvar  = "htpf";
  char* metvar = "pfmet";

  //--------------------------------------------------
  // MET/HT cuts defining the signal region
  //--------------------------------------------------

  float metcut = 300.;
  if( metcutval > 0 ) metcut = metcutval;

  float htcut = 300;
  if( htcutval > 0 )  htcut  = htcutval;

  float metcutmax = 1.e10;
  if( metcutvalmax > 0 ) metcutmax = metcutvalmax;

  float htcutmax = 1.e10;
  if( htcutvalmax > 0 )  htcutmax  = htcutvalmax;

  float x1, x2, x3, y1, y2, y3;

  x1 = 125;
  x2 = 300;
  x3 = 300;
  
  y1 = 4.5;
  y2 = 8.5;
  y3 = 8.5;

  if( metcutval == 200 ){
    //x1 = 125;
    //x2 = 300;
    //x3 = 300;
  
    //y1 = 4.5;
    y2 = 6.5;
    y3 = 6.5;
  }

  //--------------------------------------------------
  // Define binning for f(y) and g(Ht) functions
  //--------------------------------------------------

  int   nbinsx = 115;
  int   nbinsy = 111;
  float xmax   = 3000.;
  float ymax   = 60.;

  //--------------------------------------------------
  // Initialize stuff
  //--------------------------------------------------
  
  gROOT->Reset();
  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);

  //--------------------------------------------------
  // preselection and weight
  //--------------------------------------------------

  TCut sel    = selection_TCut(highpt);
  //sel+="w1==1&&w2==1";
  //sel += TCut("(w1>0&&w2>0) && (w1<3&&w2<3)");

  //TCut weight = weight_TCut();
  TCut weight = "1";
  bool isData = ( sample == "data" );
  //if( !isData ) weight = TCut("weight*(1000./204.)");
  //if( !isData ) weight = TCut("weight");
  //if( !isData ) weight = weight_TCut();

  //--------------------------------------------------
  // choose sample
  //--------------------------------------------------

  TChain *ch = new TChain();

  if     ( sample == "data"     )    ch = data;
  else if( sample == "ttpowheg" )    ch = ttpowheg;
  else if( sample == "ttall"    )    ch = ttall;
  else{
    cout << "Unrecognized sample " << sample << ", quitting" << endl;
    exit(0);
  }

  // TH1F* h = new TH1F("h","",1,0,1);

  // ch->Draw("0.5>>h",(sel+"pfmet>275&&ht>300")*weight);
  // ch->Draw("0.5>>h",(sel+"pfmet>275&&ht>300")*weight);

  //--------------------------------------------------------
  // define regions in ABCD plane y vs. HT
  // lowx and lowy regions used for extracting f(y), g(HT)
  //--------------------------------------------------------
  
  //TCut ABCD  (Form("htpf>%.0f && y>%.1f",x1,y1));                               // event falls in ABCD region
  TCut A     (Form("htpf>%.0f && htpf<%.0f && y>%.1f",x1,x2,y3));               // region A
  TCut B     (Form("htpf>%.0f && htpf<%.0f && y>%.1f && y<%.1f",x1,x2,y1,y2));  // region B
  TCut C     (Form("htpf>%.0f && y>%.1f && y<%.1f",x3,y1,y2));                  // region C
  TCut D     (Form("htpf>%.0f && y>%.1f",x3,y3));                               // region D
  TCut lowx  (Form("htpf>%.0f && htpf<%.0f && y>%.1f",x1,x2,y1));               // region for f(y)  extraction
  TCut lowy  (Form("y>%.1f && y<%.1f && htpf>%.0f",y1,y2,x1));                  // region for g(HT) extraction

  //--------------------------------------------------
  // define signal and control regions
  //--------------------------------------------------
  
  TCut sig;
  TCut control;

  if     ( sigregion == "y_ht"   ){
    sig     = D;
    sig     = "y>8.5 && htpf>300.";
    cout << "Setting signal region " << sig.GetTitle() << endl;
    control = A||B||C;
    //control = lowx || lowy;
  }
  else if( sigregion == "met_ht" ){
    //sig = TCut(Form("htpf>%.0f && pfmet>%.0f",htcut,metcut));

    sig = TCut(Form("htpf>%.0f && htpf<%.0f && pfmet>%.0f && pfmet<%.0f",htcut,htcutmax,metcut,metcutmax));

    //sig = TCut("pfmet>275 && htpf>300 && htpf<600");  // region1
    //sig = TCut("pfmet>275 && htpf>600");              // region2
    //sig = TCut("pfmet>200 && pfmet<275 && htpf>600"); // region3

    //control = (A||B||C)&&!sig;
    //lowy = lowy + !sig;
    //control = (A||B||C);
    control = (lowx||lowy) && !sig;
    lowy    = lowy && !sig;
  }

  //--------------------------------------------------
  // cout regions
  //--------------------------------------------------

  cout << endl;
  //cout << "Region ABCD    : " << ABCD.GetTitle()    << endl;
  //cout << "Region A       : " << A.GetTitle()       << endl;
  //cout << "Region B       : " << B.GetTitle()       << endl;
  //cout << "Region C       : " << C.GetTitle()       << endl;
  //cout << "Region D       : " << D.GetTitle()       << endl;
  cout << "Low x          : " << lowx.GetTitle()    << endl;
  cout << "Low y          : " << lowy.GetTitle()    << endl;
  cout << "Sig region     : " << sig.GetTitle()     << endl;
  cout << "Control region : " << control.GetTitle() << endl;
  cout << endl;

  //--------------------------------------------------
  // declare histograms
  // habcd, hx, hy are filled from the baby
  // histos with 'new' are filled with the randomly
  // generated pseudo-events
  //--------------------------------------------------

  TH2F* habcd     = new TH2F( "habcd"    , "habcd"    , 1000   ,  0 , 1500 , 1000 , 0 , 30 ); 
  TH1F* hx        = new TH1F( "hx"       , "hx"       , nbinsx , x1 , xmax ); //  <<<--- this is g(HT)
  TH1F* hy        = new TH1F( "hy"       , "hy"       , nbinsy , y1 , ymax ); //  <<<--- this is f(y)

  TH2F* habcdnew  = new TH2F( "habcdnew" , "habcdnew" , 30     ,  0 , 1500 , 30 , 0 , 30 ); 
  TH1F* hxnew     = new TH1F( "hxnew"    , "hxnew"    , nbinsx , x1 , xmax );
  TH1F* hynew     = new TH1F( "hynew"    , "hynew"    , nbinsy , y1 , ymax );

  habcd->Sumw2();
  habcdnew->Sumw2();
  hx->Sumw2();
  hy->Sumw2();
  hxnew->Sumw2();
  hynew->Sumw2();

  //-------------------------------------------------
  // fill ABCD, f(y) and g(HT) functions
  //-------------------------------------------------

  ch->Draw(Form("min(y,%.1f):min(htpf,%.0f) >> habcd" , ymax , xmax ) , (sel)           * weight);
  //ch->Draw(Form("min(htpf,%.0f)             >> hx"    , xmax        ) , (sel+lowy+!sig) * weight);
  //ch->Draw(Form("min(y,%.1f)                >> hy"    , ymax        ) , (sel+lowx+!sig) * weight);
  ch->Draw(Form("min(htpf,%.0f)             >> hx"    , xmax        ) , (sel+lowy) * weight);
  ch->Draw(Form("min(y,%.1f)                >> hy"    , ymax        ) , (sel+lowx) * weight);

  cout << "f(y)  : " << hy->Integral() << " events" << endl;
  cout << "g(HT) : " << hx->Integral() << " events" << endl;

  //ncontrol -= nsig;
  //float ncontrol = ch->GetEntries(sel+ABCD);
  //float ncontrol = ch->GetEntries(sel+(lowy||lowx));

  //int bin = hx->FindBin(1000);
  //int bin = hx->FindBin(2000);
  //hx->SetBinContent(bin,1);

  //-------------------------------------------------
  // get true yields in signal and control regions
  //-------------------------------------------------

  TH1F *htemp = new TH1F("htemp","htemp",1,0,1);
  htemp->Sumw2();

  ch->Draw("0.5>>htemp",(sel+control)*weight);
  float ncontrol    = htemp->Integral();
  //float ncontrolerr = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+sig)*weight);
  float nsig    = htemp->Integral();
  float nsigerr = htemp->GetBinError(1);

  //float ncontrol = ch->GetEntries(sel+control);
  //float nsig     = ch->GetEntries(sel+sig);

  //-----------------------------------------
  // sample from templates, get prediction
  //-----------------------------------------

  const unsigned int npoints = 1000000;

  Float_t y_;
  Float_t pfmet_;
  Float_t htpf_;

  TTree* pseudoTree = new TTree("pseudoTree","pseudoTree");
  pseudoTree->Branch("y"      , &y_     ,  "y/F");
  pseudoTree->Branch("pfmet"  , &pfmet_ ,  "pfmet/F");
  pseudoTree->Branch("htpf"   , &htpf_  ,  "htpf/F");

  for( unsigned int i = 0 ; i < npoints ; i++ ){

    htpf_  = hx->GetRandom();
    y_     = hy->GetRandom();
    pfmet_ = y_ * sqrt(htpf_);
    
    hxnew->Fill(htpf_);
    hynew->Fill(y_);
    habcdnew->Fill(htpf_,y_);

    pseudoTree->Fill();
  }

  cout << "entries " << pseudoTree->GetEntries() << endl;
  cout << "non     " << pseudoTree->GetEntries(!sig+!control);

  float npseudo_sig     = pseudoTree->GetEntries(sig);
  float npseudo_control = pseudoTree->GetEntries(control);
  float R               = npseudo_sig / npseudo_control;

  float npred = R * ncontrol;

  cout << endl;
  cout << ncontrol      << " events in control regions" << endl;
  cout << nsig          << " events in signal region"   << endl;
  cout << nsig/ncontrol << " eff sig/control"           << endl;

  cout << endl;
  cout << "R = " << npseudo_sig << " / " <<  npseudo_control << " " << R << endl;
  cout << "prediction " << npred << endl;
 
  float delta_pred = 0;
  float nprederr   = 0;

  TH1F* hxrand[nstats];
  TH1F* hyrand[nstats];

  if( doStatError ){

    cout << endl << "Doing stat error" << endl;

    for( int istat = 0 ; istat < nstats ; istat++ ){
      
      pseudoTree->Reset();
      
      hxrand[istat] = getRandHist(hx,sample,istat);
      hyrand[istat] = getRandHist(hy,sample,istat);

      for( unsigned int i = 0 ; i < npoints ; i++ ){
	
	htpf_  = hxrand[istat]->GetRandom();
	y_     = hyrand[istat]->GetRandom();
	pfmet_ = y_ * sqrt(htpf_);

	pseudoTree->Fill();
      }
            
      float stat_npseudo_sig     = (float) pseudoTree->GetEntries(sig);
      float stat_npseudo_control = (float) pseudoTree->GetEntries(control);
      float stat_R               = stat_npseudo_sig / stat_npseudo_control;
      
      float stat_npred = stat_R * ncontrol;
      
      cout << istat << " pred: " << Form("%.2f",stat_npred) << endl;

      delta_pred += pow(stat_npred - npred,2);
    }

    delta_pred /= nstats;
    delta_pred = sqrt(delta_pred);
    
    nprederr = delta_pred;
    cout << "Stat error: " << Form("%.2f",nprederr) << endl;
  }

  float op    = nsig/npred;
  float operr = op * sqrt( pow(nsigerr/nsig,2) + pow(nprederr/npred,2) );

  cout << endl;
  cout << "observed  : " << Form("%.1f%s%.1f", nsig,pm,nsigerr)                  << endl;
  cout << "predicted : " << Form("%.1f%s%.1f", npred,pm,nprederr)                << endl;
  cout << "o/p       : " << Form("%.2f%s%.2f", op , pm, operr) << endl;


  /*
  TCanvas *statcan = new TCanvas();
  statcan->cd();
  hx->Draw("hist");
  hxrand[0]->SetMarkerColor(2);
  hxrand[0]->SetLineColor(2);
  hxrand[1]->SetMarkerColor(4);
  hxrand[1]->SetLineColor(4);
  hxrand[0]->Draw("sameE1");
  hxrand[1]->Draw("sameE1");
  */

  //---------------------
  // draw stuff
  //---------------------

  TCanvas *abcd_can = new TCanvas();
  abcd_can->cd();

  habcd->GetYaxis()->SetTitle("y [ #sqrt{GeV} ]");
  habcd->GetYaxis()->SetTitleOffset(1);
  habcd->GetXaxis()->SetTitle("H_{T} [ GeV ]");
  
  habcd->SetMarkerColor(4);
  habcd->SetMarkerSize(0.3);
  habcd->Draw();
  
  //habcd->Draw("box");
  //habcd->SetLineColor(4);
  //habcd->SetMarkerColor(0);
  //habcd->SetFillColor(0);

  //TLatex *text=new TLatex();
  //text->SetTextSize(0.05);
  //text->SetTextColor(1);
      
  // drawSquare(x1,y1,x2,y2);
  // drawSquare(x3,y1,1500,y2);
  // drawSquare(x1,y3,x2,30);
  // drawSquare(x3,y3,1500,30);

  //-------------------------------------
  // draw f(y), g(HT) regions
  //-------------------------------------

  drawSquare(x1,y1,1500,y2,2);
  drawSquare(x1,y1,x2,30,8);
  
  TBox box;
  box.SetLineColor(2);
  box.SetLineWidth(3);
  box.SetFillColor(2);
  box.SetFillStyle(3004);
  box.DrawBox(x1,y1,1500,y2);
  box.SetLineColor(8);
  box.SetLineWidth(3);
  box.SetFillColor(8);
  box.SetFillStyle(3005);
  box.DrawBox(x1,y1,x2,30);

  /*
  text->DrawLatex( x1+50 , y1+1.5 , "B" );
  text->DrawLatex( x1+50 , y3+10  , "A" );
  text->DrawLatex( x3+50 , y1+1.5 , "C" );
  text->DrawLatex( x3+50 , y3+10  , "D" );
  */

  if( sigregion == "met_ht" ){
    TF1* fmet = new TF1("fmet","[0]/sqrt(x)",htcut,1500);
    fmet->SetParameter(0,metcut);
    fmet->SetLineColor(1);
    fmet->SetLineWidth(4);
    fmet->Draw("same");

    TLine line;
    line.SetLineColor(1);
    line.SetLineWidth(4);
    line.DrawLine(htcut,metcut/sqrt(htcut),htcut,30);    
    line.DrawLine(htcut,30,1500,30);    
    line.DrawLine(1500,metcut/sqrt(1500),1500,30);    

  }

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //text->DrawLatex(0.52,0.85,"CMS Preliminary");
  //text->DrawLatex(0.52,0.80,"#sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 0.98 fb^{-1}");
  //text->DrawLatex(0.52,0.75,"Events with ee/#mu#mu/e#mu");
  text->DrawLatex(0.52,0.85,"CMS Simulation");
  text->DrawLatex(0.52,0.80,"#sqrt{s} = 7 TeV, t#bar{t}#rightarrow#font[12]{l^{+}}#font[12]{l^{-}}");
  text->DrawLatex(0.52,0.75,"Events with ee/#mu#mu/e#mu");


  TCanvas *projcan = new TCanvas("projcan","",1200,600);
  projcan->Divide(2,1);

  projcan->cd(1);
  gPad->SetLogy();
  hx->Draw("hist");
  hx->GetXaxis()->SetTitle("H_{T} (GeV)");
  hxnew->Scale( hx->Integral() / hxnew->Integral() );
  hxnew->SetLineColor(2);
  hxnew->SetMarkerColor(2);
  //hxnew->Draw("sameE1");
  text->DrawLatex(0.52,0.85,"CMS Simulation");
  text->DrawLatex(0.52,0.80,"#sqrt{s} = 7 TeV, t#bar{t}#rightarrow#font[12]{l^{+}}#font[12]{l^{-}}");
  text->DrawLatex(0.52,0.75,"Events with ee/#mu#mu/e#mu");


  projcan->cd(2);
  gPad->SetLogy();
  hy->Draw("hist");
  hy->GetXaxis()->SetTitle("y (GeV^{1/2})");
  hynew->Scale( hy->Integral() / hynew->Integral() );
  hynew->SetLineColor(2);
  hynew->SetMarkerColor(2);
  //hynew->Draw("sameE1");
  text->DrawLatex(0.52,0.85,"CMS Simulation");
  text->DrawLatex(0.52,0.80,"#sqrt{s} = 7 TeV, t#bar{t}#rightarrow#font[12]{l^{+}}#font[12]{l^{-}}");
  text->DrawLatex(0.52,0.75,"Events with ee/#mu#mu/e#mu");


  TCanvas *abcdnewcan = new TCanvas();
  abcdnewcan->cd();
  habcdnew->Draw("box");
  // TCanvas *can1 = new TCanvas();
  // can1->cd();
  // hy->Draw("hist");

  // TCanvas *can2 = new TCanvas();
  // can2->cd();
  // hx->Draw("hist");

  char* sigchar = "highmet";
  if( metcutval == 200 ) sigchar = "highht";

  if( print ){
    abcd_can->Print(Form("../plots/abcdprime_098fb_%s.pdf" , sigchar));
    abcd_can->Print(Form("../plots/abcdprime_098fb_%s.eps" , sigchar));
    abcd_can->Print(Form("../plots/abcdprime_098fb_%s.C"   , sigchar));
  }

  return make_pair( op , operr );

}

void njets( char* path , bool printplot = false ){

  deleteHistos();
  initialize(path);
  
  //char* var = "njetsuncor";
  char* var = "npfjets";
  //char* var = "htuncor";

  char* sample = "data";
  TChain *ch = data;

  TH1F* data_v1  = getHist( ch , var , TCut("dilmass>81&&dilmass<101&&(leptype==0||leptype==1)&&ndavtx<4")           , "data_v1"  , 10 , 0 , 500 );
  TH1F* data_v2  = getHist( ch , var , TCut("dilmass>81&&dilmass<101&&(leptype==0||leptype==1)&&ndavtx>3&&ndavtx<6") , "data_v2"  , 10 , 0 , 500 );
  TH1F* data_v3  = getHist( ch , var , TCut("dilmass>81&&dilmass<101&&(leptype==0||leptype==1)&&ndavtx>5")           , "data_v3"  , 10 , 0 , 500 );

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
  leg->AddEntry(data_v2,"4 < nvtx #leq 5");
  leg->AddEntry(data_v3,"nvtx > 5");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  if( printplot ) c1->Print(Form("../plots/%s_%s.gif",sample,var));
}


void PU( char* path , bool printplot = false ){

  deleteHistos();
  initialize(path);
  
  char* var    = "pfmet";
  int   nbins  = 10;
  float xmax   = 300.;
  float cutval = 30;

  // char* var    = "npfjets";
  // int   nbins  = 10;
  // float xmax   = 10.;
  // float cutval = 2;

  //char* var    = "htpf";
  //int   nbins  = 20;
  //float xmax   = 1000.;
  //float cutval = 100;

  char* sample = "data";
  TChain *ch = data;

  TCut base("leptype ==2 && lep1.pt() > 20 && lep2.pt() > 20");
  TCut vtx1("ndavtx<4");
  TCut vtx2("ndavtx>3&&ndavtx<6");
  TCut vtx3("ndavtx>5");
  TCut hi(Form("%s>%.0f",var,cutval));
  TCut lo(Form("%s<%.0f",var,cutval));

  TH1F* data_v1  = getHist( ch , var , TCut(base+vtx1) , "data_v1"  , nbins , 0 , xmax );
  TH1F* data_v2  = getHist( ch , var , TCut(base+vtx2) , "data_v2"  , nbins , 0 , xmax );
  TH1F* data_v3  = getHist( ch , var , TCut(base+vtx3) , "data_v3"  , nbins , 0 , xmax );

  TH1F* nvtxlo   = getHist( ch , "ndavtx" , TCut(base+lo) , "nvtx"     , 20 , 0 , 20 );
  TH1F* nvtxhi   = getHist( ch , "ndavtx" , TCut(base+hi) , "nvtx50"   , 20 , 0 , 20 );

  data_v2->SetLineColor(2);
  data_v2->SetMarkerColor(2);
  data_v2->SetMarkerStyle(23);

  data_v3->SetLineColor(4);
  data_v3->SetMarkerColor(4);
  data_v3->SetMarkerStyle(25);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  //gPad->SetLogy();
  data_v1->GetXaxis()->SetTitle(var);
  data_v1->DrawNormalized();
  data_v2->DrawNormalized("same");
  data_v3->DrawNormalized("same");

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(data_v1,"nvtx #leq 3");
  leg->AddEntry(data_v2,"4 < nvtx #leq 5");
  leg->AddEntry(data_v3,"nvtx > 5");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  int bin = data_v1->FindBin(cutval);

  cout << "nvtx < 4      : " << data_v1->Integral(bin,100) << " / " << data_v1->Integral() << " = " << data_v1->Integral(bin,100)/data_v1->Integral() << endl;
  cout << "nvtx 4-5      : " << data_v2->Integral(bin,100) << " / " << data_v2->Integral() << " = " << data_v2->Integral(bin,100)/data_v2->Integral() << endl;
  cout << "nvtx > 5      : " << data_v3->Integral(bin,100) << " / " << data_v3->Integral() << " = " << data_v3->Integral(bin,100)/data_v3->Integral() << endl;

  TCanvas *c2 = new TCanvas();
  c2->cd();
  nvtxlo->GetXaxis()->SetTitle("nDAvtx");
  nvtxlo->DrawNormalized("hist");
  nvtxhi->SetLineColor(2);
  nvtxhi->SetMarkerColor(2);
  nvtxhi->DrawNormalized("sameE1");

  TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
  leg2->AddEntry(nvtxlo,lo.GetTitle());
  leg2->AddEntry(nvtxhi,hi.GetTitle());
  leg2->SetBorderSize(1);
  leg2->SetFillColor(0);
  leg2->Draw();

  if( printplot ){ 
    c1->Modified();
    c1->Update();
    c1->Print(Form("../plots/%s_%s.gif",sample,var));
    c2->Modified();
    c2->Update();
    c2->Print(Form("../plots/%s_%s_nvtx.gif",sample,var));
  }
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

  cout << "Trig eff correction     : 100% (ee), 90% (mm), 95% (em)" << endl;
  
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
  float x2=600;
  float x3=600;
  float x4=1500;
  
  float y1 = 4.5;
  float y2 = 8.5;
  float y3 = 8.5;
  float y4 = 30;
  
  //y3 = 13;
  //x3=600;

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
  //TH1F* hyield  = new TH1F("hyield","hyield",4,0,4);
  //hyield->Sumw2();
  //TH1F* hyield_mctot = new TH1F("hyield_mctot","hyield_mctot",4,0,4);

  for( unsigned int i = 0 ; i < mc.size() ; ++i ){
    //hmc[i] = doABCD( mc[i] , mclabels[i] , sel , weight , A , B , C , D , hyield );
    hmc[i] = doABCD( mc[i] , mclabels[i] , sel , weight , A , B , C , D );
    mctot->Add(mc[i]);
    //hyield_mctot->Add(hyield);
  }

  printLine(latex);

  //printABCDRow(hyield_mctot);
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
  text->DrawLatex(x1+50,y3+1.5,"A");
  text->DrawLatex(x3+50,y1+1.5,"C");
  text->DrawLatex(x3+50,y3+1.5,"D");

  
  text->SetNDC();
  text->SetTextSize(0.037);
  text->DrawLatex(0.35,0.85,"CMS");
  text->DrawLatex(0.35,0.80,"0.98 fb^{-1} at #sqrt{s} = 7 TeV");
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

  //bool y_vary = true;
  //const unsigned int n = 16;

  bool y_vary = false;
  const unsigned int n = 10;

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut sel    = selection_TCut(highpt);
  //TCut weight = weight_TCut();
  TCut weight = TCut("1");//weight_TCut();
  
  float x1=125;
  float x2=300;
  float x3=300;
  //float x4=1500;
  
  float y1 = 4.5;
  float y2 = 8.5;
  float y3 = 8.5;
  //float y4 = 30;

  TCut A("");
  TCut B("");
  TCut C("");
  TCut D("");
    


  float x[n];
  float xerr[n];
  float mg[n];
  float pow[n];
  float mgerr[n];
  float powerr[n];

  for( unsigned int i = 0 ; i < n ; ++i ){

    if( y_vary ){
      float ycut = 6.5 + i*0.5;

      //y2 = ycut;
      y3 = ycut;

      x[i]    = ycut;
      xerr[i] = 0.0;
    }

    else{
      float htcut = 200 + 50 * i;
      
      //x2 = htcut;
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

  if( y_vary ) gr_mg->GetXaxis()->SetTitle("y cut (GeV^{1/2})");
  else         gr_mg->GetXaxis()->SetTitle("H_{T} cut (GeV^{1/2})");
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

TH1F* getABCDYield( TChain *ch , char* sample , TCut sel , TCut weight ){


  TH1F *htempee = new TH1F("htempee","htempee",1,0,1);
  TH1F *htempmm = new TH1F("htempmm","htempmm",1,0,1);
  TH1F *htempem = new TH1F("htempem","htempem",1,0,1);

  htempee->Sumw2();
  htempmm->Sumw2();
  htempem->Sumw2();

  TCut eetype("leptype==0");
  TCut mmtype("leptype==1");
  TCut emtype("leptype==2");

  TCut eeeff("1.00");
  TCut mmeff("0.90");
  TCut emeff("0.95");

  if( strcmp(sample,"data")==0 ){
    eeeff=TCut("1.00");
    mmeff=TCut("1.00");
    emeff=TCut("1.00");
  }

  TCanvas *mytemp = new TCanvas();  
  ch->Draw("0.5>>htempee",(sel+eetype)*weight*eeeff);
  ch->Draw("0.5>>htempmm",(sel+mmtype)*weight*mmeff);
  ch->Draw("0.5>>htempem",(sel+emtype)*weight*emeff);
  delete mytemp;

  TRandom r;
  int rand = (int)r.Uniform(0,10000);

  TH1F *htot = (TH1F*) htempee->Clone(Form("hist_%i",rand));
  htot->Add(htempmm);
  htot->Add(htempem);

  delete htempee;
  delete htempmm;
  delete htempem;

  return htot;

}

//TH2F* doABCD( TChain *ch , char* label , TCut sel , TCut weight , TCut A , TCut B , TCut C , TCut D , TH1F* hyield ){
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

  //TH1F* hA = getABCDYield( ch , label , TCut(sel+A) , weight );
  //TH1F* hB = getABCDYield( ch , label , TCut(sel+B) , weight );
  //TH1F* hC = getABCDYield( ch , label , TCut(sel+C) , weight );
  //TH1F* hD = getABCDYield( ch , label , TCut(sel+D) , weight );

  float nA = hA->GetBinContent(1);
  float nB = hB->GetBinContent(1);
  float nC = hC->GetBinContent(1);
  float nD = hD->GetBinContent(1);

  float eA = hA->GetBinError(1);
  float eB = hB->GetBinError(1);
  float eC = hC->GetBinError(1);
  float eD = hD->GetBinError(1);

  // if( strcmp(label,"data") != 0 ){
  //   nA *= 0.95;
  //   nB *= 0.95;
  //   nC *= 0.95;
  //   nD *= 0.95;
  //   eA *= 0.95;
  //   eB *= 0.95;
  //   eC *= 0.95;
  //   eD *= 0.95;
  // }

  printABCDRow( label , nA , nB , nC , nD , eA , eB , eC , eD );

  delete hA;
  delete hB;
  delete hC;
  delete hD;

  // hyield->Reset();

  // hyield->SetBinContent(1,nA);
  // hyield->SetBinContent(2,nB);
  // hyield->SetBinContent(3,nC);
  // hyield->SetBinContent(4,nD);

  // hyield->SetBinError(1,eA);
  // hyield->SetBinError(2,eB);
  // hyield->SetBinError(3,eC);
  // hyield->SetBinError(4,eD);

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
    sD     << Form("%.1f %s %.1f",  D    , pm , dD       );
    sPred  << Form("%.1f %s %.1f",  pred , pm , prederr  );
    sRatio << Form("%.2f %s %.2f", ratio , pm , ratioerr );
  }
  
  cout  << delimstart << setw(width1) << sample        << setw(width2)
        << delim      << setw(width1) << sA.str()      << setw(width2)
        << delim      << setw(width1) << sB.str()      << setw(width2)
        << delim      << setw(width1) << sC.str()      << setw(width2)
        << delim      << setw(width1) << sD.str()      << setw(width2)
        << delim      << setw(width1) << sPred.str()   << setw(width2) 
        << delim      << setw(width1) << sRatio.str()  << setw(width2)  
	<< delimend << endl;
}


void printABCDHeader(){

  cout  << delimstart << setw(width1) << "sample"   << setw(width2) 
        << delim      << setw(width1) << "A"        << setw(width2) 
        << delim      << setw(width1) << "B"        << setw(width2) 
        << delim      << setw(width1) << "C"        << setw(width2) 
        << delim      << setw(width1) << "D"        << setw(width2) 
        << delim      << setw(width1) << "A \\times B / C"     << setw(width2) 
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

  bool  combine4  = false;
  bool  residual  = false;
  bool  log       = false;
  char* flavor    = "all";


  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  //----------------------------------------
  // read in selection and weight
  //----------------------------------------

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  //sel = sel + "y>8.5";
  //sel = sel + "y > 8.5 && ht > 300";

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  if( highpt ){
    //vars.push_back("mt2");        xt.push_back("MT2");			n.push_back(10); xi.push_back(0.);   xf.push_back(150.);
    vars.push_back("mt2j");        xt.push_back("MT2J");			n.push_back(18); xi.push_back(140.);   xf.push_back(500.);
    //vars.push_back("mlljj");     xt.push_back("m(#mu#mujj) (GeV)");      n.push_back(25); xi.push_back(350.); xf.push_back(1850.);
    //vars.push_back("pfmet");     xt.push_back("pfmet (GeV)");      n.push_back(20); xi.push_back(0.); xf.push_back(200.);
    // vars.push_back("pfmet");     xt.push_back("E_{T}^{miss} (GeV)");       n.push_back(30); xi.push_back(0.); xf.push_back(300.);
    // //vars.push_back("y");         xt.push_back("y #equiv MET  /  #sqrt{H_{T}} (GeV^{1/2})");    n.push_back(20); xi.push_back(0.); xf.push_back(20.);
    // vars.push_back("htpf");      xt.push_back("H_{T} (GeV)");      n.push_back(20); xi.push_back(0.); xf.push_back(1000.);
    // vars.push_back("dilmass");   xt.push_back("m(ll) (GeV)");      n.push_back(20); xi.push_back(1.); xf.push_back(301.);
    // //vars.push_back("lep1.eta()");   xt.push_back("#eta(lep1)");      n.push_back(50); xi.push_back(-3.); xf.push_back(3.);
    // //vars.push_back("lep2.eta()");   xt.push_back("#eta(lep2)");      n.push_back(50); xi.push_back(-3.); xf.push_back(3.);
    // //vars.push_back("lep1.pt()");   xt.push_back("pt(lep1)");         n.push_back(50); xi.push_back(0.); xf.push_back(100.);
    // //vars.push_back("lep2.pt()");   xt.push_back("pt(lep2)");         n.push_back(50); xi.push_back(0.); xf.push_back(100.);
    // vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(20); xi.push_back(0.); xf.push_back(300.);
    //vars.push_back("npfjets");   xt.push_back("njets");            n.push_back(10); xi.push_back(0.); xf.push_back(10.);
    //vars.push_back("ndavtx");      xt.push_back("nDAVertices");      n.push_back(20); xi.push_back(0.); xf.push_back(20.);
    //vars.push_back("htoffset");  xt.push_back("L1Offset-H_{T} (GeV)");    n.push_back(10); xi.push_back(0.); xf.push_back(1000.);
    //vars.push_back("htuncor");   xt.push_back("uncorrected H_{T} (GeV)");    n.push_back(10); xi.push_back(0.); xf.push_back(1000.);
  }else{
    vars.push_back("mt2");        xt.push_back("MT2");			n.push_back(10); xi.push_back(0.);   xf.push_back(150.);
    // vars.push_back("pfmet");     xt.push_back("pfmet (GeV)");      n.push_back(5); xi.push_back(0.); xf.push_back(250.);
    // vars.push_back("htpf");      xt.push_back("H_{T} (GeV)");      n.push_back(7); xi.push_back(0.); xf.push_back(700.);
    // vars.push_back("dilmass");   xt.push_back("m(ll) (GeV)");      n.push_back(5); xi.push_back(0.); xf.push_back(100.);
    // vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(5); xi.push_back(0.); xf.push_back(100.);
  } 

  //vars.push_back("njets");     xt.push_back("njets");            n.push_back(6);  xi.push_back(0);  xf.push_back(6);
  //vars.push_back("dilep.mass()");     xt.push_back("dilmass (GeV)");    n.push_back(100);  xi.push_back(0);  xf.push_back(200);
  
  const unsigned int nvars = vars.size();
  
  TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];
  int canCounter = -1;

  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    //if( ivar < 2 ) log = true;
    //else           log = false;

    //log = false;

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

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , !combine4 , log , flavor);

    if( printgif ) can[ivar]->Print(Form("../plots/%s.pdf",vars[ivar]));
  } 
}


void printOFRow( char* sample, float nee, float nmm, float nem, float eeerr, float mmerr, float emerr ){

  //float tot    = nee + nmm + nem;
  //float toterr = sqrt( nee*nee + nmm*nmm + nem*nem );
  
  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;
    
  if( sample == "Observed"){
    see     << nee;
    smm     << nmm;
    sem     << nem;
    //stot    << tot;
  }
  else{
    see     << Form("%.1f %s %.1f",  nee  , pm , eeerr  );
    smm     << Form("%.1f %s %.1f",  nmm  , pm , mmerr  );
    if( nem > -0.1 ) sem     << Form("%.1f %s %.1f",  nem  , pm , emerr  );
    else             sem     << "";
    //stot    << Form("%.1f %s %.1f",  tot , pm , toterr );
  }
  
  cout  << delimstart << setw(width1) << sample         << setw(width2)
        << delim      << setw(width1) << see.str()      << setw(width2)
        << delim      << setw(width1) << smm.str()      << setw(width2)
        << delim      << setw(width1) << sem.str()      << setw(width2)
      //<< delim      << setw(width1) << stot.str()     << setw(width2)  
	<< delimend << endl;

}

void printOFHeader(){
  cout  << delimstart << setw(width1) << ""      << setw(width2)
        << delim      << setw(width1) << ee      << setw(width2)
        << delim      << setw(width1) << mm      << setw(width2)
        << delim      << setw(width1) << em      << setw(width2)
      //<< delim      << setw(width1) << "total" << setw(width2)  
	<< delimend << endl;
}


void doOFSubtraction( char* path , string sample = "data", bool latex = false ){

  bool isData = sample == "data";

  //----------------------
  // R=eff(mu)/eff(e)
  //----------------------

  const float R10 = 1.28;
  const float R20 = 1.08;

  //----------------------
  // fake yields
  //----------------------

  float nee1010fake = 0.00;
  float nmm1010fake = 0.37;
  float nem1010fake = 0.89;
    
  float nee105fake  = 0.00;
  float nmm105fake  = 1.71;
  float nem105fake  = 1.07;

  //----------------------
  // trigger efficiencies
  //----------------------
  
  float trigee = 1.00;
  float trigmm = 0.90;
  float trigem = 0.95;

  if( !isData ){
    nee1010fake = 0.;
    nmm1010fake = 0.;
    nem1010fake = 0.;

    nee105fake = 0.;
    nmm105fake = 0.;
    nem105fake = 0.;

    trigee = 1;
    trigmm = 1;
    trigem = 1;
  }


  gROOT->Reset();
  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);
  initSymbols(latex);
  linelength  = (width1+width2)*4+1;
   
  TCut sel    = selection_TCut(highpt);
  //TCut weight = weight_TCut();
  TCut weight("weight * (1000./204.)");
  if( isData ) weight = TCut("1");
  //TCut weight("1");
  
  sel = sel + "y>8.5&&htpf>300";
  cout << "SELECTION: " << sel.GetTitle() << endl;

  TCut eetype("leptype==0");
  TCut mmtype("leptype==1");
  TCut emtype("leptype==2");
  TCut pt2010("lep1.pt()>20 && lep2.pt()>10");
  TCut pt1010("lep1.pt()>10 && lep2.pt()>10");
  TCut pt105 ("lep1.pt()>10 && lep2.pt()>5 ");

  TChain *ch = new TChain();
  TCut myweight;

  if( sample == "data" ){
    cout << "Doing data" << endl;
    ch = data;
    myweight = TCut("1");
  }
  else if( sample == "ttall" ){
    cout << "Doing ttbar madgraph" << endl;
    ch = ttall;
    myweight = weight;
  }
  else if( sample == "ttpowheg" ){
    cout << "Doing ttbar madgraph" << endl;
    ch = ttpowheg;
    myweight = weight;
  }
  else{
    cout << "Unrecognized sample " << sample << ", quitting" << endl;
    exit(0);
  }

  TH1F* htemp = new TH1F("htemp","htemp",1,0,1);
  htemp->Sumw2();
  
  // float nee1010 = ch->GetEntries(sel+eetype+pt1010+!pt2010);
  // float nmm1010 = ch->GetEntries(sel+mmtype+pt1010+!pt2010);
  // float nem1010 = ch->GetEntries(sel+emtype+pt1010+!pt2010);

  // float nee105 = ch->GetEntries(sel+eetype+pt105+!pt1010);
  // float nmm105 = ch->GetEntries(sel+mmtype+pt105+!pt1010);
  // float nem105 = ch->GetEntries(sel+emtype+pt105+!pt1010);

  ch->Draw("0.5>>htemp",(sel+eetype+pt1010+!pt2010)*weight);
  float nee1010    = htemp->GetBinContent(1);
  float nee1010err = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+mmtype+pt1010+!pt2010)*weight);
  float nmm1010    = htemp->GetBinContent(1);
  float nmm1010err = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+emtype+pt1010+!pt2010)*weight);
  float nem1010    = htemp->GetBinContent(1);
  float nem1010err = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+eetype+pt105+!pt1010)*weight);
  float nee105     = htemp->GetBinContent(1);
  float nee105err  = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+mmtype+pt105+!pt1010)*weight);
  float nmm105     = htemp->GetBinContent(1);
  float nmm105err  = htemp->GetBinError(1);

  ch->Draw("0.5>>htemp",(sel+emtype+pt105+!pt1010)*weight);
  float nem105     = htemp->GetBinContent(1);
  float nem105err  = htemp->GetBinError(1);


 
  if( isData ){
    TCut fakeweight("weight");

    TCanvas *ctemp = new TCanvas();
    ctemp->cd();

    datasf->Draw("0.5>>htemp",(sel+eetype+pt1010+!pt2010)*fakeweight);
    nee1010fake = htemp->Integral();

    datasf->Draw("0.5>>htemp",(sel+mmtype+pt1010+!pt2010)*fakeweight);
    nmm1010fake = htemp->Integral();

    datasf->Draw("0.5>>htemp",(sel+emtype+pt1010+!pt2010)*fakeweight);
    nem1010fake = htemp->Integral();

    datasf->Draw("0.5>>htemp",(sel+eetype+pt105+!pt1010)*fakeweight);
    nee105fake = htemp->Integral();

    datasf->Draw("0.5>>htemp",(sel+mmtype+pt105+!pt1010)*fakeweight);
    nmm105fake = htemp->Integral();

    datasf->Draw("0.5>>htemp",(sel+emtype+pt105+!pt1010)*fakeweight);
    nem105fake = htemp->Integral();

    delete ctemp;
  }


  cout << endl;
  cout << "R(10-20)    : " << R10 << endl;
  cout << "R(>20)      : " << R20 << endl;

  cout << endl;
  cout << "N(ee) 10,10 : " << Form("%.1f%s%.1f",nee1010,pm,nee1010err) << endl;
  cout << "N(mm) 10,10 : " << Form("%.1f%s%.1f",nmm1010,pm,nmm1010err) << endl;
  cout << "N(em) 10,10 : " << Form("%.1f%s%.1f",nem1010,pm,nem1010err) << endl;

  cout << endl;
  cout << "N(ee) 10,5  : " << Form("%.1f%s%.1f",nee105,pm,nee105err) << endl;
  cout << "N(mm) 10,5  : " << Form("%.1f%s%.1f",nmm105,pm,nmm105err) << endl;
  cout << "N(em) 10,5  : " << Form("%.1f%s%.1f",nem105,pm,nem105err) << endl;

  cout << endl;
  cout << "N(ee) fake 10,10 : " << nee1010fake << endl;
  cout << "N(mm) fake 10,10 : " << nmm1010fake << endl;
  cout << "N(em) fake 10,10 : " << nem1010fake << endl;

  cout << endl;
  cout << "N(ee) fake 10,5  : " << nee105fake << endl;
  cout << "N(mm) fake 10,5  : " << nmm105fake << endl;
  cout << "N(em) fake 10,5  : " << nem105fake << endl;
    
  float predee1010 = (0.5/R10) * (nem1010-nem1010fake) * (trigee/trigem);
  float predmm1010 = (0.5*R10) * (nem1010-nem1010fake) * (trigmm/trigem);
  float predmm105  = R20       * (nem105-nem105fake)   * (trigmm/trigem);

  float predee1010err = (0.5/R10) * sqrt(nem1010+pow(0.5*nem1010fake,2)) * (trigee/trigem);
  float predmm1010err = (0.5*R10) * sqrt(nem1010+pow(0.5*nem1010fake,2)) * (trigmm/trigem);
  float predmm105err  = R20       * sqrt(nem105 +pow(0.5*nem105fake,2))  * (trigmm/trigem);

  cout << endl;
  cout << "N(ee) pred 10,10 : " << predee1010 << endl;
  cout << "N(mm) pred 10,10 : " << predmm1010 << endl;
  cout << "N(mm) pred 10,5  : " << predmm105  << endl;

  // cout << endl;
  // cout << "ee obs         : " << nee1010    << endl;
  // cout << "ee fake        : " << nee1010fake << endl;
  // cout << "ee OF          : " << Form("%.1f +/- %.1f",predee1010,predee1010err) << endl;
  // cout << "ee tot pred    : " << Form("%.1f +/- %.1f",predee1010+nee1010fake,predee1010err) << endl;

  // cout << endl;
  // cout << "mm obs         : " << nmm105+nmm1010 << endl;
  // cout << "mm fake        : " << nmm105fake+nmm1010fake << endl;
  // cout << "mm OF          : " << Form("%.1f +/- %.1f",predmmtot,predmmtoterr)  << endl;
  // cout << "mm tot pred    : " << Form("%.1f +/- %.1f",predmmtot+nmm105fake+nmm1010fake,predmmtoterr)  << endl;

  // cout << endl;
  // cout << "mm obs  10,10  : " << nmm1010    << endl;
  // cout << "mm pred 10,10  : " << Form("%.1f +/- %.1f",predmm1010,predmm1010err) << endl;

  // cout << endl;
  // cout << "mm obs  10,5   : " << nmm105     << endl;
  // cout << "mm pred 10,5   : " << Form("%.1f +/- %.1f",predmm105,predmm105err)  << endl;

  float neefake = nee105fake + nee1010fake;
  float nmmfake = nmm105fake + nmm1010fake;
  float nemfake = nem105fake + nem1010fake;

  float predeetot  = predee1010;
  float predmmtot  = predmm1010 + predmm105;

  float predeetoterr  = predee1010err;
  float predmmtoterr  = sqrt( pow(predmm1010err,2) + pow(predmm105err,2) );

  cout << endl << endl;

  printLine(latex);
  printOFHeader();

  if( isData ){
    printLine(latex);
    printOFRow( "Fakes"    , neefake    , nmmfake   , nemfake , 0.5 * neefake , 0.5 * nmmfake , 0.5 * nemfake );
    
    printOFRow( "OF pred"  , predeetot , predmmtot , -1       , predeetoterr , predmmtoterr  , -1 );
    
  }

  printLine(latex);
  printOFRow( "tot pred" , predeetot+neefake , predmmtot+nmmfake , -1 , 
	      sqrt(pow(predeetoterr,2)+pow(0.5*neefake,2)) , sqrt(pow(predmmtoterr,2)+pow(0.5*nmmfake,2)) , -1 );

  printLine(latex);
  printOFRow( "Observed" , nee105 + nee1010 , nmm105+nmm1010 , nem105+nem1010 , 0 , 0 , 0 );

  printLine(latex);
  cout << endl << endl;


  float delta    = R10 * nee1010 + (1./R10) * nmm1010 + (1./R20) * nmm105 - nem1010 - nem105;
  float deltaerr = 0;
  deltaerr += pow( R10      * nee1010err , 2);
  deltaerr += pow( (1./R10) * nmm1010err , 2);
  deltaerr += pow(            nem1010err , 2);
  deltaerr += pow( (1./R20) * nmm105err  , 2);
  deltaerr += pow(            nem105err  , 2);
  deltaerr = sqrt(deltaerr);

  float R10up   = 1.05 * R10;
  float R20up   = 1.05 * R20;
  float deltaup = R10up * nee1010 + (1./R10up) * nmm1010 + (1./R20up) * nmm105 - nem1010 - nem105;

  float R10dn   = 0.95 * R10;
  float R20dn   = 0.95 * R20;
  float deltadn = R10dn * nee1010 + (1./R10dn) * nmm1010 + (1./R20dn) * nmm105 - nem1010 - nem105;

  float d1        = fabs( delta - deltaup );
  float d2        = fabs( delta - deltadn );
  float deltasyst = 0.5 * ( d1 + d2 );

  cout << endl;
  cout << Form("%.1f%s%.1f",nee1010,pm,nee1010err) << " & ";
  cout << Form("%.1f%s%.1f",nmm1010,pm,nmm1010err) << " & ";
  cout << Form("%.1f%s%.1f",nem1010,pm,nem1010err) << " & ";
  cout << Form("%.1f%s%.1f",nmm105,pm,nmm105err)   << " & ";
  cout << Form("%.1f%s%.1f",nem105,pm,nem105err)   << " & ";
  cout << Form("%.1f%s%.1f (stat)%s%.1f (syst)",delta,pm,deltaerr,pm,deltasyst) << endl;

  cout << Form("Delta  : %.1f%s%.1f (stat)%s%.1f (syst)",delta,pm,deltaerr,pm,deltasyst) << endl;
  cout << endl;


}

void makeStandardPlots( char* path , bool sigregion = false ){

  bool residual = false;
  bool log      = false;

  cout << "Plot residual? " << residual << endl;
  cout << "Do log plot?   " << log      << endl;

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  if( sigregion ){
    TCut highmet = "pfmet>275 && htpf>300";
    TCut highht = "pfmet>200 && htpf>600";

    sel = sel + ( highmet || highht );
 
    //sel = sel + "y > 8.5 && htpf > 300";
    //sel = sel + "y>13 && htpf > 300";
    //sel = sel + "pfmet>200 && htpf > 300";
    cout << "Signal region: " << sel.GetTitle() << endl;
  }

  char* filename;
  if(  highpt &&  sigregion ) filename = "datamc_highpt_sig";
  if(  highpt && !sigregion ) filename = "datamc_highpt";
  if( !highpt &&  sigregion ) filename = "datamc_lowpt_sig";
  if( !highpt && !sigregion ) filename = "datamc_lowpt";

  vector<char*> vars;
  vector<char*> xt;
  vector<int>   n;
  vector<float> xi;
  vector<float> xf;

  if( highpt ){
    if( !sigregion ){
      vars.push_back("lep1.pt()");  xt.push_back("max lepton p_{T} (GeV)");	n.push_back(10); xi.push_back(0.);   xf.push_back(200.);
      vars.push_back("lep2.pt()");  xt.push_back("min lepton p_{T} (GeV)");	n.push_back(10); xi.push_back(0.);   xf.push_back(100.);
      vars.push_back("lep1.eta()"); xt.push_back("max lepton #eta")       ;	n.push_back(10); xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("lep2.eta()"); xt.push_back("min lepton #eta");	        n.push_back(10); xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	n.push_back(10); xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	        n.push_back(10); xi.push_back(-3);   xf.push_back( 3);
      vars.push_back("dphijm");     xt.push_back("#Delta#phi(max jet,pfmet)");	n.push_back(10); xi.push_back(0);    xf.push_back(3.2);
      vars.push_back("dildphi");    xt.push_back("#Delta#phi(l_{1},l_{2})");	n.push_back(10); xi.push_back(0);    xf.push_back(3.2);
      vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(250.);
      vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(250.);
      vars.push_back("y");          xt.push_back("y (GeV^{1/2})");		n.push_back(10); xi.push_back(0.);   xf.push_back(20.);
      vars.push_back("htpf");       xt.push_back("H_{T} (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(1000.);
      vars.push_back("dilmass");    xt.push_back("m(ll) (GeV)");		n.push_back(12); xi.push_back(1.);   xf.push_back(301.);
      //vars.push_back("dilmass");    xt.push_back("m(ll) (GeV)");		n.push_back(60); xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("dilpt");      xt.push_back("p_{T}(ll) (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("npfjets");    xt.push_back("npfjets");			n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
      vars.push_back("nbtags");     xt.push_back("nbtags");			n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
      vars.push_back("ndavtx");     xt.push_back("nDAVertices");		n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
      vars.push_back("mt2jcore");   xt.push_back("MT2J");			n.push_back(10); xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("mt2");        xt.push_back("MT2");			n.push_back(10); xi.push_back(0.);   xf.push_back(150.);
      vars.push_back("meff");       xt.push_back("meff");			n.push_back(10); xi.push_back(0.);   xf.push_back(2000.);
    }else{
      vars.push_back("leptype");    xt.push_back("dilepton flavor");	        n.push_back(3);  xi.push_back(0.);   xf.push_back(3);
      vars.push_back("lep1.pt()");  xt.push_back("max lepton p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(250.);
      vars.push_back("lep2.pt()");  xt.push_back("min lepton p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(150.);
      vars.push_back("lep1.eta()"); xt.push_back("max lepton #eta")       ;	n.push_back(5);  xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("lep2.eta()"); xt.push_back("min lepton #eta");	        n.push_back(5);  xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	        n.push_back(5);  xi.push_back(-3);   xf.push_back( 3);
      vars.push_back("dphijm");     xt.push_back("#Delta#phi(max jet,pfmet)");	n.push_back(10); xi.push_back(0.);   xf.push_back(3.2);
      vars.push_back("dildphi");    xt.push_back("#Delta#phi(l_{1},l_{2})");	n.push_back(10); xi.push_back(0.);   xf.push_back(3.2);
      vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("y");          xt.push_back("y (GeV^{1/2})");		n.push_back(5);  xi.push_back(8.5);  xf.push_back(28.5);
      vars.push_back("htpf");       xt.push_back("H_{T} (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(1000.);
      vars.push_back("dilmass");    xt.push_back("m(ll) (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("dilpt");      xt.push_back("p_{T}(ll) (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("npfjets");    xt.push_back("npfjets");			n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
      vars.push_back("nbtags");     xt.push_back("nbtags");			n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
      vars.push_back("ndavtx");     xt.push_back("nDAVertices");		n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
      vars.push_back("mt2jcore");   xt.push_back("MT2J");			n.push_back(5);  xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("mt2");        xt.push_back("MT2");			n.push_back(5);  xi.push_back(0.);   xf.push_back(200.);
      vars.push_back("meff");       xt.push_back("meff");			n.push_back(5);  xi.push_back(0.);   xf.push_back(2000.);
    }
  }else{
      vars.push_back("lep1.pt()");  xt.push_back("max lepton p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(50.);
      vars.push_back("lep2.pt()");  xt.push_back("min lepton p_{T} (GeV)");	n.push_back(4);  xi.push_back(0.);   xf.push_back(20.);
      vars.push_back("lep1.eta()"); xt.push_back("max lepton #eta")       ;	n.push_back(5);  xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("lep2.eta()"); xt.push_back("min lepton #eta");	        n.push_back(5);  xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	        n.push_back(5);  xi.push_back(-3);   xf.push_back( 3);
      vars.push_back("dphijm");     xt.push_back("#phi(max jet,pfmet)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(3.2);
      vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("y");          xt.push_back("y (GeV^{1/2})");		n.push_back(5);  xi.push_back(0);    xf.push_back(20.);
      vars.push_back("htpf");       xt.push_back("H_{T} (GeV)");		n.push_back(7);  xi.push_back(0.);   xf.push_back(700.);
      vars.push_back("dilmass");    xt.push_back("m(ll) (GeV)");		n.push_back(5);  xi.push_back(1.);   xf.push_back(150.);
      vars.push_back("dilpt");      xt.push_back("p_{T}(ll) (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(150.);
      vars.push_back("npfjets");    xt.push_back("npfjets");			n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
      vars.push_back("nbtags");     xt.push_back("nbtags");			n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
      vars.push_back("ndavtx");     xt.push_back("nDAVertices");		n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
      vars.push_back("mt2jcore");   xt.push_back("MT2J");			n.push_back(5);  xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("mt2");        xt.push_back("MT2");			n.push_back(5);  xi.push_back(0.);   xf.push_back(100.);
      vars.push_back("meff");       xt.push_back("meff");			n.push_back(5);  xi.push_back(0.);   xf.push_back(1000.);
  } 
  
  const unsigned int nvars = vars.size();
  
  //TCanvas *can[nvars];
  TPad* legpad[nvars];
  TPad* plotpad[nvars];

  //int canCounter = -1;

  TCanvas* canvas = new TCanvas("canvas","canvas",1100,750);
  gStyle->SetPaperSize(22,28);
  canvas->Print(Form("../plots/%s.ps[",filename));

  TLegend *leg = getLegend( mc , mclabels , true , 0.1 , 0.3 , 0.8 , 0.7 );
  leg->SetTextSize(0.1);
  leg->SetBorderSize(1);
  
  for( unsigned int ivar = 0 ; ivar < nvars ; ++ivar ){     

    //can[ivar] = new TCanvas(Form("%s_can",vars[ivar]),Form("%s_can",vars[ivar]),1000,800);
    canvas->cd();

    legpad[ivar] = new TPad("legpad","legpad",0.8,0,1,1);
    legpad[ivar]->Draw();
    legpad[ivar]->cd();

    leg->Draw();

    //can[ivar]->cd();
    canvas->cd();

    plotpad[ivar] = new TPad("plotpad","plotpad",0,0,0.8,1);
    plotpad[ivar]->Draw();
    plotpad[ivar]->cd();
    
    plotpad[ivar]->Divide(2,2);

    plotpad[ivar]->cd(1);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==0") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "ee"  );
    plotpad[ivar]->cd(2);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==1") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "mm"  );
    plotpad[ivar]->cd(3);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==2") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "em"  );
    plotpad[ivar]->cd(4);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel)              , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , log , "all" );

    canvas->Print(Form("../plots/%s.ps",filename));
    canvas->Clear();

    //if( printgif ) can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
  } 

  canvas->Print(Form("../plots/%s.ps]",filename));
  canvas->Clear();
  
  gROOT->ProcessLine(Form(".! ps2pdf ../plots/%s.ps ../plots/%s.pdf",filename,filename));
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
  compareDataMC( mc , mclabels , data , "dilmass" , TCut(ofsel+sf) , weight , 10 , 0 , 200 , "m(ll) (GeV)" , true );

  can_of->cd(2);
  compareDataMC( mc , mclabels , data , "dilmass" , TCut(ofsel+of) , weight , 10 , 0 , 200 , "m(ll) (GeV)" , true );

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

  if( printgif ){
    can_of->Print("../plots/of.gif");
    can_of2->Print("../plots/ofsub.gif");
  }
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





void plotVictoryHist( TH1F* hpred , TH1F* hobs , TH1F* hpred_data, TH1F* hobs_data,
		      string title, string xtitle , int rebin , float metcut, char* htrange, float ymin, float ymax){

  bool plotData = true;

  hpred->SetLineColor(4);
  hobs ->SetLineColor(2);
  hpred->SetMarkerColor(4);
  hobs ->SetMarkerColor(2);
  hpred->SetMarkerStyle(23);
  hobs->SetMarkerStyle(4);
  hpred_data->SetLineColor(4);
  hobs_data ->SetLineColor(2);
  hpred_data->SetMarkerColor(4);
  hobs_data ->SetMarkerColor(2);
  hpred_data->SetMarkerStyle(23);
  hobs_data->SetMarkerStyle(4);
  hpred_data->SetMarkerSize(2);
  hobs_data->SetMarkerSize(2);
  hpred->SetFillColor(0);
  hobs ->SetFillColor(0);
  hpred->SetTitle( title.c_str() );
  hpred->GetXaxis()->SetTitle( xtitle.c_str() );
  hpred->GetXaxis()->SetLabelSize(0.04);
  hpred->GetXaxis()->SetTitleSize(0.055);
  hpred->GetYaxis()->SetTitle( "Events   " );
  hpred->GetYaxis()->SetTitleSize(0.055);
  hpred->GetYaxis()->SetTitleOffset(1);
  //hpred->GetYaxis()->SetTitle( "Entries / 2.125 GeV" );

  hpred->SetMinimum(ymin);
  hpred->SetMaximum(ymax);

  int metbin = hobs->FindBin( metcut );
  float pred = hpred->Integral( metbin , 100000000);
  float obs  = hobs->Integral(  metbin , 100000000);

  //float pred_data = hpred_data->Integral( metbin , 100000000);
  //float obs_data  = hobs_data->Integral(  metbin , 100000000);

  stringstream spred;
  spred << "Pred: " << pred ;
  stringstream sobs;
  sobs << "Obs:  " << obs ;

  if( rebin > 1 ){
    hobs->Rebin ( rebin );
    hpred->Rebin( rebin );
    hobs_data->Rebin ( rebin );
    hpred_data->Rebin( rebin );
  }
 
  if( plotData ){
    hpred->Draw("hist");
    hobs->Draw("samehist");
    hpred_data->Draw("sameE1");
    hobs_data->Draw("sameE1");
    hpred->SetLineWidth(2);
    hobs->SetLineWidth(2);
  }else{
    hpred->Draw("E1");
    hobs->Draw("sameE1");
  } 
  gPad->SetLogy(1);

  //TLatex t;
  //t.SetNDC();
  //t.DrawLatex(0.55,0.65,Form("obs/pred = %.2f",obs/pred));

  TLegend *leg = new TLegend(0.74,0.72,0.98,0.95);
  //leg->AddEntry(hpred, Form("pred>%.1f GeV = %.3f",metcut,pred) ,"l");
  //leg->AddEntry(hobs,  Form("obs>%.1f GeV = %.3f",metcut,obs)  ,"l");
  leg->AddEntry(hpred, "MC predicted" ,"l");
  leg->AddEntry(hobs,  "MC observed"  ,"l");
  if( plotData ){
    leg->AddEntry(hpred_data, "data predicted");
    leg->AddEntry(hobs_data,  "data observed" );
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.DrawLine( metcut , ymin , metcut , ymax );
  //line.DrawLine( metcut , 0.5 * TMath::Min( hobs->GetMinimum() , hpred->GetMinimum() ),
  //               metcut , 2 *   TMath::Max( hobs->GetMaximum() , hpred->GetMaximum() ) );

  if( plotData ){
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.045);
    //t->DrawLatex(0.18,0.32,"CMS");
    //t->DrawLatex(0.18,0.27,"34 pb^{-1}  at  #sqrt{s} = 7 TeV");
    text->DrawLatex(0.20,0.88,"Events with ee/#mu#mu/e#mu");
    text->DrawLatex(0.20,0.82,htrange);

    // if( issignal ){
    //   t->DrawLatex(0.20,0.82,"H_{T} > 300 GeV");
    //   //t->SetTextSize(0.05);
    //   //t->DrawLatex(0.8,0.6,"(a)");
    // }
    // else{
    //   t->DrawLatex(0.20,0.82,"125 < H_{T} < 300 GeV");
    //   //t->SetTextSize(0.05);
    //   //t->DrawLatex(0.8,0.6,"(b)");
    // }

 }
}




void printVictoryRow( string sample , TH1F* hobs , TH1F* hpred , float cut , float mcratio ){

  stringstream sobs;
  stringstream spred;
  stringstream sratio;

  int   minbin     = hobs->FindBin( cut );
  int   maxbin     = 10000;

  float pred    = hpred->Integral( minbin , maxbin);
  float prederr = calculateHistError( hpred , minbin , maxbin );

  if( sample == "data" ){
    pred    *= mcratio;
    prederr *= mcratio;
    //cout << "Scaling data " << mcratio << endl;
  }


  float obs     = hobs->Integral(  minbin , maxbin );
  float obserr  = calculateHistError( hobs , minbin , maxbin );

  float ratio    = pred > 0 ? obs / pred : 0;
  float ratioerr = obs > 0 && pred > 0 ? ratio * sqrt( pow( prederr / pred , 2 ) + pow( obserr / obs , 2 ) ) : 0;
  
  sobs   << Form( "%.2f" , obs   ) << pm << Form( "%.2f" , obserr   ); 
  spred  << Form( "%.2f" , pred  ) << pm << Form( "%.2f" , prederr  ); 
  sratio << Form( "%.2f" , ratio ) << pm << Form( "%.2f" , ratioerr ); 

  cout << delimstart << setw(width1) << sample        << setw(width2)
       << delim      << setw(width1) << spred.str()   << setw(width2)
       << delim      << setw(width1) << sobs.str()    << setw(width2)
       << delim      << setw(width1) << sratio.str()  << setw(width2)
       << delimend   << endl;

}


pair<float,float> getMCratio( TH1F* hobs , TH1F* hpred , float cut ){

  int   minbin     = hobs->FindBin( cut );
  int   maxbin     = 10000;

  float pred    = hpred->Integral( minbin , maxbin);
  float prederr = calculateHistError( hpred , minbin , maxbin );

  float obs     = hobs->Integral(  minbin , maxbin );
  float obserr  = calculateHistError( hobs , minbin , maxbin );

  float ratio   = pred > 0 ? obs / pred : 0;
  float ratioerr = obs > 0 && pred > 0 ? ratio * sqrt( pow( prederr / pred , 2 ) + pow( obserr / obs , 2 ) ) : 0;

  return make_pair( ratio , ratioerr );
}
