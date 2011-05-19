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
    data2010->Reset();
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
    t->Reset();
    wjets->Reset();
    LM0->Reset();
    LM1->Reset();
    LM2->Reset();
    LM3->Reset();
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
    t		= new TChain("t");
    wjets	= new TChain("t");
    LM0 	= new TChain("t");
    LM1 	= new TChain("t");
    LM2 	= new TChain("t");
    LM3 	= new TChain("t");
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
  data2010->Add("/tas03/home/benhoob/OSSusy/output_38X/nov5th_v6_skim/dataskim_smallTree.root");
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
  dytautau->Add(Form("%s/DYtot_smallTree.root",path));
  ww->Add(Form("%s/ww_smallTree.root",path));
  wz->Add(Form("%s/wz_smallTree.root",path));
  zz->Add(Form("%s/zz_smallTree.root",path));
  t->Add(Form("%s/tW_smallTree.root",path));
  wjets->Add(Form("%s/wjetsMG_smallTree.root",path));
  //wjets->Add(Form("%s/wjets_smallTree.root",path));
  LM0->Add(Form("%s/LM0_smallTree.root",path));
  LM1->Add(Form("%s/LM1_smallTree.root",path));
  LM2->Add(Form("%s/LM2_smallTree.root",path));
  LM3->Add(Form("%s/LM3_smallTree.root",path));
  other->Add(Form("%s/ww_smallTree.root",path));
  other->Add(Form("%s/wz_smallTree.root",path));
  other->Add(Form("%s/zz_smallTree.root",path));
  other->Add(Form("%s/tW_smallTree.root",path));


  //mc.push_back(ttall);     mclabels.push_back("ttall");
  //mc.push_back(ttdil);     mclabels.push_back("ttdil");
  //mc.push_back(ttotr);     mclabels.push_back("ttotr");
  //mc.push_back(ttpowheg);  mclabels.push_back("ttpowheg");
  //mc.push_back(ttotr);    mclabels.push_back("ttotr");
  mc.push_back(ttll);     mclabels.push_back("ttll");    mctex.push_back("$t\\bar{b}\\rightarrow\\ell^+\\ell^-$");
  mc.push_back(tttau);    mclabels.push_back("tttau");   mctex.push_back("$t\\bar{b}\\rightarrow\\ell^{\\pm}\\tau^{\\mp}$");
  mc.push_back(ttfake);   mclabels.push_back("ttfake");  mctex.push_back("$t\\bar{b}\\rightarrow$fake");
  mc.push_back(wjets);    mclabels.push_back("wjets");       mctex.push_back("$W^{\\pm}$+jets");
  //mc.push_back(zjets);    mclabels.push_back("zjets");     mctex.push_back("$Z^0$+jets");
  //mc.push_back(dydata);   mclabels.push_back("DYdata");      mctex.push_back("DYdata");
  mc.push_back(dy);       mclabels.push_back("DY");          mctex.push_back("DY");
  //mc.push_back(dytautau); mclabels.push_back("DYtautau");    mctex.push_back("DYtautau");
  mc.push_back(ww);       mclabels.push_back("WW");          mctex.push_back("W^+W^-");
  mc.push_back(wz);       mclabels.push_back("WZ");          mctex.push_back("W^{\\pm}Z^0");
  mc.push_back(zz);       mclabels.push_back("ZZ");          mctex.push_back("Z^0Z^0");
  mc.push_back(t);        mclabels.push_back("t");           mctex.push_back("single top");
  //mc.push_back(other);    mclabels.push_back("other");   mctex.push_back("WW/WZ/ZZ/t");
  mc.push_back(LM0);      mclabels.push_back("LM0");     mctex.push_back("LM0");
  mc.push_back(LM1);      mclabels.push_back("LM1");     mctex.push_back("LM1");
  mc.push_back(LM2);      mclabels.push_back("LM2");     mctex.push_back("LM2");
  mc.push_back(LM3);      mclabels.push_back("LM3");     mctex.push_back("LM3");
  
  sm->Add(Form("%s/ttall_smallTree.root",path));
  sm->Add(Form("%s/wjetsMG_smallTree.root",path));
  //sm->Add(Form("%s/ttdil_smallTree.root",path));
  //sm->Add(Form("%s/ttotr_smallTree.root",path));
  sm->Add(Form("%s/DYtot_smallTree.root",path));
  sm->Add(Form("%s/ww_smallTree.root",path));
  sm->Add(Form("%s/wz_smallTree.root",path));
  sm->Add(Form("%s/zz_smallTree.root",path));
  sm->Add(Form("%s/tW_smallTree.root",path));


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
  TCut met50("pfmet > 50");
  TCut met30("pfmet > 30");
  TCut ht100("htpf > 100");
  TCut ht200("htpf > 200");
  TCut ht250("htpf > 250");
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

  TCut highptsel = zveto + njets2 + met50 + ht100 + pt2010;
  TCut lowptsel  = zveto + njets2 + met50 + ht250 + pt105 + !pt2010 + lep2tight;

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

  //TCut weight("weight * ndavtxweight * (5.53/191.)");
  TCut weight("weight * ndavtxweight");
  //TCut weight("weight*(186./43.)*1.17*ndavtxweight");
  //TCut weight("weight");
  //TCut weight("1");

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
  TCut weight("weight*ndavtxweight");


  TH1F* hee = new TH1F("hee","",1,0,1);
  TH1F* hmm = new TH1F("hmm","",1,0,1);
  hee->Sumw2();
  hmm->Sumw2();

  TCanvas *ctemp = new TCanvas();
  sm->Draw("0.5>>hee",(eetype+zmass)*weight);
  sm->Draw("0.5>>hmm",(mmtype+zmass)*weight);
  
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
  text->SetTextColor(4);
  text->SetNDC();

  TCanvas *zcan = new TCanvas("zcan","zcan",1200,600);
  zcan->Divide(2,1);

  zcan->cd(1);
  compareDataMC( mc , mclabels , data , "dilmass" , ee , weight , 100 , 0 , 200 , "M(e^{+}e^{-}) (GeV)"     , true, false, true, true, "ee" );
  text->DrawLatex(0.2,0.75,Form("N_{DATA} = %.0f",needata));
  text->DrawLatex(0.2,0.70,Form("N_{MC}   = %.0f",neeMC));
  text->DrawLatex(0.2,0.65,Form("N_{DATA}/N_{MC} = %.2f",needata/neeMC));
  

  zcan->cd(2);
  compareDataMC( mc , mclabels , data , "dilmass" , mm , weight , 100 , 0 , 200 , "M(#mu^{+}#mu^{-}) (GeV)" , true, false, true, true, "mm" );
  text->DrawLatex(0.2,0.75,Form("N_{DATA} = %.0f",nmmdata));
  text->DrawLatex(0.2,0.70,Form("N_{MC}   = %.0f",nmmMC));
  text->DrawLatex(0.2,0.65,Form("N_{DATA}/N_{MC} = %.2f",nmmdata/nmmMC));

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

  TCut sigsel = zveto + njets2 + met50 + ht100 + pt2010 + sig;

  cout << "Using selection: " << sigsel.GetTitle() << endl;
  char* colformat = "precision=3 col=7.6::10.10:";

  data->Scan("run:lumi:event:leptype:pfmet:npfjets:htpf:dilmass:y",sigsel,colformat);
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
  
  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  
  printYields( mc , mclabels , data , sel , weight , latex );
 
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
  text->DrawLatex(0.35,0.80,"191 pb^{-1} at #sqrt{s} = 7 TeV");
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
  //float x4=1500;
  
  float y1 = 4.5;
  float y2 = 8.5;
  float y3 = 8.5;
  //float y4 = 30;

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

  for( unsigned int i = 0 ; i < n ; ++i ){

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
    vars.push_back("pfmet");     xt.push_back("pfmet (GeV)");      n.push_back(20); xi.push_back(0.); xf.push_back(200.);
    //vars.push_back("y");         xt.push_back("y #equiv MET  /  #sqrt{H_{T}} (GeV^{1/2})");    n.push_back(20); xi.push_back(0.); xf.push_back(20.);
    //vars.push_back("htpf");      xt.push_back("H_{T} (GeV)");      n.push_back(20); xi.push_back(0.); xf.push_back(1000.);
    //vars.push_back("dilmass");   xt.push_back("M(ll) (GeV)");      n.push_back(60); xi.push_back(1.); xf.push_back(301.);
    //vars.push_back("lep1.eta()");   xt.push_back("#eta(lep1)");      n.push_back(50); xi.push_back(-3.); xf.push_back(3.);
    //vars.push_back("lep2.eta()");   xt.push_back("#eta(lep2)");      n.push_back(50); xi.push_back(-3.); xf.push_back(3.);
    //vars.push_back("lep1.pt()");   xt.push_back("pt(lep1)");         n.push_back(50); xi.push_back(0.); xf.push_back(100.);
    //vars.push_back("lep2.pt()");   xt.push_back("pt(lep2)");         n.push_back(50); xi.push_back(0.); xf.push_back(100.);
    //vars.push_back("dilpt");     xt.push_back("p_{T}(ll) (GeV)");  n.push_back(20); xi.push_back(0.); xf.push_back(300.);
    //vars.push_back("npfjets");   xt.push_back("njets");            n.push_back(10); xi.push_back(0.); xf.push_back(10.);
    //vars.push_back("ndavtx");      xt.push_back("nDAVertices");      n.push_back(20); xi.push_back(0.); xf.push_back(20.);
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

  bool residual = true;
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

    compareDataMC( mc , mclabels , data , vars[ivar] , sel , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , !combine4 );

    if( printgif ) can[ivar]->Print(Form("../plots/%s.png",vars[ivar]));
  } 
}



void makeStandardPlots( char* path , bool sigregion = false ){

  bool residual = true;

  deleteHistos();

  bool highpt = false;
  if( TString(path).Contains("highpt") ) highpt = true;

  initialize(path);

  TCut sel    = selection_TCut(highpt);
  TCut weight = weight_TCut();

  if( sigregion ){
    sel = sel + "y > 8.5 && htpf > 300";
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
      vars.push_back("lep1.pt()");  xt.push_back("max lepton p_{T} (GeV)");	n.push_back(20); xi.push_back(0.);   xf.push_back(200.);
      vars.push_back("lep2.pt()");  xt.push_back("min lepton p_{T} (GeV)");	n.push_back(20); xi.push_back(0.);   xf.push_back(100.);
      //vars.push_back("lep1.eta()"); xt.push_back("max lepton #eta")       ;	n.push_back(10); xi.push_back(-3.);  xf.push_back(3.);
      //vars.push_back("lep2.eta()"); xt.push_back("min lepton #eta");	        n.push_back(10); xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	n.push_back(10); xi.push_back(0.);   xf.push_back(300.);
      //vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	        n.push_back(10); xi.push_back(-3);   xf.push_back( 3);
      //vars.push_back("dphijm");     xt.push_back("#phi(max jet,pfmet)");	n.push_back(10); xi.push_back(0.);   xf.push_back(3.2);
      vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");		n.push_back(25); xi.push_back(0.);   xf.push_back(250.);
      vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");		n.push_back(25); xi.push_back(0.);   xf.push_back(250.);
      vars.push_back("y");          xt.push_back("y (GeV^{1/2})");		n.push_back(10); xi.push_back(0.);   xf.push_back(20.);
      vars.push_back("htpf");       xt.push_back("H_{T} (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(1000.);
      vars.push_back("dilmass");    xt.push_back("M(ll) (GeV)");		n.push_back(12); xi.push_back(1.);   xf.push_back(301.);
      vars.push_back("dilmass");    xt.push_back("M(ll) (GeV)");		n.push_back(60); xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("dilpt");      xt.push_back("p_{T}(ll) (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(300.);
      vars.push_back("npfjets");    xt.push_back("npfjets");			n.push_back(10); xi.push_back(0.);   xf.push_back(10.);
      vars.push_back("nbtags");     xt.push_back("nbtags");			n.push_back(5);  xi.push_back(0.);   xf.push_back(5.);
      vars.push_back("ndavtx");     xt.push_back("nDAVertices");		n.push_back(20); xi.push_back(0.);   xf.push_back(20.);
      //vars.push_back("mt2jcore");   xt.push_back("MT2J");			n.push_back(10); xi.push_back(0.);   xf.push_back(500.);
      //vars.push_back("mt2");        xt.push_back("MT2");			n.push_back(10); xi.push_back(0.);   xf.push_back(150.);
      //vars.push_back("meff");       xt.push_back("meff");			n.push_back(10); xi.push_back(0.);   xf.push_back(2000.);
    }else{
      vars.push_back("lep1.pt()");  xt.push_back("max lepton p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(250.);
      vars.push_back("lep2.pt()");  xt.push_back("min lepton p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(150.);
      vars.push_back("lep1.eta()"); xt.push_back("max lepton #eta")       ;	n.push_back(5);  xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("lep2.eta()"); xt.push_back("min lepton #eta");	        n.push_back(5);  xi.push_back(-3.);  xf.push_back(3.);
      vars.push_back("jet.pt()");   xt.push_back("max jet p_{T} (GeV)");	n.push_back(5);  xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("jet.eta()");  xt.push_back("max jet #eta");	        n.push_back(5);  xi.push_back(-3);   xf.push_back( 3);
      vars.push_back("dphijm");     xt.push_back("#phi(max jet,pfmet)");	n.push_back(10); xi.push_back(0.);   xf.push_back(3.2);
      vars.push_back("tcmet");      xt.push_back("tcmet (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("pfmet");      xt.push_back("pfmet (GeV)");		n.push_back(10); xi.push_back(0.);   xf.push_back(500.);
      vars.push_back("y");          xt.push_back("y (GeV^{1/2})");		n.push_back(5);  xi.push_back(8.5);  xf.push_back(28.5);
      vars.push_back("htpf");       xt.push_back("H_{T} (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(1000.);
      vars.push_back("dilmass");    xt.push_back("M(ll) (GeV)");		n.push_back(5);  xi.push_back(0.);   xf.push_back(300.);
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
      vars.push_back("dilmass");    xt.push_back("M(ll) (GeV)");		n.push_back(5);  xi.push_back(1.);   xf.push_back(150.);
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

  TLegend *leg = getLegend( mc , mclabels , true , 0.2 , 0.3 , 0.6 , 0.7 );
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
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==0") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , "ee"  );
    plotpad[ivar]->cd(2);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==1") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , "mm"  );
    plotpad[ivar]->cd(3);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel+"leptype==2") , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , "em"  );
    plotpad[ivar]->cd(4);
    compareDataMC( mc , mclabels , data , vars[ivar] , TCut(sel)              , weight , n[ivar] , xi[ivar] , xf[ivar] , xt[ivar] , true , residual , false , "all" );

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




