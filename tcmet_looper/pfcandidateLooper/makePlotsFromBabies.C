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
#include <sstream>


bool printgif_ = false;
using namespace std;

void printChain(TChain* ch, TCut sel);
void plotHist(TCanvas *can ,TH1F *hA, TH1F* hB, char* title, char* var, char* leg1, char* leg2 , bool norm = false, bool text = false );
void effVsPU( char* sample, TChain *ch , vector<char*> metTypes , TCut mysel, float metcut , int rebin = 1);

void effVsPU( char* title, vector<TChain*> ch , vector<char*> metTypes , vector<char*> labels, 
              TCut mysel , float metcut , int rebin );

void ROC(TChain *chsig, TChain *chbkg, vector<char*> metTypes , 
	 TCut sel, vector<char*> labels );


void compareDataMC( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , char* var , 
		    TCut sel , int nbins ,  float xmin , float xmax , float metcut , char* xtitle = ""){

  int colors[]={2,5,7};

  assert( chmc.size() == labels.size() );
  const unsigned int nmc = chmc.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist;
  TH1F*    mchist[nmc];
  TH1F*    datahist = new TH1F("datahist","datahist",nbins,xmin,xmax);

  int bin81  = datahist->FindBin(81);
  int bin101 = datahist->FindBin(101)-1;

  //TCut weight("weight");
  TCut weight("weight*davtxweight"); 

  TLegend *leg = new TLegend(0.65,0.8,0.95,0.95);
  //TLegend *leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->AddEntry(datahist,"data","p");
  
  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){

    //char* weightchar =         "0.02 * (157.5/660.)";
    //if( imc > 0 ) weightchar = "0.02 * (1666./2000.)";
    //TCut weight(weightchar);

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),nbins,xmin,xmax);

    chmc.at(imc)->Draw(Form("TMath::Min(%s,%f)>>mc_%i",var,xmax-0.01,imc),sel*weight);

    mchist[imc]->SetFillColor( colors[imc] );

    mcstack->Add( mchist[imc] );

    if( imc == 0 ) mctothist = (TH1F*) mchist[imc]->Clone();
    else           mctothist->Add(mchist[imc]);

    leg->AddEntry(mchist[imc],labels.at(imc),"f");
    //cout << labels.at(imc) << " " << mchist[imc]->Integral(bin81,bin101) << endl;
    cout << labels.at(imc) << " " << mchist[imc]->Integral() << endl;
  }

  chdata->Draw(Form("TMath::Min(%s,%f)>>datahist",var,xmax-0.01),sel);

  //cout << "data " << datahist->Integral(bin81,bin101) << endl;
  cout << "data " << datahist->Integral() << endl;

  int mybin        = datahist->FindBin(metcut);
  float intdata    = datahist->Integral(mybin,10000);
  float intmc      = mctothist->Integral(mybin,10000);
  float totintdata = datahist->Integral();
  float totintmc   = mctothist->Integral();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextColor(4);

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(4);
  
  datahist->GetXaxis()->SetTitle(xtitle);

  float min = 0.05;
  float max = 2 * datahist->GetMaximum();

  datahist->Draw("E1");
  datahist->SetMinimum(min);
  datahist->SetMaximum(max);
  mcstack->Draw("same");
  datahist->Draw("sameE1");
  datahist->Draw("sameaxis");

  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  if( TString(var).Contains("met") ){

    t->SetTextColor(1);
    t->SetTextSize(0.035);
    t->DrawLatex(0.67,0.70,Form("N_{DATA} (tot) = %.0f",totintdata));
    t->DrawLatex(0.67,0.64,Form("N_{MC}   (tot) = %.0f",totintmc));
    t->SetTextSize(0.04);

    if( metcut > 0 ){
      t->DrawLatex(0.67,0.55,Form("N_{DATA} (>%.0f) = %.0f",metcut,intdata));
      t->DrawLatex(0.67,0.48,Form("N_{MC}   (>%.0f) = %.2f",metcut,intmc));
      
      line.SetLineColor(1);
      line.DrawLine(metcut,min,metcut,max);
    }

    // int bin60 = datahist->FindBin(60);
    // float intdata = datahist->Integral(bin60,10000);
    // float intmc   = mctothist->Integral(bin60,10000);
    
    // t->SetTextSize(0.04);
    // t->DrawLatex(0.7,0.60,Form("N_{DATA} >60 = %.0f",intdata));
    // t->DrawLatex(0.7,0.54,Form("N_{MC}   >60 = %.1f",intmc));
    // t->DrawLatex(0.7,0.450,Form("N_{DATA} >60 = %.0f",intdata));
    // t->DrawLatex(0.7,0.38,Form("N_{MC}   >60 = %.1f",intmc));
  }

  else if( TString(var).Contains("dilmass") ){
    int bin81  = datahist->FindBin(76);
    int bin101 = datahist->FindBin(106)-1;
   
    float intdata = datahist->Integral(bin81,bin101);
    float intmc   = mctothist->Integral(bin81,bin101);
    
    t->SetTextSize(0.04);
    t->DrawLatex(0.7,0.60,Form("N_{DATA}  = %.0f",intdata));
    t->DrawLatex(0.7,0.54,Form("N_{MC}    = %.1f",intmc));

    datahist->SetMinimum(min);
    datahist->SetMaximum(max);
    line.DrawLine(76,min,76,max);
    line.DrawLine(106,min,106,max);
    

  }


  gPad->SetLogy(1);


}

void makePlots( bool printgif = false ){

  printgif_ = printgif;

  gStyle->SetOptStat(0);

  
  char* path = "output/PVT/promptreco/dcsonly";

  TChain *dymm = new TChain("Events");
  dymm->Add(Form("%s/dymm_default_baby.root",path));
  
  TChain *dyee = new TChain("Events");
  dyee->Add(Form("%s/dyee_default_baby.root",path));
  
  TChain *zjets = new TChain("Events");
  zjets->Add(Form("%s/zjets_default_baby.root",path));

  TChain *tt = new TChain("Events");
  tt->Add(Form("%s/tt_default_baby.root",path));

  TChain *ww = new TChain("Events");
  ww->Add(Form("%s/ww_default_baby.root",path));
  
  TChain *h130 = new TChain("Events");
  h130->Add(Form("%s/h130_default_baby.root",path));

  TChain *chdata = new TChain("Events");
  chdata->Add(Form("%s/data_default_baby.root",path));

  TChain *dy = new TChain("Events");
  dy->Add(Form("%s/dymm_default_baby.root",path));
  dy->Add(Form("%s/dyee_default_baby.root",path));

  TCut zmass("dilep.mass()>76&&dilep.mass()<106");
  TCut jet1pv("njets30==1 && jetpv==1");
  TCut jet1("njets30==1");
  TCut jet0("njets30==0");
  TCut ee("leptype==3");
  TCut mm("leptype==0");
  TCut sf=ee||mm;
  TCut em("leptype==1||leptype==2");
  TCut dphizjet("acos(cos(dilep.phi()-jet.phi()))>2.");
  TCut lepveto("nlep==0");
  TCut muveto("softmu==0");
  TCut topveto("toptag==0");
  TCut weight("weight*davtxweight"); 

  TCut fullsel= zmass+jet1+sf+lepveto+topveto+muveto;

  


  /*
  char* path = "output/promptreco/dcsonly";
  //char* path = "output/PVT/promptreco";

  TChain *dymm = new TChain("Events");
  dymm->Add(Form("%s/dymm_default_baby.root",path));
  
  TChain *dyee = new TChain("Events");
  dyee->Add(Form("%s/dyee_default_baby.root",path));

  TChain *dy = new TChain("Events");
  dy->Add(Form("%s/dyee_default_baby.root",path));
  dy->Add(Form("%s/dymm_default_baby.root",path));
  
  // TChain *dyjets = new TChain("Events");
  // dyjets->Add(Form("%s/zjets_spring11_default_baby.root",path));
  
  TChain *ttbar = new TChain("Events");
  ttbar->Add(Form("%s/tt_default_baby.root",path));

  TChain *chdata = new TChain("Events");
  chdata->Add(Form("%s/data_default_baby.root",path));

  TCut dphizjet("acos(cos(dilep.phi()-jet.phi()))>2");
  TCut lepveto("nlep==0");
  TCut zmass("dilep.mass()>76&&dilep.mass()<106");
  TCut jet0("njets30==0");
  TCut jet1("njets30==1");
  TCut jet2("njets30>1");

  TCut ee("leptype==0");
  TCut mm("leptype==3");
  TCut em("leptype==2||leptype==1");

  TCut zee      = zmass+ee;
  TCut selem    = zmass+em;
  TCut sel      = zmass+(mm||ee);
  TCut sel0     = zmass+(mm||ee)+"njets30==0";
  TCut sel1     = zmass+(mm||ee)+"njets30==1";
  TCut sel1pv   = zmass+(mm||ee)+"njets30==1 && jetpv==1";
  TCut sel1nopv = zmass+(mm||ee)+"njets30==1 && jetpv==0";
  TCut sel2     = zmass+(mm||ee)+"njets30>1";
  TCut sel1_lepveto_dphi = sel1 + lepveto + dphizjet;
  TCut ee_jetveto = ee+jet0;
  TCut mm_jetveto = mm+jet0;
  */
  

  vector<TChain*> chmc;
  chmc.push_back(tt);  
  chmc.push_back(ww);  
  //chmc.push_back(dymm);
  //chmc.push_back(dyee);
  //chmc.push_back(dyjets);
  chmc.push_back(dy);
  
  vector<char*> labels;
  labels.push_back("t#bar{t}");
  labels.push_back("W^{+}W^{-}");
  //labels.push_back("DY#rightarrow#mu#mu");
  //labels.push_back("DY#rightarrowee");
  //labels.push_back("DY");
  labels.push_back("DY");



  //-------------------------------
  // dilepton mass
  //-------------------------------
  /*
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);
  
  c1->cd(1);
  compareDataMC( chmc , labels , chdata , "dilmass" , mm , 90 , 20 , 200 , 30 , "M_{ll} (GeV)" );

  c1->cd(2);
  compareDataMC( chmc , labels , chdata , "dilmass" , ee , 90 , 20 , 200 , 30 , "M_{ll} (GeV)" );

  c1->Print("plots/dilmass.gif");
    
  //-------------------------------
  // all jets
  //-------------------------------
  
  TCanvas *c2 = new TCanvas("c2","",1200,600);
  c2->Divide(2,1);
  
  c2->cd(2);
  compareDataMC( chmc , labels , chdata , "tcmet" , sel , 50 , 0 , 100 , 30 , "tcmet (GeV)" );
  
  c2->cd(1);
  compareDataMC( chmc , labels , chdata , "pfmet" , sel , 50 , 0 , 100 , 30 , "pfmet (GeV)" );
  
  if( printgif ) c2->Print("plots/met_datamc_allj.gif");
    
  //-------------------------------
  // 0 jets
  //-------------------------------
  
  TCanvas *c3 = new TCanvas("c3","",1200,600);
  c3->Divide(2,1);
  
  c3->cd(2);
  compareDataMC( chmc , labels , chdata , "tcmet" , sel0 , 50 , 0 , 100 , 30 , "tcmet (GeV)" );
  
  c3->cd(1);
  compareDataMC( chmc , labels , chdata , "pfmet" , sel0 , 50 , 0 , 100 , 30 , "pfmet (GeV)" );
  
  if( printgif ) c3->Print("plots/met_datamc_0j.gif");
  
  //-------------------------------
  // 1 jets
  //-------------------------------
  

  //TCut mysel = sel1 + lepveto;
  //TCut mysel = zmass+jet1+sf;
  TCut mysel = fullsel;

  TCanvas *c4 = new TCanvas("c4","",1200,600);
  c4->Divide(2,1);
  
  c4->cd(1);
  compareDataMC( chmc , labels , chdata , "tcmet" , mysel , 20 , 0 , 100 , 45 , "tcmet (GeV)" );
  
  c4->cd(2);
  compareDataMC( chmc , labels , chdata , "pfmet" , mysel , 20 , 0 , 100 , 45 , "pfmet (GeV)" );
  
  if( printgif ) c4->Print("plots/met_datamc_1j.gif");

  /*
  TCanvas *c4_pf = new TCanvas("c4","",600,600);
  c4_pf->cd();
  compareDataMC( chmc , labels , chdata , "pfmet" , mysel , 20 , 0 , 100 , 45 , "pfmet (GeV)" );
  */

  /*
  //-------------------------------
  // >=2 jets
  //-------------------------------
  
  TCanvas *c5 = new TCanvas("c5","",1200,600);
  c5->Divide(2,1);
  
  c5->cd(1);
  compareDataMC( chmc , labels , chdata , "tcmet" , sel2 , 20 , 0 , 100 , "tcmet (GeV)" );
  
  c5->cd(2);
  compareDataMC( chmc , labels , chdata , "pfmet" , sel2 , 20 , 0 , 100 , "pfmet (GeV)" );
  
  if( printgif ) c5->Print("../plots/met_datamc_2j.gif");
  
  //-------------------------------
  // e/m
  //-------------------------------
  
  TCanvas *c6 = new TCanvas("c6","",1200,600);
  c6->Divide(2,1);
  
  c6->cd(1);
  compareDataMC( chmc , labels , chdata , "tcmet" , selem , 20 , 0 , 100 , "tcmet (GeV)" );

  c6->cd(2);
  compareDataMC( chmc , labels , chdata , "pfmet" , selem , 20 , 0 , 100 , "pfmet (GeV)" );

  if( printgif ) c6->Print("../plots/met_datamc_em.gif");

  */


  //---------------------------------
  // jet-Z MET jet associated to vtx?
  //---------------------------------

  TCut fullsel_nopv = zmass+jet1+sf+lepveto+topveto+muveto+"jetpv==0";
  TCut fullsel_pv   = zmass+jet1+sf+lepveto+topveto+muveto+"jetpv==1";

  TCanvas *c7 = new TCanvas("c7","",1200,600);
  c7->Divide(2,1);
  
  c7->cd(1);
  compareDataMC( chmc , labels , chdata , "jetzmet" , fullsel_pv   , 50 , 0 , 100 , -1 , "jet-Z MET (GeV)" );

  c7->cd(2);
  compareDataMC( chmc , labels , chdata , "jetzmet" , fullsel_nopv , 50 , 0 , 100 , -1 , "jet-Z MET (GeV)" );

  if( printgif ) c7->Print("plots/jetzmet_datamc.gif");









  //TCut weight("weight");
  /*
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  
  int nbins  = 100;
  float xmin = 0.;
  float xmax = 100.;
  
  TH1F* hmc     = new TH1F("hmc",  "", nbins,xmin,xmax);
  TH1F* hdata   = new TH1F("hdata","", nbins,xmin,xmax);

  hmc->Sumw2();
  hdata->Sumw2();
 
  //char* met1 = "pfmet";
  //char* met1 = "hmetpf";
  char* met1 = "TMath::Min(pfmet,hmetpf)";

  chmc->Draw  (Form("TMath::Min(%s,99.9) >> hmc"  ,met1),(sel+type)*weight);
  chdata->Draw(Form("TMath::Min(%s,99.9) >> hdata",met1),(sel+type)       );
  
  delete ctemp;

  
  TCanvas *c1 = new TCanvas();
  c1->cd();
  plotHist( c1 , hmc , hdata , "Spring11 MC vs. data" , "tcmet (GeV)" , "MC" , "data" , false , true );
  if( printgif ) c1->Print("plots/minmet.gif");
  
  TCanvas *c2 = new TCanvas();
  c2->cd();
  
  vector<char*> metTypes;
  metTypes.push_back("pfmet");
  metTypes.push_back("pfmet");
  //metTypes.push_back("hmetpf");
  //metTypes.push_back("hmetpf");
  //metTypes.push_back("hmetpf4");
  //metTypes.push_back("hmetpf4");
  metTypes.push_back("TMath::Min(pfmet,hmetpf)");
  metTypes.push_back("TMath::Min(pfmet,hmetpf)");

  vector<char*> labels;
  labels.push_back("MC (pfmet)");
  labels.push_back("data (pfmet)");
  //labels.push_back("MC (trk-MET)");
  //labels.push_back("data (trk-MET)");
  //labels.push_back("MC (trk/neu-MET)");
  //labels.push_back("data (trk/neu-MET)");
  labels.push_back("MC (min-MET)");
  labels.push_back("data (min-MET)");


  vector<TChain*> ch;
  ch.push_back(chmc);
  ch.push_back(chdata);
  ch.push_back(chmc);
  ch.push_back(chdata);

  effVsPU( "Spring11 MC vs. data", ch , metTypes , labels , sel , 30. , 1 );
  if( printgif ) c2->Print("plots/minmet30.gif");

  TCanvas *c3 = new TCanvas();
  c3->cd();

  effVsPU( "Spring11 MC vs. data", ch , metTypes , labels , sel , 45. , 1 );
  if( printgif ) c3->Print("plots/minmet45.gif");
  */  



  /*
  vector<char*> rocMetTypes;
  rocMetTypes.push_back("pfmet");
  rocMetTypes.push_back("hmetpf");
  rocMetTypes.push_back("hmetpf4");
  rocMetTypes.push_back("TMath::Min(pfmet,hmetpf)");

  vector<char*> rocMetLabels;
  rocMetLabels.push_back("pfmet");
  rocMetLabels.push_back("trk-met");
  rocMetLabels.push_back("trk/neu-met");
  rocMetLabels.push_back("min-met");

  TCanvas *c4 = new TCanvas();
  c4->cd();

  TCut jetveto1("njets30==0&&nvtx<8&&dilep.mass()>12&&dilep.mass()<76");
  //TCut jetveto1("njets30==0&&nvtx<8");

  ROC(h130,dymm,rocMetTypes,jetveto1,rocMetLabels);

  c4->Print("plots/ROC1_m76.gif");
  //c4->Print("plots/ROC1.gif");
  
  TCanvas *c5 = new TCanvas();
  c5->cd();
  
  TCut jetveto2("njets30==0&&nvtx>7&&dilep.mass()>12&&dilep.mass()<76");
  //TCut jetveto2("njets30==0&&nvtx>7");
  
  ROC(h130,dymm,rocMetTypes,jetveto2,rocMetLabels);
  
  c5->Print("plots/ROC2_m76.gif");
  //c5->Print("plots/ROC2.gif");
  */



  
//   TCut sel("njets30==0");

//   TCanvas *ctemp = new TCanvas();
//   ctemp->cd();
  
//   int nbins  = 100;
//   float xmin = 0.;
//   float xmax = 100.;

//   TH1F* dymm_met1   = new TH1F("dymm_met1",  "",nbins,xmin,xmax);
//   TH1F* dyee_met1   = new TH1F("dyee_met1",  "",nbins,xmin,xmax);
//   TH1F* h130_met1   = new TH1F("h130_met1",  "",nbins,xmin,xmax);
  
//   TH1F* dymm_met2   = new TH1F("dymm_met2",  "",nbins,xmin,xmax);
//   TH1F* dyee_met2   = new TH1F("dyee_met2",  "",nbins,xmin,xmax);
//   TH1F* h130_met2   = new TH1F("h130_met2",  "",nbins,xmin,xmax);

//   dymm_met1->Sumw2();
//   dyee_met1->Sumw2();
//   h130_met1->Sumw2();
  
//   dymm_met2->Sumw2();
//   dyee_met2->Sumw2();
//   h130_met2->Sumw2();
  
//   char* met1 = "pfmet";
//   char* met2 = "pfmet";

//   TCut sel1("njets30==0 && nvtx<5");
//   TCut sel2("njets30==0 && nvtx>=5");

//   dymm->Draw(Form("TMath::Min(%s,199.9) >> dymm_met1",met1),sel1);
//   dyee->Draw(Form("TMath::Min(%s,199.9) >> dyee_met1",met1),sel1);
//   h130->Draw(Form("TMath::Min(%s,199.9) >> h130_met1",met1),sel1);
  
//   dymm->Draw(Form("TMath::Min(%s,199.9) >> dymm_met2",met2),sel2);
//   dyee->Draw(Form("TMath::Min(%s,199.9) >> dyee_met2",met2),sel2);
//   h130->Draw(Form("TMath::Min(%s,199.9) >> h130_met2",met2),sel2);
  
//   delete ctemp;
  
//   TCanvas *c1 = new TCanvas();
//   c1->cd();
//   //plotHist( c1 , dymm_met1 , dymm_met2 , "DYmm" , "MET" , "trk-MET" , "+ neutrals" );
//   plotHist( c1 , dymm_met1 , dymm_met2 , "DYmm" , "pfmet (GeV)" , "N_{PV} < 5" , "N_{PV} #geq 5" , true );
  
//   TCanvas *c2 = new TCanvas();
//   c2->cd();
//   plotHist( c2 , dyee_met1 , dyee_met2 , "DYee" , "MET" , "trk-MET" , "+ neutrals" , true );
  
//   TCanvas *c3 = new TCanvas();
//   c3->cd();
//   plotHist( c3 , h130_met1 , h130_met2 , "Higgs130" , "MET" , "trk-MET" , "+ neutrals" );

  //------------------------------------------
  // met flavors in 0-jet bin
  //------------------------------------------

//   TCanvas *effcan = new TCanvas();
//   effcan->cd();

//   vector<char*> metTypes1;
//   metTypes1.push_back("pfmet");
//   metTypes1.push_back("hmetpf");
//   metTypes1.push_back("hmetpf8");

//   TCut mysel1("njets30==0 && dilmass>76 && dilmass<106");

//   effVsPU( "DY#rightarrow#mu#mu" , dymm , metTypes1 , mysel1 , 45. );
//   effcan->Print("plots/dymm_met_vtx.gif");
//   effcan->Clear();

//   effVsPU( "DY#rightarrowee" , dyee , metTypes1 , mysel1 , 45. );
//   effcan->Print("plots/dyee_met_vtx.gif");
//   effcan->Clear();
  
//   effVsPU( "Higgs 130 GeV" , h130 , metTypes1 , mysel1 , 45. );
//   effcan->Print("plots/h130_met_vtx.gif");
//   effcan->Clear();

//  effVsPU( "Data" , data , metTypes1 , mysel1 , 45. );


//   TCanvas *effcan_mll = new TCanvas();
//   effcan_mll->cd();

//   TCut jetveto_dilmass("njets30==0&&dilmass<76&&dilmass>12");

//   effVsPU( "DY#rightarrow#mu#mu (12 < M_{ll} < 76 GeV)" , dymm , metTypes1 , jetveto_dilmass , 45. );
//   effcan_mll->Print("plots/dymm_met_vtx_m76.gif");
//   effcan_mll->Clear();

//   effVsPU( "DY#rightarrowee (12 < M_{ll} < 76 GeV)"     , dyee , metTypes1 , jetveto_dilmass , 45. );  
//   effcan_mll->Print("plots/dyee_met_vtx_m76.gif");
//   effcan_mll->Clear();

//   effVsPU( "Higgs 130 GeV (12 < M_{ll} < 76 GeV)"       , h130 , metTypes1 , jetveto_dilmass , 45. );
//   effcan_mll->Print("plots/h130_met_vtx_m76.gif");
//   effcan_mll->Clear();
  

  /*
  //------------------------------------------
  // met flavors in 1-jet bin
  //------------------------------------------

  TCanvas *effcan2 = new TCanvas();
  effcan2->cd();

  vector<char*> metTypes2;
  metTypes2.push_back("pfmet");
  //metTypes2.push_back("hmetpf4");
  metTypes2.push_back("jetzmet");
  metTypes2.push_back("jetzmet4");
  metTypes2.push_back("jetzmet8");

  TCut mysel2("njets30==1 && jetpv==1");

  effVsPU( "DY#rightarrow#mu#mu" , dymm , metTypes2 , mysel2 , 45. );
  //effVsPU( "DY#rightarrowee"     , dyee , metTypes2 , mysel2 , 45. , 5 );
  */
}

void ROC(TChain *chsig, TChain *chbkg, vector<char*> metTypes , 
	 TCut sel, vector<char*> labels ){
 
  
  int colors[4] = {1,2,4,8};
  int styles[4] = {2,20,22,25};

  const unsigned int npoints = 100;
  const unsigned int nMetTypes = metTypes.size();

  TGraphErrors *gr[nMetTypes];
  TH1F* hsig[nMetTypes];
  TH1F* hbkg[nMetTypes];
  TLegend *leg = new TLegend(0.5,0.2,0.7,0.5);

  float nsigtot  = chsig->GetEntries(sel);
  float nbkgtot  = chbkg->GetEntries(sel);

  for( unsigned int i = 0 ; i < nMetTypes ; ++i ){

    hsig[i] = new TH1F(Form("hsig_%i",i),Form("hsig_%i",i),npoints,0,npoints);
    hbkg[i] = new TH1F(Form("hbkg_%i",i),Form("hbkg_%i",i),npoints,0,npoints);
    
    chsig->Draw(Form("%s>>hsig_%i",metTypes.at(i),i) , sel);
    chbkg->Draw(Form("%s>>hbkg_%i",metTypes.at(i),i) , sel);

    float sig[npoints];
    float bkg[npoints];
    float sigerr[npoints];
    float bkgerr[npoints];

    bool found = false;

    for( int icut = 0 ; icut < npoints ; icut++ ){

      TCut metcut(Form("%s>%i",metTypes.at(i),icut));
      //cout << endl << Form("%s>%i",metTypes.at(i),icut) << endl;

      //float nsigpass = chsig->GetEntries(sel+metcut);
      //float nbkgpass = chbkg->GetEntries(sel+metcut);

      float nsigpass = hsig[i]->Integral(icut+1,10000);
      float nbkgpass = hbkg[i]->Integral(icut+1,10000);

      if( nsigpass < 0.1 || nbkgpass < 0.1 ){
	sig[icut]    = 0.;
	bkg[icut]    = 0.;
	sigerr[icut] = 0.;
	bkgerr[icut] = 0.;
      }else{
	sig[icut]    = nsigpass / nsigtot;
	bkg[icut]    = nbkgpass / nbkgtot;
	sigerr[icut] = sqrt( nsigpass ) / nsigtot;
	bkgerr[icut] = sqrt( nbkgpass ) / nbkgtot;
      }

      if( bkg[icut] < 0.001 &&!found ){
	cout << metTypes.at(i) << " sig " << sig[icut] << " bkg " << bkg[icut] << endl;
	found=true;
      }
      //cout << "sig " << nsigpass << " / " << nsigtot << " = " << sig[icut] << " +/- " << sigerr[icut] << endl;
      //cout << "bkg " << nbkgpass << " / " << nbkgtot << " = " << bkg[icut] << " +/- " << bkgerr[icut] << endl;
    }

    gr[i] = new TGraphErrors(npoints,bkg,sig,bkgerr,sigerr);
    gr[i]->SetLineColor(colors[i]);
    gr[i]->SetMarkerColor(colors[i]);
    gr[i]->SetMarkerStyle(styles[i]);

    if( i == 0 ){
      gr[i]->GetYaxis()->SetTitle("Sig Efficiency");
      gr[i]->GetXaxis()->SetTitle("Bkg Efficiency");
      gr[i]->SetTitle("");
      gr[i]->GetXaxis()->SetLimits(0.0001,1);
      gr[i]->GetYaxis()->SetLimits(0.0001,1);
   
      gr[i]->GetXaxis()->SetRangeUser(0.0001,1);
      gr[i]->GetYaxis()->SetRangeUser(0,1);
      gPad->SetLogx(1);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
    }else{

    }

    leg->AddEntry(gr[i],labels.at(i),"p");
  }

  gr[0]->Draw("AP");
  for( unsigned int i = 1 ; i < nMetTypes ; ++i ){
    gr[i]->Draw("sameP");
  }  

  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.DrawLine(0.001,0.001,0.001,1);

}


void effVsPU( char* title, vector<TChain*> ch , vector<char*> metTypes , vector<char*> labels, 
              TCut mysel , float metcut , int rebin ){

  int colors[4] = {4,2,1,8};
  int styles[4] = {23,25,20,28};
 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
  
  const unsigned int nMetTypes = metTypes.size();

  TH1F* hist_tot[nMetTypes];
  TH1F* hist_pass[nMetTypes];

  for( unsigned int i = 0 ; i < nMetTypes ; ++i ){
    
    hist_pass[i] = new TH1F(Form("hist_pass_%i",i), Form("hist_pass_%i",i), 15,0.5,15.5); 
    hist_tot[i]  = new TH1F(Form("hist_tot_%i",i),  Form("hist_tot_%i",i),  15,0.5,15.5);
    hist_pass[i]->Sumw2();
    hist_tot[i]->Sumw2();
    
    for( unsigned int ipu = 1 ; ipu < 16 ; ++ipu ){
      
      TCut npv   (Form("nvtx==%i" ,ipu                        ));
      TCut met45 (Form("%s>%f"    ,metTypes.at(i) , metcut ));
      TCut myweight = ("1");
      if( TString(labels.at(i)).Contains("MC") ) myweight="weight*vtxweight";
      //if( TString(labels.at(i)).Contains("MC") ) myweight="weight";
      
      float npass = ch.at(i)->GetEntries((mysel+npv+met45)*myweight);
      float ntot  = ch.at(i)->GetEntries((mysel+npv)*myweight);

      hist_pass[i]->SetBinContent( ipu , npass );
      hist_tot[i]-> SetBinContent( ipu , ntot  );

      hist_pass[i]->SetBinError( ipu , sqrt(npass) );
      hist_tot[i]-> SetBinError( ipu , sqrt(ntot)  );
      
    }

    hist_pass[i]->Rebin(rebin);
    hist_tot[i]->Rebin(rebin);
    
    hist_pass[i] -> Divide( hist_tot[i] );
    
    hist_pass[i]->SetLineColor(colors[i]);
    hist_pass[i]->SetMarkerColor(colors[i]);
    hist_pass[i]->SetMarkerStyle(styles[i]);   
    //if( i==0 ) hist_pass[i] -> SetMarkerSize(2);
    hist_pass[i]->GetXaxis()->SetTitle("N_{PV}");
    hist_pass[i]->GetYaxis()->SetTitle(Form("eff(MET > %.0f GeV)",metcut));
    hist_pass[i]->GetYaxis()->SetTitleOffset(1.2);
    hist_pass[i]->SetTitle( title );

    if( fabs( metcut - 30 ) < 0.1 ){
      hist_pass[i]->SetMinimum(0);
      hist_pass[i]->SetMinimum(0.1);
      hist_pass[i]->GetYaxis()->SetRangeUser(0,0.1);
    }
    if( fabs( metcut - 45 ) < 0.1 ){
      hist_pass[i]->SetMinimum(0);
      hist_pass[i]->SetMinimum(0.01);
      hist_pass[i]->GetYaxis()->SetRangeUser(0,0.005);
    }

    
     // if( TString(labels.at(i)).Contains("MC") ){
     //   if( i == 0 ) hist_pass[i]->Draw("hist");
     //   else         hist_pass[i]->Draw("samehist");
     //   leg->AddEntry(hist_pass[i] , labels.at(i) , "l" );    
     // }
     // else{
      if( i == 0 ) hist_pass[i]->Draw("E1");
      else         hist_pass[i]->Draw("sameE1");
      leg->AddEntry(hist_pass[i] , labels.at(i) , "p" );    
      //    }

  }

  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
  
}




// void effVsPU( char* sample, TChain *ch , vector<char*> metTypes , TCut mysel , float metcut ){

//   int colors[4] = {1,2,4,8};
//   int styles[4] = {23,25,20,28};
 
//   TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
  
//   const unsigned int nMetTypes = metTypes.size();

//   TGraphErrors* gr[nMetTypes];

//   for( unsigned int imet = 0 ; imet < nMetTypes ; ++imet ){
    
//     float eff[15];
//     float efferr[15];
//     float n[15];
//     float nerr[15];

//     for( unsigned int ipu = 1 ; ipu < 16 ; ++ipu ){
      
//       TCut npv   (Form("nvtx==%i" ,ipu                        ));
//       TCut met45 (Form("%s>%f"    ,metTypes.at(imet) , metcut ));
      
//       float npass = ch->GetEntries(mysel+npv+met45);
//       float ntot  = ch->GetEntries(mysel+npv);

//       if( ntot > 0. ){
//         eff[ipu-1]    = npass / ntot;
//         efferr[ipu-1] = sqrt(npass) / ntot;
//       }else{
//         eff[ipu-1]    = 0.;
//         efferr[ipu-1] = 0.;
//       }
//       n[ipu-1]      = ipu;
//       nerr[ipu-1]   = 0.5;
      
//     }
    
//     gr[imet] = new TGraphErrors( 15 , n , eff , nerr , efferr );
//     gr[imet]->SetLineColor(colors[imet]);
//     gr[imet]->SetMarkerColor(colors[imet]);
//     gr[imet]->SetMarkerStyle(styles[imet]);    
//     gr[imet]->GetXaxis()->SetTitle("N_{PV}");
//     gr[imet]->GetYaxis()->SetTitle(Form("eff(MET > %.0f GeV)",metcut));
//     gr[imet]->GetYaxis()->SetTitleOffset(1.2);
//     gr[imet]->SetTitle( sample );
    
    
//     if( imet == 0 ) gr[imet]->Draw("AP");
//     else            gr[imet]->Draw("sameP");
    
//     leg->AddEntry(gr[imet] , metTypes.at(imet) , "p" );
//   }

//   leg->SetBorderSize(1);
//   leg->SetFillColor(0);
//   leg->Draw();
  
// }

void effVsPU( char* sample, TChain *ch , vector<char*> metTypes , TCut mysel , float metcut , int rebin ){

  int colors[4] = {1,2,4,8};
  int styles[4] = {23,25,20,28};
 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
  
  const unsigned int nMetTypes = metTypes.size();

  TH1F* hist_tot[nMetTypes];
  TH1F* hist_pass[nMetTypes];

  for( unsigned int imet = 0 ; imet < nMetTypes ; ++imet ){
    
    hist_pass[imet] = new TH1F(Form("hist_pass_%i",imet), Form("hist_pass_%i",imet), 15,0.5,15.5); 
    hist_tot[imet]  = new TH1F(Form("hist_tot_%i",imet),  Form("hist_tot_%i",imet),  15,0.5,15.5);
    hist_pass[imet]->Sumw2();
    hist_tot[imet]->Sumw2();
    
    for( unsigned int ipu = 1 ; ipu < 16 ; ++ipu ){
      
      TCut npv   (Form("nvtx==%i" ,ipu                        ));
      TCut met45 (Form("%s>%f"    ,metTypes.at(imet) , metcut ));
      
      float npass = ch->GetEntries(mysel+npv+met45);
      float ntot  = ch->GetEntries(mysel+npv);

      hist_pass[imet]->SetBinContent( ipu , npass );
      hist_tot[imet]-> SetBinContent( ipu , ntot  );

      hist_pass[imet]->SetBinError( ipu , sqrt(npass) );
      hist_tot[imet]-> SetBinError( ipu , sqrt(ntot)  );
      
    }

    hist_pass[imet]->Rebin(rebin);
    hist_tot[imet]->Rebin(rebin);
    
    hist_pass[imet] -> Divide( hist_tot[imet] );
    
    hist_pass[imet]->SetLineColor(colors[imet]);
    hist_pass[imet]->SetMarkerColor(colors[imet]);
    hist_pass[imet]->SetMarkerStyle(styles[imet]);    
    hist_pass[imet]->GetXaxis()->SetTitle("N_{PV}");
    hist_pass[imet]->GetYaxis()->SetTitle(Form("eff(MET > %.0f GeV)",metcut));
    hist_pass[imet]->GetYaxis()->SetTitleOffset(1.2);
    hist_pass[imet]->SetTitle( sample );
    if( fabs( metcut - 30 ) < 0.1 ){
      hist_pass[imet]->SetMinimum(0);
      hist_pass[imet]->SetMinimum(0.1);
      hist_pass[imet]->GetYaxis()->SetRangeUser(0,0.1);
    }
    if( fabs( metcut - 45 ) < 0.1 ){
      hist_pass[imet]->SetMinimum(0);
      hist_pass[imet]->SetMinimum(0.01);
      hist_pass[imet]->GetYaxis()->SetRangeUser(0,0.01);
    }
    
    if( imet == 0 ) hist_pass[imet]->Draw("E1");
    else            hist_pass[imet]->Draw("sameE1");

    if( fabs( metcut - 30 ) < 0.1 ){
      hist_pass[imet]->SetMinimum(0);
      hist_pass[imet]->SetMinimum(0.1);
      hist_pass[imet]->GetYaxis()->SetRangeUser(0,0.1);
    }
    if( fabs( metcut - 45 ) < 0.1 ){
      hist_pass[imet]->SetMinimum(0);
      hist_pass[imet]->SetMinimum(0.01);
      hist_pass[imet]->GetYaxis()->SetRangeUser(0,0.01);
    }

    leg->AddEntry(hist_pass[imet] , metTypes.at(imet) , "p" );    
  }

  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
  
}



void printChain( TChain *ch , TCut sel ){

//   cout << endl << "DY -> mm" << endl;
//   printChain( dymm );

//   cout << endl << "DY -> ee" << endl;
//   printChain( dyee );

//   cout << endl << "H130" << endl;
//   printChain( h130 );

  TH1F* h_zmet     = new TH1F("h_zmet",    "",100,0,200);
  TH1F* h_hmet     = new TH1F("h_hmet",    "",100,0,200);
  TH1F* h_hmetpf   = new TH1F("h_hmetpf",  "",100,0,200);
  TH1F* h_hmetpf0  = new TH1F("h_hmetpf0", "",100,0,200);
  TH1F* h_hmetpf1  = new TH1F("h_hmetpf1", "",100,0,200);
  TH1F* h_hmetpf2  = new TH1F("h_hmetpf2", "",100,0,200);
  TH1F* h_hmetpf3  = new TH1F("h_hmetpf3", "",100,0,200);
  TH1F* h_hmetpf4  = new TH1F("h_hmetpf4", "",100,0,200);
  TH1F* h_hmetpf5  = new TH1F("h_hmetpf5", "",100,0,200);
  TH1F* h_hmetpf6  = new TH1F("h_hmetpf6", "",100,0,200);
  TH1F* h_hmetpf7  = new TH1F("h_hmetpf7", "",100,0,200);
  TH1F* h_hmetpf8  = new TH1F("h_hmetpf8", "",100,0,200);
  TH1F* h_hmetpf9  = new TH1F("h_hmetpf9", "",100,0,200);
  TH1F* h_hmetpf10 = new TH1F("h_hmetpf10","",100,0,200);

  ch->Draw("TMath::Min(zmet,199.9)    >>h_zmet",sel);
  ch->Draw("TMath::Min(hmet,199.9)    >>h_hmet",sel);
  ch->Draw("TMath::Min(hmetpf,199.9)  >>h_hmetpf",sel);
  ch->Draw("TMath::Min(hmetpf0,199.9) >>h_hmetpf0",sel);
  ch->Draw("TMath::Min(hmetpf1,199.9) >>h_hmetpf1",sel);
  ch->Draw("TMath::Min(hmetpf2,199.9) >>h_hmetpf2",sel);
  ch->Draw("TMath::Min(hmetpf3,199.9) >>h_hmetpf3",sel);
  ch->Draw("TMath::Min(hmetpf4,199.9) >>h_hmetpf4",sel);
  ch->Draw("TMath::Min(hmetpf5,199.9) >>h_hmetpf5",sel);
  ch->Draw("TMath::Min(hmetpf6,199.9) >>h_hmetpf6",sel);
  ch->Draw("TMath::Min(hmetpf7,199.9) >>h_hmetpf7",sel);
  ch->Draw("TMath::Min(hmetpf8,199.9) >>h_hmetpf8",sel);
  ch->Draw("TMath::Min(hmetpf9,199.9) >>h_hmetpf9",sel);
  ch->Draw("TMath::Min(hmetpf10,199.9)>>h_hmetpf10",sel);

  gPad->SetLogy(1);
   
  float integral = h_zmet->Integral();
  int bin30      = h_zmet->FindBin(30);

  //cout << "Efficiences met > 30 GeV" << endl;
  //cout << "zmet       " << Form("%.3f",h_zmet->Integral(bin30,10000)     / integral) << endl;
  //cout << "hmet       " << Form("%.3f",h_hmet->Integral(bin30,10000)     / integral) << endl;
  cout << "no neutrals  " << Form("%.3f",h_hmetpf->Integral(bin30,10000)   / integral) << endl;
  cout << "pt > 0 GeV   " << Form("%.3f",h_hmetpf0->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 1 GeV   " << Form("%.3f",h_hmetpf1->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 2 GeV   " << Form("%.3f",h_hmetpf2->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 3 GeV   " << Form("%.3f",h_hmetpf3->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 4 GeV   " << Form("%.3f",h_hmetpf4->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 5 GeV   " << Form("%.3f",h_hmetpf5->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 6 GeV   " << Form("%.3f",h_hmetpf6->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 7 GeV   " << Form("%.3f",h_hmetpf7->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 8 GeV   " << Form("%.3f",h_hmetpf8->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 9 GeV   " << Form("%.3f",h_hmetpf9->Integral(bin30,10000)  / integral) << endl;
  cout << "pt > 10 GeV  " << Form("%.3f",h_hmetpf10->Integral(bin30,10000) / integral) << endl;

  delete h_zmet;     
  delete h_hmet;     
  delete h_hmetpf;   
  delete h_hmetpf0;  
  delete h_hmetpf1;  
  delete h_hmetpf2;  
  delete h_hmetpf3;  
  delete h_hmetpf4;  
  delete h_hmetpf5;  
  delete h_hmetpf6;  
  delete h_hmetpf7;  
  delete h_hmetpf8;  
  delete h_hmetpf9;  
  delete h_hmetpf10; 
}

void plotHist(TCanvas *can , TH1F *hA, TH1F* hB,char* title, char* var, 
              char* leg1, char* leg2 , bool norm , bool text ){


  can->cd();

  TH1F* h1 = (TH1F*) hA->Clone();
  TH1F* h2 = (TH1F*) hB->Clone();

  if( norm ){
    h1->Scale( 1. / h1->Integral() );
    h2->Scale( 1. / h2->Integral() );
  }

  TLegend *leg = new TLegend(0.65,0.8,0.95,0.95);
  leg->AddEntry(h1,leg1,"f");
  leg->AddEntry(h2,leg2);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);

  float integral = h1->Integral();
  int bin30      = h1->FindBin(30);
  int bin45      = h1->FindBin(45);
  int bin50      = h1->FindBin(50);

  float h1_eff30 = h1->Integral(bin30,10000)/integral;
  float h2_eff30 = h2->Integral(bin30,10000)/integral;
  float h1_eff45 = h1->Integral(bin45,10000)/integral;
  float h2_eff45 = h2->Integral(bin45,10000)/integral;
  float h1_eff50 = h1->Integral(bin50,10000)/integral;
  float h2_eff50 = h2->Integral(bin50,10000)/integral;

  float h1_eff30_err = sqrt( h1->Integral(bin30,10000) )/integral;
  float h2_eff30_err = sqrt( h2->Integral(bin30,10000) )/integral;
  float h1_eff45_err = sqrt( h1->Integral(bin45,10000) )/integral;
  float h2_eff45_err = sqrt( h2->Integral(bin45,10000) )/integral;
  float h1_eff50_err = sqrt( h1->Integral(bin50,10000) )/integral;
  float h2_eff50_err = sqrt( h2->Integral(bin50,10000) )/integral;

  float n1_30 = h1->Integral(bin30,1000);
  float n1_45 = h1->Integral(bin45,1000);
  float n2_30 = h2->Integral(bin30,1000);
  float n2_45 = h2->Integral(bin45,1000);

  h1->Rebin(4);
  h2->Rebin(4);

  h1->SetMinimum(0.1);
  h1->SetMaximum(3000);
  h1->Draw("hist");
  h1->SetFillColor(kAzure-8);
  //h1->Draw("hist");
  h1->SetLineColor(4);
  h1->SetLineWidth(2);
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(var);
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(0.7);
  h2->Draw("sameE1");
  h1->SetMinimum(0.1);
  h1->SetMaximum(3000);


  if( text ){

    TLatex *t = new TLatex();
    t->SetNDC();
    
    t->SetTextSize(0.06);
    
    t->SetTextColor(4);
    t->DrawLatex(0.6,0.75,Form("N(>30) = %.1f",n1_30));
    t->SetTextColor(2);
    t->DrawLatex(0.6,0.67,Form("N(>30) = %.0f",n2_30));
    
    t->SetTextColor(4);
    t->DrawLatex(0.6,0.57,Form("N(>45) = %.1f",n1_45));
    t->SetTextColor(2);
    t->DrawLatex(0.6,0.47,Form("N(>45) = %.0f",n2_45));
    
    
//     if( strcmp( title , "Higgs130" ) == 0 ){
//       t->SetTextColor(4);
//       t->DrawLatex(0.55,0.70,Form("eff30 = %.3f",h1_eff30));
//       t->SetTextColor(2);
//       t->DrawLatex(0.55,0.62,Form("eff30 = %.3f",h2_eff30));
      
//       t->SetTextColor(4);
//       t->DrawLatex(0.55,0.52,Form("eff45 = %.3f",h1_eff45));
//       t->SetTextColor(2);
//       t->DrawLatex(0.55,0.44,Form("eff45 = %.3f",h2_eff45));
//     }
//     else{
//       t->SetTextColor(4);
//       t->DrawLatex(0.45,0.70,Form("eff30 = %.4f #pm %.4f",h1_eff30,h1_eff30_err));
//       t->SetTextColor(2);
//       t->DrawLatex(0.45,0.62,Form("eff30 = %.4f #pm %.4f",h2_eff30,h2_eff30_err));
      
//       t->SetTextColor(4);
//       t->DrawLatex(0.45,0.52,Form("eff45 = %.4f #pm %.4f",h1_eff45,h1_eff45_err));
//       t->SetTextColor(2);
//       t->DrawLatex(0.45,0.44,Form("eff45 = %.4f #pm %.4f",h2_eff45,h2_eff45_err));
//     }
  }

  gPad->SetLogy(1);

  leg->Draw();

  if( printgif_ ) can->Print(Form("plots/%s.gif",title));
}
