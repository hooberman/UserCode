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
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include <sstream>


bool printgif_ = false;
using namespace std;

void printChain(TChain* ch, TCut sel);
void plotHist(TCanvas *can ,TH1F *hA, TH1F* hB, char* title, char* var, char* leg1, char* leg2 , 
	      float metcut , bool norm = false, bool text = false );
void effVsPU( char* sample, TChain *ch , vector<char*> metTypes , TCut mysel, float metcut , int rebin = 1);

void effVsPU( char* title, vector<TChain*> ch , vector<char*> metTypes , vector<char*> labels, 
              TCut mysel , float metcut , int rebin );

void ROC(TChain *chsig, TChain *chbkg, vector<char*> metTypes , 
	 TCut sel, vector<char*> labels );


void compareMET( bool printgif = false ){

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
  
  TChain *h130 = new TChain("Events");
  h130->Add(Form("%s/h130_default_baby.root",path));

  TChain *chdata = new TChain("Events");
  chdata->Add(Form("%s/data_default_baby.root",path));

  cout << "DATA " << chdata->GetEntries() << endl;

  TChain *chmc = new TChain("Events");
  chmc->Add(Form("%s/dymm_default_baby.root",path));
  chmc->Add(Form("%s/dyee_default_baby.root",path));
  //chmc->Add(Form("%s/zjets_spring11_default_baby.root",path));
  //chmc->Add(Form("%s/tt_spring11_default_baby.root",path));
  
  cout << "MC " << chmc->GetEntries() << endl;

  TCut sel("njets30==1&&dilep.mass()>76&&dilep.mass()<106");
  //TCut sel("njets30==1&&dilep.mass()>76&&dilep.mass()<106");
  //TCut sel("njets30==0&&dilep.mass()>81&&dilep.mass()<101");
  
  TCut type("leptype==0||(leptype==3&&lep1.pt()>27)");
  //TCut type("leptype==3&&lep1.pt()>27");
  //TCut type("leptype==0");
  //TCut weight("weight*vtxweight");
  //TCut weight("weight");
  TCut weight("(7429./6555.1)*weight*davtxweight"); //powheg
  //TCut weight("(7429./6352)*10*weight*vtxweight");   //MG
  TCut seltype = sel+type;
  /*
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  
  int nbins  = 100;
  float xmin = 0.;
  float xmax = 100.;
  
  TH1F* hmc_tcmet     = new TH1F("hmc_tcmet",  "", nbins,xmin,xmax);
  TH1F* hdata_tcmet   = new TH1F("hdata_tcmet","", nbins,xmin,xmax);

  hmc_tcmet->Sumw2();
  hdata_tcmet->Sumw2();

  TH1F* hmc_pfmet     = new TH1F("hmc_pfmet",  "", nbins,xmin,xmax);
  TH1F* hdata_pfmet   = new TH1F("hdata_pfmet","", nbins,xmin,xmax);

  hmc_pfmet->Sumw2();
  hdata_pfmet->Sumw2();

  TH1F* hmc_hmetpf     = new TH1F("hmc_hmetpf",  "", nbins,xmin,xmax);
  TH1F* hdata_hmetpf   = new TH1F("hdata_hmetpf","", nbins,xmin,xmax);

  hmc_hmetpf->Sumw2();
  hdata_hmetpf->Sumw2();
 
  TH1F* hmc_hmetpf4     = new TH1F("hmc_hmetpf4",  "", nbins,xmin,xmax);
  TH1F* hdata_hmetpf4   = new TH1F("hdata_hmetpf4","", nbins,xmin,xmax);

  hmc_hmetpf4->Sumw2();
  hdata_hmetpf4->Sumw2();

  TH1F* hmc_jetzmet     = new TH1F("hmc_jetzmet",  "", nbins,xmin,xmax);
  TH1F* hdata_jetzmet   = new TH1F("hdata_jetzmet","", nbins,xmin,xmax);

  hmc_jetzmet->Sumw2();
  hdata_jetzmet->Sumw2();
 
  //char* met1 = "tcmet";
  //char* met1 = "hmetpf";
  //char* met1 = "TMath::Min(pfmet,hmetpf)";
  
  chmc->Draw  ("TMath::Min(tcmet,99.9) >> hmc_tcmet"  ,(sel+type)*weight);
  chdata->Draw("TMath::Min(tcmet,99.9) >> hdata_tcmet",(sel+type)       );

  chmc->Draw  ("TMath::Min(pfmet,99.9) >> hmc_pfmet"  ,(sel+type)*weight);
  chdata->Draw("TMath::Min(pfmet,99.9) >> hdata_pfmet",(sel+type)       );
  
  chmc->Draw  ("TMath::Min(hmetpf,99.9) >> hmc_hmetpf"  ,(sel+type)*weight);
  chdata->Draw("TMath::Min(hmetpf,99.9) >> hdata_hmetpf",(sel+type)       );
  
  chmc->Draw  ("TMath::Min(hmetpf4,99.9) >> hmc_hmetpf4"  ,(sel+type)*weight);
  chdata->Draw("TMath::Min(hmetpf4,99.9) >> hdata_hmetpf4",(sel+type)       );

  chmc->Draw  ("TMath::Min(jetzmet,99.9) >> hmc_jetzmet"  ,(sel+type)*weight);
  chdata->Draw("TMath::Min(jetzmet,99.9) >> hdata_jetzmet",(sel+type)       );
  
  cout << "Tot MC   " << hmc_tcmet->Integral() << endl;
  cout << "Tot data " << hdata_tcmet->Integral() << endl;

  TH1F* hmc_zpt     = new TH1F("hmc_zpt",  "", nbins,xmin,xmax);
  TH1F* hdata_zpt   = new TH1F("hdata_zpt","", nbins,xmin,xmax);

  hmc_zpt->Sumw2();
  hdata_zpt->Sumw2();

  chmc->Draw  ("TMath::Min(dilep.pt(),99.9) >> hmc_zpt"  ,(sel+type)*weight);
  chdata->Draw("TMath::Min(dilep.pt(),99.9) >> hdata_zpt",(sel+type)       );

  delete ctemp;
  */
  //----------------------------------
  // tcmet, pfmet plots
  //----------------------------------

  /*

  TCanvas *c1 = new TCanvas();  
  c1->cd();
  plotHist( c1 , hmc_tcmet , hdata_tcmet , "" , "tcmet (GeV)" , "MC" , "data" , 30 , false , true);
  
  if( printgif ) c1->Print("plots/tcmet_datamc_mm.gif");

  
  TCanvas *c2 = new TCanvas();
  c2->cd();
  plotHist( c2 , hmc_pfmet , hdata_pfmet , "" , "pfmet (GeV)" , "MC" , "data" , 30 , false , true );
  
  if( printgif ) c2->Print("plots/pfmet_datamce_mm.gif");
  
  
  TCanvas *c3 = new TCanvas();
  c3->cd();
  plotHist( c3 , hmc_hmetpf , hdata_hmetpf , "" , "track-MET (GeV)" , "MC" , "data" , 30 , false , true );
  
  if( printgif ) c3->Print("plots/trkmet_datamc.gif");
  
  TCanvas *c4 = new TCanvas();
  c4->cd();
  plotHist( c4 , hmc_hmetpf4 , hdata_hmetpf4 , "" , "track/neu-MET (GeV)" , "MC" , "data" , 30 , false , true );
  
  if( printgif ) c4->Print("plots/trkneumet_datamc.gif");

  
  TCanvas *c5 = new TCanvas();
  c5->cd();
  plotHist( c5 , hmc_jetzmet , hdata_jetzmet , "" , "jet/Z MET (GeV)" , "MC" , "data" , 30 , false , true );
  
  if( printgif ) c5->Print("plots/jetzmet_datamc.gif");
  */  


  /*
  //----------------------------------
  // Z pt
  //----------------------------------
  
  // TCanvas *c2 = new TCanvas("c1","",600,600);
 
  // c2->cd();
  // plotHist( c2 , hmc_zpt , hdata_zpt , "" , "Z p_{T} (GeV)" , "MC" , "data" , false , false );



  vector<TChain*> ch;
  ch.push_back(chmc);
  ch.push_back(chdata);

  vector<TChain*> ch2;
  ch2.push_back(chmc);
  ch2.push_back(chdata);
  ch2.push_back(chmc);
  ch2.push_back(chdata);

  
  //------------------------------
  // pfmet
  //------------------------------
 
  vector<char*> pfmetTypes;
  pfmetTypes.push_back("pfmet");
  pfmetTypes.push_back("pfmet");
  
  vector<char*> pflabels;
  pflabels.push_back("MC (pfmet)");
  pflabels.push_back("data (pfmet)");

  //TCanvas *c2 = new TCanvas();
  //c2->cd();
  //effVsPU( "Spring11 MC vs. data", ch , pfmetTypes , pflabels , seltype , 30. , 1 );
  //if( printgif ) c2->Print("plots/DA_pfmet30.gif");

  TCanvas *c3 = new TCanvas();
  c3->cd();
  effVsPU( "Spring11 MC vs. data", ch , pfmetTypes , pflabels , seltype , 45. , 1 );
  if( printgif ) c3->Print("plots/DA_pfmet45.gif");


  //------------------------------
  // tcmet
  //------------------------------
  
  vector<char*> tcmetTypes;
  tcmetTypes.push_back("tcmet");
  tcmetTypes.push_back("tcmet");
  
  vector<char*> tclabels;
  tclabels.push_back("MC (tcmet)");
  tclabels.push_back("data (tcmet)");

  TCanvas *c4 = new TCanvas();
  c4->cd();
  effVsPU( "Spring11 MC vs. data", ch , tcmetTypes , tclabels , seltype , 30. , 1 );
  if( printgif ) c4->Print("plots/DA_tcmet30.gif");

  TCanvas *c5 = new TCanvas();
  c5->cd();
  effVsPU( "Spring11 MC vs. data", ch , tcmetTypes , tclabels , seltype , 45. , 1 );
  if( printgif ) c5->Print("plots/DA_tcmet45.gif");
  */

  //------------------------------
  // pfmet+trk-MET
  //------------------------------
  /*
  vector<char*> pftrkmetTypes;
  pftrkmetTypes.push_back("pfmet");
  pftrkmetTypes.push_back("pfmet");
  //pftrkmetTypes.push_back("hmetpf");
  //pftrkmetTypes.push_back("hmetpf");
  //pftrkmetTypes.push_back("hmetpf4");
  //pftrkmetTypes.push_back("hmetpf4");
  pftrkmetTypes.push_back("TMath::Min(pfmet,hmetpf)");
  pftrkmetTypes.push_back("TMath::Min(pfmet,hmetpf)");
  
  vector<char*> pftrklabels;
  pftrklabels.push_back("MC (pfmet)");
  pftrklabels.push_back("data (pfmet)");
  //pftrklabels.push_back("MC (trk-MET)");
  //pftrklabels.push_back("data (trk-MET)");
  //pftrklabels.push_back("MC (trk/neu-MET)");
  //pftrklabels.push_back("data (trk/neu-MET)");
  pftrklabels.push_back("MC (min-MET)");
  pftrklabels.push_back("data (min-MET)");

  TCanvas *c6 = new TCanvas();
  c6->cd();
  effVsPU( "Spring11 MC vs. data", ch2 , pftrkmetTypes , pftrklabels , seltype , 30. , 1 );
  if( printgif ) c6->Print("plots/DA_pftrkmet30.gif");

  TCanvas *c7 = new TCanvas();
  c7->cd();
  effVsPU( "Spring11 MC vs. data", ch2 , pftrkmetTypes , pftrklabels , seltype , 45. , 1 );
  if( printgif ) c7->Print("plots/DA_pftrkmet45.gif");
  */

 
  /*

  vector<char*> metTypes;
  //metTypes.push_back("pfmet");
  //metTypes.push_back("pfmet");
  metTypes.push_back("tcmet");
  metTypes.push_back("tcmet");
  //metTypes.push_back("hmetpf");
  //metTypes.push_back("hmetpf");
  //metTypes.push_back("hmetpf4");
  //metTypes.push_back("hmetpf4");
  //metTypes.push_back("TMath::Min(pfmet,hmetpf)");
  //metTypes.push_back("TMath::Min(pfmet,hmetpf)");
  //metTypes.push_back("jetzmet");
  //metTypes.push_back("jetzmet");

  vector<char*> labels;
  //labels.push_back("MC (pfmet)");
  //labels.push_back("data (pfmet)");
  labels.push_back("MC (tcmet)");
  labels.push_back("data (tcmet)");
  //labels.push_back("MC (trk-MET)");
  //labels.push_back("data (trk-MET)");
  //labels.push_back("MC (trk/neu-MET)");
  //labels.push_back("data (trk/neu-MET)");
  //labels.push_back("MC (min-MET)");
  //labels.push_back("data (min-MET)");
  //labels.push_back("MC (jet/Z-MET)");
  //labels.push_back("data (jet/Z-MET)");

  */


  //---------------------------------------
  // ROC curves pfmet, track-met, min-MET
  //---------------------------------------
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

  //TCut jetveto1("njets30==0&&nvtx<8&&dilep.mass()>12&&dilep.mass()<76");
  TCut jetveto1("njets30==0&&nvtx<8");

  ROC(h130,dymm,rocMetTypes,jetveto1,rocMetLabels);

  c4->Print("plots/ROC1_m76.gif");
  //c4->Print("plots/ROC1.gif");

  TCanvas *c5 = new TCanvas();
  c5->cd();
  
  //TCut jetveto2("njets30==0&&nvtx>7&&dilep.mass()>12&&dilep.mass()<76");
  TCut jetveto2("njets30==0&&nvtx>7");
  
  ROC(h130,dymm,rocMetTypes,jetveto2,rocMetLabels);
  
  c5->Print("plots/ROC2_m76.gif");
  //c5->Print("plots/ROC2.gif");
  */



  /*
  vector<char*> rocMetTypes;
  rocMetTypes.push_back("pfmet");
  rocMetTypes.push_back("jetzmet");
  rocMetTypes.push_back("TMath::Min(pfmet,jetzmet)");

  vector<char*> rocMetLabels;
  rocMetLabels.push_back("pfmet");
  rocMetLabels.push_back("jet-Z MET");
  rocMetLabels.push_back("min-MET");

  TCanvas *c6 = new TCanvas();
  c6->cd();

  TCut jet1_loPU("njets30==1&&jetpv==1&&nvtx<8&&dilep.mass()>76&&dilep.mass()<106");

  ROC(h130,dymm,rocMetTypes,jet1_loPU,rocMetLabels);

  c6->Print("plots/jetzmet_loPU.gif");
  
  TCanvas *c7 = new TCanvas();
  c7->cd();
  
  TCut jet1_hiPU("njets30==1&&jetpv==1&&nvtx>7&&dilep.mass()>76&&dilep.mass()<106");
  
  ROC(h130,dymm,rocMetTypes,jet1_hiPU,rocMetLabels);
  
  c7->Print("plots/jetzmet_hiPU.gif");
  */


  /*
  
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
  

*/
  //------------------------------------------
  // met flavors in 1-jet bin
  //------------------------------------------

  TCanvas *effcan2 = new TCanvas("effcan2","",1200,600);
  effcan2->Divide(2,1);

  vector<char*> metTypes2;
  metTypes2.push_back("pfmet");
  //metTypes2.push_back("hmetpf4");
  //metTypes2.push_back("jetzmet");
  //metTypes2.push_back("jetzmet4");
  //metTypes2.push_back("jetzmet8");
  metTypes2.push_back("TMath::Min(pfmet,jetzmet)");

  TCut dphizjet("acos(cos(dilep.phi()-jet.phi()))>2");
  TCut lepveto("nlep==0");
  TCut zmass("dilep.mass()>76&&dilep.mass()<106");
  TCut ee("leptype==0");
  TCut mm("leptype==3");

  TCut sel1          = zmass+(mm||ee)+lepveto+dphizjet+"njets30==1 && jetpv==1";

  //TCut mysel2("njets30==1 && (leptype==0||leptype==3) && jetpv==1 && dilep.mass()>76 && dilep.mass()<106");

  effcan2->cd(1);
  effVsPU( "DY" , chmc , metTypes2 , sel1 , 30. );

  effcan2->cd(2);
  effVsPU( "data" , chdata , metTypes2 , sel1 , 30. );




  //effVsPU( "DY#rightarrowee"     , dyee , metTypes2 , mysel2 , 45. , 5 );

  
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
  TLegend *leg = new TLegend(0.7,0.2,0.9,0.5);

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
 
  TLegend *leg = new TLegend(0.2,0.6,0.6,0.9);
  
  const unsigned int nMetTypes = metTypes.size();

  TH1F* hist_tot[nMetTypes];
  TH1F* hist_pass[nMetTypes];

  TH1F* htot_temp  = new TH1F("htot_temp","",1,0,1);
  TH1F* hpass_temp = new TH1F("hpass_temp","",1,0,1);

  htot_temp->Sumw2();
  hpass_temp->Sumw2();

  int npasstot = 0;

  for( unsigned int i = 0 ; i < nMetTypes ; ++i ){
    
    hist_pass[i] = new TH1F(Form("hist_pass_%i",i), Form("hist_pass_%i",i), 15,0.5,15.5); 
    hist_tot[i]  = new TH1F(Form("hist_tot_%i",i),  Form("hist_tot_%i",i),  15,0.5,15.5);
    hist_pass[i]->Sumw2();
    hist_tot[i]->Sumw2();
    
    for( unsigned int ipu = 1 ; ipu < 16 ; ++ipu ){
      
      //TCut npv   (Form("nvtx==%i" ,ipu                        ));
      TCut npv   (Form("ndavtx==%i" ,ipu                        ));
      TCut met45 (Form("%s>%f"    ,metTypes.at(i) , metcut ));
      TCut myweight("1");
      TCut mymcweight("weight*davtxweight");
      if( TString(labels.at(i)).Contains("MC") ) myweight=mymcweight;
      //if( TString(labels.at(i)).Contains("MC") ) myweight="weight";
      
      //float npass = ch.at(i)->GetEntries((mysel+npv+met45)*myweight);
      //float ntot  = ch.at(i)->GetEntries((mysel+npv)*myweight);
      //float npasserr = sqrt(npass);
      //float ntoterr  = sqrt(ntot);

      ch.at(i)->Draw("0.5>>hpass_temp" , (mysel+npv+met45) * myweight);
      ch.at(i)->Draw("0.5>>htot_temp"  , (mysel+npv)       * myweight);
      
      float npass = hpass_temp->GetBinContent(1);
      float ntot  = htot_temp->GetBinContent(1);

      npasstot += npass;

      float npasserr = hpass_temp->GetBinError(1);
      float ntoterr  = htot_temp->GetBinError(1);

      hist_pass[i]->SetBinContent( ipu , npass );
      hist_tot[i]-> SetBinContent( ipu , ntot  );

      hist_pass[i]->SetBinError( ipu , npasserr );
      hist_tot[i]-> SetBinError( ipu , ntoterr  );

      //cout << hist_pass[i]->GetBinContent(ipu) << " +/- " << hist_pass[i]->GetBinError(ipu) << endl;
      //cout << hist_tot[i]->GetBinContent(ipu)  << " +/- " << hist_tot[i]->GetBinError(ipu) << endl;
      
    }

    cout << "npass all " << npasstot << endl;

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
      hist_pass[i]->GetYaxis()->SetLabelSize(0.03);
      hist_pass[i]->SetMinimum(0);
      hist_pass[i]->SetMinimum(0.01);
      hist_pass[i]->GetYaxis()->SetRangeUser(0,0.015);
    }
  }

  for( unsigned int i = 0 ; i < nMetTypes ; ++i ){
    
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

  float npassall = 0.;

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

      npassall += npass;

      hist_pass[imet]->SetBinContent( ipu , npass );
      hist_tot[imet]-> SetBinContent( ipu , ntot  );

      hist_pass[imet]->SetBinError( ipu , sqrt(npass) );
      hist_tot[imet]-> SetBinError( ipu , sqrt(ntot)  );
      
    }

    cout << "npass all " << npassall << endl;

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
              char* leg1, char* leg2 , float metcut , bool norm , bool text ){


  //can->cd();

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

  int mybin = h1->FindBin(metcut);

  float intdata = h2->Integral(mybin,10000);
  float intmc   = h1->Integral(mybin,10000);
  float totintdata = h2->Integral();
  float totintmc   = h1->Integral();

  h1->Rebin(2);
  h2->Rebin(2);

  h1->SetMinimum(0.1);
  h1->SetMaximum(3000);
  h1->Draw("hist");
  h1->SetFillColor(5);
  //h1->Draw("hist");
  h1->SetLineColor(1);
  h1->SetLineWidth(1);
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(var);
  h2->SetLineColor(1);
  h2->SetMarkerColor(1);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(0.7);
  h2->Draw("sameE1");
  h1->SetMinimum(0.1);
  h1->SetMaximum(5000);

  if( text ){

    TLatex *t = new TLatex();
    t->SetNDC();

    
    t->SetTextSize(0.035);
    t->DrawLatex(0.67,0.70,Form("N_{DATA} (tot) = %.0f",totintdata));
    t->DrawLatex(0.67,0.64,Form("N_{MC}   (tot) = %.0f",totintmc));
    t->SetTextSize(0.04);
    t->DrawLatex(0.67,0.55,Form("N_{DATA} (>%.0f) = %.0f",metcut,intdata));
    t->DrawLatex(0.67,0.48,Form("N_{MC}   (>%.0f) = %.1f",metcut,intmc));

    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawLine(metcut,0.1,metcut,5000);

  }

  if( text ){

    // TLatex *t = new TLatex();
    // t->SetNDC();
    
    // t->SetTextSize(0.06);
    
    // t->SetTextColor(4);
    // t->DrawLatex(0.6,0.75,Form("N(>30) = %.1f",n1_30));
    // t->SetTextColor(2);
    // t->DrawLatex(0.6,0.67,Form("N(>30) = %.0f",n2_30));
    
    // t->SetTextColor(4);
    // t->DrawLatex(0.6,0.57,Form("N(>45) = %.1f",n1_45));
    // t->SetTextColor(2);
    // t->DrawLatex(0.6,0.47,Form("N(>45) = %.0f",n2_45));
    
    
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

  //if( printgif_ ) can->Print(Form("plots/%s.gif",title));
}
