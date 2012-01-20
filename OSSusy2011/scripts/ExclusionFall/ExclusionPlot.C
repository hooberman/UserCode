#include <sstream>
#include <string>
#include <iomanip>
#include "ExclusionPlot.hh"
 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TMarker.h"
#include <vector>
#include "TMath.h"

const float m0max    = 2000; // x-axis range
const float m12max   =  700; // y-axis range
const float spacing  =  250; // spacing between constant squark/gluino lines
const int   tanBeta_ =   10; // tan beta 

void ExclusionPlot(){
  gStyle->SetPalette(1);

  CommandMSUGRA("ExclusionPlot.root");
}

void CommandMSUGRA(TString plotName_){
  


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 
  gStyle->SetTextFont(42);
  
  //convert tanb value to string
  std::stringstream tmp;
  tmp << tanBeta_;
  TString tanb( tmp.str() );
  
  // Output file
  cout << " create " << plotName_ << endl;
  TFile* output = new TFile( plotName_, "RECREATE" );
  if ( !output || output->IsZombie() ) { std::cout << " zombie alarm output is a zombie " << std::endl; }
  
  //-----------------------------------
  //set old exclusion Limits
  //-----------------------------------

  TGraph* LEP_ch      = set_lep_ch(tanBeta_);         //chargino curve
  TGraph* LEP_sl      = set_lep_sl(tanBeta_);         //slepton curve
  TGraph* TEV_sg_cdf  = set_tev_sg_cdf(tanBeta_);     //squark gluino cdf
  TGraph* TEV_sg_d0   = set_tev_sg_d0(tanBeta_);      //squark gluino d0
  TGraph* stau        = set_tev_stau(tanBeta_);       //stau 
  TGraph* NoEWSB      = set_NoEWSB(tanBeta_);         //noEWSB
  //TGraph* TEV_tlp_cdf = set_tev_tlp_cdf(tanBeta_);    //trilepton cdf
  //TGraph* TEV_tlp_d0  = set_tev_tlp_d0(tanBeta_);     //trilepton d0
  TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1(tanBeta_); //D0 sneutrino
  TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2(tanBeta_); //D0 sneutrino

  //-----------------------------------
  // constant sqquark and gluino lines
  //-----------------------------------

  const unsigned int nslines = 5;
  const unsigned int nglines = 4;

  TF1* lnsq[nslines];
  TF1* lngl[nglines];

  TLatex* sq_text[nslines];
  TLatex* gl_text[nglines];

  for(unsigned int i = 0; i < nslines; i++){
    lnsq[i]    = constant_squark(tanBeta_,i);
    sq_text[i] = constant_squark_text(i,*lnsq[i],tanBeta_);
  }

  for(unsigned int i = 0; i < nglines; i++){
    lngl[i]    = constant_gluino(tanBeta_,i);
    gl_text[i] = constant_gluino_text(i,*lngl[i]);
  }

  //-----------------------------------
  // Legends
  //-----------------------------------

  TLegend* legst     = makeStauLegend(0.05,tanBeta_);
  TLegend* legexp    = makeExpLegend( *TEV_sg_cdf,*TEV_sg_d0,*LEP_ch,*LEP_sl,*TEV_sn_d0_1,0.035,tanBeta_);
  TLegend* legNoEWSB = makeNoEWSBLegend(0.05,tanBeta_);

  //-----------------------------------
  // make Canvas
  //-----------------------------------

  TCanvas* cvsSys = new TCanvas("cvsnm","cvsnm",0,0,800,600);
  gStyle->SetOptTitle(0);
  cvsSys->SetFillColor(0);
  cvsSys->GetPad(0)->SetRightMargin(0.07);
  cvsSys->Range(-120.5298,26.16437,736.0927,500);
  cvsSys->SetFillColor(0);
  cvsSys->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderSize(2);
  cvsSys->GetPad(0)->SetLeftMargin(0.1407035);
  cvsSys->GetPad(0)->SetTopMargin(0.08);
  cvsSys->GetPad(0)->SetBottomMargin(0.13);
  cvsSys->SetTitle("tan#beta="+tanb);
 
  output->cd();

  //-----------------------------------
  // get contours
  //-----------------------------------
  
  TGraph* observedLimit           = observedLimit_OS2011();
  TGraph* expectedLimit           = expectedLimit_OS2011();
  TGraph* expectedLimitM1         = expectedLimitM1_OS2011();
  TGraph* expectedLimitP1         = expectedLimitP1_OS2011();
  TGraph* observedLimitTheoryUp   = observedLimitTheoryUp_OS2011();
  TGraph* observedLimitTheoryDown = observedLimitTheoryDown_OS2011();
  TGraph* expectedLimitTheoryUp   = expectedLimitTheoryUp_OS2011();
  TGraph* expectedLimitTheoryDown = expectedLimitTheoryDown_OS2011();

  //-----------------------------------
  // format contours
  //-----------------------------------

  observedLimit->SetLineWidth(2);
  observedLimit->SetLineColor(2);

  observedLimitTheoryUp->SetLineWidth(2);
  observedLimitTheoryUp->SetLineColor(2);
  observedLimitTheoryUp->SetLineStyle(2);

  observedLimitTheoryDown->SetLineWidth(2);
  observedLimitTheoryDown->SetLineColor(2);
  observedLimitTheoryDown->SetLineStyle(2);

  expectedLimit->SetLineWidth(2);
  expectedLimit->SetLineColor(4);

  expectedLimitTheoryUp->SetLineWidth(2);
  expectedLimitTheoryUp->SetLineColor(4);
  expectedLimitTheoryUp->SetLineStyle(2);

  expectedLimitTheoryDown->SetLineWidth(2);
  expectedLimitTheoryDown->SetLineColor(4);
  expectedLimitTheoryDown->SetLineStyle(2);

  expectedLimitM1->SetLineWidth(2);
  expectedLimitM1->SetLineColor(4);
  expectedLimitM1->SetLineStyle(3);

  expectedLimitP1->SetLineWidth(2);
  expectedLimitP1->SetLineColor(4);
  expectedLimitP1->SetLineStyle(3);



  //-----------------------------------
  // draw dummy histogram
  //-----------------------------------

  double m0min = 0;
  TH2D* hist = new TH2D("h","h",100,m0min,m0max,100,120,m12max);
  hist->Draw();  
  hist->GetXaxis()->SetTitle("m_{0} (GeV/c^{2})");
  hist->GetYaxis()->SetTitle("m_{1/2} (GeV/c^{2})");
  hist->GetYaxis()->SetTitleOffset(1.);
  hist->GetXaxis()->SetNdivisions(506);
  hist->GetYaxis()->SetNdivisions(506);

  //-----------------------------------
  // make legend
  //-----------------------------------
  
  TLegend* myleg = new TLegend(0.3,0.75,0.54,0.9,NULL,"brNDC");
  myleg->SetFillColor(0); 
  myleg->SetShadowColor(0);
  myleg->SetTextSize(0.03);
  myleg->SetBorderSize(0);

  myleg->AddEntry(observedLimit,"NLO observed limit","l");
      
  //-----------------------------------
  // constant squark/gluino contours
  //-----------------------------------

  for (unsigned int it=1;it<nslines;it++) {   
    lnsq[it]->Draw("same");
    sq_text[it]->Draw();
  }

  for (unsigned int it=1;it<nglines;it++) {   
    lngl[it]->Draw("same");   
    gl_text[it]->Draw();
  }

  //-----------------------------------
  // draw exclusion contours
  //-----------------------------------

  observedLimit->Draw("c");
  observedLimitTheoryUp->Draw("samec");
  observedLimitTheoryDown->Draw("samec");
  expectedLimit->Draw("samec");
  expectedLimitTheoryUp->Draw("samec");
  expectedLimitTheoryDown->Draw("samec");
  expectedLimitM1->Draw("samec");
  expectedLimitP1->Draw("samec");

  //-----------------------------------
  // draw prior excluded regions
  //-----------------------------------
    
  LEP_ch->Draw("fsame");
  LEP_sl->Draw("fsame");
  TEV_sg_cdf->Draw("fsame");
  TEV_sg_d0->Draw("same");  
  TEV_sg_d0->Draw("fsame");

  //-----------------------------------
  // other labels
  //-----------------------------------

  Double_t xpos  = 0;
  Double_t xposi = 0;
  Double_t ypos  = 0;
  TLatex* lumilabel = new TLatex(1000,m12max+15,"#sqrt{s} = 7 TeV, #scale[0.6]{#int} L dt = 4.7 fb^{-1}");

  lumilabel->SetTextSize(0.05);
  lumilabel->Draw("same");

  TLatex* cmslabel = new TLatex(10.,m12max+15,"CMS Preliminary");
  cmslabel->SetTextSize(0.05);
  cmslabel->Draw("same");

  TLatex* cmssmpars = new TLatex(120,540,"tan#beta = "+tanb+",  A_{0} = 0,  #mu > 0");
  cmssmpars->SetTextSize(0.045);
  cmssmpars->Draw("same");

  stau->Draw("fsame");
  NoEWSB->Draw("fsame");

  legexp->Draw();
  legst->Draw();
  myleg->Draw();
  legNoEWSB->Draw();

  //-----------------------------------
  // LM points
  //-----------------------------------

  TMarker* LM0 = new TMarker(200.,160.,20);
  TMarker* LM1 = new TMarker(60.,250.,20);
  TMarker* LM3 = new TMarker(330.,240.,20);
  TMarker* LM6 = new TMarker(80.,400.,20);
    
  LM0->SetMarkerSize(1.2);
  LM1->SetMarkerSize(1.2);
    
  TLatex* tLM0 = new TLatex(205.,160.," LM0");
  TLatex* tLM1 = new TLatex(80.,245.,"LM1");
  TLatex* tLM3 = new TLatex(350.,235.,"LM3");
  TLatex* tLM6 = new TLatex(100.,395.,"LM6");

  tLM0->SetTextSize(0.035);
  tLM1->SetTextSize(0.035);
  tLM3->SetTextSize(0.035);
  tLM6->SetTextSize(0.035);
  
  LM1->Draw("same");
  tLM1->Draw("same");
  LM3->Draw("same");
  tLM3->Draw("same");
  LM6->Draw("same");
  tLM6->Draw("same");

  //-----------------------------------
  // write output
  //-----------------------------------
 
  hist->Draw("sameaxis");
  cvsSys->RedrawAxis();
  cvsSys->Update();
  cvsSys->Write();  
  cvsSys->SaveAs("RA6_ExclusionLimit_tanb"+tanb+".eps");
  cvsSys->SaveAs("RA6_ExclusionLimit_tanb"+tanb+".pdf");
  cvsSys->SaveAs("RA6_ExclusionLimit_tanb"+tanb+".png");
  
  output->Write();
  
}


void setPlottingStyle(TH1F& hsig){
  
  hsig.SetStats(kFALSE);
  
  hsig.SetAxisRange(80,500,"Y");
  hsig.SetAxisRange(0,520,"X");
  hsig.SetAxisRange(200,520,"X");

  hsig.GetXaxis()->SetTitle("m_{0} (GeV)");
  hsig.GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hsig.GetYaxis()->SetTitleOffset(0.8);
  hsig.GetYaxis()->SetTitleSize(0.06);
  hsig.GetYaxis()->SetLabelSize(0.06);
  hsig.GetXaxis()->SetTitleOffset(0.9);
  hsig.GetXaxis()->SetTitleSize(0.06);
  hsig.GetXaxis()->SetLabelSize(0.06);

  hsig.SetLineWidth(1);  
  hsig.SetLineColor(kBlue);  
  
}




TGraph* set_sneutrino_d0_1(Int_t tanBeta){

  tanBeta = tanBeta;

  double sn_m0[14]= {0,  0, 48, 55, 80, 90,100,105,109,105,100, 72, 55,0};
  double sn_m12[14]={0,140,210,220,237,241,242,241,230,220,210,170,150,0};

  TGraph* sn_d0_gr = new TGraph(14,sn_m0,sn_m12);

  sn_d0_gr->SetFillColor(kGreen+3);
  sn_d0_gr->SetFillStyle(1001);

  return sn_d0_gr;
}

TGraph* set_sneutrino_d0_2(Int_t tanBeta){

  tanBeta = tanBeta;

  double sn_m0[9]= {0, 45, 75,115,130,150,163,185,0};
  double sn_m12[9]={0,140,170,213,202,183,168,140,0};

  TGraph* sn_d0_gr_2 = new TGraph(9,sn_m0,sn_m12);

  sn_d0_gr_2->SetFillColor(kGreen+3);
  sn_d0_gr_2->SetFillStyle(1001);

  return sn_d0_gr_2;
}

TGraph* set_lep_ch(Int_t tanBeta){
  if(tanBeta ==  3) return set_lep_ch_tanBeta3();
  if(tanBeta == 10) return set_lep_ch_tanBeta10();
  if(tanBeta == 50) return set_lep_ch_tanBeta50();

  return 0;
}

TGraph* set_lep_ch_tanBeta10(){

  const unsigned int npoints = 50;

  double ch_m0[npoints];
  double ch_m12[npoints];

  //----------------------------
  // updated, from Frederic
  //----------------------------
  int i = -1;

  ch_m0[++i] =   0.   ; ch_m12[i] = 163.0;
  ch_m0[++i] = 100.   ; ch_m12[i] = 162.0;
  ch_m0[++i] = 150.   ; ch_m12[i] = 161.3;  
  ch_m0[++i] = 200.   ; ch_m12[i] = 160.7;
  ch_m0[++i] = 250.   ; ch_m12[i] = 160.0;
  ch_m0[++i] = 300.   ; ch_m12[i] = 159.2;
  ch_m0[++i] = 350.   ; ch_m12[i] = 158.4;
  ch_m0[++i] = 400.   ; ch_m12[i] = 157.6;
  ch_m0[++i] = 450.   ; ch_m12[i] = 156.7;
  ch_m0[++i] = 500.   ; ch_m12[i] = 155.8;
  ch_m0[++i] = 550.   ; ch_m12[i] = 154.9;
  ch_m0[++i] = 600.   ; ch_m12[i] = 154.0;
  ch_m0[++i] = 650.   ; ch_m12[i] = 153.1;
  ch_m0[++i] = 700.   ; ch_m12[i] = 152.3;
  ch_m0[++i] = 750.   ; ch_m12[i] = 151.5;
  ch_m0[++i] = 800.   ; ch_m12[i] = 150.7;
  ch_m0[++i] = 850.   ; ch_m12[i] = 150.0;
  ch_m0[++i] = 900.   ; ch_m12[i] = 149.3;
  ch_m0[++i] = 950.   ; ch_m12[i] = 148.7;
  ch_m0[++i] = 1000.  ; ch_m12[i] = 148.2;
  ch_m0[++i] = 1050.  ; ch_m12[i] = 147.7;
  ch_m0[++i] = 1100.  ; ch_m12[i] = 147.3;
  ch_m0[++i] = 1150.  ; ch_m12[i] = 147.1;
  ch_m0[++i] = 1200.  ; ch_m12[i] = 146.8;
  ch_m0[++i] = 1250.  ; ch_m12[i] = 146.7;
  ch_m0[++i] = 1300.  ; ch_m12[i] = 146.6;
  ch_m0[++i] = 1350.  ; ch_m12[i] = 146.6;
  ch_m0[++i] = 1400.  ; ch_m12[i] = 146.7;
  ch_m0[++i] = 1450.  ; ch_m12[i] = 146.9;
  ch_m0[++i] = 1500.  ; ch_m12[i] = 147.1;
  ch_m0[++i] = 1550.  ; ch_m12[i] = 147.5;
  ch_m0[++i] = 1600.  ; ch_m12[i] = 148.1;
  ch_m0[++i] = 1650.  ; ch_m12[i] = 148.7;
  ch_m0[++i] = 1700.  ; ch_m12[i] = 149.5;
  ch_m0[++i] = 1750.  ; ch_m12[i] = 150.5;
  ch_m0[++i] = 1800.  ; ch_m12[i] = 151.7;
  ch_m0[++i] = 1850.  ; ch_m12[i] = 153.1;
  ch_m0[++i] = 1900.  ; ch_m12[i] = 154.7;
  ch_m0[++i] = 1950.  ; ch_m12[i] = 156.6;
  ch_m0[++i] = 2000.  ; ch_m12[i] = 159.1;
  ch_m0[++i] = 2050.  ; ch_m12[i] = 161.5;
  ch_m0[++i] = 2050.  ; ch_m12[i] = 0;
  ch_m0[++i] = 0.     ; ch_m12[i] = 0;


  //---------------------------------------------
  // old region, extended linearly to m0 1 TeV
  //---------------------------------------------
  /*
  ch_m0[0] = 0;
  ch_m0[1] = 100;
  ch_m0[2] = 200;
  ch_m0[3] = 300;
  ch_m0[4] = 400;
  ch_m0[5] = 500;
  ch_m0[6] = 600;
  ch_m0[7] = 700;
  ch_m0[8] = 800; 
  ch_m0[9] = 1200; 
  ch_m0[10] = 1200;
  ch_m0[11] = 0;

  ch_m12[0] = 163;
  ch_m12[1] = 162;
  ch_m12[2] = 161;
  ch_m12[3] = 160;
  ch_m12[4] = 159;
  ch_m12[5] = 158;
  ch_m12[6] = 157;
  ch_m12[7] = 156;
  ch_m12[8] = 155.4;
  ch_m12[9] = 151.6;
  ch_m12[10] = 0;
  ch_m12[11] = 0;
  */
  
  TGraph* ch_gr = new TGraph(npoints,ch_m0,ch_m12);

  ch_gr->SetFillColor(8);
  ch_gr->SetLineColor(8);
  //  ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}



TGraph* set_lep_ch_tanBeta3(){

  double ch_m0[17];
  double ch_m12[17];

  ch_m0[0] = 0;
  ch_m0[1] = 100;
  ch_m0[2] = 150;
  ch_m0[3] = 200;
  ch_m0[4] = 250;
  ch_m0[5] = 300;
  ch_m0[6] = 350;
  ch_m0[7] = 400;
  ch_m0[8] = 450;
  ch_m0[9] = 500;
  ch_m0[10] = 550;
  ch_m0[11] = 600;
  ch_m0[12] = 650;
  ch_m0[13] = 700;
  ch_m0[14] = 750;
  ch_m0[15] = 750;
  ch_m0[16] = 0;
  
  ch_m12[0] = 170;
  ch_m12[1] = 168;
  ch_m12[2] = 167;
  ch_m12[3] = 165;
  ch_m12[4] = 163;
  ch_m12[5] = 161;
  ch_m12[6] = 158;
  ch_m12[7] = 156;
  ch_m12[8] = 154;
  ch_m12[9] = 152;
  ch_m12[10] = 150;
  ch_m12[11] = 148;
  ch_m12[12] = 147;
  ch_m12[13] = 145;
  ch_m12[14] = 144;
  ch_m12[15] = 0;
  ch_m12[16] = 0;
  
  TGraph* ch_gr = new TGraph(17,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  // ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}


TGraph* set_lep_ch_tanBeta50(){

  double ch_m0[21];
  double ch_m12[21];

  ch_m0[0] = 200;
  ch_m0[1] = 250;
  ch_m0[2] = 300;
  ch_m0[3] = 350;
  ch_m0[4] = 400;
  ch_m0[5] = 450;
  ch_m0[6] = 500;
  ch_m0[7] = 550;
  ch_m0[8] = 600;
  ch_m0[9] = 650;
  ch_m0[10] = 700;
  ch_m0[11] = 750;
  ch_m0[12] = 800;
  ch_m0[13] =850;
  ch_m0[14] = 900;
  ch_m0[15] = 950;
  ch_m0[16] = 1000;
  ch_m0[17] = 1050;
  ch_m0[18] = 1100;
  ch_m0[19] = 1100;
  ch_m0[20] = 200;
 
  ch_m12[0] = 157;
  ch_m12[1] = 156;
  ch_m12[2] = 156;
  ch_m12[3] = 155;
  ch_m12[4] = 155;
  ch_m12[5] = 154;
  ch_m12[6] = 154;
  ch_m12[7] = 153;
  ch_m12[8] = 153;
  ch_m12[9] = 152;
  ch_m12[10] = 152;
  ch_m12[11] = 152;
  ch_m12[12] = 152;
  ch_m12[13] = 152;
  ch_m12[14] = 152;
  ch_m12[15] = 153;
  ch_m12[16] = 153;
  ch_m12[17] = 153;
  ch_m12[18] = 154;
  ch_m12[19] = 0;
  ch_m12[20] = 0;
  
  
  TGraph* ch_gr = new TGraph(21,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}




TGraph* set_lep_sl(Int_t tanBeta){

  // CMS SUSY Summer2010 implementation
  //  double sl_m0[] =  {0,  0, 30, 50, 60, 75, 80,90,100};
  //  double sl_m12[] = {0,245,240,220,200,150,100,50,0}; 
  
  //contour from D0 trilepton paper (PLB 680 (2009) 34-43)

  double *sl_m0 = 0;
  double *sl_m12 = 0;
  int n = 0;

  double sl_m0_3[] ={0,  0, 10, 20, 30, 40, 50, 60, 70, 77,88,95};
  double sl_m12_3[]={0,245,242,239,232,222,209,189,165,140,60,0};
  int n_3 = 12;

  double sl_m0_10[]={ 0,  0, 11, 20, 24, 49, 70, 82,88,90};
  double sl_m12_10[]={0,240,237,233,230,200,150,100,50,0};
  int n_10 = 10;

  if (tanBeta==3){
    sl_m0 = sl_m0_3;
    sl_m12 = sl_m12_3;
    n = n_3;
  }
  //CMS PTDR-II
  //* Selectron_R line mass=99, ISASUGRA7.69, A0=0, m_top=175, tan(beta]=10
  if (tanBeta==10 || tanBeta==50){
    sl_m0 = sl_m0_10;
    sl_m12 = sl_m12_10;
    n = n_10;
  }

  TGraph* lep_sl = new TGraph(n,sl_m0,sl_m12);

  lep_sl->SetFillColor(5);
  lep_sl->SetLineColor(5);
  lep_sl->SetFillStyle(1001);
  
  return lep_sl;
}


TGraph* set_tev_sg_cdf(Int_t tanBeta){

  tanBeta = tanBeta;

  //  double sg_m0[] =  {0,  0, 20, 50,100,150,200,250,300,350,400,450,500,550,600,600};
  //  double sg_m12[] = {0,160,169,170,160,155,150,122,116,112,110,106,105,100, 98,  0};
  //  int np=16;
  //New CHF from CDF plot in ICHEP2010 talk (E. Halkiadakis)
  double sg_m0[]= {0,  0, 30, 75,150,185,225,310,360,400,430,500,600,600};
  double sg_m12[]={0,162,168,170,160,150,130,120,109,108,100, 96, 95,  0};
  int np=14;

  TGraph* sg_gr = new TGraph(np,sg_m0,sg_m12);

  //  gStyle->SetHatchesLineWidth(3);

  sg_gr->SetFillColor(2);
  sg_gr->SetLineColor(2);
  //  sg_gr->SetLineWidth(3);
  sg_gr->SetFillStyle(1001); 

  return sg_gr;

}

TGraph* set_tev_sg_d0(Int_t tanBeta){

  tanBeta = tanBeta;

  //  double sgd_m0[] = {0, 0,  50, 100,150,200,250,300,350,400,450,500,550,600,600};
  //  double sgd_m12[] = {0,168,167,162,157,145,125,120,110,108,95, 94 ,94 ,93,0};
  //  int np=15;
  double sgd_m0[]= {0,  0, 30, 80,150,240,320,400,500,600,600,0};
  double sgd_m12[]={0,167,166,162,156,138,121,109,105,105,  0,0};
  int npd=12;

  TGraph* sgd_gr = new TGraph(npd,sgd_m0,sgd_m12);

  gStyle->SetHatchesLineWidth(3);

  sgd_gr->SetFillColor(kMagenta+3);
  sgd_gr->SetLineColor(kMagenta+3);
  sgd_gr->SetLineWidth(3);
  sgd_gr->SetFillStyle(3335);

  return sgd_gr;

}

TGraph* set_tev_tlp_cdf(Int_t tanBeta){
  double tlp1_m0[] = {   0, 20, 40, 60, 70, 80, 90, 80, 70, 60};
  double tlp1_m12[] = {170,185,200,215,220,215,210,190,175,160};
  TGraph* tlp1_gr = new TGraph(10,tlp1_m0,tlp1_m12);

  tlp1_gr->SetFillColor(4);
  tlp1_gr->SetLineColor(4);
  tlp1_gr->SetFillStyle(1001);

  return tlp1_gr;
}

TGraph* set_tev_tlp_d0(Int_t tanBeta){
  double tlp2_m0[] = {  70, 80, 90,100,105,110,120,130,140};
  double tlp2_m12[] = {160,172,184,196,205,195,185,173,160};
  TGraph* tlp2_gr = new TGraph(9,tlp2_m0,tlp2_m12);

  tlp2_gr->SetFillColor(4);
  tlp2_gr->SetFillStyle(1001); 

  return tlp2_gr;

}

//------------------------------
// new (Frederic)
//------------------------------

TGraph* set_tev_stau(Int_t tanBeta){

    double st_m0_tanBeta10[] =  {0,   10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 , 130,  147, 0, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 6500};
    double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518, 559, 600., 682., 750.,750, 842., 921., 999., 1076, 1152, 1228, 1304, 1378, 1453, 1527, 1600, 1673, 1746, 1818, 1890, 1962, 2034, 2105, 2175, 2246, 2316, 2386, 2456, 2526, 2595}; 


    double st_m0_tanBeta40[] = {240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880};
    double st_m12_tanBeta40[] = {186, 256, 329, 400, 470, 537, 603, 666, 727, 787, 845, 902, 958, 1013, 1067, 1121, 1174, 1226, 1278, 1330, 1381, 1431, 1481, 1531, 1581, 1630, 1679, 1728, 1779, 1825, 1874, 1920, 1971};

    TGraph* st_gr_tanBeta10 = new TGraph(15,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(10,st_m0_tanBeta40,st_m12_tanBeta40);

    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);


    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 40)return st_gr_tanBeta40;

    return 0;
}
//-----------------------
// new (Sanjay)
//-----------------------

TF1* constant_squark(int tanBeta,int i){
//---lines of constant gluino/squark. 
// Min squark mass from 1st and 2nd generations using fit for tanbeta = 10.


  double coef1[] = {2.67058e+04, 6.39642e+04, 1.16565e+05, 1.95737e+05, 2.86190e+05};
  double coef2[] = {1.98772e-01, 2.11242e-01, 2.17734e-01, 2.39535e-01, 2.39768e-01};
  double coef3[] = {2.67058e+04, 6.39641e+04, 1.16565e+05, 1.95736e+05, 2.86189e+05};

  
  char hname[200];

  sprintf(hname,"lnsq_%i",i); 
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]+[2])",0,1000);
  lnsq->SetParameter(0,coef1[i-1]);
  lnsq->SetParameter(1,coef2[i-1]);
  lnsq->SetParameter(2,coef3[i-1]);
  lnsq->SetLineWidth(1);
  lnsq->SetLineColor(kGray);

  return lnsq;
}

TGraph* set_NoEWSB(Int_t tanBeta){

    double st_m0_tanBeta10[]  = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta10[] = {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0}; 

//// Needs to be modified
    double st_m0_tanBeta40[] = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta40[]= {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0};

    TGraph* st_gr_tanBeta10 = new TGraph(25,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(25,st_m0_tanBeta40,st_m12_tanBeta40);

    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);

    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 40)return st_gr_tanBeta40;

    return 0;
}

//-------------------
// new (Sanjay)
//-------------------

TF1* constant_gluino(int tanBeta,int i){
//---lines of constant gluino/squark
  char hname[200];
  sprintf(hname,"lngl_%i",i); 

  double coef1[] = {201.77, 311.027, 431.582, 553.895, 676.137};
  double coef2[] = {-0.0146608, -0.01677, -0.022244, -0.0271851, -0.0292212};
    
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,2000);
  lngl->SetParameter(0,coef1[i-1]);
  lngl->SetParameter(1,coef2[i-1]);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kGray);

  return lngl;
}

//-------------------
// old (2010)
//-------------------
/*
TF1* constant_gluino(int tanBeta,int i){

  tanBeta = tanBeta;

  //---lines of constant gluino/squark
  
  //double coef1 = 0.35;
  //double coef2[] = {5,5,4.6,4.1};

  char hname[200];

  sprintf(hname,"lngl_%i",i); 
    
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,m0max);
  lngl->SetParameter(0,(500+spacing*(i-1))/2.4);
  lngl->SetParameter(1,-40./1400);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kGray);

  return lngl;
}
*/

TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta_){
  char legnm[200];

  //sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV/c^{2}}",500+150*(it-1));
  sprintf(legnm,"#font[92]{#tilde{q}(%i)}",500+(int)spacing*(it-1));
  Double_t place_x = 160;
  Double_t angle   = -8.;
  if(tanBeta_ == 50)            place_x = 290;
  if(tanBeta_ == 3 && it == 1 ){
    place_x = 220;
    angle = -15.;
  }

  double x = place_x+10*(it-1);
  double y = lnsq.Eval(place_x+10*(it-1))+5;

  if( it == 1 ){
    x     = 200;
    y     = 230;
    angle = -20;
  }
  else if( it == 2 ){
    x     = 225;
    y     = 355;
    angle = -15;
  }
  else if( it == 3 ){
    x     = 250;
    y     = 480;
    angle = -10;
  }

  TLatex* t3 = new TLatex(x,y,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(angle);
  t3->SetTextColor(kGray+2);

  
  return t3;
}

TLatex* constant_gluino_text(Int_t it,TF1& lngl){
  char legnm[200];

  //sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV/c^{2}}",500+150*(it-1));
  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)}",500+(int)spacing*(it-1));
  //TLatex* t4 = new TLatex(400,18+lngl.Eval(480),legnm);

  double x = 855;
  double y = 25+lngl.Eval(x+80);

  if( it == 1 ){
    x  = 855;
    y -= 35.;
  }

  TLatex* t4 = new TLatex(x,y,legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}



TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_){

  txtsz = txtsz;

  Double_t ypos_1 = 0.78;
  Double_t ypos_2 = 0.80;
  Double_t xpos_1 = 0.17;
  Double_t xpos_2 = 0.19;
  if(tanBeta_ == 50){
    xpos_1 = 0.17;
    xpos_2 = 0.18;
    ypos_1 = 0.76;
    ypos_2 = 0.78;

  }
  TLegend* legst = new TLegend(xpos_1,ypos_1,xpos_2,ypos_2);
  legst->SetHeader("#tilde{#tau} = LSP");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(85);

  return legst;
}

TLegend* makeNoEWSBLegend(Double_t txtsz,Int_t tanBeta_){

  txtsz = txtsz;

  Double_t ypos_1 = 0.10;
  Double_t ypos_2 = 0.20;
  Double_t xpos_1 = 0.82;
  Double_t xpos_2 = 0.92;
  if(tanBeta_ == 40){
    xpos_1 = 0.10;
    xpos_2 = 0.20;
    ypos_1 = 0.85;
    ypos_2 = 0.95;

  }

  TLegend* legst = new TLegend(xpos_1,ypos_1,xpos_2,ypos_2);
  legst->SetHeader("No EWSB");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(30);

  return legst;
}


TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph& tev_sn,Double_t txtsz,Int_t tanbeta){

  //TLegend* legexp = new TLegend(0.61,0.65,0.91,0.9,NULL,"brNDC");
  TLegend* legexp = new TLegend(0.57,0.65,0.91,0.9,NULL,"brNDC");

  legexp->SetFillColor(0);
  legexp->SetShadowColor(0);
  legexp->SetTextSize(txtsz);
  legexp->SetBorderSize(0);

  sg_gr.SetLineColor(1);
  
  //legexp->AddEntry(&sg_gr,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, 2 fb^{-1}}","f"); 
  legexp->AddEntry(&sg_gr,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0, 2 fb^{-1}}","f");   

  //  sgd_gr.SetLineColor(1);
  //  sgd_gr.SetLineWidth(1);

  //legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, 2.1 fb^{-1}}","f");  
  legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0, 2.1 fb^{-1}}","f");  

  ch_gr.SetLineColor(1);
  legexp->AddEntry(&ch_gr,"LEP2   #tilde{#chi}_{1}^{#pm}","f");  
  
  sl_gr.SetLineColor(1);
  if(tanbeta != 50) legexp->AddEntry(&sl_gr,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); 
  if(tanbeta == 3) legexp->AddEntry(&tev_sn,"D0  #chi^{#pm}_{1}, #chi^{0}_{2}","f");  
 

  return legexp;

}

