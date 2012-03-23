#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TCanvas.h"
#include "TPaveText.h"
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
#include "TFile.h"
#include "TColor.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>
#include "tdrstyle_SUSY.C"

void doPlot( TCanvas* can, TH1F* hist_VV , TH1F* hist_OF , TH1F* hist_Z , TH1F* hist_QCD , TH1F* hist_data , TH1F* hist_LM , bool residual = false , bool print = false );

TH1F* smoothHist( TH1F* hin , int n ){
						
  if( n%2 == 0 ) cout << "Error! n = " << n << " must be odd" << endl;

  TH1F* hout = (TH1F*) hin->Clone("hout");

  int diff     = (n-1) / 2;
  int firstbin = (n+1) / 2;
  int lastbin  = hin->GetXaxis()->GetNbins() - firstbin;

  for( int ibin = firstbin ; ibin <= lastbin ; ++ibin ){
									
    float ave = 0;

    for( int jbin = ibin - diff ; jbin <= ibin + diff ; ++jbin ){
      ave += hin->GetBinContent(jbin) / (float) n;
    }

    hout->SetBinContent(ibin,ave);

  }
  
  return hout;

}

void makeZPlot_2jets( bool print = false ){

  setTDRStyle();

  //------------------------------------------------------
  // read in data histogram
  //------------------------------------------------------

  TFile* file1     = new TFile("data_47fb_2jets.root");
  TH1F*  hist_data = (TH1F*) file1->Get("hdata");

  //------------------------------------------------------
  // read in histograms: photon pred, OF, WZ/ZZ
  //------------------------------------------------------

  TFile *file2 = new TFile("histsave_nj2_4p7fb.root");

  TH1F* hist_OF_2011A      = (TH1F*) file2->Get("ttbarpred_2011A_all");
  TH1F* hist_OF_2011B      = (TH1F*) file2->Get("ttbarpred_2011B_all");
  TH1F* hist_photon_2011A  = (TH1F*) file2->Get("phopred_2011A_all");
  TH1F* hist_photon_2011B  = (TH1F*) file2->Get("phopred_2011B_all");
  TH1F* hist_VV_2011A      = (TH1F*) file2->Get("dibosonpred_2011A_all");
  TH1F* hist_VV_2011B      = (TH1F*) file2->Get("dibosonpred_2011B_all");

  TH1F* hist_OF = (TH1F*) hist_OF_2011A->Clone("OF");
  hist_OF->Add(hist_OF_2011B);

  TH1F* hist_photon = (TH1F*) hist_photon_2011A->Clone("photon");
  hist_photon->Add(hist_photon_2011B);

  TH1F* hist_VV = (TH1F*) hist_VV_2011A->Clone("VV");
  hist_VV->Add(hist_VV_2011B);
  //hist_VV->Add(hist_OF,-1);

  //------------------------------------------------------
  // read in QCD histogram
  //------------------------------------------------------

  TFile* file3  = new TFile("data-10gev-nominal-1.2ht-6percent-final.root");
  TH1F*  hist_QCD_temp = (TH1F*) file3->Get("metEstimated2JetsIncl");
  hist_QCD_temp->Rebin(5);
  TH1F* hist_QCD = (TH1F*) hist_photon->Clone("QCD");
  for( unsigned int ibin = 1 ; ibin <= hist_photon->GetNbinsX() ; ibin++ )
    hist_QCD->SetBinContent( ibin , hist_QCD_temp->GetBinContent(ibin) );
  
  //TH1F*  hist_QCD = (TH1F*) hist_photon->Clone("QCD");

  //------------------------------------------------------
  // get LM4/LM8
  //------------------------------------------------------

  TFile* fileLM = new TFile("LM4_histos.root");
  TH1F* hist_LM4 = (TH1F*) fileLM->Get("LM4_2jets");

  //------------------------------------------------------
  // make the plot
  //------------------------------------------------------

  TCanvas *main_canvas = new TCanvas("main_canvas", "main_canvas",0,0,600,750);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  main_canvas->SetHighLightColor(2);
  main_canvas->Range(0,0,1,1);
  main_canvas->SetFillColor(0);
  main_canvas->SetBorderMode(0);
  main_canvas->SetBorderSize(0);
  main_canvas->SetTickx(1);
  main_canvas->SetTicky(1);
  main_canvas->SetLeftMargin(0.13);
  main_canvas->SetRightMargin(0.05);
  main_canvas->SetTopMargin(0.08);
  main_canvas->SetBottomMargin(0.13);
  main_canvas->SetFrameFillStyle(0);
  main_canvas->SetFrameBorderMode(0);

  doPlot( main_canvas , hist_VV , hist_OF , hist_photon , hist_QCD , hist_data , hist_LM4 , true , print );

}


void doPlot( TCanvas *can , TH1F* hist_VV , TH1F* hist_OF , TH1F* hist_photon , TH1F* hist_QCD , TH1F* hist_data , TH1F* hist_LM , bool residual , bool print ){

  //-----------------------------------------
  // dummy check
  //-----------------------------------------
  
  int bin = hist_VV->FindBin(100);
  cout << "data   : " << hist_data->Integral(bin,10000)   << endl;
  cout << "QCD    : " << hist_QCD->Integral(bin,10000)    << endl;
  cout << "gjets  : " << hist_photon->Integral(bin,10000) << endl;
  cout << "OF     : " << hist_OF->Integral(bin,10000)     << endl;
  cout << "VV     : " << hist_VV->Integral(bin,10000)     << endl;
  cout << "LM8    : " << hist_LM->Integral(bin,10000)     << endl;

  //-----------------------------------------
  // create a TPad for the main plot
  //-----------------------------------------

  can->cd(); 

  if( residual ){

    TPad *plotpad = new TPad("plotpad", "plotpad",0,0.2,1,1);
    plotpad->Draw();
    plotpad->cd();
    plotpad->SetLogy();
    plotpad->Range(0,0,1,1);
    plotpad->SetFillColor(0);
    plotpad->SetBorderMode(0);
    plotpad->SetBorderSize(0);
    plotpad->SetTickx(1);
    plotpad->SetTicky(1);
    plotpad->SetLeftMargin(0.13);
    plotpad->SetRightMargin(0.05);
    plotpad->SetTopMargin(0.08);
    plotpad->SetBottomMargin(0.13);
    plotpad->SetFrameFillStyle(0);
    plotpad->SetFrameBorderMode(0);

  }else{
    gPad->SetLogy();
  }

  //-----------------------------------------
  // histogram formatting
  //-----------------------------------------

  hist_data->SetMinimum(0.1);
  //hist_data->SetMaximum(29035);
  hist_data->SetLineStyle(0);
  hist_data->SetMarkerStyle(20);
  hist_data->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  hist_data->GetXaxis()->CenterTitle(true);
  hist_data->GetXaxis()->SetLabelFont(42);
  hist_data->GetXaxis()->SetLabelOffset(0.007);
  hist_data->GetXaxis()->SetTitleSize(0.05);
  hist_data->GetXaxis()->SetTitleFont(42);
  hist_data->GetYaxis()->SetTitle("events / 5 GeV");
  hist_data->GetYaxis()->CenterTitle(true);
  hist_data->GetYaxis()->SetLabelFont(42);
  hist_data->GetYaxis()->SetLabelOffset(0.007);
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetTitleFont(42);
  hist_data->GetZaxis()->SetLabelFont(42);
  hist_data->GetZaxis()->SetLabelOffset(0.007);
  hist_data->GetZaxis()->SetTitleSize(0.05);
  hist_data->GetZaxis()->SetTitleFont(42);
  hist_data->Draw("E1");

  hist_VV->SetFillColor(kGreen-10);
  //hist_OF->SetFillColor(kMagenta+2);
  hist_OF->SetFillColor(kMagenta-8);
  //hist_VV->SetFillStyle(3006);
  hist_photon->SetLineColor(kRed);
  hist_QCD->SetLineColor(kBlue);
  hist_photon->SetFillColor(0);
  hist_VV->SetLineColor(1);
  hist_OF->SetLineColor(1);
  hist_VV->SetLineWidth(1);
  hist_OF->SetLineWidth(1);
  hist_photon->SetLineWidth(2);
  hist_QCD->SetLineWidth(2);

  hist_LM->SetLineColor(kOrange+1);
  hist_LM->SetLineWidth(2);

  TH1F* hist_LM_smooth = smoothHist(hist_LM,5);
  //hist_LM_smooth->SetLineColor(1);

  //-----------------------------------------
  // make a stack of the prediction histos
  //-----------------------------------------

  THStack* mcstack = new THStack("mcstack","mcstack");
  mcstack->Add(hist_OF,"hist");
  mcstack->Add(hist_VV,"hist");
  hist_photon->Add(hist_OF);
  hist_photon->Add(hist_VV);

  //-----------------------------------------
  // make a hist of the total prediction
  // this is used only for the ratio plot
  //-----------------------------------------

  TH1F* totpred = (TH1F*) hist_photon->Clone("totpred");
  //totpred->Add(hist_QCD);
  //totpred->Scale(0.5);

  //-----------------------------------------
  // draw the histos
  //-----------------------------------------

  mcstack->Draw("same");
  hist_data->Draw("sameE1");
  hist_photon->Draw("samehist");
  //hist_QCD->Draw("samehist");
  //hist_LM->Draw("same");
  hist_LM_smooth->Draw("same");
  hist_data->Draw("sameaxis");
  hist_data->Draw("sameE1");


  //-----------------------------------------
  // make the legend
  //-----------------------------------------

  TLegend *leg = new TLegend(0.55,0.62,0.9,0.9);
  leg->AddEntry(hist_data,"data","lp");
  leg->AddEntry(hist_photon    ,"total bkg (#gamma+jets)","l");
  //leg->AddEntry(hist_QCD    ,"total bkg (QCD)","l");
  leg->AddEntry(hist_VV   ,"WZ/ZZ prediction","f");
  leg->AddEntry(hist_OF   ,"OF prediction","f");
  leg->AddEntry(hist_LM   ,"LM4","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  //-----------------------------------------
  // make the CMS preliminary text
  //-----------------------------------------

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.042);
  text->DrawLatex(0.60,0.48,"ee/#mu#mu + #geq2 jets");
  text->DrawLatex(0.14,0.95,"CMS Preliminary       #sqrt{s} = 7 TeV,   L_{int} = 4.7 fb^{-1}");

  //-----------------------------------------
  // make a TPad for the ratio histogram
  //-----------------------------------------

  if( residual ){

    //TPad* respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
    
    //plotpad->cd();
    //plotpad->Modified();
    can->cd();    

    // ------------>Primitives in pad: coverpad
    TPad* coverpad = new TPad("coverpad", "coverpad",0.122,0.16,1,0.303);
    coverpad->Draw();
    coverpad->cd();
    coverpad->Range(0,0,1,1);
    coverpad->SetFillColor(0);
    coverpad->SetBorderMode(0);
    coverpad->SetBorderSize(0);
    coverpad->SetTickx(1);
    coverpad->SetTicky(1);
    coverpad->SetLeftMargin(0.13);
    coverpad->SetRightMargin(0.05);
    coverpad->SetTopMargin(0.08);
    coverpad->SetBottomMargin(0.13);
    coverpad->SetFrameFillStyle(0);
    coverpad->SetFrameBorderMode(0);
    coverpad->Modified();
    can->cd();

    TPad* bottompad = new TPad("bottompad", "Ratio Pad",0,0,1,0.289);
    bottompad->Draw();
    bottompad->cd();
    bottompad->Range(-63.41463,-1,424.3902,2.333333);
    bottompad->SetFillColor(0);
    bottompad->SetBorderMode(0);
    bottompad->SetBorderSize(2);
    bottompad->SetTickx(1);
    bottompad->SetTicky(1);
    bottompad->SetLeftMargin(0.13);
    bottompad->SetRightMargin(0.05);
    bottompad->SetBottomMargin(0.3);
    bottompad->SetFrameFillStyle(0);
    bottompad->SetFrameBorderMode(0);
    bottompad->SetFrameFillStyle(0);
    bottompad->SetFrameBorderMode(0);


    //-----------------------------------------
    // make the histogram of syst uncertainty
    // in the background prediction
    //-----------------------------------------

    TH1F* systUp = (TH1F*) hist_data->Clone("systUp");
    TH1F* systDn = (TH1F*) hist_data->Clone("systDn");

    for( unsigned int ibin = 1 ; ibin <= hist_data->GetXaxis()->GetNbins() ; ibin++ ){
      if( totpred->GetBinContent(ibin) > 0 ){
	systUp->SetBinContent( ibin , 1 + totpred->GetBinError(ibin) / totpred->GetBinContent(ibin) );
	systDn->SetBinContent( ibin , 1 - totpred->GetBinError(ibin) / totpred->GetBinContent(ibin) );
      }else{
	systUp->SetBinContent( ibin , 0 );
	systDn->SetBinContent( ibin , 0 );
      }
    }

    //-----------------------------------------
    // make the ratio: observed/predicted
    // we set the errors to 0 because we want
    // the error bars to reflect the observed
    // statistical uncertainty only
    //-----------------------------------------

    for( unsigned int ibin = 1 ; ibin <= totpred->GetXaxis()->GetNbins() ; ibin++ )
      totpred->SetBinError(ibin,0);
   
    TH1F* ratio = (TH1F*) hist_data->Clone("data/pred");
    ratio->Divide(totpred);

    //-----------------------------------------
    // histogram formatting
    //-----------------------------------------

    ratio->GetYaxis()->SetRangeUser(0,2);
    ratio->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    ratio->GetXaxis()->CenterTitle(true);
    ratio->GetXaxis()->SetLabelFont(42);
    ratio->GetXaxis()->SetLabelOffset(0.007);
    ratio->GetXaxis()->SetLabelSize(0.1104);
    ratio->GetXaxis()->SetTitleSize(0.138);
    ratio->GetXaxis()->SetTitleFont(42);
    ratio->GetYaxis()->SetTitle("ratio");
    ratio->GetYaxis()->CenterTitle(true);
    ratio->GetYaxis()->SetNdivisions(-502);
    ratio->GetYaxis()->SetLabelFont(42);
    ratio->GetYaxis()->SetLabelOffset(0.007);
    ratio->GetYaxis()->SetLabelSize(0.1104);
    ratio->GetYaxis()->SetTitleSize(0.138);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetYaxis()->SetTitleFont(42);
    ratio->GetZaxis()->SetLabelFont(42);
    ratio->GetZaxis()->SetLabelOffset(0.007);
    ratio->GetZaxis()->SetTitleSize(0.05);
    ratio->GetZaxis()->SetTitleFont(42);


    //-----------------------------------------
    // draw histograms
    //-----------------------------------------

    systUp->SetFillColor(5);
    systDn->SetFillColor(5);

    ratio->Draw();
    //systUp->Draw("samehist");
    //systDn->Draw("samehist");

    TLine line;
    ratio->Draw("same");
    line.DrawLine(0,1,350,1);
    ratio->GetYaxis()->SetTitle("data/pred");
  }

  //-----------------------------------------
  // print file 
  //-----------------------------------------

  if(print ){
    can->Print("metplot_2jets.gif");
    can->Print("metplot_2jets.png");
    can->Print("metplot_2jets.pdf");
    can->Print("metplot_2jets.C");
    can->Print("metplot_2jets.ps");
    gROOT->ProcessLine(".! ps2pdf metplot_2jets.ps metplot_2jets_ppt.pdf");
  }

}
