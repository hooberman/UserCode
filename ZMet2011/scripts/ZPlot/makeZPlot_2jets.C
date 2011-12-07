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
#include "TFile.h"
#include "TColor.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>
#include "tdrstyle_SUSY.C"

void doPlot( TCanvas* can, TH1F* hist_VV , TH1F* hist_OF , TH1F* hist_Z , TH1F* hist_QCD , TH1F* hist_data , bool residual = false , bool print = false );

void makeZPlot_2jets( bool print = false ){

  setTDRStyle();

  //------------------------------------------------------
  // read in data histogram
  //------------------------------------------------------

  TFile* file1     = new TFile("data_35fb_2jets.root");
  TH1F*  hist_data = (TH1F*) file1->Get("hdata");

  //------------------------------------------------------
  // read in histograms: photon pred, OF, WZ/ZZ
  //------------------------------------------------------

  TFile *file2 = new TFile("histsave_2jets.root");

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
  // make the plot
  //------------------------------------------------------

  TCanvas *can = new TCanvas();
  doPlot( can , hist_VV , hist_OF , hist_photon , hist_QCD , hist_data , true , print );

}


void doPlot( TCanvas *can , TH1F* hist_VV , TH1F* hist_OF , TH1F* hist_photon , TH1F* hist_QCD , TH1F* hist_data , bool residual , bool print ){

  //-----------------------------------------
  // dummy check
  //-----------------------------------------
  
  int bin = hist_VV->FindBin(100);
  cout << "data   : " << hist_data->Integral(bin,10000)   << endl;
  cout << "QCD    : " << hist_QCD->Integral(bin,10000)    << endl;
  cout << "gjets  : " << hist_photon->Integral(bin,10000) << endl;
  cout << "OF     : " << hist_OF->Integral(bin,10000)     << endl;
  cout << "VV     : " << hist_VV->Integral(bin,10000)     << endl;

  //-----------------------------------------
  // create a TPad for the main plot
  //-----------------------------------------

  can->cd(); 

  if( residual ){
    TPad* plotpad = new TPad("plotpad","plotpad",0.0,0.0,1.0,0.8);
    plotpad->Draw();
    plotpad->cd();
    plotpad->SetLogy();
    plotpad->SetRightMargin(0.05);
  }else{
    gPad->SetLogy();
  }

  //-----------------------------------------
  // histogram formatting
  //-----------------------------------------

  hist_data->SetLineColor(1);
  hist_data->SetMarkerColor(1);
  hist_data->SetMarkerSize(1);
  hist_data->GetXaxis()->SetTitle("E_{T}^{miss} (GeV/c)");
  hist_data->SetMinimum(0.05);
  hist_data->Draw("E1");
  hist_data->SetMinimum(0.1);
  hist_VV->SetFillColor(kViolet-1);
  hist_OF->SetFillColor(kGreen-2);
  //hist_photon->SetFillColor(kRed-4);
  hist_photon->SetLineColor(kRed);
  hist_QCD->SetLineColor(kBlue);
  hist_photon->SetFillColor(0);
  hist_VV->SetLineColor(1);
  hist_OF->SetLineColor(1);
  //hist_photon->SetLineColor(1);
  hist_VV->SetLineWidth(1);
  hist_OF->SetLineWidth(1);
  hist_photon->SetLineWidth(2);
  hist_QCD->SetLineWidth(2);

  //-----------------------------------------
  // make a stack of the prediction histos
  //-----------------------------------------

  THStack* mcstack = new THStack("mcstack","mcstack");
  mcstack->Add(hist_OF,"hist");
  mcstack->Add(hist_VV,"hist");
  //mcstack->Add(hist_photon,"hist");
  hist_photon->Add(hist_OF);
  hist_photon->Add(hist_VV);
  //hist_QCD->Add(hist_OF);
  //hist_QCD->Add(hist_VV);

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
  hist_data->Draw("sameaxis");
  hist_data->Draw("sameE1");

  //-----------------------------------------
  // make the legend
  //-----------------------------------------

  TLegend *leg = new TLegend(0.55,0.45,0.95,0.7);
  leg->AddEntry(hist_data,"data","lp");
  leg->AddEntry(hist_photon    ,"total bkg (#gamma+jets)","l");
  //leg->AddEntry(hist_QCD    ,"total bkg (QCD)","l");
  leg->AddEntry(hist_VV   ,"WZ/ZZ prediction","f");
  leg->AddEntry(hist_OF   ,"OF prediction","f");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  //-----------------------------------------
  // make the CMS preliminary text
  //-----------------------------------------

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.05);
  text->DrawLatex(0.42,0.88,"CMS Preliminary");
  text->DrawLatex(0.42,0.81,"#sqrt{s} = 7 TeV, #scale[0.6]{#int}L dt = 3.5 fb^{-1}");
  text->DrawLatex(0.42,0.74,"Events with ee/#mu#mu + #geq2 jets");

  //-----------------------------------------
  // make a TPad for the ratio histogram
  //-----------------------------------------

  if( residual ){
    can->cd();
    TPad* respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
    respad->Draw();
    respad->cd();
    respad->SetRightMargin(0.05);

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
   
    TH1F* ratio = (TH1F*) hist_data->Clone("ratio");
    ratio->Divide(totpred);

    //-----------------------------------------
    // histogram formatting
    //-----------------------------------------

    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetYaxis()->SetLabelSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetTitleSize(0.24);
    ratio->GetYaxis()->SetTitleOffset(0.25);
    ratio->GetYaxis()->SetRangeUser(0,2);
    ratio->GetXaxis()->SetTitle("");

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
    //can->Print("metplot_2jets.ps");
    //gROOT->ProcessLine(".! ps2pdf metplot_2jets.ps");
  }

}
