#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TPostScript.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TLine.h"
#include "TBox.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TCut.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TEllipse.h"
#include "TPaveStats.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;

TH1F* getPullHist(TH1F* h1, TH1F* h2);
bool drawlog(char* prefix);
inline double fround(double n, unsigned d);
void setTDRStyle();
string getDataset( string rel, string sample );
float chi2prob( TH1F* h , float xmin , float xmax);
void formatHist(TH1F* h1, TH1F* h2, const char* sample, char* var);

int main(int argc, char* argv[]){

  if( argc != 5){
    cout << "USAGE: tcmetValidation SAMPLE RELEASE1 RELEASE2 DATA/MC" << endl;
    exit(0);
  }

  setTDRStyle();
  gStyle->SetOptStat(0);

  string sample = argv[1];
  string rel1   = argv[2];
  string rel2   = argv[3];
  string datamc = argv[4];

  bool isData = ( datamc == "DATA" ) ? true : false;

  string dataset1 = getDataset( rel1 , sample );
  string dataset2 = getDataset( rel2 , sample );

  //ps output file name
  //char* psFileName  = Form( "%s/tcmetValidation.ps" , rel2.c_str() );
  ofstream ofile( Form("%s/webpage/%s.html" , rel2.c_str() , sample.c_str() ) );

  //2 files to compare
  TFile *f1 = TFile::Open( Form( "%s/%s.root" , rel1.c_str() , sample.c_str() ) );
  TFile *f2 = TFile::Open( Form( "%s/%s.root" , rel2.c_str() , sample.c_str() ) );

  //in root file, path to tcmet histos
  char* tcmetpath = "DQMData/RecoMETV/MET_Global/tcMet/METTask";
  
  //list variables for comparison
  vector<char*> vars;
  vector<float> xmin;
  vector<float> xmax;
  vector<int>   rebin;
  vector<float> prob;

  //--------------------------------
  //event vars
  //-------------------------------- 
  vars.push_back("MET");                             xmin.push_back(0);    xmax.push_back(100);  rebin.push_back(1);
  vars.push_back("SumET");                           xmin.push_back(0);    xmax.push_back(2000); rebin.push_back(20);
  vars.push_back("MEx");                             xmin.push_back(-250); xmax.push_back(250);  rebin.push_back(10);
  vars.push_back("MEy");                             xmin.push_back(-250); xmax.push_back(250);  rebin.push_back(10);
  vars.push_back("METPhi");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  vars.push_back("METSig");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  vars.push_back("Nevents");                         xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);

  if( !isData ){
    vars.push_back("METPhiResolution_GenMETTrue");   xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
    vars.push_back("METResolution_GenMETTrue");      xmin.push_back(-100); xmax.push_back(100);  rebin.push_back(2);
    vars.push_back("METPhiResolution_GenMETCalo");   xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
    vars.push_back("METResolution_GenMETCalo");      xmin.push_back(-200); xmax.push_back(200);  rebin.push_back(5);
  }

  //vars.push_back("fracTracks");                    xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("dMET");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("dMETx");                         xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("dMETy");                         xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("dMUx");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("dMUy");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);

  //--------------------------------
  //muon vars
  //--------------------------------
  vars.push_back("muonPt");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  vars.push_back("muonEta");                         xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(5);
  vars.push_back("muonD0");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  vars.push_back("muonNhits");                       xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  vars.push_back("muonNormalizedChi2");              xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  vars.push_back("muonSAhits");                      xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);  
  vars.push_back("MExCorrection");                   xmin.push_back(-100); xmax.push_back(100);  rebin.push_back(5);
  vars.push_back("MEyCorrection");                   xmin.push_back(-100); xmax.push_back(100);  rebin.push_back(5);
  vars.push_back("CorrectionFlag");                  xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  vars.push_back("nMus");                            xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);

  //--------------------------------
  //electron vars
  //--------------------------------
  vars.push_back("electronEta");                     xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(5);
  vars.push_back("electronPt");                      xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  vars.push_back("electronHoverE");                  xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);

  //--------------------------------
  //track vars
  //--------------------------------
  vars.push_back("trackPt");                         xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(4);
  vars.push_back("trackEta");                        xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  vars.push_back("trackNhits");                      xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  vars.push_back("trackNormalizedChi2");             xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);

  //vars.push_back("trackD0PVTX");                     xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  //vars.push_back("trackDZPVTX");                     xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  //vars.push_back("trackD0BS");                       xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(2);
  //vars.push_back("trackPtErr");                    xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("nEls");                          xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("nMusAsPis");                     xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("trackQuality");                  xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);
  //vars.push_back("trackAlgo");                     xmin.push_back(-999); xmax.push_back(-999); rebin.push_back(1);



    
  //make TLegend
  TH1F* h1dummy = new TH1F("h1dummy","",1,0,1);
  h1dummy->SetLineColor(2);
  h1dummy->SetMarkerColor(2);
  h1dummy->SetMarkerSize(0);
  TH1F* h2dummy = new TH1F("h2dummy","",1,0,1);
  h2dummy->SetLineColor(4);
  h2dummy->SetMarkerColor(4);
  h2dummy->SetMarkerSize(0);
  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry( h1dummy , rel1.c_str() );
  leg->AddEntry( h2dummy , rel2.c_str() );
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  //declare canvases, histos, etc
  const int nvar = vars.size();
  TH1F    *h1[nvar];
  TH1F    *h2[nvar];
  TPad    *mainpad[nvar];
  TPad    *pullpad[nvar];
  TLine line;
  line.SetLineStyle(2);
  TLatex t;
  t.SetNDC();

  //TCanvas *canvas = new TCanvas("canvas","canvas",800,800);
  //canvas->Print(Form("%s[",psFileName));
  //canvas->Clear();
  
  TCanvas *canvas[nvar];

  //loop over variables
  for(int ivar=0;ivar<nvar;ivar++){
    
    cout << vars.at(ivar) << endl;
    //make canvas and pad
    //canvas->cd();
    canvas[ivar]  = new TCanvas(Form("%s_can",vars.at(ivar)) , Form("%s_can",vars.at(ivar)) , 800 , 800);
    canvas[ivar]  -> cd();
    mainpad[ivar] = new TPad(Form("%s_mainpad",vars.at(ivar)),Form("%s_mainpad",vars.at(ivar)),0,0,1,0.8);
    mainpad[ivar] -> Draw();
    mainpad[ivar] -> cd();
    
    //get histos
    h1[ivar] = (TH1F*) f1->Get(Form("%s_%s",tcmetpath,vars.at(ivar)));
    h2[ivar] = (TH1F*) f2->Get(Form("%s_%s",tcmetpath,vars.at(ivar)));

//     if( h1[ivar]->Integral() != h2[ivar]->Integral() ){
//       cout << vars.at(ivar) << " h1 " << h1[ivar]->Integral() << " h2 " << h2[ivar]->Integral() << endl;
//       h2[ivar]->Scale( h1[ivar]->Integral() / h2[ivar]->Integral() );
//     }
    if( h1[ivar]->GetEntries() != h2[ivar]->GetEntries() ){
      cout << vars.at(ivar) << " h1 " << h1[ivar]->GetEntries() << " h2 " << h2[ivar]->GetEntries() << endl;
      h2[ivar]->Scale( h1[ivar]->GetEntries() / h2[ivar]->GetEntries() );
    }


    //format and draw histos
    //if(drawlog(vars.at(ivar))) mainpad[ivar]->SetLogy(1); 
    if( rebin.at(ivar) > 1 ) {
      h1[ivar] -> Rebin( rebin.at(ivar) );
      h2[ivar] -> Rebin( rebin.at(ivar) );
    }
    if( xmin.at(ivar) != -999 ){
      h1[ivar]->GetXaxis()->SetRangeUser( xmin.at(ivar) , xmax.at(ivar) );
      if( strcmp( sample.c_str() , "ttbar" ) == 0 && strcmp( vars.at(ivar) , "MET" ) == 0 ){
        h1[ivar]->Rebin(2);
        h2[ivar]->Rebin(2);
        h1[ivar]->GetXaxis()->SetRangeUser( 0 , 300 );
      }
    }


    mainpad[ivar]->SetLogy(1); 
    h1[ivar] -> SetLineColor(2);
    h2[ivar] -> SetLineColor(4);
    h2[ivar] -> SetMarkerColor(4);
    h2[ivar] -> SetMarkerSize(0.5);
    h2[ivar] -> SetMarkerStyle(20);
    h1[ivar] -> Draw();
    h1[ivar] -> SetTitle("");
    h1[ivar] -> GetXaxis() -> SetTitle(vars.at(ivar));
    h2[ivar] -> Draw("sameE1");

    formatHist( h1[ivar] , h2[ivar] , sample.c_str() , vars.at(ivar) );
    //leg->Draw();
    
    //make canvas and pad
    //canvas->cd();
    canvas[ivar]  -> cd();
    pullpad[ivar] = new TPad(Form("%s_pullpad",vars.at(ivar)),Form("%s_pullpad",vars.at(ivar)),0,0.8,1,1.);
    pullpad[ivar] -> Draw();
    pullpad[ivar] -> cd();

    //format and draw pull hist
    TH1F* hpull = getPullHist(h1[ivar],h2[ivar]);
    hpull->Draw("E1");
    hpull->GetXaxis()->SetLabelSize(0);
    hpull->GetXaxis()->SetTitleSize(0);
    hpull->GetYaxis()->SetTitleSize(0.16);
    hpull->GetYaxis()->SetLabelSize(0.16);
    hpull->SetLineColor(4);
    hpull->SetMarkerColor(4);
    hpull->SetMarkerSize(0.5);
    hpull->SetMarkerStyle(20);
    //hpull->GetYaxis()->SetRangeUser(-4,4);
    hpull->GetYaxis()->SetTitle("Pull");
    hpull->SetTitle("");
    //hpull->GetYaxis()->SetNdivisions(5);
    prob.push_back( chi2prob( hpull , xmin.at(ivar) , xmax.at(ivar) ) );

    //draw guidelines
    if( xmin.at(ivar) != -999){
      line.SetLineStyle(1);
      line.DrawLine(xmin.at(ivar), 0, xmax.at(ivar), 0);
      line.SetLineStyle(2);
      //line.DrawLine(xmin.at(ivar), 1, xmax.at(ivar), 1);
      //line.DrawLine(xmin.at(ivar),-1, xmax.at(ivar),-1);
    }else{
      line.SetLineStyle(1);
      line.DrawLine( hpull->GetXaxis()->GetXmin() , 0 , hpull->GetXaxis()->GetXmax() , 0);
      line.SetLineStyle(2);
      //line.DrawLine( hpull->GetXaxis()->GetXmin() ,  1 , hpull->GetXaxis()->GetXmax() ,  1);
      //line.DrawLine( hpull->GetXaxis()->GetXmin() , -1 , hpull->GetXaxis()->GetXmax() , -1);
    }

    //canvas->Update();
    //canvas->Print(psFileName);
    canvas[ivar] -> Modified();
    canvas[ivar] -> Update();
    canvas[ivar] -> Print( Form( "%s/webpage/plots/%s_%s.gif" , rel2.c_str() , sample.c_str() , vars.at(ivar) ));
    
  }
  
  TCanvas *cprob = new TCanvas("cprob","",800,600);
  cprob->cd();
  gPad->SetLogx(1);
  gStyle->SetPadLeftMargin(0.2);
 
  //make summary histo
  TH1F* hprob = new TH1F("hprob","",vars.size(),0,vars.size());

  for( int ivar = 1 ; ivar <= vars.size() ; ivar++ ){
    hprob->SetBinContent(           vars.size() - ivar + 1 , prob.at(ivar - 1) );
    hprob->GetXaxis()->SetBinLabel( vars.size() - ivar + 1 , vars.at(ivar - 1) );
    //cout << vars.at(ivar-1) << setw(20) << prob.at(ivar-1) << endl;
  }

  hprob->GetXaxis()->SetLabelSize(0.029);
  hprob->GetYaxis()->SetLabelSize(0.035);
  hprob->SetBarWidth(0.8);
  hprob->SetBarOffset(0.1);
  hprob->SetFillColor(4);
  hprob->GetYaxis()->SetTitle("#chi^{2}/ndf Probability");
  hprob->GetYaxis()->SetTitleSize(0.04);
  hprob->GetYaxis()->SetNdivisions(7);
  hprob->GetXaxis()->SetLabelOffset(0.001);
  hprob->SetMaximum(1.5);
  hprob->SetMinimum(1e-10);
  TH1 *hprobcopy = hprob->DrawCopy("hbar0");

  cprob -> Print( Form( "%s/webpage/plots/%s_%s.gif" , rel2.c_str() , sample.c_str() , "summary" ));
  
  ofile << "<table><tr>" << endl;
  ofile << "<td><img SRC=" << "plots/" << sample << "_summary.gif> </td>" << endl;
  ofile << "</tr></table>" << endl;

  for( int ivar = 0 ; ivar < vars.size() ; ivar++ ){
    ofile << " <H3> " << "<FONT color=#ff0000>" << rel1 << "  :  " << dataset1 << " </FONT> " << " </H3>" << endl;
    ofile << " <H3> " << "<FONT color=#0000FF>" << rel2 << "  :  " << dataset2 << " </FONT> " << " </H3>" << endl;
    ofile << "<table><tr>" << endl;
    ofile << "<td><img SRC=" << "plots/" << sample << "_" << vars.at(ivar) << ".gif> </td>" << endl;
    ofile << "</tr></table>" << endl;
  }

  //canvas->Update();
  //canvas->Print(Form("%s]",psFileName));
 
  return 0;
}

float chi2prob( TH1F* h , float xmin , float xmax){

  int ndf    = 0;
  float chi2 = 0.;

  int ibinmin = 1;
  int ibinmax = h->GetXaxis()->GetNbins() + 1;

  if( xmin != -999){
    ibinmin = h->FindBin(xmin);
    ibinmax = h->FindBin(xmax);
  }

  for(int ibin = ibinmin ; ibin < ibinmax ; ibin++){
    //if( fabs( h->GetBinContent( ibin ) ) 1.e-10 ){
      ndf++;
      if( h->GetBinError(ibin) > 0 ){
        chi2 += pow( h->GetBinContent(ibin) / h->GetBinError(ibin) , 2);
      }else{
        chi2 += pow( h->GetBinContent(ibin) , 2);
      }

      //}
    //cout << ibin << " " << h->GetBinCenter(ibin) << " " << h->GetBinContent(ibin) << endl;
  }

  return TMath::Prob( chi2 , ndf );
}

// pull = ( h1-h2 ) / h1
TH1F* getPullHist(TH1F* h1, TH1F* h2){
  
  TH1F* hout = (TH1F*) h1->Clone(Form("%s_clone",h1->GetName()));
  
  for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){
  
    float val = h2->GetBinContent(ibin) - h1->GetBinContent(ibin);
    float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
    if(fabs(err) < 1.e-10)  err = sqrt(h2->GetBinContent(ibin) + h1->GetBinContent(ibin));
    
    float denom = fabs( h1->GetBinContent(ibin) ) > 0. ? h1->GetBinContent(ibin) : 1;

    hout -> SetBinContent( ibin, val / denom );
    hout -> SetBinError(   ibin, err / denom );
  }
  return hout;
}

// pull = (h1-h2)/sigma
// TH1F* getPullHist(TH1F* h1, TH1F* h2){
  
//   TH1F* hout = (TH1F*) h1->Clone(Form("%s_clone",h1->GetName()));
  
//   for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){
  
//     float val = h2->GetBinContent(ibin) - h1->GetBinContent(ibin);
//     float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
//     if(fabs(err) < 1.e-10)  err = sqrt(h2->GetBinContent(ibin) + h1->GetBinContent(ibin));
    
//     hout -> SetBinContent(ibin,fabs(err) > 0 ? val/err : val);
//     hout -> SetBinError(ibin,1);
//   }
//   return hout;
// }


bool drawlog(char* prefix){
  
  if(strcmp(prefix,"MET") == 0)      return true;
  
  return false; 
}

inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void setTDRStyle() {
  
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(800); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  //tdrStyle->SetPadTopMargin(0.1);
  //tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetPadLeftMargin(0.2);
  //tdrStyle->SetPadRightMargin(0.16);

  // For the Global title:

  //tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  //tdrStyle->SetTitleXOffset(0.9);
  //tdrStyle->SetTitleYOffset(1.25);
  tdrStyle->SetTitleXOffset(1);
  tdrStyle->SetTitleYOffset(1);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->SetPalette(1);
  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}



string getDataset( string rel , string sample ){
  
  ifstream ifile( Form( "%s/.%sdataset" , rel.c_str() , sample.c_str() ) );
  string dataset;

  if( ifile.is_open() ){
    string dataset;
    ifile >> dataset;
    return dataset;
  }else{
    return "ERROR CANNOT FIND DATASETNAME";
  }

}




void formatHist(TH1F* h1, TH1F* h2, const char* sample, char* var){

  TLine line;
  TLatex t;
  t.SetNDC();

  //line.SetLineColor(4);
  line.SetLineWidth(2);
  
  if( strcmp(var,"MET")==0 && strcmp(sample,"ttbar")!=0 ){
    float metcut = 30;
    if( strcmp(sample,"data") == 0 ) metcut=15;
    if( strcmp(sample,"qcd")  == 0 ) metcut=40;
    if( strcmp(sample,"zmm")  == 0 ) metcut=20;
    if( strcmp(sample,"zee")  == 0 ) metcut=20;
    
    line.DrawLine( metcut  , h1->GetMinimum() , metcut  , 2 * h1->GetMaximum() );

    stringstream s1;
    s1 << "N(met>" << metcut << " GeV) " << h1->Integral( h1->FindBin(metcut) , 1000000);
    stringstream s2;
    s2 << "N(met>" << metcut << " GeV) " << h2->Integral( h2->FindBin(metcut) , 1000000);

    t.SetTextColor(2);
    t.DrawLatex( 0.5 , 0.7 , s1.str().c_str() );
    t.SetTextColor(4);
    t.DrawLatex( 0.5 , 0.6 , s2.str().c_str() );
  }

  if( strcmp(var,"METResolution_GenMETTrue")==0 && strcmp(sample,"ttbar")==0 ){
    
    stringstream s1;
    s1 << "mean " << fround(h1->GetMean(1),1) << " RMS " << fround(h1->GetRMS(1),1);
    stringstream s2;
    s2 << "mean " << fround(h2->GetMean(1),1) << " RMS " << fround(h2->GetRMS(1),1);

    t.SetTextSize(0.04);
    t.SetTextColor(2);
    t.DrawLatex( 0.25 , 0.85 , s1.str().c_str() );
    t.SetTextColor(4);
    t.DrawLatex( 0.25 , 0.8  , s2.str().c_str() );

  }

}

