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

int main(int argc, char* argv[]){

  if( argc != 4){
    cout << "USAGE: tcmetValidation SAMPLE RELEASE1 RELEASE2" << endl;
    exit(0);
  }

  setTDRStyle();
  gStyle->SetOptStat(0);

  string sample = argv[1];
  string rel1   = argv[2];
  string rel2   = argv[3];

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
  vars.push_back("MET");
  vars.push_back("SumET");
  vars.push_back("MEx");
  vars.push_back("MEy");
  //vars.push_back("dMET");
  //vars.push_back("dMETx");
  //vars.push_back("dMETy");
  //vars.push_back("dMUx");
  //vars.push_back("dMUy");
  vars.push_back("CorrectionFlag");
  vars.push_back("METPhi");
  vars.push_back("METPhiResolution_GenMETTrue");
  vars.push_back("METResolution_GenMETTrue");
  vars.push_back("Nevents");
  vars.push_back("electronHoverE");
  //vars.push_back("fracTracks");
  vars.push_back("muonEta");
  vars.push_back("muonNormalizedChi2");
  //vars.push_back("muonSAhits");
  //vars.push_back("nMus");
  //vars.push_back("trackAlgo");
  vars.push_back("trackEta");
  vars.push_back("trackNormalizedChi2");
  //vars.push_back("trackPtErr");
  vars.push_back("METPhiResolution_GenMETCalo");
  vars.push_back("METResolution_GenMETCalo");
  vars.push_back("METSig");
  vars.push_back("MExCorrection");
  vars.push_back("MEyCorrection");
  vars.push_back("electronEta");
  vars.push_back("electronPt");
  vars.push_back("muonD0");
  vars.push_back("muonNhits");
  vars.push_back("muonPt");
  //vars.push_back("nEls");
  //vars.push_back("nMusAsPis");
  vars.push_back("trackD0");
  vars.push_back("trackNhits");
  vars.push_back("trackPt");
  //vars.push_back("trackQuality");
  

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

    h1[ivar] -> Rebin(5);
    h2[ivar] -> Rebin(5);

    //format and draw histos
    //if(drawlog(vars.at(ivar))) mainpad[ivar]->SetLogy(1); 
    mainpad[ivar]->SetLogy(1); 
    h1[ivar] -> SetLineColor(2);
    h2[ivar] -> SetLineColor(4);
    h1[ivar] -> Draw();
    h1[ivar] -> SetTitle("");
    h1[ivar] -> GetXaxis() -> SetTitle(vars.at(ivar));
    h2[ivar] -> Draw("same");
    leg->Draw();
    
    //make canvas and pad
    //canvas->cd();
    canvas[ivar]  -> cd();
    pullpad[ivar] = new TPad(Form("%s_pullpad",vars.at(ivar)),Form("%s_pullpad",vars.at(ivar)),0,0.8,1,1.);
    pullpad[ivar] -> Draw();
    pullpad[ivar] -> cd();

    //format and draw pull hist
    TH1F* hpull = getPullHist(h1[ivar],h2[ivar]);
    hpull->Draw();
    hpull->GetXaxis()->SetLabelSize(0);
    hpull->GetXaxis()->SetTitleSize(0);
    hpull->GetYaxis()->SetTitleSize(0.16);
    hpull->GetYaxis()->SetLabelSize(0.16);
    hpull->GetYaxis()->SetRangeUser(-4,4);
    hpull->GetYaxis()->SetTitle("Pull");
    hpull->SetTitle("");
    hpull->GetYaxis()->SetNdivisions(5);

    //draw guidelines
    line.DrawLine(hpull->GetXaxis()->GetXmin(), 0, hpull->GetXaxis()->GetXmax(), 0);
    line.DrawLine(hpull->GetXaxis()->GetXmin(), 1, hpull->GetXaxis()->GetXmax(), 1);
    line.DrawLine(hpull->GetXaxis()->GetXmin(),-1, hpull->GetXaxis()->GetXmax(),-1);

    //canvas->Update();
    //canvas->Print(psFileName);
    canvas[ivar] -> Modified();
    canvas[ivar] -> Update();
    canvas[ivar] -> Print( Form( "%s/webpage/plots/%s_%s.gif" , rel2.c_str() , sample.c_str() , vars.at(ivar) ));

    //ofile << "<H2> Selector : " << aSel->GetName() << " : </H2>" << endl;
    //ofile << "<H3> Kaon Fake Rate : </H3>" << endl;
    //ofile << "<font color="#ff0000">" << dataset1 << "</font>" << endl; 
    ofile << "<H3> " << rel1 << "  :  " << dataset1 << " </H3>" << endl;
    ofile << "<H3> " << rel2 << "  :  " << dataset2 << " </H3>" << endl;
    ofile << "<table><tr>" << endl;
    ofile << "<td><img SRC=" << "plots/" << sample << "_" << vars.at(ivar) << ".gif> </td>" << endl;
    ofile << "</tr></table>" << endl;
    
  }
  
  //canvas->Update();
  //canvas->Print(Form("%s]",psFileName));
 

  return 0;
}

TH1F* getPullHist(TH1F* h1, TH1F* h2){
  
  TH1F* hout = (TH1F*) h1->Clone(Form("%s_clone",h1->GetName()));
  
  for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){
  
    float val = h2->GetBinContent(ibin) - h1->GetBinContent(ibin);
    float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
    if(fabs(err) < 1.e-10)  err = sqrt(h2->GetBinContent(ibin) + h1->GetBinContent(ibin));
    
    hout -> SetBinContent(ibin,fabs(err) > 0 ? val/err : val);
    //hout -> SetBinError(ibin,1);
  }
  
  return hout;
}


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
  //tdrStyle->SetPadLeftMargin(0.13);
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




