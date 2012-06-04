#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

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

using namespace std;

void extractK( bool exclusive = false , bool printplot = false ){

  char* iter = (char*) "V00-00-12";

  TChain *data = new TChain("T1");
  data->Add(Form("../output/%s/data_baby.root",iter));

  TChain *mc = new TChain("T1");
  mc->Add(Form("../output/%s/ttbar_baby.root",iter));
  mc->Add(Form("../output/%s/zjets_baby.root",iter));

  TCut pfleptons("pflep1.pt() > 20 && pflep2.pt() > 20");
  TCut transveto("el1tv==0 && el2tv==0");
  TCut Zmass("dilmasspf>81 && dilmasspf<101");
  TCut njets2("njets>=2");
  TCut em("leptype==2");

  TCut sel;
  sel += njets2;
  sel += pfleptons;
  sel += transveto;
  sel += em;

  TCut weight("vtxweight * weight");

  cout << "OF entries (total)  " << data->GetEntries(sel+em) << endl;
  cout << "OF entries (Z mass) " << data->GetEntries(sel+Zmass+em) << endl;

  TCut metcuts[6];
  if( exclusive ){
    metcuts[0] = TCut("pfmet>0   && pfmet<30" );
    metcuts[1] = TCut("pfmet>30  && pfmet<60" );
    metcuts[2] = TCut("pfmet>60  && pfmet<100");
    metcuts[3] = TCut("pfmet>100 && pfmet<200");
    metcuts[4] = TCut("pfmet>200 && pfmet<300");
    metcuts[5] = TCut("pfmet>300");
  }else{
    metcuts[0] = TCut("pfmet>0");
    metcuts[1] = TCut("pfmet>30");
    metcuts[2] = TCut("pfmet>60");
    metcuts[3] = TCut("pfmet>100");
    metcuts[4] = TCut("pfmet>200");
    metcuts[5] = TCut("pfmet>300");
  }


  float K[6];
  float Kerr[6];

  float KMC[6];
  float KMCerr[6];

  TH1F* htot = new TH1F("htot","",1,0,1);
  TH1F* hZ   = new TH1F("hZ"  ,"",1,0,1);
  htot->Sumw2();
  hZ->Sumw2();

  for( int i = 0 ; i < 6 ; ++i ){
    float tot = (float) data->GetEntries(sel+em+metcuts[i]);
    float   Z = (float) data->GetEntries(sel+em+metcuts[i]+Zmass);
    
    K[i]    = Z/tot;
    Kerr[i] = sqrt(Z)/tot;

    // float totMC = (float) mc->GetEntries(sel+em+metcuts[i]);
    // float   ZMC = (float) mc->GetEntries(sel+em+metcuts[i]+Zmass);
    
    // KMC[i]    = ZMC/totMC;
    // KMCerr[i] = sqrt(ZMC)/totMC;

    mc->Draw("0.5>>htot" , (sel+em+metcuts[i]      ) * weight );
    mc->Draw("0.5>>hZ"   , (sel+em+metcuts[i]+Zmass) * weight );

    float totMC    = htot->GetBinContent(1);
    float   ZMC    = hZ->GetBinContent(1);
    float   ZMCerr = hZ->GetBinError(1);

    KMC[i]    = ZMC    / totMC;
    KMCerr[i] = ZMCerr / totMC;

    cout << endl;
    cout << "--------------------------------------------------------------" << endl;
    cout << metcuts[i].GetTitle() << endl << endl;
    cout << "data  : " << endl;
    cout << "total : " << tot << endl;
    cout << "Z     : " <<   Z << endl;
    cout << "K     : " << Form("%.2f +/- %.3f",K[i],Kerr[i]) << endl << endl;
    cout << "MC    : " << endl;
    cout << "total : " << totMC << endl;
    cout << "Z     : " <<   ZMC << endl;
    cout << "K     : " << Form("%.2f +/- %.3f",KMC[i],KMCerr[i]) << endl;
    cout << "--------------------------------------------------------------" << endl;
    cout << endl;
  }

  float x[6]    = {1,2,3,4,5,6};
  float xerr[6] = {0,0,0,0,0,0};

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",6,0.5,6.5,10,0,0.2);
  hdummy->Draw();
  hdummy->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  hdummy->GetYaxis()->SetTitle("K");
  if( exclusive ){
    hdummy->GetXaxis()->SetBinLabel(1,"0-30");
    hdummy->GetXaxis()->SetBinLabel(2,"30-60");
    hdummy->GetXaxis()->SetBinLabel(3,"60-100");
    hdummy->GetXaxis()->SetBinLabel(4,"100-200");
    hdummy->GetXaxis()->SetBinLabel(5,"200-300");
    hdummy->GetXaxis()->SetBinLabel(6,">300");
  }
  else{
    hdummy->GetXaxis()->SetBinLabel(1,">0");
    hdummy->GetXaxis()->SetBinLabel(2,">30");
    hdummy->GetXaxis()->SetBinLabel(3,">60");
    hdummy->GetXaxis()->SetBinLabel(4,">100");
    hdummy->GetXaxis()->SetBinLabel(5,">200");
    hdummy->GetXaxis()->SetBinLabel(6,">300");
  }

  TGraphErrors *gr   = new TGraphErrors(6,x,K,xerr,Kerr);
  TGraphErrors *grMC = new TGraphErrors(6,x,KMC,xerr,KMCerr);
  gr->GetXaxis()->SetTitle("E_{T}^{miss} bin");
  gr->GetYaxis()->SetTitle("K");
  gr->SetMaximum(0.2);
  grMC->SetLineColor(2);
  grMC->SetMarkerColor(2);
  grMC->SetMarkerStyle(25);

  gr->Draw("sameP");
  grMC->Draw("sameP");

  TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->AddEntry(gr,"data","lp");
  leg->AddEntry(grMC,"MC","lp");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  if( printplot ){
    if(exclusive ) c1->Print("../plots/extractK_exclusive.pdf");
    else           c1->Print("../plots/extractK_inclusive.pdf");
  }

}
