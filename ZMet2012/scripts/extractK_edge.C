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

void extractK_edge( bool exclusive = false , bool highMet = false , bool central = false , bool printplot = false );

void doPlots( bool printplot = false ){

  extractK_edge( false , false , false , printplot );
  extractK_edge( false , false , true  , printplot );
  extractK_edge( false , true  , false , printplot );
  extractK_edge( false , true  , true  , printplot );

  extractK_edge( true  , false , false , printplot );
  extractK_edge( true  , false , true  , printplot );
  extractK_edge( true  , true  , false , printplot );
  extractK_edge( true  , true  , true  , printplot );

  //extractK_edge( true  , printplot , false );
  //extractK_edge( false , printplot , true  );
  //extractK( false , printplot , false );

}

void extractK_edge( bool exclusive , bool highMet , bool central , bool printplot ){

  char* sigChar = "lowMet";
  if( highMet ) sigChar = "highMet";

  char* lepChar = "inclusiveLeptons";
  if( central ) lepChar = "centralLeptons";

  char* exChar = "inclusive";
  if( exclusive ) exChar = "exclusive";

  cout << "Signal region       : " << sigChar << endl;
  cout << "Leptons             : " << lepChar << endl;
  cout << "Inclusive/exclusive : " << exChar  << endl;

  char* iter = (char*) "V00-01-04";

  //char* suffix = "";
  char* suffix = "_2jets";

  TChain *data = new TChain("T1");
  data->Add(Form("../output/%s/data_ALL_53X_baby%s.root",iter,suffix));

  TChain *mc = new TChain("T1");

  mc->Add(Form("../output/%s/ttbar_53X_baby%s.root"       ,iter,suffix));
  mc->Add(Form("../output/%s/zjets_53X_baby%s.root"       ,iter,suffix));
  mc->Add(Form("../output/%s/ww_53X_baby%s.root"          ,iter,suffix));
  mc->Add(Form("../output/%s/wz_53X_baby%s.root"          ,iter,suffix));
  mc->Add(Form("../output/%s/zz_53X_baby%s.root"          ,iter,suffix));
  mc->Add(Form("../output/%s/t_53X_baby%s.root"           ,iter,suffix));
  mc->Add(Form("../output/%s/ttZ_53X_baby%s.root"         ,iter,suffix));
  mc->Add(Form("../output/%s/ttW_53X_baby%s.root"         ,iter,suffix));
  mc->Add(Form("../output/%s/VVV_53X_baby%s.root"         ,iter,suffix));
  //mc->Add(Form("../output/%s/zjets_full_53X_baby%s.root"  ,iter,suffix));

  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut filters("csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1");
  TCut runrange("isdata==0 || (run<197556 || run>198913)");
  TCut pt2020("lep1.pt()>20 && lep2.pt()>20");
  TCut njets2("njets40>=2");
  TCut njets3("njets40 >= 3");
  TCut ht100("ht40 >= 100.0");
  TCut eta14("abs(lep1.eta())<1.4 && abs(lep2.eta())<1.4");

  TCut Zmass("dilmass>81 && dilmass<101");


  TCut sel;
  sel += em;
  sel += filters;
  sel += runrange;
  sel += pt2020;
  
  if( highMet ){
    sel += njets2;
    sel += ht100;
  }

  else{
    sel += njets3;
  }

  if( central ){
    sel += eta14;
  }

  TCut weight("vtxweight * weight");

  cout << "Using selection : " << sel.GetTitle() << endl;
  cout << "Using weight    : " << weight.GetTitle() << endl;

  cout << "OF entries (total)  " << data->GetEntries(sel+em) << endl;
  cout << "OF entries (Z mass) " << data->GetEntries(sel+Zmass+em) << endl;
  cout << "K                   " << (float) data->GetEntries(sel+Zmass+em) / (float) data->GetEntries(sel+em) << endl;

  int mynbins = 6;

  const unsigned int nbins = mynbins;

  TCut metcuts[nbins];

  if( exclusive ){
    metcuts[0] = TCut("pfmet>0   && pfmet<30" );
    metcuts[1] = TCut("pfmet>30  && pfmet<60" );
    metcuts[2] = TCut("pfmet>60  && pfmet<100");
    metcuts[3] = TCut("pfmet>100 && pfmet<150");
    metcuts[4] = TCut("pfmet>150 && pfmet<300");
    metcuts[5] = TCut("pfmet>300");
  }else{
    metcuts[0] = TCut("pfmet>0");
    metcuts[1] = TCut("pfmet>30");
    metcuts[2] = TCut("pfmet>60");
    metcuts[3] = TCut("pfmet>100");
    metcuts[4] = TCut("pfmet>150");
    metcuts[5] = TCut("pfmet>300");
  }

  float K[nbins];
  float Kerr[nbins];

  float KMC[nbins];
  float KMCerr[nbins];

  TH1F* htot = new TH1F("htot","",1,0,1);
  TH1F* hZ   = new TH1F("hZ"  ,"",1,0,1);
  htot->Sumw2();
  hZ->Sumw2();

  for( int i = 0 ; i < nbins ; ++i ){
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

    float   totMC     = htot->GetBinContent(1);
    float   ZMC       = hZ->GetBinContent(1);
    float   totMCerr  = htot->GetBinError(1);
    float   ZMCerr    = hZ->GetBinError(1);

    KMC[i]    = ZMC    / totMC;
    //KMCerr[i] = ZMCerr / totMC;
    KMCerr[i] = KMC[i] * sqrt(pow( totMCerr/totMC , 2 ) + pow( ZMCerr/ZMC , 2 ) );

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

  float x[nbins];
  float xerr[nbins];
  
  for(int i = 0 ; i < nbins ; ++i ){
    x[i]    = i+1;
    xerr[i] = 0;
  }

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetGridy();

  TH2F* hdummy = new TH2F("hdummy","",nbins,0.5,nbins+0.5,10,0,0.2);
  hdummy->Draw();
  hdummy->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  hdummy->GetYaxis()->SetTitle("K");

  if( exclusive ){
    hdummy->GetXaxis()->SetBinLabel(1,"0-30");
    hdummy->GetXaxis()->SetBinLabel(2,"30-60");
    hdummy->GetXaxis()->SetBinLabel(3,"60-100");
    hdummy->GetXaxis()->SetBinLabel(4,"100-150");
    hdummy->GetXaxis()->SetBinLabel(5,"150-300");
    hdummy->GetXaxis()->SetBinLabel(6,">300");
  }
  else{
    hdummy->GetXaxis()->SetBinLabel(1,">0");
    hdummy->GetXaxis()->SetBinLabel(2,">30");
    hdummy->GetXaxis()->SetBinLabel(3,">60");
    hdummy->GetXaxis()->SetBinLabel(4,">100");
    hdummy->GetXaxis()->SetBinLabel(5,">150");
    hdummy->GetXaxis()->SetBinLabel(6,">300");
  }


  TGraphErrors *gr   = new TGraphErrors(nbins,x,K,xerr,Kerr);
  TGraphErrors *grMC = new TGraphErrors(nbins,x,KMC,xerr,KMCerr);
  gr->GetXaxis()->SetTitle("E_{T}^{miss} bin");
  gr->GetYaxis()->SetTitle("K");
  gr->SetMaximum(0.2);
  grMC->SetLineColor(2);
  grMC->SetMarkerColor(2);
  grMC->SetMarkerStyle(25);

  gr->Draw("sameP");
  grMC->Draw("sameP");

  TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->AddEntry(gr  ,"data","lp");
  leg->AddEntry(grMC,"MC"  ,"lp");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  if( printplot ){
    c1->Print( Form("../plots/extractK_edge_%s_%s_%s.pdf",sigChar,lepChar,exChar) );
  }

}
