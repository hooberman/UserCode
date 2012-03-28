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
#include "TMath.h"
#include <sstream>
#include <iomanip>


using namespace std;

void LMyields(){

  TChain *LM = new TChain("T1");
  //LM->Add("../output/V00-02-03/LM4_baby.root");
  //LM->Add("../output/V00-02-03/LM8_baby.root");
  //LM->Add("../output/V00-02-04/LM4v2_baby.root");
  LM->Add("../output/V00-02-04/LM8v2_baby.root");


  TCut sel("dilmass>81 && dilmass<101 && leptype<2");
  TCut weight("4.7 * weight");
  TCut ee("leptype==0");
  TCut mm("leptype==1");
  TCut em("leptype==2");
  TCut eeweight("1.00");
  TCut mmweight("0.90");
  TCut emweight("0.95");

  // TCut njets2("njets>=2");
  // TCut njetsup2("njetsup>=2");
  // TCut njetsdn2("njetsdn>=2");

  TCut njets2("njets>=3");
  TCut njetsup2("njetsup>=3");
  TCut njetsdn2("njetsdn>=3");

  vector<TCut> metcuts;

  //----------------------------------
  // inclusive bins
  //----------------------------------

  metcuts.push_back("pfmet>30");
  metcuts.push_back("pfmet>60");
  metcuts.push_back("pfmet>100");
  metcuts.push_back("pfmet>200");
  metcuts.push_back("pfmet>300");

  //----------------------------------
  // exclusive bins
  //----------------------------------

  // metcuts.push_back("pfmet>30 && pfmet<60");
  // metcuts.push_back("pfmet>60 && pfmet<100");
  // metcuts.push_back("pfmet>100 && pfmet<200");
  // metcuts.push_back("pfmet>200 && pfmet<300");
  // metcuts.push_back("pfmet>300");

  TH1F* hee = new TH1F("hee","",1,0,1);
  TH1F* hmm = new TH1F("hmm","",1,0,1);
  TH1F* hem = new TH1F("hem","",1,0,1);
  hee->Sumw2();
  hmm->Sumw2();
  hem->Sumw2();

  const unsigned int n = metcuts.size();

  float yield[n];
  float yieldsc[n];
  float yieldscerr[n];
  float yieldup[n];
  float yielddn[n];

  for( unsigned int i = 0 ; i < n ; ++i ){

    //--------------------------------
    // nominal
    //--------------------------------

    LM->Draw("0.5>>hee",(sel+metcuts.at(i)+njets2+ee)*weight*eeweight);
    LM->Draw("0.5>>hmm",(sel+metcuts.at(i)+njets2+mm)*weight*mmweight);

    float tot    = hee->GetBinContent(1) + hmm->GetBinContent(1);
    float toterr = sqrt( pow( hee->GetBinError(1) , 2 ) + pow( hmm->GetBinError(1) , 2 ) );

    //--------------------------------
    // OF contamination
    //--------------------------------

    LM->Draw("0.5>>hem",(metcuts.at(i)+njets2+em)*weight*emweight);

    float oftot = 0.16 * hem->GetBinContent(1);

    //--------------------------------
    // JES up
    //--------------------------------

    TString metupstring(metcuts.at(i).GetTitle());
    metupstring.ReplaceAll("pfmet","pfmetup");
    TCut metup(metupstring);

    LM->Draw("0.5>>hee",(sel+metup+ee+njetsup2)*weight*eeweight);
    LM->Draw("0.5>>hmm",(sel+metup+mm+njetsup2)*weight*mmweight);

    float totup    = hee->GetBinContent(1) + hmm->GetBinContent(1);
    float toterrup = sqrt( pow( hee->GetBinError(1) , 2 ) + pow( hmm->GetBinError(1) , 2 ) );

    //--------------------------------
    // JES up
    //--------------------------------

    TString metdnstring(metcuts.at(i).GetTitle());
    metdnstring.ReplaceAll("pfmet","pfmetdn");
    TCut metdn(metdnstring);

    LM->Draw("0.5>>hee",(sel+metdn+ee+njetsdn2)*weight*eeweight);
    LM->Draw("0.5>>hmm",(sel+metdn+mm+njetsdn2)*weight*mmweight);

    float totdn    = hee->GetBinContent(1) + hmm->GetBinContent(1);
    float toterrdn = sqrt( pow( hee->GetBinError(1) , 2 ) + pow( hmm->GetBinError(1) , 2 ) );

    float dup   = totup/tot-1;
    float ddn   = 1-totdn/tot;
    float djes  = 0.5 * (dup+ddn);
    

    cout << endl << endl;
    cout << "----------------------------------------------" << endl;
    cout << metcuts.at(i).GetTitle() << " "  << Form("%.1f +/- %.1f",tot,toterr) << endl;
    cout << metup.GetTitle()         << " "  << Form("%.1f +/- %.1f",totup,toterrup) << endl;
    cout << metdn.GetTitle()         << " "  << Form("%.1f +/- %.1f",totdn,toterrdn) << endl;
    cout << "JES uncertainty " << Form("%.2f",djes) << endl;

    float staterr = toterr/tot;
    float tot_uncertainty = sqrt( 0.06*0.06 + 0.05*0.05 + djes*djes + staterr*staterr);
    
    cout << "Yield    " << Form("%.1f $\\pm$ %.1f",tot,tot_uncertainty*tot) << endl;
    cout << "Yield SC " << Form("%.1f $\\pm$ %.1f",tot-oftot,tot_uncertainty*(tot-oftot)) << endl;

    yield[i] = tot;
    yieldup[i] = totup;
    yielddn[i] = totdn;

    yieldsc[i]    = tot-oftot;
    yieldscerr[i] = tot_uncertainty*(tot-oftot);
  }
  
  cout << "LM_yield[nbins] = { " << Form("%.1f",yield[0]) << " , " << Form("%.1f",yield[1]) << " , " << Form("%.1f",yield[2])
       << " , " << Form("%.1f",yield[3]) << " , " << Form("%.1f",yield[4]) << " }; " << endl; 
  cout << "LM_yield_JESup[nbins] = { " << Form("%.1f",yieldup[0]) << " , " << Form("%.1f",yieldup[1]) << " , " << Form("%.1f",yieldup[2])
       << " , " << Form("%.1f",yieldup[3]) << " , " << Form("%.1f",yieldup[4]) << " }; " << endl; 
  cout << "LM_yield_JESdown[nbins] = { " << Form("%.1f",yielddn[0]) << " , " << Form("%.1f",yielddn[1]) << " , " << Form("%.1f",yielddn[2])
       << " , " << Form("%.1f",yielddn[3]) << " , " << Form("%.1f",yielddn[4]) << " }; " << endl; 




  cout << endl << endl;

  cout << "LMX  & " 
       << Form("$%.0f \\pm %.1f$",yieldsc[0],yieldscerr[0]) << "  &  " 
       << Form("$%.0f \\pm %.1f$",yieldsc[1],yieldscerr[1]) << "  &  " 
       << Form("$%.0f \\pm %.1f$",yieldsc[2],yieldscerr[2]) << "  &  " 
       << Form("$%.0f \\pm %.1f$",yieldsc[3],yieldscerr[3]) << "  &  " 
       << Form("$%.0f \\pm %.1f$",yieldsc[4],yieldscerr[4]) << endl;


  TH1F* hsignom  = new TH1F("hsignom","",3,100,400);
  TH1F* hsigup   = new TH1F("hsigup" ,"",3,100,400);
  TH1F* hsigdn   = new TH1F("hsigdn" ,"",3,100,400);
  TH1F* hbkg     = new TH1F("hbkg"   ,"",3,100,400);

  hbkg->SetBinContent(1,213);
  hbkg->SetBinContent(2,14);
  hbkg->SetBinContent(3,2.3);
  hbkg->SetBinError(1,22);
  hbkg->SetBinError(2,2.5);
  hbkg->SetBinError(3,0.7);

  hsignom->SetBinContent(1,yield[2]);
  hsigup->SetBinContent (1,yieldup[2]);
  hsigdn->SetBinContent (1,yielddn[2]);
  hsignom->SetBinContent(2,yield[3]);
  hsigup->SetBinContent (2,yieldup[3]);
  hsigdn->SetBinContent (2,yielddn[3]);
  hsignom->SetBinContent(3,yield[4]);
  hsigup->SetBinContent (3,yieldup[4]);
  hsigdn->SetBinContent (3,yielddn[4]);

  TCanvas *can = new TCanvas();
  can->cd();
  gPad->SetLogy();

  hsignom->SetLineWidth(2);
  hsigup->SetLineColor(2);
  hsigup->SetLineStyle(2);
  hsigdn->SetLineColor(4);
  hsigdn->SetLineStyle(2);

  hbkg->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
  hbkg->Draw();
  hsignom->Draw("same");
  hsigup->Draw("same");
  hsigdn->Draw("same");

  TLegend *leg = new TLegend(0.4,0.7,0.9,0.9);
  leg->AddEntry(hbkg,"predicted bkg","lp");
  leg->AddEntry(hsignom,"expected signal (nominal)","l");
  leg->AddEntry(hsigup,"expected signal (JES +7.5%)","l");
  leg->AddEntry(hsigdn,"expected signal (JES -7.5%)","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  can->Print("LMshapes_LM8.pdf");
}
