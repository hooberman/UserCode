#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include <iomanip>
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "THStack.h"
#include "TProfile.h"
#include <sstream>

bool  drawpull   = true;
bool  printtext  = false; 
float metval1    = 30;
float metval2    = 60;
float xt         = 0.5;
float yt         = 0.65;

using namespace std;

float calculateHistError( TH1F* h , int minbin , int maxbin );
TH1F* getCloneHist(TH1F* hin, int color);
TH1F* getPullHist(TH1F* h1, TH1F* h2);
inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void compareMET( bool printgif = false ){

  char* iter = "V01-02";

  int colors[5]={2,5,7,3,7};
 
  const float maxmet = 150;
  
  vector<char*> mcfilenames;
  vector<char*> mcleg;
  vector<int>   rebin;
 
  
  //ZJets MC with PhotonJet templates
  //TFile *f = TFile::Open(Form("output/%s/ZJets_PhotonJetTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate_njets1.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate_reweight_njets1.root",iter));
  TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate_reweight.root",iter));
  string sample = "Z+jets MC";
  string filename="ZJets";
  

  /*
  //LM4 MC with PhotonJet templates
  TFile *f = TFile::Open(Form("output/%s/babylooper_LM4_PhotonJetTemplate.root",iter));
  string sample = "LM4 MC";
  string filename="LM4";
  */
  
  //dilepton data with EG templates
  //TFile *f = TFile::Open(Form("output/%s/lepdata_EGTemplate.root",iter));
  //mcfilenames.push_back(Form("output/%s/TTbar_PhotonJetTemplate.root",iter));
  //mcfilenames.push_back(Form("output/%s/ZJets_PhotonJetTemplate.root",iter));
  
  /*
  TFile *f = TFile::Open(Form("output/%s/babylooper_lepdata_EGTemplate.root",iter));
  mcfilenames.push_back(Form("output/%s/babylooper_TTbar_PhotonJetTemplate.root",iter));
  mcfilenames.push_back(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));

  mcleg.push_back("t#bar{t}");
  mcleg.push_back("Z+jets");
  string sample = "Z+jets DATA";
  string filename = "lep";
  */

  vector<string> predhist;
  vector<string> obshist;
  vector<string> xtitle;
  vector<string> title;
  
  //histos to make
  //predhist.push_back("metPredicted");            obshist.push_back("metObserved");           title.push_back("");   rebin.push_back(5);
  //predhist.push_back("metPredicted_sf");         obshist.push_back("metObserved_sf");        title.push_back("(SF)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_df");         obshist.push_back("metObserved_df");        title.push_back("(DF)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_ptlt40");     obshist.push_back("metObserved_ptlt40");    title.push_back("(p_{T} < 40 GeV)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_pt40_60");    obshist.push_back("metObserved_pt40_60");   title.push_back("(40 < p_{T} < 60 GeV)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_ptgt60");     obshist.push_back("metObserved_ptgt60");    title.push_back("(p_{T} > 60 GeV)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_ptlt50");     obshist.push_back("metObserved_ptlt50");    title.push_back("(p_{T} < 50 GeV)"); rebin.push_back(5); 
  //predhist.push_back("metPredicted_ptgt50");     obshist.push_back("metObserved_ptgt50");    title.push_back("(p_{T} > 50 GeV)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_ee");         obshist.push_back("metObserved_ee");        title.push_back("(ee)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_mm");         obshist.push_back("metObserved_mm");        title.push_back("(#mu#mu)"); rebin.push_back(5);
  predhist.push_back("metPredicted_njets0");     obshist.push_back("metObserved_njets0");    title.push_back("(nJets = 1)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_njets1");     obshist.push_back("metObserved_njets1");    title.push_back("(nJets = 2)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_njets2");     obshist.push_back("metObserved_njets2");    title.push_back("(nJets #geq 3)"); rebin.push_back(5);
  //predhist.push_back("metPredicted_njets4");     obshist.push_back("metObserved_njets4");    title.push_back("(nJets=4)"); rebin.push_back(20);
  //predhist.push_back("metParPredicted");         obshist.push_back("metParObserved");        title.push_back(""); rebin.push_back(5);
  //predhist.push_back("metPerpPredicted");        obshist.push_back("metPerpObserved");       title.push_back(""); rebin.push_back(5);

  const unsigned int nhist    = predhist.size();

  const unsigned int nMC = mcfilenames.size();
  TH1F* metPredicted[nhist];
  TH1F* metObserved[nhist];
  TH1F* metObserved_MC[nhist][nMC];
  TFile* mcfiles[nMC];
  for( unsigned int iMC = 0 ; iMC < nMC ; ++iMC ){
    cout << "Opening MC file " << mcfilenames.at(iMC) << endl;
    mcfiles[iMC] = TFile::Open( mcfilenames.at(iMC) );
  }
  THStack* MCStack[nhist];

  TCanvas *can[nhist];
  TPad *pullpad[nhist];
  TPad *mainpad[nhist];

  TLine line;

  for( unsigned int i = 0 ; i < nhist ; i++ ){


    metPredicted[i] = (TH1F*) f->Get( predhist[i].c_str() );
    metObserved[i]  = (TH1F*) f->Get( obshist[i].c_str() );

    for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
      metObserved_MC[i][iMC] = (TH1F*) mcfiles[iMC]->Get( obshist[i].c_str() );
    }

    if( metPredicted[i] == 0 ){
      cout << "ERROR CAN'T FIND " << predhist[i] << endl;
      exit(0);
    }
    if( metObserved[i] == 0 ){
      cout << "ERROR CAN'T FIND " << obshist[i] << endl;
      exit(0);
    }

    can[i]=new TCanvas(Form("can_%i",i),"",800,600);
    can[i]->cd();

    if(drawpull)    mainpad[i] = new TPad(Form("mainpad_%i",i),Form("mainpad_%i",i),0,0,1,0.6);
    else            mainpad[i] = new TPad(Form("mainpad_%i",i),Form("mainpad_%i",i),0,0,1,1);

    mainpad[i] -> Draw();
    mainpad[i] -> cd();
    mainpad[i]->SetLogy(1);
    
    metPredicted[i]->Rebin( rebin.at(i) );
    metObserved[i]->Rebin( rebin.at(i) );
    //metPredicted[i]->Scale( metObserved[i]->Integral() / metPredicted[i]->Integral() );

    //if( i == 1 ) metPredicted[i]->GetXaxis()->SetRangeUser(-maxmet,maxmet);
    //else         metPredicted[i]->GetXaxis()->SetRangeUser(0,maxmet);
    
    metPredicted[i]->GetXaxis()->SetRangeUser(0,maxmet);
    metObserved[i]->GetXaxis()->SetRangeUser(0,maxmet);
   
    //metObserved[i]->GetXaxis()->SetTitle( xtitle[i].c_str() );
    metObserved[i]->GetXaxis()->SetTitle( "tcmet (GeV)" );
    metObserved[i]->SetTitle( Form("%s %s",sample.c_str(),title[i].c_str()) );
    metPredicted[i]->SetLineColor(4);
    metPredicted[i]->SetLineWidth(2);
    metObserved[i]->SetLineWidth(2);
    metPredicted[i]->SetMarkerColor(4);
    metPredicted[i]->SetMarkerSize(0);
    metObserved[i]->SetLineColor(1);
    metObserved[i]->SetMarkerColor(1);
    metObserved[i]->SetMinimum(0.002);
    metObserved[i]->SetMarkerSize(1);

    float max = metPredicted[i]->GetMaximum();
    if( metObserved[i]->GetMaximum() > max ) max = metObserved[i]->GetMaximum();
    metPredicted[i]->SetMaximum( 1.5 * max );
    
    //metPredicted[i]->Draw("hist");
    metObserved[i]->Draw("E1");
    MCStack[i] = new THStack(Form("stack_%i",i),Form("stack_%i",i));

    float mcyield1[nMC];
    float mcyield2[nMC];

    int bin1  = metPredicted[i]->FindBin(metval1);
    int bin2  = metPredicted[i]->FindBin(metval2);
    int maxbin = metPredicted[i]->GetXaxis()->GetNbins() + 1;


    for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
      metObserved_MC[i][iMC]->Rebin( rebin.at(i) );
      //metObserved_MC[i][iMC]->SetLineColor( colors[iMC] );
      mcyield1[iMC]= metObserved_MC[i][iMC]->Integral(bin1,maxbin);
      mcyield2[iMC]= metObserved_MC[i][iMC]->Integral(bin2,maxbin);
      metObserved_MC[i][iMC]->SetFillColor( colors[iMC] );
    
      MCStack[i]->Add( metObserved_MC[i][iMC] );
    }

    MCStack[i]->Draw("samehist");
    metPredicted[i]->Draw("samehist");
    metObserved[i]->Draw("sameE1");
    metObserved[i]->Draw("axissame");
  
    float npred1    = metPredicted[i]->Integral( bin1 , maxbin );
    float nprederr1 = calculateHistError( metPredicted[i] , bin1 , maxbin );

    float npred2    = metPredicted[i]->Integral( bin2 , maxbin );
    float nprederr2 = calculateHistError( metPredicted[i] , bin2 , maxbin );

    float nobs1    = metObserved[i]->Integral( bin1 , maxbin );
    float nobserr1 = calculateHistError( metObserved[i] , bin1 , maxbin );

    float nobs2  = metObserved[i]->Integral( bin2 , maxbin );
    float nobserr2 = calculateHistError( metObserved[i] , bin2 , maxbin );    

    TH1F* met_df = (TH1F*) f->Get("metObserved_df");
    met_df->Rebin( rebin.at(i) );
    float nem1   = met_df->Integral(bin1,maxbin);
    float nem2   = met_df->Integral(bin2,maxbin);

    cout << endl;
    cout << "|" << setw(15) << ""                << setw(2)
         << "|" << setw(15) << "N(met>30) GeV"   << setw(2)
         << "|" << setw(15) << "N(met>60) GeV"   << setw(2) << "|" << endl;
    cout << "|" << setw(15) << "data"            << setw(2)
         << "|" << setw(15) << nobs1             << setw(2)
         << "|" << setw(15) << nobs2             << setw(2) << "|" << endl;
    cout << "|" << setw(15) << "pred"            << setw(2)
         << "|" << setw(15) << fround(npred1,4)  << setw(2)
         << "|" << setw(15) << fround(npred2,4)  << setw(2) << "|" << endl;
    for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
      cout << "|" << setw(15) << mcleg.at(iMC)       << setw(2)
           << "|" << setw(15) << fround(mcyield1[iMC],2)  << setw(2)
           << "|" << setw(15) << fround(mcyield2[iMC],2)  << setw(2) << "|" << endl;
    }
    cout << "|" << setw(15) << "emu"             << setw(2)
         << "|" << setw(15) << nem1              << setw(2)
         << "|" << setw(15) << nem2              << setw(2) << "|" << endl;
    cout << endl;

    cout << "Predicted N(met > " << metval1 << ") " << npred1 << " +/- " << nprederr1 << endl;
    cout << "Observed  N(met > " << metval1 << ") " << nobs1  << " +/- " << nobserr1  << endl;
    cout << "Predicted N(met > " << metval2 << ") " << npred2 << " +/- " << nprederr2 << endl;
    cout << "Observed  N(met > " << metval2 << ") " << nobs2  << " +/- " << nobserr2  << endl;
    

    stringstream s1;
    stringstream s2;
    stringstream s3;
    stringstream s4;
    //s1 << "N(met > " << metval1 << ") " << fround(npred1,2) << endl; //" #pm " << fround(nprederr1,2) << endl;
    //s2 << "N(met > " << metval1 << ") " << fround(nobs1,2)  << endl; //" #pm " << fround(nobserr1,2)  << endl;
    //s3 << "N(met > " << metval2 << ") " << fround(npred2,2) << endl; //" #pm " << fround(nprederr2,2) << endl;
    //s4 << "N(met > " << metval2 << ") " << fround(nobs2,2)  << endl; //" #pm " << fround(nobserr2,2)  << endl;
    s1 << "N(met > " << metval1 << ") " << fround(npred1,2) << " #pm " << fround(nprederr1,2) << endl;
    s2 << "N(met > " << metval1 << ") " << fround(nobs1,2)  << " #pm " << fround(nobserr1,2)  << endl;
    s3 << "N(met > " << metval2 << ") " << fround(npred2,2) << " #pm " << fround(nprederr2,2) << endl;
    s4 << "N(met > " << metval2 << ") " << fround(nobs2,2)  << " #pm " << fround(nobserr2,2)  << endl;



    if( printtext ){
      TLatex * t = new TLatex();
      t->SetTextSize(0.04);
      t->SetNDC();
      t->SetTextColor(1);
      t->DrawLatex(xt,yt+0.15,s2.str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(xt,yt+0.1,s1.str().c_str());
      t->SetTextColor(1);
      t->DrawLatex(xt,yt+0.05,s4.str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(xt,yt,s3.str().c_str());
    }

    TLegend *leg = new TLegend(0.75,0.65,0.9,0.95);
    leg->AddEntry(metObserved[i],"Data","p");  
    leg->AddEntry(metPredicted[i],"Predicted");
    for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
      leg->AddEntry(metObserved_MC[i][iMC],mcleg[iMC],"f");
    }
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    leg->Draw();

    if( drawpull ){

      can[i]  -> cd();
      pullpad[i] = new TPad(Form("pullpad_%i",i),Form("pullpad_%i",i),0,0.6,1,1.);
      pullpad[i] -> Draw();
      pullpad[i] -> cd();
      
      //format and draw pull hist
      TH1F* hpull = getPullHist( metPredicted[i] , metObserved[i] );
      hpull->Draw("E1");
      hpull->GetXaxis()->SetLabelSize(0);
      hpull->GetXaxis()->SetTitleSize(0);
      hpull->GetYaxis()->SetTitleSize(0.075);
      hpull->GetYaxis()->SetLabelSize(0.075);
      hpull->SetLineColor(1);
      hpull->SetMarkerColor(1);
      hpull->SetMarkerSize(1);
      hpull->SetMarkerStyle(20);
      hpull->GetYaxis()->SetRangeUser(-1,1);
      hpull->GetYaxis()->SetTitle("(data-pred)/pred");
      hpull->GetYaxis()->SetTitleOffset(0.35);
      hpull->GetYaxis()->SetTitleSize(0.12);
      hpull->SetTitle("");
      hpull->GetYaxis()->SetNdivisions(5);
      //TF1* fpol = new TF1("fpol","pol1",0,50);
      //hpull->Fit(fpol,"R");
      
      
      line.SetLineStyle(1);
      line.DrawLine(0, 0, maxmet, 0);
      line.SetLineStyle(2);
      //line.DrawLine(0, 1, maxmet, 1);
      //line.DrawLine(0,-1, maxmet,-1);
    }
    
    if( printgif ) can[i]->Print(Form("plots/%s_%s.gif",filename.c_str(),predhist[i].c_str()));
  }

}

float calculateHistError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}


TH1F* getCloneHist(TH1F* hin, int color){

  TH1F* hout=new TH1F(hin->GetName(),hin->GetTitle(),
		      hin->GetNbinsX(),hin->GetXaxis()->GetXmin(),hin->GetXaxis()->GetXmax());

  for(int ibin=1;ibin<=hin->GetNbinsX();ibin++)
    hout->SetBinContent(ibin,hin->GetBinContent(ibin));

    hout->SetLineColor(1);
    hout->SetMarkerColor(color);
    hout->SetFillColor(color);

  return hout;
}


// pull = ( h1-h2 ) / h1
TH1F* getPullHist(TH1F* h1, TH1F* h2){
  
  TH1F* hout = (TH1F*) h1->Clone(Form("%s_clone",h1->GetName()));
  
  for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){
  
    float val = h2->GetBinContent(ibin) - h1->GetBinContent(ibin);
    float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
    if(fabs(err) < 1.e-10)  err = sqrt(h2->GetBinContent(ibin) + h1->GetBinContent(ibin));
    
    //float denom = fabs( h1->GetBinContent(ibin) ) > 0. ? h1->GetBinContent(ibin) : 1;
    float denom = h1->GetBinContent(ibin);

    //cout << "bin " << h1->GetBinCenter(ibin) << " h1 " << h1->GetBinContent(ibin) << " h2 " << h2->GetBinContent(ibin) << endl;
    if( h1->GetBinContent(ibin) > 0 && h2->GetBinContent(ibin) > 0 ){
      hout -> SetBinContent( ibin, val / denom );
      hout -> SetBinError(   ibin, err / denom );
    }
    else{
      //cout << "NO ENTRIES " << endl;
      hout -> SetBinContent( ibin , -9999 );
      hout -> SetBinError( ibin , 1 );
    }
  }

  return hout;
}

// TH1F* getPullHist(TH1F* h1, TH1F* h2){
  
//   TH1F* hout = (TH1F*) h1->Clone(Form("%s_clone",h1->GetName()));
  
//   for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){
  
//     float val = h2->GetBinContent(ibin) - h1->GetBinContent(ibin);
//     float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
//     //if(fabs(err) < 1.e-10)  err = sqrt(h2->GetBinContent(ibin) + h1->GetBinContent(ibin));
    
//     //hout -> SetBinContent(ibin,fabs(err) > 0 ? val/err : val);
//     if( fabs( val ) > 1.e-10 ){
//       hout -> SetBinContent(ibin, val/err);
//       hout -> SetBinError(ibin,1);
//     }
//   }
//   return hout;
// }
