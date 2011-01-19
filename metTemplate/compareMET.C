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
#include "Tools/histtools.cc"

//-------------------------------------------------------

//physics analysis params
const float metval1    = 30;
const float metval2    = 60;
const float metval3    = 120;
const float R          = 1.111;
const bool  useK       = true;
const float K          = 0.156;
const bool  addttbar   = true;
const char* iter       = "v5";

//plotting/text params
const int   colors[5]           = {9,8,5,2,4};
const bool  drawpull            = true;
const bool  printtext           = false; 
const bool  drawIntegralDist  	= false; 
const float xt         		= 0.5;
const float yt         		= 0.65;
const float maxmet     		= -1;
const int   nprec      		= 3;
const int   width1     		= 20;
const int   width2     		= 5;     

//-------------------------------------------------------

const bool makeLatexPlot = false;        //plot in latex style
char* pm         = " +/- ";
char* delim      = "|";
char* delimstart = "|";
char* delimend   = "|";
char* ee         = "ee";
char* mm         = "mm";
char* em         = "em";

//-------------------------------------------------------

void printLine(){

  if( makeLatexPlot ){
    cout << "\\hline" << endl;
  }
  else{
    cout << "---------------------------------------------------"
         << "--------------------------------------------------" << endl;
  }
}

//-------------------------------------------------------

float addquad(float x,float y){ return sqrt(x*x+y*y); }

//-------------------------------------------------------

using namespace std;

float calculateHistError( TH1F* h , int minbin , int maxbin );
TH1F* getCloneHist(TH1F* hin, int color);
TH1F* getPullHist(TH1F* h1, TH1F* h2);
inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

void compareMET( bool printgif = false ){

  if( makeLatexPlot ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }
  

   
  vector<char*> mcfilenames;
  vector<char*> mcleg;
  vector<int>   rebin;


 
  /*
  //ZJets MC with PhotonJet templates
  //TFile *f = TFile::Open(Form("output/%s/ZJets_PhotonJetTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate_njets1.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate_reweight_njets1.root",iter));
  TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_ZJets_PhotonJetTemplate_reweight.root",iter));
  string sample = "Z+jets MC";
  string filename="ZJets";
  */

  /*
  //LM4 MC with PhotonJet templates
  TFile *f = TFile::Open(Form("output/%s/babylooper_LM4_PhotonJetTemplate.root",iter));
  string sample = "LM4 MC";
  string filename="LM4";
  */
  
  //-------------------------------
  //dilepton data with EG templates
  //-------------------------------

  TFile *f = TFile::Open(Form("output/%s/babylooper_lepdata_skim_EGStitchedTemplate.root",iter));
  //TFile *f = TFile::Open(Form("output/%s/babylooper_lepdata_skim_EGTemplate.root",iter));  
  //TFile *f = TFile::Open(Form("output/%s/babylooper_lepdata_skim_victorTemplate.root",iter));  
  //TFile *f = TFile::Open(Form("output/%s/babylooper_lepdata_victorTemplate.root",iter)); 
  //mcfilenames.push_back(Form("output/%s/babylooper_WW_PhotonJetTemplate.root",iter));        mcleg.push_back("WW");
  //mcfilenames.push_back(Form("output/%s/babylooper_WZ_PhotonJetTemplate.root",iter));        mcleg.push_back("WZ");
  //mcfilenames.push_back(Form("output/%s/babylooper_ZZ_PhotonJetTemplate.root",iter));        mcleg.push_back("ZZ");

  //-------------------------------  
  //mc files
  //-------------------------------

  mcfilenames.push_back(Form("output/%s/babylooper_VV_PhotonJetTemplate.root",iter));        mcleg.push_back("VV");
  //mcfilenames.push_back(Form("output/%s/babylooper_tW_PhotonJetTemplate.root",iter));        mcleg.push_back("single top"); 
  mcfilenames.push_back(Form("output/%s/babylooper_TTbar_PhotonJetTemplate.root",iter));     mcleg.push_back("t#bar{t}");
  mcfilenames.push_back(Form("output/%s/babylooper_ZJets_PhotonJetTemplate.root",iter));     mcleg.push_back("Z+jets");
  //mcfilenames.push_back(Form("output/%s/babylooper_WJets_PhotonJetTemplate.root",iter));     mcleg.push_back("W+jets");
  
  //---------------------------------
  //Read in em met hist from ttbar MC
  //---------------------------------

  TFile* fttmc           = TFile::Open(Form("output/%s/babylooper_TTbar_PhotonJetTemplate.root",iter));
  TH1F*  httmc           = (TH1F*) fttmc->Get("metObserved_df");

  //TH1F*  httdata         = (TH1F*) f->Get("metObserved_df");
  //TH1F*  httdata_nozveto = (TH1F*) f->Get("metObserved_df_nozveto");
  //float nem = httdata->GetEntries();
  //cout << "nem " << nem << endl;
  //httmc->Scale( nem / httmc->Integral() );
  //httmc->Rebin(5);
  //httmc->SetLineColor(2);

  string sample = "";
  //string sample = "Z+jets DATA";
  string filename = "lep";
  
  assert( mcfilenames.size() == mcleg.size() );

  vector<string> predhist;
  vector<string> obshist;
  vector<string> xtitle;
  vector<string> title;

  //--------------
  //histos to make
  //--------------

  predhist.push_back("metPredicted");        obshist.push_back("metObserved");      title.push_back("(ee+#mu#mu)");  rebin.push_back(5);
  predhist.push_back("metPredicted_ee");     obshist.push_back("metObserved_ee");   title.push_back("(ee)");         rebin.push_back(5);
  predhist.push_back("metPredicted_mm");     obshist.push_back("metObserved_mm");   title.push_back("(#mu#mu)");     rebin.push_back(5);
  
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
  THStack* MCStack_integral[nhist];

  TCanvas* can[nhist];
  TPad*    pullpad[nhist];
  TPad*    mainpad[nhist];

  TLine line;
  
  for( unsigned int i = 0 ; i < nhist ; i++ ){

    //-------------------------------------
    //get predicted and observed met histos
    //-------------------------------------

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

    //-------------------------------------
    //draw predicted and observed met histos
    //-------------------------------------

    can[i]=new TCanvas(Form("can_%i",i),"",800,600);
    can[i]->cd();
  
    if(drawpull)    mainpad[i] = new TPad(Form("mainpad_%i",i),Form("mainpad_%i",i),0,0,1,0.6);
    else            mainpad[i] = new TPad(Form("mainpad_%i",i),Form("mainpad_%i",i),0,0,1,1);

    mainpad[i] -> Draw();
    mainpad[i] -> cd();
    mainpad[i]->SetLogy(1);
    
    metPredicted[i]->Rebin( rebin.at(i) );
    metObserved[i]->Rebin( rebin.at(i) );
    
    if( maxmet > 0 ){
      metPredicted[i]->GetXaxis()->SetRangeUser(0,maxmet);
      metObserved[i]->GetXaxis()->SetRangeUser(0,maxmet);
    }

    metObserved[i]->GetXaxis()->SetTitle( "pfmet (GeV)" );
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
    MCStack[i]          = new THStack(Form("stack_%i",i),Form("stack_%i",i));
    MCStack_integral[i] = new THStack(Form("stackintegral_%i",i),Form("stackintegral_%i",i));

    int bin1   = metPredicted[i]->FindBin(metval1);
    int bin2   = metPredicted[i]->FindBin(metval2);
    int bin3   = metPredicted[i]->FindBin(metval3);
    int maxbin = metPredicted[i]->GetXaxis()->GetNbins() + 1;
  
    //-----------------------------------
    //get MC yields, save in MCStack
    //-----------------------------------

    float mcyield1[nMC];
    float mcyield2[nMC];
    float mcyield3[nMC];
  
    for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
      metObserved_MC[i][iMC]->Rebin( rebin.at(i) );
      //metObserved_MC[i][iMC]->SetLineColor( colors[iMC] );
      mcyield1[iMC]= metObserved_MC[i][iMC]->Integral(bin1,maxbin);
      mcyield2[iMC]= metObserved_MC[i][iMC]->Integral(bin2,maxbin);
      mcyield3[iMC]= metObserved_MC[i][iMC]->Integral(bin3,maxbin);
      metObserved_MC[i][iMC]->SetFillColor( colors[iMC] );
      MCStack[i]->Add( metObserved_MC[i][iMC] );
      MCStack_integral[i]->Add( (TH1F*) cumulate( *metObserved_MC[i][iMC] , false).Clone() );
    }

    //-----------------------------------
    //get yields predicted from templates
    //-----------------------------------

    float npred1    = metPredicted[i]->Integral( bin1 , maxbin );
    float nprederr1 = calculateHistError( metPredicted[i] , bin1 , maxbin );

    float npred2    = metPredicted[i]->Integral( bin2 , maxbin );
    float nprederr2 = calculateHistError( metPredicted[i] , bin2 , maxbin );

    float npred3    = metPredicted[i]->Integral( bin3 , maxbin );
    float nprederr3 = calculateHistError( metPredicted[i] , bin3 , maxbin );

    //-----------------------------------
    //get observed yields
    //-----------------------------------

    float nobs1    = metObserved[i]->Integral( bin1 , maxbin );
    float nobserr1 = calculateHistError( metObserved[i] , bin1 , maxbin );

    float nobs2  = metObserved[i]->Integral( bin2 , maxbin );
    float nobserr2 = calculateHistError( metObserved[i] , bin2 , maxbin );    

    float nobs3  = metObserved[i]->Integral( bin3 , maxbin );
    float nobserr3 = calculateHistError( metObserved[i] , bin3 , maxbin ); 

    //---------------------------
    //get yields from em met hist
    //---------------------------

    char* dfhistname = "metObserved_df";
    if( useK ) dfhistname = "metObserved_df_nozveto";

    TH1F* met_df = (TH1F*) f->Get( dfhistname )->Clone("metObserved_df_clone");
    met_df->Rebin( rebin.at(i) );
    float nem1   = met_df->Integral(bin1,maxbin);
    float nem2   = met_df->Integral(bin2,maxbin);
    float nem3   = met_df->Integral(bin3,maxbin);

    float nemerr1 = sqrt( nem1 );
    float nemerr2 = sqrt( nem2 );
    float nemerr3 = sqrt( nem3 );

    float scale = 1;

    //------------------
    //scale using R, K
    //------------------

    if( TString(predhist[i]).Contains("ee") ){
      if( useK ) scale = K/(2*R);      
      else       scale = 1/(2*R);
      cout << "Scaling ee by: " << scale << endl;
    }
    else if( TString(predhist[i]).Contains("mm") ){
      if( useK ) scale = (R*K)/2;      
      else       scale = R/2;
      cout << "Scaling mm by: " << scale << endl;
    }
    else{
      if( useK ) scale = K * ( R/2 +  1/(2*R) );
      else       scale = R/2 +  1/(2*R);
      cout << "Scaling ee+mm by: " << scale << endl;
    }

    nem1    *= scale;
    nem2    *= scale;
    nem3    *= scale;
    nemerr1 *= scale;
    nemerr2 *= scale;
    nemerr3 *= scale;

    TH1F* httmc_clone = (TH1F*) httmc->Clone();
    httmc_clone->Rebin( rebin.at(i) );

    httmc_clone->SetLineColor(2);
    httmc_clone->SetMarkerColor(2);
    httmc_clone->SetLineWidth(2);

    float ntt = httmc_clone->Integral( bin2 , maxbin );
    httmc_clone->Scale( nem2 / ntt ); 

    cout << "MC ttbar yield > 60 GeV " << httmc_clone->Integral( bin2 , maxbin ) << endl;
    //------------------
    //print output table
    //------------------

    cout << endl;
    printLine();
    cout << delimstart << setw(width1) << ""                << setw(width2)
         << delim      << setw(width1) << "N(met>30)  GeV"  << setw(width2)
         << delim      << setw(width1) << "N(met>60)  GeV"  << setw(width2) 
         << delim      << setw(width1) << "N(met>120) GeV"  << setw(width2) << delimend << endl;
     
    cout << delimstart << setw(width1) << "Z pred"                                  << setw(width2)
         << delim      << setw(width1) << Form("%.2f %s %.2f",npred1,pm,nprederr1)  << setw(width2)
         << delim      << setw(width1) << Form("%.2f %s %.2f",npred2,pm,nprederr2)  << setw(width2) 
         << delim      << setw(width1) << Form("%.2f %s %.2f",npred3,pm,nprederr3)  << setw(width2) << delimend << endl;
    
    cout << delimstart << setw(width1) << "OFOS"                                << setw(width2)
         << delim      << setw(width1) << Form("%.2f %s %.2f",nem1,pm,nemerr1)  << setw(width2)
         << delim      << setw(width1) << Form("%.2f %s %.2f",nem2,pm,nemerr2)  << setw(width2) 
         << delim      << setw(width1) << Form("%.2f %s %.2f",nem3,pm,nemerr3)  << setw(width2) << delimend << endl;
    
    printLine();
    cout << delimstart << setw(width1) << "Z pred + OFOS"                                                 << setw(width2)
         << delim      << setw(width1) << Form("%.2f %s %.2f",nem1+npred1,pm,addquad(nemerr1,nprederr1))  << setw(width2)
         << delim      << setw(width1) << Form("%.2f %s %.2f",nem2+npred2,pm,addquad(nemerr2,nprederr2))  << setw(width2) 
         << delim      << setw(width1) << Form("%.2f %s %.2f",nem3+npred3,pm,addquad(nemerr3,nprederr3))  << setw(width2) << delimend << endl;
    
    printLine();
    cout << delimstart << setw(width1) << "data"            << setw(width2)
         << delim      << setw(width1) << nobs1             << setw(width2)
         << delim      << setw(width1) << nobs2             << setw(width2) 
         << delim      << setw(width1) << nobs3             << setw(width2) << delimend << endl;
   
    printLine();
 
    cout << endl;
    
    stringstream s1;
    stringstream s2;
    stringstream s3;
    stringstream s4;
    stringstream s5;
    stringstream s6;
  
    s1 << "N(met > " << metval1 << ") " << fround(npred1,nprec) << endl; //" #pm " << fround(nprederr1,2) << endl;
    s2 << "N(met > " << metval1 << ") " << fround(nobs1,nprec)  << endl; //" #pm " << fround(nobserr1,2)  << endl;
    s3 << "N(met > " << metval2 << ") " << fround(npred2,nprec) << endl; //" #pm " << fround(nprederr2,2) << endl;
    s4 << "N(met > " << metval2 << ") " << fround(nobs2,nprec)  << endl; //" #pm " << fround(nobserr2,2)  << endl;
    s5 << "N(met > " << metval3 << ") " << fround(npred3,nprec) << endl; //" #pm " << fround(nprederr2,2) << endl;
    s6 << "N(met > " << metval3 << ") " << fround(nobs3,nprec)  << endl; //" #pm " << fround(nobserr2,2)  << endl;
    
    if( addttbar ) metPredicted[i]->Add(httmc_clone);

    if( drawIntegralDist ){
      bool cumulateAscending = false;

      TH1F *metPredicted_integral    = (TH1F*) cumulate(*metPredicted[i] , cumulateAscending).Clone();
      TH1F *metObserved_integral     = (TH1F*) cumulate(*metObserved[i]  , cumulateAscending).Clone();
      TH1F *httmc_integral           = (TH1F*) cumulate(*httmc_clone     , cumulateAscending).Clone();

      MCStack_integral[i]->Draw("hist");
      metPredicted_integral->Draw("samehist");
      if( addttbar ) httmc_integral->Draw("samehist");
      metObserved_integral->Draw("sameE1");
      metObserved_integral->Draw("axissame");
    }else{
      MCStack[i]->Draw("samehist");
      metPredicted[i]->Draw("samehist");
      if( addttbar ) httmc_clone->Draw("same");
      //metPredicted[i]->Draw("sameE1");
      metObserved[i]->Draw("sameE1");
      metObserved[i]->Draw("axissame");
    }

    if( printtext ){
      TLatex * t = new TLatex();
      t->SetTextSize(0.04);
      t->SetNDC();
      
      t->SetTextColor(1);
      t->DrawLatex(xt,yt+0.15,s2.str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(xt,yt+0.10,s1.str().c_str());
      
      t->SetTextColor(1);
      t->DrawLatex(xt,yt+0.05,s4.str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(xt,yt,     s3.str().c_str());
      
      t->SetTextColor(1);
      t->DrawLatex(xt,yt-0.05,s6.str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(xt,yt-0.10,s5.str().c_str());
   
    }
  
    TLegend *leg = new TLegend(0.65,0.45,0.95,0.9);
    leg->AddEntry(metObserved[i],"Observed","p");  
    leg->AddEntry(metPredicted[i],"Z pred + OFOS");
    if(addttbar) leg->AddEntry(httmc_clone,    "OFOS","l");
    for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
      leg->AddEntry(metObserved_MC[i][iMC],mcleg[iMC],"f");
    }
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    leg->Draw();

    TLatex *t = new TLatex();
    t->SetNDC();
    t->DrawLatex(0.35,0.85, "CMS");
    t->DrawLatex(0.35,0.8,  "34.0 pb^{-1} at #sqrt{s} = 7 TeV");

    if( TString(predhist[i]).Contains("ee") )
      t->DrawLatex(0.35,0.75, "Events with ee");
    else if( TString(predhist[i]).Contains("mm") )
      t->DrawLatex(0.35,0.75, "Events with #mu#mu");
    else
      t->DrawLatex(0.35,0.75, "Events with ee/#mu#mu");


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
      hpull->GetYaxis()->SetTitle("(data-Z pred)/(Z pred)");
      hpull->GetYaxis()->SetTitleOffset(0.35);
      hpull->GetYaxis()->SetTitleSize(0.09);
      hpull->SetTitle("");
      hpull->GetYaxis()->SetNdivisions(5);
      //TF1* fpol = new TF1("fpol","pol1",0,50);
      //hpull->Fit(fpol,"R");
      
      line.SetLineWidth(3);
      line.SetLineStyle(1);
      line.DrawLine(0, 0, maxmet, 0);
      line.SetLineStyle(2);
      gPad->SetGridy(1);
      //line.DrawLine(0, 1, maxmet, 1);
      //line.DrawLine(0,-1, maxmet,-1);
    }
    
    if( printgif ) can[i]->Print(Form("plots/%s_%s.png",filename.c_str(),predhist[i].c_str()));
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




//     s1 << "N(met > " << metval1 << ") " << fround(npred1,2) << " #pm " << fround(nprederr1,2) << endl;
//     s2 << "N(met > " << metval1 << ") " << fround(nobs1,2)  << " #pm " << fround(nobserr1,2)  << endl;
//     s3 << "N(met > " << metval2 << ") " << fround(npred2,2) << " #pm " << fround(nprederr2,2) << endl;
//     s4 << "N(met > " << metval2 << ") " << fround(nobs2,2)  << " #pm " << fround(nobserr2,2)  << endl;

      
    
//     for(unsigned int iMC = 0 ; iMC < nMC ; ++ iMC ){
//       cout << delimstart << setw(width1) << mcleg.at(iMC)       << setw(width2)
//            << delim      << setw(width1) << fround(mcyield1[iMC],2)  << setw(width2)
//            << delim      << setw(width1) << fround(mcyield2[iMC],2)  << setw(width2) 
//            << delim      << setw(width1) << fround(mcyield3[iMC],2)  << setw(width2) << delim << endl;
//     }

//ttbar MC normalized to em data yield 
//    cout << delimstart << setw(width1) << "tt"              << setw(width2)
//         << delim      << setw(width1) << Form("%.2f %s %.2f",ntt1,ntterr1)  << setw(width2)
//         << delim      << setw(width1) << Form("%.2f %s %.2f",ntt2,ntterr2)  << setw(width2) 
//         << delim      << setw(width1) << Form("%.2f %s %.2f",ntt3,ntterr3)  << setw(width2) << delim << endl;
  






    //--------------------------------------------
    //get yields from ttbar MC, normalized to data
    //--------------------------------------------

//     float ntt1      = httmc_clone->Integral( bin1 , maxbin );
//     float ntterr1   = calculateHistError( httmc_clone , bin1 , maxbin );

//     float ntt2      = httmc_clone->Integral( bin2 , maxbin );
//     float ntterr2   = calculateHistError( httmc_clone , bin2 , maxbin );

//     float ntt3      = httmc_clone->Integral( bin3 , maxbin );
//     float ntterr3   = calculateHistError( httmc_clone , bin3 , maxbin );
