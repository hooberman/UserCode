#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TLatex.h"
#include "TLine.h"
//#include "histtools.h"

using namespace std;

//-------------------------------------------------------

//const float scale_factor = 2;
const int   width1      = 20;
const int   width2      = 4;
//const char* iter        = "../output/V00-00-10/highpt/";
const char* iter        = "../output/V00-01-00/highpt/";
const bool  plotData    = true;
const int   nprec       = 2;
const char* data        = "data";
const bool  issignal    = false; //signal/control?

const char* metvar   = "y";
const char* dilptvar = "dilpt/sqrt(htpf)";

//-------------------------------------------------------

const bool  makeLatexPlot = true;        //plot in latex style
char* pm         = " +/- ";
char* delim      = "|";
char* delimstart = "|";
char* delimend   = "|";
char* ee         = "ee";
char* mm         = "mm";
char* em         = "em";

string plottitle  = "control region (125 < H_{T} < 300 GeV)";

pair<float,float> getMCratio( TH1F* hobs , TH1F* hpred , float cut );
void plotHist( TH1F* hpred , TH1F* hobs , TH1F* hpred_data , TH1F* hobs_data,
               string title, string xtitle , int rebin, float metcut);

inline double fround(double n, double d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}
void printRow( string sample , TH1F* hobs , TH1F* hpred , float cut , float mcratio = 1 );
void printLine(){

  if( makeLatexPlot ){
    cout << "\\hline" << endl;
  }
  else{
    cout << "-------------------------------------------------------------------------------------------------" << endl;
        
  }
}
void printHeader(){
  
  cout << delimstart << setw(width1) << "Sample"        << setw(width2)
       << delim      << setw(width1) << "Predicted"     << setw(width2)
       << delim      << setw(width1) << "Observed"      << setw(width2)
       << delim      << setw(width1) << "Obs/pred"      << setw(width2)
       << delimend   << endl;
}


void deleteHistos() {
   // Delete all existing histograms in memory
   TObject* obj;
   TList* list = gDirectory->GetList() ;
   TIterator* myiter = list->MakeIterator();
   while ((obj=myiter->Next())) {
     if (obj->IsA()->InheritsFrom(TH1::Class()) ||
         obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
   }
}


pair<float,float> victory( float ycut = 8.5 , float htcut = 300.0 , bool printgif =false );

void run(){

  /*
  pair<float,float> mcratio_pair = victory();
  float mcratio    = mcratio_pair.first;
  float mcratioerr = mcratio_pair.second;

  cout << "MC ratio " << mcratio << " +/- " << mcratioerr << endl;  
  */

  //bool y_vary = true;
  //const unsigned int n = 16;

  bool y_vary = false;
  const unsigned int n = 10;

  float ycut  = 8.5;
  //float ycut  = 13.0;
  float htcut = 300.0;


  float x[n];
  float xerr[n];
  float op[n];
  float operr[n];

  for( unsigned int i = 0 ; i < n ; ++i ){

    deleteHistos();

    if( y_vary ){
      ycut = 6.5 + i * 0.5;
      x[i]    = ycut;
      xerr[i] = 0.0;
    }

    else{
      htcut = 200 + 50 * i;
      x[i]    = htcut;
      xerr[i] = 0.0;
    }

    pair<float,float> mcratio_pair = victory( ycut , htcut , false );
    float mcratio    = mcratio_pair.first;
    float mcratioerr = mcratio_pair.second;

    op[i]    = mcratio;
    operr[i] = mcratioerr;

    cout << "y > " << ycut << " ht > " << htcut << " O/P " << Form("%.2f +/- %.2f",op[i],operr[i]) << endl;

  }



  TGraphErrors *gr  = new TGraphErrors(n,x,op ,xerr,operr);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  gPad->SetGridx(1);
  gPad->SetGridy(1);

  gr->SetLineColor(2);
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(20);

  if( y_vary ) gr->GetXaxis()->SetTitle("y cut (GeV^{1/2})");
  else         gr->GetXaxis()->SetTitle("H_{T} cut (GeV)");
  gr->GetYaxis()->SetTitle("observed / predicted");
  gr->Draw("AP");

}


pair<float,float> victory( float ycut  , float htcut  , bool printgif  ){

  cout << "y > " << ycut << " ht > " << htcut << endl;


  deleteHistos();

  if( issignal ) plottitle  = Form("signal region (H_{T} > %.0f GeV)",htcut);

  if( makeLatexPlot ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }

  int nbins  = 40;
  float xmin = 0;
  float xmax = 21.25;
  
  //int nbins  = 400;
  //float xmin = 0;
  //float xmax = 20.;

  //declare MC samples to include
  vector <string> samples;
  //samples.push_back("ttall");
  samples.push_back("ttpowheg");
  //samples.push_back("ttdil");
  //samples.push_back("ttotr");
  // samples.push_back("DYtot"); 
  // samples.push_back("wjets");
  // samples.push_back("ww");
  // samples.push_back("wz");
  // samples.push_back("zz");
  // samples.push_back("tW");
  //samples.push_back("LM1");


  //samples.push_back("DYee");
  //samples.push_back("DYmm");
  //samples.push_back("DYtautau");

  const unsigned int nsamples = samples.size();


  //-------------------------------
  //create MC chains and histos
  //-------------------------------
  TChain *ch[nsamples];
  TChain *ch_mctot = new TChain("t");
  
  TH1F* hmety[nsamples];
  TH1F* hdilpty[nsamples];
  TH1F* hdilpt[nsamples];

  TH1F* hmety_mctot     = new TH1F("hmety_mctot"    , "hmety_mctot"   , nbins,xmin,xmax);     
  TH1F* hdilpty_mctot   = new TH1F("hdilpty_mctot"  , "hdilpty_mctot" , nbins,xmin,xmax);     
  TH1F* hdilpt_mctot    = new TH1F("hdilpt_mctot"   , "hdilpt_mctot"  , 300,0,300);           

  hmety_mctot->Sumw2();    
  hdilpty_mctot->Sumw2();    
  hdilpt_mctot->Sumw2();       

  for( unsigned int i = 0 ; i < nsamples ; ++i ){
    ch[i] = new TChain("t");
    ch[i]->Add(Form( "%s/%s_smallTree.root" , iter , samples.at(i).c_str() ));
    ch_mctot->Add(Form( "%s/%s_smallTree.root" , iter , samples.at(i).c_str() ));

    //cout << Form( "%s/%s_smallTree.root" , iter , samples.at(i).c_str()) << endl;

    hmety[i]    =  new TH1F(Form("%s_hmety"    , samples.at(i).c_str()) , 
                            Form("%s_hmety"    , samples.at(i).c_str()) , nbins,xmin,xmax);

    hdilpty[i]  =  new TH1F(Form("%s_hdilpty"  , samples.at(i).c_str()) , 
                            Form("%s_hmety"    , samples.at(i).c_str()) , nbins,xmin,xmax); 
    
    hdilpt[i]   =  new TH1F(Form("%s_hdilpt"   , samples.at(i).c_str()) , 
                            Form("%s_hmety"    , samples.at(i).c_str()) , 300,0,300);           
    
    hdilpty[i]->Sumw2();
    hmety[i]->Sumw2();
    hdilpt[i]->Sumw2();
  }

  //create data chain and histos
  TChain *chdata = new TChain("t");
  if( plotData )
    chdata->Add(Form("%s/%s_smallTree.root",iter,data));
 
  TH1F* hmety_data          =  new TH1F("hmety_data",   "", nbins,xmin,xmax);
  TH1F* hdilpty_data        =  new TH1F("hdilpty_data", "", nbins,xmin,xmax); 
  TH1F* hdilpt_data         =  new TH1F("hdilpt_data",  "", 300,0,300); 

  hdilpty_data->Sumw2();
  hmety_data->Sumw2();

  //declare cuts for event selection
  char* jetcutstring = Form("htpf>125&&htpf<%.0f&&npfjets>1",htcut);
  if( issignal ) jetcutstring = Form("htpf>%.0f&&npfjets>1",htcut);
  TCut jetcut(jetcutstring);
  TCut zcut("passz==0");  
  TCut metcut("pfmet>50");
  //TCut weight("weight * ndavtxweight");
  TCut weight("weight");
  //TCut weight("weight * (1000./204.)");
  //TCut weight("1");

  //declare cuts for event selection
  //char* jetcutstring = "sumjetpt>125&&sumjetpt<300&&njets>1";
  //if( issignal ) jetcutstring = "sumjetpt>300&&njets>1";
  //TCut jetcut(jetcutstring);
  //TCut zcut("passz==0");  
  //TCut metcut("tcmet>50");
  //TCut weight("weight");

  TCut fcut       = ( jetcut + zcut + metcut ) * weight;
  TCut fcut_data  = ( jetcut + zcut + metcut );

  cout << "-------------------------------------------------------" << endl;
  cout << "Using selection: " << fcut.GetTitle() << endl;
  cout << "-------------------------------------------------------" << endl;


  TCanvas *ctemp =new TCanvas();
  ctemp->cd();

  //total MC
  //ch_mctot->Draw(Form("TMath::Min(%s,19.9)>>hmety_mctot"   , metvar   ) ,  fcut );
  //ch_mctot->Draw(Form("TMath::Min(%s,19.9)>>hdilpty_mctot" , dilptvar ) ,  fcut );
  ch_mctot->Draw("TMath::Min(dilpt,299.9)>>hdilpt_mctot"                ,  fcut );

  //apply dilepton pt scale factor
  //float k = hdilpt_mctot->Integral() /  hdilpt_mctot->Integral( hdilpt_mctot->FindBin(50) , 100000 );
  float k = hdilpt_mctot->Integral() /  hdilpt_mctot->Integral( hdilpt_mctot->FindBin(50) , 100000 );
  cout << "Scaling factor " << k << " " << hdilpt_mctot->Integral() << " / " << hdilpt_mctot->Integral( hdilpt_mctot->FindBin(50) , 100000 ) << endl;
  hdilpty_mctot->Scale( k );

  //data
  if( plotData ){
    chdata->Draw(Form("TMath::Min(%s,19.9)>>hmety_data"   , metvar   ) ,  fcut_data );
    chdata->Draw(Form("TMath::Min(%s,19.9)>>hdilpty_data" , dilptvar ) ,  fcut_data );
    chdata->Draw("TMath::Min(dilpt,299.9)>>hdilpt_data",                  fcut_data );
  }
  delete ctemp;

  float kdata = 1;
  if( plotData ){
    kdata = k;
    cout << "Using data scale factor from MC " << k << endl;
    hdilpty_data->Scale( kdata );
  }

  
  //MC
  TH1F* hmety_temp     = new TH1F("hmety_temp"    , "" , nbins,xmin,xmax);     
  TH1F* hdilpty_temp   = new TH1F("hdilpty_temp"  , "" , nbins,xmin,xmax);     
  TH1F* hdilpt_temp    = new TH1F("hdilpt_temp"   , "" , 300 , 0 , 300);     

  hmety_temp->Sumw2();
  hdilpty_temp->Sumw2();

  TCut fcut_MC  = ( jetcut + zcut + metcut );

  TCut fcut_MC_clone = fcut_MC;

  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    fcut_MC = fcut_MC_clone;
    if( samples[i] == "DYtot" ) fcut_MC = fcut_MC + "ntaus==2";

    ch[i]->Draw(Form("TMath::Min(%s,19.9)>>hmety_temp"    , metvar   )  , fcut_MC * weight );
    ch[i]->Draw(Form("TMath::Min(%s,19.9)>>hdilpty_temp"  , dilptvar )  , fcut_MC * weight );
    ch[i]->Draw("TMath::Min(dilpt,299.9)>>hdilpt_temp",                   fcut_MC * weight );

    hmety[i]   = (TH1F*) hmety_temp->Clone(   Form("%s_hmety"    , samples.at(i).c_str()) );
    hdilpty[i] = (TH1F*) hdilpty_temp->Clone( Form("%s_hdilpty"  , samples.at(i).c_str()) );
    hdilpt[i]  = (TH1F*) hdilpt_temp->Clone(  Form("%s_hdilpt"   , samples.at(i).c_str()) );

    hdilpty[i]->Scale( k );
    hmety_temp->Reset();
    hdilpty_temp->Reset();

    hmety_mctot->Add(hmety[i]);
    hdilpty_mctot->Add(hdilpty[i]);
  }


  //dineutrino pt
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd();  
  
  TPad *plotpad = new TPad("plotpad","",0,0,1,0.9);
  plotpad->Draw();
  plotpad->cd();


  plotHist( hdilpty_mctot , hmety_mctot , hdilpty_data , hmety_data , 
            plottitle , "y   (GeV^{1/2})" , 4 , ycut );

  //if( printgif ) c1->Print("plots/victory.gif");

  
  c1->cd();
  TPad *titlepad = new TPad("titlepad","",0,0.9,1,1.0);
  titlepad->Draw(); 
  titlepad->cd();
  
  TLatex *title = new TLatex();
  title->SetTextSize(0.5);
  title->SetNDC();
  title->DrawLatex(0.18,0.15,"CMS");
  title->DrawLatex(0.58,0.15,"L_{int} = 191 pb^{-1},  #sqrt{s} = 7 TeV");

  cout << endl << endl;

  printLine();
  printHeader();
  printLine();

  for( unsigned int i = 0 ; i < nsamples ; ++i ){
    printRow( samples.at(i) , hmety[i] , hdilpty[i] , ycut );
  }


  printLine();
  printRow( "total MC" , hmety_mctot , hdilpty_mctot , ycut );
  printLine();

  //------------------------------
  // scale data by MC obs/pred
  //------------------------------

  pair<float,float> mcratio_pair = getMCratio(  hmety_mctot , hdilpty_mctot , ycut );  
  float mcratio    = mcratio_pair.first;
  float mcratioerr = mcratio_pair.second;

  if( plotData ){
    printRow( "data" , hmety_data , hdilpty_data , ycut , mcratio );
    printLine();
  }

  if( printgif ){
    
    if( issignal ){
      c1->Print("../plots/victory_signal.pdf");
      c1->Print("../plots/victory_signal.png");
    }
    else{
      c1->Print("../plots/victory_control.pdf");
      c1->Print("../plots/victory_control.png");
    }
  }


  return make_pair( mcratio , mcratioerr );

 
}


float calculateHistError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}


pair<float,float> getMCratio( TH1F* hobs , TH1F* hpred , float cut ){

  int   minbin     = hobs->FindBin( cut );
  int   maxbin     = 10000;

  float pred    = hpred->Integral( minbin , maxbin);
  float prederr = calculateHistError( hpred , minbin , maxbin );

  float obs     = hobs->Integral(  minbin , maxbin );
  float obserr  = calculateHistError( hobs , minbin , maxbin );

  float ratio   = pred > 0 ? obs / pred : 0;
  float ratioerr = obs > 0 && pred > 0 ? ratio * sqrt( pow( prederr / pred , 2 ) + pow( obserr / obs , 2 ) ) : 0;

  return make_pair( ratio , ratioerr );
}

void printRow( string sample , TH1F* hobs , TH1F* hpred , float cut , float mcratio ){

  stringstream sobs;
  stringstream spred;
  stringstream sratio;

  int   minbin     = hobs->FindBin( cut );
  int   maxbin     = 10000;

  float pred    = hpred->Integral( minbin , maxbin);
  float prederr = calculateHistError( hpred , minbin , maxbin );

  if( sample == "data" ){
    pred    *= mcratio;
    prederr *= mcratio;
    cout << "Scaling data " << mcratio << endl;
  }


  float obs     = hobs->Integral(  minbin , maxbin );
  float obserr  = calculateHistError( hobs , minbin , maxbin );

  float ratio    = pred > 0 ? obs / pred : 0;
  float ratioerr = obs > 0 && pred > 0 ? ratio * sqrt( pow( prederr / pred , 2 ) + pow( obserr / obs , 2 ) ) : 0;
  
  sobs   << Form( "%.1f" , obs   ) << pm << Form( "%.1f" , obserr   ); 
  spred  << Form( "%.1f" , pred  ) << pm << Form( "%.1f" , prederr  ); 
  sratio << Form( "%.2f" , ratio ) << pm << Form( "%.2f" , ratioerr ); 

  cout << delimstart << setw(width1) << sample        << setw(width2)
       << delim      << setw(width1) << spred.str()   << setw(width2)
       << delim      << setw(width1) << sobs.str()    << setw(width2)
       << delim      << setw(width1) << sratio.str()  << setw(width2)
       << delimend   << endl;

}




void plotHist( TH1F* hpred , TH1F* hobs , TH1F* hpred_data, TH1F* hobs_data,
               string title, string xtitle , int rebin , float metcut){

  hpred->SetLineColor(4);
  hobs ->SetLineColor(2);
  hpred->SetMarkerColor(4);
  hobs ->SetMarkerColor(2);
  hpred->SetMarkerStyle(23);
  hobs->SetMarkerStyle(4);
  hpred_data->SetLineColor(4);
  hobs_data ->SetLineColor(2);
  hpred_data->SetMarkerColor(4);
  hobs_data ->SetMarkerColor(2);
  hpred_data->SetMarkerStyle(23);
  hobs_data->SetMarkerStyle(4);
  hpred_data->SetMarkerSize(2);
  hobs_data->SetMarkerSize(2);
  hpred->SetFillColor(0);
  hobs ->SetFillColor(0);
  hpred->SetTitle( title.c_str() );
  hpred->GetXaxis()->SetTitle( xtitle.c_str() );
  hpred->GetXaxis()->SetLabelSize(0.04);
  hpred->GetXaxis()->SetTitleSize(0.055);
  hpred->GetYaxis()->SetTitle( "Events   " );
  hpred->GetYaxis()->SetTitleSize(0.055);
  hpred->GetYaxis()->SetTitleOffset(1);
  //hpred->GetYaxis()->SetTitle( "Entries / 2.125 GeV" );

  hpred->SetMinimum(0.01);
  hpred->SetMaximum(10000.);

  int metbin = hobs->FindBin( metcut );
  float pred = hpred->Integral( metbin , 100000000);
  float obs  = hobs->Integral(  metbin , 100000000);

  float pred_data = hpred_data->Integral( metbin , 100000000);
  float obs_data  = hobs_data->Integral(  metbin , 100000000);

  //cout << "Predicted (met>" << metcut << ") " << pred << endl;
  //cout << "Observed  (met>" << metcut << ") " << obs  << endl;
  //cout << "obs/pred            " << obs/pred << endl;

  //int width1 = 10;
  //int width2 = 2;
  /*
  cout << endl;
  cout << "|" << setw(width1) << ""          << setw(width2)
       << "|" << setw(width1) << "Predicted" << setw(width2)
       << "|" << setw(width1) << "Observed"  << setw(width2)
       << "|" << setw(width1) << "Obs/Pred"  << setw(width2) << "|" << endl;

  cout << "|" << setw(width1) << "MC"      << setw(width2)
       << "|" << setw(width1) << fround(pred,nprec)      << setw(width2)
       << "|" << setw(width1) << fround(obs,nprec)       << setw(width2)
       << "|" << setw(width1) << fround(obs/pred,nprec)  << setw(width2) << "|" << endl;

  if( plotData ){
    cout << "|" << setw(width1) << "data"      << setw(width2)
         << "|" << setw(width1) << fround(pred_data,nprec)      << setw(width2)
         << "|" << setw(width1) << fround(obs_data,nprec)       << setw(width2)
         << "|" << setw(width1) << fround(obs_data/pred_data,nprec)  << setw(width2) << "|" << endl;
  }
  cout << endl;
  */

  stringstream spred;
  spred << "Pred: " << pred ;
  stringstream sobs;
  sobs << "Obs:  " << obs ;

  if( rebin > 1 ){
    hobs->Rebin ( rebin );
    hpred->Rebin( rebin );
    hobs_data->Rebin ( rebin );
    hpred_data->Rebin( rebin );
  }
 
  if( plotData ){
    hpred->Draw("hist");
    hobs->Draw("samehist");
    hpred_data->Draw("sameE1");
    hobs_data->Draw("sameE1");
    hpred->SetLineWidth(2);
    hobs->SetLineWidth(2);
  }else{
    hpred->Draw("E1");
    hobs->Draw("sameE1");
  } 
  gPad->SetLogy(1);

  //TLatex t;
  //t.SetNDC();
  //t.DrawLatex(0.55,0.65,Form("obs/pred = %.2f",obs/pred));

  TLegend *leg = new TLegend(0.64,0.72,0.98,0.95);
  //leg->AddEntry(hpred, Form("pred>%.1f GeV = %.3f",metcut,pred) ,"l");
  //leg->AddEntry(hobs,  Form("obs>%.1f GeV = %.3f",metcut,obs)  ,"l");
  leg->AddEntry(hpred, "MC predicted" ,"l");
  leg->AddEntry(hobs,  "MC observed"  ,"l");
  if( plotData ){
    leg->AddEntry(hpred_data, "data predicted");
    leg->AddEntry(hobs_data,  "data observed" );
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.DrawLine( metcut , 0.01 , metcut , 10000 );
  //line.DrawLine( metcut , 0.5 * TMath::Min( hobs->GetMinimum() , hpred->GetMinimum() ),
  //               metcut , 2 *   TMath::Max( hobs->GetMaximum() , hpred->GetMaximum() ) );

  if( plotData ){
    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextSize(0.045);
    //t->DrawLatex(0.18,0.32,"CMS");
    //t->DrawLatex(0.18,0.27,"34 pb^{-1}  at  #sqrt{s} = 7 TeV");
    t->DrawLatex(0.20,0.88,"Events with ee/#mu#mu/e#mu");
    
    if( issignal ){
      t->DrawLatex(0.20,0.82,"H_{T} > 300 GeV");
      //t->SetTextSize(0.05);
      //t->DrawLatex(0.8,0.6,"(a)");
    }
    else{
      t->DrawLatex(0.20,0.82,"125 < H_{T} < 300 GeV");
      //t->SetTextSize(0.05);
      //t->DrawLatex(0.8,0.6,"(b)");
    }

 }
}
