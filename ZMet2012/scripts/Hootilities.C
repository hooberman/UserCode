#include "Hootilities.h"
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
#include "TF1.h"

float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

void printLine( bool latex ){

  if( latex ){
    cout << "\\hline" << endl;
  }
  else{
    for( int i = 0 ; i < linelength ; ++i ) cout << "-";
    cout << endl;
  }
}

void printHeader(){

  cout << delimstart << setw(width1) << "Sample"    << setw(width2)
       << delim      << setw(width1) << ee          << setw(width2)
       << delim      << setw(width1) << mm          << setw(width2)
       << delim      << setw(width1) << em          << setw(width2)
       << delim      << setw(width1) << "tot"       << setw(width2) 
       << delimend   << endl;

}


void print( TH1F* h , string label , bool correlatedError = false ){

  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  if( label == "data" ){
    see  << Form( "%.0f" , h->GetBinContent(1) );
    smm  << Form( "%.0f" , h->GetBinContent(2) );
    sem  << Form( "%.0f" , h->GetBinContent(3) );
    stot << Form( "%.0f" , h->Integral()       );
  }else{
    //see  << Form( "%.1f" , h->GetBinContent(1) );
    //smm  << Form( "%.1f" , h->GetBinContent(2) );
    //sem  << Form( "%.1f" , h->GetBinContent(3) );
    //stot << Form( "%.1f" , h->Integral()       );
    
    see  << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
    smm  << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
    sem  << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );
    
    float error = 0;
    if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
    else                  error = histError(h,1,4);
    
    stot << Form( "%.1f" , h->Integral()       ) << pm << Form( "%.1f" , error  );
  }

  cout << delimstart << setw(width1) << label      << setw(width2)
       << delim      << setw(width1) << see.str()  << setw(width2)
       << delim      << setw(width1) << smm.str()  << setw(width2)
       << delim      << setw(width1) << sem.str()  << setw(width2)
       << delim      << setw(width1) << stot.str() << setw(width2)
       << delimend   << endl;
  
  
}


#include <TList.h>
#include <TIterator.h>

void deleteHistos() {
   // Delete all existing histograms in memory
   TObject* obj;
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();
   while ((obj=iter->Next())) {
     if (obj->IsA()->InheritsFrom(TH1::Class()) ||
         obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
   }
}

void initSymbols( bool latex ){

  //-------------------------------------------------------
  // table format
  //-------------------------------------------------------

  width1      = 15;
  width2      = 4;
  linelength  = (width1+width2)*5+1;

  //-------------------------------------------------------
  // symbols
  //-------------------------------------------------------
  
  if( latex ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "$ee$";
    mm         = "$\\mu\\mu$";
    em         = "$e\\mu$";
  }else{
    pm         = " +/- ";
    delim      = "|";
    delimstart = "|";
    delimend   = "|";
    ee         = "ee";
    mm         = "mm";
    em         = "em";
  }

}

//TH1F* doDYestimate( TChain *ch , TCut sel ){
void doDYestimate( TChain *ch , TCut sel , TH1F* hyield ){

  //-------------------
  // set parameters
  //-------------------

  bool  verbose = true;
  float R       = 0.13;
  float R_err   = 0.07;
  float k       = 1.11;

  if( verbose ){
    cout << "R : " << R << " +/- " << R_err << endl;
    cout << "k : " << k << endl;
  }

  //-------------------
  // invert Z-veto
  //-------------------

  TString selstring(sel.GetTitle());
  selstring.ReplaceAll("passz == 0","dilmass>76&&dilmass<106");
  TCut newsel = TCut(selstring);

  if( verbose ){
    cout << "Pre  : " << sel.GetTitle() << endl;
    cout << "Post : " << newsel.GetTitle() << endl;
  }

  //-------------------
  // get data yields
  //-------------------

  float ndata_ee_in = ch->GetEntries(newsel+"leptype==0");
  float ndata_mm_in = ch->GetEntries(newsel+"leptype==1");
  float ndata_em_in = ch->GetEntries(newsel+"leptype==2");

  if( verbose ){
    cout << "nee     " << ndata_ee_in << endl;
    cout << "nmm     " << ndata_mm_in << endl;
    cout << "nem     " << ndata_em_in << endl;
  }

  //-------------------
  // do DY estimate
  //-------------------

  float neepred     = R     * ( ndata_ee_in - (0.5/k)  * ndata_em_in );
  float nmmpred     = R     * ( ndata_mm_in - k/2.     * ndata_em_in );

  float neeprederr  = 0.;
  float nmmprederr  = 0.;

  neeprederr  += pow( R_err * ( ndata_ee_in - (0.5/k) * ndata_em_in ) , 2 );
  nmmprederr  += pow( R_err * ( ndata_mm_in - k/2.    * ndata_em_in ) , 2 );

  neeprederr  += pow( R * sqrt( ndata_ee_in + pow(0.5/k,2)  * ndata_em_in ) , 2 );
  nmmprederr  += pow( R * sqrt( ndata_mm_in + pow(k/2.,2)   * ndata_em_in ) , 2 );

  neeprederr = sqrt(neeprederr);
  nmmprederr = sqrt(nmmprederr);

  float ntotpred    = neepred + nmmpred;
  float ntotprederr = sqrt( pow(neeprederr,2) + pow(nmmprederr,2) );

  if( verbose ){
    cout << "nee  pred " << neepred  << " +/- " << neeprederr  << endl;
    cout << "nmm  pred " << nmmpred  << " +/- " << nmmprederr  << endl;
    cout << "ntot pred " << ntotpred << " +/- " << ntotprederr << endl;
  }

  //-------------------
  // set hist contents
  //-------------------

  hyield->SetBinContent( 1 , neepred );
  hyield->SetBinContent( 2 , nmmpred );
  hyield->SetBinContent( 3 , 0 );
  hyield->SetBinContent( 4 , 0 );

  hyield->SetBinError( 1 , neeprederr );
  hyield->SetBinError( 2 , nmmprederr );
  hyield->SetBinError( 3 , 0 );
  hyield->SetBinError( 4 , 0 );

}

void printYields( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , TCut sel , TCut weight , bool latex ){

  initSymbols( latex );

  TCanvas *ctemp = new TCanvas();

  printLine(latex);
  printHeader();
  printLine(latex);

  TH1F* hyield = new TH1F("hyield","yield",4,0,4);
  TH1F* hmctot = new TH1F("hmctot","hmctot",4,0,4);
  hyield->Sumw2();
  hmctot->Sumw2();

  //----------------------
  // print SM MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < chmc.size() ; imc++){

    bool correlatedError = false;

    if( TString(labels[imc]).Contains("LM") || TString(labels[imc]).Contains("T2tt") ) continue;

    // data-driven DY estimate
    if( strcmp(labels[imc],"DYdata")   == 0 ){
      //hyield = doDYestimate( chmc[imc] , sel );
      doDYestimate( chmc[imc] , sel , hyield );
    }

    // fake estimate
    else if( strcmp(labels[imc],"single fakes")   == 0 || strcmp(labels[imc],"double fakes") == 0){

      //correlatedError = true;

      TString weightstring(weight.GetTitle());
      weightstring.ReplaceAll("ndavtxweight","1");
      TCut newweight = TCut(weightstring);

      chmc[imc]->Draw("leptype>>hyield",sel*newweight);

      // SF --> SF - 2 X DF
      if( strcmp(labels[imc],"single fakes")   == 0 ){
	
	TH1F *hyielddf = new TH1F("hyielddf","hyielddf",4,0,4);

	chmc[imc+1]->Draw("leptype>>hyielddf",sel*newweight);
	hyield->Add(hyielddf,-2);
      }      

      hyield->SetBinError(1,0.5*hyield->GetBinContent(1));
      hyield->SetBinError(2,0.5*hyield->GetBinContent(2));
      hyield->SetBinError(3,0.5*hyield->GetBinContent(3));
      
    }

    //vanilla MC
    else{
      chmc[imc]->Draw("leptype>>hyield",sel*weight);

      //do efficiency correction
      /*
      //ee
      hyield->SetBinContent  ( 1 , hyield->GetBinContent(1) * 0.95);
      hyield->SetBinError    ( 1 , hyield->GetBinError(1)   * 0.95);

      //mumu
      hyield->SetBinContent  ( 2 , hyield->GetBinContent(2) * 0.88);
      hyield->SetBinError    ( 2 , hyield->GetBinError(2)   * 0.88);

      //emu
      hyield->SetBinContent  ( 3 , hyield->GetBinContent(3) * 0.92);
      hyield->SetBinError    ( 3 , hyield->GetBinError(3)   * 0.92);
      */
      //cout << "printYieldTable apply trigger efficiencies 0.95, 0.88, 0.92" << endl;
    }

    if( imc == 0 ) hmctot = (TH1F*) hyield->Clone();
    else           hmctot->Add(hyield);
    
    print( hyield , labels[imc] , correlatedError );

    //hyield->Reset();
  }

  printLine(latex);

  //-------------------------------
  // print sum of SM MC samples
  //-------------------------------

  print( hmctot , "tot SM MC" );

  printLine(latex);
 
  chdata->Draw("leptype>>hyield",sel);

  print( hyield , "data" );
    
  printLine(latex);

  //----------------------
  // print SUSY MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < chmc.size() ; imc++){

    if( !TString(labels[imc]).Contains("LM") )   continue;
    if( !TString(labels[imc]).Contains("T2tt") ) continue;

    chmc[imc]->Draw("leptype>>hyield",sel*weight);

    //do efficiency correction
    //hyield->SetBinContent  ( 2 , hyield->GetBinContent(2) * 0.90);
    //hyield->SetBinContent  ( 3 , hyield->GetBinContent(3) * 0.95);
    //hyield->SetBinError    ( 2 , hyield->GetBinError(2)   * 0.90);
    //hyield->SetBinError    ( 3 , hyield->GetBinError(3)   * 0.95);
    
    print( hyield , labels[imc] );

  }

  printLine(latex);


  
  delete ctemp;
}


TLegend *getLegend( vector<TChain*> chmc , vector<char*> labels , bool overlayData, float x1, float y1, float x2, float y2){

  //int colors[]={6,2,7,4,5,8,9,15,12};
  int colors[]={4,7,2,5,8,9,15,12};
  
  TLegend *leg = new TLegend(x1,y1,x2,y2);

  TH1F*    datahist = new TH1F("datahist","datahist",1,0,1);
  datahist->Sumw2();

  if( overlayData ) leg->AddEntry(datahist,"data");

  const int nmc = chmc.size();
  TH1F*    mchist[nmc];

  //-----------------
  // SM samples
  //-----------------

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( TString(t).Contains("LM") )   continue;
    if( TString(t).Contains("T2tt") ) continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    if( TString( labels.at(imc) ).Contains("LM") || TString( labels.at(imc) ).Contains("LM") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    if( strcmp("tt",t)      == 0 ) t = "t#bar{t}";
    if( strcmp("ttll",t)    == 0 ) t = "t#bar{t} #rightarrow ll";
    if( strcmp("tttau",t)   == 0 ) t = "t#bar{t} #rightarrow l#tau/#tau#tau";
    if( strcmp("ttfake",t)  == 0 ) t = "t#bar{t} #rightarrow fake";
    if( strcmp("t",t)       == 0 ) t = "single top";
    if( strcmp("wjets",t)   == 0 ) t = "W+jets";
    if( strcmp("zjets",t)   == 0 ) t = "Z+jets";
    if( strcmp("ww",t)      == 0 ) t = "WW";
    if( strcmp("wz",t)      == 0 ) t = "WZ";
    if( strcmp("zz",t)      == 0 ) t = "ZZ";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  //-----------------
  // LM samples
  //-----------------

  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    char* t = labels.at(imc);

    if( ! (TString(t).Contains("LM")||TString(t).Contains("T2tt")) )   continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    if( TString( labels.at(imc) ).Contains("LM") || TString( labels.at(imc) ).Contains("T2tt") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    if( strcmp("ttall",t) == 0 ) t = "t#bar{t}";
    if( strcmp("t",t)     == 0 ) t = "single top";
    if( strcmp("wjets",t) == 0 ) t = "W+jets";
    if( strcmp("WW",t)    == 0 ) t = "W^{+}W^{-}";
    if( strcmp("WZ",t)    == 0 ) t = "W^{#pm}Z^{0}";
    if( strcmp("ZZ",t)    == 0 ) t = "Z^{0}Z^{0}";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],t,"f");
    
  }

  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  
  return leg;

}

void compareDataMC( vector<TChain*> chmc , vector<char*> labels , TChain* chdata , char* var , 
		    TCut sel , TCut weight , int nbins ,  float xmin , float xmax ,  
		    char* xtitle , bool overlayData , bool residual , bool drawLegend , bool log , bool normalize , char* flavor ){

  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();

  if( residual ){
    fullpad = new TPad("fullpad","fullpad",0,0,1,1);
    fullpad->Draw();
    fullpad->cd();

    plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
    plotpad->Draw();
    plotpad->cd();
    if( log ) plotpad->SetLogy();
  }
  else{
    if( log ) gPad->SetLogy();
  }

  TString tvar(var);
  tvar.ReplaceAll("()","");
  tvar.ReplaceAll("(","");
  tvar.ReplaceAll(")","");
  tvar.ReplaceAll(".","");
  tvar.ReplaceAll("*","");
  tvar.ReplaceAll(" ","");
  const char* myvar = tvar;

  cout << "Plotting var " << myvar << " flavor " << flavor << endl;

  //int colors[]={6,2,7,4,5,8,9,15,12};
  //int colors[]={kRed+2,5,7,5,5,8,9,15,12};
  int colors[]={4,7,2,5,8,9,15,12};

  assert( chmc.size() == labels.size() );
  const unsigned int nmc = chmc.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist = new TH1F();
  TH1F*    mchist[nmc];
  TH1F*    datahist = new TH1F(Form("%s_datahist_%s",myvar,flavor),Form("%s_datahist_%s",myvar,flavor),nbins,xmin,xmax);

  float trigeff = 1.0;
  // if     ( TString(flavor).Contains("ee")  ) trigeff = 0.95;
  // else if( TString(flavor).Contains("mm")  ) trigeff = 0.88;
  // else if( TString(flavor).Contains("em")  ) trigeff = 0.92;
  // else if( TString(flavor).Contains("all") ) trigeff = 0.91;

  cout << "compareDataMC : apply trigeff " << trigeff << endl;

  if     ( TString(flavor).Contains("ee")  ) sel+="leptype==0";
  else if( TString(flavor).Contains("mm")  ) sel+="leptype==1";
  else if( TString(flavor).Contains("em")  ) sel+="leptype==2";
  //else if( TString(flavor).Contains("all") ) 

  TCut trigweight(Form("%.2f",trigeff));

  TH1F* htemp = new TH1F("htemp","htemp",1,0,1);

  float nmctot = 0.0;

  for( int imc = nmc-1 ; imc > -1 ; imc-- ){
    chmc.at(imc)->Draw("0.5>>htemp",sel*weight*trigweight);
    nmctot += htemp->Integral();
  }

  chdata->Draw("0.5>>htemp",sel);
  float ndata = htemp->Integral();

  float SF = ndata/nmctot;
  if( normalize ){
    cout << "Data, MC, SF " << ndata << ", " << nmctot << ", " << SF << endl;
    cout << "Scaling MC by " << SF << "!!!!!!" << endl;
  }

  float mcintegral = 0.0;
  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc-1 ; imc > -1 ; imc-- ){

    mchist[imc] = new TH1F(Form("%s_mc_%i_%s",myvar,imc,flavor),Form("%s_mc_%i_%s",myvar,imc,flavor),nbins,xmin,xmax);
    mchist[imc]->Sumw2();

    chmc.at(imc)->Draw(Form("TMath::Min(%s,%f)>>%s_mc_%i_%s",var,xmax-0.01,myvar,imc,flavor),sel*weight*trigweight);

    if( normalize ) mchist[imc]->Scale(SF);

    if( TString( labels.at(imc) ).Contains("LM") || TString( labels.at(imc) ).Contains("T2tt") ){
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineStyle(2);
    }else{
      mchist[imc]->SetFillColor( colors[imc] );
    }

    // if( strcmp(labels[imc],"ttfake")  == 0 || strcmp(labels[imc],"wjets")  == 0 ){
    //   cout << "Scaling " << labels[imc] << " by 3.8" << endl;
    //   mchist[imc]->Scale(3.8);
    // }

    mcstack->Add( mchist[imc] );

    //if( imc == 0 ) mctothist = (TH1F*) mchist[imc]->Clone();
    if( imc == nmc-1 ) mctothist = (TH1F*) mchist[imc]->Clone();
    else               mctothist->Add(mchist[imc]);

    cout << "MC yield " << labels[imc] << " " << Form("%.2f",mchist[imc]->Integral()) << endl;
    mcintegral += mchist[imc]->Integral();
  }

  cout << "MC total yield " << Form("%.2f",mcintegral) << endl;

  chdata->Draw(Form("TMath::Min(%s,%f)>>%s_datahist_%s",var,xmax-0.01,myvar,flavor),sel);

  if( overlayData ){

    float max = datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin());
    if( mctothist->GetMaximum() > max ) max = mctothist->GetMaximum();
    if( log ) datahist->SetMaximum( 15 * max );
    else      datahist->SetMaximum( 1.4 * max );

    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->GetXaxis()->SetTitleSize(0.05);
    datahist->GetXaxis()->SetLabelSize(0.04);
    datahist->GetXaxis()->SetTitleOffset(1.2);
    datahist->Draw("E1");
    mcstack->Draw("samehist");
    datahist->Draw("sameE1");
    datahist->Draw("sameaxis");
    
    if(!log) datahist->GetYaxis()->SetRangeUser(0.,1.4*max);
    
    cout << "data yield " << datahist->Integral() << endl;

  }
  else{
    mctothist->GetXaxis()->SetTitle(xtitle);
    mctothist->Draw();
    mcstack->Draw("same");
    mctothist->Draw("sameaxis");
  }

  if( drawLegend ){
    TLegend* myleg = getLegend( chmc , labels , overlayData );
    myleg->Draw();
  }

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  text->DrawLatex(0.2,0.88,"CMS Preliminary");
  //text->DrawLatex(0.2,0.83,"0.98 fb^{-1} at #sqrt{s} = 7 TeV");
  text->DrawLatex(0.2,0.83,"#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 5.1 fb^{-1}");

  if     ( TString(flavor).Contains("ee")  ) text->DrawLatex(0.2,0.78,"Events with ee");
  else if( TString(flavor).Contains("mm")  ) text->DrawLatex(0.2,0.78,"Events with #mu#mu");
  else if( TString(flavor).Contains("em")  ) text->DrawLatex(0.2,0.78,"Events with e#mu");
  else if( TString(flavor).Contains("all") ) text->DrawLatex(0.2,0.78,"Events with ee/#mu#mu/e#mu");

  if( residual ){
    fullpad->cd();

    respad = new TPad("respad","respad",0,0.8,1,1);
    respad->Draw();
    respad->cd();

    gPad->SetGridy();

    TH1F* ratio = (TH1F*) datahist->Clone(Form("%s_ratio",datahist->GetName()));
    ratio->Divide(mctothist);

    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetLabelSize(0.2);
    //ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetRangeUser(0.6,1.4);
    ratio->GetYaxis()->SetTitle("data/MC  ");
    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetXaxis()->SetTitleSize(0);
    ratio->SetMarkerSize(0.7);

    TF1* fpol1 = new TF1("fpol1","pol1",datahist->GetXaxis()->GetXmin(),datahist->GetXaxis()->GetXmax());
    ratio->Fit(fpol1,"R");

    ratio->Draw();

    TLine line;
    line.SetLineWidth(1);
    line.DrawLine(datahist->GetXaxis()->GetXmin(),1,datahist->GetXaxis()->GetXmax(),1);

  }






}
