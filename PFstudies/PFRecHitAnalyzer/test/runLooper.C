#include "TChain.h"
#include "looper.C"
#include "TCanvas.h"
#include "TGraph.h"
//#include "TMultiGraph.h"

//#include "looper.h"

void runLooper(char* prefix , bool isData = false){

  TChain* ch = new TChain("Events");

  if( strcmp( prefix , "zmm" ) == 0 ){
    ch->Add("/tas05/disk00/benhoob/tcmetTestFiles/output/PFstudies_zmm.root");
  }

  if( strcmp( prefix , "qcd" ) == 0 ){
    ch->Add("/tas05/disk00/benhoob/tcmetTestFiles/output/PFstudies_qcd.root");
  }

  if( strcmp( prefix , "ttbar" ) == 0 ){
    ch->Add("/tas05/disk00/benhoob/tcmetTestFiles/output/PFstudies_ttbar_clusters.root");
  }
  
  looper* myLooper = new looper();
  
  cout << "Running on sample " << prefix << endl;

  //tune threshold so that sumet matches calosumet for each subdetector
  float eb_threshold  = 0.18;
  float ee_threshold  = 0.6;
  float hb_threshold  = 0.8;
  float he_threshold  = 1.5;
  float hfh_threshold = 1.7;
  float hfe_threshold = 8.;
  float hfshort_threshold = 0.;
  float hflong_threshold  = 0.;
  
  //calotower thresholds
//   float eb_threshold  = 0.2;
//   float ee_threshold  = 0.45;
//   float hb_threshold  = 0.7;
//   float he_threshold  = 0.8;
//   float hfh_threshold = 0.;
//   float hfe_threshold = 0.;
//   float hfshort_threshold = 0.85;
//   float hflong_threshold  = 0.5;

  //pfrechit seed thresholds
//   float eb_threshold  = 0.23;
//   float ee_threshold  = 0.6;
//   float hb_threshold  = 0.8;
//   float he_threshold  = 0.8;
//   float hfh_threshold = 1.4;
//   float hfe_threshold = 1.4;

  
  looper::metStruct mystruct = myLooper->ScanChain(ch, prefix, "_matchsumet" ,isData, -1,  
                                                   eb_threshold, ee_threshold, hb_threshold, 
                                                   he_threshold, hfh_threshold, hfe_threshold,
                                                   hfshort_threshold, hflong_threshold);
  
  /*
  const unsigned int n = 10;
  float x[n];
  float dmet_rms[n];
  float metxy_rms[n];
  float subdet_metxy_rms[n];

  const char* subdet = "hfe";
  const char* suffix = "";

  for( unsigned int i = 0 ; i < n ; ++i ){

    if( strcmp( subdet , "eb" ) == 0 ){
      eb_threshold = i * 0.05;
      x[i] = eb_threshold;
      suffix = Form("_eb%f",eb_threshold);
    }

    else if( strcmp( subdet , "ee" ) == 0 ){
      ee_threshold = i * 0.1;
      x[i] = ee_threshold;
      suffix = Form("_ee%f",ee_threshold);
    }

    else if( strcmp( subdet , "hb" ) == 0 ){
      hb_threshold = i * 0.2;
      x[i] = hb_threshold;
      suffix = Form("_hb%f",hb_threshold);
    }

    else if( strcmp( subdet , "he" ) == 0 ){
      he_threshold = i * 0.5;
      x[i] = he_threshold;
      suffix = Form("_he%f",he_threshold);
    }

    else if( strcmp( subdet , "hfh" ) == 0 ){
      hfh_threshold = i * 0.5;
      x[i] = hfh_threshold;
      suffix = Form("_hfh%f",hfh_threshold);
    }
    
    else if( strcmp( subdet , "hfe" ) == 0 ){
      hfe_threshold = i * 1.;
      x[i] = hfe_threshold;
      suffix = Form("_hfe%f",hfe_threshold);
    }

    looper::metStruct mystruct = myLooper->ScanChain(ch, prefix, suffix ,isData, -1,  
                                                     eb_threshold, ee_threshold, hb_threshold, 
                                                     he_threshold, hfh_threshold, hfe_threshold,
                                                     hfshort_threshold, hflong_threshold);
    
    
    cout << "dmet mean  : " << fround(mystruct.dmet_mean,2) << endl;
    cout << "dmet rms   : " << fround(mystruct.dmet_rms,2) << endl;
    cout << "met rms    : " << fround(mystruct.met_rms,2) << endl;
    cout << "EB met_rms : " << fround(mystruct.ebmet_rms,2) << endl;
 
    dmet_rms[i]         = mystruct.dmet_rms;
    metxy_rms[i]        = mystruct.met_rms;
    
    if     ( strcmp( subdet , "eb" ) == 0  )  subdet_metxy_rms[i] = mystruct.ebmet_rms;
    else if( strcmp( subdet , "ee" ) == 0  )  subdet_metxy_rms[i] = mystruct.eemet_rms;
    else if( strcmp( subdet , "hb" ) == 0  )  subdet_metxy_rms[i] = mystruct.hbmet_rms;
    else if( strcmp( subdet , "he" ) == 0  )  subdet_metxy_rms[i] = mystruct.hemet_rms;
    else if( strcmp( subdet , "hfh" ) == 0 )  subdet_metxy_rms[i] = mystruct.hfhmet_rms;
    else if( strcmp( subdet , "hfe" ) == 0 )  subdet_metxy_rms[i] = mystruct.hfemet_rms;
  }

  /*
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd();
  
  TGraph *gdmet_rms         = new TGraph( n , x , dmet_rms );
  TGraph *gmetxy_rms        = new TGraph( n , x , metxy_rms );
  TGraph *gsubdet_metxy_rms = new TGraph( n , x , subdet_metxy_rms );
  
  gdmet_rms->SetLineColor(2);
  gdmet_rms->SetMarkerColor(2);
  gmetxy_rms->SetLineColor(4);
  gmetxy_rms->SetMarkerColor(4);

  gdmet_rms->GetYaxis()->SetTitle("RMS (GeV)");
  gdmet_rms->GetXaxis()->SetTitle(Form("%s rechit threshold (GeV)",subdet));
  gdmet_rms->SetTitle("QCD MC");
  gdmet_rms->GetYaxis()->SetRangeUser(6,10);
 
  gdmet_rms->Draw("AP");
  gmetxy_rms->Draw("sameP");
  gsubdet_metxy_rms->Draw("sameP");

  TLegend *leg = new TLegend(0.8,0.8,1,1);
  leg->AddEntry( gdmet_rms , "met-genmet","p");
  leg->AddEntry( gmetxy_rms , "met(x/y)","p");
  leg->AddEntry( gsubdet_metxy_rms , "subdet met(x/y)","p");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();
  */
}

