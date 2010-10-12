{
  bool useFit = false;

  gStyle->SetOptFit(0);

  //TFile *f1 = TFile::Open("crabjobs/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/res/merged.root");
  //TFile *f2 = TFile::Open("crabjobs/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1_TEC5_BP10/res/merged.root");

  //TFile *f1 = TFile::Open("crabjobs/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1_PEAK/res/merged.root");
  //TFile *f2 = TFile::Open("crabjobs/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/res/merged.root");

  //TFile *f1 = TFile::Open("crabjobs/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1_PEAK/res/merged.root");
  //TFile *f2 = TFile::Open("crabjobs/MinimumBias_Run2010A-PromptReco-v4_GR10_P_V7_peakLA/res/merged.root");
  //TFile *f2 = TFile::Open("MinimumBias_Run2010A-PromptReco-v4_GR10_P_V7_peakLA_dwvsdtantheta/res/merged.root");

  TFile *f1 = TFile::Open("MinimumBias_Commissioning10-GOODCOLL-Jun14thSkim_v1/res/merged.root");
  //TFile *f2 = TFile::Open("MinimumBias_Run2010A-PromptReco-v4_GR10_P_V7_peakLA/res/merged.root");
  TFile *f2 = TFile::Open("MinimumBias_Run2010A-PromptReco-v4_GR10_P_V7_peakLA_dwvsdtantheta/res/merged.root");

  const unsigned int nrings = 7;

  TH1F *hdw1[nrings];
  TH1F *hdw2[nrings];

  float mean1[nrings];
  float mean2[nrings];
  float diff[nrings];

  float mean1err[nrings];
  float mean2err[nrings];
  float differr[nrings];

  TLatex t;
  t.SetNDC();

  TCanvas *c1 = new TCanvas("c1","",1600,800);
  c1->Divide(4,2);

  for( unsigned int ring = 0 ; ring < nrings ; ++ring ){
  
    c1->cd( ring + 1 );

    //hdw1[ring] = (TH1F*) f1->Get(Form("PeakDecoResiduals/PeakDecoResiduals/hdw_TEC_wheel5_ring%i_stereo0",ring+1));
    //hdw2[ring] = (TH1F*) f2->Get(Form("PeakDecoResiduals/histos/hdw_TEC_cut_wheel5_ring%i_stereo0",ring+1));

    //hdw1[ring] = (TH1F*) f1->Get(Form("PeakDecoResiduals/PeakDecoResiduals/hdw_TEC_wheel0_ring%i_stereo1",ring+1));
    //hdw2[ring] = (TH1F*) f2->Get(Form("PeakDecoResiduals/PeakDecoResiduals/hdw_TEC_wheel0_ring%i_stereo1",ring+1));

    //hdw1[ring] = (TH1F*) f1->Get(Form("PeakDecoResiduals/PeakDecoResiduals/hdw_TEC_wheel0_ring%i_stereo0",ring+1));
    //hdw2[ring] = (TH1F*) f2->Get(Form("PeakDecoResiduals/histos/hdw_TEC_wheel0_ring%i_stereo0",ring+1));
 
    hdw1[ring] = (TH1F*) f1->Get(Form("PeakDecoResiduals/histos/hdw_TEC_wheel5_ring%i_stereo0",ring+1));
    hdw2[ring] = (TH1F*) f2->Get(Form("PeakDecoResiduals/histos/hdw_TEC_wheel5_ring%i_stereo0",ring+1));


    hdw1[ring]->SetLineColor(2);
    hdw2[ring]->SetLineColor(4);

    hdw1[ring]->Rebin(5);
    hdw2[ring]->Rebin(5);

    float integral1 = hdw1[ring]->Integral();
    float integral2 = hdw2[ring]->Integral();

    if( integral1 > 0 && integral2 > 0 ){
      hdw1[ring]->Scale( 1. / integral1 );
      hdw2[ring]->Scale( 1. / integral2 );
    }

    //TF1 *g1=new TF1("g1","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",-5000,5000);
    //g1->SetParameters(hdw1[ring]->GetMaximum()/3.,0,500,hdw1[ring]->GetMaximum()/3.,1000,hdw1[ring]->GetMaximum()/3.,2000);
    //g1->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
    
    //TF1 *g2=new TF1("g2","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",-5000,5000);
    //g2->SetParameters(hdw1[ring]->GetMaximum()/3.,0,500,hdw1[ring]->GetMaximum()/3.,1000,hdw1[ring]->GetMaximum()/3.,2000);
    //g2->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
    
    TF1* g1 = new TF1("g1","gaus",-500,500);
    TF1* g2 = new TF1("g2","gaus",-500,500);

    g1->SetLineColor(2);
    g2->SetLineColor(4);

    hdw1[ring]->Fit( g1 , "R" );
    hdw2[ring]->Fit( g2 , "R" );

    if( useFit ){
      mean1[ring]     = g1->GetParameter(1);
      mean2[ring]     = g2->GetParameter(1);
      //mean1err[ring]  = g1->GetParError(1);
      //mean2err[ring]  = g2->GetParError(1);
      mean1err[ring]  = 0.;
      mean2err[ring]  = 0.;
      //mean1err[ring]  = 0.;
      //if( hdw1[ring]->GetEntries() > 0 ) mean1err[ring] = hdw1[ring]->GetRMS(1) / sqrt( hdw1[ring]->GetEntries() ) ;
      //mean2err[ring]  = 0.;
      //if( hdw2[ring]->GetEntries() > 0 ) mean2err[ring] = hdw2[ring]->GetRMS(1) / sqrt( hdw2[ring]->GetEntries() ) ;
    }else{
      mean1[ring]     = hdw1[ring]->GetMean(1);
      mean2[ring]     = hdw2[ring]->GetMean(1);
      mean1err[ring]  = 0.;
      if( hdw1[ring]->GetEntries() > 0 ) mean1err[ring] = hdw1[ring]->GetRMS(1) / sqrt( hdw1[ring]->GetEntries() ) ;
      mean2err[ring]  = 0.;
      if( hdw2[ring]->GetEntries() > 0 ) mean2err[ring] = hdw2[ring]->GetRMS(1) / sqrt( hdw2[ring]->GetEntries() ) ;
      
    }

   
    diff[ring]      = mean2[ring] - mean1[ring];
    differr[ring]   = sqrt( pow( mean1err[ring] , 2 ) + pow( mean2err[ring] , 2 ) );

    cout << "mean1 " << mean1[ring] << " +/- " << mean1err[ring] << endl;   
    cout << "mean2 " << mean2[ring] << " +/- " << mean2err[ring] << endl;
    cout << "diff  " << diff[ring]  << " +/- " << differr[ring]  << endl;

    hdw1[ring]->Draw();
    hdw1[ring]->GetXaxis()->SetNdivisions(5);
    hdw2[ring]->Draw("same");

 

    t.SetTextColor(2);
    t.DrawLatex(0.15,0.8,Form("<#DeltaW> = %.1f",mean1[ring]));

    t.SetTextColor(4);
    t.DrawLatex(0.15,0.7, Form("<#DeltaW> = %.1f",mean2[ring]));
  }

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->cd();

  float ringnum[nrings];
  float ringnumerr[nrings];
  float dwexp[nrings]={-15,-15,-15,-15,-25,-25,-25};
  float dwexperr[nrings];

  for( unsigned int i = 0 ; i < nrings ; ++i ){
    ringnum[i] = i + 1;
    ringnumerr[i] = 0;
    dwexperr[i]=0;
  }



  TGraphErrors *gr    = new TGraphErrors(nrings,ringnum,diff,ringnumerr,differr);
  TGraphErrors *grexp = new TGraphErrors(nrings,ringnum,dwexp,ringnumerr,dwexperr);
  //TGraphErrors *gr = new TGraphErrors(nrings,ringnum,diff);
  gr->GetXaxis()->SetTitle("Ring Number");
  gr->GetYaxis()->SetTitle("#DeltaW shift");
  //gr->SetTitle("TEClateBP = 0.1, wheel5");
  gr->SetTitle("");
  gr->SetTitle("");
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);
  gr->Draw("AP");

  grexp->SetLineColor(4);
  grexp->SetMarkerColor(4);
  grexp->SetMarkerStyle(29);
  grexp->SetMarkerSize(1.5);
  grexp->Draw("sameP");

  gPad->SetGridy(1);


}
