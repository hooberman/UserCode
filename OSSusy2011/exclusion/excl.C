{

  gROOT->ProcessLine(".L OSstuff.hh");

  TFile *f = TFile::Open("exclusion_Spring11.root"); // new scan
  //TFile *f = TFile::Open("exclusion_Fall10_tcmet_JPT.root"); // new scan
  //TFile *f = TFile::Open("exclusion_Fall10_pfmet_pfjets.root"); // new scan


  
  TH2F* hobs = (TH2F*) f->Get("hexcl_NLO_obs");
  TH2F* hexp = (TH2F*) f->Get("hexcl_NLO_exp");

  /*
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  gPad->SetGridx();
  gPad->SetGridy();

  hexp->GetXaxis()->SetRangeUser(0,1500);
  hexp->GetXaxis()->SetTitle("m_{0} [GeV]");
  hexp->GetYaxis()->SetTitle("m_{1/2} [GeV]");
				 
  hexp->Draw("colz");


  TGraphErrors *gexp = getNLOexpTanbeta10();
  gexp->SetMarkerColor(4);
  gexp->Draw("sameP");
  gexp->Draw("samec");
  */


  TCanvas *c2 = new TCanvas("c2","",1200,600);
  gPad->SetGridx();
  gPad->SetGridy();

  hobs->GetXaxis()->SetRangeUser(0,1500);
  //hobs->GetXaxis()->SetRangeUser(0,1000);
  //hobs->GetYaxis()->SetRangeUser(100,450);
  hobs->GetXaxis()->SetTitle("m_{0} [GeV]");
  hobs->GetYaxis()->SetTitle("m_{1/2} [GeV]");

  hobs->Draw("colz");

  //TGraphErrors *gobs = getNLOobsTanbeta10_smooth();
  //TGraphErrors *gobs = getNLOobsTanbeta10_funky();
  TGraphErrors *gobs = getNLOobsTanbeta10_intermediate();
  //TGraphErrors *gobs = getNLOobsTanbeta10();
  gobs->SetMarkerColor(4);
  gobs->Draw("sameP");
  gobs->Draw("samec");




  //TGraphErrors *gr = getNLOobsTanbeta10();





}
