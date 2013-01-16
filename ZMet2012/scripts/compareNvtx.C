{

  TChain *data = new TChain("T1");
  data->Add("../output/V00-02-02/data_53X_2012A_baby.root");
  data->Add("../output/V00-02-02/data_53X_2012B_baby.root");
  data->Add("../output/V00-02-02/data_53X_2012C_baby.root");
  data->Add("../output/V00-02-02/data_53X_2012D_baby.root");

  TChain *mc   = new TChain("T1");
  mc->Add("../output/V00-02-02/zjets_53X_slim_baby.root");
  //mc->Add("../output/V00-00-12/ttbar_baby.root");
  
  TH1F* hdata        = new TH1F("hdata"        ,"",40,0,40);
  TH1F* hmc          = new TH1F("hmc"          ,"",40,0,40);
  TH1F* hmc_weighted = new TH1F("hmc_weighted" ,"",40,0,40);

  data->Draw("nvtx>>hdata");
  mc->Draw("nvtx>>hmc");
  mc->Draw("nvtx>>hmc_weighted","vtxweight");

  cout << "MC (unweighted) " << hmc->Integral()                            << endl;
  cout << "MC (weighted)   " << hmc_weighted->Integral()                   << endl;
  cout << "MC (ratio)      " << hmc_weighted->Integral() / hmc->Integral() << endl;
  
  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.05);

  hdata->GetXaxis()->SetTitle("N_{VTX}");
  hdata->Sumw2();
  hdata->DrawNormalized("E1");
  hmc->SetLineColor(2);
  hmc_weighted->SetLineColor(4);
  hmc->DrawNormalized("samehist");
  hmc_weighted->DrawNormalized("samehist");

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(hdata,"data","lp");
  leg->AddEntry(hmc,"MC","l");
  leg->AddEntry(hmc_weighted,"MC (weighted)","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();




}
