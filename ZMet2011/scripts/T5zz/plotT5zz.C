{

  TFile *f = TFile::Open("histos.root");

  TH2F* h_T5zz     = (TH2F*) f->Get("h_T5zz");
  TH2F* h_T5zzl    = (TH2F*) f->Get("h_T5zzl");
  TH2F* h_T5zzgmsb = (TH2F*) f->Get("h_T5zzgmsb");

  TGraph* g_T5zz     = (TGraph*) f->Get("g_T5zz");
  TGraph* g_T5zzl    = (TGraph*) f->Get("g_T5zzl");
  TGraph* g_T5zzgmsb = (TGraph*) f->Get("g_T5zzgmsb");

  TCanvas* c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_T5zz->SetMinimum(0.001);
  h_T5zz->SetMaximum(10);
  h_T5zz->Draw("colz");
  g_T5zz->Draw("same");

  TCanvas* c2 = new TCanvas();
  c2->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_T5zzl->SetMinimum(0.001);
  h_T5zzl->SetMaximum(10);
  h_T5zzl->Draw("colz");
  g_T5zzl->Draw("same");

  TCanvas* c3 = new TCanvas();
  c3->cd();
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_T5zzgmsb->SetMinimum(0.001);
  h_T5zzgmsb->SetMaximum(10);
  h_T5zzgmsb->Draw("colz");
  g_T5zzgmsb->Draw("same");






}
