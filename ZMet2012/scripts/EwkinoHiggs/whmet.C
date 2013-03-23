{

  //---------------------------------
  // get total cross section
  //---------------------------------

  TFile* f = TFile::Open("C1N2_8TeV_finer.root");
  TH1F*  h = (TH1F*) f->Get("C1N2_8TeV_NLO");

  //---------------------------------
  // define branching fractions
  //---------------------------------

  float Wlv = 2 * 0.108;
  float Wjj = 1 - (3*0.108);

  // see Table 1 of http://arxiv.org/pdf/1203.3456.pdf
  float Hbb = 0.577;
  float Htt = 0.063;
  float HWW = 0.215;
  float HZZ = 0.026;
  float Hgg = 0.0023;

  //---------------------------------
  // create and scale histograms
  //---------------------------------

  TH1F* h_Wlv_Hbb = (TH1F*) h->Clone("h_Wlv_Hbb");
  TH1F* h_Wlv_Hgg = (TH1F*) h->Clone("h_Wlv_Hgg");
  TH1F* h_Wjj_Hgg = (TH1F*) h->Clone("h_Wjj_Hgg");

  TH1F* h_Wlv_Htt = (TH1F*) h->Clone("h_Wlv_Htt");
  TH1F* h_Wlv_HWW = (TH1F*) h->Clone("h_Wlv_HWW");
  TH1F* h_Wlv_HZZ = (TH1F*) h->Clone("h_Wlv_HZZ");
  
  h_Wlv_Hbb->Scale( Wlv * Hbb );
  h_Wlv_Hgg->Scale( Wlv * Hgg );
  h_Wjj_Hgg->Scale( Wjj * Hgg );

  h_Wlv_Htt->Scale( Wlv * Htt );
  h_Wlv_HWW->Scale( Wlv * HWW );
  h_Wlv_HZZ->Scale( Wlv * HZZ );

  h_Wlv_Hbb->SetLineColor(2);
  h_Wlv_Hgg->SetLineColor(4);
  h_Wjj_Hgg->SetLineColor(7);
  h_Wlv_Htt->SetLineColor(6);
  h_Wlv_HWW->SetLineColor(8);
  h_Wlv_HZZ->SetLineColor(9);

  h->SetLineWidth(4);
  h_Wlv_Hbb->SetLineWidth(2);
  h_Wlv_Hgg->SetLineWidth(2);
  h_Wjj_Hgg->SetLineWidth(2);
  h_Wlv_Htt->SetLineWidth(2);
  h_Wlv_HWW->SetLineWidth(2);
  h_Wlv_HZZ->SetLineWidth(2);

  //---------------------------------
  // make plots
  //---------------------------------

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetLogy();

  h->GetXaxis()->SetTitle("m_{#chi_{1}^{#pm}} = m_{#chi_{2}^{0}} [GeV]");
  h->GetYaxis()->SetTitle("#sigma(pp#rightarrow #chi_{1}^{#pm} #chi_{2}^{0}) #times BF [pb]");
  h->Draw("c");
  h->SetMinimum(0.0001);
  h->SetMaximum(100);
  h->GetXaxis()->SetRangeUser(125,400);
  h_Wlv_Hbb->Draw("samec");
  h_Wlv_Hgg->Draw("samec");
  h_Wjj_Hgg->Draw("samec");

  h_Wlv_Htt->Draw("samec");
  h_Wlv_HWW->Draw("samec");
  h_Wlv_HZZ->Draw("samec");

  TLatex *t = new TLatex();
  t->SetNDC();

  t->SetTextSize(0.05);
  t->DrawLatex(0.2,0.92,"8 TeV NLO cross sections");
  t->DrawLatex(0.3,0.80,"#chi_{1}^{#pm} #rightarrow W #chi_{1}^{0}");
  t->DrawLatex(0.3,0.72,"#chi_{2}^{0} #rightarrow h #chi_{1}^{0}");

  TLegend *leg = new TLegend(0.6,0.57,0.85,0.858);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h,"total","l");
  leg->AddEntry(h_Wlv_Hbb,"W(#font[12]{l}#nu)H(b#bar{b})","l");
  leg->AddEntry(h_Wlv_HWW,"W(#font[12]{l}#nu)H(WW)","l");
  leg->AddEntry(h_Wlv_Htt,"W(#font[12]{l}#nu)H(#tau#tau)","l");
  leg->AddEntry(h_Wlv_HZZ,"W(#font[12]{l}#nu)H(ZZ)","l");
  leg->AddEntry(h_Wjj_Hgg,"W(jj)H(#gamma#gamma)","l");
  leg->AddEntry(h_Wlv_Hgg,"W(#font[12]{l}#nu)H(#gamma#gamma)","l");

  leg->Draw();

  c1->Print("whmet.pdf");


}
