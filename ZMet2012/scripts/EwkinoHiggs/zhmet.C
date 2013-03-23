{

  //---------------------------------
  // get total cross section
  //---------------------------------

  TH1F*  h = new TH1F("htot","",16,100,420);

  h->SetBinContent(1 ,  7.288);
  h->SetBinContent(2 ,  3.764);
  h->SetBinContent(3 ,  2.141);
  h->SetBinContent(4 ,  1.304);
  h->SetBinContent(5 ,  0.837);
  h->SetBinContent(6 ,  0.558);
  h->SetBinContent(7 ,  0.382);
  h->SetBinContent(8 ,  0.271);
  h->SetBinContent(9 ,  0.195);
  h->SetBinContent(10,  0.142);
  h->SetBinContent(11,  0.106);
  h->SetBinContent(12,  0.0798);
  h->SetBinContent(13,  0.0608);
  h->SetBinContent(14,  0.0468);
  h->SetBinContent(15,  0.0366);
  h->SetBinContent(16,  0.0287);

  h->Scale(0.5);

  //---------------------------------
  // define branching fractions
  //---------------------------------

  float Zll = 0.066;
  float Zvv = 0.200;
  float Zjj = 0.700;

  // see Table 1 of http://arxiv.org/pdf/1203.3456.pdf
  float Hbb = 0.577;
  float Htt = 0.063;
  float HWW = 0.215;
  float HZZ = 0.026;
  float Hgg = 0.0023;

  //---------------------------------
  // create and scale histograms
  //---------------------------------

  TH1F* h_Zll_Hbb = (TH1F*) h->Clone("h_Zll_Hbb");
  TH1F* h_Zll_Hgg = (TH1F*) h->Clone("h_Zll_Hgg");
  TH1F* h_Zjj_Hgg = (TH1F*) h->Clone("h_Zjj_Hgg");
  TH1F* h_Zvv_Hgg = (TH1F*) h->Clone("h_Zvv_Hgg");

  TH1F* h_Zll_Htt = (TH1F*) h->Clone("h_Zll_Htt");
  TH1F* h_Zll_HWW = (TH1F*) h->Clone("h_Zll_HWW");
  TH1F* h_Zll_HZZ = (TH1F*) h->Clone("h_Zll_HZZ");
  
  h_Zll_Hbb->Scale( Zll * Hbb );
  h_Zll_Hgg->Scale( Zll * Hgg );
  h_Zjj_Hgg->Scale( Zjj * Hgg );
  h_Zvv_Hgg->Scale( Zvv * Hgg );

  h_Zll_Htt->Scale( Zll * Htt );
  h_Zll_HWW->Scale( Zll * HWW );
  h_Zll_HZZ->Scale( Zll * HZZ );

  h_Zll_Hbb->SetLineColor(2);
  h_Zll_Hgg->SetLineColor(4);
  h_Zjj_Hgg->SetLineColor(kOrange);
  h_Zvv_Hgg->SetLineColor(7);
  h_Zll_Htt->SetLineColor(6);
  h_Zll_HWW->SetLineColor(8);
  h_Zll_HZZ->SetLineColor(9);

  h->SetLineWidth(4);
  h_Zll_Hbb->SetLineWidth(2);
  h_Zll_Hgg->SetLineWidth(2);
  h_Zjj_Hgg->SetLineWidth(2);
  h_Zvv_Hgg->SetLineWidth(2);
  h_Zll_Htt->SetLineWidth(2);
  h_Zll_HWW->SetLineWidth(2);
  h_Zll_HZZ->SetLineWidth(2);

  //---------------------------------
  // make plots
  //---------------------------------

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();

  gPad->SetTopMargin(0.1);
  gPad->SetLogy();

  h->GetXaxis()->SetTitle("higgsino mass parameter #mu [GeV]");
  h->GetYaxis()->SetTitle("(1/2) #times #sigma(GMSB) #times BF [pb]");
  h->Draw("c");
  h->SetMinimum(0.0001);
  h->SetMaximum(100);
  h->GetXaxis()->SetRangeUser(125,400);
  h_Zll_Hbb->Draw("samec");
  h_Zll_Hgg->Draw("samec");
  h_Zjj_Hgg->Draw("samec");
  h_Zvv_Hgg->Draw("samec");

  h_Zll_Htt->Draw("samec");
  h_Zll_HWW->Draw("samec");
  h_Zll_HZZ->Draw("samec");

  TLatex *t = new TLatex();
  t->SetNDC();

  t->SetTextSize(0.05);
  t->DrawLatex(0.2,0.92,"8 TeV NLO cross sections");
  //t->DrawLatex(0.3,0.80,"#chi_{1}^{#pm} #rightarrow W #chi_{1}^{0}");
  //t->DrawLatex(0.3,0.72,"#chi_{2}^{0} #rightarrow h #chi_{1}^{0}");

  TLegend *leg = new TLegend(0.6,0.57,0.85,0.858);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h,"total","l");
  leg->AddEntry(h_Zll_Hbb,"Z(#font[12]{ll})H(b#bar{b})","l");
  leg->AddEntry(h_Zll_HWW,"Z(#font[12]{ll})H(WW)","l");
  leg->AddEntry(h_Zll_Htt,"Z(#font[12]{ll})H(#tau#tau)","l");
  leg->AddEntry(h_Zll_HZZ,"Z(#font[12]{ll})H(ZZ)","l");
  leg->AddEntry(h_Zjj_Hgg,"Z(jj)H(#gamma#gamma)","l");
  leg->AddEntry(h_Zvv_Hgg,"Z(#nu#nu)H(#gamma#gamma)","l");
  leg->AddEntry(h_Zll_Hgg,"Z(#font[12]{ll})H(#gamma#gamma)","l");

  leg->Draw();

  c1->Print("zhmet.pdf");
  

}
