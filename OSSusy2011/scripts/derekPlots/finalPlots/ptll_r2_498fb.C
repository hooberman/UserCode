{
//=========Macro generated from canvas: c1_n5/c1_n5
//=========  (Wed Dec  7 16:51:06 2011) by ROOT version5.27/06b
   TCanvas *c1_n5 = new TCanvas("c1_n5", "c1_n5",0,0,700,500);
   gStyle->SetOptStat(0);
   c1_n5->Range(3.124997,-1.5,471.875,3.5);
   c1_n5->SetFillColor(0);
   c1_n5->SetBorderMode(0);
   c1_n5->SetBorderSize(2);
   c1_n5->SetLogy();
   c1_n5->SetGridx();
   c1_n5->SetGridy();
   c1_n5->SetFrameBorderMode(0);
   c1_n5->SetFrameBorderMode(0);
   Double_t xAxis10[5] = {50, 125, 200, 275, 425}; 
   
   TH1F *h_r2_met = new TH1F("h_r2_met","SR2: 600 GeV < H_{T}",4, xAxis10);
   h_r2_met->SetBinContent(1,55);
   h_r2_met->SetBinContent(2,28);
   h_r2_met->SetBinContent(3,18);
   h_r2_met->SetBinContent(4,5.5);
   h_r2_met->SetBinError(1,7.416198);
   h_r2_met->SetBinError(2,5.291503);
   h_r2_met->SetBinError(3,4.242641);
   h_r2_met->SetBinError(4,1.658312);
   h_r2_met->SetMinimum(0.1);
   h_r2_met->SetMaximum(1000);
   h_r2_met->SetEntries(113);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff0000");
   h_r2_met->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   h_r2_met->SetMarkerColor(ci);
   h_r2_met->SetMarkerStyle(24);
   h_r2_met->GetXaxis()->SetTitle("E_{ T }^{ MISS } ( GeV )");
   h_r2_met->GetYaxis()->SetTitle("Events / 25 GeV");
   h_r2_met->Draw("E1");
   Double_t xAxis11[5] = {50, 125, 200, 275, 425}; 
   
   TH1F *h_r2_pt = new TH1F("h_r2_pt","h_pt",4, xAxis11);
   h_r2_pt->SetBinContent(0,51.33334);
   h_r2_pt->SetBinContent(1,46.25075);
   h_r2_pt->SetBinContent(2,40.63969);
   h_r2_pt->SetBinContent(3,14.04667);
   h_r2_pt->SetBinContent(4,5.280644);
   h_r2_pt->SetBinError(0,8.935986);
   h_r2_pt->SetBinError(1,10.08115);
   h_r2_pt->SetBinError(2,7.777778);
   h_r2_pt->SetBinError(3,4.115613);
   h_r2_pt->SetBinError(4,3.465855);
   h_r2_pt->SetBinError(5,6.93171);
   h_r2_pt->SetEntries(118);

   ci = TColor::GetColor("#0000ff");
   h_r2_pt->SetFillColor(ci);
   //h_r2_pt->SetFillStyle(3002);

   ci = TColor::GetColor("#0000ff");
   h_r2_pt->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   h_r2_pt->SetMarkerColor(ci);
   h_r2_pt->SetMarkerStyle(23);
   h_r2_pt->Draw("SAME E2");
   Double_t xAxis12[5] = {50, 125, 200, 275, 425}; 
   
   TH1F *h_r2_met = new TH1F("h_r2_met","SR2: 600 GeV < H_{T}",4, xAxis12);
   h_r2_met->SetBinContent(1,55);
   h_r2_met->SetBinContent(2,28);
   h_r2_met->SetBinContent(3,18);
   h_r2_met->SetBinContent(4,5.5);
   h_r2_met->SetBinError(1,7.416198);
   h_r2_met->SetBinError(2,5.291503);
   h_r2_met->SetBinError(3,4.242641);
   h_r2_met->SetBinError(4,1.658312);
   h_r2_met->SetMinimum(0.1);
   h_r2_met->SetMaximum(1000);
   h_r2_met->SetEntries(113);

   ci = TColor::GetColor("#ff0000");
   h_r2_met->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   h_r2_met->SetMarkerColor(ci);
   h_r2_met->SetMarkerStyle(24);
   h_r2_met->GetXaxis()->SetTitle("E_{ T }^{ MISS } ( GeV )");
   h_r2_met->GetYaxis()->SetTitle("Events / 25 GeV");
   h_r2_met->Draw("E1 SAME");
   TLatex *   tex = new TLatex(0.15,0.85,"CMS");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.775,"#sqrt{s} = 7 TeV ,");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.275,0.775,"#scale[0.6]{#int}Ldt = 4.7 fb^{-1}");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   TArrow *arrow = new TArrow(275,0,275,1000,0.05,"<");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineStyle(9);
   arrow->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.9339831,0.3058621,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   text = pt->AddText("SR2: 600 GeV < H_{T}");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.8,0.8,0.99,0.99,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("h_r2_met","Observed","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h_r2_pt","Predicted","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1);
   leg->Draw();
   c1_n5->Modified();
   c1_n5->cd();
   c1_n5->SetSelected(c1_n5);

   TH1F* hmet = (TH1F*) h_r2_met->Clone("hmet");
   TH1F* hpt  = (TH1F*) h_r2_pt->Clone("hmet");

   hmet->SetLineWidth(3);
   hpt->SetFillColor(kCyan-9);

   gStyle->SetErrorX(0.5);

   TCanvas *c1 = new TCanvas();
   c1->cd();
   gPad->SetLogy();
   gPad->SetTopMargin(0.1);

   hpt->GetYaxis()->SetTitle("events / 75 GeV");
   hpt->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
   hpt->SetLineColor(0);
   hpt->SetMaximum(1000);
   hpt->SetMinimum(1);
   hpt->Draw("E2");
   hmet->Draw("same");

   TLine line;
   line.SetLineStyle(2);
   line.SetLineWidth(2);
   line.DrawLine(200,1,200,1000);
   line.DrawLine(275,1,275,1000);

   TLegend *leg = new TLegend(0.2,0.2,0.45,0.4);
   leg->AddEntry(hpt,"predicted","pf");
   leg->AddEntry(hmet,"observed","lp");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->Draw();

   TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextSize(0.04);
   t->DrawLatex(0.18,0.93,"CMS                                #sqrt{s} = 7 TeV,  #scale[0.6]{#int} Ldt = 4.98 fb^{-1}");
   t->DrawLatex(0.7,0.8,"H_{T} > 600 GeV");


   c1->Print("ptll_SR2_498fb.pdf");


}
