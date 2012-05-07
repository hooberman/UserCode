{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Wed Dec  7 16:51:06 2011) by ROOT version5.27/06b
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",0,0,700,500);
   gStyle->SetOptStat(0);
   c1_n2->Range(3.124997,-1.625,471.875,4.625);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetLogy();
   c1_n2->SetGridx();
   c1_n2->SetGridy();
   c1_n2->SetFrameBorderMode(0);
   c1_n2->SetFrameBorderMode(0);
   Double_t xAxis1[11] = {50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 425}; 
   
   TH1F *h_r1r2_met = new TH1F("h_r1r2_met","High E_{T}^{MISS}: 300 GeV < H_{T}",10, xAxis1);
   h_r1r2_met->SetBinContent(1,234);
   h_r1r2_met->SetBinContent(2,428);
   h_r1r2_met->SetBinContent(3,259);
   h_r1r2_met->SetBinContent(4,163);
   h_r1r2_met->SetBinContent(5,106);
   h_r1r2_met->SetBinContent(6,51);
   h_r1r2_met->SetBinContent(7,41);
   h_r1r2_met->SetBinContent(8,24);
   h_r1r2_met->SetBinContent(9,14);
   h_r1r2_met->SetBinContent(10,5);
   h_r1r2_met->SetBinError(1,15.29706);
   h_r1r2_met->SetBinError(2,20.68816);
   h_r1r2_met->SetBinError(3,16.09348);
   h_r1r2_met->SetBinError(4,12.76715);
   h_r1r2_met->SetBinError(5,10.29563);
   h_r1r2_met->SetBinError(6,7.141428);
   h_r1r2_met->SetBinError(7,6.403124);
   h_r1r2_met->SetBinError(8,4.898979);
   h_r1r2_met->SetBinError(9,3.741657);
   h_r1r2_met->SetBinError(10,0.9128709);
   h_r1r2_met->SetMinimum(0.1);
   h_r1r2_met->SetMaximum(10000);
   h_r1r2_met->SetEntries(1351);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff0000");
   h_r1r2_met->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   h_r1r2_met->SetMarkerColor(ci);
   h_r1r2_met->SetMarkerStyle(24);
   h_r1r2_met->GetXaxis()->SetTitle("E_{ T }^{ MISS } ( GeV )");
   h_r1r2_met->GetYaxis()->SetTitle("Events / 25 GeV");
   h_r1r2_met->Draw("E1");
   Double_t xAxis2[11] = {50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 425}; 
   
   TH1F *h_r1r2_pt = new TH1F("h_r1r2_pt","h_pt",10, xAxis2);
   h_r1r2_pt->SetBinContent(0,690.6212);
   h_r1r2_pt->SetBinContent(1,242.4681);
   h_r1r2_pt->SetBinContent(2,389.7155);
   h_r1r2_pt->SetBinContent(3,263.0812);
   h_r1r2_pt->SetBinContent(4,150.7529);
   h_r1r2_pt->SetBinContent(5,123.7933);
   h_r1r2_pt->SetBinContent(6,81.33868);
   h_r1r2_pt->SetBinContent(7,41.67567);
   h_r1r2_pt->SetBinContent(8,22.00697);
   h_r1r2_pt->SetBinContent(9,15.79362);
   h_r1r2_pt->SetBinContent(10,3.53613);
   h_r1r2_pt->SetBinError(0,33.69886);
   h_r1r2_pt->SetBinError(1,28.99819);
   h_r1r2_pt->SetBinError(2,25.52696);
   h_r1r2_pt->SetBinError(3,19.93652);
   h_r1r2_pt->SetBinError(4,14.70739);
   h_r1r2_pt->SetBinError(5,12.84268);
   h_r1r2_pt->SetBinError(6,10.13637);
   h_r1r2_pt->SetBinError(7,7.712621);
   h_r1r2_pt->SetBinError(8,5.453646);
   h_r1r2_pt->SetBinError(9,4.650885);
   h_r1r2_pt->SetBinError(10,1.995671);
   h_r1r2_pt->SetBinError(11,11.97402);
   h_r1r2_pt->SetEntries(1362);

   ci = TColor::GetColor("#0000ff");
   h_r1r2_pt->SetFillColor(ci);
   //h_r1r2_pt->SetFillStyle(3002);

   ci = TColor::GetColor("#0000ff");
   h_r1r2_pt->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   h_r1r2_pt->SetMarkerColor(ci);
   h_r1r2_pt->SetMarkerStyle(23);
   h_r1r2_pt->Draw("SAME E2");
   Double_t xAxis3[11] = {50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 425}; 
   
   TH1F *h_r1r2_met = new TH1F("h_r1r2_met","High E_{T}^{MISS}: 300 GeV < H_{T}",10, xAxis3);
   h_r1r2_met->SetBinContent(1,234);
   h_r1r2_met->SetBinContent(2,428);
   h_r1r2_met->SetBinContent(3,259);
   h_r1r2_met->SetBinContent(4,163);
   h_r1r2_met->SetBinContent(5,106);
   h_r1r2_met->SetBinContent(6,51);
   h_r1r2_met->SetBinContent(7,41);
   h_r1r2_met->SetBinContent(8,24);
   h_r1r2_met->SetBinContent(9,14);
   h_r1r2_met->SetBinContent(10,5);
   h_r1r2_met->SetBinError(1,15.29706);
   h_r1r2_met->SetBinError(2,20.68816);
   h_r1r2_met->SetBinError(3,16.09348);
   h_r1r2_met->SetBinError(4,12.76715);
   h_r1r2_met->SetBinError(5,10.29563);
   h_r1r2_met->SetBinError(6,7.141428);
   h_r1r2_met->SetBinError(7,6.403124);
   h_r1r2_met->SetBinError(8,4.898979);
   h_r1r2_met->SetBinError(9,3.741657);
   h_r1r2_met->SetBinError(10,0.9128709);
   h_r1r2_met->SetMinimum(0.1);
   h_r1r2_met->SetMaximum(10000);
   h_r1r2_met->SetEntries(1351);

   ci = TColor::GetColor("#ff0000");
   h_r1r2_met->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   h_r1r2_met->SetMarkerColor(ci);
   h_r1r2_met->SetMarkerStyle(24);
   h_r1r2_met->GetXaxis()->SetTitle("E_{ T }^{ MISS } ( GeV )");
   h_r1r2_met->GetYaxis()->SetTitle("Events / 25 GeV");
   h_r1r2_met->Draw("E1 SAME");
   TLatex *   tex = new TLatex(0.15,0.85,"CMS Preliminary");
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
   TArrow *arrow = new TArrow(275,0,275,10000,0.05,"<");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineStyle(9);
   arrow->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.9111017,0.3949425,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   text = pt->AddText("High E_{T}^{MISS}: 300 GeV < H_{T}");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.8,0.8,0.99,0.99,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("h_r1r2_met","Observed","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h_r1r2_pt","Predicted","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1);
   leg->Draw();
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);


   TH1F* hmet = (TH1F*) h_r1r2_met->Clone("hmet");
   TH1F* hpt  = (TH1F*) h_r1r2_pt->Clone("hmet");

   //hmet->SetLineWidth(3);
   //hpt->SetFillColor(7);

   hmet->SetLineWidth(3);
   hpt->SetFillColor(kCyan-9);


   gStyle->SetErrorX(0.5);

   TCanvas *c1 = new TCanvas();
   c1->cd();
   gPad->SetLogy();
   gPad->SetTopMargin(0.1);

   hpt->GetYaxis()->SetTitle("events / 25 GeV");
   hpt->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
   hpt->SetLineColor(0);
   hpt->SetMaximum(1000);
   hpt->SetMinimum(1);
   hpt->Draw("E2");
   hmet->Draw("same");

   TLine line;
   line.SetLineStyle(2);
   line.SetLineWidth(2);
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
   t->DrawLatex(0.18,0.93,"CMS Preliminary          #sqrt{s} = 7 TeV,  #scale[0.6]{#int} Ldt = 4.98 fb^{-1}");
   t->DrawLatex(0.7,0.8,"H_{T} > 300 GeV");


   c1->Print("ptll_HighMet_498fb.pdf");


}
