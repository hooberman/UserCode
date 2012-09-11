{
//=========Macro generated from canvas: c1_n7/c1_n7
//=========  (Sun Jan  8 17:13:00 2012) by ROOT version5.27/06b
   TCanvas *c1_n7 = new TCanvas("c1_n7", "c1_n7",0,0,700,500);
   gStyle->SetOptStat(0);
   c1_n7->Range(3.124997,-1.75,471.875,5.75);
   c1_n7->SetFillColor(0);
   c1_n7->SetBorderMode(0);
   c1_n7->SetBorderSize(2);
   c1_n7->SetLogy();
   c1_n7->SetGridx();
   c1_n7->SetGridy();
   c1_n7->SetTickx(1);
   c1_n7->SetTicky(1);
   c1_n7->SetFrameBorderMode(0);
   c1_n7->SetFrameBorderMode(0);
   Double_t xAxis16[11] = {50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 425}; 
   
   TH1F *h_r4_met = new TH1F("h_r4_met","SR4: 125 GeV < H_{T} < 300 GeV",10, xAxis16);
   h_r4_met->SetBinContent(1,1722);
   h_r4_met->SetBinContent(2,2264);
   h_r4_met->SetBinContent(3,1321);
   h_r4_met->SetBinContent(4,623);
   h_r4_met->SetBinContent(5,226);
   h_r4_met->SetBinContent(6,89);
   h_r4_met->SetBinContent(7,34);
   h_r4_met->SetBinContent(8,16);
   h_r4_met->SetBinContent(9,3);
   h_r4_met->SetBinContent(10,1);
   h_r4_met->SetBinError(1,41.49699);
   h_r4_met->SetBinError(2,47.58151);
   h_r4_met->SetBinError(3,36.34556);
   h_r4_met->SetBinError(4,24.95997);
   h_r4_met->SetBinError(5,15.0333);
   h_r4_met->SetBinError(6,9.433981);
   h_r4_met->SetBinError(7,5.830952);
   h_r4_met->SetBinError(8,4);
   h_r4_met->SetBinError(9,1.732051);
   h_r4_met->SetBinError(10,0.4082483);
   h_r4_met->SetMinimum(0.1);
   h_r4_met->SetMaximum(100000);
   h_r4_met->SetEntries(6305);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff0000");
   h_r4_met->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   h_r4_met->SetMarkerColor(ci);
   h_r4_met->SetMarkerStyle(24);
   h_r4_met->GetXaxis()->SetTitle("E_{ T }^{ MISS } ( GeV )");
   h_r4_met->GetYaxis()->SetTitle("Events / 25 GeV");
   h_r4_met->Draw("E1");
   Double_t xAxis17[11] = {50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 425}; 
   
   TH1F *h_r4_pt = new TH1F("h_r4_pt","h_pt",10, xAxis17);
   h_r4_pt->SetBinContent(0,4523.22);
   h_r4_pt->SetBinContent(1,1749.142);
   h_r4_pt->SetBinContent(2,2447.733);
   h_r4_pt->SetBinContent(3,1365.489);
   h_r4_pt->SetBinContent(4,614.1578);
   h_r4_pt->SetBinContent(5,248.0721);
   h_r4_pt->SetBinContent(6,97.12337);
   h_r4_pt->SetBinContent(7,29.2541);
   h_r4_pt->SetBinContent(8,16.65287);
   h_r4_pt->SetBinContent(9,10.78248);
   h_r4_pt->SetBinContent(10,2.051857);
   h_r4_pt->SetBinError(0,92.71696);
   h_r4_pt->SetBinError(1,82.05285);
   h_r4_pt->SetBinError(2,64.6733);
   h_r4_pt->SetBinError(3,43.46319);
   h_r4_pt->SetBinError(4,27.34361);
   h_r4_pt->SetBinError(5,17.9294);
   h_r4_pt->SetBinError(6,12.46249);
   h_r4_pt->SetBinError(7,7.836014);
   h_r4_pt->SetBinError(8,5.375461);
   h_r4_pt->SetBinError(9,5.375461);
   h_r4_pt->SetBinError(10,1.252771);
   h_r4_pt->SetBinError(11,7.516628);
   h_r4_pt->SetEntries(6316);

   ci = TColor::GetColor("#0000ff");
   h_r4_pt->SetFillColor(ci);
   //h_r4_pt->SetFillStyle(3002);

   ci = TColor::GetColor("#0000ff");
   h_r4_pt->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   h_r4_pt->SetMarkerColor(ci);
   h_r4_pt->SetMarkerStyle(23);
   h_r4_pt->Draw("SAME E2");
   Double_t xAxis18[11] = {50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 425}; 
   
   TH1F *h_r4_met = new TH1F("h_r4_met","SR4: 125 GeV < H_{T} < 300 GeV",10, xAxis18);
   h_r4_met->SetBinContent(1,1722);
   h_r4_met->SetBinContent(2,2264);
   h_r4_met->SetBinContent(3,1321);
   h_r4_met->SetBinContent(4,623);
   h_r4_met->SetBinContent(5,226);
   h_r4_met->SetBinContent(6,89);
   h_r4_met->SetBinContent(7,34);
   h_r4_met->SetBinContent(8,16);
   h_r4_met->SetBinContent(9,3);
   h_r4_met->SetBinContent(10,1);
   h_r4_met->SetBinError(1,41.49699);
   h_r4_met->SetBinError(2,47.58151);
   h_r4_met->SetBinError(3,36.34556);
   h_r4_met->SetBinError(4,24.95997);
   h_r4_met->SetBinError(5,15.0333);
   h_r4_met->SetBinError(6,9.433981);
   h_r4_met->SetBinError(7,5.830952);
   h_r4_met->SetBinError(8,4);
   h_r4_met->SetBinError(9,1.732051);
   h_r4_met->SetBinError(10,0.4082483);
   h_r4_met->SetMinimum(0.1);
   h_r4_met->SetMaximum(100000);
   h_r4_met->SetEntries(6305);

   ci = TColor::GetColor("#ff0000");
   h_r4_met->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   h_r4_met->SetMarkerColor(ci);
   h_r4_met->SetMarkerStyle(24);
   h_r4_met->GetXaxis()->SetTitle("E_{ T }^{ MISS } ( GeV )");
   h_r4_met->GetYaxis()->SetTitle("Events / 25 GeV");
   h_r4_met->Draw("E1 SAME");
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
      tex = new TLatex(0.275,0.775,"#scale[0.6]{#int}Ldt = 4.98 fb^{-1}");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   TArrow *arrow = new TArrow(275,0,275,100000,0.05,"<");
   arrow->SetFillColor(1);
   arrow->SetFillStyle(1001);
   arrow->SetLineStyle(9);
   arrow->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.9339831,0.4667816,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   text = pt->AddText("SR4: 125 GeV < H_{T} < 300 GeV");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.8,0.8,0.99,0.99,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("h_r4_met","Observed","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h_r4_pt","Predicted","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1);
   leg->Draw();
   c1_n7->Modified();
   c1_n7->cd();
   c1_n7->SetSelected(c1_n7);


   TH1F* hmet = (TH1F*) h_r4_met->Clone("hmet");
   TH1F* hpt  = (TH1F*) h_r4_pt->Clone("hmet");

   //hmet->SetLineWidth(3);
   //hpt->SetFillColor(7);
   hmet->SetLineWidth(3);
   hpt->SetFillColor(kCyan-9);


   gStyle->SetErrorX(0.5);

   /*
   TCanvas *c1 = new TCanvas();
   c1->cd();
   gPad->SetLogy();
   gPad->SetTopMargin(0.1);

   hpt->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
   hpt->GetYaxis()->SetTitle("events / 25 GeV");
   hpt->SetLineColor(0);
   hpt->SetMaximum(10000);
   hpt->SetMinimum(0.5);
   hpt->Draw("E2");
   hmet->Draw("same");

   TLine line;
   line.SetLineStyle(2);
   line.SetLineWidth(2);
   line.DrawLine(275,0.5,275,10000);

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
   t->DrawLatex(0.68,0.8,"H_{T} 125 - 300 GeV");
*/

   TCanvas *c1 = new TCanvas();
   c1->cd();

   TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextSize(0.04);
   t->DrawLatex(0.18,0.95,"CMS                                #sqrt{s} = 7 TeV,  #scale[0.6]{#int} Ldt = 4.98 fb^{-1}");

   TPad *mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
   mainpad->Draw();
   mainpad->cd();

   mainpad->SetLogy();
   mainpad->SetTopMargin(0.1);
   mainpad->SetRightMargin(0.04);

   // gPad->SetLogy();
   // gPad->SetTopMargin(0.1);
   // gPad->SetRightMargin(0.04);

   hpt->GetYaxis()->SetTitle("events / 25 GeV");
   hpt->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
   hpt->SetLineColor(0);
   hpt->SetMaximum(10000);
   hpt->SetMinimum(0.5);
   hpt->Draw("E2");
   hmet->Draw("same");

   TLine line;
   line.SetLineStyle(2);
   line.SetLineWidth(2);
   //line.DrawLine(200,1,200,1000);
   line.DrawLine(275,0.5,275,10000);

   TLegend *leg = new TLegend(0.2,0.2,0.45,0.4);
   leg->AddEntry(hpt,"predicted","pf");
   leg->AddEntry(hmet,"observed","lp");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->Draw();

   t->DrawLatex(0.7,0.8,"H_{T} 125-300 GeV");

   c1->cd();

   TPad *respad = new TPad("respad","respad",0.0,0.74,1.0,0.92);
   respad->Draw();
   respad->cd();
   respad->SetRightMargin(0.04);
   respad->SetTopMargin(0.1);
   respad->SetGridy();

   TH1F* hratio = (TH1F*) hmet->Clone("hratio");
   hratio->Divide(hpt);
   hratio->GetYaxis()->SetRangeUser(0,2);
   hratio->GetYaxis()->SetNdivisions(5);
   hratio->GetYaxis()->SetLabelSize(0.2);
   hratio->GetXaxis()->SetTitleSize(0.0);
   hratio->GetXaxis()->SetLabelSize(0.0);
   hratio->GetYaxis()->SetTitleSize(0.24);
   hratio->GetYaxis()->SetTitleOffset(0.3);
   hratio->GetYaxis()->SetTitle("ratio");
   //hratio->SetMarkerColor(1);
   //hratio->SetLineColor(1);
   hratio->Draw();

   TLine line;
   line.SetLineWidth(2);
   line.DrawLine(50,1,425,1);
   hratio->Draw("same");


   cout << "met : " << 6*hmet->GetBinContent(10) << " +/- " << 6*hmet->GetBinError(10) << endl;
   cout << "met : " << 6*hpt->GetBinContent(10)  << " +/- " << 6*hpt->GetBinError(10)  << endl;

   c1->Print("ptll_SR4_498fb.pdf");

}

