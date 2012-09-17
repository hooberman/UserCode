{
//=========Macro generated from canvas: canv_pt/Ratios for pt
//=========  (Wed Jan 25 13:39:24 2012) by ROOT version5.27/06b
   TCanvas *canv_pt = new TCanvas("canv_pt", "Ratios for pt",62,64,500,500);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canv_pt->SetHighLightColor(2);
   canv_pt->Range(-0.8780477,-0.1004286,79.60976,0.672099);
   canv_pt->SetFillColor(0);
   canv_pt->SetBorderMode(0);
   canv_pt->SetBorderSize(2);
   canv_pt->SetTickx(1);
   canv_pt->SetTicky(1);
   canv_pt->SetLeftMargin(0.16);
   canv_pt->SetRightMargin(0.02);
   canv_pt->SetTopMargin(0.05);
   canv_pt->SetBottomMargin(0.13);
   canv_pt->SetFrameFillStyle(0);
   canv_pt->SetFrameBorderMode(0);
   canv_pt->SetFrameFillStyle(0);
   canv_pt->SetFrameBorderMode(0);
   
   TH1F *hframe__1 = new TH1F("hframe__1","pt",1000,12,78);
   hframe__1->SetMinimum(0);
   hframe__1->SetMaximum(0.6334726);
   hframe__1->SetDirectory(0);
   hframe__1->SetStats(0);
   hframe__1->SetLineStyle(0);
   hframe__1->SetMarkerStyle(20);
   hframe__1->GetXaxis()->SetTitle("p^{#tau_{h}, vis}_{T} [GeV]");
   hframe__1->GetXaxis()->SetLabelFont(42);
   hframe__1->GetXaxis()->SetLabelOffset(0.007);
   hframe__1->GetXaxis()->SetLabelSize(0.05);
   hframe__1->GetXaxis()->SetTitleSize(0.06);
   hframe__1->GetXaxis()->SetTitleOffset(0.9);
   hframe__1->GetXaxis()->SetTitleFont(42);
   hframe__1->GetYaxis()->SetTitle("Efficiency");
   hframe__1->GetYaxis()->SetLabelFont(42);
   hframe__1->GetYaxis()->SetLabelOffset(0.007);
   hframe__1->GetYaxis()->SetLabelSize(0.05);
   hframe__1->GetYaxis()->SetTitleSize(0.06);
   hframe__1->GetYaxis()->SetTitleOffset(1.25);
   hframe__1->GetYaxis()->SetTitleFont(42);
   hframe__1->GetZaxis()->SetLabelFont(42);
   hframe__1->GetZaxis()->SetLabelOffset(0.007);
   hframe__1->GetZaxis()->SetLabelSize(0.05);
   hframe__1->GetZaxis()->SetTitleSize(0.06);
   hframe__1->GetZaxis()->SetTitleFont(42);
   hframe__1->Draw(" ");
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("pt");
   multigraph->SetTitle("pt;p^{#tau_{h}, vis}_{T} [GeV];Efficiency");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(13);
   grae->SetName("hpsTau-pt-Efficiency-LM13Truth_");
   grae->SetTitle("#splitline{hpsTau-pt-Efficiency-LM13Truth}{}");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);
   grae->SetLineColor(2);
   grae->SetMarkerColor(2);
   grae->SetMarkerStyle(21);
   grae->SetPoint(0,16.25,0.376415);
   grae->SetPointError(0,1.25,1.25,0.008632,0.008659);
   grae->SetPoint(1,18.75,0.390507);
   grae->SetPointError(1,1.25,1.25,0.009157,0.00917);
   grae->SetPoint(2,21.25,0.422905);
   grae->SetPointError(2,1.25,1.25,0.010548,0.010548);
   grae->SetPoint(3,23.75,0.438852);
   grae->SetPointError(3,1.25,1.25,0.011388,0.011368);
   grae->SetPoint(4,26.25,0.416334);
   grae->SetPointError(4,1.25,1.25,0.011586,0.011585);
   grae->SetPoint(5,28.75,0.400166);
   grae->SetPointError(5,1.25,1.25,0.012535,0.012556);
   grae->SetPoint(6,31.25,0.448596);
   grae->SetPointError(6,1.25,1.25,0.014393,0.014366);
   grae->SetPoint(7,33.75,0.464301);
   grae->SetPointError(7,1.25,1.25,0.017248,0.017295);
   grae->SetPoint(8,37.5,0.449827);
   grae->SetPointError(8,2.5,2.5,0.012574,0.012588);
   grae->SetPoint(9,42.5,0.43194);
   grae->SetPointError(9,2.5,2.5,0.014432,0.0145);
   grae->SetPoint(10,47.5,0.451513);
   grae->SetPointError(10,2.5,2.5,0.016759,0.016776);
   grae->SetPoint(11,55,0.452391);
   grae->SetPointError(11,5,5,0.013218,0.013248);
   grae->SetPoint(12,67.5,0.413535);
   grae->SetPointError(12,7.5,7.5,0.014262,0.014374);
   
   TH1F *Graph1 = new TH1F("Graph1","#splitline{hpsTau-pt-Efficiency-LM13Truth}{}",100,9,81);
   Graph1->SetMinimum(0.3564017);
   Graph1->SetMaximum(0.4929773);
   Graph1->SetDirectory(0);
   Graph1->SetStats(0);
   Graph1->SetLineStyle(0);
   Graph1->SetMarkerStyle(20);
   Graph1->GetXaxis()->SetLabelFont(42);
   Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph1->GetXaxis()->SetTitleOffset(0.9);
   Graph1->GetXaxis()->SetTitleFont(42);
   Graph1->GetYaxis()->SetLabelFont(42);
   Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph1->GetYaxis()->SetTitleOffset(1.25);
   Graph1->GetYaxis()->SetTitleFont(42);
   Graph1->GetZaxis()->SetLabelFont(42);
   Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph1);
   
   
   TF1 *correction_pt = new TF1("correction_pt","[1]+0.5*[2]*(1-[1])*(TMath::Erf((x)/[0])+1)",9,81);
   correction_pt->SetFillColor(19);
   correction_pt->SetFillStyle(0);
   correction_pt->SetMarkerStyle(20);
   correction_pt->SetLineColor(2);
   correction_pt->SetLineWidth(1);
   correction_pt->SetChisquare(20.83675);
   correction_pt->SetNDF(10);
   correction_pt->GetXaxis()->SetLabelFont(42);
   correction_pt->GetXaxis()->SetLabelOffset(0.007);
   correction_pt->GetXaxis()->SetLabelSize(0.05);
   correction_pt->GetXaxis()->SetTitleSize(0.06);
   correction_pt->GetXaxis()->SetTitleOffset(0.9);
   correction_pt->GetXaxis()->SetTitleFont(42);
   correction_pt->GetYaxis()->SetLabelFont(42);
   correction_pt->GetYaxis()->SetLabelOffset(0.007);
   correction_pt->GetYaxis()->SetLabelSize(0.05);
   correction_pt->GetYaxis()->SetTitleSize(0.06);
   correction_pt->GetYaxis()->SetTitleOffset(1.25);
   correction_pt->GetYaxis()->SetTitleFont(42);
   correction_pt->SetParameter(0,16.69773);
   correction_pt->SetParError(0,3.645578);
   correction_pt->SetParLimits(0,0,0);
   correction_pt->SetParameter(1,-0.3176219);
   correction_pt->SetParError(1,0.4317906);
   correction_pt->SetParLimits(1,0,0);
   correction_pt->SetParameter(2,0.5747245);
   correction_pt->SetParError(2,0.1371748);
   correction_pt->SetParLimits(2,0,0);
   grae->GetListOfFunctions()->Add(correction_pt);
   
   // TPaveStats *ptstats = new TPaveStats(0.71,0.87,0.98,0.995,"brNDC");
   // ptstats->SetName("stats");
   // ptstats->SetBorderSize(1);
   // ptstats->SetTextAlign(12);
   // ptstats->SetTextFont(42);
   // TText *text = ptstats->AddText("#chi^{2} / ndf = 20.84 / 10");
   // text = ptstats->AddText("p0       =  16.7 #pm 3.646 ");
   // text = ptstats->AddText("p1       = -0.3176 #pm 0.4318 ");
   // text = ptstats->AddText("p2       = 0.5747 #pm 0.1372 ");
   // ptstats->SetOptStat(0);
   // ptstats->SetOptFit(0);
   // ptstats->Draw();
   // grae->GetListOfFunctions()->Add(ptstats);
   // ptstats->SetParent(grae->GetListOfFunctions());
   multigraph->Add(grae,"");
   multigraph->Draw("P");
   // multigraph->GetXaxis()->SetTitle("p^{#tau_{h}, vis}_{T} [GeV]");
   // multigraph->GetXaxis()->SetLabelFont(42);
   // multigraph->GetXaxis()->SetLabelOffset(0.007);
   // multigraph->GetXaxis()->SetLabelSize(0.05);
   // multigraph->GetXaxis()->SetTitleSize(0.06);
   // multigraph->GetXaxis()->SetTitleOffset(0.9);
   // multigraph->GetXaxis()->SetTitleFont(42);
   // multigraph->GetYaxis()->SetTitle("Efficiency");
   // multigraph->GetYaxis()->SetLabelFont(42);
   // multigraph->GetYaxis()->SetLabelOffset(0.007);
   // multigraph->GetYaxis()->SetLabelSize(0.05);
   // multigraph->GetYaxis()->SetTitleSize(0.06);
   // multigraph->GetYaxis()->SetTitleOffset(1.25);
   // multigraph->GetYaxis()->SetTitleFont(42);
   canv_pt->Modified();
   canv_pt->cd();
   canv_pt->SetSelected(canv_pt);


   TCanvas *c2 = new TCanvas();
   c2->cd();

   TGraphAsymmErrors *gr = (TGraphAsymmErrors*) grae->Clone();

   TF1 *fit = new TF1("fit","[0]*TMath::Erf((x-10.0)/[2])+[1]*(1.0-TMath::Erf((x-10.0)/[2]))",10,250);
   fit->SetParameter(0,1);
   fit->SetParameter(1,0.2);
   fit->SetParameter(2,20);

   fit->SetLineColor(2);
   fit->SetLineWidth(2);
   gr->Fit(fit);

   gPad->SetRightMargin(0.1);
   gPad->SetTopMargin(0.1);
   gPad->SetGridx();
   gPad->SetGridy();

   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);

   gr->GetXaxis()->SetRangeUser(15,75);
   gr->SetMinimum(0.2);
   gr->SetMaximum(0.5);
   gr->SetMarkerStyle(20);
   gr->GetXaxis()->SetTitle("generated visible tau p_{T} (GeV)");
   gr->GetYaxis()->SetTitle("efficiency");
   gr->Draw("AP");

   
   TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextSize(0.05);
   t->DrawLatex(0.25,0.92,"CMS Simulation, #sqrt{s} = 7 TeV");

   c2->Print("LM6_tauefficiency.pdf");

}
