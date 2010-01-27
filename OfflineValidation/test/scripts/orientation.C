{
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);

int dx=900;
int dy=600;
TCanvas *c1=new TCanvas("c1","",dx,dy);
c1->Divide(3,2);
TCanvas *c2=new TCanvas("c2","",dx,dy);
c2->Divide(3,2);
TCanvas *c3=new TCanvas("c3","",dx,dy);
c3->Divide(3,2);

TFile f("rootfiles/peak.root");

f.cd();
TrackerOfflineValidation.cd();
MyStrip.cd();
TOBBarrel_4.cd();
TOBHalfBarrel_1.cd();

c1->cd(1);
h_uOrientation_TOBLayer_0->Draw();
h_uOrientation_TOBLayer_0->SetTitle("uOrientation TOB Layer 1");
c1->cd(2);
h_uOrientation_TOBLayer_1->Draw();
h_uOrientation_TOBLayer_1->SetTitle("uOrientation TOB Layer 2");
c1->cd(3);
h_uOrientation_TOBLayer_2->Draw();
h_uOrientation_TOBLayer_2->SetTitle("uOrientation TOB Layer 3");
c1->cd(4);
h_uOrientation_TOBLayer_3->Draw();
h_uOrientation_TOBLayer_3->SetTitle("uOrientation TOB Layer 4");
c1->cd(5);
h_uOrientation_TOBLayer_4->Draw();
h_uOrientation_TOBLayer_4->SetTitle("uOrientation TOB Layer 5");
c1->cd(6);
h_uOrientation_TOBLayer_5->Draw();
h_uOrientation_TOBLayer_5->SetTitle("uOrientation TOB Layer 6");

c2->cd(1);
h_vOrientation_TOBLayer_0->Draw();
h_vOrientation_TOBLayer_0->SetTitle("vOrientation TOB Layer 1");
c2->cd(2);
h_vOrientation_TOBLayer_1->Draw();
h_vOrientation_TOBLayer_1->SetTitle("vOrientation TOB Layer 2");
c2->cd(3);
h_vOrientation_TOBLayer_2->Draw();
h_vOrientation_TOBLayer_2->SetTitle("vOrientation TOB Layer 3");
c2->cd(4);
h_vOrientation_TOBLayer_3->Draw();
h_vOrientation_TOBLayer_3->SetTitle("vOrientation TOB Layer 4");
c2->cd(5);
h_vOrientation_TOBLayer_4->Draw();
h_vOrientation_TOBLayer_4->SetTitle("vOrientation TOB Layer 5");
c2->cd(6);
h_vOrientation_TOBLayer_5->Draw();
h_vOrientation_TOBLayer_5->SetTitle("vOrientation TOB Layer 6");

c3->cd(1);
h_wOrientation_TOBLayer_0->Draw();
h_wOrientation_TOBLayer_0->SetTitle("wOrientation TOB Layer 1");
c3->cd(2);
h_wOrientation_TOBLayer_1->Draw();
h_wOrientation_TOBLayer_1->SetTitle("wOrientation TOB Layer 2");
c3->cd(3);
h_wOrientation_TOBLayer_2->Draw();
h_wOrientation_TOBLayer_2->SetTitle("wOrientation TOB Layer 3");
c3->cd(4);
h_wOrientation_TOBLayer_3->Draw();
h_wOrientation_TOBLayer_3->SetTitle("wOrientation TOB Layer 4");
c3->cd(5);
h_wOrientation_TOBLayer_4->Draw();
h_wOrientation_TOBLayer_4->SetTitle("wOrientation TOB Layer 5");
c3->cd(6);
h_wOrientation_TOBLayer_5->Draw();
h_wOrientation_TOBLayer_5->SetTitle("wOrientation TOB Layer 6");

f.cd();
TrackerOfflineValidation.cd();
MyStrip.cd();
TOBBarrel_4.cd();
TOBHalfBarrel_2.cd();

c1->cd(1);
h_uOrientation_TOBLayer_0->Draw("same");
h_uOrientation_TOBLayer_0->SetLineColor(2);
c1->cd(2);
h_uOrientation_TOBLayer_1->Draw("same");
h_uOrientation_TOBLayer_1->SetLineColor(2);
c1->cd(3);
h_uOrientation_TOBLayer_2->Draw("same");
h_uOrientation_TOBLayer_2->SetLineColor(2);
c1->cd(4);
h_uOrientation_TOBLayer_3->Draw("same");
h_uOrientation_TOBLayer_3->SetLineColor(2);
c1->cd(5);
h_uOrientation_TOBLayer_4->Draw("same");
h_uOrientation_TOBLayer_4->SetLineColor(2);
c1->cd(6);
h_uOrientation_TOBLayer_5->Draw("same");
h_uOrientation_TOBLayer_5->SetLineColor(2);

c2->cd(1);
h_vOrientation_TOBLayer_0->Draw("same");
h_vOrientation_TOBLayer_0->SetLineColor(2);
c2->cd(2);
h_vOrientation_TOBLayer_1->Draw("same");
h_vOrientation_TOBLayer_1->SetLineColor(2);
c2->cd(3);
h_vOrientation_TOBLayer_2->Draw("same");
h_vOrientation_TOBLayer_2->SetLineColor(2);
c2->cd(4);
h_vOrientation_TOBLayer_3->Draw("same");
h_vOrientation_TOBLayer_3->SetLineColor(2);
c2->cd(5);
h_vOrientation_TOBLayer_4->Draw("same");
h_vOrientation_TOBLayer_4->SetLineColor(2);
c2->cd(6);
h_vOrientation_TOBLayer_5->Draw("same");
h_vOrientation_TOBLayer_5->SetLineColor(2);


c3->cd(1);
h_wOrientation_TOBLayer_0->Draw("same");
h_wOrientation_TOBLayer_0->SetLineColor(2);
c3->cd(2);
h_wOrientation_TOBLayer_1->Draw("same");
h_wOrientation_TOBLayer_1->SetLineColor(2);
c3->cd(3);
h_wOrientation_TOBLayer_2->Draw("same");
h_wOrientation_TOBLayer_2->SetLineColor(2);
c3->cd(4);
h_wOrientation_TOBLayer_3->Draw("same");
h_wOrientation_TOBLayer_3->SetLineColor(2);
c3->cd(5);
h_wOrientation_TOBLayer_4->Draw("same");
h_wOrientation_TOBLayer_4->SetLineColor(2);
c3->cd(6);
h_wOrientation_TOBLayer_5->Draw("same");
h_wOrientation_TOBLayer_5->SetLineColor(2);



}
