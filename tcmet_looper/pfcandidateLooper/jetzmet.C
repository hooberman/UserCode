{

  bool printgif = true;

  gStyle->SetOptStat(0);

  char* path = "output/PVT/promptreco/dcsonly";

  TChain *dy = new TChain("Events");
  dy->Add(Form("%s/dyee_default_baby.root",path));
  dy->Add(Form("%s/dymm_default_baby.root",path));

  TChain *chdata = new TChain("Events");
  chdata->Add(Form("%s/data_default_baby.root",path));

  TChain *h130 = new TChain("Events");
  h130->Add(Form("%s/h130_default_baby.root",path));
  
  TCut zmass("dilep.mass()>76&&dilep.mass()<106");
  TCut jet0("njets30==0");
  TCut jet1("njets30==1");
  TCut jet2("njets30>1");
  TCut ee("leptype==0");
  TCut mm("leptype==3");
  TCut em("leptype==2||leptype==1");
  TCut dphizjet("acos(cos(dilep.phi()-jet.phi()))>2");
  TCut lepveto("nlep==0");

  //TCut sel1          = zmass+(mm||ee)+"njets30==1 && jetpv==1";
  TCut sel1          = (zmass+(mm||ee)+lepveto+dphizjet+"njets30==1 && jetpv==1")*"weight*davtxweight";
  TCut sel1_em       = zmass+em+"njets30==1 && jetpv==1";
  TCut higgssel1     = (mm||ee)+"njets30==1 && jetpv==1";

  TCanvas *ctemp = new TCanvas();
  chdata->Draw("TMath::Min(jetzmet,99.9)>>hj(50,0,100)" , sel1);
  chdata->Draw("TMath::Min(pfmet,99.9)>>hpf(50,0,100)"  , sel1);
  chdata->Draw("TMath::Min(TMath::Min(pfmet,jetzmet),99.9)>>hmin(50,0,100)"  , sel1);
  dy->Draw("TMath::Min(pfmet,99.9)>>hpf_dy(50,0,100)"  , sel1);
  dy->Draw("TMath::Min(TMath::Min(pfmet,jetzmet),99.9)>>hmin_dy(50,0,100)"  , sel1);
  h130->Draw("TMath::Min(pfmet,99.9)>>hpf_higgs(50,0,100)"  , higgssel1);
  h130->Draw("TMath::Min(TMath::Min(pfmet,jetzmet),99.9)>>hmin_higgs(50,0,100)"  , higgssel1);
  chdata->Draw("TMath::Min(pfmet,99.9):TMath::Min(jetzmet,99.9)>>hscatter(50,0,100,50,0,100)" , sel1);
  chdata->Draw("TMath::Min(pfmet,99.9):TMath::Min(jetzmet,99.9)>>hscatter_em(50,0,100,50,0,100)" , sel1_em);
  h130->Draw("TMath::Min(pfmet,99.9):TMath::Min(jetzmet,99.9)>>hscatter_higgs(50,0,100,50,0,100)" , higgssel1);
  delete ctemp;

  TLine line;
  line.SetLineColor(2);

  //----------------------------------
  // pfmet and jet-Z MET in data
  //----------------------------------

  /*
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);

  c1->cd(1);
  gPad->SetLogy();

  hpf->Draw();
  hpf->GetXaxis()->SetTitle("MET (GeV)");
  hj->SetMarkerColor(2);
  hj->SetLineColor(2);
  hj->Draw("sameE1");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(hpf,"pfmet","l");
  leg->AddEntry(hj,"jet-Z MET","p");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  c1->cd(2);

  hscatter->Draw();
  hscatter->GetXaxis()->SetTitle("jet-Z MET (GeV)");
  hscatter->GetYaxis()->SetTitle("pfmet (GeV)");
  line.DrawLine(0,45,100,45);
  line.DrawLine(45,0,45,100);


  if( printgif ) c1->Print("plots/met_data.gif");
  */
  
  
  //-------------------------------------
  // pfmet and min-MET in data and DY MC
  //-------------------------------------
  
  TCanvas *c2 = new TCanvas("c2","",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  gPad->SetLogy();
  
  hpf->Draw();

  hpf->Draw();
  hpf->GetXaxis()->SetTitle("MET (GeV)");
  hmin->SetMarkerColor(2);
  hmin->SetLineColor(2);
  hmin->Draw("sameE1");

  TLegend *leg2 = new TLegend(0.7,0.7,0.9,0.9);
  leg2->AddEntry(hpf,"pfmet","l");
  leg2->AddEntry(hmin,"min-MET","p");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(1);
  leg2->Draw();

  c2->cd(2);
  gPad->SetLogy();
  
  hpf_dy->Draw();

  hpf_dy->Draw();
  hpf_dy->GetXaxis()->SetTitle("MET (GeV)");
  hmin_dy->SetMarkerColor(2);
  hmin_dy->SetLineColor(2);
  hmin_dy->Draw("sameE1");
  leg2->Draw();
  
  if( printgif ) c2->Print("plots/minmet.gif");
    
  //--------------------------------------------------
  // pfmet vs. min-MET, scatter plot for Higgs MC
  //--------------------------------------------------
  /*
  TCanvas *c3 = new TCanvas("c3","",1200,600);
  c3->Divide(2,1);

  c3->cd(1);
  gPad->SetLogy();
  
  hpf_higgs->Draw();

  hpf_higgs->Draw();
  hpf_higgs->GetXaxis()->SetTitle("MET (GeV)");
  hmin_higgs->SetMarkerColor(2);
  hmin_higgs->SetLineColor(2);
  hmin_higgs->Draw("samesE1");
  leg2->Draw();

  c3->cd(2);
  //gStyle->SetOptStat(0);
  hscatter_higgs->Draw();
  hscatter_higgs->GetXaxis()->SetTitle("jet-Z MET (GeV)");
  hscatter_higgs->GetYaxis()->SetTitle("pfmet (GeV)");

  if( printgif ) c3->Print("plots/higgs_scatter.gif");
  */


  /*
  //-------------------------------------------------
  // pfmet vs. jet-Z MET in data, SF vs. DF
  //-------------------------------------------------

  TCanvas *c4 = new TCanvas("c5","",1200,600);
  c4->Divide(2,1);

 
  c4->cd(2);
  hscatter->Draw();
  hscatter->GetXaxis()->SetTitle("jet-Z MET (GeV)");
  hscatter->GetYaxis()->SetTitle("pfmet (GeV)");
  line.DrawLine(0,45,100,45);
  line.DrawLine(45,0,45,100);

  c4->cd(1);
  hscatter_em->Draw();
  hscatter_em->GetXaxis()->SetTitle("jet-Z MET (GeV)");
  hscatter_em->GetYaxis()->SetTitle("pfmet (GeV)");
  line.DrawLine(0,45,100,45);
  line.DrawLine(45,0,45,100);

  if( printgif ) c4->Print("plots/scatter.gif");
  */
}
