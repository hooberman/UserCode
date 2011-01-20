{

  TChain *chdata = new TChain("T1");
  chdata->Add("output/oct15th_v3/lepdata_baby.root");

  TChain *chmc = new TChain("T1");
  chmc->Add("output/oct15th_v3/ZJets_baby.root");
  
  chdata->Draw("dilmass>>hbb(40,70,110)",  "leptype==0&&ecaltype==1");
  chmc  ->Draw("dilmass>>hbbmc(40,70,110)","(leptype==0&&ecaltype==1)*weight");
  
  chdata->Draw("dilmasscor>>hbbcor(40,70,110)",  "leptype==0&&ecaltype==1");
  chmc  ->Draw("dilmasscor>>hbbmccor(40,70,110)","(leptype==0&&ecaltype==1)*weight");

  chdata->Draw("dilmass>>hee(40,70,110)",  "leptype==0&&ecaltype==2");
  chmc  ->Draw("dilmass>>heemc(40,70,110)","(leptype==0&&ecaltype==2)*weight");
  
  chdata->Draw("dilmasscor>>heecor(40,70,110)",  "leptype==0&&ecaltype==2");
  chmc  ->Draw("dilmasscor>>heemccor(40,70,110)","(leptype==0&&ecaltype==2)*weight");

  chdata->Draw("dilmass>>heb(40,70,110)",  "leptype==0&&ecaltype==3");
  chmc  ->Draw("dilmass>>hebmc(40,70,110)","(leptype==0&&ecaltype==3)*weight");
  
  chdata->Draw("dilmasscor>>hebcor(40,70,110)",  "leptype==0&&ecaltype==3");
  chmc  ->Draw("dilmasscor>>hebmccor(40,70,110)","(leptype==0&&ecaltype==3)*weight");

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(3,2);

  c1->cd(1);
  hbbmc->Rebin(2);
  hbb->Rebin(2);
  hbbmc->SetFillColor(5);
  hbbmc->Draw("hist");
  hbb->Draw("sameE1");
  hbbmc->GetXaxis()->SetTitle("dilepton mass (GeV)");
  hbbmc->SetTitle("EB-EB (un-corrected)");

  c1->cd(4);
  hbbmccor->Rebin(2);
  hbbcor->Rebin(2);
  hbbmccor->SetFillColor(5);
  hbbmccor->Draw("hist");
  hbbcor->Draw("sameE1");
  hbbmccor->GetXaxis()->SetTitle("dilepton mass (GeV)");
  hbbmccor->SetTitle("EB-EB (corrected)");

  c1->cd(3);
  heemc->Rebin(2);
  hee->Rebin(2);
  heemc->SetFillColor(5);
  heemc->Draw("hist");
  hee->Draw("sameE1");
  heemc->GetXaxis()->SetTitle("dilepton mass (GeV)");
  heemc->SetTitle("EE-EE (un-corrected)");

  c1->cd(6);
  heemccor->Rebin(2);
  heecor->Rebin(2);
  heemccor->SetFillColor(5);
  heemccor->Draw("hist");
  heecor->Draw("sameE1");
  heemccor->GetXaxis()->SetTitle("dilepton mass (GeV)");
  heemccor->SetTitle("EE-EE (corrected)");

  c1->cd(2);
  hebmc->Rebin(2);
  heb->Rebin(2);
  hebmc->SetFillColor(5);
  hebmc->Draw("hist");
  heb->Draw("sameE1");
  hebmc->GetXaxis()->SetTitle("dilepton mass (GeV)");
  hebmc->SetTitle("EE-EB (un-corrected)");

  c1->cd(5);
  hebmccor->Rebin(2);
  hebcor->Rebin(2);
  hebmccor->SetFillColor(5);
  hebmccor->Draw("hist");
  hebcor->Draw("sameE1");
  hebmccor->GetXaxis()->SetTitle("dilepton mass (GeV)");
  hebmccor->SetTitle("EE-EB (corrected)");
}
