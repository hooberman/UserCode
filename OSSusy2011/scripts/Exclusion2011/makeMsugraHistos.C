{

  TChain *ch = new TChain("t");
  ch->Add("../../output/V00-02-06/highpt/LMscan_smallTree.root");

  const int   nm0points    = 100;
  const float m0min        = 20.;
  const float m0max        = 2020.;
  const int   nm12points   = 38;
  const float m12min       = 20.;
  const float m12max       = 780.;
  
  TH2F* msugra        = new TH2F("msugra"    ,"msugra"    ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_jup    = new TH2F("msugra_jup","msugra_jup",nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_jdn    = new TH2F("msugra_jdn","msugra_jdn",nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_kup    = new TH2F("msugra_kup","msugra_kup",nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_kdn    = new TH2F("msugra_kdn","msugra_kdn",nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  msugra->Sumw2();
  msugra_jup->Sumw2();
  msugra_jdn->Sumw2();
  msugra_kup->Sumw2();
  msugra_kdn->Sumw2();

  TCut pass("pass&&!passz");
  //TCut sig("pfmet>275&&ht>300");   TCut sigup("pfmetUp>275&&htUp>300");  TCut sigdn("pfmetDown>275&&htDown>300");
  TCut sig("pfmet>275&&ht>600");   TCut sigup("pfmetUp>275&&htUp>600");  TCut sigdn("pfmetDown>275&&htDown>600");
  TCut ee("leptype==0");
  TCut mm("leptype==1");
  TCut em("leptype==2");

  float nee = ch->GetEntries(pass+sig+ee);
  float nmm = ch->GetEntries(pass+sig+mm);
  float nem = ch->GetEntries(pass+sig+em);
  float ntot = nee + nmm + nem;
  float fee = nee / ntot;
  float fmm = nmm / ntot;
  float fem = nem / ntot;
  float scale = fee * 1.05 + fmm * 1.12 + fem * 1.08;

  TCut scaleweight(Form("%.3f",scale));
  TCut weight  ("weight * ndavtxweight * trgeff");
  TCut weightup("weight * ndavtxweight * trgeff * ksusyup/ksusy");
  TCut weightdn("weight * ndavtxweight * trgeff * ksusydn/ksusy");

  ch->Draw("m12:m0>>msugra"     , (pass+sig)   * weight   * scaleweight );
  ch->Draw("m12:m0>>msugra_jup" , (pass+sigup) * weight   * scaleweight );
  ch->Draw("m12:m0>>msugra_jdn" , (pass+sigdn) * weight   * scaleweight );
  ch->Draw("m12:m0>>msugra_kup" , (pass+sig)   * weightup * scaleweight );
  ch->Draw("m12:m0>>msugra_kdn" , (pass+sig)   * weightdn * scaleweight );

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();

  gPad->SetRightMargin(0.2);
  msugra->Draw("colz");

  //int LMbin = msugra->FindBin(60,260); //LM6
  int LMbin = msugra->FindBin(80,400); //LM6
  //int LMbin = msugra->FindBin(80,400);

  cout << "LM yields bin " << LMbin << endl;
  cout << "nominal    " << Form("%.1f +/- %.1f",msugra->GetBinContent(LMbin),msugra->GetBinError(LMbin)) << endl;
  cout << "JES up     " << msugra_jup->GetBinContent(LMbin) << endl;
  cout << "JES dn     " << msugra_jdn->GetBinContent(LMbin) << endl;
  cout << "k up       " << msugra_kup->GetBinContent(LMbin) << endl;
  cout << "k dn       " << msugra_kdn->GetBinContent(LMbin) << endl;

  TFile *file = TFile::Open("histos.root","RECREATE");
  file->cd();
  msugra->Write();
  msugra_jup->Write();
  msugra_jdn->Write();
  msugra_kup->Write();
  msugra_kdn->Write();
  file->Close();
  
}
