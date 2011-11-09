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

  TCut pass("pass&&!passz");
  //TCut sig("pfmet>275&&ht>300");   TCut sigup("pfmetUp>275&&htUp>300");  TCut sigdn("pfmetDown>275&&htDown>300");
  TCut sig("pfmet>275&&ht>600");   TCut sigup("pfmetUp>275&&htUp>600");  TCut sigdn("pfmetDown>275&&htDown>600");


  TCut weight("weight * ndavtxweight * trgeff");
  TCut weightup("weight * ndavtxweight * trgeff * ksusyup/ksusy");
  TCut weightdn("weight * ndavtxweight * trgeff * ksusydn/ksusy");

  ch->Draw("m12:m0>>msugra"     , (pass+sig)   * weight   );
  ch->Draw("m12:m0>>msugra_jup" , (pass+sigup) * weight   );
  ch->Draw("m12:m0>>msugra_jdn" , (pass+sigdn) * weight   );
  ch->Draw("m12:m0>>msugra_kup" , (pass+sig)   * weightup );
  ch->Draw("m12:m0>>msugra_kdn" , (pass+sig)   * weightdn );

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();

  gPad->SetRightMargin(0.2);
  msugra->Draw("colz");

  int LM6bin = msugra->FindBin(80,400);

  cout << "LM6 yields bin " << LM6bin << endl;
  cout << "nominal    " << msugra->GetBinContent(LM6bin) << endl;
  cout << "JES up     " << msugra_jup->GetBinContent(LM6bin) << endl;
  cout << "JES dn     " << msugra_jdn->GetBinContent(LM6bin) << endl;
  cout << "k up       " << msugra_kup->GetBinContent(LM6bin) << endl;
  cout << "k dn       " << msugra_kdn->GetBinContent(LM6bin) << endl;

  TFile *file = TFile::Open("histos.root","RECREATE");
  file->cd();
  msugra->Write();
  msugra_jup->Write();
  msugra_jdn->Write();
  msugra_kup->Write();
  msugra_kdn->Write();
  file->Close();
  
}
