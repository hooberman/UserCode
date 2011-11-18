{


  //----------------------------
  // load TChain
  //----------------------------

  TChain *ch = new TChain("t");
  ch->Add("../../output/V00-02-06/highpt/LMscan_smallTree.root");

  //----------------------------
  // make histograms
  //----------------------------

  const int   nm0points    = 100;
  const float m0min        = 20.;
  const float m0max        = 2020.;
  const int   nm12points   = 38;
  const float m12min       = 20.;
  const float m12max       = 780.;
  
  TH2F* msugra        = new TH2F("msugra"      ,"msugra"      ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_ptll   = new TH2F("msugra_ptll" ,"msugra_ptll" ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_jup    = new TH2F("msugra_jup"  ,"msugra_jup"  ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_jdn    = new TH2F("msugra_jdn"  ,"msugra_jdn"  ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_kup    = new TH2F("msugra_kup"  ,"msugra_kup"  ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);
  TH2F* msugra_kdn    = new TH2F("msugra_kdn"  ,"msugra_kdn"  ,nm0points,m0min-10,m0max-10,nm12points,m12min-10,m12max-10);

  msugra->Sumw2();
  msugra_jup->Sumw2();
  msugra_jdn->Sumw2();
  msugra_kup->Sumw2();
  msugra_kdn->Sumw2();
  msugra_ptll->Sumw2();

  //----------------------------
  // selection
  //----------------------------

  TCut pass("pass&&!passz");
  float KKC = 1;

  //--- high MET
  //TCut sig("pfmet>275&&ht>300");   TCut sigup("pfmetUp>275&&htUp>300");  TCut sigdn("pfmetDown>275&&htDown>300"); KKC = 2.475;
  
  //--- tight SR2
  //TCut sig("pfmet>275&&ht>600");   TCut sigup("pfmetUp>275&&htUp>600");  TCut sigdn("pfmetDown>275&&htDown>600"); KKC = 2.464;

  //--- high HT
  TCut sig("pfmet>200&&ht>600");   TCut sigup("pfmetUp>200&&htUp>600");  TCut sigdn("pfmetDown>200&&htDown>600"); KKC = 2.156;

  TCut ptll("dilpt>200&&ht>600 && ((leptype<2&&pfmet>75) || (leptype==2&&pfmet>50))"); 
  TCut ee("leptype==0");
  TCut mm("leptype==1");
  TCut em("leptype==2");

  //----------------------------
  // Z peak scale factor
  //----------------------------

  float nee   = ch->GetEntries(pass+sig+ee);
  float nmm   = ch->GetEntries(pass+sig+mm);
  float nem   = ch->GetEntries(pass+sig+em);
  float ntot  = nee + nmm + nem;
  float fee   = nee / ntot;
  float fmm   = nmm / ntot;
  float fem   = nem / ntot;
  float scale = fee * 1.05 + fmm * 1.12 + fem * 1.08;

  //----------------------------
  // Z peak scale factor (ptll)
  //----------------------------

  float neep   = ch->GetEntries(pass+ptll+ee);
  float nmmp   = ch->GetEntries(pass+ptll+mm);
  float nemp   = ch->GetEntries(pass+ptll+em);
  float ntotp  = neep + nmmp + nemp;
  float feep   = neep / ntotp;
  float fmmp   = nmmp / ntotp;
  float femp   = nemp / ntotp;
  float scalep = feep * 1.05 + fmmp * 1.12 + femp * 1.08;

  //----------------------------
  // weights
  //----------------------------

  TCut  scaleweight(Form("%.3f",scale));
  TCut  scalepweight(Form("%.3f",scalep));
  TCut weight  ("weight * ndavtxweight * trgeff");
  TCut weightup("weight * ndavtxweight * trgeff * ksusyup/ksusy");
  TCut weightdn("weight * ndavtxweight * trgeff * ksusydn/ksusy");

  //----------------------------
  // draw histos
  //----------------------------

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  ch->Draw("m12:m0>>msugra"     , (pass+sig)   * weight   * scaleweight  );
  ch->Draw("m12:m0>>msugra_jup" , (pass+sigup) * weight   * scaleweight  );
  ch->Draw("m12:m0>>msugra_jdn" , (pass+sigdn) * weight   * scaleweight  );
  ch->Draw("m12:m0>>msugra_kup" , (pass+sig)   * weightup * scaleweight  );
  ch->Draw("m12:m0>>msugra_kdn" , (pass+sig)   * weightdn * scaleweight  );
  ch->Draw("m12:m0>>msugra_ptll", (pass+ptll)  * weight   * scalepweight );
  delete ctemp;

  cout << "Scaling ptll hist by " << KKC << endl;
  msugra_ptll->Scale(KKC);

  //----------------------------
  // display yield
  //----------------------------

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();

  gPad->SetRightMargin(0.2);
  msugra->Draw("colz");

  //-------------------------------
  // show yields for sample point
  //-------------------------------

  //int LMbin = msugra->FindBin(60,260); //LM1
  int LMbin = msugra->FindBin(80,400); //LM6
  //int LMbin = msugra->FindBin(80,400);

  cout << "LM yields bin " << LMbin << endl;
  cout << "nominal    " << Form("%.1f +/- %.1f",msugra->GetBinContent(LMbin),msugra->GetBinError(LMbin)) << endl;
  cout << "JES up     " << msugra_jup->GetBinContent(LMbin) << endl;
  cout << "JES dn     " << msugra_jdn->GetBinContent(LMbin) << endl;
  cout << "k up       " << msugra_kup->GetBinContent(LMbin) << endl;
  cout << "k dn       " << msugra_kdn->GetBinContent(LMbin) << endl;

  //-------------------------------
  // write out histograms
  //-------------------------------

  TFile *file = TFile::Open("histos.root","RECREATE");
  file->cd();
  msugra->Write();
  msugra_jup->Write();
  msugra_jdn->Write();
  msugra_kup->Write();
  msugra_kdn->Write();
  msugra_ptll->Write();
  file->Close();
  
}
