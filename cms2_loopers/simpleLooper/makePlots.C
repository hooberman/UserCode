{

  TChain *ch = new TChain("Events");
  ch->Add("output/data_default_baby.root");

  TH1F* hmm  = new TH1F("hmm" ,"hmm",50,0,50);
  TH1F* hmmp = new TH1F("hmmp","hmmp",50,0,50);

  ch->Draw("ht>>hmm","mm==1");
  ch->Draw("ht>>hmmp","mm==1&&mmht==1");

  TGraphAsymmErrors *grmm = new TGraphAsymmErrors();
  grmm->BayesDivide(hmm,hmmp);

  grmm->Draw("AP");

}
