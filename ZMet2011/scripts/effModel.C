{


  TChain *chreco = new TChain("T1");
  chreco->Add("../output/V00-02-05/LM4v2_baby.root");

  TChain *chgen = new TChain("T1");
  chgen->Add("../output/V00-02-05/LM4v2_gen_baby.root");

  TH1F* h = new TH1F("h","h",1,0,1);
  h->Sumw2();

  TCut sel("dilmass>81 && dilmass<101 && njets>=2 && leptype<2");
  TCut met100("pfmet>100");
  TCut met200("pfmet>200");
  TCut met300("pfmet>300");

  TCut recoweight("trgeff * weight * 4.7");
  //TCut recoweight("trgeff");

  chreco->Draw("0.5>>h",sel*recoweight);
  cout << "No MET cut: reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;
  chreco->Draw("0.5>>h",(sel+met100)*recoweight);
  cout << "MET > 100: reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;
  chreco->Draw("0.5>>h",(sel+met200)*recoweight);
  cout << "MET > 200: reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;
  chreco->Draw("0.5>>h",(sel+met300)*recoweight);
  cout << "MET > 300: reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;

  TCut effweight("weight * 4.7");
  //TCut effweight("1");
  TCut eff0("eff0");
  TCut eff100("eff100");
  TCut eff200("eff200");
  TCut eff300("eff300");

  chgen->Draw("0.5>>h",eff0 * effweight);
  cout << "No MET cut: gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;
  chgen->Draw("0.5>>h",eff100 * effweight);
  cout << "MET > 100: gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;
  chgen->Draw("0.5>>h",eff200 * effweight);
  cout << "MET > 200: gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;
  chgen->Draw("0.5>>h",eff300 * effweight);
  cout << "MET > 300: gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;




}


/*
No MET cut: reco  136.599
No MET cut: gen   142.292
MET > 100: reco  103.422
MET > 100: gen   107.065
MET > 200: reco  58.4184
MET > 200: gen   60.3771
MET > 300: reco  25.7122
MET > 300: gen   26.2945

No MET cut: reco  136.599
No MET cut: gen   74.8801
MET > 100: reco  103.422
MET > 100: gen   54.8818
MET > 200: reco  58.4184
MET > 200: gen   29.2828
MET > 300: reco  25.7122
MET > 300: gen   12.803

*/

