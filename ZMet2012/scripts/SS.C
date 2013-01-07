{

  TChain* ch = new TChain("T1");

  ch->Add("../output/V00-02-01/data_53X_2012A_baby.root");
  ch->Add("../output/V00-02-01/data_53X_2012B_baby.root");
  ch->Add("../output/V00-02-01/data_53X_2012C_baby.root");
  ch->Add("../output/V00-02-01/data_53X_2012D_baby.root");

  //TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  //TCut pt2010("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut ptcuts("lep1.pt()>30.0 && lep2.pt()>30.0");
  //TCut ptcuts("lep1.pt()>20.0 && lep2.pt()>20.0");
  //TCut ptcuts("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut mumu("leptype==1");
  TCut SS("id1 * id2 > 0");
  TCut SSmu("ssmu1==1 && ssmu2==1");

  cout << "entries                     " << ch->GetEntries(ptcuts)              << endl;
  cout << "entries dimuons             " << ch->GetEntries(ptcuts+mumu)         << endl;
  cout << "entries SS dimuons          " << ch->GetEntries(ptcuts+mumu+SS)      << endl;
  cout << "entries SS dimuons (SS mu)  " << ch->GetEntries(ptcuts+mumu+SS+SSmu) << endl;

  TCut sel = ptcuts + mumu + SS + SSmu;

  cout << "Use selection " << sel.GetTitle() << endl;
  
  TH1F* hmass = new TH1F("hmass","",100,0,200);

  TCanvas *ctemp = new TCanvas();
  ch->Draw("min(dilmass,199.9)>>hmass",sel);
  delete ctemp;

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->cd();

  hmass->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
  hmass->GetYaxis()->SetTitle("entries / GeV");
  //hmass->Draw("hist");
  //hmass->Draw("sameE1");
  hmass->Draw("E1");


}
