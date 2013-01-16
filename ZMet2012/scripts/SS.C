{
  gStyle->SetOptStat(111111);


  TChain* ch = new TChain("T1");

  ch->Add("../output/V00-02-01/data_53X_2012A_baby_samesign.root");
  ch->Add("../output/V00-02-01/data_53X_2012B_baby_samesign.root");
  ch->Add("../output/V00-02-01/data_53X_2012C_baby_samesign.root");
  ch->Add("../output/V00-02-01/data_53X_2012D_baby_samesign.root");

  //TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  //TCut pt2010("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut ptcuts("lep1.pt()>20.0 && lep2.pt()>20.0");
  //TCut ptcuts("lep1.pt()>20.0 && lep2.pt()>20.0");
  //TCut ptcuts("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut mumu("leptype==1");
  TCut SS("id1 * id2 > 0");
  TCut SSmu("ssmu1==1 && ssmu2==1");
  TCut njets2("njets >= 2");
  TCut bveto("nbcsvm == 0");
  TCut mjj90("mjj>70 && mjj<110");
  TCut mjj105("mjj>85 && mjj<125");
  TCut nlep3("lep3.pt()>20");
  TCut met50("pfmet>50.0");

  TCut sel       = ptcuts + mumu + SS + SSmu + nlep3 + met50;
  char* filename = "../plots/SS_baseline_temp.pdf"; 

  cout << "entries                     " << ch->GetEntries(ptcuts)              << endl;
  cout << "entries dimuons             " << ch->GetEntries(ptcuts+mumu)         << endl;
  cout << "entries SS dimuons          " << ch->GetEntries(ptcuts+mumu+SS)      << endl;
  cout << "entries SS dimuons (SS mu)  " << ch->GetEntries(ptcuts+mumu+SS+SSmu) << endl;
  cout << "entries (full selection)    " << ch->GetEntries(sel)                 << endl;


  cout << "Use selection " << sel.GetTitle() << endl;
  
  //TH1F* hmass = new TH1F("hmass","hmass",100,0,200);
  TH1F* hmass = new TH1F("hmass","",100,-10,10);
  hmass->SetLineColor(2);
  hmass->SetMarkerColor(2);

  TH2F* hmass2 = new TH2F("hmass2","hmass2",100,0,200,100,0,200);

  char* mass12 = "sqrt(  pow(lep1->mass(),2) + pow(lep2->mass(),2) + 2 * lep1->E() * lep2->E() - 2 * lep1->Px() * lep2->Px() - 2 * lep1->Py() * lep2->Py() - 2 * lep1->Pz() * lep2->Pz() )";
  char* mass13 = "sqrt(  pow(lep1->mass(),2) + pow(lep3->mass(),2) + 2 * lep1->E() * lep3->E() - 2 * lep1->Px() * lep3->Px() - 2 * lep1->Py() * lep3->Py() - 2 * lep1->Pz() * lep3->Pz() )";
  char* mass23 = "sqrt(  pow(lep2->mass(),2) + pow(lep3->mass(),2) + 2 * lep2->E() * lep3->E() - 2 * lep2->Px() * lep3->Px() - 2 * lep2->Py() * lep3->Py() - 2 * lep2->Pz() * lep3->Pz() )";

  TCanvas *ctemp = new TCanvas();
  ch->Draw(Form("dilmass-%s>>hmass",mass12),sel);

  ch->SetScanField(100);
  ch->Scan(Form("id1:id2:dilmass:%s:%s:%s:lep4.pt()",mass12,mass13,mass23),sel);
  //ch->Draw("dilmass>>hmass",sel);
  //ch->Draw(Form("dilmass:%s>>hmass2",myvar),sel);
  //ch->Draw(Form("%s:%s>>hmass2",mass13,mass23),sel);
  delete ctemp;

  TLine line;
  line.SetLineStyle(2);


  TCanvas *c1 = new TCanvas("c1","c1",1600,600);
  c1->Divide(2,1);

  c1->cd(1);
  //gPad->SetGridx(1);
  hmass->SetMinimum(0);
  hmass->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
  hmass->GetYaxis()->SetTitle("entries / 2 GeV");
  hmass->Draw("E1");
  line.DrawLine(105,0,105,hmass->GetMaximum());

  c1->cd(2);
  //gPad->SetGridx(1);
  TH1F* hmass_clone = (TH1F*) hmass->Clone("hmass_clone");
  hmass_clone->SetMinimum(0);
  hmass_clone->GetXaxis()->SetRangeUser(90,120);
  hmass_clone->Draw("E1");
  line.DrawLine(105,0,105,hmass_clone->GetMaximum());


  /*
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->cd();

  hmass2->Draw();
  */

  //c1->Print(filename);

}
