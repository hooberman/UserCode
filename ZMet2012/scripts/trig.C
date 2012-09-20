{

  TChain *ch = new TChain("T1");
  //ch->Add("../output/V00-00-24/dataskim2010_all_baby_2jets.root");
  ch->Add("../output/V00-01-06/data_53X_baby_2jets.root");

  //TCut presel("pfmet>100 && njets40>=3 && lep1.pt()>20.0 && lep2.pt()>20.0 && dilmass>15.0 && dilmass<70.0");
  TCut presel("pfmet>150 && njets40>=2 && ht40>=100.0 && dilmass>15.0 && dilmass<70.0");

  //TCut leptons("");

  TCut eetrig("leptype==0 && ee==1");  
  TCut alleetrig("leptype==0 && (ee==1||el27==1)");

  TCut mmtrig("leptype==1 && mm==1");  
  TCut allmmtrig("leptype==1 && (mm==1||mu==1||mu21==1||mmtk==1)");

  TCut emtrig("leptype==2 && (em==1||me==1)");  
  TCut allemtrig("leptype==2 && (em==1||me==1||mu==1||mu21==1||el27==1)");

  TH1F* hee    = new TH1F("hee"   ,"",15,0,300);
  TH1F* heeall = new TH1F("heeall","",15,0,300);

  TH1F* hmm    = new TH1F("hmm"   ,"",15,0,300);
  TH1F* hmmall = new TH1F("hmmall","",15,0,300);

  TH1F* hem    = new TH1F("hem"   ,"",15,0,300);
  TH1F* hemall = new TH1F("hemall","",15,0,300);

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  ch->Draw("min(dilmass,299)>>hee"   ,presel+eetrig);
  ch->Draw("min(dilmass,299)>>heeall",presel+alleetrig);

  ch->Draw("min(dilmass,299)>>hmm"   ,presel+mmtrig);
  ch->Draw("min(dilmass,299)>>hmmall",presel+allmmtrig);

  ch->Draw("min(dilmass,299)>>hem"   ,presel+emtrig);
  ch->Draw("min(dilmass,299)>>hemall",presel+allemtrig);
  delete ctemp;

  TCanvas *c1 = new TCanvas("c1","",1500,500);
  c1->Divide(3,1);

  c1->cd(1);
  heeall->GetXaxis()->SetTitle("M(ee) [GeV]");
  heeall->SetLineColor(2);
  heeall->SetMarkerColor(2);
  heeall->Draw("E1");
  hee->Draw("samehist");
  heeall->Draw("sameE1");

  TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
  leg->AddEntry(hee,"dilepton","l");
  leg->AddEntry(heeall,"dilepton OR single lepton");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  //leg->Draw();

  c1->cd(2);
  hmmall->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
  hmmall->SetLineColor(2);
  hmmall->SetMarkerColor(2);
  hmmall->Draw("E1");
  hmm->Draw("samehist");
  hmmall->Draw("sameE1");

  c1->cd(3);
  hemall->GetXaxis()->SetTitle("M(e#mu) [GeV]");
  hemall->SetLineColor(2);
  hemall->SetMarkerColor(2);
  hemall->Draw("E1");
  hem->Draw("samehist");
  hemall->Draw("sameE1");
  
  cout << endl;
  cout << "ee       " << hee->GetEntries()    << endl;
  cout << "ee all   " << heeall->GetEntries() << endl;
  cout << "increase " << heeall->GetEntries() / (float) hee->GetEntries() << endl;
  cout << endl;
  cout << "mm       " << hmm->GetEntries()    << endl;
  cout << "mm all   " << hmmall->GetEntries() << endl;
  cout << "increase " << hmmall->GetEntries() / (float) hmm->GetEntries() << endl;
  cout << endl;
  cout << "em       " << hem->GetEntries()    << endl;
  cout << "em all   " << hemall->GetEntries() << endl;
  cout << "increase " << hemall->GetEntries() / (float) hem->GetEntries() << endl;

  // cout << "mm     " << hmm->GetEntries()    << endl;
  // cout << "mm all " << hmmall->GetEntries() << endl;
  // cout << "em     " << hem->GetEntries()    << endl;
  // cout << "em all " << hemall->GetEntries() << endl;

  c1->Print("../plots/trig.pdf");





}
