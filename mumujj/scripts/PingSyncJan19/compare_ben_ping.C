{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //--------------------------------
  // get Ping's histogram
  //--------------------------------

  TFile* fping = TFile::Open("mass_2011_1b-v1.root");

  TH1F* hping2011a = (TH1F*) fping->Get("h2011a");
  TH1F* hping2011b = (TH1F*) fping->Get("h2011b");

  hping2011a->SetMarkerSize(1);
  hping2011a->SetMarkerStyle(20);
  hping2011b->SetMarkerSize(1);
  hping2011b->SetMarkerStyle(20);

  cout << "Ping entries 2011a " << hping2011a->GetEntries() << endl;
  cout << "Ping entries 2011b " << hping2011b->GetEntries() << endl;

  //--------------------------------
  // get Ben's histogram
  //--------------------------------

  TChain *data2011a = new TChain("t");
  TChain *data2011b = new TChain("t");

  char* path = "../../output/V00-01-04";

  data2011a->Add(Form("%s/datamay10_smallTree.root",path));
  data2011a->Add(Form("%s/dataPRv4_smallTree.root",path));
  data2011a->Add(Form("%s/dataaug05_smallTree.root",path));
  data2011a->Add(Form("%s/dataPRv6_smallTree.root",path));

  data2011b->Add(Form("%s/data2011B_smallTree.root",path));

  //--------------------------------
  // compare histograms
  //--------------------------------

  TCut njets2("njets>=2");
  TCut hlt("hlt1==1 || hlt2==1");
  TCut leadjet120("jet1.pt()>120");
  TCut mll50("dilmass>50");
  TCut zveto("dilmass<70 || dilmass>100");
  TCut pt2520("lep1.pt()>25 && lep2.pt()>20");
  TCut pt3020("lep1.pt()>30 && lep2.pt()>20");
  TCut ngoodlep2("ngoodlep>=2");
  //TCut nbtags1("nbtags2024>=1");
  TCut nbtags1("nbtags2024c>=1");

  // TCut pt3020("lep1.pt()>30 && lep2.pt()>20");
  // TCut nbtags1("nbtags33>=1");
  // TCut iso1("iso1<0.15");
  // TCut iso2("iso2<0.15");
  // TCut id1("passid1==1");
  // TCut id2("passid2==1");
  // TCut mu2("ngoodmu>=2");
  // TCut bump("mmjj>850 && mmjj<1150");
  // TCut left("mmjj>650 && mmjj<850");
  // TCut right("mmjj>1150");

  TCut sel;
  //sel += pt3020;
  sel += ngoodlep2;
  sel += pt2520;
  sel += hlt;
  sel += mll50;
  sel += zveto;
  sel += njets2;
  sel += leadjet120;
  sel += nbtags1;
  //sel += "mmjj > 600 && mmjj < 1400";
  //sel += "mmjj > 850 && mmjj < 1150";

  cout << "Using selection   " << sel.GetTitle()          << endl;
  cout << "Ben entries 2011a " << data2011a->GetEntries(sel) << endl;
  cout << "Ben entries 2011b " << data2011b->GetEntries(sel) << endl;

  TH1F* hben2011a = new TH1F("hben2011a","hben2011a",30,50,1850);
  TH1F* hben2011b = new TH1F("hben2011b","hben2011b",30,50,1850);

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  data2011a->Draw("min(mmjj,1849)>>hben2011a",sel);
  data2011b->Draw("min(mmjj,1849)>>hben2011b",sel);
  delete ctemp;

  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);

  c1->cd(1);
  gPad->SetLogy();
  gStyle->SetErrorX(0.5);
  gPad->SetGridx();
  gPad->SetGridy();
  hping2011a->SetFillColor(4);
  hping2011a->SetMarkerColor(4);
  hping2011a->Draw("E2");
  hping2011a->GetXaxis()->SetTitle("M(#mu#mujj) (GeV)");
  hben2011a->SetLineColor(2);
  hben2011a->SetMarkerColor(2);
  hben2011a->Draw("sameE1");

  TLegend *leg = new TLegend(0.7,0.8,0.99,1.0);
  leg->AddEntry(hping2011a,"Ping (2011A)","f");
  leg->AddEntry(hben2011a,"Ben (2011A)","lp");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  c1->cd(2);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetErrorX(0.5);
  hping2011b->SetFillColor(4);
  hping2011b->SetMarkerColor(4);
  hping2011b->Draw("E2");
  hping2011b->GetXaxis()->SetTitle("M(#mu#mujj) (GeV)");
  hben2011b->SetLineColor(2);
  hben2011b->SetMarkerColor(2);
  hben2011b->Draw("sameE1");

  TLegend *leg = new TLegend(0.7,0.8,0.99,1.0);
  leg->AddEntry(hping2011b,"Ping (2011B)","f");
  leg->AddEntry(hben2011b,"Ben (2011B)","lp");
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->Draw();

  TFile *f = TFile::Open("ben_mmjj.root","RECREATE");
  f->cd();
  hben2011a->Write();
  hben2011b->Write();
  f->Close();


}
