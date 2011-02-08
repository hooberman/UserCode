{

  #include <iomanip>

  bool isData = true;

  char* file = "output/v7/ZJets_baby.root";
  if( isData ) file = "output/v7/lepdata_skim_baby.root";

  cout << "Adding " << file << endl;

  TFile *f = TFile::Open( file );

  TH1F* hee = new TH1F("hee","",50,0,100);
  TH1F* hmm = new TH1F("hmm","",50,0,100);
  TH1F* hem = new TH1F("hem","",50,0,100);

  TCut sel("njets>1");
  TCut ee("leptype==0&&jetptll-ptll>-5&&jetptlt-ptlt>-5&&dilmass>81&&dilmass<101");
  TCut mm("leptype==1&&nmatchedpfmuons==2&&dilmasspf>81&&dilmasspf<101");
  TCut em("leptype==2&&nmatchedpfmuons==1&&dilmass>81&&dilmass<101");
  TCut weight("weight");

  TCanvas *c1 = new TCanvas();
  c1->cd();

  TPad *plotpad = new TPad("plotpad","plotpad",0.0,0.0,1.0,0.8);
  plotpad->Draw();
  plotpad->cd();

  T1->Draw("TMath::Min(pfmet,99.9)>>hee",(sel+ee)*weight);
  T1->Draw("TMath::Min(pfmet,99.9)>>hmm",(sel+mm)*weight);
  T1->Draw("TMath::Min(pfmet,99.9)>>hem",(sel+em)*weight);

  hee->Sumw2();
  hmm->Sumw2();
  hem->Sumw2();

  hmm->Draw();
  hmm->GetXaxis()->SetTitle("pfmet (GeV)");
  hee->SetLineColor(2);
  hee->SetMarkerColor(2);
  hem->SetLineColor(4);
  hem->SetMarkerColor(4);
  hem->SetMarkerStyle(4);
  hee->Draw("samehist");
  hem->Draw("same");
  gPad->SetLogy(1);

  TLegend *leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->AddEntry(hee,"ee","l");
  leg->AddEntry(hmm,"#mu#mu");
  leg->AddEntry(hem,"e#mu");
  leg->SetFillColor(0);
  leg->Draw();

//   TLatex *t = new TLatex();
//   t->SetNDC();
//   if( isData ){
//     t->DrawLatex(0.4,0.6,  "CMS");
//     //t->DrawLatex(0.4,0.53, "Selected Z+#geq2jet Events (DATA)");
//     t->DrawLatex(0.4,0.53, "Selected Z Events (DATA)");
//     t->DrawLatex(0.4,0.46, "34.0 pb^{-1} at #sqrt{s} = 7 TeV");
//   }else{
//     t->DrawLatex(0.4,0.6, "CMS");
//     t->DrawLatex(0.4,0.53,"Selected Z+0jet Events (Z+jets MC)");
//   }

  c1->cd();

  TPad *pullpad = new TPad("pullpad","pullpad",0.0,0.8,1.0,1.0);
  pullpad->Draw();
  pullpad->cd();

  TH1F* hratio = (TH1F*) hmm->Clone();
  hratio->Divide(hee);
  hratio->Draw();
  gPad->SetGridy();
  hratio->SetMinimum(0.8);
  hratio->SetMaximum(1.6);
  hratio->GetYaxis()->SetLabelSize(0.15);
  hratio->GetYaxis()->SetNdivisions(7);
  hratio->GetYaxis()->SetTitle("#mu#mu / ee  ");
  hratio->GetYaxis()->SetTitleSize(0.2);
  hratio->GetYaxis()->SetTitleOffset(0.2);


  int bin30 = hee->FindBin(30);
  int bin60 = hee->FindBin(60);
  int bin120 = hee->FindBin(120);
  
  cout << "ee tot         " << hee->Integral() << endl;
  cout << "ee met>30  GeV " << hee->Integral(bin30,1000) << endl;
  cout << "ee met>60  GeV " << hee->Integral(bin60,1000) << endl;
  cout << "ee met>120 GeV " << hee->Integral(bin120,1000) << endl;

  cout << "mm tot         " << hmm->Integral() << endl;
  cout << "mm met>30  GeV " << hmm->Integral(bin30,1000) << endl;
  cout << "mm met>60  GeV " << hmm->Integral(bin60,1000) << endl;
  cout << "mm met>120 GeV " << hmm->Integral(bin120,1000) << endl;

  int width1 = 15;
  int width2 =  5;

  cout << "|" << setw(width1) << ""               << setw(width2)
       << "|" << setw(width1) << "N(met>30 GeV)"  << setw(width2)
       << "|" << setw(width1) << "N(met>60 GeV)"  << setw(width2)
       << "|" << setw(width1) << "N(met>120 GeV)" << setw(width2) << "|" << endl;

  cout << "|" << setw(width1) << "ee"                       << setw(width2)
       << "|" << setw(width1) << hee->Integral(bin30,1000)  << setw(width2)
       << "|" << setw(width1) << hee->Integral(bin60,1000)  << setw(width2)
       << "|" << setw(width1) << hee->Integral(bin120,1000) << setw(width2) << "|" << endl;

  cout << "|" << setw(width1) << "mm"                       << setw(width2)
       << "|" << setw(width1) << hmm->Integral(bin30,1000)  << setw(width2)
       << "|" << setw(width1) << hmm->Integral(bin60,1000)  << setw(width2)
       << "|" << setw(width1) << hmm->Integral(bin120,1000) << setw(width2) << "|" << endl;

  cout << "|" << setw(width1) << "em"                       << setw(width2)
       << "|" << setw(width1) << hem->Integral(bin30,1000)  << setw(width2)
       << "|" << setw(width1) << hem->Integral(bin60,1000)  << setw(width2)
       << "|" << setw(width1) << hem->Integral(bin120,1000) << setw(width2) << "|" << endl;
}



