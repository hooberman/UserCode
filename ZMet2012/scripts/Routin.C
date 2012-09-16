{

  TChain *ch = new TChain("T1");
  ch->Add("../output/V00-00-21/dataskim2010_baby_2jets.root");

  TCut sel("dilmass>12 && pfmet<60 && njets40>=2");
  TCut peak("dilmass>81 && dilmass<101");
  TCut lo("dilmass>12 && dilmass<70");
  TCut ee("leptype==0 && ee==1");
  TCut mm("leptype==1 && mm==1");
  TCut em("leptype==2 && (em==1 || me==1)");

  TH1F* hee = new TH1F("hee","",100,1,201);
  TH1F* hmm = new TH1F("hmm","",100,1,201);
  TH1F* hsf = new TH1F("hsf","",100,1,201);
  TH1F* hem = new TH1F("hem","",100,1,201);

  TCanvas *ctemp = new TCanvas();
  ctemp->cd();
  ch->Draw("dilmass>>hee",sel+ee);
  ch->Draw("dilmass>>hmm",sel+mm);
  ch->Draw("dilmass>>hsf",sel+(ee||mm));
  ch->Draw("dilmass>>hem",sel+em);
  delete ctemp;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  c1->cd(1);
  gPad->SetLogy();
  hsf->SetMinimum(1);
  hsf->SetMaximum(1e5);
  hsf->GetXaxis()->SetTitle("dilepton mass [GeV]");
  hsf->Draw();
  hem->SetLineColor(4);
  hem->Draw("same");

  TLine line;
  line.SetLineColor(2);
  line.SetLineStyle(2);

  line.DrawLine(12 ,1,12,1e5);
  line.DrawLine(70 ,1,70,1e5);
  line.DrawLine(81 ,1,81,1e5);
  line.DrawLine(101,1,101,1e5);

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->AddEntry(hsf,"SF (ee+#mu#mu)","l");
  leg->AddEntry(hem,"OF (e#mu)","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  float nsf_peak = ch->GetEntries(sel+(ee||mm)+peak);
  float nsf_lo   = ch->GetEntries(sel+(ee||mm)+lo);

  float nof_peak = ch->GetEntries(sel+(em)+peak);
  float nof_lo   = ch->GetEntries(sel+(em)+lo);

  float Rtot = ( nsf_lo - nof_lo ) / ( nsf_peak - nof_peak );

  cout << "R " << Rtot << endl;

  const unsigned int nbins = 6;

  TCut metcuts[nbins];
  metcuts[0] = TCut("pfmet>0  && pfmet<10");
  metcuts[1] = TCut("pfmet>10 && pfmet<20");
  metcuts[2] = TCut("pfmet>20 && pfmet<30");
  metcuts[3] = TCut("pfmet>30 && pfmet<40");
  metcuts[4] = TCut("pfmet>40 && pfmet<50");
  metcuts[5] = TCut("pfmet>50 && pfmet<60");
  // metcuts[6] = TCut("pfmet>60 && pfmet<70");
  // metcuts[7] = TCut("pfmet>70 && pfmet<80");
  // metcuts[8] = TCut("pfmet>80 && pfmet<90");
  // metcuts[9] = TCut("pfmet>90 && pfmet<100");

  float R[nbins];
  float Rerr[nbins];

  float x[nbins];
  float xerr[nbins];


  for( int i = 0 ; i < nbins ; ++i ){

    x[i]     = 5.0 + 10.0 * i;
    xerr[i]  = 0.0;

    nsf_peak = ch->GetEntries(sel+(ee||mm)+metcuts[i]+peak);
    nsf_lo   = ch->GetEntries(sel+(ee||mm)+metcuts[i]+lo);
    
    nof_peak = ch->GetEntries(sel+(em)+metcuts[i]+peak);
    nof_lo   = ch->GetEntries(sel+(em)+metcuts[i]+lo);
    
    R[i]     = ( nsf_lo - nof_lo ) / ( nsf_peak - nof_peak );
    Rerr[i]  = ( sqrt( nsf_lo + nof_lo ) ) / ( nsf_peak - nof_peak );

  }

  TGraphErrors* gr = new TGraphErrors(nbins,x,R,xerr,Rerr);

  TCanvas *c3 = new TCanvas();
  c3->cd();
  gr->GetXaxis()->SetTitle("E_{T}^{miss} range [GeV]");
  gr->GetYaxis()->SetTitle("R_{low/in}");
  gr->Draw("AP");
  

  // float Ree      = nee_lo / nee_peak;

  // float nmm_peak = ch->GetEntries(sel+mm+peak);
  // float nmm_lo   = ch->GetEntries(sel+mm+lo);
  // float Rmm      = nmm_lo / nmm_peak;

  // cout << "Ree " << Ree << endl;
  // cout << "Rmm " << Rmm << endl;


  // TCanvas *c2 = new TCanvas();
  // c2->cd();
  // hee->DrawNormalized("hist");
  // hmm->SetLineColor(2);
  // hmm->DrawNormalized("samehist");







}
