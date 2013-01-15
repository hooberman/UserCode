{

  char* iter = "V00-01-04";
  //char* iter = "temp";

  TCut runrange("isdata==0 || (run<197556 || run>198913)");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut ee("leptype==0 && (ee==1 || isdata==0)");
  TCut mm("leptype==1 && (mm==1 || isdata==0)");
  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut Zmass("dilmass>81 && dilmass<101");
  TCut sf = ee||mm;
  TCut njets2("njets>=2");
  TCut bveto("nbcsvm==0");
  TCut mjj("mjj>70.0 && mjj<110.0");
  TCut nlep2("nlep==2");

  TCut sel;

  sel += runrange;
  sel += filters;
  sel += (ee||mm||em);
  sel += pt2020;
  sel += Zmass;
  sel += sf;
  sel += njets2;
  sel += bveto;
  sel += mjj;
  sel += nlep2;

  TCut weight("weight * 9.2");

  cout << "Baseline selection    : " << sel.GetTitle()    << endl;
  cout << "weight                : " << weight.GetTitle() << endl;

  TChain *gmsb_524 = new TChain("T1");
  gmsb_524->Add("../output/V00-01-08/gmsb_baby_oldIso.root");
  
  TChain *gmsb_526 = new TChain("T1");
  gmsb_526->Add("../output/V00-01-09/gmsb_526_baby_oldIso.root");
  
  const unsigned int npoints = 15;
  
  int muvalues[npoints] = {130,150,170,190,210,230,250,270,290,310,330,350,370,390,410};

  TH1F*    h524[npoints];
  TH1F*    h526[npoints];;
  TCanvas* can[npoints];
  TLegend* leg[npoints];
  TCut     mucut[npoints];

  float ptbins[]  = {0., 30., 60., 80., 100., 120., 150., 200., 300.};
  int   nptbin  = 8;

  TLatex *t = new TLatex();
  t->SetNDC();

  for( int i = 0 ; i < npoints ; ++i ){

    mucut[i] = TCut(Form("mg==%i",muvalues[i]));
    cout << "mu " << muvalues[i] << " " << mucut[i].GetTitle() << endl;

    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),600,600);
    can[i]->cd();

    h524[i] = new TH1F(Form("h524_%i",i),Form("h524_%i",i),nptbin,ptbins);
    h526[i] = new TH1F(Form("h526_%i",i),Form("h526_%i",i),nptbin,ptbins);

    h524[i]->Sumw2();
    h526[i]->Sumw2();

    gmsb_524->Draw(Form("min(pfmet,299)>>h524_%i",i),(sel+mucut[i])*weight);
    gmsb_526->Draw(Form("min(pfmet,299)>>h526_%i",i),(sel+mucut[i])*weight);

    h524[i]->SetLineColor(2);
    h526[i]->SetLineColor(4);
    h524[i]->SetMarkerColor(2);
    h526[i]->SetMarkerColor(4);
    h524[i]->SetMarkerStyle(21);
    h526[i]->SetMarkerStyle(25);

    h524[i]->SetMinimum(0);

    float max = h524[i]->GetMaximum();
    if( h526[i]->GetMaximum() > max ) max = h526[i]->GetMaximum();
    h524[i]->SetMaximum(1.1*max);

    h524[i]->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");

    h524[i]->Draw();
    h526[i]->Draw("same");
    
    leg[i] = new TLegend(0.7,0.7,0.9,0.9);
    leg[i]->AddEntry(h524[i],"524");
    leg[i]->AddEntry(h526[i],"526");
    leg[i]->SetBorderSize(0);
    leg[i]->SetFillColor(0);
    leg[i]->Draw();

    t->DrawLatex(0.6,0.5,Form("#mu = %i GeV",muvalues[i]));

    can[i]->Print(Form("../plots/GMSB_%i.pdf",muvalues[i]));
  }


}
