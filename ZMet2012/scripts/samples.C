{

  char* iter = "V00-01-04";
  //char* iter = "temp";

  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut pt2010("lep1.pt()>20.0 && lep2.pt()>10.0");
  TCut pfleptons("pflep1.pt() > 20 && pflep2.pt() > 20");
  TCut transveto("el1tv==0 && el2tv==0");
  TCut Zmasspf("dilmasspf>81 && dilmasspf<101");
  //TCut Zmass("dilmass>81 && dilmass<101");
  TCut Zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets>=2");
  TCut njets40_2("njets40>=2");
  TCut njets40_3("njets40>=3");
  TCut ht40_100("ht40>=100.0");
  TCut met100("pfmet>100.0");
  TCut met150("pfmet>150.0");
  //TCut ee("leptype==0 && jetptll-ptll>-5 && jetptlt-ptlt>-5");
  TCut ee("leptype==0 && (ee==1 || isdata==0)");
  TCut mm("leptype==1 && (mm==1 || isdata==0)");
  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut sf = ee||mm;
  TCut bveto("nbcsvm==0");
  TCut nlep2("nlep==2");
  TCut mjj("mjj>70.0 && mjj<110.0");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut runrange("isdata==0 || (run<197556 || run>198913)");

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

  //------------------------
  // high-MET SR (Aachen)
  //------------------------
  /*
  sel += runrange;
  sel += filters;
  sel += (ee||mm||em);
  sel += pt2010;
  sel += njets40_2;
  sel += ht40_100;
  sel += met150;
  */

  //------------------------
  // low-MET SR (ETH)
  //------------------------
  /*
  sel += runrange;
  sel += filters;
  sel += (ee||mm||em);
  sel += pt2020;
  sel += njets40_3;
  sel += met100;
  */

  //sel += Zmass;
  //sel += (ee||mm);

  // sel += "njets40>=2 && ht40>=100.0";
  // sel += pfleptons;
  // sel += transveto;
  // sel += Zmasspf;
  // sel += bveto;
  // sel += nlep2;
  // sel += mjj;

  TCut weight("weight * 9.2 * trgeff * vtxweight");

  cout << "Baby location         : " << iter              << endl;
  cout << "Baseline selection    : " << sel.GetTitle()    << endl;
  cout << "ee selection          : " << ee.GetTitle()     << endl;
  cout << "mm selection          : " << mm.GetTitle()     << endl;
  cout << "em selection          : " << em.GetTitle()     << endl;
  cout << "weight                : " << weight.GetTitle() << endl;

  char* colformat = "precision=4 col=7.6::10.10:";

  TChain *tt = new TChain("T1");
  tt->Add(Form("../output/%s/ttbar_53X_baby.root",iter));

  TChain *zjets = new TChain("T1");
  zjets->Add(Form("../output/%s/zjets_full_53X_baby_2jets.root",iter));

  TChain *wz = new TChain("T1");
  wz->Add(Form("../output/%s/wz_53X_baby.root",iter));

  TChain *wzmg = new TChain("T1");
  wzmg->Add(Form("../output/%s/wzmg_baby.root",iter));

  TChain *zz = new TChain("T1");
  zz->Add(Form("../output/%s/zz_53X_baby.root",iter));

  TChain *data = new TChain("T1");
  //data->Add(Form("../output/%s/data_baby.root",iter));
  //data->Add(Form("../output/%s/data_ALL_53X_baby_2jets.root",iter));
  //data->Add(Form("../output/%s/data_53X_baby_2jets.root",iter));
  data->Add(Form("../output/%s/data_53X_baby_2jets.root",iter));
  data->Add(Form("../output/%s/data_2012C_53X_baby_2jets.root",iter));

  TChain *dataskim = new TChain("T1");
  dataskim->Add(Form("../output/%s/dataskim2010_baby.root",iter));

  TChain *wzsms = new TChain("T1");
  wzsms->Add(Form("../output/%s/wzsms_baby_oldIso.root",iter));

  iter = "V00-01-05";
  cout << "SWITCHING TO " << iter << endl;

  TChain *ttw = new TChain("T1");
  ttw->Add(Form("../output/%s/ttW_53X_baby.root",iter));

  TChain *ttz = new TChain("T1");
  ttz->Add(Form("../output/%s/ttZ_53X_baby.root",iter));

  TChain *ttv = new TChain("T1");
  ttv->Add(Form("../output/%s/ttW_53X_baby.root",iter));
  ttv->Add(Form("../output/%s/ttZ_53X_baby.root",iter));

  TChain *vvv = new TChain("T1");
  vvv->Add(Form("../output/%s/VVV_53X_baby.root",iter));

  TChain *rare = new TChain("T1");
  rare->Add(Form("../output/%s/ttZ_53X_baby.root",iter));
  rare->Add(Form("../output/%s/VVV_53X_baby.root",iter));

  TChain *wzsms = new TChain("T1");
  wzsms->Add(Form("../output/%s/wzsms_baby_oldIso.root",iter));
  

}
