{

  char* iter = "V00-00-17";
  //char* iter = "temp";

  TCut pfleptons("pflep1.pt() > 20 && pflep2.pt() > 20");
  TCut transveto("el1tv==0 && el2tv==0");
  TCut Zmass("dilmasspf>81 && dilmasspf<101");
  TCut njets2("njets>=2");
  //TCut ee("leptype==0 && jetptll-ptll>-5 && jetptlt-ptlt>-5");
  TCut ee("leptype==0 && (ee==1 || isdata==0)");
  TCut mm("leptype==1 && (mm==1 || isdata==0)");
  TCut em("leptype==2 && (em==1 || me==1 || isdata==0)");
  TCut sf = ee||mm;
  TCut bveto("nbm==0");
  TCut nlep2("nlep==2");
  TCut mjj("mjj>70.0 && mjj<110.0");

  TCut sel;
  sel += njets2;
  sel += pfleptons;
  //sel += transveto;
  sel += Zmass;
  sel += (ee||mm||em);
  sel += bveto;
  sel += nlep2;
  sel += mjj;

  TCut weight("weight * 5.10 * vtxweight");

  cout << "Baby location         : " << iter              << endl;
  cout << "Baseline selection    : " << sel.GetTitle()    << endl;
  cout << "ee selection          : " << ee.GetTitle()     << endl;
  cout << "mm selection          : " << mm.GetTitle()     << endl;
  cout << "em selection          : " << em.GetTitle()     << endl;
  cout << "weight                : " << weight.GetTitle() << endl;

  char* colformat = "precision=4 col=7.6::10.10:";

  TChain *tt = new TChain("T1");
  tt->Add(Form("../output/%s/ttbar_baby.root",iter));

  TChain *zjets = new TChain("T1");
  zjets->Add(Form("../output/%s/zjets_baby.root",iter));

  TChain *wz = new TChain("T1");
  wz->Add(Form("../output/%s/wz_baby.root",iter));

  TChain *wzmg = new TChain("T1");
  wzmg->Add(Form("../output/%s/wzmg_baby.root",iter));

  TChain *zz = new TChain("T1");
  zz->Add(Form("../output/%s/zz_baby.root",iter));

  TChain *data = new TChain("T1");
  data->Add(Form("../output/%s/data_baby.root",iter));

  TChain *dataskim = new TChain("T1");
  dataskim->Add(Form("../output/%s/dataskim_baby.root",iter));


}
