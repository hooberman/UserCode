{
  char* path = "output/promptreco/dcsonly";

  TChain *dymm = new TChain("Events");
  dymm->Add(Form("%s/dymm_default_baby.root",path));
  
  TChain *dyee = new TChain("Events");
  dyee->Add(Form("%s/dyee_default_baby.root",path));
  
  TChain *zjets = new TChain("Events");
  zjets->Add(Form("%s/zjets_default_baby.root",path));

  TChain *tt = new TChain("Events");
  tt->Add(Form("%s/tt_default_baby.root",path));

  TChain *ww = new TChain("Events");
  ww->Add(Form("%s/ww_default_baby.root",path));
  
  TChain *h130 = new TChain("Events");
  h130->Add(Form("%s/h130_default_baby.root",path));

  TChain *chdata = new TChain("Events");
  chdata->Add(Form("%s/data_default_baby.root",path));

  TChain *dy = new TChain("Events");
  dy->Add(Form("%s/dymm_default_baby.root",path));
  dy->Add(Form("%s/dyee_default_baby.root",path));
  //dy->Add(Form("%s/dymm_default_tenpercent_baby.root",path));
  //dy->Add(Form("%s/dyee_default_tenpercent_baby.root",path));

  TCut zmass("dilep.mass()>76&&dilep.mass()<106");
  TCut jet1pv("njets30==1 && jetpv==1");
  TCut jet1("njets30==1");
  TCut jet0("njets30==0");
  TCut ee("leptype==3");
  TCut mm("leptype==0");
  TCut sf=ee||mm;
  TCut em("leptype==1||leptype==2");
  TCut dphizjet("acos(cos(dilep.phi()-jet.phi()))>2.");
  TCut lepveto("nlep==0");
  TCut muveto("softmu==0");
  TCut topveto("toptag==0");
  TCut weight("weight*davtxweight"); 

  //TCut fullsel= zmass+jet1+sf+lepveto+topveto+muveto;
  TCut fullsel= zmass+jet1+sf;

}
