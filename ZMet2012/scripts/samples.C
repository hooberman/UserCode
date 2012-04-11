{

  char* iter = "V00-00-00";
  //char* iter = "temp";

  TCut zmass("dilmass>81 && dilmass<101");
  TCut ee     = "leptype==0";
  TCut mm     = "leptype==1";
  TCut em     = "leptype==2";
  TCut sel    = zmass+(mm||ee);
  TCut selall = zmass+(mm||ee||em);

  TCut weight("weight * 0.032");

  TChain *tt = new TChain("T1");
  tt->Add(Form("../output/%s/ttbar_baby.root",iter));

  TChain *zjets = new TChain("T1");
  zjets->Add(Form("../output/%s/zjets_baby.root",iter));

  TChain *data = new TChain("T1");
  data->Add(Form("../output/%s/data_baby.root",iter));


}
