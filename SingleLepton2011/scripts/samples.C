{

  char* path = "../output/V00-00-01";
  cout << "Loading samples at " << path << endl;

  //--------------------------------------------------
  // declare samples
  //--------------------------------------------------

  TChain *data = new TChain("t");
  data->Add(Form("%s/data_smallTree.root",path));

  TChain *ttall = new TChain("t");
  ttall->Add(Form("%s/ttall_smallTree.root",path));

  TChain *ttl = new TChain("t");
  ttl->Add(Form("%s/ttl_smallTree.root",path));

  TChain *ttll = new TChain("t");
  ttll->Add(Form("%s/ttll_smallTree.root",path));

  TChain *ttltau = new TChain("t");
  ttltau->Add(Form("%s/ttltau_smallTree.root",path));

  TChain *tttau = new TChain("t");
  tttau->Add(Form("%s/tttau_smallTree.root",path));

  TChain *tttautau = new TChain("t");
  tttautau->Add(Form("%s/tttautau_smallTree.root",path));

  TChain *ttotr = new TChain("t");
  ttotr->Add(Form("%s/ttotr_smallTree.root",path));

  TChain *wjets = new TChain("t");
  wjets->Add(Form("%s/wjets_smallTree.root",path));
  
  //--------------------------------------------------
  // declare selection
  //--------------------------------------------------

  TCut njets3("njets >= 3");
  TCut met60("pfmet > 60");
  TCut ht300("ht > 300");

  TCut weight("weight*0.976");

  TCut  sel       = njets3 + met60 + ht300;
  char* colformat = "precision=3 col=7.6::10.10:";

  cout << "sel       : " << sel.GetTitle()       << endl;
  cout << "weight    : " << weight.GetTitle()    << endl;
  cout << "colformat : " << colformat            << endl;
  
  cout << endl << endl;

}
