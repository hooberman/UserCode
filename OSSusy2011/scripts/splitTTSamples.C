{

  //char* path = "../output/temp/highpt";
  char* path = "../output/V00-02-01/highpt";

  char* infilename     = Form("%s/ttall_smallTree.root",path);
  char* ttllfilename   = Form("%s/ttll_smallTree.root",path);
  char* tttaufilename  = Form("%s/tttau_smallTree.root",path);
  char* ttfakefilename = Form("%s/ttfake_smallTree.root",path);

  char* ll   = "(w1>0&&w2>0) && (w1<3&&w2<3)";
  char* tau  = "(w1>0&&w2>0) && (w1>2||w2>2)";
  char* fake = "!(w1>0&&w2>0)";

  cout << "Reading in : " << infilename     << endl;
  cout << "Writing to : " << ttllfilename   << " " << ll   << endl;
  cout << "Writing to : " << tttaufilename  << " " << tau  << endl;
  cout << "Writing to : " << ttfakefilename << " " << fake << endl;
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  //-------------------
  // ttll
  //-------------------

  TFile *llfile = TFile::Open(ttllfilename, "RECREATE");
  assert( llfile != 0 );
  TTree* lltree = chain->CopyTree( ll );
  lltree->Write();
  llfile->Close();

  //-------------------
  // tttau
  //-------------------

  TFile *taufile = TFile::Open(tttaufilename, "RECREATE");
  assert( taufile != 0 );
  TTree* tautree = chain->CopyTree( tau );
  tautree->Write();
  taufile->Close();

  //-------------------
  // ttll
  //-------------------

  TFile *fakefile = TFile::Open(ttfakefilename, "RECREATE");
  assert( fakefile != 0 );
  TTree* faketree = chain->CopyTree( fake );
  faketree->Write();
  fakefile->Close();





}
