{

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  char* path                    = "../output/V00-00-09";
  char* infilename		= Form("%s/Zjets_smallTree*root",path);

  //--------------------------------------------------
  // list of output files
  //--------------------------------------------------

  char* outfilename		= Form("%s/Zjets_2jets_smallTree.root",path);    // 0 W leptons

  //--------------------------------------------------
  // list of cuts definining output files
  //--------------------------------------------------
  
  char* njets2                  = "njets>=2 || njetsuncor>=2 || njetsplain>=2";

  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------

  cout << "Reading in : " << infilename                            << endl;
  cout << "Writing to : " << outfilename    << " " << njets2       << endl;

  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  //--------------------------------------------------------------
  // 2-lepton (2nd lepton = tau->had->multi-prong)
  //--------------------------------------------------------------

  TFile *outfile = TFile::Open(outfilename, "RECREATE");
  assert( outfile != 0 );
  TTree* outtree = chain->CopyTree( njets2 );
  outtree->Write();
  outfile->Close();






}
