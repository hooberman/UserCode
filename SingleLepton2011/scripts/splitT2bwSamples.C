{

  char* path = "../output/V00-01-03";
  //char* path = "../output/temp";

  //-------------------------------------
  // declare input file, tree, etc.
  //-------------------------------------

  char* infilename = Form("%s/T2bw_smallTree.root" , path);
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  //-------------------------------------
  // which grid points to skim?
  //-------------------------------------
  
  vector<string> filenames;
  vector<int>   mg;
  vector<int>   ml;
  vector<double> x;
  
  mg.push_back(250);   ml.push_back(50);   x.push_back(0.5);
  mg.push_back(250);   ml.push_back(50);   x.push_back(0.75);

  const unsigned int npoints = mg.size();

  //-------------------------------------
  // skim grid points
  //-------------------------------------

  TFile *file[npoints];
  TTree *tree[npoints];

  for( unsigned int i = 0 ; i < npoints ; ++i ){

    char* filename = Form("%s/T2bw_%i_%i_%.2f_smallTree.root",path,mg.at(i),ml.at(i),x.at(i));
    TCut cut(Form("mg==%i && ml==%i && abs(x-%.2f)<0.001",mg.at(i),ml.at(i),x.at(i)));

    cout << "Skimming: " << filename << " selection: " << cut.GetTitle() << endl;
    
    file[i] = TFile::Open(filename, "RECREATE");
    assert( file[i] != 0 );
    tree[i] = chain->CopyTree( cut );
    tree[i]->Write();
    file[i]->Close();
    
  }

}
