{

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  char* path                    = "../output/V00-01-03";
  char* infilename		= Form("%s/ttall_smallTree*root",path);

  // char* ttlfilename		= Form("%s/ttl_smallTree.root",path);
  // char* ttllfilename		= Form("%s/ttll_smallTree.root",path);
  // char* ttltaufilename		= Form("%s/ttltau_smallTree.root",path);
  // char* tttaufilename		= Form("%s/tttau_smallTree.root",path);
  // char* tttautaufilename	= Form("%s/tttautau_smallTree.root",path);
  // char* ttotrfilename		= Form("%s/ttotr_smallTree.root",path);

  //--------------------------------------------------
  // list of output files
  //--------------------------------------------------

  char* fakefilename		= Form("%s/ttfake_smallTree.root",path);    // 0 W leptons
  char* slfilename		= Form("%s/ttsl_smallTree.root",path);      // 1 W lepton
  char* dlfilename		= Form("%s/ttdl_smallTree.root",path);      // 2 W leptons
  char* dllostfilename		= Form("%s/ttdllost_smallTree.root",path);  // 2nd lepton outside tracker acceptance
  char* dllepfilename		= Form("%s/ttdllep_smallTree.root",path);   // 2nd lepton = e, mu
  char* dltauhfilename		= Form("%s/ttdltauh_smallTree.root",path);  // 2nd lepton = tau->had
  char* dltaulfilename		= Form("%s/ttdltaul_smallTree.root",path);  // 2nd lepton = tau->e/mu
  char* dltauh1filename		= Form("%s/ttdltauh1_smallTree.root",path); // 2nd lepton = tau->had (1-prong)
  char* dltauhmfilename		= Form("%s/ttdltauhm_smallTree.root",path); // 2nd lepton = tau->had (multi-prong)

  // char* ttl			= "(w1==1||w1==2) && nleps==1 && ntaus==0";
  // char* ttll			= "(w1==1||w1==2) && nleps==2 && ntaus==0";
  // char* ttltau			= "(w1==1||w1==2) && nleps==2 && ntaus==1";
  // char* tttau			= "(w1==3||w1==4) && nleps==1 && ntaus==1";
  // char* tttautau		= "(w1==3||w1==4) && nleps==2 && ntaus==2";
  // TCut  ttotrcut		= !( TCut(ttl) || TCut(ttll) || TCut(ttltau) || TCut(tttau) || TCut(tttautau) );
  // char* ttotr			= ttotrcut.GetTitle();

  //--------------------------------------------------
  // list of cuts definining output files
  //--------------------------------------------------
  
  char* fake                    = "nleps==0";
  char* sl                      = "nleps==1";
  char* dl                      = "nleps==2";
  char* dllost                  = "nleps==2 && abs(mclep2.eta())>2.5";
  char* dllep                   = "nleps==2 && abs(mclep2.eta())<2.5 && abs(mcid2)<14";
  char* dltauh                  = "nleps==2 && abs(mclep2.eta())<2.5 && abs(mcid2)>14 && mcdecay2==1";
  char* dltaul                  = "nleps==2 && abs(mclep2.eta())<2.5 && abs(mcid2)>14 && mcdecay2==2";
  char* dltauh1                 = "nleps==2 && abs(mclep2.eta())<2.5 && abs(mcid2)>14 && mcdecay2==1 && mcndec2==1";
  char* dltauhm                 = "nleps==2 && abs(mclep2.eta())<2.5 && abs(mcid2)>14 && mcdecay2==1 && mcndec2 >1";

  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------

  cout << "Reading in : " << infilename                            << endl;
  // cout << "Writing to : " << ttlfilename        << " " << ttl      << endl;
  // cout << "Writing to : " << ttllfilename       << " " << ttll     << endl;
  // cout << "Writing to : " << ttltaufilename     << " " << ttltau   << endl;
  // cout << "Writing to : " << tttaufilename      << " " << tttau    << endl;
  // cout << "Writing to : " << tttautaufilename   << " " << tttautau << endl;
  // cout << "Writing to : " << ttotrfilename      << " " << ttotr    << endl;
  cout << "Writing to : " << fakefilename       << " " << fake      << endl;
  cout << "Writing to : " << slfilename         << " " << sl        << endl;
  cout << "Writing to : " << dlfilename         << " " << dl        << endl;
  cout << "Writing to : " << dllostfilename     << " " << dllost    << endl;
  cout << "Writing to : " << dllepfilename      << " " << dllep     << endl;
  cout << "Writing to : " << dltauhfilename     << " " << dltauh    << endl;
  cout << "Writing to : " << dltaulfilename     << " " << dltaul    << endl;
  cout << "Writing to : " << dltauh1filename    << " " << dltauh1   << endl;
  cout << "Writing to : " << dltauhmfilename    << " " << dltauhm   << endl;

  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  /*
  //-------------------
  // ttl
  //-------------------

  TFile *ttlfile = TFile::Open(ttlfilename, "RECREATE");
  assert( ttlfile != 0 );
  TTree* ttltree = chain->CopyTree( ttl );
  ttltree->Write();
  ttlfile->Close();

  //-------------------
  // ttll
  //-------------------

  TFile *ttllfile = TFile::Open(ttllfilename, "RECREATE");
  assert( ttllfile != 0 );
  TTree* ttlltree = chain->CopyTree( ttll );
  ttlltree->Write();
  ttllfile->Close();

  //-------------------
  // ttltau
  //-------------------

  TFile *ttltaufile = TFile::Open(ttltaufilename, "RECREATE");
  assert( ttltaufile != 0 );
  TTree* ttltautree = chain->CopyTree( ttltau );
  ttltautree->Write();
  ttltaufile->Close();

  //-------------------
  // tttau
  //-------------------

  TFile *tttaufile = TFile::Open(tttaufilename, "RECREATE");
  assert( tttaufile != 0 );
  TTree* tttautree = chain->CopyTree( tttau );
  tttautree->Write();
  tttaufile->Close();

  //-------------------
  // tttautau
  //-------------------

  TFile *tttautaufile = TFile::Open(tttautaufilename, "RECREATE");
  assert( tttautaufile != 0 );
  TTree* tttautautree = chain->CopyTree( tttautau );
  tttautautree->Write();
  tttautaufile->Close();

  //-------------------
  // ttotr
  //-------------------

  TFile *ttotrfile = TFile::Open(ttotrfilename, "RECREATE");
  assert( ttotrfile != 0 );
  TTree* ttotrtree = chain->CopyTree( ttotr );
  ttotrtree->Write();
  ttotrfile->Close();
  */

  //-------------------
  // fake
  //-------------------

  TFile *fakefile = TFile::Open(fakefilename, "RECREATE");
  assert( fakefile != 0 );
  TTree* faketree = chain->CopyTree( fake );
  faketree->Write();
  fakefile->Close();

  //-------------------
  // 1-lepton
  //-------------------

  TFile *slfile = TFile::Open(slfilename, "RECREATE");
  assert( slfile != 0 );
  TTree* sltree = chain->CopyTree( sl );
  sltree->Write();
  slfile->Close();

  //-------------------
  // 2-lepton
  //-------------------

  TFile *dlfile = TFile::Open(dlfilename, "RECREATE");
  assert( dlfile != 0 );
  TTree* dltree = chain->CopyTree( dl );
  dltree->Write();
  dlfile->Close();

  //--------------------------------------------------------------
  // 2-lepton (2nd lepton lost)
  //--------------------------------------------------------------

  TFile *dllostfile = TFile::Open(dllostfilename, "RECREATE");
  assert( dllostfile != 0 );
  TTree* dllosttree = chain->CopyTree( dllost );
  dllosttree->Write();
  dllostfile->Close();

  //--------------------------------------------------------------
  // 2-lepton (2nd lepton = e, mu)
  //--------------------------------------------------------------

  TFile *dllepfile = TFile::Open(dllepfilename, "RECREATE");
  assert( dllepfile != 0 );
  TTree* dlleptree = chain->CopyTree( dllep );
  dlleptree->Write();
  dllepfile->Close();

  //--------------------------------------------------------------
  // 2-lepton (2nd lepton = tau->had)
  //--------------------------------------------------------------

  TFile *dltauhfile = TFile::Open(dltauhfilename, "RECREATE");
  assert( dltauhfile != 0 );
  TTree* dltauhtree = chain->CopyTree( dltauh );
  dltauhtree->Write();
  dltauhfile->Close();

  //--------------------------------------------------------------
  // 2-lepton (2nd lepton = tau->e,mu)
  //--------------------------------------------------------------

  TFile *dltaulfile = TFile::Open(dltaulfilename, "RECREATE");
  assert( dltaulfile != 0 );
  TTree* dltaultree = chain->CopyTree( dltaul );
  dltaultree->Write();
  dltaulfile->Close();
  //--------------------------------------------------------------
  // 2-lepton (2nd lepton = tau->had->1-prong)
  //--------------------------------------------------------------

  TFile *dltauh1file = TFile::Open(dltauh1filename, "RECREATE");
  assert( dltauh1file != 0 );
  TTree* dltauh1tree = chain->CopyTree( dltauh1 );
  dltauh1tree->Write();
  dltauh1file->Close();

  //--------------------------------------------------------------
  // 2-lepton (2nd lepton = tau->had->multi-prong)
  //--------------------------------------------------------------

  TFile *dltauhmfile = TFile::Open(dltauhmfilename, "RECREATE");
  assert( dltauhmfile != 0 );
  TTree* dltauhmtree = chain->CopyTree( dltauhm );
  dltauhmtree->Write();
  dltauhmfile->Close();






}
