{
  //--------------------------------------------------
  // Flavor history information
  //--------------------------------------------------
  // from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFlavorHistory
  // The hierarchy is:
  // 1. W+bb with >= 2 jets from the ME (dr > 0.5)
  // 2. W+b or W+bb with 1 jet from the ME
  // 3. W+cc from the ME (dr > 0.5)
  // 4. W+c or W+cc with 1 jet from the ME
  // 5. W+bb with 1 jet from the parton shower (dr == 0.0)
  // 6. W+cc with 1 jet from the parton shower (dr == 0.0)
  // 7. W+bb with >= 2 partons but 1 jet from the ME (dr == 0.0)
  // 8. W+cc with >= 2 partons but 1 jet from the ME (dr == 0.0)
  // 9. W+bb with >= 2 partons but 2 jets from the PS (dr > 0.5)
  // 10. W+cc with >= 2 partons but 2 jets from the PS (dr > 0.5)
  // 11. Veto of all the previous (W+ light jets)
  // Note that paths 1-6, plus sample 11, are the paths you will use in your analysis. 
  // Paths 7-10 are not included in any of the samples as this is the double counting that we are trying to eliminate.
  // The physics content of the samples for your analysis should break down as follows:
  // W+bb :
  // From the "V+QQ" sample: Paths 1 and 2 (for the matrix element contribution).
  // From the "W+jets" sample: Path 5 (for the gluon splitting contribution).
  // W+cc :
  // From the "V+QQ" sample: Paths 3 and 4 (for the matrix element contribution).
  // From the "W+Jets" sample: Path 6 (for the gluon splitting contribution).
  // W+c
  // From the "W+c" sample: Path 4 (for true flavor excitation).
  // W+light flavor:
  // From the "W+jets" sample: Path 11 (veto of all heavy flavor content).
  // The remainder of the paths (7-10) are defined for studies and completeness. They contain samples of phase space that are unreliable (such as using the matrix element calculation in the "colinear" case, or using the parton shower in the case of widely separated partons), and hence the events should be discarded for the analysis. Note that the normalization still must be taken into account when stitching samples back together.

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  char* path                    = "../output/V00-03-00_flavtest";
  char* infilename		= Form("%s/wjets_smallTree*root",path);

  //--------------------------------------------------
  // list of output files
  //--------------------------------------------------

  char* wbxfilename		= Form("%s/wbx_smallTree.root",path);    // W + bx
  char* wcxfilename		= Form("%s/wcx_smallTree.root",path);    // W + cx
  char* wlffilename		= Form("%s/wlf_smallTree.root",path);    // W + light flavor

  //--------------------------------------------------
  // list of cuts definining output files
  //--------------------------------------------------
  
  char* wbx                     = "wflav==1 || wflav==2 || wflav==5";
  char* wcx                     = "wflav==3 || wflav==4 || wflav==6";
  char* wlf                     = "wflav==11";

  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------

  cout << "Reading in : " << infilename                            << endl;
  cout << "Writing to : " << wbxfilename        << " " << wbx      << endl;
  cout << "Writing to : " << wcxfilename        << " " << wcx      << endl;
  cout << "Writing to : " << wlffilename        << " " << wlf      << endl;

  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  //-------------------
  // W + bx
  //-------------------

  TFile *wbxfile = TFile::Open(wbxfilename, "RECREATE");
  assert( wbxfile != 0 );
  TTree* wbxtree = chain->CopyTree( wbx );
  wbxtree->Write();
  wbxfile->Close();

  //-------------------
  // W + cx
  //-------------------

  TFile *wcxfile = TFile::Open(wcxfilename, "RECREATE");
  assert( wcxfile != 0 );
  TTree* wcxtree = chain->CopyTree( wcx );
  wcxtree->Write();
  wcxfile->Close();

  //-------------------
  // W + lf
  //-------------------

  TFile *wlffile = TFile::Open(wlffilename, "RECREATE");
  assert( wlffile != 0 );
  TTree* wlftree = chain->CopyTree( wlf );
  wlftree->Write();
  wlffile->Close();


}
