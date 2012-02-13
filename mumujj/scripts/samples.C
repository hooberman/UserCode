{
  //--------------------------------
  // get Ben's histogram
  //--------------------------------

  char* path = "../output/V00-01-05";

  TChain *data2011a = new TChain("t");
  TChain *data2011b = new TChain("t");
  TChain *data2011  = new TChain("t");

  data2011a->Add(Form("%s/datamay10_smallTree.root",path));
  data2011a->Add(Form("%s/dataPRv4_smallTree.root",path));
  data2011a->Add(Form("%s/dataaug05_smallTree.root",path));
  data2011a->Add(Form("%s/dataPRv6_smallTree.root",path));
  data2011b->Add(Form("%s/data2011B_smallTree.root",path));

  data2011->Add(Form("%s/datamay10_smallTree.root",path));
  data2011->Add(Form("%s/dataPRv4_smallTree.root",path));
  data2011->Add(Form("%s/dataaug05_smallTree.root",path));
  data2011->Add(Form("%s/dataPRv6_smallTree.root",path));
  data2011->Add(Form("%s/data2011B_smallTree.root",path));

  //--------------------------------
  // compare histograms
  //--------------------------------

  TCut njets2("njets>=2");
  TCut hlt("hlt1==1 || hlt2==1");
  TCut leadjet120("jet1.pt()>120");
  TCut mll50("dilmass>50");
  TCut zveto("dilmass<70 || dilmass>100");
  TCut pt2520("lep1.pt()>25 && lep2.pt()>20");
  TCut pt3020("lep1.pt()>30 && lep2.pt()>20");
  TCut ngoodlep2("ngoodlep>=2");
  //TCut nbtags1("nbtags2024>=1");
  TCut nbtags1("nbtags2024c>=1");
  TCut trg("mtrg==1");
  TCut mm("leptype==1");

  // TCut pt3020("lep1.pt()>30 && lep2.pt()>20");
  // TCut nbtags1("nbtags33>=1");
  // TCut iso1("iso1<0.15");
  // TCut iso2("iso2<0.15");
  // TCut id1("passid1==1");
  // TCut id2("passid2==1");
  // TCut mu2("ngoodmu>=2");
  // TCut bump("mmjj>850 && mmjj<1150");
  // TCut left("mmjj>650 && mmjj<850");
  // TCut right("mmjj>1150");

  TCut sel;
  //sel += pt3020;
  sel += ngoodlep2;
  sel += pt2520;
  sel += hlt;
  sel += mll50;
  sel += zveto;
  sel += njets2;
  sel += leadjet120;
  sel += nbtags1;
  //sel += trg;
  //sel += mm;
  //sel += "mmjj > 600 && mmjj < 1400";
  //sel += "mmjj > 850 && mmjj < 1150";

  cout << "Using selection   " << sel.GetTitle()             << endl;
  cout << "Ben entries 2011a " << data2011a->GetEntries(sel) << endl;
  cout << "Ben entries 2011b " << data2011b->GetEntries(sel) << endl;
  cout << "Ben entries 2011  " << data2011->GetEntries(sel)  << endl;






}
