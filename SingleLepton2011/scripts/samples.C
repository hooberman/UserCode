{

  char* path = "../output/V00-00-00";

  //--------------------------------------------------
  // declare samples
  //--------------------------------------------------

  cout << "Loading samples at " << path << endl;

  TChain *data = new TChain("t");
  data->Add(Form("%s/data_smallTree.root",path));

  TChain *datasf = new TChain("t");
  datasf->Add(Form("%s/data_smallTree_singleFake.root",path));

  TChain *datadf = new TChain("t");
  datadf->Add(Form("%s/data_smallTree_doubleFake.root",path));

  TChain *ttdil = new TChain("t");
  ttdil->Add(Form("%s/ttdil_smallTree.root",path));

  TChain *ttpowheg = new TChain("t");
  ttpowheg->Add(Form("%s/ttpowheg_smallTree.root",path));

  TChain *ttall = new TChain("t");
  ttall->Add(Form("%s/ttall_smallTree.root",path));

  TChain *ttotr = new TChain("t");
  ttotr->Add(Form("%s/ttotr_smallTree.root",path));

  TChain *dy = new TChain("t");
  dy->Add(Form("%s/DYtot_smallTree.root",path));
  
  TChain *ww = new TChain("t");
  ww->Add(Form("%s/ww_smallTree.root",path));

  TChain *wz = new TChain("t");
  wz->Add(Form("%s/wz_smallTree.root",path));

  TChain *zz = new TChain("t");
  zz->Add(Form("%s/zz_smallTree.root",path));

  TChain *zjets = new TChain("t");
  zjets->Add(Form("%s/Zjets_smallTree.root",path));

  TChain *t = new TChain("t");
  t->Add(Form("%s/tW_smallTree.root",path));

  TChain *wjets = new TChain("t");
  wjets->Add(Form("%s/wjets_smallTree.root",path));

  TChain *dyeemm = new TChain("t");
  dyeemm->Add(Form("%s/DYee_smallTree.root",path));
  dyeemm->Add(Form("%s/DYmm_smallTree.root",path));

  TChain *dyee = new TChain("t");
  dyee->Add(Form("%s/DYee_smallTree.root",path));

  TChain *dymm = new TChain("t");
  dymm->Add(Form("%s/DYmm_smallTree.root",path));

  TChain *dytt = new TChain("t");
  dytt->Add(Form("%s/DYtautau_smallTree.root",path));

  TChain *LM0 = new TChain("t");
  LM0->Add(Form("%s/LM0_smallTree.root",path));

  TChain *LM1 = new TChain("t");
  LM1->Add(Form("%s/LM1_smallTree.root",path));

  TChain *LM2 = new TChain("t");
  LM2->Add(Form("%s/LM2_smallTree.root",path));

  TChain *LM3 = new TChain("t");
  LM3->Add(Form("%s/LM3_smallTree.root",path));

  TChain *LM4 = new TChain("t");
  LM4->Add(Form("%s/LM4_smallTree.root",path));

  TChain *LM5 = new TChain("t");
  LM5->Add(Form("%s/LM5_smallTree.root",path));

  TChain *LM6 = new TChain("t");
  LM6->Add(Form("%s/LM6_smallTree.root",path));

  TChain *t1lh = new TChain("t");
  t1lh->Add(Form("%s/T1lh_smallTree.root",path));

  TChain *T2tt = new TChain("t");
  T2tt->Add(Form("%s/T2tt_smallTree.root",path));

  TChain *mc = new TChain("t");
  mc->Add(Form("%s/ttdil_smallTree.root",path));
  mc->Add(Form("%s/ttotr_smallTree.root",path));
  mc->Add(Form("%s/DYtot_smallTree.root",path));
  mc->Add(Form("%s/ww_smallTree.root",path));
  mc->Add(Form("%s/wz_smallTree.root",path));
  mc->Add(Form("%s/zz_smallTree.root",path));
  mc->Add(Form("%s/tW_smallTree.root",path));
  mc->Add(Form("%s/wjets_smallTree.root",path));
                          
  //--------------------------------------------------
  // declare selection
  //--------------------------------------------------

  TCut zveto("passz == 0");
  TCut njets2("npfjets >= 2");
  TCut njets2("npfjets >= 3");
  TCut met50("pfmet > 50");
  TCut met100("pfmet > 100");
  TCut ht100("htpf > 100");
  TCut ht200("htpf > 200");
  TCut pt1010("lep1.pt()>10 && lep2.pt()>10");
  TCut pt2010("lep1.pt()>20 && lep2.pt()>10");
  TCut pt105 ("lep1.pt()>10 && lep2.pt()>5" );
  TCut pt55  ("lep1.pt()>5 && lep2.pt()>5" );
  TCut A("ht>125 && ht<300 && y>8.5");
  TCut B("ht>125 && ht<300 && y>4.5 && y<8.5");
  TCut C("ht>300 && y>4.5 && y<8.5");
  TCut D("ht>300 && y>8.5");
  TCut lep2tight("lep2.pt()>10 || isont2<0.15");

  //TCut weight("weight*ndavtxweight");
  TCut weight("weight*ndavtxweight*trgeff*0.976");

  TCut  sel  = njets3 + met25;
  char* colformat = "precision=3 col=7.6::10.10:";

  cout << "sel       : " << sel.GetTitle()       << endl;
  cout << "weight    : " << weight.GetTitle()    << endl;
  cout << "colformat : " << colformat            << endl;
  cout << endl << endl;

}
