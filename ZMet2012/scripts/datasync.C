{

  char* iter = "V00-01-00";
  //char* iter = "temp";

  TCut pfleptons("pflep1.pt() > 20 && pflep2.pt() > 20");
  TCut transveto("el1tv==0 && el2tv==0");
  TCut Zmasspf("dilmasspf>81 && dilmasspf<101");
  TCut Zmass("dilmass>81 && dilmass<101");
  TCut njets2("njets>=2");

  TCut pt2020("lep1.pt()>20.0 && lep2.pt()>20.0");
  TCut nlep2("lep3.pt()<20.0");

  TCut ee("leptype==0");
  TCut mm("leptype==1");
  TCut em("leptype==2");
  TCut eetrg("ee==1");
  TCut mmtrg("mm==1");
  TCut emtrg("em==1 || me==1");
  TCut njets0("njets==0");
  TCut njets1("njets==1");
  TCut njets2("njets>=2");
  TCut filters("csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1");

  //---------------------------------------------------
  // RelValZEE
  //---------------------------------------------------

  TChain *chmcee = new TChain("T1");
  chmcee->Add(Form("../output/%s/RelValZEE_baby_pt2020.root",iter));

  cout << "RelValZEE" << endl;
  cout << "baseline (no filters)    " << chmcee->GetEntries(ee+pt2020+nlep2)                << endl;
  cout << "baseline                 " << chmcee->GetEntries(ee+pt2020+nlep2+filters)        << endl;
  cout << "baseline + 0 jets        " << chmcee->GetEntries(ee+pt2020+nlep2+filters+njets0) << endl;
  cout << "baseline + 1 jet         " << chmcee->GetEntries(ee+pt2020+nlep2+filters+njets1) << endl;
  cout << "baseline + >=2 jets      " << chmcee->GetEntries(ee+pt2020+nlep2+filters+njets2) << endl;
  cout << "baseline + pfmet>10      " << chmcee->GetEntries(ee+pt2020+nlep2+filters+"pfmet>10.0") << endl;
  cout << "baseline + T1 pfmet>10   " << chmcee->GetEntries(ee+pt2020+nlep2+filters+"pfmett1>10.0") << endl;

  //---------------------------------------------------
  // RelValZMM data
  //---------------------------------------------------

  TChain *chmcmm = new TChain("T1");
  chmcmm->Add(Form("../output/%s/RelValZMM_baby_pt2020.root",iter));

  cout << "RelValZMM" << endl;
  cout << "baseline (no filters)    " << chmcmm->GetEntries(mm+pt2020+nlep2)                << endl;
  cout << "baseline                 " << chmcmm->GetEntries(mm+pt2020+nlep2+filters)        << endl;
  cout << "baseline + 0 jets        " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+njets0) << endl;
  cout << "baseline + 1 jet         " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+njets1) << endl;
  cout << "baseline + >=2 jets      " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+njets2) << endl;
  cout << "baseline + pfmet>10      " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+"pfmet>10.0") << endl;
  cout << "baseline + T1 pfmet>10   " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+"pfmett1>10.0") << endl;

  // cout << "baseline + 0 b           " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+"nbcsvm==0") << endl;
  // cout << "baseline + trigger       " << chmcmm->GetEntries(mm+pt2020+nlep2+filters+"mm==1")      << endl;

  //---------------------------------------------------
  // DoubleElectron data
  //---------------------------------------------------

  TChain *chee = new TChain("T1");
  chee->Add(Form("../output/%s/DoubleElectron_199752_baby_pt2020.root",iter));

  cout << "DoubleEle" << endl;
  cout << "baseline (no filters)    " << chee->GetEntries(ee+pt2020+nlep2)                << endl;
  cout << "baseline                 " << chee->GetEntries(ee+pt2020+nlep2+filters)        << endl;
  cout << "baseline + 0 b           " << chee->GetEntries(ee+pt2020+nlep2+filters+"nbcsvm==0") << endl;
  cout << "baseline + 0 jets        " << chee->GetEntries(ee+pt2020+nlep2+filters+njets0) << endl;
  cout << "baseline + 1 jet         " << chee->GetEntries(ee+pt2020+nlep2+filters+njets1) << endl;
  cout << "baseline + >=2 jets      " << chee->GetEntries(ee+pt2020+nlep2+filters+njets2) << endl;
  cout << "baseline + pfmet>10      " << chee->GetEntries(ee+pt2020+nlep2+filters+"pfmet>10.0") << endl;
  cout << "baseline + trigger       " << chee->GetEntries(ee+pt2020+nlep2+filters+"ee==1")      << endl;

  //---------------------------------------------------
  // DoubleMu data
  //---------------------------------------------------

  TChain *chmm = new TChain("T1");
  chmm->Add(Form("../output/%s/DoubleMu_199752_baby_pt2020.root",iter));

  cout << "DoubleMu" << endl;
  cout << "baseline (no filters)    " << chmm->GetEntries(mm+pt2020+nlep2)                << endl;
  cout << "baseline                 " << chmm->GetEntries(mm+pt2020+nlep2+filters)        << endl;
  cout << "baseline + 0 b           " << chmm->GetEntries(mm+pt2020+nlep2+filters+"nbcsvm==0") << endl;
  cout << "baseline + 0 jets        " << chmm->GetEntries(mm+pt2020+nlep2+filters+njets0) << endl;
  cout << "baseline + 1 jet         " << chmm->GetEntries(mm+pt2020+nlep2+filters+njets1) << endl;
  cout << "baseline + >=2 jets      " << chmm->GetEntries(mm+pt2020+nlep2+filters+njets2) << endl;
  cout << "baseline + pfmet>10      " << chmm->GetEntries(mm+pt2020+nlep2+filters+"pfmet>10.0") << endl;
  cout << "baseline + trigger       " << chmm->GetEntries(mm+pt2020+nlep2+filters+"mm==1")      << endl;

  //---------------------------------------------------
  // MuEG data
  //---------------------------------------------------

  TChain *chem = new TChain("T1");
  chem->Add(Form("../output/%s/MuEG_199752_baby_pt2020.root",iter));

  cout << "MuEG" << endl;
  cout << "baseline (no filters)    " << chem->GetEntries(em+pt2020+nlep2)                << endl;
  cout << "baseline                 " << chem->GetEntries(em+pt2020+nlep2+filters)        << endl;
  cout << "baseline + 0 b           " << chem->GetEntries(em+pt2020+nlep2+filters+"nbcsvm==0") << endl;
  cout << "baseline + 0 jets        " << chem->GetEntries(em+pt2020+nlep2+filters+njets0) << endl;
  cout << "baseline + 1 jet         " << chem->GetEntries(em+pt2020+nlep2+filters+njets1) << endl;
  cout << "baseline + >=2 jets      " << chem->GetEntries(em+pt2020+nlep2+filters+njets2) << endl;
  cout << "baseline + pfmet>10      " << chem->GetEntries(em+pt2020+nlep2+filters+"pfmet>10.0") << endl;
  cout << "baseline + trigger       " << chem->GetEntries(em+pt2020+nlep2+filters+"(em==1||me==1)")      << endl;


}
