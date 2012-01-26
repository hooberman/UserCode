{

  char* prefix = "/hadoop/cms/store/user/benhoob/CMS2_V04-02-20-04/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1";

  const unsigned int nfiles = 2060;

  TChain *ch[nfiles];

  TCut sel("evt_run==167830 && evt_lumiBlock==44 && evt_event==40088398");


  for( unsigned int i = 2040 ; i < nfiles ; ++i ){

    cout << Form("Adding %s/ntuple_%i_*root",prefix,i) << endl;

    ch[i] = new TChain("Events");
    ch[i]->Add(Form("%s/ntuple_%i_*root",prefix,i));
    //ch[i]->Add(Form("%s/merged_ntuple_%i.root",prefix,i));
    //ch[i]->Add(Form("%s/skimmed_ntuple_%i.root",prefix,i));

    int nentries = ch[i]->GetEntries(sel);


    char* line="";
    if( nentries > 0 ) line = "        <<<------------------------------";
    cout << "File " << i << " entries " << nentries << line << endl;



  }

}
