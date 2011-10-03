{

  //char* prefix = "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_4_V04-02-20-01/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/CMSSW_4_2_4_V04-02-20-01_merged/V04-02-20-01/";
  //char* prefix = "/nfs-4/userdata/spadhi/TAS/tanbeta10/dilepSim/";
  //char* prefix = "cms2_data/mSUGRA_m0-20to2000_m12-20to760_tanb-10andA0-0_7TeV-Pythia6Z_Summer11-PU_S4_START42_V11_FastSim-v1_AODSIM/V04-02-20-01/";
  //char* prefix = "/hadoop/cms/store/user/warren/CMS2_V04-02-20-03/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2/";
  //char* prefix = "/nfs-7/userdata/warren/SMS-T1Lh_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v2";
  //char* prefix = "/hadoop/cms/store/user/benhoob/CMS2_V04-02-29/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/";
  char* prefix = "/hadoop/cms/store/user/benhoob/CMS2_V04-02-20-04/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1";

  const unsigned int nfiles = 2060;

  TChain *ch[nfiles];

  //TCut LM1a("sparm_m0==60&&sparm_m12==240");
  //TCut LM1b("sparm_m0==60&&sparm_m12==260");
  //TCut LM6("sparm_m0==80&&sparm_m12==400");
  //TCut LM6("sparm_m0>1000.0");
  //TCut sel("sparm_mG==400 && sparm_mL==150");
  //TCut sel("sparm_mG<250 && sparm_mL<75");

  TCut sel1("sparm_mG==400 && sparm_mL==100");
  TCut sel2("sparm_mG==800 && sparm_mL==100");
  TCut sel3("sparm_mG==800 && sparm_mL==400");
  TCut sel=sel1||sel2||sel3;

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
