{
  const char* skimname = "met50skim";
  
  //Load CORE stuff
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  // gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  // gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  // gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/SimpleFakeRate.cc+");
  // gROOT->ProcessLine(".L ../CORE/mcSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/MT2/MT2.cc+");
  // gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
  // gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");
  gROOT->ProcessLine(".L ntupleFilter.cc+");

  vector<char*> samples;
  vector<char*> prefix;

  prefix.push_back("/nfs-3/userdata/cms2"); samples.push_back("DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/");
  prefix.push_back("/nfs-3/userdata/cms2"); samples.push_back("DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/");
  prefix.push_back("/nfs-3/userdata/cms2"); samples.push_back("DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/");
  prefix.push_back("/nfs-3/userdata/cms2"); samples.push_back("DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/");
  prefix.push_back("/nfs-7/userdata/cms2"); samples.push_back("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/");
  prefix.push_back("/nfs-7/userdata/cms2"); samples.push_back("DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/");
  prefix.push_back("/nfs-7/userdata/cms2"); samples.push_back("DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/");

  const unsigned int nsamples = samples.size();

  for( unsigned int i = 0 ; i < nsamples ; i++ ){

    cout << endl << endl;
    cout << "Doing sample    :  " << samples.at(i) << endl;
    gROOT->ProcessLine(Form(".! mkdir -p output/%s/%s" , samples.at(i) , skimname ));

    //char* infile  = Form( "%s/%s/merged_ntuple.root"      , prefix.at(i) , samples.at(i)            );
    char* infile  = Form( "%s/%s/merged*root"      , prefix.at(i) , samples.at(i)            );
    char* outfile = Form( "output/%s/%s/merged_ntuple.root" , samples.at(i) , skimname );

    cout << "Reading in file :  " << infile  << endl;
    cout << "Writing to file :  " << outfile << endl;
    cout << endl << endl;

    ntupleFilter( infile , outfile );
  }






}



  // Load CORE stuff
  // gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");  
  // gROOT->ProcessLine(".L ../CORE/susySelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/mcSUSYkfactor.cc+");
  // gROOT->ProcessLine(".L ../CORE/triggerSuperModel.cc+");
  // gROOT->ProcessLine(".L ../CORE/jetSelections.cc+");
  // gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");


      /*
  //TAS
  //samples.push_back("DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //TAS
  //samples.push_back("DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //UAF
  samples.push_back("DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //UAF
  samples.push_back("DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //TAS
  //samples.push_back("DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //UAF
  samples.push_back("DYToTauTau_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v2/V04-01-01");

  //TAS
  //samples.push_back("DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");
  */
