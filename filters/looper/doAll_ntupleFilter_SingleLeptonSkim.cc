{
  const char* skimname = "SingleLeptonSkim";
  
  //Load CORE stuff
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/utilities.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ntupleFilter_SingleLeptonSkim.cc++");

  vector<string> samples;
  vector<string> prefix;

  //----------------------
  // W+jets
  //----------------------

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton");

  //----------------------
  // QCD
  //----------------------

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31");

  //----------------------
  // QCD
  //----------------------

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29");

  prefix.push_back("/nfs-3/userdata/cms2");
  samples.push_back("DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29");

  prefix.push_back("/nfs-3/userdata/cms2");  
  samples.push_back("DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29");

  prefix.push_back("/nfs-3/userdata/cms2");
  samples.push_back("DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29");

  prefix.push_back("/nfs-3/userdata/cms2");
  samples.push_back("DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29");

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29");





  const unsigned int nsamples = samples.size();

  for( unsigned int i = 0 ; i < nsamples ; i++ ){

    cout << endl << endl;
    cout << "Doing sample    :  " << samples.at(i) << endl;
    gROOT->ProcessLine(Form(".! mkdir -p ../output/%s/%s" , samples.at(i).c_str() , skimname ));

    //char* infile  = Form( "%s/%s/merged_ntuple.root"      , prefix.at(i).c_str() , samples.at(i).c_str() );
    char* infile  = Form( "%s/%s/merged*root"      , prefix.at(i).c_str() , samples.at(i).c_str()      );
    char* outfile = Form( "../output/%s/%s/merged_ntuple.root" , samples.at(i).c_str() , skimname );

    cout << "Reading in file :  " << infile  << endl;
    cout << "Writing to file :  " << outfile << endl;
    cout << endl << endl;
    
    ntupleFilter( infile , outfile );
  }

}
