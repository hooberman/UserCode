{
  const char* skimname = "met50skim";
  
  //Load CORE stuff
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L ntupleFilter.cc+");

  vector<char*> samples;

  //DONE
  //samples.push_back("DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //TAS
  samples.push_back("DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //UAF
  samples.push_back("DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //UAF
  samples.push_back("DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //TAS
  samples.push_back("DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  //UAF
  samples.push_back("DYToTauTau_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v2/V04-01-01");

  //TAS
  samples.push_back("DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01");

  const unsigned int nsamples = samples.size();

  for( unsigned int i = 0 ; i < nsamples ; i++ ){
    cout << "Doing sample " << samples.at(i) << endl;
    gROOT->ProcessLine(Form(".! mkdir -p output/%s/%s" , samples.at(i) , skimname ));


    char* infile  = Form( "cms2/%s/merged*root"      , samples.at(i)            );
    char* outfile = Form( "output/%s/%s/merged_ntuple.root" , samples.at(i) , skimname );
    
    cout << "Reading in file " << infile  << endl;
    cout << "Writing to file " << outfile << endl;

    ntupleFilter( infile , outfile );
  }






}
