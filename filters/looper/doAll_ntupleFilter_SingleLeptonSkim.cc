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

  prefix.push_back("/nfs-7/userdata/cms2");
  samples.push_back("WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton");

  const unsigned int nsamples = samples.size();

  for( unsigned int i = 0 ; i < nsamples ; i++ ){

    cout << endl << endl;
    cout << "Doing sample    :  " << samples.at(i) << endl;
    gROOT->ProcessLine(Form(".! mkdir -p ../output/%s/%s" , samples.at(i).c_str() , skimname ));

    char* infile  = Form( "%s/%s/merged_ntuple.root"      , prefix.at(i).c_str() , samples.at(i).c_str() );
    //char* infile  = Form( "%s/%s/merged*root"      , prefix.at(i).c_str() , samples.at(i).c_str()      );
    char* outfile = Form( "../output/%s/%s/merged_ntuple.root" , samples.at(i).c_str() , skimname );

    cout << "Reading in file :  " << infile  << endl;
    cout << "Writing to file :  " << outfile << endl;
    cout << endl << endl;
    
    ntupleFilter( infile , outfile );
  }

}
