{

  gROOT->ProcessLine(".L trainMVA_smurf.C+");
  gROOT->ProcessLine(".L evaluateMVA_smurf.C+");

  gROOT->ProcessLine(".! rm SmurfBabies/tas-2020/*root");
  gROOT->ProcessLine(".! cp data/*root SmurfBabies/tas-2020/.");
  
  //----------------------------------------
  // choose which MVA types to add to babies
  //----------------------------------------

  //char* mva = "Fisher";
  char* mva = "BDT,MLPBNN";

  //-------------------------------------------------
  // choose which Higgs mass MVA's to add to babies
  //-------------------------------------------------
  
  vector<int> mH;
  mH.push_back(130);
  mH.push_back(160);
  mH.push_back(200);
  mH.push_back(250);
  
  //------------------------------------------
  // run!
  //------------------------------------------

  for( unsigned int i = 0 ; i < mH.size() ; ++i ){

    //train MVA
    trainMVA_smurf(mH.at(i),mva);

    //copy MVA training files and root output
    gROOT->ProcessLine( Form( ".! cp -r weights/* SmurfTraining/hww%i_ww/weights/." , mH.at(i) ) );
    gROOT->ProcessLine( Form( ".! cp trainMVA_smurf.root SmurfTraining/hww%i_ww/." , mH.at(i) ) );

    //evaluate MVA
    evaluateMVA_smurf(mH.at(i),mva);
    gROOT->ProcessLine( Form( ".! cp  SmurfTraining/hww%i_ww/output/*.root SmurfBabies/tas-2020/." , mH.at(i) ) );

  }


}
