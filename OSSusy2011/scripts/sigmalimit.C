{


  bool doLimits = true;

  //-------------------------
  // xsec
  //-------------------------

  float   LM1_sigma_2010    = 4.89 * 1.42;
  float   LM1_sigma_highmet = 4.89 * 1.42;
  float   LM1_sigma_highht  = 4.89 * 1.47;

  float   LM3_sigma_2010    = 3.44 * 1.51;
  float   LM3_sigma_highmet = 3.44 * 1.49;
  float   LM3_sigma_highht  = 3.44 * 1.56;

  float   LM6_sigma_2010    = 0.31 * 1.51;
  float   LM6_sigma_highmet = 0.31 * 1.52;
  float   LM6_sigma_highht  = 0.31 * 1.56;

  //-------------------------
  // acceptances
  //-------------------------

  float   LM1_acc_2010    = 0.038 * (114514./219190.);
  float   LM1_acc_highmet = 0.024 * (114514./219190.);
  float   LM1_acc_highht  = 0.018 * (114514./219190.);

  float   LM3_acc_2010    = 0.021 * (123627./220000.);
  float   LM3_acc_highmet = 0.012 * (123627./220000.);
  float   LM3_acc_highht  = 0.012 * (123627./220000.);

  float   LM6_acc_2010    = 0.063 * (121523./220000.);
  float   LM6_acc_highmet = 0.051 * (121523./220000.);
  float   LM6_acc_highht  = 0.045 * (121523./220000.);


  //-------------------------
  // yields
  //-------------------------

  int   nobs_2010      = 45;
  int   nobs_highmet   = 8;
  int   nobs_highht    = 4;

  //-------------------------
  // backgrounds
  //-------------------------

  float bkg_2010       = 42.;
  float bkgerr_2010    = 8.5.;
  float bkg_highmet    = 4.2;
  float bkgerr_highmet = 1.3.;
  float bkg_highht     = 5.1.;
  float bkgerr_highht  = 1.7.;

  //-------------------------
  // efficiencies
  //-------------------------

  float LM1_2010    = 0.50;
  float LM1_highmet = 0.46;
  float LM1_highht  = 0.42;

  float LM3_2010    = 0.46;
  float LM3_highmet = 0.43;
  float LM3_highht  = 0.40;

  float LM6_2010    = 0.52;
  float LM6_highmet = 0.54;
  float LM6_highht  = 0.52;

  //-------------------------
  // JES uncertainties
  //-------------------------

  float LM1_JES_2010    = 0.09;
  float LM1_JES_highmet = 0.22;
  float LM1_JES_highht  = 0.30;

  float LM3_JES_2010    = 0.10;
  float LM3_JES_highmet = 0.27;
  float LM3_JES_highht  = 0.32;

  float LM6_JES_2010    = 0.06;
  float LM6_JES_highmet = 0.10;
  float LM6_JES_highht  = 0.14;

  //-------------------------
  // tot uncertainties
  //-------------------------

  float systerr2 = 0.02*0.02 + 0.02*0.02 + 0.02*0.02;

  float LM1_2010_err    = LM1_2010    * sqrt( LM1_JES_2010    * LM1_JES_2010    + systerr2 );
  float LM1_highmet_err = LM1_highmet * sqrt( LM1_JES_highmet * LM1_JES_highmet + systerr2 );
  float LM1_highht_err  = LM1_highht  * sqrt( LM1_JES_highht  * LM1_JES_highht  + systerr2 );
  
  float LM3_2010_err    = LM3_2010    * sqrt( LM3_JES_2010    * LM3_JES_2010    + systerr2 );
  float LM3_highmet_err = LM3_highmet * sqrt( LM3_JES_highmet * LM3_JES_highmet + systerr2 );
  float LM3_highht_err  = LM3_highht  * sqrt( LM3_JES_highht  * LM3_JES_highht  + systerr2 );
  
  float LM6_2010_err    = LM6_2010    * sqrt( LM6_JES_2010    * LM6_JES_2010    + systerr2 );
  float LM6_highmet_err = LM6_highmet * sqrt( LM6_JES_highmet * LM6_JES_highmet + systerr2 );
  float LM6_highht_err  = LM6_highht  * sqrt( LM6_JES_highht  * LM6_JES_highht  + systerr2 );
  
  
  char* pm = "$\\pm$";

  cout << Form("LM1  &  %.2f %s %.2f  &    %.2f %s %.2f  &    %.2f %s %.2f \\\\",
	       LM1_2010,pm,LM1_2010_err,LM1_highmet,pm,LM1_highmet_err,LM1_highht,pm,LM1_highht_err) << endl;
  cout << Form("LM3  &  %.2f %s %.2f  &    %.2f %s %.2f  &    %.2f %s %.2f \\\\",
	       LM3_2010,pm,LM3_2010_err,LM3_highmet,pm,LM3_highmet_err,LM3_highht,pm,LM3_highht_err) << endl;
  cout << Form("LM6  &  %.2f %s %.2f  &    %.2f %s %.2f  &    %.2f %s %.2f \\\\",
	       LM6_2010,pm,LM6_2010_err,LM6_highmet,pm,LM6_highmet_err,LM6_highht,pm,LM6_highht_err) << endl;
  


  if( doLimits ){

    gROOT->ProcessLine(".L ~/ntupling/CMS2/NtupleMacros/Statistics/cl95cms.c+");

    //sigma95 = CL95(ilum, slum, eff, seff, bck, sbck, n, gauss = false, nuisanceModel = 0)

    cout << "LM1 2010" << endl;
    float LM1_2010_UL = CL95( 976 , 59 , LM1_2010 , LM1_2010_err , bkg_2010 , bkgerr_2010 , nobs_2010 , kFALSE ,1 );
    cout << "LM1 highmet" << endl;
    float LM1_highmet_UL = CL95( 976 , 59 , LM1_highmet , LM1_highmet_err , bkg_highmet , bkgerr_highmet , nobs_highmet , kFALSE ,1 );
    cout << "LM1 highht" << endl;
    float LM1_highht_UL = CL95( 976 , 59 , LM1_highht , LM1_highht_err , bkg_highht , bkgerr_highht , nobs_highht , kFALSE ,1 );

    cout << endl;

    cout << "LM3 2010" << endl;
    float LM3_2010_UL = CL95( 976 , 59 , LM3_2010 , LM3_2010_err , bkg_2010 , bkgerr_2010 , nobs_2010 , kFALSE ,1 );
    cout << "LM3 highmet" << endl;
    float LM3_highmet_UL = CL95( 976 , 59 , LM3_highmet , LM3_highmet_err , bkg_highmet , bkgerr_highmet , nobs_highmet , kFALSE ,1 );
    cout << "LM3 highht" << endl;
    float LM3_highht_UL = CL95( 976 , 59 , LM3_highht , LM3_highht_err , bkg_highht , bkgerr_highht , nobs_highht , kFALSE ,1 );

    cout << endl;

    cout << "LM6 2010" << endl;
    float LM6_2010_UL = CL95( 976 , 59 , LM6_2010 , LM6_2010_err , bkg_2010 , bkgerr_2010 , nobs_2010 , kFALSE ,1 );
    cout << "LM6 highmet" << endl;
    float LM6_highmet_UL = CL95( 976 , 59 , LM6_highmet , LM6_highmet_err , bkg_highmet , bkgerr_highmet , nobs_highmet , kFALSE ,1 );
    cout << "LM6 highht" << endl;
    float LM6_highht_UL = CL95( 976 , 59 , LM6_highht , LM6_highht_err , bkg_highht , bkgerr_highht , nobs_highht , kFALSE ,1 );

    cout << endl;
    cout << "LM1 2010    " << Form("%.3f",LM1_2010_UL)    << endl;
    cout << "LM1 highmet " << Form("%.3f",LM1_highmet_UL) << endl;
    cout << "LM1 highht  " << Form("%.3f",LM1_highht_UL)  << endl;
    cout << endl;
    cout << "LM3 2010    " << Form("%.3f",LM3_2010_UL)    << endl;
    cout << "LM3 highmet " << Form("%.3f",LM3_highmet_UL) << endl;
    cout << "LM3 highht  " << Form("%.3f",LM3_highht_UL)  << endl;
    cout << endl;
    cout << "LM6 2010    " << Form("%.3f",LM6_2010_UL)    << endl;
    cout << "LM6 highmet " << Form("%.3f",LM6_highmet_UL) << endl;
    cout << "LM6 highht  " << Form("%.3f",LM6_highht_UL)  << endl;
    cout << endl;
  }


    cout << endl;
    cout << "LM1 2010    sigma X acc   " << Form("%.3f",LM1_sigma_2010    * LM1_acc_2010    )    << endl;
    cout << "LM1 highmet sigma X acc   " << Form("%.3f",LM1_sigma_highmet * LM1_acc_highmet )    << endl;
    cout << "LM1 highht  sigma X acc   " << Form("%.3f",LM1_sigma_highht  * LM1_acc_highht  )    << endl;
    cout << endl;
    cout << "LM3 2010    sigma X acc   " << Form("%.3f",LM3_sigma_2010    * LM3_acc_2010    )    << endl;
    cout << "LM3 highmet sigma X acc   " << Form("%.3f",LM3_sigma_highmet * LM3_acc_highmet )    << endl;
    cout << "LM3 highht  sigma X acc   " << Form("%.3f",LM3_sigma_highht  * LM3_acc_highht  )    << endl;
    cout << endl;
    cout << "LM6 2010    sigma X acc   " << Form("%.3f",LM6_sigma_2010    * LM6_acc_2010    )    << endl;
    cout << "LM6 highmet sigma X acc   " << Form("%.3f",LM6_sigma_highmet * LM6_acc_highmet )    << endl;
    cout << "LM6 highht  sigma X acc   " << Form("%.3f",LM6_sigma_highht  * LM6_acc_highht  )    << endl;


}
