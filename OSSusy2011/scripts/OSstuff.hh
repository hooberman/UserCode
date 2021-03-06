
TGraphErrors* getObserved_NLOunc( ){

  Int_t n=67;
  Double_t x[n],y[n];
  Double_t xerr[n],yerr[n];
  for( unsigned int ierr = 0 ; ierr < n ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
  }

  Int_t i = -1;

  x[++i]=100;
  y[i]=350.;

  x[++i]=300;
  y[i]=340.;

  x[++i]=500;
  y[i]=320.;

  x[++i]=300;
  y[i]=340.;

  x[++i]=800;
  y[i]=320.;

  x[++i]=840;
  y[i]=100.;


  /*
  x[++i]=0;
  y[i]=255.;

  x[++i]=20;
  y[i]=255.;

  x[++i]=30.;
  y[i]=265.;

  x[++i]=40.;
  y[i]=265.;

  x[++i]=70.;
  y[i]=265.;
  
  x[++i]=80.;
  y[i]=275.;

  x[++i]=100.;
  y[i]=275.;

  x[++i]=120.;
  y[i]=275.;

  x[++i]=130.;
  y[i]=275.;

  x[++i]=140.;
  y[i]=275.;

  x[++i]=145.;
  y[i]=265.;

  x[++i]=150.;
  y[i]=255.;

  x[++i]=160.;
  y[i]=235.;

  x[++i]=170.;
  y[i]=225.;

  x[++i]=180.;
  y[i]=215.;

  x[++i]=190.;
  y[i]=208.;

  x[++i]=200.;
  y[i]=205.;

  //x[++i]=210.;
  //y[i]=195.;

  //x[++i]=220.;
  //y[i]=170.;

  x[++i]=230.;
  y[i]=185.;

  x[++i]=240.;
  y[i]=175.;

  x[++i]=250.;
  y[i]=175.;

  //x[++i]=260.;
  //y[i]=155.;

  x[++i]=270.;
  y[i]=165.;

  x[++i]=280.;
  y[i]=162.;

  x[++i]=290.;
  y[i]=165.;

  x[++i]=300.;
  y[i]=165.;

  x[++i]=310.;
  y[i]=160.;

  x[++i]=300.;
  y[i]=155.;

  x[++i]=290.;
  y[i]=155.;

  x[++i]=280.;
  y[i]=160.;

  x[++i]=275.;
  y[i]=150.;

  x[++i]=280.;
  y[i]=135.;

  x[++i]=290.;
  y[i]=125.;;

  x[++i]=300.;
  y[i]=125.;;

  x[++i]=310.;
  y[i]=120.;

  x[++i]=320.;
  y[i]=125.;

  x[++i]=330.;
  y[i]=125.;

  x[++i]=340.;
  y[i]=135.;

  x[++i]=350.;
  y[i]=135.;

  x[++i]=360.;
  y[i]=135.;

  x[++i]=370.;
  y[i]=135.;

  x[++i]=380.;
  y[i]=135.;

  x[++i]=390.;
  y[i]=135.;

  x[++i]=400.;
  y[i]=145.;

  x[++i]=420.;
  y[i]=145.;

  x[++i]=430.;
  y[i]=145.;

  x[++i]=440.;
  y[i]=145.;

  x[++i]=450.;
  y[i]=145.;

  x[++i]=460.;
  y[i]=145.;

  x[++i]=465.;
  y[i]=140.;

  x[++i]=460.;
  y[i]=135.;

  x[++i]=450.;
  y[i]=135.;

  x[++i]=440.;
  y[i]=135.;

  x[++i]=430.;
  y[i]=135.;

  x[++i]=420.;
  y[i]=135.;

  x[++i]=410.;
  y[i]=125.;

  x[++i]=400.;
  y[i]=125.;

  x[++i]=390.;
  y[i]=125.;

  x[++i]=380.;
  y[i]=125.;

  x[++i]=370.;
  y[i]=125.;

  x[++i]=360.;
  y[i]=115.;

  x[++i]=350.;
  y[i]=115.;

  x[++i]=340.;
  y[i]=115.;

  x[++i]=330.;
  y[i]=115.;

  x[++i]=320.;
  y[i]=115.;

  x[++i]=310.;
  y[i]=115.;

  x[++i]=300.;
  y[i]=100.;
  */

  TGraphErrors* gr = new TGraphErrors(n,x,y,xerr,yerr);
  gr->SetMarkerColor(kWhite);
  return gr;

}



TGraphErrors* getExpected_NLOunc( ){


  Int_t nexp=71;
  Double_t xexp[nexp],yexp[nexp];
  Double_t xerr[nexp],yerr[nexp];
  for( unsigned int ierr = 0 ; ierr < nexp ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
  }

  Int_t i = -1;

  xexp[++i]=0;
  yexp[i]=255.;

  xexp[++i]=20;
  yexp[i]=245.;

  xexp[++i]=30.;
  yexp[i]=255.;

  xexp[++i]=40.;
  yexp[i]=255.;

  xexp[++i]=50.;
  yexp[i]=265.;

  xexp[++i]=60.;
  yexp[i]=265.;

  xexp[++i]=70.;
  yexp[i]=265.;
  
  xexp[++i]=80.;
  yexp[i]=273.;

  xexp[++i]=90.;
  yexp[i]=265.;

  xexp[++i]=100.;
  yexp[i]=265.;

  xexp[++i]=110.;
  yexp[i]=265.;

  xexp[++i]=120.;
  yexp[i]=275.;

  xexp[++i]=130.;
  yexp[i]=275.;

  xexp[++i]=140.;
  yexp[i]=235.;

  xexp[++i]=150.;
  yexp[i]=235.;

  xexp[++i]=160.;
  yexp[i]=235.;

  xexp[++i]=170.;
  yexp[i]=225.;

  xexp[++i]=180.;
  yexp[i]=205.;

  xexp[++i]=190.;
  yexp[i]=205.;

  xexp[++i]=200.;
  yexp[i]=205.;

  xexp[++i]=210.;
  yexp[i]=195.;

  xexp[++i]=220.;
  yexp[i]=155.;

  xexp[++i]=230.;
  yexp[i]=155.;

  xexp[++i]=240.;
  yexp[i]=165.;

  xexp[++i]=250.;
  yexp[i]=155.;

  xexp[++i]=260.;
  yexp[i]=155.;

  xexp[++i]=270.;
  yexp[i]=155.;

  //
  xexp[++i]=280.;
  yexp[i]=160.;

  xexp[++i]=290.;
  yexp[i]=165.;

  xexp[++i]=300.;
  yexp[i]=165.;

  xexp[++i]=310.;
  yexp[i]=160.;

  xexp[++i]=300.;
  yexp[i]=155.;

  xexp[++i]=290.;
  yexp[i]=159.;

  xexp[++i]=280.;
  yexp[i]=160.;

  xexp[++i]=275.;
  yexp[i]=150.;

  xexp[++i]=255.;
  yexp[i]=135.;

  xexp[++i]=290.;
  yexp[i]=125.;;

  xexp[++i]=300.;
  yexp[i]=125.;;

  xexp[++i]=310.;
  yexp[i]=120.;

  xexp[++i]=320.;
  yexp[i]=125.;

  xexp[++i]=330.;
  yexp[i]=125.;

  xexp[++i]=340.;
  yexp[i]=125.;

  xexp[++i]=350.;
  yexp[i]=125.;

  xexp[++i]=360.;
  yexp[i]=135.;

  xexp[++i]=370.;
  yexp[i]=135.;

  xexp[++i]=380.;
  yexp[i]=135.;

  xexp[++i]=390.;
  yexp[i]=135.;

  xexp[++i]=400.;
  yexp[i]=145.;

  xexp[++i]=420.;
  yexp[i]=145.;

  xexp[++i]=430.;
  yexp[i]=145.;

  xexp[++i]=440.;
  yexp[i]=145.;

  xexp[++i]=450.;
  yexp[i]=145.;

  xexp[++i]=460.;
  yexp[i]=145.;

  xexp[++i]=460.;
  yexp[i]=140.;

  xexp[++i]=460.;
  yexp[i]=135.;

  xexp[++i]=450.;
  yexp[i]=135.;

  xexp[++i]=440.;
  yexp[i]=135.;

  xexp[++i]=430.;
  yexp[i]=135.;

  xexp[++i]=420.;
  yexp[i]=135.;

  xexp[++i]=410.;
  yexp[i]=130.;

  xexp[++i]=400.;
  yexp[i]=130.;

  xexp[++i]=390.;
  yexp[i]=125.;

  xexp[++i]=380.;
  yexp[i]=125.;

  xexp[++i]=370.;
  yexp[i]=125.;

  xexp[++i]=360.;
  yexp[i]=115.;

  xexp[++i]=350.;
  yexp[i]=115.;

  xexp[++i]=340.;
  yexp[i]=115.;

  xexp[++i]=330.;
  yexp[i]=115.;

  xexp[++i]=320.;
  yexp[i]=115.;

  xexp[++i]=310.;
  yexp[i]=115.;

  xexp[++i]=300.;
  yexp[i]=100.;

  //cout << i<< endl;

  TGraphErrors* grexp = new TGraphErrors(nexp,xexp,yexp,xerr,yerr);
  grexp->SetMarkerColor(kWhite);
  return grexp;

}

TGraphErrors* getLO( ){
  
  Int_t nLO=39;
  Double_t xLO[nLO],yLO[nLO];
  Double_t xerr[nLO],yerr[nLO];
  for( unsigned int ierr = 0 ; ierr < 39 ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xLO[ierr] = 0.;
    yLO[ierr] = 0.;
  }

  Int_t i = -1;

  xLO[++i]=0;
  yLO[i]=245.;

  xLO[++i]=20;
  yLO[i]=245.;

  xLO[++i]=30.;
  yLO[i]=245.;

  xLO[++i]=40.;
  yLO[i]=255.;

  xLO[++i]=50.;
  yLO[i]=255.;

  xLO[++i]=60.;
  yLO[i]=255.;
  
  xLO[++i]=70.;
  yLO[i]=263.;
  
  xLO[++i]=80.;
  yLO[i]=265.;

  xLO[++i]=100.;
  yLO[i]=265.;

  xLO[++i]=120.;
  yLO[i]=265.;

  xLO[++i]=130.;
  yLO[i]=245.;

  xLO[++i]=140.;
  yLO[i]=235.;

  xLO[++i]=150.;
  yLO[i]=225.;

  xLO[++i]=160.;
  yLO[i]=215.;

  xLO[++i]=170.;
  yLO[i]=205.;

  xLO[++i]=180.;
  yLO[i]=195.;

  xLO[++i]=190.;
  yLO[i]=195.;

  xLO[++i]=200.;
  yLO[i]=180.;

  xLO[++i]=210.;
  yLO[i]=175.;

  xLO[++i]=220.;
  yLO[i]=155.;

  xLO[++i]=230.;
  yLO[i]=150.;

  xLO[++i]=240.;
  yLO[i]=145.;

  xLO[++i]=250.;
  yLO[i]=155.;

  xLO[++i]=255.;
  yLO[i]=140.;

  xLO[++i]=255.;
  yLO[i]=135.;

  xLO[++i]=240.;
  yLO[i]=125.;

  xLO[++i]=230.;
  yLO[i]=115.;

  xLO[++i]=240.;
  yLO[i]=115.;

  xLO[++i]=250.;
  yLO[i]=115.;

  xLO[++i]=260.;
  yLO[i]=115.;

  xLO[++i]=270.;
  yLO[i]=125.;

  xLO[++i]=280.;
  yLO[i]=125.;

  xLO[++i]=290.;
  yLO[i]=125.;;

  xLO[++i]=300.;
  yLO[i]=125.;;

  xLO[++i]=310.;
  yLO[i]=120.;

  xLO[++i]=320.;
  yLO[i]=125.;

  xLO[++i]=320.;
  yLO[i]=120.;

  xLO[++i]=310.;
  yLO[i]=110.;

  xLO[++i]=305.;
  yLO[i]=100.;
  

  //cout << i<< endl;
  TGraphErrors* grLO  = new TGraphErrors(nLO,xLO, yLO,xerr,yerr);
  grLO->SetMarkerColor(kGreen+2);
  grLO->SetMarkerStyle(21);
  return grLO;



  
//     Double_t xxLO[9],yyLO[9];
//     xxLO[0]=345.;
//     yyLO[0]=120.;

//     xxLO[1]=350.;
//     yyLO[1]=125.;

//     xxLO[2]=360.;
//     yyLO[2]=135.;

//     xxLO[3]=370.;
//     yyLO[3]=135.;

//     xxLO[4]=375.;
//     yyLO[4]=130.;

//     xxLO[5]=370.;
//     yyLO[5]=125.;

//     xxLO[6]=360.;
//     yyLO[6]=115.;

//     xxLO[7]=350.;
//     yyLO[7]=115.;

//     xxLO[8]=345.;
//     yyLO[8]=120.;
 
//     grLO2 = new TGraph(9,  xxLO,yyLO);

//     // Draw the curves
//     grLO.Draw("C");
//     grLO2.Draw("C");

}


TGraphErrors* getLO_signalCont( ){

    TGraphErrors *gre = getLO();
    return gre;
  
}


TGraphErrors* getNLOobsTanbeta10( ){
  
  const unsigned int nNLO=50;

  Double_t xNLO[nNLO],yNLO[nNLO],xerr[nNLO],yerr[nNLO];
  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Int_t i = -1;

  //------------------
  // new scan
  //------------------

  xNLO[++i]=0.;
  yNLO[i]=460.;

  xNLO[++i]=80;
  yNLO[i]=410.;

  xNLO[++i]=160;
  yNLO[i]=330.;

  xNLO[++i]=200;
  yNLO[i]=300.;

  xNLO[++i]=350;
  yNLO[i]=310.;

  xNLO[++i]=450;
  yNLO[i]=300.;

  // xNLO[++i]=520;
  // yNLO[i]=300.;

  // xNLO[++i]=560;
  // yNLO[i]=280.;

  xNLO[++i]=600;
  yNLO[i]=270.;

  xNLO[++i]=700;
  yNLO[i]=250.;

  xNLO[++i]=850;
  yNLO[i]=230.;

  xNLO[++i]=900;
  yNLO[i]=220.;

  xNLO[++i]=1000;
  yNLO[i]=200.;

  xNLO[++i]=1200;
  yNLO[i]=140.;

  xNLO[++i]=1400;
  yNLO[i]=100.;



  //------------------
  // new scan
  //------------------

  /*
  xNLO[++i]=50;
  yNLO[i]=500.;

  xNLO[++i]=80;
  yNLO[i]=420.;

  xNLO[++i]=120;
  yNLO[i]=370.;

  xNLO[++i]=150;
  yNLO[i]=350.;

  xNLO[++i]=300;
  yNLO[i]=330.;

  xNLO[++i]=500;
  yNLO[i]=310.;

  xNLO[++i]=700;
  yNLO[i]=280.;

  xNLO[++i]=800;
  yNLO[i]=250.;

  xNLO[++i]=850;
  yNLO[i]=210.;

  xNLO[++i]=900;
  yNLO[i]=100.;

  xNLO[++i]=1000;
  yNLO[i]=0.;
  */

  /*
  xNLO[++i]=0;
  yNLO[i]=410.;

  xNLO[++i]=80;
  yNLO[i]=410.;

  xNLO[++i]=120;
  yNLO[i]=370.;

  xNLO[++i]=140;
  yNLO[i]=350.;

  xNLO[++i]=160;
  yNLO[i]=350.;

  xNLO[++i]=180;
  yNLO[i]=350.;

  xNLO[++i]=200;
  yNLO[i]=350.;

  xNLO[++i]=220;
  yNLO[i]=350.;

  xNLO[++i]=240;
  yNLO[i]=350.;

  xNLO[++i]=260;
  yNLO[i]=350.;

  xNLO[++i]=280;
  yNLO[i]=350.;

  xNLO[++i]=300;
  yNLO[i]=330.;

  xNLO[++i]=350;
  yNLO[i]=310.;

  xNLO[++i]=370;
  yNLO[i]=310.;

  xNLO[++i]=390;
  yNLO[i]=310.;

  xNLO[++i]=410;
  yNLO[i]=310.;

  xNLO[++i]=430;
  yNLO[i]=310.;

  xNLO[++i]=450;
  yNLO[i]=310.;

  xNLO[++i]=470;
  yNLO[i]=310.;

  xNLO[++i]=500;
  yNLO[i]=310.;

  xNLO[++i]=550;
  yNLO[i]=270.;

  xNLO[++i]=570;
  yNLO[i]=270.;

  xNLO[++i]=590;
  yNLO[i]=270.;

  xNLO[++i]=610;
  yNLO[i]=270.;

  xNLO[++i]=630;
  yNLO[i]=270.;

  xNLO[++i]=650;
  yNLO[i]=270.;

  xNLO[++i]=670;
  yNLO[i]=270.;

  xNLO[++i]=690;
  yNLO[i]=270.;

  xNLO[++i]=710;
  yNLO[i]=270.;

  xNLO[++i]=730;
  yNLO[i]=270.;

  xNLO[++i]=750;
  yNLO[i]=270.;

  xNLO[++i]=800;
  yNLO[i]=250.;

  xNLO[++i]=830;
  yNLO[i]=230.;

  xNLO[++i]=850;
  yNLO[i]=210.;

  xNLO[++i]=870;
  yNLO[i]=130.;

  xNLO[++i]=890;
  yNLO[i]=110.;

  xNLO[++i]=910;
  yNLO[i]=110.;

  xNLO[++i]=920;
  yNLO[i]=110.;

  xNLO[++i]=922;
  yNLO[i]=110.;

  xNLO[++i]=924;
  yNLO[i]=110.;

  xNLO[++i]=926;
  yNLO[i]=110.;

  xNLO[++i]=928;
  yNLO[i]=110.;

  xNLO[++i]=930;
  yNLO[i]=110.;

  xNLO[++i]=950;
  yNLO[i]=0.;

  */

  TGraphErrors* grtb10  = new TGraphErrors(nNLO,xNLO, yNLO,xerr,yerr);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  return grtb10;

}


TGraphErrors* getNLOexpTanbeta10( ){
  
  const unsigned int nNLO=50;

  Double_t xNLO[nNLO],yNLO[nNLO],xerr[nNLO],yerr[nNLO];
  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Int_t i = -1;


  //------------------
  // take obs, subtract 15
  //------------------


  xNLO[++i]=0.;
  yNLO[i]=460.;

  xNLO[++i]=80;
  yNLO[i]=410.;

  xNLO[++i]=160;
  yNLO[i]=330.;

  xNLO[++i]=200;
  yNLO[i]=300.;

  xNLO[++i]=350;
  yNLO[i]=310.;

  xNLO[++i]=450;
  yNLO[i]=300.;

  // xNLO[++i]=520;
  // yNLO[i]=300.;

  // xNLO[++i]=560;
  // yNLO[i]=280.;

  xNLO[++i]=600;
  yNLO[i]=270.;

  xNLO[++i]=700;
  yNLO[i]=250.;

  xNLO[++i]=850;
  yNLO[i]=230.;

  xNLO[++i]=900;
  yNLO[i]=220.;

  xNLO[++i]=1000;
  yNLO[i]=200.;

  xNLO[++i]=1200;
  yNLO[i]=140.;

  xNLO[++i]=1400;
  yNLO[i]=100.;



  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    if( yNLO[ierr] > 0 ) yNLO[ierr] -= 15.0;

  }

    /*
  //------------
  // old scan
  //------------

  xNLO[++i]=0;
  yNLO[i]=440.;

  xNLO[++i]=80;
  yNLO[i]=390.;

  xNLO[++i]=120;
  yNLO[i]=350.;

  xNLO[++i]=170;
  yNLO[i]=330.;

  xNLO[++i]=180;
  yNLO[i]=300.;

  xNLO[++i]=200;
  yNLO[i]=280.;

  xNLO[++i]=220;
  yNLO[i]=280.;

  xNLO[++i]=300;
  yNLO[i]=280.;

  xNLO[++i]=400;
  yNLO[i]=300.;

  xNLO[++i]=500;
  yNLO[i]=300.;

  xNLO[++i]=550;
  yNLO[i]=270.;

  xNLO[++i]=600;
  yNLO[i]=260.;

  xNLO[++i]=700;
  yNLO[i]=240.;

  xNLO[++i]=800;
  yNLO[i]=240.;

  xNLO[++i]=850;
  yNLO[i]=230.;

  xNLO[++i]=900;
  yNLO[i]=220.;

  xNLO[++i]=1000;
  yNLO[i]=200.;

  xNLO[++i]=1200;
  yNLO[i]=160.;

  xNLO[++i]=1400;
  yNLO[i]=80.;
    */

  //------------
  // new scan
  //------------
  /*
  xNLO[++i]=50;
  yNLO[i]=450.;

  xNLO[++i]=80;
  yNLO[i]=390.;

  xNLO[++i]=120;
  yNLO[i]=350.;

  xNLO[++i]=180;
  yNLO[i]=320.;

  xNLO[++i]=300;
  yNLO[i]=300.;

  xNLO[++i]=500;
  yNLO[i]=290.;

  xNLO[++i]=600;
  yNLO[i]=270.;

  xNLO[++i]=700;
  yNLO[i]=250.;

  xNLO[++i]=800;
  yNLO[i]=220.;

  xNLO[++i]=850;
  yNLO[i]=190.;

  xNLO[++i]=900;
  yNLO[i]=100.;

  xNLO[++i]=1000;
  yNLO[i]=0.;
  */
  /*
  xNLO[++i]=0;
  yNLO[i]=410.;

  xNLO[++i]=80;
  yNLO[i]=410.;

  xNLO[++i]=120;
  yNLO[i]=370.;

  xNLO[++i]=140;
  yNLO[i]=350.;

  xNLO[++i]=160;
  yNLO[i]=350.;

  xNLO[++i]=180;
  yNLO[i]=350.;

  xNLO[++i]=200;
  yNLO[i]=350.;

  xNLO[++i]=220;
  yNLO[i]=350.;

  xNLO[++i]=240;
  yNLO[i]=350.;

  xNLO[++i]=260;
  yNLO[i]=350.;

  xNLO[++i]=280;
  yNLO[i]=350.;

  xNLO[++i]=300;
  yNLO[i]=330.;

  xNLO[++i]=350;
  yNLO[i]=310.;

  xNLO[++i]=370;
  yNLO[i]=310.;

  xNLO[++i]=390;
  yNLO[i]=310.;

  xNLO[++i]=410;
  yNLO[i]=310.;

  xNLO[++i]=430;
  yNLO[i]=310.;

  xNLO[++i]=450;
  yNLO[i]=310.;

  xNLO[++i]=470;
  yNLO[i]=310.;

  xNLO[++i]=500;
  yNLO[i]=310.;

  xNLO[++i]=550;
  yNLO[i]=270.;

  xNLO[++i]=570;
  yNLO[i]=270.;

  xNLO[++i]=590;
  yNLO[i]=270.;

  xNLO[++i]=610;
  yNLO[i]=270.;

  xNLO[++i]=630;
  yNLO[i]=270.;

  xNLO[++i]=650;
  yNLO[i]=270.;

  xNLO[++i]=670;
  yNLO[i]=270.;

  xNLO[++i]=690;
  yNLO[i]=270.;

  xNLO[++i]=710;
  yNLO[i]=270.;

  xNLO[++i]=730;
  yNLO[i]=270.;

  xNLO[++i]=750;
  yNLO[i]=270.;

  xNLO[++i]=800;
  yNLO[i]=250.;

  xNLO[++i]=830;
  yNLO[i]=230.;

  xNLO[++i]=850;
  yNLO[i]=210.;

  xNLO[++i]=870;
  yNLO[i]=130.;

  xNLO[++i]=890;
  yNLO[i]=110.;

  xNLO[++i]=910;
  yNLO[i]=110.;

  xNLO[++i]=920;
  yNLO[i]=110.;

  xNLO[++i]=922;
  yNLO[i]=110.;

  xNLO[++i]=924;
  yNLO[i]=110.;

  xNLO[++i]=926;
  yNLO[i]=110.;

  xNLO[++i]=928;
  yNLO[i]=110.;

  xNLO[++i]=930;
  yNLO[i]=110.;

  xNLO[++i]=950;
  yNLO[i]=0.;

  */

  TGraphErrors* grtb10  = new TGraphErrors(nNLO,xNLO, yNLO,xerr,yerr);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  return grtb10;

}



