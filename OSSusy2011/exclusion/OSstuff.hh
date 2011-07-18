const unsigned int nNLO=20;
Double_t xNLO[100],yNLO[100],xerr[100],yerr[100];
  
void initialize(){

  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    //xerr[ierr] = 0.;
    //yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Int_t i = -1;

  xNLO[++i]=0.;    yNLO[i]=460.; //0
  xNLO[++i]=80;    yNLO[i]=410.; //1
  xNLO[++i]=100;   yNLO[i]=365.; //2
  xNLO[++i]=120;   yNLO[i]=360.; //3
  xNLO[++i]=150;   yNLO[i]=355.; //4
  xNLO[++i]=180;   yNLO[i]=325.; //5
  xNLO[++i]=200;   yNLO[i]=290.; //6
  xNLO[++i]=210;   yNLO[i]=285.; //7
  xNLO[++i]=270;   yNLO[i]=300.; //8
  xNLO[++i]=350;   yNLO[i]=310.; //9
  xNLO[++i]=450;   yNLO[i]=310.; //10
  xNLO[++i]=520;   yNLO[i]=320.; //11
  xNLO[++i]=550;   yNLO[i]=280.; //12
  xNLO[++i]=600;   yNLO[i]=265.; //13
  xNLO[++i]=700;   yNLO[i]=250.; //14
  xNLO[++i]=900;   yNLO[i]=230.; //15
  xNLO[++i]=1000;  yNLO[i]=220.; //16
  xNLO[++i]=1400;  yNLO[i]=180.; //17
  xNLO[++i]=1800;  yNLO[i]=156.; //18
  xNLO[++i]=1900;  yNLO[i]=150.; //19
}

TGraphErrors* getNLOobsTanbeta10( ){

  initialize();

  TGraphErrors* gr  = new TGraphErrors(nNLO,xNLO, yNLO,xerr,yerr);
  gr->SetMarkerColor(kWhite);
  return gr;

}

TGraphErrors* getNLOexpTanbeta10( ){
  
  initialize();

  for( int j = 0 ; j <= 100 ; ++j ){
    xNLO[j] -= 10;
    yNLO[j] -= 10;
  }

  TGraphErrors* gr  = new TGraphErrors(nNLO,xNLO, yNLO,xerr,yerr);
  gr->SetMarkerColor(kWhite);
  return gr;

}

TGraphErrors* getNLOexpUpTanbeta10( ){
  
  initialize();

  for( int j = 0 ; j <= 100 ; ++j ){
    xNLO[j] += 10;
    yNLO[j] += 10;
  }

  TGraphErrors* gr  = new TGraphErrors(nNLO,xNLO, yNLO,xerr,yerr);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(21);

  return gr;

}

TGraphErrors* getNLOexpTanbeta10_band( ){
  
  initialize();
  
  Double_t x[100],y[100],ex[100],ey[100];
  
  for( unsigned int j = 0 ; j < nNLO ; ++j ){
    x[j]  = xNLO[j] + 15;
    y[j]  = yNLO[j] + 15;
    ex[j] = 0;
    ey[j] = 0;
  }
  
  for( unsigned int j = 0 ; j < nNLO ; ++j ){
    x[j+nNLO]  = xNLO[nNLO-j-1] - 30;
    y[j+nNLO]  = yNLO[nNLO-j-1] - 30;
    ex[j+nNLO] = 0;
    ey[j+nNLO] = 0;
  }
  
  TGraphErrors* gr  = new TGraphErrors(2*nNLO,x, y,ex,ey);
  gr->SetMarkerColor(4);
  gr->SetFillColor(4);
  gr->SetMarkerStyle(21);
  gr->SetFillStyle(3002);
  gr->SetLineColor(4);
  gr->SetLineWidth(2);

  return gr;

}

TGraphErrors* getNLOexpDownTanbeta10( ){
  
  initialize();

  for( int j = 0 ; j <= 100 ; ++j ){
    xNLO[j] -= 30;
    yNLO[j] -= 30;
  }


  TGraphErrors* gr  = new TGraphErrors(nNLO,xNLO, yNLO,xerr,yerr);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(21);

  return gr;

}


TGraphErrors* getNLOobsTanbeta10_2010( ){
  
  Int_t n=47;
  Double_t x[n],y[n],xerr2[n],yerr2[n];
  for( int ierr = 0 ; ierr < n ; ++ierr ){
    xerr2[ierr] = 0.;
    yerr2[ierr] = 0.;
    x[ierr] = 0.;
    y[ierr] = 0.;
  }

  Int_t i = -1;

  x[++i]=0;
  y[i]=225.;

  x[++i]=20;
  y[i]=225.;

  x[++i]=30.;
  y[i]=235.;

  x[++i]=40.;
  y[i]=235.;

  x[++i]=50.;
  y[i]=235.;

  //x[++i]=60.;
  //y[i]=225.;

  x[++i]=70.;
  y[i]=235.;
  
  x[++i]=80.;
  y[i]=245.;

  x[++i]=90.;
  y[i]=245.;

  x[++i]=100.;
  y[i]=245.;

  x[++i]=115.;
  y[i]=240.;

  x[++i]=110.;
  y[i]=235.;

  //x[++i]=105.;
  //y[i]=230.;

  x[++i]=105.;
  y[i]=213.;

  x[++i]=105.;
  y[i]=205.;

  x[++i]=95.;
  y[i]=200.;

  x[++i]=90.;
  y[i]=185.;

  x[++i]=70.;
  y[i]=155.;

  x[++i]=65.;
  y[i]=145.;

  x[++i]=75.;
  y[i]=155.;

  x[++i]=90.;
  y[i]=170.;

  x[++i]=100.;
  y[i]=185.;

  x[++i]=110.;
  y[i]=195.;

  x[++i]=120.;
  y[i]=205.;

  x[++i]=130.;
  y[i]=215.;

  x[++i]=140.;
  y[i]=215.;

  x[++i]=150.;
  y[i]=215.;

  x[++i]=160.;
  y[i]=185.;

  x[++i]=170.;
  y[i]=175.;

  x[++i]=170.;
  y[i]=155.;

  x[++i]=160.;
  y[i]=140.;

  x[++i]=170.;
  y[i]=135.;

  x[++i]=180.;
  y[i]=125.;

  x[++i]=190.;
  y[i]=135.;
  
  x[++i]=200.;
  y[i]=145.;

  x[++i]=210.;
  y[i]=155.;

  x[++i]=220.;
  y[i]=165.;

  x[++i]=230.;
  y[i]=160.;

  x[++i]=240.;
  y[i]=145.;

  x[++i]=250.;
  y[i]=155.;
  
  x[++i]=260.;
  y[i]=170.;

  x[++i]=270.;
  y[i]=170.;

  x[++i]=280.;
  y[i]=170.;

  x[++i]=280.;
  y[i]=160.;

  x[++i]=270.;
  y[i]=150.;

  x[++i]=260.;
  y[i]=145.;

  x[++i]=230.;
  y[i]=125.;
  


  TGraphErrors* grtb10  = new TGraphErrors(n,x, y,xerr2,yerr2);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  grtb10->SetLineColor(6);
  grtb10->SetLineWidth(3);
  return grtb10;

}



