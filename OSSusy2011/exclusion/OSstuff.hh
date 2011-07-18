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

  xNLO[++i]=0.;
  yNLO[i]=460.;

  xNLO[++i]=80;
  yNLO[i]=410.;

  xNLO[++i]=100;
  yNLO[i]=365.;

  xNLO[++i]=120;
  yNLO[i]=360.;

  xNLO[++i]=150;
  yNLO[i]=355.;

  xNLO[++i]=180;
  yNLO[i]=325.;

  xNLO[++i]=200;
  yNLO[i]=290.;

  xNLO[++i]=210;
  yNLO[i]=285.;

  xNLO[++i]=270;
  yNLO[i]=300.;

  xNLO[++i]=350;
  yNLO[i]=310.;

  xNLO[++i]=450;
  yNLO[i]=310.;

  xNLO[++i]=520;
  yNLO[i]=320.;

  xNLO[++i]=550;
  yNLO[i]=280.;

  xNLO[++i]=600;
  yNLO[i]=265.;

  xNLO[++i]=700;
  yNLO[i]=250.;

  xNLO[++i]=900;
  yNLO[i]=230.;

  xNLO[++i]=1000;
  yNLO[i]=220.;

  xNLO[++i]=1400;
  yNLO[i]=180.;

  xNLO[++i]=1800;
  yNLO[i]=156.;

  xNLO[++i]=1900;
  yNLO[i]=150.;
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
    x[j]  = xNLO[j] + 10;
    y[j]  = yNLO[j] + 10;
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
