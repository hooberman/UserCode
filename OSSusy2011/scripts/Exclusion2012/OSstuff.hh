

TGraph *observedLimit_OS2011(){

  const unsigned int no = 21;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 125;  yo[i] = 435;
  xo[++i] = 150;  yo[i] = 420;
  xo[++i] = 180;  yo[i] = 380;
  xo[++i] = 190;  yo[i] = 350;
  xo[++i] = 195;  yo[i] = 330;
  xo[++i] = 210;  yo[i] = 320;
  xo[++i] = 250;  yo[i] = 330;
  xo[++i] = 300;  yo[i] = 340;
  xo[++i] = 450;  yo[i] = 360;
  xo[++i] = 580;  yo[i] = 360;
  xo[++i] = 600;  yo[i] = 350;
  xo[++i] = 620;  yo[i] = 300;
  xo[++i] = 650;  yo[i] = 285;
  xo[++i] = 700;  yo[i] = 275;
  xo[++i] = 830;  yo[i] = 260; //
  xo[++i] = 1000; yo[i] = 230;
  xo[++i] = 1400; yo[i] = 200; //
  xo[++i] = 2000; yo[i] = 190; //
  xo[++i] = 3000; yo[i] = 185; //

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph *expectedLimit_OS2011(){

  const unsigned int no = 19;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 125;  yo[i] = 435;
  xo[++i] = 150;  yo[i] = 420;
  xo[++i] = 180;  yo[i] = 380;
  xo[++i] = 190;  yo[i] = 350;
  xo[++i] = 195;  yo[i] = 330;
  xo[++i] = 210;  yo[i] = 320;
  xo[++i] = 250;  yo[i] = 330;
  xo[++i] = 300;  yo[i] = 340;
  xo[++i] = 450;  yo[i] = 350;
  xo[++i] = 570;  yo[i] = 360;
  xo[++i] = 590;  yo[i] = 350;
  xo[++i] = 580;  yo[i] = 300;
  xo[++i] = 800;  yo[i] = 250;
  xo[++i] = 1000; yo[i] = 220;
  xo[++i] = 1400; yo[i] = 190;//
  xo[++i] = 2000; yo[i] = 185;//
  xo[++i] = 3000; yo[i] = 185;//

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;

}
  
TGraph *expectedLimitP1_OS2011(){

  const unsigned int no = 23;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 410;
  xo[++i] = 80;   yo[i] = 410;
  xo[++i] = 90;   yo[i] = 400;
  xo[++i] = 100;  yo[i] = 390;
  xo[++i] = 100;  yo[i] = 370;
  xo[++i] = 120;  yo[i] = 365;
  xo[++i] = 145;  yo[i] = 360;
  xo[++i] = 180;  yo[i] = 350;
  xo[++i] = 185;  yo[i] = 315;
  xo[++i] = 190;  yo[i] = 280;
  xo[++i] = 210;  yo[i] = 275;
  xo[++i] = 240;  yo[i] = 280;
  xo[++i] = 400;  yo[i] = 320;
  xo[++i] = 460;  yo[i] = 320;
  xo[++i] = 480;  yo[i] = 270;
  xo[++i] = 550;  yo[i] = 245;
  xo[++i] = 600;  yo[i] = 240;
  xo[++i] = 750;  yo[i] = 230;
  xo[++i] = 900;  yo[i] = 210;
  xo[++i] = 1200; yo[i] = 150;
  xo[++i] = 1600; yo[i] = 135;
  xo[++i] = 2000; yo[i] = 130;
  xo[++i] = 3000; yo[i] = 120;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}
  
TGraph *expectedLimitM1_OS2011(){

  const unsigned int no = 19;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 200;  yo[i] = 420;
  xo[++i] = 220;  yo[i] = 400;
  xo[++i] = 225;  yo[i] = 360;
  xo[++i] = 245;  yo[i] = 357;
  xo[++i] = 265;  yo[i] = 355;
  xo[++i] = 300;  yo[i] = 360;
  xo[++i] = 400;  yo[i] = 370;
  xo[++i] = 450;  yo[i] = 375;
  xo[++i] = 550;  yo[i] = 377;
  xo[++i] = 640;  yo[i] = 380;
  xo[++i] = 730;  yo[i] = 300;
  xo[++i] = 900;  yo[i] = 270;
  xo[++i] = 1000; yo[i] = 250;
  xo[++i] = 1200; yo[i] = 240;
  xo[++i] = 1400; yo[i] = 230;
  xo[++i] = 2000; yo[i] = 220;
  xo[++i] = 3000; yo[i] = 210;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph *observedLimitTheoryUp_OS2011(){


  const unsigned int no = 19;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 200;  yo[i] = 420;
  xo[++i] = 220;  yo[i] = 400;
  xo[++i] = 225;  yo[i] = 350;
  xo[++i] = 245;  yo[i] = 347;
  xo[++i] = 265;  yo[i] = 345;
  xo[++i] = 300;  yo[i] = 350;
  xo[++i] = 400;  yo[i] = 370;
  xo[++i] = 450;  yo[i] = 375;
  xo[++i] = 550;  yo[i] = 372;
  xo[++i] = 640;  yo[i] = 370;
  xo[++i] = 730;  yo[i] = 300;
  xo[++i] = 900;  yo[i] = 270;
  xo[++i] = 1000; yo[i] = 250;
  xo[++i] = 1200; yo[i] = 230;
  xo[++i] = 1500; yo[i] = 210;
  xo[++i] = 2000; yo[i] = 200;
  xo[++i] = 3000; yo[i] = 200;

  for( int ip = 0 ; ip < 18 ; ip++ ) yo[ip] += 5;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;

}

TGraph *observedLimitTheoryDown_OS2011(){

  const unsigned int no = 21;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 440;
  xo[++i] = 100;  yo[i] = 440;
  xo[++i] = 125;  yo[i] = 425;
  xo[++i] = 150;  yo[i] = 400;
  xo[++i] = 180;  yo[i] = 360;
  xo[++i] = 190;  yo[i] = 330;
  xo[++i] = 195;  yo[i] = 310;
  xo[++i] = 210;  yo[i] = 300;
  xo[++i] = 250;  yo[i] = 310;
  xo[++i] = 300;  yo[i] = 320;
  xo[++i] = 450;  yo[i] = 340;
  xo[++i] = 550;  yo[i] = 340;
  xo[++i] = 570;  yo[i] = 330;
  xo[++i] = 620;  yo[i] = 280;
  xo[++i] = 630;  yo[i] = 275;
  xo[++i] = 700;  yo[i] = 265;
  xo[++i] = 830;  yo[i] = 250;
  xo[++i] = 1000; yo[i] = 220;
  xo[++i] = 1400; yo[i] = 185;
  xo[++i] = 2000; yo[i] = 175;
  xo[++i] = 3000; yo[i] = 175;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
  
}


TGraph *expectedLimitTheoryUp_OS2011(){

  const unsigned int no = 19;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 450;
  xo[++i] = 100;  yo[i] = 450;
  xo[++i] = 180;  yo[i] = 420;
  xo[++i] = 200;  yo[i] = 400;
  xo[++i] = 205;  yo[i] = 350;
  xo[++i] = 225;  yo[i] = 347;
  xo[++i] = 245;  yo[i] = 345;
  xo[++i] = 300;  yo[i] = 350;
  xo[++i] = 400;  yo[i] = 370;
  xo[++i] = 450;  yo[i] = 375;
  xo[++i] = 550;  yo[i] = 372;
  xo[++i] = 620;  yo[i] = 350;
  xo[++i] = 650;  yo[i] = 310;
  xo[++i] = 800;  yo[i] = 260;
  xo[++i] = 1000; yo[i] = 240;
  xo[++i] = 1200; yo[i] = 220;
  xo[++i] = 1400; yo[i] = 200;
  xo[++i] = 2000; yo[i] = 190;
  xo[++i] = 3000; yo[i] = 190;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph* expectedLimitTheoryDown_OS2011(){
  const unsigned int no = 18;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 430;
  xo[++i] = 100;  yo[i] = 420;
  xo[++i] = 110;  yo[i] = 410;
  xo[++i] = 120;  yo[i] = 390;
  xo[++i] = 170;  yo[i] = 360;
  xo[++i] = 190;  yo[i] = 320;
  xo[++i] = 230;  yo[i] = 315;
  xo[++i] = 320;  yo[i] = 324;
  xo[++i] = 400;  yo[i] = 330;
  xo[++i] = 490;  yo[i] = 330;
  xo[++i] = 560;  yo[i] = 270;
  xo[++i] = 620;  yo[i] = 260;
  xo[++i] = 800;  yo[i] = 230;
  xo[++i] = 1000; yo[i] = 200;
  xo[++i] = 1200; yo[i] = 180;
  xo[++i] = 1400; yo[i] = 170;
  xo[++i] = 2000; yo[i] = 170;
  xo[++i] = 3000; yo[i] = 170;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;

}


TGraphErrors* getNLOobsTanbeta10_2010( ){
  
  const Int_t n=47;
  Double_t x[n];// ,y[n],xerr2[n],yerr2[n];
  Double_t y[n];
  Double_t xerr2[n];
  Double_t yerr2[n];

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

  x[++i]=60.;
  y[i]=225.;

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

  x[++i]=105.;
  y[i]=230.;

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
  
  TGraphErrors* grtb10  = new TGraphErrors(n,x, y,xerr2,yerr2);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  grtb10->SetLineColor(6);
  grtb10->SetLineWidth(3);
  return grtb10;

}




  
TGraph *expectedLimitBand_OS2011(){

  const unsigned int no1 = 24;
  float xo1[no1];
  float yo1[no1];

  int i = -1;

  xo1[++i] = 60;   yo1[i] = 410;
  xo1[++i] = 80;   yo1[i] = 410;
  xo1[++i] = 90;   yo1[i] = 400;
  xo1[++i] = 100;  yo1[i] = 390;
  xo1[++i] = 100;  yo1[i] = 370;

  xo1[++i] = 120;  yo1[i] = 365;
  xo1[++i] = 145;  yo1[i] = 360;
  xo1[++i] = 180;  yo1[i] = 350;
  xo1[++i] = 185;  yo1[i] = 315;
  xo1[++i] = 190;  yo1[i] = 280;

  xo1[++i] = 210;  yo1[i] = 275;
  xo1[++i] = 240;  yo1[i] = 280;
  xo1[++i] = 400;  yo1[i] = 320;
  xo1[++i] = 460;  yo1[i] = 320;
  xo1[++i] = 480;  yo1[i] = 270;

  xo1[++i] = 550;  yo1[i] = 245;
  xo1[++i] = 600;  yo1[i] = 240;
  xo1[++i] = 750;  yo1[i] = 230;
  xo1[++i] = 900;  yo1[i] = 210;
  xo1[++i] = 1200; yo1[i] = 150;

  xo1[++i] = 1600; yo1[i] = 135;
  xo1[++i] = 2000; yo1[i] = 130;
  xo1[++i] = 3000; yo1[i] = 120;

  xo1[++i] = 6000; yo1[i] = 270;


  const unsigned int no2 = 19;
  float xo2[no2];
  float yo2[no2];

  i = -1;

  xo2[++i] = 60;   yo2[i] = 450;
  xo2[++i] = 100;  yo2[i] = 450;
  xo2[++i] = 200;  yo2[i] = 420;
  xo2[++i] = 220;  yo2[i] = 400;
  xo2[++i] = 225;  yo2[i] = 360;

  xo2[++i] = 245;  yo2[i] = 357;
  xo2[++i] = 265;  yo2[i] = 355;
  xo2[++i] = 300;  yo2[i] = 360;
  xo2[++i] = 400;  yo2[i] = 370;
  xo2[++i] = 450;  yo2[i] = 375;

  xo2[++i] = 550;  yo2[i] = 377;
  xo2[++i] = 640;  yo2[i] = 380;
  xo2[++i] = 730;  yo2[i] = 300;
  xo2[++i] = 900;  yo2[i] = 270;
  xo2[++i] = 1000; yo2[i] = 250;

  xo2[++i] = 1200; yo2[i] = 240; 
  xo2[++i] = 1400; yo2[i] = 230; 
  xo2[++i] = 2000; yo2[i] = 220; 
  xo2[++i] = 3000; yo2[i] = 210;

  /*
  const unsigned int no1 = 21;
  float xo1[no1];
  float yo1[no1];

  int i = -1;

  xo1[++i] = 60;   yo1[i] = 410;
  xo1[++i] = 80;   yo1[i] = 410;
  xo1[++i] = 90;   yo1[i] = 400;
  xo1[++i] = 100;  yo1[i] = 390;
  xo1[++i] = 100;  yo1[i] = 370;
  xo1[++i] = 120;  yo1[i] = 365;
  xo1[++i] = 145;  yo1[i] = 360;
  xo1[++i] = 180;  yo1[i] = 350;
  xo1[++i] = 185;  yo1[i] = 315;
  xo1[++i] = 190;  yo1[i] = 280;
  xo1[++i] = 210;  yo1[i] = 275;
  xo1[++i] = 240;  yo1[i] = 280;
  xo1[++i] = 400;  yo1[i] = 320;
  xo1[++i] = 460;  yo1[i] = 320;
  xo1[++i] = 480;  yo1[i] = 270;
  xo1[++i] = 550;  yo1[i] = 245;
  xo1[++i] = 600;  yo1[i] = 240;
  xo1[++i] = 800;  yo1[i] = 210;
  xo1[++i] = 1200; yo1[i] = 150;
  xo1[++i] = 2000; yo1[i] = 140;

  xo1[++i] = 200000; yo1[i] = 170;

  const unsigned int no2 = 18;
  float xo2[no2];
  float yo2[no2];

  i = -1;

  xo2[++i] = 60;   yo2[i] = 450;
  xo2[++i] = 100;  yo2[i] = 450;
  xo2[++i] = 200;  yo2[i] = 420;
  xo2[++i] = 220;  yo2[i] = 400;
  xo2[++i] = 225;  yo2[i] = 350;
  xo2[++i] = 245;  yo2[i] = 347;
  xo2[++i] = 265;  yo2[i] = 345;
  xo2[++i] = 300;  yo2[i] = 350;
  xo2[++i] = 400;  yo2[i] = 370;
  xo2[++i] = 450;  yo2[i] = 375;
  xo2[++i] = 550;  yo2[i] = 372;
  xo2[++i] = 640;  yo2[i] = 370;
  xo2[++i] = 730;  yo2[i] = 300;
  xo2[++i] = 900;  yo2[i] = 270;
  xo2[++i] = 1000; yo2[i] = 250;
  xo2[++i] = 1200; yo2[i] = 230;
  xo2[++i] = 1400; yo2[i] = 220;
  xo2[++i] = 2000; yo2[i] = 210;
  */

  const unsigned int no = no1 + no2;
  float xo[no];
  float yo[no];

  for( unsigned int ip = 0 ; ip < no1 ; ip++ ){
    xo[ip] = xo1[ip];
    yo[ip] = yo1[ip];
  }

  for( unsigned int ip = 0 ; ip < no2 ; ip++ ){
    xo[ip+no1] = xo2[no2-1-ip];
    yo[ip+no1] = yo2[no2-1-ip];
  }
  

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);
  gro->SetFillColor(2);
  
  return gro;
}
