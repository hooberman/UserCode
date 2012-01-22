
TGraph *observedLimit_OS2011(){

  const unsigned int no = 18;
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
  xo[++i] = 590;  yo[i] = 300;
  xo[++i] = 800;  yo[i] = 260;
  xo[++i] = 1000; yo[i] = 230;
  xo[++i] = 1400; yo[i] = 190;
  xo[++i] = 2000; yo[i] = 180;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph *expectedLimit_OS2011(){

  const unsigned int no = 18;
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
  xo[++i] = 1400; yo[i] = 180;
  xo[++i] = 2000; yo[i] = 165;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;

}
  
TGraph *expectedLimitP1_OS2011(){

  const unsigned int no = 20;
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
  xo[++i] = 800;  yo[i] = 210;
  xo[++i] = 1200; yo[i] = 150;
  xo[++i] = 2000; yo[i] = 140;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}
  
TGraph *expectedLimitM1_OS2011(){

  const unsigned int no = 18;
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
  xo[++i] = 1400; yo[i] = 220;
  xo[++i] = 2000; yo[i] = 210;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph *observedLimitTheoryUp_OS2011(){
  const unsigned int no = 18;
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
  xo[++i] = 1400; yo[i] = 220;
  xo[++i] = 2000; yo[i] = 210;

  for( int ip = 0 ; ip < 18 ; ip++ ) yo[ip] += 5;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;

}

TGraph *observedLimitTheoryDown_OS2011(){

  const unsigned int no = 17;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 430;
  xo[++i] = 100;  yo[i] = 420;
  xo[++i] = 110;  yo[i] = 410;
  xo[++i] = 120;  yo[i] = 390;
  xo[++i] = 180;  yo[i] = 360;
  xo[++i] = 190;  yo[i] = 310;
  xo[++i] = 230;  yo[i] = 305;
  xo[++i] = 320;  yo[i] = 310;
  xo[++i] = 400;  yo[i] = 320;
  xo[++i] = 500;  yo[i] = 330;
  xo[++i] = 550;  yo[i] = 270;
  xo[++i] = 600;  yo[i] = 265;
  xo[++i] = 800;  yo[i] = 240;
  xo[++i] = 1000; yo[i] = 220;
  xo[++i] = 1200; yo[i] = 190;
  xo[++i] = 1400; yo[i] = 170;
  xo[++i] = 2000; yo[i] = 170;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}


TGraph *expectedLimitTheoryUp_OS2011(){

  const unsigned int no = 18;
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
  xo[++i] = 2000; yo[i] = 200;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;
}

TGraph* expectedLimitTheoryDown_OS2011(){
  const unsigned int no = 17;
  float xo[no];
  float yo[no];

  int i = -1;

  xo[++i] = 60;   yo[i] = 430;
  xo[++i] = 100;  yo[i] = 420;
  xo[++i] = 110;  yo[i] = 410;
  xo[++i] = 120;  yo[i] = 390;
  xo[++i] = 170;  yo[i] = 360;
  xo[++i] = 190;  yo[i] = 310;
  xo[++i] = 230;  yo[i] = 305;
  xo[++i] = 320;  yo[i] = 310;
  xo[++i] = 400;  yo[i] = 320;
  xo[++i] = 490;  yo[i] = 320;
  xo[++i] = 540;  yo[i] = 270;
  xo[++i] = 600;  yo[i] = 260;
  xo[++i] = 800;  yo[i] = 230;
  xo[++i] = 1000; yo[i] = 200;
  xo[++i] = 1200; yo[i] = 180;
  xo[++i] = 1400; yo[i] = 170;
  xo[++i] = 2000; yo[i] = 170;

  TGraph* gro = new TGraph(no,xo,yo);
  gro->SetLineWidth(5);
  gro->SetMarkerColor(2);

  return gro;

}

