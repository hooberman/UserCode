TGraph * Left25_observed(){
  double xx_obs[] = {100 , 100 , 100 , 125 , 150 , 175 , 200 , 225 , 250 , 275 , 300 , 325 , 350 , 375 , 400 , 425 , 450 , 475 , 475 , 475 , 500 , 500 , 500 , 500 , 500 , 500 , 500 , 475 , 450 , 425 , 400 , 375 , 350 , 325 , 300 , 275 , 250 , 225 , 200 , 175 , 150 , 125 , 100 , 100 , 100 , 125};
  double yy_obs[] = {0 , 25 , 50 , 75 , 100 , 125 , 150 , 175 , 200 , 200 , 225 , 250 , 250 , 250 , 250 , 250 , 250 , 225 , 200 , 175 , 150 , 125 , 100 , 75 , 50 , 25 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , 25 , 50 , 75};
  TGraph* observed = new TGraph(46, xx_obs, yy_obs);
  return observed;
}

TGraph * Left25_observedup(){
  double xx_theoryup[] = {100 , 100 , 100 , 121.464 , 146.464 , 171.464 , 196.464 , 221.464 , 247.764 , 272.764 , 296.464 , 322.764 , 350 , 375 , 400 , 425 , 452.236 , 479.472 , 480 , 479.472 , 504.472 , 505 , 505 , 505 , 505 , 505 , 504.16 , 476.213 , 451.213 , 426.213 , 401.213 , 376.213 , 351.213 , 326.213 , 301.213 , 276.213 , 251.213 , 226.213 , 201.213 , 176.213 , 151.213 , 125.981 , 100 , 100 , 100};
  double yy_theoryup[] = {0 , 23.7873 , 52.2361 , 78.5355 , 103.536 , 128.536 , 153.536 , 178.536 , 204.472 , 204.472 , 228.536 , 254.472 , 255 , 255 , 255 , 255 , 254.472 , 227.236 , 200 , 177.236 , 152.236 , 125 , 100 , 75 , 50 , 25 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , 23.7873 , 52.2361};
  TGraph* observed_theoryup = new TGraph(45, xx_theoryup, yy_theoryup);
  return observed_theoryup;
}

TGraph * Left25_observeddown(){
  double xx_theorydown[] = {100 , 100 , 100 , 128.536 , 153.536 , 178.536 , 203.536 , 228.536 , 252.236 , 277.236 , 303.536 , 327.236 , 350 , 375 , 400 , 425 , 447.764 , 470.528 , 470 , 470.528 , 495.528 , 495 , 495 , 495 , 495 , 495 , 495.84 , 473.787 , 448.787 , 423.787 , 398.787 , 373.787 , 348.787 , 323.787 , 298.787 , 273.787 , 248.787 , 223.787 , 198.787 , 173.787 , 148.787 , 124.019 , 100 , 100 , 100};
  double yy_theorydown[] = {0 , 26.2127 , 47.7639 , 71.4645 , 96.4645 , 121.464 , 146.464 , 171.464 , 195.528 , 195.528 , 221.464 , 245.528 , 245 , 245 , 245 , 245 , 245.528 , 222.764 , 200 , 172.764 , 147.764 , 125 , 100 , 75 , 50 , 25 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , -50 , 26.2127 , 47.7639};
  TGraph* observed_theorydown = new TGraph(45, xx_theorydown, yy_theorydown);
  return observed_theorydown;
}
