#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;



void makeThisCard(ofstream* dofile, char* cardname, int nobs, float nbkg, float sbkg){

  ofstream* ofile = new ofstream();

  ofile->open(Form("%s.txt",cardname));

  *ofile <<  "imax 1  number of channels"                                                    << endl;
  *ofile <<  "jmax 1  number of backgrounds"                                                 << endl;
  *ofile <<  "kmax 2  number of nuisance parameters (sources of systematical uncertainties)" << endl;
  *ofile <<  "------------"                                                                  << endl;
  *ofile <<  "bin         1"                                                                 << endl;
  *ofile <<  Form("observation %i",nobs)                                                     << endl;
  *ofile <<  "------------"                                                                  << endl;
  *ofile <<  "bin             1      1"                                                      << endl;
  *ofile <<  "process       signal  background"                                              << endl; 
  *ofile <<  "process         0      1"                                                      << endl;
  *ofile <<  Form("rate            1    %.1f",nbkg)                << endl;
  *ofile <<  "------------" << endl;
  *ofile <<  "deltaS  lnN   1.10    -     uncertainty on signal" << endl;
  *ofile <<  Form("deltaB  lnN     -   %.2f    uncertainty on background",sbkg) << endl; 

  ofile->close();

  *dofile << Form("../test/lands.exe -d %s.txt -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --freq --nToysForCLsb 10000 --nToysForCLb 5000 --seed 1234 -n %s",cardname,cardname) << endl;

}




void makeCards(){

  ofstream* dofile = new ofstream();
  
  dofile->open("doLimits_T2tt.sh");

  // makeThisCard(dofile,"LM150" , 240 , 240.0 , 1.20);
  // makeThisCard(dofile,"LM200" ,  87 ,  87.1 , 1.31);
  // makeThisCard(dofile,"LM250" ,  32 ,  31.5 , 1.32);
  // makeThisCard(dofile,"LM300" ,  11 ,  11.0 , 1.38);

  // makeThisCard(dofile,"HM150" ,  30 ,  29.5 , 1.29);
  // makeThisCard(dofile,"HM200" ,  17 ,  16.6 , 1.33);
  // makeThisCard(dofile,"HM250" ,  10 ,   9.7 , 1.36);
  // makeThisCard(dofile,"HM300" ,   5 ,   4.8 , 1.38);

  // makeThisCard(dofile,"BDT1L" , 741 , 740.9 , 1.12);
  // makeThisCard(dofile,"BDT1T" , 120 , 119.4 , 1.18);
  // makeThisCard(dofile,"BDT2"  ,  81 ,  80.8 , 1.21);
  // makeThisCard(dofile,"BDT3"  ,  13 ,  13.0 , 1.36);
  // makeThisCard(dofile,"BDT4"  ,   3 ,   2.7 , 1.49);
  // makeThisCard(dofile,"BDT5"  ,  84 ,  83.6 , 1.41);

  // MT > 150 GeV numbers

  // http://www.t2.ucsd.edu/tastwiki/pub/CMS/StopPredictions/PredResults_CNCLM_MT150.pdf
  makeThisCard(dofile,"LM150_MT150" , 151 , 151.1 , 1.28);
  makeThisCard(dofile,"LM200_MT150" ,  56 ,  56.2 , 1.33);
  makeThisCard(dofile,"LM250_MT150" ,  23 ,  23.1 , 1.34);
  makeThisCard(dofile,"LM300_MT150" ,   9 ,   9.2 , 1.39);

  // http://www.t2.ucsd.edu/tastwiki/pub/CMS/StopPredictions/PredResults_CNCHM_MT150.pdf
  makeThisCard(dofile,"HM150_MT150" ,  18 ,  18.2 , 1.34);
  makeThisCard(dofile,"HM200_MT150" ,  12 ,  11.8 , 1.36);
  makeThisCard(dofile,"HM250_MT150" ,   8 ,   7.5 , 1.39);
  makeThisCard(dofile,"HM300_MT150" ,   4 ,   4.1 , 1.42);

  // http://www.t2.ucsd.edu/tastwiki/pub/CMS/StopPredictions/PredResults_BDT_T2tt_MT150.pdf
  makeThisCard(dofile,"BDT1L_MT150" , 450 , 449.8 , 1.19);
  makeThisCard(dofile,"BDT1T_MT150" ,  83 ,  83.0 , 1.23);
  makeThisCard(dofile,"BDT2_MT150"  ,  62 ,  62.1 , 1.25);
  makeThisCard(dofile,"BDT3_MT150"  ,  11 ,  10.7 , 1.36);
  makeThisCard(dofile,"BDT4_MT150"  ,   3 ,   2.5 , 1.49);
  makeThisCard(dofile,"BDT5_MT150"  ,  56 ,  55.5 , 1.33);





  /*

  dofile->open("doLimits_T2bw.sh");
  
  // from http://www.t2.ucsd.edu/tastwiki/pub/CMS/20130323OsChats/PredResults_CNCLM_T2bw.pdf
  // makeThisCard(dofile,"T2bw_LM100"  ,  1670 ,   1670 , 1.09);
  // makeThisCard(dofile,"T2bw_LM150"  ,   519 ,    519 , 1.13);
  // makeThisCard(dofile,"T2bw_LM200"  ,   180 ,    180 , 1.17);
  // makeThisCard(dofile,"T2bw_LM250"  ,    66 ,     66 , 1.23);
  // makeThisCard(dofile,"T2bw_LM300"  ,    24 ,     24 , 1.29);
  // makeThisCard(dofile,"T2bw_LM350"  ,    10 ,    9.8 , 1.36);
  // makeThisCard(dofile,"T2bw_LM400"  ,     4 ,    4.2 , 1.48);

  // http://www.t2.ucsd.edu/tastwiki/pub/CMS/20130323OsChats/PredResults_CNCHM_T2bw.pdf
  makeThisCard(dofile,"T2bw_HM100"  ,    78 ,   78.2 , 1.18);
  makeThisCard(dofile,"T2bw_HM150"  ,    38 ,   37.7 , 1.22);
  makeThisCard(dofile,"T2bw_HM200"  ,    19 ,   18.5 , 1.27);
  makeThisCard(dofile,"T2bw_HM250"  ,    10 ,    9.5 , 1.30);

  // http://www.t2.ucsd.edu/tastwiki/pub/CMS/20130323OsChats/PredResults_BDT_T2bw_0p5.pdf
  // makeThisCard(dofile,"T2bw_BDT1"   ,    76 ,   75.6 , 1.17);
  // makeThisCard(dofile,"T2bw_BDT2"   ,    50 ,   49.9 , 1.20);
  // makeThisCard(dofile,"T2bw_BDT3"   ,    15 ,   15.3 , 1.29);

  // makeThisCard(dofile,"T2bw_LM100_unc15"  ,  1670 ,   1670 , 1.15);
  // makeThisCard(dofile,"T2bw_LM150_unc15"  ,   519 ,    519 , 1.15);
  // makeThisCard(dofile,"T2bw_LM200_unc15"  ,   180 ,    180 , 1.15);
  // makeThisCard(dofile,"T2bw_LM250_unc15"  ,    66 ,     66 , 1.15);
  // makeThisCard(dofile,"T2bw_LM300_unc15"  ,    24 ,     24 , 1.15);
  // makeThisCard(dofile,"T2bw_LM350_unc15"  ,    10 ,    9.8 , 1.15);
  // makeThisCard(dofile,"T2bw_LM400_unc15"  ,     4 ,    4.2 , 1.15);

  // makeThisCard(dofile,"T2bw_HM100_unc15"  ,    32 ,   31.6 , 1.15);
  // makeThisCard(dofile,"T2bw_HM150_unc15"  ,    13 ,   13.2 , 1.15);
  // makeThisCard(dofile,"T2bw_HM200_unc15"  ,     6 ,    5.8 , 1.15);
  // makeThisCard(dofile,"T2bw_HM250_unc15"  ,     3 ,    2.6 , 1.15);
  */  



}
