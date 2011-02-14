#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCut.h"
#include "TRandom3.h"
#include <iomanip>

using namespace std;


void efficiency(){

  int mH = 130;
  //int mH = 160;
  //int mH = 200;

  TCut met_projpt = "(event_type != 2 && met_projpt>35. ) || (event_type == 2 && met_projpt > 20.)";
  TCut pt2020     = "lephard_pt > 20 && lepsoft_pt > 20";
  TCut pt2010     = "lephard_pt > 20 && lepsoft_pt > 10";
  TCut jetveto    = "jets_num==0 && extralep_num==0 && lowptbtags_num==0 && softmu_num == 0";
  TCut mll12      = "dil_mass > 12.";
  TCut h130       = "dil_dphi < 1.05 && dil_mass < 45. && lephard_pt > 25 && lepsoft_pt > 20";
  TCut h160       = "dil_dphi < 1.05 && dil_mass < 50. && lephard_pt > 30 && lepsoft_pt > 25";
  TCut h200       = "dil_dphi < 1.75 && dil_mass < 90. && lephard_pt > 40 && lepsoft_pt > 25";
  //TCut weight     = "event_scale1fb * 0.0355";
  TCut weight     = "1";
  TCut sel        = pt2010 + met_projpt + jetveto + mll12;
  TCut hsel;

  char* iter = "v2";

  TChain *bkg = new TChain("Events");
  //bkg->Add(Form("babies/%s/WWTo2L2Nu_PU_testFinal_baby.root",iter));
//   bkg->Add(Form("babies/%s/GluGluToWWTo4L_PU_testFinal_baby.root",iter));
//   bkg->Add(Form("babies/%s/WZ_PU_testFinal_baby.root",iter));
//   bkg->Add(Form("babies/%s/ZZ_PU_testFinal_baby.root",iter));
//   bkg->Add(Form("babies/%s/TTJets_PU_testFinal_baby.root",iter));
//   bkg->Add(Form("babies/%s/tW_PU_testFinal_baby.root",iter));
  bkg->Add(Form("babies/%s/WJetsToLNu_PU_testFinal_baby.root",iter));
//   bkg->Add(Form("babies/%s/DYToMuMuM20_PU_testFinal_baby.root",iter) );
//   bkg->Add(Form("babies/%s/DYToMuMuM10To20_PU_testFinal_baby.root",iter) );
//   bkg->Add(Form("babies/%s/DYToEEM20_PU_testFinal_baby.root",iter) );
//   bkg->Add(Form("babies/%s/DYToEEM10To20_PU_testFinal_baby.root",iter) );
//   bkg->Add(Form("babies/%s/DYToTauTauM20_PU_testFinal_baby.root",iter) );
//   bkg->Add(Form("babies/%s/DYToTauTauM10To20_PU_testFinal_baby.root",iter) );

  TChain* sig = new TChain("Events");

  if( mH == 130 ){
     sig->Add(Form("babies/%s/HToWWTo2L2NuM130_PU_testFinal_baby.root",iter));
     sig->Add(Form("babies/%s/HToWWToLNuTauNuM130_PU_testFinal_baby.root",iter));
     sig->Add(Form("babies/%s/HToWWTo2Tau2NuM130_PU_testFinal_baby.root",iter));
     hsel = sel + h130;
   }
   else if( mH == 160 ){
     sig->Add(Form("babies/%s/HToWWTo2L2NuM160_PU_testFinal_baby.root",iter));
     sig->Add(Form("babies/%s/HToWWToLNuTauNuM160_PU_testFinal_baby.root",iter));
     sig->Add(Form("babies/%s/HToWWTo2Tau2NuM160_PU_testFinal_baby.root",iter));
     hsel = sel + h160;
   }
   else if( mH == 200 ){
     sig->Add(Form("babies/%s/HToWWTo2L2NuM200_PU_testFinal_baby.root",iter));
     sig->Add(Form("babies/%s/HToWWToLNuTauNuM200_PU_testFinal_baby.root",iter));
     sig->Add(Form("babies/%s/HToWWTo2Tau2NuM200_PU_testFinal_baby.root",iter));
     hsel = sel + h200;
  }
   else{
     std::cout << "Error, unrecognized higgs mass " << mH << " GeV, quitting" << std::endl;
     exit(0);
   }


  TH1F* htotbkg  = new TH1F("htotbkg","",1,0,1);
  TH1F* htotsig  = new TH1F("htotsig","",1,0,1);
  TH1F* hpassbkg = new TH1F("hpassbkg","",1,0,1);
  TH1F* hpasssig = new TH1F("hpasssig","",1,0,1);

  bkg->Draw("0.5>>htotbkg", sel * weight );
  sig->Draw("0.5>>htotsig", sel * weight );

  bkg->Draw("0.5>>hpassbkg", hsel * weight);
  sig->Draw("0.5>>hpasssig", hsel * weight);

  cout << "total WW      " << htotbkg->Integral() << endl;
  cout << "total H->WW   " << htotsig->Integral() << endl;
  cout << "pass WW       " << hpassbkg->Integral() << endl;
  cout << "pass H->WW    " << hpasssig->Integral() << endl;
  cout << "H->WW eff     " << hpasssig->Integral() / htotsig->Integral() << endl;
  cout << "WW eff        " << hpassbkg->Integral() / htotbkg->Integral() << endl;

}
