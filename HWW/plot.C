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
#include "TCanvas.h"
#include <iomanip>

using namespace std;


void plot(){

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
  TCut weight     = "event_scale1fb * 0.0355";
  TCut sel        = pt2010 + met_projpt + jetveto + mll12;
  TCut hsel;

  char* iter = "v2";

//   TChain *bkg = new TChain("Events");
// //   bkg->Add(Form("babies/%s/WWTo2L2Nu_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/GluGluToWWTo4L_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/WZ_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/ZZ_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/TTJets_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/tW_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/WJetsToLNu_PU_testFinal_baby.root",iter));
// //   bkg->Add(Form("babies/%s/DYToMuMuM20_PU_testFinal_baby.root",iter) );
// //   bkg->Add(Form("babies/%s/DYToMuMuM10To20_PU_testFinal_baby.root",iter) );
// //   bkg->Add(Form("babies/%s/DYToEEM20_PU_testFinal_baby.root",iter) );
// //   bkg->Add(Form("babies/%s/DYToEEM10To20_PU_testFinal_baby.root",iter) );
// //   bkg->Add(Form("babies/%s/DYToTauTauM20_PU_testFinal_baby.root",iter) );
// //   bkg->Add(Form("babies/%s/DYToTauTauM10To20_PU_testFinal_baby.root",iter) );

  TChain *wwbkg = new TChain("Events");
  wwbkg->Add(Form("babies/%s/WWTo2L2Nu_PU_testFinal_baby.root",iter));

  TChain *wjetsbkg = new TChain("Events");
  wjetsbkg->Add(Form("babies/%s/WJetsToLNu_PU_testFinal_baby.root",iter));

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


  TCanvas *c1 = new TCanvas();
  c1->cd();

  TH1F* hdphi_ww     = new TH1F("hdphi_ww","",    20,0,3.14);
  TH1F* hdphi_wjets  = new TH1F("hdphi_wjets","", 20,0,3.14);
  TH1F* hdphi_h130   = new TH1F("hdphi_h130","",  20,0,3.14);

  wwbkg->Draw   ("dil_dphi>>hdphi_ww", sel * weight );
  wjetsbkg->Draw("dil_dphi>>hdphi_wjets", sel * weight );
  sig->Draw     ("dil_dphi>>hdphi_h130", sel * weight );

  cout << "WW yield "     << hdphi_ww->Integral() << endl;
  cout << "W+jets yield " << hdphi_wjets->Integral() << endl;
  cout << "Higgs yield  " << hdphi_h130->Integral() << endl;

  hdphi_wjets->SetLineColor(2);
  hdphi_h130->SetLineColor(4);


  hdphi_ww->Draw();
  hdphi_wjets->Draw("same");
  hdphi_h130->Draw("same");


  TCanvas *c2 = new TCanvas();
  c2->cd();

  TH1F* hmll_ww     = new TH1F("hmll_ww","",    20,0,100);
  TH1F* hmll_wjets  = new TH1F("hmll_wjets","", 20,0,100);
  TH1F* hmll_h130   = new TH1F("hmll_h130","",  20,0,100);

  wwbkg->Draw   ("dil_mass>>hmll_ww", sel * weight );
  wjetsbkg->Draw("dil_mass>>hmll_wjets", sel * weight );
  sig->Draw     ("dil_mass>>hmll_h130", sel * weight );

  cout << "WW yield "     << hmll_ww->Integral() << endl;
  cout << "W+jets yield " << hmll_wjets->Integral() << endl;
  cout << "Higgs yield  " << hmll_h130->Integral() << endl;

  hmll_wjets->SetLineColor(2);
  hmll_h130->SetLineColor(4);


  hmll_ww->Draw();
  hmll_wjets->Draw("same");
  hmll_h130->Draw("same");

}
