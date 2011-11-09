{

  //---------------------------------------
  // load TChain
  //---------------------------------------
  
  TChain *LM = new TChain("t");
  //LM->Add("../output/V00-02-06/highpt/LM6_smallTree.root");
  LM->Add("/tas/benhoob/home/OSSusy2011/output/JulyPAS/highpt/LM1_smallTree.root");

  //---------------------------------------
  // selection
  //---------------------------------------

  TCut myweight("weight * 3.5");
  TCut presel("pfmet>50 && njets>=2 && ht>100 && !passz");
  TCut SR1("ht>300&&ht<600&&pfmet>275");
  TCut SR2("ht>600&&pfmet>275");
  TCut SR3("ht>600&&pfmet>200&&pfmet<275");
  TCut SF("leptype==0 || leptype==1");
  TCut OF("leptype==2");

  //---------------------------------------
  // preselection and SR1,SR2,SR3 yields
  //---------------------------------------

  TH1F* h = new TH1F("h","h",1,0,1);
  h->Sumw2();

  cout << endl;
  LM->Draw("0.5>>h",presel*myweight);
  cout << "presel   : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  cout << endl;
  LM->Draw("0.5>>h",(presel+SR1)*myweight);
  cout << "SR1      : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  LM->Draw("0.5>>h",(presel+SR2)*myweight);
  cout << "SR2      : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  LM->Draw("0.5>>h",(presel+SR3)*myweight);
  cout << "SR3      : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  LM->Draw("0.5>>h",(presel+(SR1||SR2||SR3))*myweight);
  cout << "allSR    : " << Form("%.1f +/- %.1f",h->GetBinContent(1),h->GetBinError(1)) << endl; 

  //---------------------------------------
  // yields in 6 bins for shape analysis
  //---------------------------------------

  TH1F* h1SF = new TH1F("h1SF","h1SF",1,0,1);
  TH1F* h1OF = new TH1F("h1OF","h1OF",1,0,1);
  TH1F* h2SF = new TH1F("h2SF","h2SF",1,0,1);
  TH1F* h2OF = new TH1F("h2OF","h2OF",1,0,1);
  TH1F* h3SF = new TH1F("h3SF","h3SF",1,0,1);
  TH1F* h3OF = new TH1F("h3OF","h3OF",1,0,1);

  h1SF->Sumw2();
  h1OF->Sumw2();
  h2SF->Sumw2();
  h2OF->Sumw2();
  h3SF->Sumw2();
  h3OF->Sumw2();

  //------nominal-------//
  LM->Draw("0.5>>h1SF",(presel+SR1+SF)*myweight);
  LM->Draw("0.5>>h1OF",(presel+SR1+OF)*myweight);
  LM->Draw("0.5>>h2SF",(presel+SR2+SF)*myweight);
  LM->Draw("0.5>>h2OF",(presel+SR2+OF)*myweight);
  LM->Draw("0.5>>h3SF",(presel+SR3+SF)*myweight);
  LM->Draw("0.5>>h3OF",(presel+SR3+OF)*myweight);
  
  cout << endl << endl;
  cout << "float   LM_yield[nbins]           = {  " 
       << Form("%.1f",h1SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h1OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3OF->GetBinContent(1)) << "  }" << endl;
  //cout << endl << endl;

  //------JESup-------//
  TString preselupstring(presel);
  TString SR1upstring(SR1);
  TString SR2upstring(SR2);
  TString SR3upstring(SR3);

  preselupstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  preselupstring.ReplaceAll("ht"     ,  "htUp");
  preselupstring.ReplaceAll("njets"  ,  "njetsUp");

  SR1upstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  SR1upstring.ReplaceAll("ht"     ,  "htUp");
  SR1upstring.ReplaceAll("njets"  ,  "njetsUp");

  SR2upstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  SR2upstring.ReplaceAll("ht"     ,  "htUp");
  SR2upstring.ReplaceAll("njets"  ,  "njetsUp");

  SR3upstring.ReplaceAll("pfmet"  ,  "pfmetUp");
  SR3upstring.ReplaceAll("ht"     ,  "htUp");
  SR3upstring.ReplaceAll("njets"  ,  "njetsUp");

  TCut preselup(preselupstring);
  TCut SR1up(SR1upstring);
  TCut SR2up(SR2upstring);
  TCut SR3up(SR3upstring);

  LM->Draw("0.5>>h1SF",(preselup+SR1up+SF)*myweight);
  LM->Draw("0.5>>h1OF",(preselup+SR1up+OF)*myweight);
  LM->Draw("0.5>>h2SF",(preselup+SR2up+SF)*myweight);
  LM->Draw("0.5>>h2OF",(preselup+SR2up+OF)*myweight);
  LM->Draw("0.5>>h3SF",(preselup+SR3up+SF)*myweight);
  LM->Draw("0.5>>h3OF",(preselup+SR3up+OF)*myweight);
  
  //cout << endl << endl;
  cout << "float   LM_yield_JESup[nbins]      = {  " 
       << Form("%.1f",h1SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h1OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3OF->GetBinContent(1)) << "  }" << endl;
  //cout << endl << endl;

  //------JESdn-------//
  TString preseldnstring(presel);
  TString SR1dnstring(SR1);
  TString SR2dnstring(SR2);
  TString SR3dnstring(SR3);

  preseldnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  preseldnstring.ReplaceAll("ht"     ,  "htDown");
  preseldnstring.ReplaceAll("njets"  ,  "njetsDown");

  SR1dnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  SR1dnstring.ReplaceAll("ht"     ,  "htDown");
  SR1dnstring.ReplaceAll("njets"  ,  "njetsDown");

  SR2dnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  SR2dnstring.ReplaceAll("ht"     ,  "htDown");
  SR2dnstring.ReplaceAll("njets"  ,  "njetsDown");

  SR3dnstring.ReplaceAll("pfmet"  ,  "pfmetDown");
  SR3dnstring.ReplaceAll("ht"     ,  "htDown");
  SR3dnstring.ReplaceAll("njets"  ,  "njetsDown");

  TCut preseldn(preseldnstring);
  TCut SR1dn(SR1dnstring);
  TCut SR2dn(SR2dnstring);
  TCut SR3dn(SR3dnstring);

  LM->Draw("0.5>>h1SF",(preseldn+SR1dn+SF)*myweight);
  LM->Draw("0.5>>h1OF",(preseldn+SR1dn+OF)*myweight);
  LM->Draw("0.5>>h2SF",(preseldn+SR2dn+SF)*myweight);
  LM->Draw("0.5>>h2OF",(preseldn+SR2dn+OF)*myweight);
  LM->Draw("0.5>>h3SF",(preseldn+SR3dn+SF)*myweight);
  LM->Draw("0.5>>h3OF",(preseldn+SR3dn+OF)*myweight);
  
  //cout << endl << endl;
  cout << "float   LM_yield_JESdn[nbins]      = {  " 
       << Form("%.1f",h1SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h1OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h2OF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3SF->GetBinContent(1)) << "   ,   "
       << Form("%.1f",h3OF->GetBinContent(1)) << "  }" << endl;
  cout << endl << endl;








}




  // cout << endl;
  // cout << "SR1 SF   : " << Form("%.1f +/- %.1f",h1SF->GetBinContent(1),h1SF->GetBinError(1)) << endl;
  // cout << "SR1 OF   : " << Form("%.1f +/- %.1f",h1OF->GetBinContent(1),h1OF->GetBinError(1)) << endl;
  // cout << "SR2 SF   : " << Form("%.1f +/- %.1f",h2SF->GetBinContent(1),h2SF->GetBinError(1)) << endl;
  // cout << "SR2 OF   : " << Form("%.1f +/- %.1f",h2OF->GetBinContent(1),h2OF->GetBinError(1)) << endl;
  // cout << "SR3 SF   : " << Form("%.1f +/- %.1f",h3SF->GetBinContent(1),h3SF->GetBinError(1)) << endl;
  // cout << "SR3 OF   : " << Form("%.1f +/- %.1f",h3OF->GetBinContent(1),h3OF->GetBinError(1)) << endl;
