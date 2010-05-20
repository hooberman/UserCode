void tcmetValidation(){
  
  //ps output file name
  char* psFileName  = "tcmetValidation.ps";

  //2 files to compare
  TFile *f1=TFile::Open("/tmp/benhoob/CMSSW_3_6_0_pre1_test.root");
  TFile *f2=TFile::Open("/tmp/benhoob/CMSSW_3_6_0_pre1.root");

  //labels for TLegend
  string leg1 = "CMSSW_3_5_8";
  string leg2 = "CMSSW_3_6_1";

  //in root file, path to tcmet histos
  char* tcmetpath = "DQMData/RecoMETV/MET_Global/tcMet/METTask";
  
  //list variables for comparison
  vector<char*> vars;
  vars.push_back("MET");
  vars.push_back("SumET");
  vars.push_back("MEx");
  vars.push_back("MEy");
  /*
  vars.push_back("dMET");
  vars.push_back("dMETx");
  vars.push_back("dMETy");
  vars.push_back("dMUx");
  vars.push_back("dMUy");
  vars.push_back("CorrectionFlag");
  vars.push_back("METPhi");
  vars.push_back("METPhiResolution_GenMETTrue");
  vars.push_back("METResolution_GenMETTrue");
  vars.push_back("Nevents");
  vars.push_back("electronHoverE");
  vars.push_back("fracTracks");
  vars.push_back("muonEta");
  vars.push_back("muonNormalizedChi2");
  vars.push_back("muonSAhits");
  vars.push_back("nMus");
  vars.push_back("trackAlgo");
  vars.push_back("trackEta");
  vars.push_back("trackNormalizedChi2");
  vars.push_back("trackPtErr");
  vars.push_back("METPhiResolution_GenMETCalo");
  vars.push_back("METResolution_GenMETCalo");
  vars.push_back("METSig");
  vars.push_back("MExCorrection");
  vars.push_back("MEyCorrection");
  vars.push_back("electronEta");
  vars.push_back("electronPt");
  vars.push_back("muonD0");
  vars.push_back("muonNhits");
  vars.push_back("muonPt");
  vars.push_back("nEls");
  vars.push_back("nMusAsPis");
  vars.push_back("trackD0");
  vars.push_back("trackNhits");
  vars.push_back("trackPt");
  vars.push_back("trackQuality");
  */

  //make TLegend
  TH1F* h1dummy = new TH1F("h1dummy","",1,0,1);
  h1dummy->SetLineColor(2);
  h1dummy->SetMarkerColor(2);
  TH1F* h2dummy = new TH1F("h2dummy","",1,0,1);
  h2dummy->SetLineColor(4);
  h2dummy->SetMarkerColor(4);
  TLegend *leg = new TLegend(0.8,0.8,0.9,0.9);
  leg->AddEntry(h1dummy , leg1.c_str());
  leg->AddEntry(h2dummy , leg2.c_str());
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  //declare canvases, histos, etc
  const int nvar = vars.size();
  TH1F    *h1[nvar];
  TH1F    *h2[nvar];
  TPad    *mainpad[nvar];
  TPad    *pullpad[nvar];
  TLine line;
  line.SetLineStyle(2);
  TLatex t;
  t.SetNDC();

  TCanvas *canvas = new TCanvas("canvas","canvas",800,800);
  canvas->Print(Form("%s[",psFileName));
  canvas->Clear();
  
  //loop over variables
  for(int ivar=0;ivar<nvar;ivar++){
    
    //make canvas and pad
    canvas->cd();
    mainpad[ivar] = new TPad(Form("%s_mainpad",vars.at(ivar)),Form("%s_mainpad",vars.at(ivar)),0,0,1,0.8);
    mainpad[ivar] -> Draw();
    mainpad[ivar] -> cd();
    
    //get histos
    h1[ivar] = (TH1F*) f1->Get(Form("%s_%s",tcmetpath,vars.at(ivar)));
    h2[ivar] = (TH1F*) f2->Get(Form("%s_%s",tcmetpath,vars.at(ivar)));

    h1[ivar] -> Rebin(5);
    h2[ivar] -> Rebin(5);

    //format and draw histos
    if(drawlog(vars.at(ivar))) mainpad[ivar]->SetLogy(1); 
    h1[ivar] -> SetLineColor(2);
    h2[ivar] -> SetLineColor(4);
    h1[ivar] -> Draw();
    h1[ivar] -> GetXaxis() -> SetTitle(vars.at(ivar));
    h2[ivar] -> Draw("same");
    leg->Draw();
    
    //make canvas and pad
    canvas->cd();
    pullpad[ivar] = new TPad(Form("%s_pullpad",vars.at(ivar)),Form("%s_pullpad",vars.at(ivar)),0,0.8,1,1.);
    pullpad[ivar] -> Draw();
    pullpad[ivar] -> cd();

    //format and draw pull hist
    TH1F* hpull = getPullHist(h1[ivar],h2[ivar]);
    hpull->Draw();
    hpull->GetXaxis()->SetLabelSize(0);
    hpull->GetXaxis()->SetTitleSize(0);
    hpull->GetYaxis()->SetTitleSize(0.16);
    hpull->GetYaxis()->SetLabelSize(0.16);
    hpull->GetYaxis()->SetRangeUser(-4,4);
    hpull->GetYaxis()->SetTitle("Pull");
    hpull->GetYaxis()->SetNdivisions(5);

    //draw guidelines
    line.DrawLine(hpull->GetXaxis()->GetXmin(), 0, hpull->GetXaxis()->GetXmax(), 0);
    line.DrawLine(hpull->GetXaxis()->GetXmin(), 1, hpull->GetXaxis()->GetXmax(), 1);
    line.DrawLine(hpull->GetXaxis()->GetXmin(),-1, hpull->GetXaxis()->GetXmax(),-1);

    canvas->Update();
    canvas->Print(psFileName);
      
  }
  
  canvas->Update();
  canvas->Print(Form("%s]",psFileName));
  
  

}

TH1F* getPullHist(TH1F* h1, TH1F* h2){
  
  TH1F* hout = (TH1F*) h1->Clone(Form("%s_clone",h1->GetName()));
  
  for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){
  
    float val = h2->GetBinContent(ibin) - h1->GetBinContent(ibin);
    float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
    if(fabs(err) < 1.e-10)  err = sqrt(h2->GetBinContent(ibin) + h1->GetBinContent(ibin));
    
    hout -> SetBinContent(ibin,fabs(err) > 0 ? val/err : val);
    //hout -> SetBinError(ibin,1);
  }
  
  return hout;
}

bool drawlog(char* prefix){
  
  if(strcmp(prefix,"MET") == 0)      return true;
  
  return false; 
}

inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

