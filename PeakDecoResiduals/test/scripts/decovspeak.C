TGraphErrors* averageTGraph(TGraphErrors* g1, TGraphErrors *g2);

inline double fround(double n, unsigned d){
  return floor(n * pow(10., d) + .5) / pow(10., d);
}

bool draw(char* var, vector<char*> cantitles){
  
  for(int i=0;i<cantitles.size();i++) {
    if(strcmp(var,cantitles.at(i)) == 0) return true;
  }
  return false;
}

void decovspeak(){

  //Config params-----------------------------------------------
  bool printeps=false;
  bool printgif=true;
  bool fit=false;

  string dir="crabjobs/trial3";

  vector<char*> cantitles;
  cantitles.push_back("du_dw");                    //delta u, delta w TH1s
  //cantitles.push_back("dtquantities");             //DT time TH1s
  //cantitles.push_back("charge_nstrips");           //charge, nstrips TH1s
  //cantitles.push_back("deltawEC");                 //delta w for EC TH1s
  //cantitles.push_back("duvsdtantheta");            //du vs. delta tan theta TH2
  //cantitles.push_back("dwvsdttime");               //dw vs. dt time TH2
  //cantitles.push_back("chargevsdttime");           //charge vs dt time TH2
  cantitles.push_back("duvsdtanthetasplit");       //v+/- split du vs. delta tan(theta) TH2s
  //cantitles.push_back("duvstrkthetasplit");        //v+/- split du vs. tan(theta_trk)   TH2
  //cantitles.push_back("dwvsdt_tgraph");            //dw vs. dt time TGraphs
  //cantitles.push_back("chargevsdttime_tgraph");    //charge vs. dt time TGraphs
  cantitles.push_back("duvsdtantheta_tgraph");     //v+/- split du vs. delta tan(theta) TGraphs
  //cantitles.push_back("duvstrktheta_tgraph");      //v+/- split du vs. tan(theta_trk)   TGraphs
  //cantitles.push_back("duvsdtanthetaall_tgraph");  //du vs. delta tan(theta)         
  //cantitles.push_back("duvsdtanthetaproj");        //v+/- split du vs. delta tan(theta) & projections
  //cantitles.push_back("duvstrkthetaproj");         //v+/- split du vs. tan(thetatrk) & projections 
  //cantitles.push_back("duvstrkthetaprojy");        //v+/v- split du vs. dtantheta & Y projections
  //cantitles.push_back("dtanladt");                   //delta tan(LA) vs. DT time
  //cantitles.push_back("deltatanLAdt");             //time dependence of delta tan(LA)
    
  const unsigned int ncan = cantitles.size();
  TCanvas *can[ncan];
  int idx = 0;

  //open root files---------------------------------------------
  string decofile    = dir+"/root/deco.root";
  string peakfile    = dir+"/root/peak.root";
  TFile *fdeco = TFile::Open(decofile.c_str());
  TFile *fpeak = TFile::Open(peakfile.c_str());
  
  //load macros and set style-----------------------------------
  gROOT->LoadMacro("scripts/cloneHist.C");
  gROOT->LoadMacro("scripts/setStats.C");
  gROOT->LoadMacro("scripts/plotHists.C"); 
  gROOT->LoadMacro("scripts/plotHistsTH2.C");
  gROOT->ProcessLine(".L scripts/getTGraphFromTH2.C+");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat("mr");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  //make legend----------------------------------------
  int dx=1200;
  int dy=600;
  TH1F *h1=new TH1F("h1","h1",1,0,1);
  h1->SetLineColor(2);
  h1->SetMarkerColor(2);
  TH1F *h2=new TH1F("h2","h2",1,0,1);
  h2->SetMarkerColor(4);

  TLegend *leg1=new TLegend(0.15,0.55,0.35,0.85);
  leg1->AddEntry(h1,"PEAK");
  leg1->AddEntry(h2,"DECO");
  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);

  //get histos-----------------------------------------------
  fpeak->cd();
  PeakDecoResiduals->cd();
  PeakDecoResiduals->cd();
  TH1F *dupeak =            (TH1F*)du_TOB->Clone();
  TH1F *dwpeak =            (TH1F*)dw_TOB->Clone();
  TH1F *dttimepeak =        (TH1F*)dttime_TOB->Clone();
  TH1F *dttimeerrpeak =     (TH1F*)dttimeerr_TOB->Clone();
  TH1F *ndtpeak =           (TH1F*)ndt_TOB->Clone();
  TH1F *chargepeak =        (TH1F*)charge_TOB->Clone(); 
  TH1F *nstripspeak =       (TH1F*)nstrips_TOB->Clone(); 
  TH2F* duthetapeak =       (TH2F*)duvsdtantheta_TOB->Clone();
  TH2F* duvsdtpeak =        (TH2F*)duvsdttime_TOB->Clone();
  TH2F* dwvsdtpeak =        (TH2F*)dwvsdttime_TOB->Clone();
  TH2F* chvsdtpeak =        (TH2F*)chargevsdttime_TOB->Clone();
  TH2F* duthetapeakp =      (TH2F*)duvsdtantheta_TOB_vp_wp->Clone();
  duthetapeakp       -> Add((TH2F*)duvsdtantheta_TOB_vp_wm->Clone());
  TH2F* duthetapeakm =      (TH2F*)duvsdtantheta_TOB_vm_wp->Clone();
  duthetapeakm       -> Add((TH2F*)duvsdtantheta_TOB_vm_wm->Clone());
  TH2F* dutrkpeakp   =      (TH2F*)duvstrktheta_TOB_vp_wp->Clone();
  dutrkpeakp         -> Add((TH2F*)duvstrktheta_TOB_vp_wm->Clone());
  TH2F* dutrkpeakm   =      (TH2F*)duvstrktheta_TOB_vm_wp->Clone();
  dutrkpeakm         -> Add((TH2F*)duvstrktheta_TOB_vm_wm->Clone());
  
  TH2F* h_vp_peak[8];
  TH2F* h_vm_peak[8];
  
  for(int ih=0;ih<8;ih++){
    h_vp_peak[ih] = (TH2F*) fpeak->
      Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_vp_dt_TOB_%i",ih));
    h_vm_peak[ih] = (TH2F*) fpeak->
      Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_vm_dt_TOB_%i",ih));
  }

  fdeco->cd();
  PeakDecoResiduals->cd();
  PeakDecoResiduals->cd();
  TH1F *dudeco =            (TH1F*)du_TOB->Clone();
  TH1F *dwdeco =            (TH1F*)dw_TOB->Clone();
  TH1F *dttimedeco =        (TH1F*)dttime_TOB->Clone();
  TH1F *dttimeerrdeco =     (TH1F*)dttimeerr_TOB->Clone();
  TH1F *ndtdeco =           (TH1F*)ndt_TOB->Clone();
  TH1F *chargedeco =        (TH1F*)charge_TOB->Clone(); 
  TH1F *nstripsdeco =       (TH1F*)nstrips_TOB->Clone(); 
  TH2F* duthetadeco =       (TH2F*)duvsdtantheta_TOB->Clone();
  TH2F* duvsdtdeco =        (TH2F*)duvsdttime_TOB->Clone();
  TH2F* dwvsdtdeco =        (TH2F*)dwvsdttime_TOB->Clone();
  TH2F* chvsdtdeco =        (TH2F*)chargevsdttime_TOB->Clone();
  TH2F* duthetadecop =      (TH2F*)duvsdtantheta_TOB_vp_wp->Clone();
  duthetadecop       -> Add((TH2F*)duvsdtantheta_TOB_vp_wm->Clone());
  TH2F* duthetadecom =      (TH2F*)duvsdtantheta_TOB_vm_wp->Clone();
  duthetadecom       -> Add((TH2F*)duvsdtantheta_TOB_vm_wm->Clone());
  TH2F* dutrkdecop   =      (TH2F*)duvstrktheta_TOB_vp_wp->Clone();
  dutrkdecop         -> Add((TH2F*)duvstrktheta_TOB_vp_wm->Clone());
  TH2F* dutrkdecom   =      (TH2F*)duvstrktheta_TOB_vm_wp->Clone();
  dutrkdecom         -> Add((TH2F*)duvstrktheta_TOB_vm_wm->Clone());

  TH2F* h_vp_deco[8];
  TH2F* h_vm_deco[8];
  
  for(int ih=0;ih<8;ih++){
    h_vp_deco[ih] = (TH2F*) fdeco->
      Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_vp_dt_TOB_%i",ih));
    h_vm_deco[ih] = (TH2F*) fdeco->
      Get(Form("PeakDecoResiduals/PeakDecoResiduals/duvsdtantheta_vm_dt_TOB_%i",ih));

  }

  //delta u, delta w TH1s---------------------------------------------------

  if(draw("du_dw",cantitles)){

    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),dx,dy);
    can[idx]->Divide(2,1);
    
    can[idx]->cd(1);
    setStats(dupeak,dudeco,0.5,0.65,0.15,false);  
    plotHists(dupeak,dudeco,"TOB","#Delta u [#mum]",-500,500,1,3);
    leg1->Draw();
    can[idx]->cd(2);
    setStats(dwpeak,dwdeco,0.5,0.65,0.15,false);  
    plotHists(dwpeak,dwdeco,"TOB","#Delta w [#mum]",-1000,1000,1,3);
    leg1->Draw();

    idx++;
  }
  
  //DT quantities------------------------------------------------------------

  if(draw("dtquantities",cantitles)){
 
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),dx,2*dy);
    can[idx]->Divide(2,2);
    
    can[idx]->cd(1);
    setStats(dttimepeak,dttimedeco,0.5,0.65,0.15,false);  
    plotHists(dttimepeak,dttimedeco,"TOB","DT time (ns)",-50,50,1,0);
    leg1->Draw();
  
    can[idx]->cd(2);
    setStats(dttimeerrpeak,dttimeerrdeco,0.5,0.65,0.15,false);  
    plotHists(dttimeerrpeak,dttimeerrdeco,"TOB","#delta DT time (ns)",0,100,1,0);
  
    TLegend *leg2=new TLegend(0.4,0.55,0.6,0.85);
    leg2->AddEntry(h1,"PEAK");
    leg2->AddEntry(h2,"DECO");
    leg2->SetBorderSize(1);
    leg2->SetFillColor(0);
    leg2->Draw();
    
    can[idx]->cd(3);
    setStats(ndtpeak,ndtdeco,0.5,0.65,0.15,false);  
    plotHists(ndtpeak,ndtdeco,"TOB","n DT Hits",0,50,1,0);
    leg1->Draw();
    
    idx++;
  }
	
  //charge, nstrips------------------------------------------------------------

  if(draw("charge_nstrips",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),dx,dy);
    can[idx]->Divide(2,1);

    can[idx]->cd(1);
    setStats(chargepeak,chargedeco,0.5,0.65,0.15,false);  
    plotHists(chargepeak,chargedeco,"TOB","Charge",0,1000,1,0);
    can[idx]->cd(2);
    setStats(nstripspeak,nstripsdeco,0.5,0.65,0.15,false);  
    plotHists(nstripspeak,nstripsdeco,"TOB","NStrips",0,10,1,0);

    idx++;
  }
  

  //du vs. delta tan theta TH2-------------------------------------------------

  if(draw("",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),dx,dy);
    can[idx]->Divide(2,1);
    
    plotHistsTH2(duthetapeak,
		 "TOB PEAK","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],1);
    plotHistsTH2(duthetadeco,
		 "TOB DECO","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],2);
  
    idx++;
  }
  
  //dw vs. dt time TH2---------------------------------------------------------

  if(draw("",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),dx,dy);
    can[idx]->Divide(2,1);
  
    plotHistsTH2(dwvsdtpeak,
		 "TOB PEAK","DT Time [ns]",
		 "#Delta w [#mum]",-25,25,-500,500,0,can[idx],1);
    plotHistsTH2(dwvsdtdeco,
		 "TOB DECO","DT Time [ns]",
		 "#Delta w [#mum]",-25,25,-500,500,0,can[idx],2);

    idx++;
  }

  //charge vs. dt time TH2------------------------------------------------------
  
  if(draw("",cantitles)){
     
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),dx,dy);
    can[idx]->Divide(2,1);

    plotHistsTH2(chvsdtpeak,
		 "TOB PEAK","DT Time [ns]",
		 "Charge",-25,25,0,500,0,can[idx],1);
    plotHistsTH2(chvsdtdeco,
		 "TOB DECO","DT Time [ns]",
		 "Charge",-25,25,0,500,0,can[idx],2);

    idx++;
  }
  
  //du vs. delta tan(theta) split v+/v- histos-----------------------------------

  if(draw("",cantitles)){

    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,900);
    can[idx]->Divide(2,2);

    plotHistsTH2(duthetapeakp,
		 "TOB PEAK (v+)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],1);
    plotHistsTH2(duthetadecop,
		 "TOB DECO (v+)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],2);
    plotHistsTH2(duthetapeakm,
		 "TOB PEAK (v-)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],3);
    plotHistsTH2(duthetadecom,
		 "TOB DECO (v-)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],4);
    idx++;
  }

  //du vs. tan(theta_trk) split v+/v- histos-------------------------------------
  
  if(draw("",cantitles)){
  
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,900);
    can[idx]->Divide(2,2);
 
    plotHistsTH2(dutrkpeakp,
		 "TOB PEAK (v+)","tan(#theta_{trk})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],1);
    plotHistsTH2(dutrkdecop,
		 "TOB DECO (v+)","tan(#theta_{trk})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],2);
    plotHistsTH2(dutrkpeakm,
		 "TOB PEAK (v-)","tan(#theta_{trk})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],3);
    plotHistsTH2(dutrkdecom,
		 "TOB DECO (v-)","tan(#theta_{trk})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[idx],4);

    idx++;
  }

  //divide du vs. DT time histos into DT slices
  const int nslices=8;
  const float xmin[]={-20.,-12.,-8.,-4.,0.,4., 8.,12.};
  const float xmax[]={-12., -8.,-4., 0.,4.,8.,12.,20.};

  float xpeak[nslices];
  float xpeakerr[nslices];
  float xdeco[nslices];
  float xdecoerr[nslices];
  float ypeak[nslices];
  float ypeakerr[nslices];
  float ydeco[nslices];
  float ydecoerr[nslices];

  float xchpeak[nslices];
  float xchpeakerr[nslices];
  float xchdeco[nslices];
  float xchdecoerr[nslices];
  float ychpeak[nslices];
  float ychpeakerr[nslices];
  float ychdeco[nslices];
  float ychdecoerr[nslices];

  TH2 *dwvsdtpeakslice[nslices];
  TH2 *dwvsdtdecoslice[nslices];
  TH1 *dwpeakproj[nslices];
  TH1 *dwdecoproj[nslices];
  TH1 *dtpeakproj[nslices];
  TH1 *dtdecoproj[nslices];

  TH2 *chvsdtpeakslice[nslices];
  TH2 *chvsdtdecoslice[nslices];
  TH1 *chpeakproj[nslices];
  TH1 *chdecoproj[nslices];
  TH1 *dtpeakproj2[nslices];
  TH1 *dtdecoproj2[nslices];

  //   string projtitle[nslices]={"TOB (-20 ns < dt < -12 ns)",
  // 			     "TOB (-12 ns < dt < -8 ns)",
  // 			     "TOB (-8 ns < dt < -4 ns)",
  // 			     "TOB (-4 ns < dt < 0 ns)",
  // 			     "TOB (0 ns < dt < 4 ns)",
  // 			     "TOB (4 ns < dt < 8 ns)",
  // 			     "TOB (8 ns < dt < 12 ns)",
  // 			     "TOB (12 ns < dt < 20 ns)"};
  

  vector<float> dtbins;
  dtbins.push_back(-20.);
  dtbins.push_back(-12.);
  dtbins.push_back(-8.);
  dtbins.push_back(-4.);
  dtbins.push_back(0.);
  dtbins.push_back(4.);
  dtbins.push_back(8.);
  dtbins.push_back(12.);
  dtbins.push_back(20.);

  TLine *line=new TLine;
  
  //TGraph <dw> vs. <dt time>------------------------------------------------------

  if(draw("",cantitles)){
   
    TCanvas *can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),600,450);
    can[idx]->cd();

    dwvsdtpeak->SetName("dwvsdtpeak");
    TGraphErrors *gpeak = getTGraphFromTH2(dwvsdtpeak,dtbins,0);
    gpeak->SetMarkerColor(2);
    gpeak->SetLineColor(2);
    gpeak->SetMarkerStyle(8);
    gpeak->SetMarkerSize(1);
    TF1 *fp=new TF1("fp","pol1");
    fp->SetLineColor(2);
    gpeak->Fit(fp);
    stringstream speak;
    speak<<"Slope = "<<fround(fp->GetParameter(1),2)<<" #mum / ns"<<endl;
    
    dwvsdtdeco->SetName("dwvsdtdeco");
    TGraphErrors *gdeco = getTGraphFromTH2(dwvsdtdeco,dtbins,0);
    gdeco->SetMarkerColor(4);
    gdeco->SetLineColor(4);
    gdeco->SetMarkerStyle(8);
    gdeco->SetMarkerSize(1);
    TF1 *fd=new TF1("fd","pol1");
    fd->SetLineColor(4);
    gdeco->Fit(fd);
    stringstream sdeco;
    sdeco<<"Slope = "<<fround(fd->GetParameter(1),2)<<" #mum / ns"<<endl;
    
    TMultiGraph *mg=new TMultiGraph();
    mg->Add(gpeak);
    mg->Add(gdeco);
    mg->SetTitle("TOB");
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("<DT Time> [ns]");
    mg->GetYaxis()->SetTitle("<#Deltaw> [#mum]");
    mg->SetTitle("TOB");
    
    TLatex *t=new TLatex();
    t->SetNDC();
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.75,speak.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.65,sdeco.str().c_str());

    idx++;
  }
  
  //TGraph <charge> vs. <dt time>--------------------------------------------------

  if(draw("",cantitles)){

    TCanvas *can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),800,600);
    can[idx]->cd();

    chvsdtdeco->SetName("chvsdtpeak");
    TGraphErrors *gchpeak = getTGraphFromTH2(chvsdtpeak,dtbins,0);
    gchpeak->SetMarkerColor(2);
    gchpeak->SetLineColor(2);
    gchpeak->SetMarkerStyle(8);
    gchpeak->SetMarkerSize(1);
    
    chvsdtdeco->SetName("chvsdtdeco");
    TGraphErrors *gchdeco = getTGraphFromTH2(chvsdtdeco,dtbins,0);
    gchdeco->SetMarkerColor(4);
    gchdeco->SetLineColor(4);
    gchdeco->SetMarkerStyle(8);
    gchdeco->SetMarkerSize(1);
    
    TMultiGraph *mgch=new TMultiGraph();
    mgch->Add(gchpeak);
    mgch->Add(gchdeco);
    mgch->SetTitle("TOB");
    mgch->Draw("AP");
    mgch->GetXaxis()->SetTitle("<DT Time> [ns]");
    mgch->GetYaxis()->SetTitle("<Cluster Charge>");
    mgch->SetTitle("TOB");

    idx++;
  }
 

  //TGraph du vs. delta tan theta----------------------------------------
  if(draw("duvsdtantheta_tgraph",cantitles)){
    
    bool addbpcorr = false; //add tgraph for BP-corrected deco
    float size = 0.1;
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
  
    can[idx]->cd(1);
    vector<float> thetabins;
//     thetabins.push_back(-1.);
//     thetabins.push_back(-0.75);
//     thetabins.push_back(-0.5);
//     thetabins.push_back(-0.25);
//     thetabins.push_back(0.);
//     thetabins.push_back(0.25);
//     thetabins.push_back(0.5);
//     thetabins.push_back(0.75);
//     thetabins.push_back(1.);

    thetabins.push_back(-0.9);
    thetabins.push_back(-0.7);
    thetabins.push_back(-0.5);
    thetabins.push_back(-0.3);
    thetabins.push_back(-0.1);
    thetabins.push_back(0.1);
    thetabins.push_back(0.3);
    thetabins.push_back(0.5);
    thetabins.push_back(0.7);
    thetabins.push_back(0.9);

//     thetabins.push_back(-0.9);
//     thetabins.push_back(-0.8);
//     thetabins.push_back(-0.7);
//     thetabins.push_back(-0.6);
//     thetabins.push_back(-0.5);
//     thetabins.push_back(-0.4);
//     thetabins.push_back(-0.3);
//     thetabins.push_back(-0.2);
//     thetabins.push_back(-0.1);
//     thetabins.push_back(0.);
//     thetabins.push_back(0.1);
//     thetabins.push_back(0.2);
//     thetabins.push_back(0.3);
//     thetabins.push_back(0.4);
//     thetabins.push_back(0.5);
//     thetabins.push_back(0.6);
//     thetabins.push_back(0.7);
//     thetabins.push_back(0.8);
//     thetabins.push_back(0.9);
    
    duthetapeakp->SetName("duthetapeakp");
    TGraphErrors *gduthetapeakp = getTGraphFromTH2(duthetapeakp,thetabins,0);
    gduthetapeakp->SetMarkerColor(2);
    gduthetapeakp->SetLineColor(2);
    gduthetapeakp->SetMarkerStyle(8);
    gduthetapeakp->SetMarkerSize(size);
    TF1 *fduthetapeakp=new TF1("fduthetapeakp","pol1");
    fduthetapeakp->SetLineColor(2);
    gduthetapeakp->Fit(fduthetapeakp);
    stringstream sduthetapeakp;
    sduthetapeakp<<"m = "<<fround(fduthetapeakp->GetParameter(1),3)
		 <<" #pm " <<fround(fduthetapeakp->GetParError(1),3)
		 <<"   b = "<<fround(fduthetapeakp->GetParameter(0),3)
		 <<" #pm "<<fround(fduthetapeakp->GetParError(0),3)
		 <<endl;
    
    duthetadecop->SetName("duthetadecop");
    TGraphErrors *gduthetadecop = getTGraphFromTH2(duthetadecop,thetabins,0);
    gduthetadecop->SetMarkerColor(4);
    gduthetadecop->SetLineColor(4);
    gduthetadecop->SetMarkerStyle(8);
    gduthetadecop->SetMarkerSize(size);
    TF1 *fduthetadecop=new TF1("fduthetadecop","pol1");
    fduthetadecop->SetLineColor(4);
    gduthetadecop->Fit(fduthetadecop);
    stringstream sduthetadecop;
    sduthetadecop<<"m = "<<fround(fduthetadecop->GetParameter(1),3)
		 <<" #pm " <<fround(fduthetadecop->GetParError(1),3)
		 <<"   b = "<<fround(fduthetadecop->GetParameter(0),3)
		 <<" #pm "<<fround(fduthetadecop->GetParError(0),3)
		 <<endl;

    //add deco with BP correction------------------------------------------
    if(addbpcorr){
      TFile fdecobp("crabjobs/trial17/root/deco.root");
      fdecobp.cd();
      PeakDecoResiduals.cd();
      PeakDecoResiduals.cd();
      TH2F* duthetadecobp_vp = (TH2F*)duvsdtantheta_vp_TOB->Clone();
      TH2F* duthetadecobp_vm = (TH2F*)duvsdtantheta_vm_TOB->Clone();
      
      duthetadecobp_vp->SetName("duthetadecobp_vp");
      TGraphErrors *gduthetadecobp_vp = getTGraphFromTH2(duthetadecobp_vp,thetabins,0);
      gduthetadecobp_vp->SetMarkerColor(6);
      gduthetadecobp_vp->SetLineColor(6);
      gduthetadecobp_vp->SetMarkerStyle(8);
      gduthetadecobp_vp->SetMarkerSize(size);
      TF1 *fduthetadecobp_vp=new TF1("fduthetadecobp_vp","pol1");
      fduthetadecobp_vp->SetLineColor(6);
      gduthetadecobp_vp->Fit(fduthetadecobp_vp);
      stringstream sduthetadecobp_vp;
    
      sduthetadecobp_vp<<"m = "<<fround(fduthetadecobp_vp->GetParameter(1),3)
		       <<" #pm " <<fround(fduthetadecobp_vp->GetParError(1),3)
		       <<"   b = "<<fround(fduthetadecobp_vp->GetParameter(0),3)
		       <<" #pm "<<fround(fduthetadecobp_vp->GetParError(0),3)
		       <<endl;
    }

    TMultiGraph *mgdtp=new TMultiGraph();
    mgdtp->Add(gduthetapeakp);
    mgdtp->Add(gduthetadecop);
    if(addbpcorr) mgdtp->Add(gduthetadecobp_vp);
    mgdtp->SetTitle("TOB (v+)");
    mgdtp->Draw("AP");
    mgdtp->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
    mgdtp->GetYaxis()->SetTitle("<#Deltau> [#mum]");
    mgdtp->SetTitle("TOB");
    //mgdtp->GetXaxis()->SetLimits(-1,1);
    mgdtp->GetYaxis()->SetLimits(-20,20);
    mgdtp->SetMinimum(-20);
    mgdtp->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    //t->DrawLatex(0.5,0.8,"y = mx + b");
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.85,sduthetapeakp.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.75,sduthetadecop.str().c_str());
    if(addbpcorr){
      t->SetTextColor(6);
      t->DrawLatex(0.15,0.65,sduthetadecobp_vp.str().c_str());
    }

    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);
    
    can[idx]->cd(2);
    duthetapeakm->SetName("duthetapeakm");
    TGraphErrors *gduthetapeakm = getTGraphFromTH2(duthetapeakm,thetabins,0);
    gduthetapeakm->SetMarkerColor(2);
    gduthetapeakm->SetLineColor(2);
    gduthetapeakm->SetMarkerStyle(8);
    gduthetapeakm->SetMarkerSize(size);
    TF1 *fduthetapeakm=new TF1("fduthetapeakm","pol1");
    fduthetapeakm->SetLineColor(2);
    gduthetapeakm->Fit(fduthetapeakm);
    stringstream sduthetapeakm;
    sduthetapeakm<<"m = "<<fround(fduthetapeakm->GetParameter(1),3)
		 <<" #pm " <<fround(fduthetapeakm->GetParError(1),3)
		 <<"   b = "<<fround(fduthetapeakm->GetParameter(0),3)
		 <<" #pm "<<fround(fduthetapeakm->GetParError(0),3)
		 <<endl;

    duthetadecom->SetName("duthetadecom");
    TGraphErrors *gduthetadecom = getTGraphFromTH2(duthetadecom,thetabins,0);
    gduthetadecom->SetMarkerColor(4);
    gduthetadecom->SetLineColor(4);
    gduthetadecom->SetMarkerStyle(8);
    gduthetadecom->SetMarkerSize(size);
    TF1 *fduthetadecom=new TF1("fduthetadecom","pol1");
    fduthetadecom->SetLineColor(4);
    gduthetadecom->Fit(fduthetadecom);
    stringstream sduthetadecom;
    
    sduthetadecom<<"m = "<<fround(fduthetadecom->GetParameter(1),3)
		 <<" #pm " <<fround(fduthetadecom->GetParError(1),3)
		 <<"   b = "<<fround(fduthetadecom->GetParameter(0),3)
		 <<" #pm "<<fround(fduthetadecom->GetParError(0),3)
		 <<endl;
    
    //add deco with BP correction------------------------------------------
    if(addbpcorr){
      duthetadecobp_vm->SetName("duthetadecobp_vm");
      TGraphErrors *gduthetadecobp_vm = getTGraphFromTH2(duthetadecobp_vm,thetabins,0);
      gduthetadecobp_vm->SetMarkerColor(6);
      gduthetadecobp_vm->SetLineColor(6);
      gduthetadecobp_vm->SetMarkerStyle(8);
      gduthetadecobp_vm->SetMarkerSize(size);
      TF1 *fduthetadecobp_vm=new TF1("fduthetadecobp_vm","pol1");
      fduthetadecobp_vm->SetLineColor(6);
      gduthetadecobp_vm->Fit(fduthetadecobp_vm);
      stringstream sduthetadecobp_vm;
      
      sduthetadecobp_vm<<"m = "<<fround(fduthetadecobp_vm->GetParameter(1),3)
		       <<" #pm " <<fround(fduthetadecobp_vm->GetParError(1),3)
		       <<"   b = "<<fround(fduthetadecobp_vm->GetParameter(0),3)
		       <<" #pm "<<fround(fduthetadecobp_vm->GetParError(0),3)
		       <<endl;
    }

    TMultiGraph *mgdtm=new TMultiGraph();
    mgdtm->Add(gduthetapeakm);
    mgdtm->Add(gduthetadecom);
    if(addbpcorr) mgdtm->Add(gduthetadecobp_vm);
    mgdtm->SetTitle("TOB (v-)");
    mgdtm->Draw("AP");
    mgdtm->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
    mgdtm->GetYaxis()->SetTitle("<#Deltau> [#mum]");
    mgdtm->SetTitle("TOB");
    //mgdtm->GetXaxis()->SetLimits(-1,1);
    mgdtm->GetYaxis()->SetLimits(-20,20);
    mgdtm->SetMinimum(-20);
    mgdtm->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    //t->DrawLatex(0.5,0.8,"y = mx + b");
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.85,sduthetapeakm.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.75,sduthetadecom.str().c_str());
    if(addbpcorr){
      t->SetTextColor(6);
      t->DrawLatex(0.15,0.65,sduthetadecobp_vm.str().c_str());
    }

    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);

    idx++;
  }
  /*
  //TGraph du vs. trk theta----------------------------------------
  if(drawcan[12]){
    can[12]=new TCanvas(cantitles[12].c_str(),cantitles[12].c_str(),1200,450);
    can[12]->Divide(2,1);
  
    can[12]->cd(1);
    dutrkpeakp->SetName("dutrkpeakp");
    TGraphErrors *gdutrkpeakp = getTGraphFromTH2(dutrkpeakp,thetabins,0);
    gdutrkpeakp->SetMarkerColor(2);
    gdutrkpeakp->SetLineColor(2);
    gdutrkpeakp->SetMarkerStyle(8);
    gdutrkpeakp->SetMarkerSize(1);
    TF1 *fdutrkpeakp=new TF1("fdutrkpeakp","pol1");
    fdutrkpeakp->SetLineColor(2);
    gdutrkpeakp->Fit(fdutrkpeakp);
    stringstream sdutrkpeakp;
    sdutrkpeakp<<"#Delta_{W} = "<<fround(fdutrkpeakp->GetParameter(1),1)
	       <<" #epsilon = "<<fround(fdutrkpeakp->GetParameter(0),1)
	       <<" #Deltatan(#theta) = "<<fround(fdutrkpeakp->GetParameter(0)/fdutrkpeakp->GetParameter(1),2)<<endl;
    
    dutrkdecop->SetName("dutrkdecop");
    TGraphErrors *gdutrkdecop = getTGraphFromTH2(dutrkdecop,thetabins,0);
    gdutrkdecop->SetMarkerColor(4);
    gdutrkdecop->SetLineColor(4);
    gdutrkdecop->SetMarkerStyle(8);
    gdutrkdecop->SetMarkerSize(1);
    TF1 *fdutrkdecop=new TF1("fdutrkdecop","pol1");
    fdutrkdecop->SetLineColor(4);
    gdutrkdecop->Fit(fdutrkdecop);
    stringstream sdutrkdecop;
    sdutrkdecop<<"#Delta_{W} = "<<fround(fdutrkdecop->GetParameter(1),1)
	       <<" #epsilon = "<<fround(fdutrkdecop->GetParameter(0),1)
	       <<" #Deltatan(#theta) = "<<fround(fdutrkdecop->GetParameter(0)/fdutrkdecop->GetParameter(1),2)<<endl;
    
    TMultiGraph *mgdtp=new TMultiGraph();
    mgdtp->Add(gdutrkpeakp);
    mgdtp->Add(gdutrkdecop);
    mgdtp->SetTitle("TOB (v+)");
    mgdtp->Draw("AP");
    mgdtp->GetXaxis()->SetTitle("<tan(#theta_{trk})>");
    mgdtp->GetYaxis()->SetTitle("<#Delta_{u}> [#mum]");
    mgdtp->SetTitle("TOB");
    mgdtp->GetXaxis()->SetLimits(-1,1);
    mgdtp->SetMinimum(-20);
    mgdtp->SetMaximum(20);
 
    TLatex *t=new TLatex();
    t->SetNDC();
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.85,sdutrkpeakp.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.75,sdutrkdecop.str().c_str());
    
    line->DrawLine(-1,0,1,0);
    line->DrawLine(0,-20,0,20);
    
    can[12]->cd(2);
    dutrkpeakm->SetName("dutrkpeakm");
    TGraphErrors *gdutrkpeakm = getTGraphFromTH2(dutrkpeakm,thetabins,0);
    gdutrkpeakm->SetMarkerColor(2);
    gdutrkpeakm->SetLineColor(2);
    gdutrkpeakm->SetMarkerStyle(8);
    gdutrkpeakm->SetMarkerSize(1);
    TF1 *fdutrkpeakm=new TF1("fdutrkpeakm","pol1");
    fdutrkpeakm->SetLineColor(2);
    gdutrkpeakm->Fit(fdutrkpeakm);
    stringstream sdutrkpeakm;
    sdutrkpeakm<<"#Delta_{W} = "<<fround(fdutrkpeakm->GetParameter(1),1)
	       <<" #epsilon = "<<fround(fdutrkpeakm->GetParameter(0),1)
	       <<" #Deltatan(#theta) = "<<fround(fdutrkpeakm->GetParameter(0)/fdutrkpeakm->GetParameter(1),2)<<endl;
    
    dutrkdecom->SetName("dutrkdecom");
    TGraphErrors *gdutrkdecom = getTGraphFromTH2(dutrkdecom,thetabins,0);
    gdutrkdecom->SetMarkerColor(4);
    gdutrkdecom->SetLineColor(4);
    gdutrkdecom->SetMarkerStyle(8);
    gdutrkdecom->SetMarkerSize(1);
    TF1 *fdutrkdecom=new TF1("fdutrkdecom","pol1");
    fdutrkdecom->SetLineColor(4);
    gdutrkdecom->Fit(fdutrkdecom);
    stringstream sdutrkdecom;
    sdutrkdecom<<"#Delta_{W} = "<<fround(fdutrkdecom->GetParameter(1),1)
	       <<" #epsilon = "<<fround(fdutrkdecom->GetParameter(0),1)
	       <<" #Deltatan(#theta) = "<<fround(fdutrkdecom->GetParameter(0)/fdutrkdecom->GetParameter(1),2)<<endl;
    
    
    TMultiGraph *mgdtm=new TMultiGraph();
    mgdtm->Add(gdutrkpeakm);
    mgdtm->Add(gdutrkdecom);
    mgdtm->SetTitle("TOB (v-)");
    mgdtm->Draw("AP");
    mgdtm->GetXaxis()->SetTitle("<tan(#theta_{trk})>");
    mgdtm->GetYaxis()->SetTitle("<#Delta_{u}> [#mum]");
    mgdtm->SetTitle("TOB");
    mgdtm->GetXaxis()->SetLimits(-1,1);
    mgdtm->SetMinimum(-20);
    mgdtm->SetMaximum(20);
  
    TLatex *t=new TLatex();
    t->SetNDC();
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.75,sdutrkpeakm.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.65,sdutrkdecom.str().c_str());
    
    line->DrawLine(-1,0,1,0);
    line->DrawLine(0,-20,0,20);
  }


  //TGraph du vs. delta tan theta----------------------------------------
  if(drawcan[13]){
    can[13]=new TCanvas(cantitles[13].c_str(),cantitles[13].c_str(),600,450);
    can[13]->cd();
    
    vector<float> thetabins;
//     thetabins.push_back(-1.);
//     thetabins.push_back(-0.75);
//     thetabins.push_back(-0.5);
//     thetabins.push_back(-0.25);
//     thetabins.push_back(0.);
//     thetabins.push_back(0.25);
//     thetabins.push_back(0.5);
//     thetabins.push_back(0.75);
//     thetabins.push_back(1.);

    thetabins.push_back(-0.9);
    thetabins.push_back(-0.7);
    thetabins.push_back(-0.5);
    thetabins.push_back(-0.3);
    thetabins.push_back(-0.1);
    thetabins.push_back(0.1);
    thetabins.push_back(0.3);
    thetabins.push_back(0.5);
    thetabins.push_back(0.7);
    thetabins.push_back(0.9);
    
    duthetapeak->SetName("duthetapeak");
    TGraphErrors *gduthetapeak = getTGraphFromTH2(duthetapeak,thetabins,0);
    gduthetapeak->SetMarkerColor(2);
    gduthetapeak->SetLineColor(2);
    gduthetapeak->SetMarkerStyle(8);
    gduthetapeak->SetMarkerSize(1);
    TF1 *fduthetapeak=new TF1("fduthetapeak","pol1");
    fduthetapeak->SetLineColor(2);
    gduthetapeak->Fit(fduthetapeak);
    stringstream sduthetapeak;
    sduthetapeak<<"Slope = #Deltaw = "<<fround(fdtppeak->GetParameter(1),1)<<" #mum"<<endl;
    
    duthetadeco->SetName("duthetadeco");
    TGraphErrors *gduthetadeco = getTGraphFromTH2(duthetadeco,thetabins,0);
    gduthetadeco->SetMarkerColor(4);
    gduthetadeco->SetLineColor(4);
    gduthetadeco->SetMarkerStyle(8);
    gduthetadeco->SetMarkerSize(1);
    TF1 *fduthetadeco=new TF1("fduthetadeco","pol1");
    fduthetadeco->SetLineColor(4);
    gduthetadeco->Fit(fduthetadeco);
    stringstream sduthetadeco;
    sduthetadeco<<"Slope = #Deltaw = "<<fround(fduthetadeco->GetParameter(1),1)<<" #mum"<<endl;
      
    TMultiGraph *mgdtp=new TMultiGraph();
    mgdtp->Add(gduthetapeak);
    mgdtp->Add(gduthetadeco);
    mgdtp->SetTitle("TOB");
    mgdtp->Draw("AP");
    mgdtp->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
    mgdtp->GetYaxis()->SetTitle("<#Deltau> [#mum]");
    mgdtp->SetTitle("TOB");
    mgdtp->GetXaxis()->SetLimits(-1,1);
    mgdtp->GetYaxis()->SetLimits(-20,20);
    mgdtp->SetMinimum(-20);
    mgdtp->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.75,sduthetapeak.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.65,sduthetadeco.str().c_str());
    
    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);
    
  
  }

  if(drawcan[14]){
    can[14]=new TCanvas(cantitles[14].c_str(),cantitles[14].c_str(),1200,900);
    can[14]->Divide(2,2);

    plotHistsTH2(duthetapeakp,
		 "TOB PEAK (v+)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[14],1);
    plotHistsTH2(duthetapeakm,
		 "TOB PEAK (v-)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[14],2);
    
    can[14]->cd(3);
    TH1F* duthetapeakpproj        = (TH1F*) duthetapeakp->ProjectionY();

    duthetapeakpproj->SetLineWidth(3);
    duthetapeakpproj->Draw();
    duthetapeakpproj->SetTitle("y-projection");
    duthetapeakpproj->GetXaxis()->SetTitle("#Deltau [#mum]");

    
    TH1F* duthetapeakpproj        = (TH1F*)duthetapeakp->ProjectionX();
    //TH1F* duthetapeakpproj10      = (TH1F*)duthetapeakp->ProjectionX("",99,102)->Clone();
    //TH1F* duthetapeakpproj25      = (TH1F*)duthetapeakp->ProjectionX("",96,98)->Clone();
    //duthetapeakpproj25            -> Add((TH1F*)duthetapeakp->ProjectionX("",103,105)->Clone());
    TH1F* duthetapeakpproj25      = (TH1F*)duthetapeakp->ProjectionX("",96,105)->Clone();
    TH1F* duthetapeakpproj50      = (TH1F*)duthetapeakp->ProjectionX("",91,95)->Clone();
    duthetapeakpproj50            -> Add((TH1F*)duthetapeakp->ProjectionX("",106,110)->Clone());
    TH1F* duthetapeakpproj100     = (TH1F*)duthetapeakp->ProjectionX("",81,90)->Clone();
    duthetapeakpproj100           -> Add((TH1F*)duthetapeakp->ProjectionX("",111,120)->Clone());
    TH1F* duthetapeakpproj500     = (TH1F*)duthetapeakp->ProjectionX("",0,80)->Clone();
    duthetapeakpproj500           -> Add((TH1F*)duthetapeakp->ProjectionX("",121,201)->Clone());
    
    duthetapeakpproj->SetLineWidth(3);
    duthetapeakpproj->Draw();
    duthetapeakpproj->SetTitle("x-projection");
    duthetapeakpproj->GetXaxis()->SetTitle("tan(#theta_{trk})-tan(#theta_{L})");

    //duthetapeakpproj10->SetLineColor(2);
    duthetapeakpproj25->SetLineColor(4);
    duthetapeakpproj50->SetLineColor(7);
    duthetapeakpproj100->SetLineColor(6);
    duthetapeakpproj500->SetLineColor(2);

    //duthetapeakpproj10->SetMarkerColor(2);
    duthetapeakpproj25->SetMarkerColor(4);
    duthetapeakpproj50->SetMarkerColor(7);
    duthetapeakpproj100->SetMarkerColor(6);
    duthetapeakpproj500->SetMarkerColor(2);

    //duthetapeakpproj10 ->Draw("same");
    duthetapeakpproj25 ->Draw("same");
    duthetapeakpproj50 ->Draw("same");
    duthetapeakpproj100->Draw("same");
    duthetapeakpproj500->Draw("same");


    TLegend *legproj=new TLegend(0.65,0.5,0.9,0.9);
    legproj->AddEntry(duthetapeakpproj,   "total");
    //legproj->AddEntry(duthetapeakpproj500,"|#Deltau| > 100 #mum","f");
    //legproj->AddEntry(duthetapeakpproj100,"#Deltau| = 50-100 #mum","f");
    //legproj->AddEntry(duthetapeakpproj50, "|#Deltau| = 10-50 #mum","f");   
    //legproj->AddEntry(duthetapeakpproj10, "|#Deltau| < 10 #mum","f");
    legproj->AddEntry(duthetapeakpproj500,"|#Deltau| > 100 #mum");
    legproj->AddEntry(duthetapeakpproj100,"#Deltau| = 50-100 #mum");
    legproj->AddEntry(duthetapeakpproj50, "|#Deltau| = 25-50 #mum");   
    legproj->AddEntry(duthetapeakpproj25, "|#Deltau| < 25 #mum");   
    //legproj->AddEntry(duthetapeakpproj25, "|#Deltau| = 10-25 #mum");   
    //legproj->AddEntry(duthetapeakpproj10, "|#Deltau| < 10 #mum");

    legproj->SetBorderSize(1);
    legproj->SetFillColor(0);
    legproj->Draw();
    
    


    //     THStack* stackp=new THStack("stackp","stackp");
    //     duthetapeakpproj10->SetFillColor(2);
    //     duthetapeakpproj50->SetFillColor(4);
    //     duthetapeakpproj100->SetFillColor(6);
    //     duthetapeakpproj500->SetFillColor(5);

    //     stackp->Add( duthetapeakpproj10 );
    //     stackp->Add( duthetapeakpproj50 );
    //     stackp->Add( duthetapeakpproj100 ); 
    //     stackp->Add( duthetapeakpproj500 );
    
    //    stackp->Draw("same");
    //    duthetapeakpproj->Draw("same");

    


    can[14]->cd(4);
    
    TH1F* duthetapeakmproj        = (TH1F*)duthetapeakm->ProjectionX();
    //TH1F* duthetapeakmproj10      = (TH1F*)duthetapeakm->ProjectionX("",99,102)->Clone();
    //TH1F* duthetapeakmproj25      = (TH1F*)duthetapeakm->ProjectionX("",96,98)->Clone();
    //duthetapeakmproj25            -> Add((TH1F*)duthetapeakm->ProjectionX("",103,105)->Clone());
    TH1F* duthetapeakmproj25      = (TH1F*)duthetapeakm->ProjectionX("",96,105)->Clone();
    TH1F* duthetapeakmproj50      = (TH1F*)duthetapeakm->ProjectionX("",91,95)->Clone();
    duthetapeakmproj50            -> Add((TH1F*)duthetapeakm->ProjectionX("",106,110)->Clone());
    TH1F* duthetapeakmproj100     = (TH1F*)duthetapeakm->ProjectionX("",81,90)->Clone();
    duthetapeakmproj100           -> Add((TH1F*)duthetapeakm->ProjectionX("",111,120)->Clone());
    TH1F* duthetapeakmproj500     = (TH1F*)duthetapeakm->ProjectionX("",0,80)->Clone();
    duthetapeakmproj500           -> Add((TH1F*)duthetapeakm->ProjectionX("",121,201)->Clone());
    
    duthetapeakmproj->SetLineWidth(3);
    duthetapeakmproj->Draw();
    duthetapeakmproj->SetTitle("x-projection");
    duthetapeakmproj->GetXaxis()->SetTitle("tan(#theta_{trk})-tan(#theta_{L})");

    //duthetapeakmproj10->SetLineColor(2);
    duthetapeakmproj25->SetLineColor(4);
    duthetapeakmproj50->SetLineColor(7);
    duthetapeakmproj100->SetLineColor(6);
    duthetapeakmproj500->SetLineColor(2);

    //duthetapeakmproj10->SetMarkerColor(2);
    duthetapeakmproj25->SetMarkerColor(4);
    duthetapeakmproj50->SetMarkerColor(7);
    duthetapeakmproj100->SetMarkerColor(6);
    duthetapeakmproj500->SetMarkerColor(2);

    //duthetapeakmproj10 ->Draw("same");
    duthetapeakmproj25 ->Draw("same");
    duthetapeakmproj50 ->Draw("same");
    duthetapeakmproj100->Draw("same");
    duthetapeakmproj500->Draw("same");
    
    //     THStack* stackm=new THStack("stackm","stackm");
    //     duthetapeakmproj10->SetFillColor(2);
    //     duthetapeakmproj50->SetFillColor(4);
    //     duthetapeakmproj100->SetFillColor(6);
    //     duthetapeakmproj500->SetFillColor(5);

    //     stackm->Add( duthetapeakmproj10 );
    //     stackm->Add( duthetapeakmproj50 );
    //     stackm->Add( duthetapeakmproj100 ); 
    //     stackm->Add( duthetapeakmproj500 );
    
    //    stackm->Draw("same");
    //    duthetapeakmproj->Draw("same");
    //legproj->Draw();
  }

  if(drawcan[15]){
    can[15]=new TCanvas(cantitles[15].c_str(),cantitles[15].c_str(),1200,900);
    can[15]->Divide(2,2);

    plotHistsTH2(dutrkpeakp,
		 "TOB PEAK (v+)","tan(#theta_{trk})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[15],1);
    plotHistsTH2(dutrkpeakm,
		 "TOB PEAK (v-)","tan(#theta_{trk})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[15],2);
    
    can[15]->cd(3);

    TH1F* dutrkpeakpproj        = (TH1F*)dutrkpeakp->ProjectionX();
    //TH1F* dutrkpeakpproj10      = (TH1F*)dutrkpeakp->ProjectionX("",99,102)->Clone();
    //TH1F* dutrkpeakpproj25      = (TH1F*)dutrkpeakp->ProjectionX("",96,98)->Clone();
    //dutrkpeakpproj25            -> Add((TH1F*)dutrkpeakp->ProjectionX("",103,105)->Clone());
    TH1F* dutrkpeakpproj25      = (TH1F*)dutrkpeakp->ProjectionX("",96,105)->Clone();
    TH1F* dutrkpeakpproj50      = (TH1F*)dutrkpeakp->ProjectionX("",91,95)->Clone();
    dutrkpeakpproj50            -> Add((TH1F*)dutrkpeakp->ProjectionX("",106,110)->Clone());
    TH1F* dutrkpeakpproj100     = (TH1F*)dutrkpeakp->ProjectionX("",81,90)->Clone();
    dutrkpeakpproj100           -> Add((TH1F*)dutrkpeakp->ProjectionX("",111,120)->Clone());
    TH1F* dutrkpeakpproj500     = (TH1F*)dutrkpeakp->ProjectionX("",0,80)->Clone();
    dutrkpeakpproj500           -> Add((TH1F*)dutrkpeakp->ProjectionX("",121,201)->Clone());
    
    dutrkpeakpproj->SetLineWidth(3);
    dutrkpeakpproj->SetMaximum(1.05*dutrkpeakpproj->GetMaximum());
    dutrkpeakpproj->Draw();
    dutrkpeakpproj->SetTitle("x-projection");
    dutrkpeakpproj->GetXaxis()->SetTitle("tan(#theta_{trk})");

    //dutrkpeakpproj10->SetLineColor(2);
    dutrkpeakpproj25->SetLineColor(4);
    dutrkpeakpproj50->SetLineColor(7);
    dutrkpeakpproj100->SetLineColor(6);
    dutrkpeakpproj500->SetLineColor(2);

    //dutrkpeakpproj10->SetMarkerColor(2);
    dutrkpeakpproj25->SetMarkerColor(4);
    dutrkpeakpproj50->SetMarkerColor(7);
    dutrkpeakpproj100->SetMarkerColor(6);
    dutrkpeakpproj500->SetMarkerColor(2);

    //dutrkpeakpproj10 ->Draw("same");
    dutrkpeakpproj25 ->Draw("same");
    dutrkpeakpproj50 ->Draw("same");
    dutrkpeakpproj100->Draw("same");
    dutrkpeakpproj500->Draw("same");

    TLine *line=new TLine();
    line->SetLineStyle(2);
    line->DrawLine(0.,0.,0.,dutrkpeakpproj->GetMaximum());
    line->DrawLine(-0.1,0.,-0.1,dutrkpeakpproj->GetMaximum());

    //     THStack* stackp=new THStack("stackp","stackp");
    //     dutrkpeakpproj10->SetFillColor(2);
    //     dutrkpeakpproj50->SetFillColor(4);
    //     dutrkpeakpproj100->SetFillColor(6);
    //     dutrkpeakpproj500->SetFillColor(5);

    //     stackp->Add( dutrkpeakpproj10 );
    //     stackp->Add( dutrkpeakpproj50 );
    //     stackp->Add( dutrkpeakpproj100 ); 
    //     stackp->Add( dutrkpeakpproj500 );
    
    //    stackp->Draw("same");
    //    dutrkpeakpproj->Draw("same");

    TLegend *legproj=new TLegend(0.65,0.5,0.9,0.9);
    legproj->AddEntry(dutrkpeakpproj,   "total");
    //legproj->AddEntry(dutrkpeakpproj500,"|#Deltau| > 100 #mum","f");
    //legproj->AddEntry(dutrkpeakpproj100,"#Deltau| = 50-100 #mum","f");
    //legproj->AddEntry(dutrkpeakpproj50, "|#Deltau| = 10-50 #mum","f");   
    //legproj->AddEntry(dutrkpeakpproj10, "|#Deltau| < 10 #mum","f");
    legproj->AddEntry(dutrkpeakpproj500,"|#Deltau| > 100 #mum");
    legproj->AddEntry(dutrkpeakpproj100,"#Deltau| = 50-100 #mum");
    legproj->AddEntry(dutrkpeakpproj50, "|#Deltau| = 25-50 #mum");   
    legproj->AddEntry(dutrkpeakpproj25, "|#Deltau| < 25 #mum");   
    //legproj->AddEntry(dutrkpeakpproj25, "|#Deltau| = 10-25 #mum");   
    //legproj->AddEntry(dutrkpeakpproj10, "|#Deltau| < 10 #mum");
    
    legproj->SetBorderSize(1);
    legproj->SetFillColor(0);
    legproj->Draw();

    can[15]->cd(4);

    TH1F* dutrkpeakmproj        = (TH1F*)dutrkpeakm->ProjectionX();
    //TH1F* dutrkpeakmproj10      = (TH1F*)dutrkpeakm->ProjectionX("",99,102)->Clone();
    //TH1F* dutrkpeakmproj25      = (TH1F*)dutrkpeakm->ProjectionX("",96,98)->Clone();
    //dutrkpeakmproj25            -> Add((TH1F*)dutrkpeakm->ProjectionX("",103,105)->Clone());
    TH1F* dutrkpeakmproj25      = (TH1F*)dutrkpeakm->ProjectionX("",96,105)->Clone();
    TH1F* dutrkpeakmproj50      = (TH1F*)dutrkpeakm->ProjectionX("",91,95)->Clone();
    dutrkpeakmproj50            -> Add((TH1F*)dutrkpeakm->ProjectionX("",106,110)->Clone());
    TH1F* dutrkpeakmproj100     = (TH1F*)dutrkpeakm->ProjectionX("",81,90)->Clone();
    dutrkpeakmproj100           -> Add((TH1F*)dutrkpeakm->ProjectionX("",111,120)->Clone());
    TH1F* dutrkpeakmproj500     = (TH1F*)dutrkpeakm->ProjectionX("",0,80)->Clone();
    dutrkpeakmproj500           -> Add((TH1F*)dutrkpeakm->ProjectionX("",121,201)->Clone());
    
    dutrkpeakmproj->SetLineWidth(3);
    dutrkpeakmproj->SetMaximum(1.05*dutrkpeakmproj->GetMaximum());
    dutrkpeakmproj->Draw();
    dutrkpeakmproj->SetTitle("x-projection");
    dutrkpeakmproj->GetXaxis()->SetTitle("tan(#theta_{trk})");

    //dutrkpeakmproj10->SetLineColor(2);
    dutrkpeakmproj25->SetLineColor(4);
    dutrkpeakmproj50->SetLineColor(7);
    dutrkpeakmproj100->SetLineColor(6);
    dutrkpeakmproj500->SetLineColor(2);

    //dutrkpeakmproj10->SetMarkerColor(2);
    dutrkpeakmproj25->SetMarkerColor(4);
    dutrkpeakmproj50->SetMarkerColor(7);
    dutrkpeakmproj100->SetMarkerColor(6);
    dutrkpeakmproj500->SetMarkerColor(2);

    //dutrkpeakmproj10 ->Draw("same");
    dutrkpeakmproj25 ->Draw("same");
    dutrkpeakmproj50 ->Draw("same");
    dutrkpeakmproj100->Draw("same");
    dutrkpeakmproj500->Draw("same");

    line->DrawLine(0.,0.,0.,dutrkpeakmproj->GetMaximum());
    line->DrawLine(0.1,0.,0.1,dutrkpeakmproj->GetMaximum());
    //     THStack* stackm=new THStack("stackm","stackm");
    //     dutrkpeakmproj10->SetFillColor(2);
    //     dutrkpeakmproj50->SetFillColor(4);
    //     dutrkpeakmproj100->SetFillColor(6);
    //     dutrkpeakmproj500->SetFillColor(5);

    //     stackm->Add( dutrkpeakmproj10 );
    //     stackm->Add( dutrkpeakmproj50 );
    //     stackm->Add( dutrkpeakmproj100 ); 
    //     stackm->Add( dutrkpeakmproj500 );
    
    //    stackm->Draw("same");
    //    dutrkpeakmproj->Draw("same");
    legproj->Draw();
  }
 
  //du vs. delta tan(theta) split v+/v- histos (TEC)
  if(drawcan[16]){
    can[16]=new TCanvas(cantitles[16].c_str(),cantitles[16].c_str(),1200,900);
    can[16]->Divide(2,2);
    
    plotHistsTH2(duthetapeak_TECthick,
		 "TEC thick PEAK","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[16],1);
    plotHistsTH2(duthetadeco_TECthick,
		 "TEC thick DECO","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[16],2);
    plotHistsTH2(duthetapeak_TECthin,
 		 "TEC thin PEAK","tan(#theta_{trk})-tan(#theta_{L})",
 		 "#Delta u [#mum]",-2,2,-200,200,0,can[16],3);
    plotHistsTH2(duthetadeco_TECthin,
 		 "TEC thin DECO","tan(#theta_{trk})-tan(#theta_{L})",
 		 "#Delta u [#mum]",-2,2,-200,200,0,can[16],4);
  }

  if(drawcan[17]){
    can[17]=new TCanvas(cantitles[17].c_str(),cantitles[17].c_str(),1200,900);
    can[17]->Divide(2,2);

    plotHistsTH2(duthetapeakp,
		 "TOB PEAK (v+)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[17],1);
    plotHistsTH2(duthetapeakm,
		 "TOB PEAK (v-)","tan(#theta_{trk})-tan(#theta_{L})",
		 "#Delta u [#mum]",-2,2,-200,200,0,can[17],2);
    
    can[17]->cd(3);
    TH1F* duthetapeakpproj        = (TH1F*) duthetapeakp->ProjectionY()->Clone();
    TH1F* duthetapeakpproj1       = (TH1F*) duthetapeakp->ProjectionY("p1",49,52)->Clone();
    TH1F* duthetapeakpproj2       = (TH1F*) duthetapeakp->ProjectionY("p2",56,59)->Clone();
    TH1F* duthetapeakpproj3       = (TH1F*) duthetapeakp->ProjectionY("p3",42,45)->Clone();

    duthetapeakpproj->Scale (1./duthetapeakpproj ->Integral());
    duthetapeakpproj1->Scale(1./duthetapeakpproj1->Integral());
    duthetapeakpproj2->Scale(1./duthetapeakpproj2->Integral());
    duthetapeakpproj3->Scale(1./duthetapeakpproj3->Integral());

    stringstream sp1;
    stringstream sp2;
    stringstream sp3;

    //     TF1 *fp1=new TF1("fp1","gaus",-200,200);
    //     TF1 *fp2=new TF1("fp2","gaus",-200,200);
    //     TF1 *fp3=new TF1("fp3","gaus",-200,200);

    //     fp2->SetLineColor(2);
    //     fp3->SetLineColor(4);
   
    //     duthetapeakpproj1->Fit(fp1);
    //     duthetapeakpproj2->Fit(fp2);
    //     duthetapeakpproj3->Fit(fp3);

    //     sp1<<"#sigma = "<<fround(fp1->GetParameter(2),1)<<" #mum"<<endl;
    //     sp2<<"#sigma = "<<fround(fp2->GetParameter(2),1)<<" #mum"<<endl;
    //     sp3<<"#sigma = "<<fround(fp3->GetParameter(2),1)<<" #mum"<<endl;

    duthetapeakpproj1->Draw();
    duthetapeakpproj1->GetXaxis()->SetRangeUser(-200,200);
    duthetapeakpproj2->GetXaxis()->SetRangeUser(-200,200);
    duthetapeakpproj3->GetXaxis()->SetRangeUser(-200,200);
    duthetapeakpproj1->SetTitle("y-projection");
    duthetapeakpproj1->GetXaxis()->SetTitle("#Deltau [#mum]");

    duthetapeakpproj2->SetLineColor(2);
    duthetapeakpproj3->SetLineColor(4);
    duthetapeakpproj2->Draw("same");
    duthetapeakpproj3->Draw("same");
    duthetapeakpproj1->SetMaximum(1.05*duthetapeakpproj2->GetMaximum());
    
    sp1<<"RMS = "<<fround(duthetapeakpproj1->GetRMS(1),1)<<" #mum"<<endl;
    sp2<<"RMS = "<<fround(duthetapeakpproj2->GetRMS(1),1)<<" #mum"<<endl;
    sp3<<"RMS = "<<fround(duthetapeakpproj3->GetRMS(1),1)<<" #mum"<<endl;

    TLatex *t=new TLatex();
    t->SetNDC();
    t->DrawLatex(0.15,0.8,sp1.str().c_str());
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.7,sp2.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.6,sp3.str().c_str());

    TLegend *leg17=new TLegend(0.6,0.6,0.9,0.9);
    leg17->AddEntry(duthetapeakpproj1,"-0.08 < #Deltatan(#theta) <  0.08");
    leg17->AddEntry(duthetapeakpproj2," 0.20 < #Deltatan(#theta) <  0.36");
    leg17->AddEntry(duthetapeakpproj3,"-0.36 < #Deltatan(#theta) < -0.20");
    leg17->SetFillColor(0);
    leg17->SetBorderSize(1);
    leg17->Draw();


    can[17]->cd(4);
    TH1F* duthetapeakmproj        = (TH1F*) duthetapeakm->ProjectionY()->Clone();
    TH1F* duthetapeakmproj1       = (TH1F*) duthetapeakm->ProjectionY("p1",49,52)->Clone();
    TH1F* duthetapeakmproj2       = (TH1F*) duthetapeakm->ProjectionY("p2",56,59)->Clone();
    TH1F* duthetapeakmproj3       = (TH1F*) duthetapeakm->ProjectionY("p3",42,45)->Clone();

    duthetapeakmproj->Scale (1./duthetapeakmproj ->Integral());
    duthetapeakmproj1->Scale(1./duthetapeakmproj1->Integral());
    duthetapeakmproj2->Scale(1./duthetapeakmproj2->Integral());
    duthetapeakmproj3->Scale(1./duthetapeakmproj3->Integral());

    stringstream sm1;
    stringstream sm2;
    stringstream sm3;

    //     TF1 *fp1=new TF1("fp1","gaus",-200,200);
    //     TF1 *fp2=new TF1("fp2","gaus",-200,200);
    //     TF1 *fp3=new TF1("fp3","gaus",-200,200);

    //     fp2->SetLineColor(2);
    //     fp3->SetLineColor(4);
   
    //     duthetapeakmproj1->Fit(fp1);
    //     duthetapeakmproj2->Fit(fp2);
    //     duthetapeakmproj3->Fit(fp3);

    //     sm1<<"#sigma = "<<fround(fp1->GetParameter(2),1)<<" #mum"<<endl;
    //     sm2<<"#sigma = "<<fround(fp2->GetParameter(2),1)<<" #mum"<<endl;
    //     sm3<<"#sigma = "<<fround(fp3->GetParameter(2),1)<<" #mum"<<endl;

    duthetapeakmproj1->Draw();
    duthetapeakmproj1->GetXaxis()->SetRangeUser(-200,200);
    duthetapeakmproj2->GetXaxis()->SetRangeUser(-200,200);
    duthetapeakmproj3->GetXaxis()->SetRangeUser(-200,200);
    duthetapeakmproj1->SetTitle("y-projection");
    duthetapeakmproj1->GetXaxis()->SetTitle("#Deltau [#mum]");

    duthetapeakmproj2->SetLineColor(2);
    duthetapeakmproj3->SetLineColor(4);
    duthetapeakmproj2->Draw("same");
    duthetapeakmproj3->Draw("same");
    duthetapeakmproj1->SetMaximum(1.05*duthetapeakmproj3->GetMaximum());
    
    sm1<<"RMS = "<<fround(duthetapeakmproj1->GetRMS(1),1)<<" #mum"<<endl;
    sm2<<"RMS = "<<fround(duthetapeakmproj2->GetRMS(1),1)<<" #mum"<<endl;
    sm3<<"RMS = "<<fround(duthetapeakmproj3->GetRMS(1),1)<<" #mum"<<endl;

    t->SetNDC();
    t->SetTextColor(1);
    t->DrawLatex(0.15,0.8,sm1.str().c_str());
    t->SetTextColor(2);
    t->DrawLatex(0.15,0.7,sm2.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.15,0.6,sm3.str().c_str());
    
    leg17->Draw();

  }

  if(drawcan[18]){
    can[18]=new TCanvas(cantitles[18].c_str(),cantitles[18].c_str(),1200,900);
    can[18]->Divide(2,2);
    
    TH2F* hp[8];
    fpeak.cd();
    PeakDecoResiduals.cd();
    PeakDecoResiduals.cd();
    hp[0]=(TH2F*) duvsdtantheta_vp_dt_TOB_1->Clone();

    TH2F* hm[8];


    can[18]->cd(1);
    hp[0]->Draw("colz");
    

  }

  if(drawcan[19]){
    can[19]=new TCanvas(cantitles[19].c_str(),cantitles[19].c_str(),1200,900);
    can[19]->Divide(2,2);
  
    can[19]->cd(1);
    setStats(dwpeak_TECthick,dwdeco_TECthick,0.5,0.65,0.15,false);  
    plotHists(dwpeak_TECthick,dwdeco_TECthick,"TEC thick","#Delta w [#mum]",-1000,1000,1,3);
  
    can[19]->cd(2);
    setStats(dwpeak_TECthin,dwdeco_TECthin,0.5,0.65,0.15,false);  
    plotHists(dwpeak_TECthin,dwdeco_TECthin,"TEC thin","#Delta w [#mum]",-1000,1000,1,3);

    can[19]->cd(3);
    setStats(dw2peak_TECthick,dw2deco_TECthick,0.5,0.65,0.15,false);  
    plotHists(dw2peak_TECthick,dw2deco_TECthick,"TEC thick","#Delta w [#mum]",-1000,1000,1,3);
  
    can[19]->cd(4);
    setStats(dw2peak_TECthin,dw2deco_TECthin,0.5,0.65,0.15,false);  
    plotHists(dw2peak_TECthin,dw2deco_TECthin,"TEC thin","#Delta w [#mum]",-1000,1000,1,3);
  }
  */
 
  
  if(draw("dtanladt",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
   
    //    can[20]->cd(1);
    //    h_vp_peak[0]->Draw("colz");
    //    can[idx]->cd(2);
    //    h_vm_peak[0]->Draw("colz");
    //    can[idx]->cd(3);
    //    h_vp_deco[0]->Draw("colz");
    //    can[idx]->cd(4);
    //    h_vm_deco[0]->Draw("colz");

    vector<float> thetabins;
    thetabins.push_back(-1.);
    thetabins.push_back(-0.75);
    thetabins.push_back(-0.5);
    thetabins.push_back(-0.25);
    thetabins.push_back(0.);
    thetabins.push_back(0.25);
    thetabins.push_back(0.5);
    thetabins.push_back(0.75);
    thetabins.push_back(1.);

    TGraphErrors* g_vp_peak[8];
    TGraphErrors* g_vm_peak[8];
    TGraphErrors* g_vp_deco[8];
    TGraphErrors* g_vm_deco[8];
    stringstream  s_vp_peak[8];
    stringstream  s_vm_peak[8];
    stringstream  s_vp_deco[8];
    stringstream  s_vm_deco[8];
    TF1*          f_vp_peak[8];
    TF1*          f_vm_peak[8];
    TF1*          f_vp_deco[8];
    TF1*          f_vm_deco[8];
    float         m_vp_peak[8];
    float         m_vm_peak[8];
    float         m_vp_deco[8];
    float         m_vm_deco[8];
    float         b_vp_peak[8];
    float         b_vm_peak[8];
    float         b_vp_deco[8];
    float         b_vm_deco[8];
    float         merr_vp_peak[8];
    float         merr_vm_peak[8];
    float         merr_vp_deco[8];
    float         merr_vm_deco[8];
    float         berr_vp_peak[8];
    float         berr_vm_peak[8];
    float         berr_vp_deco[8];
    float         berr_vm_deco[8];
    float         dtanla_vp_peak[8];
    float         dtanla_vm_peak[8];
    float         dtanla_vp_deco[8];
    float         dtanla_vm_deco[8];
    float         dtanlaerr_vp_peak[8];
    float         dtanlaerr_vm_peak[8];
    float         dtanlaerr_vp_deco[8];
    float         dtanlaerr_vm_deco[8];

    float xlim[]={-20,-12,-8,-4,0,4,8,12,20};
   
    float xpeak[8];
    float xerrpeak[8];
    float xdeco[8];
    float xerrdeco[8];

    for(int i=0;i<8;i++){
      dttimepeak->GetXaxis()->SetRangeUser(xlim[i],xlim[i+1]);
      dttimedeco->GetXaxis()->SetRangeUser(xlim[i],xlim[i+1]);

      xpeak[i]   =dttimepeak->GetMean(1);
      xerrpeak[i]=dttimepeak->GetRMS(1);
      xdeco[i]   =dttimedeco->GetMean(1);
      xerrdeco[i]=dttimedeco->GetRMS(1);

      //cout<<dttimepeak->GetMean(1)<<" "<<dttimepeak->GetRMS(1)<<endl;

    }
   
    for(int i=0;i<8;i++){
   
      //duthetapeakp->SetName("duthetapeakp");
      g_vp_peak[i] = getTGraphFromTH2(h_vp_peak[i],thetabins,0);
      g_vp_peak[i]->SetMarkerColor(2);
      g_vp_peak[i]->SetLineColor(2);
      g_vp_peak[i]->SetMarkerStyle(8);
      g_vp_peak[i]->SetMarkerSize(1);
      f_vp_peak[i]=new TF1(Form("f_vp_peak_%i",i),"pol1");
      f_vp_peak[i]->SetLineColor(2);
      g_vp_peak[i]->Fit(f_vp_peak[i]);
      s_vp_peak[i]<<"m = "<<fround(f_vp_peak[i]->GetParameter(1),1)
		  <<"   b = "<<fround(f_vp_peak[i]->GetParameter(0),1)<<endl;
     
      b_vp_peak[i]         = f_vp_peak[i]->GetParameter(0);
      m_vp_peak[i]         = f_vp_peak[i]->GetParameter(1);
      berr_vp_peak[i]      = f_vp_peak[i]->GetParError(0);
      merr_vp_peak[i]      = f_vp_peak[i]->GetParError(1);
      dtanla_vp_peak[i]    = b_vp_peak[i]/(250.-m_vp_peak[i]);
      dtanlaerr_vp_peak[i] = berr_vp_peak[i]/(250.-m_vp_peak[i]);
     
      //can[idx]->cd(1);
      //g_vp_peak[i]->Draw("AP");

      g_vm_peak[i] = getTGraphFromTH2(h_vm_peak[i],thetabins,0,1);
      g_vm_peak[i]->SetMarkerColor(2);
      g_vm_peak[i]->SetLineColor(2);
      g_vm_peak[i]->SetMarkerStyle(8);
      g_vm_peak[i]->SetMarkerSize(1);
      f_vm_peak[i]=new TF1(Form("f_vm_peak_%i",i),"pol1");
      f_vm_peak[i]->SetLineColor(2);
      g_vm_peak[i]->Fit(f_vm_peak[i]);
      s_vm_peak[i]<<"m = "<<fround(f_vm_peak[i]->GetParameter(1),1)
		  <<"   b = "<<fround(f_vm_peak[i]->GetParameter(0),1)<<endl;
     
      b_vm_peak[i]         = f_vm_peak[i]->GetParameter(0);
      m_vm_peak[i]         = f_vm_peak[i]->GetParameter(1);
      berr_vm_peak[i]      = f_vm_peak[i]->GetParError(0);
      merr_vm_peak[i]      = f_vm_peak[i]->GetParError(1);
      dtanla_vm_peak[i]    = b_vm_peak[i]/(250.-m_vm_peak[i]);
      dtanlaerr_vm_peak[i] = berr_vm_peak[i]/(250.-m_vm_peak[i]);

      //can[idx]->cd(2);
      //g_vm_peak[i]->Draw("AP");

     
      g_vp_deco[i] = getTGraphFromTH2(h_vp_deco[i],thetabins,0);
      g_vp_deco[i]->SetMarkerColor(2);
      g_vp_deco[i]->SetLineColor(2);
      g_vp_deco[i]->SetMarkerStyle(8);
      g_vp_deco[i]->SetMarkerSize(1);
      f_vp_deco[i]=new TF1(Form("f_vp_deco_%i",i),"pol1");
      f_vp_deco[i]->SetLineColor(2);
      g_vp_deco[i]->Fit(f_vp_deco[i]);
      s_vp_deco[i]<<"m = "<<fround(f_vp_deco[i]->GetParameter(1),1)
		  <<"   b = "<<fround(f_vp_deco[i]->GetParameter(0),1)<<endl;
     
      b_vp_deco[i]         = f_vp_deco[i]->GetParameter(0);
      m_vp_deco[i]         = f_vp_deco[i]->GetParameter(1);
      berr_vp_deco[i]      = f_vp_deco[i]->GetParError(0);
      merr_vp_deco[i]      = f_vp_deco[i]->GetParError(1);
      dtanla_vp_deco[i]    = b_vp_deco[i]/(250.-m_vp_deco[i]);
      dtanlaerr_vp_deco[i] = berr_vp_deco[i]/(250.-m_vp_deco[i]);

      //can[idx]->cd(3);
      //g_vp_deco[i]->Draw("AP");

      g_vm_deco[i] = getTGraphFromTH2(h_vm_deco[i],thetabins,0,1);
      g_vm_deco[i]->SetMarkerColor(2);
      g_vm_deco[i]->SetLineColor(2);
      g_vm_deco[i]->SetMarkerStyle(8);
      g_vm_deco[i]->SetMarkerSize(1);
      f_vm_deco[i]=new TF1(Form("f_vm_deco_%i",i),"pol1");
      f_vm_deco[i]->SetLineColor(2);
      g_vm_deco[i]->Fit(f_vm_deco[i]);
      s_vm_deco[i]<<"m = "<<fround(f_vm_deco[i]->GetParameter(1),1)
		  <<"   b = "<<fround(f_vm_deco[i]->GetParameter(0),1)<<endl;
     
      b_vm_deco[i]         = f_vm_deco[i]->GetParameter(0);
      m_vm_deco[i]         = f_vm_deco[i]->GetParameter(1);
      berr_vm_deco[i]      = f_vm_deco[i]->GetParError(0);
      merr_vm_deco[i]      = f_vm_deco[i]->GetParError(1);
      dtanla_vm_deco[i]    = b_vm_deco[i]/(250.-m_vm_deco[i]);
      dtanlaerr_vm_deco[i] = berr_vm_deco[i]/(250.-m_vm_deco[i]);

      //can[idx]->cd(4);
      //g_vm_deco[i]->Draw("AP");

    }

    float dtanla_peak[8];
    float dtanla_deco[8];
    float dtanlaerr_peak[8];
    float dtanlaerr_deco[8];


    for(int i = 0 ; i < 8 ; i ++){

      //cout<<" dtanla peak (V+) "<<dtanla_vp_peak[i]
      //	  <<" dtanla peak (V-) "<<dtanla_vm_peak[i]
      //	  <<" dtanla deco (V+) "<<dtanla_vp_deco[i]
      //	  <<" dtanla deco (V-) "<<dtanla_vm_deco[i]<<endl;

      dtanla_peak[i]    = 0.5*(dtanla_vp_peak[i]-dtanla_vm_peak[i]);  
      dtanla_deco[i]    = 0.5*(dtanla_vp_deco[i]-dtanla_vm_deco[i]);  
      dtanlaerr_peak[i] = (1./sqrt(2))*(dtanlaerr_vp_peak[i]+dtanlaerr_vm_peak[i]);
      dtanlaerr_deco[i] = (1./sqrt(2))*(dtanlaerr_vp_deco[i]+dtanlaerr_vm_deco[i]);

      //cout<<" dtanla peak "<<dtanla_peak[i]
      //	  <<" dtanla deco "<<dtanla_deco[i]<<endl;

    }

    cout<<"Peak tan(LA) : "<<dtanla_peak[3]<<" +/- "<<dtanlaerr_peak[3]<<endl;
    cout<<"Deco tan(LA) : "<<dtanla_deco[3]<<" +/- "<<dtanlaerr_deco[3]<<endl;
    

    float x[8]    = {-16,-10,-6,-2,2,6,10,16};
    float xerr[8] = {0,0,0,0,0,0,0,0};

    //TGraphErrors* gtot_vp_peak=new TGraphErrors(8,xpeak,b_vp_peak,xerrpeak,berr_vp_peak);
    //TGraphErrors* gtot_vm_peak=new TGraphErrors(8,xpeak,b_vm_peak,xerrpeak,berr_vm_peak);
    //TGraphErrors* gtot_vp_deco=new TGraphErrors(8,xdeco,b_vp_deco,xerrdeco,berr_vp_deco);
    //TGraphErrors* gtot_vm_deco=new TGraphErrors(8,xdeco,b_vm_deco,xerrdeco,berr_vm_deco);

    TGraphErrors* gtot_vp_peak=new TGraphErrors(8,xpeak,dtanla_vp_peak,xerrpeak,dtanlaerr_vp_peak);
    TGraphErrors* gtot_vm_peak=new TGraphErrors(8,xpeak,dtanla_vm_peak,xerrpeak,dtanlaerr_vm_peak);
    TGraphErrors* gtot_vp_deco=new TGraphErrors(8,xdeco,dtanla_vp_deco,xerrdeco,dtanlaerr_vp_deco);
    TGraphErrors* gtot_vm_deco=new TGraphErrors(8,xdeco,dtanla_vm_deco,xerrdeco,dtanlaerr_vm_deco);
 
    gtot_vp_peak->SetLineColor(2);
    gtot_vp_peak->SetMarkerColor(2);
    gtot_vm_peak->SetLineColor(2);
    gtot_vm_peak->SetMarkerColor(2);

    gtot_vp_deco->SetLineColor(4);
    gtot_vp_deco->SetMarkerColor(4);
    gtot_vm_deco->SetLineColor(4);
    gtot_vm_deco->SetMarkerColor(4);

    gtot_vp_peak->SetMarkerStyle(20);
    gtot_vp_peak->SetMarkerSize(0.5);
    gtot_vm_peak->SetMarkerStyle(20);
    gtot_vm_peak->SetMarkerSize(0.5);

    gtot_vp_deco->SetMarkerStyle(20);
    gtot_vp_deco->SetMarkerSize(0.5);
    gtot_vm_deco->SetMarkerStyle(20);
    gtot_vm_deco->SetMarkerSize(0.5);

    can[idx]->cd(1);
    TMultiGraph *gtot_vp=new TMultiGraph();
    gtot_vp->Add(gtot_vp_peak);
    gtot_vp->Add(gtot_vp_deco);
    gtot_vp->SetTitle("TOB (V+)");
    gtot_vp->Draw("AP");
    gtot_vp->GetXaxis()->SetTitle("DT Time (ns)");
    gtot_vp->GetYaxis()->SetTitle("#Deltatan(#theta_{LA})");
    gtot_vp->GetXaxis()->SetTitleSize(0.05);
    gtot_vp->GetYaxis()->SetTitleSize(0.05);
    gtot_vp->Draw("AP");

    can[idx]->cd(2);
    TMultiGraph *gtot_vm=new TMultiGraph();
    gtot_vm->Add(gtot_vm_peak);
    gtot_vm->Add(gtot_vm_deco);
    gtot_vm->SetTitle("TOB (V-)");
    gtot_vm->Draw("AP");
    gtot_vm->GetXaxis()->SetTitle("DT Time (ns)");
    gtot_vm->GetYaxis()->SetTitle("#Deltatan(#theta_{LA})");
    gtot_vm->GetXaxis()->SetTitleSize(0.05);
    gtot_vm->GetYaxis()->SetTitleSize(0.05);
    gtot_vm->Draw("AP");
     
    idx++;

    TCanvas *ctemp = new TCanvas("ctemp","",800,600);
    ctemp->cd();
    TGraphErrors *gpeak = averageTGraph(gtot_vp_peak,gtot_vm_peak);
    TGraphErrors *gdeco = averageTGraph(gtot_vp_deco,gtot_vm_deco);
    
    gpeak->SetLineColor(6);
    gpeak->SetMarkerColor(6);
    gdeco->SetLineColor(8);
    gdeco->SetMarkerColor(8);
    gpeak->SetMarkerStyle(20);
    gpeak->SetMarkerSize(0.5);
    gdeco->SetMarkerStyle(20);
    gdeco->SetMarkerSize(0.5);
   
    gpeak->SetTitle("TOB");
    gpeak->GetXaxis()->SetTitle("DT Time (ns)");
    gpeak->GetYaxis()->SetTitle("#Deltatan(#theta_{LA})");
    gpeak->Draw("AP");
    gdeco->Draw("sameP");
    gpeak->GetYaxis()->SetRangeUser(-0.11,-0.05);
    
 }

  
  for(int ican=0;ican<ncan;ican++){
    can[ican]->Modified();
    can[ican]->Update();
    if(printgif) can[ican]->Print(Form("%s/plots/%s_TOB.gif",dir,cantitles[ican]));
  
  }
}


TGraphErrors* averageTGraph(TGraphErrors* g1, TGraphErrors *g2){

  Double_t* x1  = g1->GetX();
  Double_t* y1  = g1->GetY();
  Double_t* ex1 = g1->GetEX();
  Double_t* ey1 = g1->GetEY();
  Int_t n1      = g1->GetN();
  
  Double_t* x2  = g2->GetX();
  Double_t* y2  = g2->GetY();
  Double_t* ex2 = g2->GetEX();
  Double_t* ey2 = g2->GetEY();
  Int_t n2      = g2->GetN();

  assert(n1 == n2);
  
  Double_t* x  = new Double_t[n1];
  Double_t* y  = new Double_t[n1];
  Double_t* ex = new Double_t[n1];
  Double_t* ey = new Double_t[n1];

  for(int i = 0 ; i < n1 ; i++){

    ex[i] = sqrt( 1. / ( 1. / pow(ex1[i],2) + 1. / pow(ex2[i],2) ) );
    x[i]  = (x1[i] / pow(ex1[i],2) + x2[i] / pow(ex2[i],2) ) * pow(ex[i],2);

    ey[i] = sqrt( 1. / ( 1. / pow(ey1[i],2) + 1. / pow(ey2[i],2) ) );
    y[i]  = (y1[i] / pow(ey1[i],2) + y2[i] / pow(ey2[i],2) ) * pow(ey[i],2) - 0.09095;

    ex[i] = 0.;
  }

  TGraphErrors *g = new TGraphErrors(n1,x,y,ex,ey);
  return g;

}















//   fpeak.cd();
//   TrackerOfflineValidation->cd();
//   MyStrip->cd();
//   TH1F *dupeak =            cloneHist(h_Xprime_TOBBarrel_3);
//   TH1F *dwpeak =            cloneHist(h_resXprimeOverTheta_TOBBarrel_3);
//   TH1F *dttimepeak =        (TH1F*)h_dttime_TOBBarrel_3->Clone();
//   TH1F *dttimeerrpeak =     (TH1F*)h_dttimeerr_TOBBarrel_3->Clone();
//   TH1F *ndtpeak =           (TH1F*)h_ndt_TOBBarrel_3->Clone();
//   TH1F *chargepeak =        (TH1F*)h_charge_TOBBarrel_3->Clone(); 
//   TH1F *nstripspeak =       (TH1F*)h_nstrips_TOBBarrel_3->Clone(); 
//   TH1F *dwpeak_EC4thick =   (TH1F*)h_resXprimeOverThetaThick_TECEndcap_4->Clone(); 
//   TH1F *dwpeak_EC5thick =   (TH1F*)h_resXprimeOverThetaThick_TECEndcap_5->Clone(); 
//   TH1F *dwpeak_EC4thin =    (TH1F*)h_resXprimeOverThetaThin_TECEndcap_4->Clone(); 
//   TH1F *dwpeak_EC5thin =    (TH1F*)h_resXprimeOverThetaThin_TECEndcap_5->Clone(); 

//   fdeco.cd();
//   TrackerOfflineValidation->cd();
//   MyStrip->cd();
//   TH1F *dudeco =            cloneHist(h_Xprime_TOBBarrel_3);
//   TH1F *dwdeco =            cloneHist(h_resXprimeOverTheta_TOBBarrel_3);
//   TH1F *dttimedeco =        (TH1F*)h_dttime_TOBBarrel_3->Clone();
//   TH1F *dttimeerrdeco =     (TH1F*)h_dttimeerr_TOBBarrel_3->Clone();
//   TH1F *ndtdeco =           (TH1F*)h_ndt_TOBBarrel_3->Clone();
//   TH1F *chargedeco =        (TH1F*)h_charge_TOBBarrel_3->Clone(); 
//   TH1F *nstripsdeco =       (TH1F*)h_nstrips_TOBBarrel_3->Clone(); 
//   TH1F *dwdeco_EC4thick =   (TH1F*)h_resXprimeOverThetaThick_TECEndcap_4->Clone(); 
//   TH1F *dwdeco_EC5thick =   (TH1F*)h_resXprimeOverThetaThick_TECEndcap_5->Clone(); 
//   TH1F *dwdeco_EC4thin =    (TH1F*)h_resXprimeOverThetaThin_TECEndcap_4->Clone(); 
//   TH1F *dwdeco_EC5thin =    (TH1F*)h_resXprimeOverThetaThin_TECEndcap_5->Clone(); 
