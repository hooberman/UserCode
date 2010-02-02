{

//Config params-----------------------------------------------
bool printeps=false;
bool printgif=true;
bool fit=false;
//string dir        = "peaksqlite";  //peak sqlite file
//string dir        = "nosqlite";  //no sqlite file
//string dir        = "decosqlite";  //deco sqlite file
string dir="crabjobs/trial1";

//open root files---------------------------------------------
string decofile    = dir+"/root/deco.root";
string peakfile    = dir+"/root/peak.root";
string decofileTH2 = dir+"/root/deco_TH2.root";
string peakfileTH2 = dir+"/root/peak_TH2.root";

TFile fdeco(decofile.c_str());
TFile fpeak(peakfile.c_str());
TFile fdecoTH2(decofileTH2.c_str());
TFile fpeakTH2(peakfileTH2.c_str());

//load macros and set style-----------------------------------
gROOT->LoadMacro("scripts/cloneHist.C");
gROOT->LoadMacro("scripts/setStats.C");
gROOT->LoadMacro("scripts/plotHists.C"); 
gROOT->LoadMacro("scripts/plotHistsTH2.C");
gROOT->LoadMacro("scripts/suppressHist.C");
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat("mr");
gStyle->SetOptFit(0);

//make canvases/legend----------------------------------------
int dx=1200;
int dy=600;
TCanvas *c1=new TCanvas("c1","",dx,dy);
c1->Divide(2,1);
TCanvas *c2=new TCanvas("c2","",dx,dy);
c2->Divide(2,1);
TCanvas *c3=new TCanvas("c3","",dx,dy);
c3->Divide(2,1);
TCanvas *c4=new TCanvas("c4","",dx,2*dy);
c4->Divide(2,2);
TCanvas *c5=new TCanvas("c5","",1200,800);
c5->Divide(4,2);
//TCanvas *c6=new TCanvas("c6","",1200,800);
//c6->Divide(4,2);
TCanvas *c7=new TCanvas("c7","",800,600);
// TCanvas *c8=new TCanvas("c8","",dx,dy);
// c8->Divide(2,1);
// TCanvas *c9=new TCanvas("c9","",dx,dy);
// c9->Divide(2,1);
TCanvas *c10=new TCanvas("c10","",dx,dy);
c10->Divide(2,1);
TCanvas *c11=new TCanvas("c11","",dx,dy);
c11->Divide(2,1);
TCanvas *c12=new TCanvas("c12","",1200,800);
c12->Divide(4,2);
TCanvas *c13=new TCanvas("c13","",800,600);

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

//get 1D histos-----------------------------------------------
fpeak.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TH1F *dupeak=        cloneHist(h_Xprime_TECBarrel_3);
TH1F *dwpeak=        cloneHist(h_resXprimeOverTheta_TECBarrel_3);
TH1F *dttimepeak=    (TH1F*)h_dttime_TECBarrel_3->Clone();
TH1F *dttimeerrpeak= (TH1F*)h_dttimeerr_TECBarrel_3->Clone();
TH1F *ndtpeak=       (TH1F*)h_ndt_TECBarrel_3->Clone();
TH1F *chargepeak=    (TH1F*)h_charge_TECBarrel_3->Clone(); 
TH1F *nstripspeak=   (TH1F*)h_nstrips_TECBarrel_3->Clone(); 

fdeco.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TH1F *dudeco=        cloneHist(h_Xprime_TECBarrel_3);
TH1F *dwdeco=        cloneHist(h_resXprimeOverTheta_TECBarrel_3);
TH1F *dttimedeco=    (TH1F*)h_dttime_TECBarrel_3->Clone();
TH1F *dttimeerrdeco= (TH1F*)h_dttimeerr_TECBarrel_3->Clone();
TH1F *ndtdeco=       (TH1F*)h_ndt_TECBarrel_3->Clone();
TH1F *chargedeco=    (TH1F*)h_charge_TECBarrel_3->Clone(); 
TH1F *nstripsdeco=   (TH1F*)h_nstrips_TECBarrel_3->Clone(); 

//draw 1D histo-----------------------------------------------
c1->cd(1);
setStats(dupeak,dudeco,0.5,0.65,0.15,false);  
plotHists(dupeak,dudeco,"TEC","#Delta u [#mum]",-500,500,1,3);
leg1->Draw();
c1->cd(2);
setStats(dwpeak,dwdeco,0.5,0.65,0.15,false);  
plotHists(dwpeak,dwdeco,"TEC","#Delta w [#mum]",-1000,1000,1,3);
leg1->Draw();

c4->cd(1);
setStats(dttimepeak,dttimedeco,0.5,0.65,0.15,false);  
plotHists(dttimepeak,dttimedeco,"TEC","DT time (ns)",-50,50,1,0);
leg1->Draw();
c4->cd(2);
setStats(dttimeerrpeak,dttimeerrdeco,0.5,0.65,0.15,false);  
plotHists(dttimeerrpeak,dttimeerrdeco,"TEC","#delta DT time (ns)",0,100,1,0);
TLegend *leg2=new TLegend(0.4,0.55,0.6,0.85);
leg2->AddEntry(h1,"PEAK");
leg2->AddEntry(h2,"DECO");
leg2->SetBorderSize(1);
leg2->SetFillColor(0);
leg2->Draw();
c4->cd(3);
setStats(ndtpeak,ndtdeco,0.5,0.65,0.15,false);  
plotHists(ndtpeak,ndtdeco,"TEC","n DT Hits",0,50,1,0);
leg1->Draw();

c10->cd(1);
setStats(chargepeak,chargedeco,0.5,0.65,0.15,false);  
plotHists(chargepeak,chargedeco,"TEC","Charge",0,1000,1,0);
c10->cd(2);
setStats(nstripspeak,nstripsdeco,0.5,0.65,0.15,false);  
plotHists(nstripspeak,nstripsdeco,"TEC","NStrips",0,10,1,0);

//get 2D histos-----------------------------------------------
//gStyle->SetOptStat(0);

fpeakTH2.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TH2F* duthetapeak =  (TH2F*)h_resXprimeVsTheta_TECBarrel_3->Clone();
//TH2F* dutrkpeak =  (TH2F*)h_duvstrkangle_TECBarrel_3->Clone();
//TH2F* dulorpeak =  (TH2F*)h_duvslorangle_TECBarrel_3->Clone();
TH2F* duvsdtpeak =   (TH2F*)h_duvsdttime_TECBarrel_3->Clone();
TH2F* dwvsdtpeak =   (TH2F*)h_dwvsdttime_TECBarrel_3->Clone();
TH2F* chvsdtpeak = (TH2F*)h_chargevsdttime_TECBarrel_3->Clone();

fdecoTH2.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TH2F* duthetadeco =  (TH2F*)h_resXprimeVsTheta_TECBarrel_3->Clone();
//TH2F* dutrkdeco =  (TH2F*)h_duvstrkangle_TECBarrel_3->Clone();
//TH2F* dulordeco =  (TH2F*)h_duvslorangle_TECBarrel_3->Clone();
TH2F* duvsdtdeco =   (TH2F*)h_duvsdttime_TECBarrel_3->Clone();
TH2F* dwvsdtdeco =   (TH2F*)h_dwvsdttime_TECBarrel_3->Clone();
TH2F* chvsdtdeco = (TH2F*)h_chargevsdttime_TECBarrel_3->Clone();

//draw 2D histo-----------------------------------------------
plotHistsTH2(duthetapeak,
	     "TEC PEAK","tan(#theta_{trk})-tan(#theta_{L})",
	     "#Delta u [#mum]",-2,2,-200,200,0,c2,1);
plotHistsTH2(duthetadeco,
	     "TEC DECO","tan(#theta_{trk})-tan(#theta_{L})",
	     "#Delta u [#mum]",-2,2,-200,200,0,c2,2);
plotHistsTH2(dwvsdtpeak,
	     "TEC PEAK","DT Time [ns]",
	     "#Delta w [#mum]",-25,25,-500,500,0,c3,1);
plotHistsTH2(dwvsdtdeco,
	     "TEC DECO","DT Time [ns]",
	     "#Delta w [#mum]",-25,25,-500,500,0,c3,2);
plotHistsTH2(chvsdtpeak,
	     "TEC PEAK","DT Time [ns]",
	     "Charge",-25,25,0,500,0,c11,1);
plotHistsTH2(chvsdtdeco,
	     "TEC DECO","DT Time [ns]",
	     "Charge",-25,25,0,500,0,c11,2);
// plotHistsTH2(dutrkpeak,
// 	     "TEC PEAK","tan(#theta_{trk})",
// 	     "#Delta u [#mum]",-2,2,-200,200,0,c8,1);
// plotHistsTH2(dutrkdeco,
// 	     "TEC DECO","tan(#theta_{trk})",
// 	     "#Delta u [#mum]",-2,2,-200,200,0,c8,2);
// plotHistsTH2(dulorpeak,
// 	     "TEC PEAK","tan(#theta_{LA})",
// 	     "#Delta u [#mum]",-0.2,0.2,-200,200,0,c9,1);
// plotHistsTH2(dulordeco,
// 	     "TEC DECO","tan(#theta_{LA})",
// 	     "#Delta u [#mum]",-0.2,0.2,-200,200,0,c9,2);


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



string projtitle[nslices]={"TEC (-20 ns < dt < -12 ns)",
			   "TEC (-12 ns < dt < -8 ns)",
			   "TEC (-8 ns < dt < -4 ns)",
			   "TEC (-4 ns < dt < 0 ns)",
			   "TEC (0 ns < dt < 4 ns)",
			   "TEC (4 ns < dt < 8 ns)",
			   "TEC (8 ns < dt < 12 ns)",
			   "TEC (12 ns < dt < 20 ns)"};
 

for(int islice=0;islice<nslices;islice++){

  //delta w vs. dt time---------------------------------------------------------
  dwvsdtpeakslice[islice]=suppressHist((TH2F*)dwvsdtpeak->Clone(),xmin[islice],xmax[islice]);
  dwvsdtdecoslice[islice]=suppressHist((TH2F*)dwvsdtdeco->Clone(),xmin[islice],xmax[islice]);
   
  dwpeakproj[islice]=(TH1F*)dwvsdtpeakslice[islice]->ProjectionY()->Clone();
  dwdecoproj[islice]=(TH1F*)dwvsdtdecoslice[islice]->ProjectionY()->Clone();
  dtpeakproj[islice]=(TH1F*)dwvsdtpeakslice[islice]->ProjectionX()->Clone();
  dtdecoproj[islice]=(TH1F*)dwvsdtdecoslice[islice]->ProjectionX()->Clone();

  c5->cd(islice+1);
  
  if(fit){
    //triple gaussian
    TF1 *ftgp=new TF1("ftgp",
		      "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))");
    TF1 *ftgd=new TF1("ftgd",
		      "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))");
    
    ftgp->SetParameters(0.1,0.,100,0.1,200,0.1,300);
    ftgd->SetParameters(0.1,0.,100,0.1,200,0.1,300);
    
    ftgp->SetLineWidth(1);
    ftgd->SetLineWidth(1);
    ftgp->SetLineColor(2);
    ftgd->SetLineColor(4);
    
    dwpeakproj[islice]->Sumw2();
    dwdecoproj[islice]->Sumw2();
    
    dwpeakproj[islice]->Scale(1./dwpeakproj[islice]->GetEntries());
    dwdecoproj[islice]->Scale(1./dwdecoproj[islice]->GetEntries());
    dwpeakproj[islice]->Fit(ftgp);
    dwdecoproj[islice]->Fit(ftgd);
    
  }

  setStats(dwpeakproj[islice],dwdecoproj[islice],0.5,0.65,0.15,false);  
  plotHists(dwpeakproj[islice],dwdecoproj[islice],
	    projtitle[islice].c_str(),"#Delta w [#mum]",-1000,1000,2,0);

  if(fit){
    ypeak[islice]   = ftgp->GetParameter(1);
    ydeco[islice]   = ftgd->GetParameter(1);
    ypeakerr[islice]= ftgp->GetParError(1);
    ydecoerr[islice]= ftgd->GetParError(1);
  }else{
    ypeak[islice]   = dwpeakproj[islice]->GetMean(1);
    ydeco[islice]   = dwdecoproj[islice]->GetMean(1);
    ypeakerr[islice]= dwpeakproj[islice]->GetEntries() > 0 ? 
      dwpeakproj[islice]->GetRMS(1)/sqrt(dwpeakproj[islice]->GetEntries()) : 0.;
    ydecoerr[islice]= dwdecoproj[islice]->GetEntries() > 0 ?
      dwdecoproj[islice]->GetRMS(1)/sqrt(dwdecoproj[islice]->GetEntries()) : 0.;
  }
  
  //c6->cd(islice+1);
  //setStats(dtpeakproj[islice],dtdecoproj[islice],0.5,0.65,0.15,false);  
  //plotHists(dtpeakproj[islice],dtdecoproj[islice],"TEC","DT Time [ns]",-20,20,1,0);

  xpeak[islice]    = dtpeakproj[islice]->GetMean(1);
  xdeco[islice]    = dtdecoproj[islice]->GetMean(1);
  xpeakerr[islice] = dtpeakproj[islice]->GetEntries() > 0 ?
    dtpeakproj[islice]->GetRMS(1)/sqrt(dtpeakproj[islice]->GetEntries()) : 0.;
  xdecoerr[islice] = dtpeakproj[islice]->GetEntries() > 0 ?
    dtpeakproj[islice]->GetRMS(1)/sqrt(dtpeakproj[islice]->GetEntries()) : 0.;

  //charge w vs. dt time---------------------------------------------------------
  chvsdtpeakslice[islice]=suppressHist((TH2F*)chvsdtpeak->Clone(),xmin[islice],xmax[islice]);
  chvsdtdecoslice[islice]=suppressHist((TH2F*)chvsdtdeco->Clone(),xmin[islice],xmax[islice]);
   
  chpeakproj[islice]=(TH1F*)chvsdtpeakslice[islice]->ProjectionY()->Clone();
  chdecoproj[islice]=(TH1F*)chvsdtdecoslice[islice]->ProjectionY()->Clone();
  dtpeakproj2[islice]=(TH1F*)chvsdtpeakslice[islice]->ProjectionX()->Clone();
  dtdecoproj2[islice]=(TH1F*)chvsdtdecoslice[islice]->ProjectionX()->Clone();

  c12->cd(islice+1);

  if(fit){
    //landau
    TF1 *flanp=new TF1("flanp","landau");
    TF1 *fland=new TF1("fland","landau");
    
    flanp->SetLineWidth(1);
    fland->SetLineWidth(1);
    flanp->SetLineColor(2);
    fland->SetLineColor(4);
    
    chpeakproj[islice]->Sumw2();
    chdecoproj[islice]->Sumw2();
    
    chpeakproj[islice]->Scale(1./chpeakproj[islice]->GetEntries());
    chdecoproj[islice]->Scale(1./chdecoproj[islice]->GetEntries());
    chpeakproj[islice]->Fit(flanp);
    chdecoproj[islice]->Fit(fland);
  }    

  setStats(chpeakproj[islice],chdecoproj[islice],0.5,0.65,0.15,false);  
  plotHists(chpeakproj[islice],chdecoproj[islice],
	    projtitle[islice].c_str(),"Charge",0,500,1,0);
  
  if(fit){
    ychpeak[islice]   = flanp->GetParameter(1);
    ychdeco[islice]   = fland->GetParameter(1);
    ychpeakerr[islice]= flanp->GetParError(1);
    ychdecoerr[islice]= fland->GetParError(1);
  }else{
    ychpeak[islice]   = chpeakproj[islice]->GetMean(1);
    ychdeco[islice]   = chdecoproj[islice]->GetMean(1);
    ychpeakerr[islice]= chpeakproj[islice]->GetEntries() > 0 ? 
      chpeakproj[islice]->GetRMS(1)/sqrt(chpeakproj[islice]->GetEntries()) : 0.;
    ychdecoerr[islice]= chdecoproj[islice]->GetEntries() > 0 ?
      chdecoproj[islice]->GetRMS(1)/sqrt(chdecoproj[islice]->GetEntries()) : 0.;
  }
  
  //c6->cd(islice+1);
  //setStats(dtpeakproj[islice],dtdecoproj[islice],0.5,0.65,0.15,false);  
  //plotHists(dtpeakproj[islice],dtdecoproj[islice],"TEC","DT Time [ns]",-20,20,1,0);

  xchpeak[islice]    = dtpeakproj2[islice]->GetMean(1);
  xchdeco[islice]    = dtdecoproj2[islice]->GetMean(1);
  xchpeakerr[islice] = dtpeakproj2[islice]->GetEntries() > 0 ?
    dtpeakproj2[islice]->GetRMS(1)/sqrt(dtpeakproj2[islice]->GetEntries()) : 0.;
  xchdecoerr[islice] = dtpeakproj2[islice]->GetEntries() > 0 ?
    dtpeakproj2[islice]->GetRMS(1)/sqrt(dtpeakproj2[islice]->GetEntries()) : 0.;

}

//tgraph dw vs dt time-------------------------------
TGraphErrors *gpeak=new TGraphErrors(nslices,xpeak,ypeak,xpeakerr,ypeakerr);
gpeak->SetMarkerStyle(8);
gpeak->SetMarkerSize(1);
gpeak->SetMarkerColor(2);
gpeak->SetLineColor(2);
TF1 *fp=new TF1("fp","pol1");
fp->SetLineColor(2);
gpeak->Fit(fp);
stringstream speak;
speak<<"Slope = "<<fround(fp->GetParameter(1),2)<<" #mum / ns"<<endl;

TGraphErrors *gdeco=new TGraphErrors(nslices,xdeco,ydeco,xdecoerr,ydecoerr);
gdeco->SetMarkerStyle(8);
gdeco->SetMarkerSize(1);
gdeco->SetMarkerColor(4);
gdeco->SetLineColor(4);
TF1 *fd=new TF1("fd","pol1");
fd->SetLineColor(4);
gdeco->Fit(fd);
stringstream sdeco;
sdeco<<"Slope = "<<fround(fd->GetParameter(1),2)<<" #mum / ns"<<endl;

c7->cd();
TMultiGraph *mg=new TMultiGraph();
mg->Add(gpeak);
mg->Add(gdeco);
mg->SetTitle("TEC");
mg->Draw("AP");
mg->GetXaxis()->SetTitle("<DT Time> [ns]");
mg->GetYaxis()->SetTitle("<#Delta_{W}> [#mum]");
mg->SetTitle("TEC");

TLatex *t=new TLatex();
t->SetNDC();
t->SetTextColor(2);
t->DrawLatex(0.15,0.75,speak.str().c_str());
t->SetTextColor(4);
t->DrawLatex(0.15,0.65,sdeco.str().c_str());


//tgraph charge vs dt time
TGraphErrors *gchpeak=new TGraphErrors(nslices,xchpeak,ychpeak,xchpeakerr,ychpeakerr);
gchpeak->SetMarkerStyle(8);
gchpeak->SetMarkerSize(1);
gchpeak->SetMarkerColor(2);
gchpeak->SetLineColor(2);
//TF1 *fchp=new TF1("fchp","pol1");
//fchp->SetLineColor(2);
//gchpeak->Fit(fchp);
//stringstream schpeak;
//schpeak<<"Slope = "<<fround(fchp->GetParameter(1),2)<<" #mum / ns"<<endl;

TGraphErrors *gchdeco=new TGraphErrors(nslices,xchdeco,ychdeco,xchdecoerr,ychdecoerr);
gchdeco->SetMarkerStyle(8);
gchdeco->SetMarkerSize(1);
gchdeco->SetMarkerColor(4);
gchdeco->SetLineColor(4);
//TF1 *fchd=new TF1("fchd","pol1");
//fchd->SetLineColor(4);
//gchdeco->Fit(fchd);
//stringstream schdeco;
//schdeco<<"Slope = "<<fround(fchd->GetParameter(1),2)<<" #mum / ns"<<endl;

c13->cd();
TMultiGraph *mgch=new TMultiGraph();
mgch->Add(gchpeak);
mgch->Add(gchdeco);
mgch->SetTitle("TEC");
mgch->Draw("AP");
mgch->GetXaxis()->SetTitle("<DT Time> [ns]");
mgch->GetYaxis()->SetTitle("<Charge>");
mgch->SetTitle("TEC");

// TLatex *t=new TLatex();
// t->SetNDC();
// t->SetTextColor(2);
// t->DrawLatex(0.15,0.75,speak.str().c_str());
// t->SetTextColor(4);
// t->DrawLatex(0.15,0.65,sdeco.str().c_str());


c1->Modified();
c1->Update();
c2->Modified();
c2->Update();
c3->Modified();
c3->Update();
c4->Modified();
c4->Update();
c5->Modified();
c5->Update();
c7->Modified();
c7->Update();
c10->Modified();
c10->Update();
c11->Modified();
c11->Update();
c12->Modified();
c12->Update();
c13->Modified();
c13->Update();

if(printgif){
  c1->Print (Form("%s/plots/deltau_deltaw_TEC.gif",dir));
  c2->Print (Form("%s/plots/duvstantheta_TEC.gif",dir));
  c3->Print (Form("%s/plots/duvsdt_TEC.gif",dir));
  c4->Print (Form("%s/plots/dttime_TEC.gif",dir));
  c5->Print (Form("%s/plots/dwslice_TEC.gif",dir));
  c7->Print (Form("%s/plots/dwvsdt_TEC.gif",dir));
  c10->Print(Form("%s/plots/charge_nstrips_TEC.gif",dir));
  c11->Print(Form("%s/plots/chargevsdt_TEC.gif",dir));
  c12->Print(Form("%s/plots/chargeslice_TEC.gif",dir));
  c13->Print(Form("%s/plots/chargevsdtslice_TEC.gif",dir));

  //c6->Print(Form("localjobs/%s/plots/dtslice_TEC.gif",dir));
}


}
