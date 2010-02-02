{

bool printeps=false;
bool printgif=false;

gROOT->LoadMacro("scripts/setStats.C");
gROOT->LoadMacro("scripts/plotHistsTH2.C");
 
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetPalette(1);

TFile fdeco("rootfiles/deco_TH2.root");
TrackerOfflineValidation->cd();
MyStrip->cd();
TFile fpeak("rootfiles/peak_TH2.root");
TrackerOfflineValidation->cd();
MyStrip->cd();

int dx=900;
int dy=600;
TCanvas *c1=new TCanvas("c1","",dx,dy);
c1->Divide(3,2);
TCanvas *c2=new TCanvas("c2","",dx,dy);
c2->Divide(3,2);
TCanvas *c3=new TCanvas("c3","",dx,dy);
c3->Divide(3,2);
TCanvas *c4=new TCanvas("c4","",dx,dy);
c4->Divide(3,2);
TCanvas *c5=new TCanvas("c5","",800,400);
c5->Divide(2,1);

const int nplots=6;

TH2* hwvstheta_TOB1_peak[nplots];
TH2* hwvstheta_TOB1_deco[nplots];
TH2* hwvstheta_TOB2_peak[nplots];
TH2* hwvstheta_TOB2_deco[nplots];

char* hb1peak[]={"TOB- L1 PEAK",
		 "TOB- L2 PEAK",
		 "TOB- L3 PEAK",
		 "TOB- L4 PEAK",
		 "TOB- L5 PEAK",
		 "TOB- L6 PEAK"};

char* hb2peak[]={"TOB+ L1 PEAK",
		 "TOB+ L2 PEAK",
		 "TOB+ L3 PEAK",
		 "TOB+ L4 PEAK",
		 "TOB+ L5 PEAK",
		 "TOB+ L6 PEAK"};

char* hb1deco[]={"TOB- L1 DECO",
		 "TOB- L2 DECO",
		 "TOB- L3 DECO",
		 "TOB- L4 DECO",
		 "TOB- L5 DECO",
		 "TOB- L6 DECO"};

char* hb2deco[]={"TOB+ L1 DECO",
		 "TOB+ L2 DECO",
		 "TOB+ L3 DECO",
		 "TOB+ L4 DECO",
		 "TOB+ L5 DECO",
		 "TOB+ L6 DECO"};

fpeak.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TH2F* hpeak=(TH2F*)h_resXprimeVsTheta_TOBBarrel_3->Clone();
TOBBarrel_4->cd();
TOBHalfBarrel_1->cd();

hwvstheta_TOB1_peak[0]=(TH2F*)h_resXprimeVsTheta_TOBLayer_0->Clone();
hwvstheta_TOB1_peak[1]=(TH2F*)h_resXprimeVsTheta_TOBLayer_1->Clone();
hwvstheta_TOB1_peak[2]=(TH2F*)h_resXprimeVsTheta_TOBLayer_2->Clone();
hwvstheta_TOB1_peak[3]=(TH2F*)h_resXprimeVsTheta_TOBLayer_3->Clone();
hwvstheta_TOB1_peak[4]=(TH2F*)h_resXprimeVsTheta_TOBLayer_4->Clone();
hwvstheta_TOB1_peak[5]=(TH2F*)h_resXprimeVsTheta_TOBLayer_5->Clone();

fpeak.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TOBBarrel_4->cd();
TOBHalfBarrel_2->cd();

hwvstheta_TOB2_peak[0]=(TH2F*)h_resXprimeVsTheta_TOBLayer_0->Clone();
hwvstheta_TOB2_peak[1]=(TH2F*)h_resXprimeVsTheta_TOBLayer_1->Clone();
hwvstheta_TOB2_peak[2]=(TH2F*)h_resXprimeVsTheta_TOBLayer_2->Clone();
hwvstheta_TOB2_peak[3]=(TH2F*)h_resXprimeVsTheta_TOBLayer_3->Clone();
hwvstheta_TOB2_peak[4]=(TH2F*)h_resXprimeVsTheta_TOBLayer_4->Clone();
hwvstheta_TOB2_peak[5]=(TH2F*)h_resXprimeVsTheta_TOBLayer_5->Clone();

fdeco.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TH2F* hdeco=(TH2F*)h_resXprimeVsTheta_TOBBarrel_3->Clone();
TOBBarrel_4->cd();
TOBHalfBarrel_1->cd();

hwvstheta_TOB1_deco[0]=(TH2F*)h_resXprimeVsTheta_TOBLayer_0->Clone();
hwvstheta_TOB1_deco[1]=(TH2F*)h_resXprimeVsTheta_TOBLayer_1->Clone();
hwvstheta_TOB1_deco[2]=(TH2F*)h_resXprimeVsTheta_TOBLayer_2->Clone();
hwvstheta_TOB1_deco[3]=(TH2F*)h_resXprimeVsTheta_TOBLayer_3->Clone();
hwvstheta_TOB1_deco[4]=(TH2F*)h_resXprimeVsTheta_TOBLayer_4->Clone();
hwvstheta_TOB1_deco[5]=(TH2F*)h_resXprimeVsTheta_TOBLayer_5->Clone();

fdeco.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TOBBarrel_4->cd();
TOBHalfBarrel_2->cd();

hwvstheta_TOB2_deco[0]=(TH2F*)h_resXprimeVsTheta_TOBLayer_0->Clone();
hwvstheta_TOB2_deco[1]=(TH2F*)h_resXprimeVsTheta_TOBLayer_1->Clone();
hwvstheta_TOB2_deco[2]=(TH2F*)h_resXprimeVsTheta_TOBLayer_2->Clone();
hwvstheta_TOB2_deco[3]=(TH2F*)h_resXprimeVsTheta_TOBLayer_3->Clone();
hwvstheta_TOB2_deco[4]=(TH2F*)h_resXprimeVsTheta_TOBLayer_4->Clone();
hwvstheta_TOB2_deco[5]=(TH2F*)h_resXprimeVsTheta_TOBLayer_5->Clone();


for(int i=0;i<6;i++){

  plotHistsTH2(hwvstheta_TOB1_peak[i],hb1peak[i],"tan(#theta_{trk})-tan(#theta_{L})","#Delta u [#mum]",-5,5,-1000,1000,0,c1,i+1);

  plotHistsTH2(hwvstheta_TOB1_deco[i],hb1deco[i],"tan(#theta_{trk})-tan(#theta_{L})","#Delta u [#mum]",-5,5,-1000,1000,0,c2,i+1);

  plotHistsTH2(hwvstheta_TOB2_peak[i],hb2peak[i],"tan(#theta_{trk})-tan(#theta_{L})","#Delta u [#mum]",-5,5,-1000,1000,0,c3,i+1);

  plotHistsTH2(hwvstheta_TOB2_deco[i],hb2deco[i],"tan(#theta_{trk})-tan(#theta_{L})","#Delta u [#mum]",-5,5,-1000,1000,0,c4,i+1);

}

plotHistsTH2(hpeak,"TOB PEAK","tan(#theta_{trk})-tan(#theta_{L})","#Delta u [#mum]",-2,2,-200,200,0,c5,1);
plotHistsTH2(hdeco,"TOB DECO","tan(#theta_{trk})-tan(#theta_{L})","#Delta u [#mum]",-2,2,-200,200,0,c5,2);



if(printeps){
   c1->Print("plots/dutantheta_HB1_peak.eps");
   c2->Print("plots/dutantheta_HB1_deco.eps");
   c3->Print("plots/dutantheta_HB2_peak.eps");
   c4->Print("plots/dutantheta_HB2_deco.eps");
}
if(printgif){
   c1->Print("plots/dutantheta_HB1_peak.gif");
   c2->Print("plots/dutantheta_HB1_deco.gif");
   c3->Print("plots/dutantheta_HB2_peak.gif");
   c4->Print("plots/dutantheta_HB2_deco.gif");
}

}














//hpeak[i] = static_cast<TH1D *>((fpeak)->Get(Form("%s",variables[i])));
//hdeco[i] = static_cast<TH1D *>((fdeco)->Get(Form("%s",variables[i])));

// const string variables[nplots]={"h_Xprime_TIBBarrel_0",
// 				"h_Xprime_TECEndcap_4",
// 				"h_Xprime_TIDEndcap_1",
// 				"h_Xprime_TOBBarrel_3",
// 				"h_Xprime_TECEndcap_5",
// 				"h_Xprime_TIDEndcap_2"
// }


/*
TH1* cloneHist(TH1* hist){

  stringstream name;
  name<<hist->GetName()<<"_clone";

  TH1F* h=new TH1F(name.str().c_str(),hist->GetTitle(),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  for(int ibin=1;ibin<=hist->GetNbinsX();ibin++){
    h->SetBinContent(ibin,hist->GetBinContent(ibin));
  }
  return h;
  
}

void setStats(TH1* s,TH1* r, double startingY, double startingX = .1,bool fit){
  if (startingY<0){
    s->SetStats(0);
    r->SetStats(0);
  } else {
    //gStyle->SetOptStat(1001);
    
    if (fit){
      s->Fit("gaus");
      TF1* f1 = (TF1*) s->GetListOfFunctions()->FindObject("gaus");
      f1->SetLineColor(2);
      f1->SetLineWidth(1);
    }
    s->Draw();
    gPad->Update(); 
    TPaveStats* st1 = (TPaveStats*) s->GetListOfFunctions()->FindObject("stats");
    if (fit) {st1->SetOptFit(0010);    st1->SetOptStat(1001);}
    st1->SetX1NDC(startingX);
    st1->SetX2NDC(startingX+0.30);
    st1->SetY1NDC(startingY+0.20);
    st1->SetY2NDC(startingY+0.35);
    st1->SetTextColor(2);
    if (fit) {
      r->Fit("gaus");
      TF1* f2 = (TF1*) r->GetListOfFunctions()->FindObject("gaus");
      f2->SetLineColor(4);
      f2->SetLineWidth(1);    
    }
    r->Draw("sames");
    gPad->Update(); 
    TPaveStats* st2 = (TPaveStats*) r->GetListOfFunctions()->FindObject("stats");
    if (fit) {st2->SetOptFit(0010);    st2->SetOptStat(1001);}
    st2->SetX1NDC(startingX);
    st2->SetX2NDC(startingX+0.30);
    st2->SetY1NDC(startingY);
    st2->SetY2NDC(startingY+0.15);
    st2->SetTextColor(4);
  }
}

void plotHists(TH1* h1, TH1* h2, char* title, char* xtitle,float xmin, float xmax,int rebin){
  if(rebin>1){
    h1->Rebin(rebin);
    h2->Rebin(rebin);
  }
  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(xtitle);
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h1->Draw();
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->Draw("sames");
  float max=h1->GetMaximum();
  if(h2->GetMaximum()>max)max=h2->GetMaximum();
  h1->SetMaximum(1.05*max);
}

void run(){
*/
