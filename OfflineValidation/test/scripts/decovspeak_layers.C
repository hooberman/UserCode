{
string dir="Mon_Nov_23_17:27:01_CST_2009";
string decofile="localjobs/"+dir+"/root/deco.root";
string peakfile="localjobs/"+dir+"/root/peak.root";

bool printeps=false;
bool printgif=false;

gROOT->LoadMacro("scripts/cloneHist.C");
gROOT->LoadMacro("scripts/setStats.C");
gROOT->LoadMacro("scripts/plotHists.C");
 
gROOT->SetStyle("Plain");
gStyle->SetOptStat("mr");
gStyle->SetOptFit(0);


TFile fdeco(decofile.c_str());
TrackerOfflineValidation->cd();
MyStrip->cd();
TFile fpeak(peakfile.c_str());
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
TCanvas *c5=new TCanvas("c5","",dx,dy);
c5->Divide(3,2);
TCanvas *c6=new TCanvas("c6","",dx,dy);
c6->Divide(3,2);
TCanvas *c7=new TCanvas("c7","",800,400);
c7->Divide(2,1);

const int nplots=6;

TH1* hxprime_peak[nplots];
TH1* hxprime_deco[nplots];

TH1* hwprime_peak[nplots];
TH1* hwprime_deco[nplots];

TH1* hwprime_TOB1_peak[nplots];
TH1* hwprime_TOB1_deco[nplots];

TH1* hwprime_TOB2_peak[nplots];
TH1* hwprime_TOB2_deco[nplots];

TH1* hxprime_TOB1_peak[nplots];
TH1* hxprime_TOB1_deco[nplots];

TH1* hxprime_TOB2_peak[nplots];
TH1* hxprime_TOB2_deco[nplots];

TH1* du;
TH1* dw;

char* toblayer1[]={"TOB- Layer 1",
		   "TOB- Layer 2)",
		   "TOB- Layer 3",
		   "TOB- Layer 4",
		   "TOB- Layer 5",
		   "TOB- Layer 6"};

char* toblayer2[]={"TOB+ Layer 1",
		   "TOB+ Layer 2",
		   "TOB+ Layer 3",
		   "TOB+ Layer 4",
		   "TOB+ Layer 5",
		   "TOB+ Layer 6"};


char* det[]={"TIB",
	     "TEC+",
	     "TID+",
	     "TOB",
	     "TEC-",
	     "TID-"};

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

fpeak.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();

hxprime_peak[0]=cloneHist(h_Xprime_TIBBarrel_0);
hxprime_peak[1]=cloneHist(h_Xprime_TECEndcap_4);
hxprime_peak[2]=cloneHist(h_Xprime_TIDEndcap_1);
hxprime_peak[3]=cloneHist(h_Xprime_TOBBarrel_3);
hxprime_peak[4]=cloneHist(h_Xprime_TECEndcap_5);
hxprime_peak[5]=cloneHist(h_Xprime_TIDEndcap_2);

hwprime_peak[0]=cloneHist(h_resXprimeOverTheta_TIBBarrel_0);
hwprime_peak[1]=cloneHist(h_resXprimeOverTheta_TECEndcap_4);
hwprime_peak[2]=cloneHist(h_resXprimeOverTheta_TIDEndcap_1);
hwprime_peak[3]=cloneHist(h_resXprimeOverTheta_TOBBarrel_3);
hwprime_peak[4]=cloneHist(h_resXprimeOverTheta_TECEndcap_5);
hwprime_peak[5]=cloneHist(h_resXprimeOverTheta_TIDEndcap_2);

TOBBarrel_4->cd();
TOBHalfBarrel_1->cd();

hwprime_TOB1_peak[0]=cloneHist(h_resXprimeOverTheta_TOBLayer_0);
hwprime_TOB1_peak[1]=cloneHist(h_resXprimeOverTheta_TOBLayer_1);
hwprime_TOB1_peak[2]=cloneHist(h_resXprimeOverTheta_TOBLayer_2);
hwprime_TOB1_peak[3]=cloneHist(h_resXprimeOverTheta_TOBLayer_3);
hwprime_TOB1_peak[4]=cloneHist(h_resXprimeOverTheta_TOBLayer_4);
hwprime_TOB1_peak[5]=cloneHist(h_resXprimeOverTheta_TOBLayer_5);

hxprime_TOB1_peak[0]=cloneHist(h_Xprime_TOBLayer_0);
hxprime_TOB1_peak[1]=cloneHist(h_Xprime_TOBLayer_1);
hxprime_TOB1_peak[2]=cloneHist(h_Xprime_TOBLayer_2);
hxprime_TOB1_peak[3]=cloneHist(h_Xprime_TOBLayer_3);
hxprime_TOB1_peak[4]=cloneHist(h_Xprime_TOBLayer_4);
hxprime_TOB1_peak[5]=cloneHist(h_Xprime_TOBLayer_5);

fpeak.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TOBBarrel_4->cd();
TOBHalfBarrel_2->cd();
hwprime_TOB2_peak[0]=cloneHist(h_resXprimeOverTheta_TOBLayer_0);
hwprime_TOB2_peak[1]=cloneHist(h_resXprimeOverTheta_TOBLayer_1);
hwprime_TOB2_peak[2]=cloneHist(h_resXprimeOverTheta_TOBLayer_2);
hwprime_TOB2_peak[3]=cloneHist(h_resXprimeOverTheta_TOBLayer_3);
hwprime_TOB2_peak[4]=cloneHist(h_resXprimeOverTheta_TOBLayer_4);
hwprime_TOB2_peak[5]=cloneHist(h_resXprimeOverTheta_TOBLayer_5);

hxprime_TOB2_peak[0]=cloneHist(h_Xprime_TOBLayer_0);
hxprime_TOB2_peak[1]=cloneHist(h_Xprime_TOBLayer_1);
hxprime_TOB2_peak[2]=cloneHist(h_Xprime_TOBLayer_2);
hxprime_TOB2_peak[3]=cloneHist(h_Xprime_TOBLayer_3);
hxprime_TOB2_peak[4]=cloneHist(h_Xprime_TOBLayer_4);
hxprime_TOB2_peak[5]=cloneHist(h_Xprime_TOBLayer_5);

fdeco.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();

hxprime_deco[0]=cloneHist(h_Xprime_TIBBarrel_0);
hxprime_deco[1]=cloneHist(h_Xprime_TECEndcap_4);
hxprime_deco[2]=cloneHist(h_Xprime_TIDEndcap_1);
hxprime_deco[3]=cloneHist(h_Xprime_TOBBarrel_3);
hxprime_deco[4]=cloneHist(h_Xprime_TECEndcap_5);
hxprime_deco[5]=cloneHist(h_Xprime_TIDEndcap_2);

hwprime_deco[0]=cloneHist(h_resXprimeOverTheta_TIBBarrel_0);
hwprime_deco[1]=cloneHist(h_resXprimeOverTheta_TECEndcap_4);
hwprime_deco[2]=cloneHist(h_resXprimeOverTheta_TIDEndcap_1);
hwprime_deco[3]=cloneHist(h_resXprimeOverTheta_TOBBarrel_3);
hwprime_deco[4]=cloneHist(h_resXprimeOverTheta_TECEndcap_5);
hwprime_deco[5]=cloneHist(h_resXprimeOverTheta_TIDEndcap_2);


TOBBarrel_4->cd();
TOBHalfBarrel_1->cd();

hwprime_TOB1_deco[0]=cloneHist(h_resXprimeOverTheta_TOBLayer_0);
hwprime_TOB1_deco[1]=cloneHist(h_resXprimeOverTheta_TOBLayer_1);
hwprime_TOB1_deco[2]=cloneHist(h_resXprimeOverTheta_TOBLayer_2);
hwprime_TOB1_deco[3]=cloneHist(h_resXprimeOverTheta_TOBLayer_3);
hwprime_TOB1_deco[4]=cloneHist(h_resXprimeOverTheta_TOBLayer_4);
hwprime_TOB1_deco[5]=cloneHist(h_resXprimeOverTheta_TOBLayer_5);

hxprime_TOB1_deco[0]=cloneHist(h_Xprime_TOBLayer_0);
hxprime_TOB1_deco[1]=cloneHist(h_Xprime_TOBLayer_1);
hxprime_TOB1_deco[2]=cloneHist(h_Xprime_TOBLayer_2);
hxprime_TOB1_deco[3]=cloneHist(h_Xprime_TOBLayer_3);
hxprime_TOB1_deco[4]=cloneHist(h_Xprime_TOBLayer_4);
hxprime_TOB1_deco[5]=cloneHist(h_Xprime_TOBLayer_5);

fdeco.cd();
TrackerOfflineValidation->cd();
MyStrip->cd();
TOBBarrel_4->cd();
TOBHalfBarrel_2->cd();
hwprime_TOB2_deco[0]=cloneHist(h_resXprimeOverTheta_TOBLayer_0);
hwprime_TOB2_deco[1]=cloneHist(h_resXprimeOverTheta_TOBLayer_1);
hwprime_TOB2_deco[2]=cloneHist(h_resXprimeOverTheta_TOBLayer_2);
hwprime_TOB2_deco[3]=cloneHist(h_resXprimeOverTheta_TOBLayer_3);
hwprime_TOB2_deco[4]=cloneHist(h_resXprimeOverTheta_TOBLayer_4);
hwprime_TOB2_deco[5]=cloneHist(h_resXprimeOverTheta_TOBLayer_5);

hxprime_TOB2_deco[0]=cloneHist(h_Xprime_TOBLayer_0);
hxprime_TOB2_deco[1]=cloneHist(h_Xprime_TOBLayer_1);
hxprime_TOB2_deco[2]=cloneHist(h_Xprime_TOBLayer_2);
hxprime_TOB2_deco[3]=cloneHist(h_Xprime_TOBLayer_3);
hxprime_TOB2_deco[4]=cloneHist(h_Xprime_TOBLayer_4);
hxprime_TOB2_deco[5]=cloneHist(h_Xprime_TOBLayer_5);

for(int i=0;i<6;i++){

  c1->cd(i+1);
  setStats(hxprime_peak[i],hxprime_deco[i],0.5,0.65,0.15,false);  
  plotHists(hxprime_peak[i],hxprime_deco[i],det[i],"#Delta u [#mum]",-1000,1000,1,0);
  leg1->Draw();

  c2->cd(i+1);  
  setStats(hwprime_peak[i],hwprime_deco[i],0.5,0.65,0.15,false);  
  plotHists(hwprime_peak[i],hwprime_deco[i],det[i],"#Delta w [#mum]",-1000,1000,1,0);
  leg1->Draw();

  c3->cd(i+1);
  setStats(hwprime_TOB1_peak[i],hwprime_TOB1_deco[i],0.5,0.65,0.15,false);  
  plotHists(hwprime_TOB1_peak[i],hwprime_TOB1_deco[i],toblayer1[i],"#Delta w [#mum]",-1000,1000,1,3);
  leg1->Draw();

  c4->cd(i+1);
  setStats(hwprime_TOB2_peak[i],hwprime_TOB2_deco[i],0.5,0.65,0.15,false);  
  plotHists(hwprime_TOB2_peak[i],hwprime_TOB2_deco[i],toblayer2[i],"#Delta w [#mum]",-1000,1000,1,3);
  leg1->Draw();

  c5->cd(i+1);
  setStats(hxprime_TOB1_peak[i],hxprime_TOB1_deco[i],0.5,0.65,0.15,false); 
  plotHists(hxprime_TOB1_peak[i],hxprime_TOB1_deco[i],toblayer1[i],"#Delta u [#mum]",-1000,1000,1,3);
  leg1->Draw();

  c6->cd(i+1);
  setStats(hxprime_TOB2_peak[i],hxprime_TOB2_deco[i],0.5,0.65,0.15,false); 
  plotHists(hxprime_TOB2_peak[i],hxprime_TOB2_deco[i],toblayer2[i],"#Delta u [#mum]",-1000,1000,1,3);
  leg1->Draw();
}

c7->cd(1);
setStats(hxprime_peak[3],hxprime_deco[3],0.5,0.65,0.15,false);  
plotHists(hxprime_peak[3],hxprime_deco[3],"TOB","#Delta u [#mum]",-500,500,1,3);
leg1->Draw();
c7->cd(2);
setStats(hwprime_peak[3],hwprime_deco[3],0.5,0.65,0.15,false);  
plotHists(hwprime_peak[3],hwprime_deco[3],"TOB","#Delta w [#mum]",-1000,1000,1,4);
leg1->Draw();

if(printeps){
  c1->Print("plots/deltau.eps");
  c2->Print("plots/deltaw.eps");
  c3->Print("plots/deltaw_HB1.eps");
  c4->Print("plots/deltaw_HB2.eps");
  c5->Print("plots/deltau_HB1.eps");
  c6->Print("plots/deltau_HB2.eps");
}
if(printgif){
  c1->Print("plots/deltau.gif");
  c2->Print("plots/deltaw.gif");
  c3->Print("plots/deltaw_HB1.gif");
  c4->Print("plots/deltaw_HB2.gif");
  c5->Print("plots/deltau_HB1.gif");
  c6->Print("plots/deltau_HB2.gif");
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
