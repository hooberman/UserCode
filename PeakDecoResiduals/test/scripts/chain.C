bool draw(char* var, vector<char*> cantitles){
  
  for(int i=0;i<cantitles.size();i++) {
    if(strcmp(var,cantitles.at(i)) == 0) return true;
  }
  return false;
}

int getDTSlice(float t){

  int slice = -1;
  if(t > -20  && t < -12)  slice = 0;
  if(t > -12  && t < -8 )  slice = 1;
  if(t > -8   && t < -4 )  slice = 2;
  if(t > -4   && t < 0  )  slice = 3;
  if(t >  0   && t < 4  )  slice = 4;
  if(t >  4   && t < 8  )  slice = 5;
  if(t >  8   && t < 12 )  slice = 6;
  if(t >  12  && t < 20 )  slice = 7;

  return slice;
}

void chain(){
   
  enum subdetenum      { PixelBarrel = 0, PixelEndcap = 1, TIB = 2, TID = 3,
			 TOB = 4, TECthick = 5, TECthin = 6};

  int detlayers[] = {0,0,4,0,6,0,0};

  bool writeTFile   = false;
  bool printgif     = false;
  int  mysubdet     = TOB; 
  int  ndiv         = 1;
  const int nlayers = detlayers[mysubdet];


  //Configure variables to plot--------------------------
  vector<char*> cantitles;
  //cantitles.push_back("du_dw");
  //cantitles.push_back("nstripsvstantrk");
  cantitles.push_back("nstripsvstantrktgraph");
  cantitles.push_back("nstripsvstantrktgraphdt");
  cantitles.push_back("tanladt");
  //cantitles.push_back("duvsdtantheta_tgraph");
  //cantitles.push_back("duvsdtantheta_layers_tgraph");
  
  const int ncan = cantitles.size();

  TCanvas *can[ncan];
  bool drawcan[ncan];

  //drawcan[0]  = true;  //v+/- split du vs. delta tan(theta) TGraphs
  //drawcan[0]  = true;  //v+/- split du vs. delta tan(theta) TGraphs
  int idx = 0;

  TFile *ofile;
  if(writeTFile){
    ofile=TFile::Open("histos.root","RECREATE");
    ofile->cd();
  }

  //Load stuff--------------------------------------------
  gROOT->LoadMacro("scripts/fround.C");
  gROOT->LoadMacro("scripts/cloneHist.C");
  gROOT->LoadMacro("scripts/setStats.C");
  gROOT->LoadMacro("scripts/plotHists.C"); 
  gROOT->LoadMacro("scripts/plotHistsTH2.C");
  gROOT->ProcessLine(".L scripts/getTGraphFromTH2.C+");
  //gROOT->SetStyle("Plain");
  //gStyle->SetPalette(1);

  //dummy legend------------------------------------------
  TH1F *h1=new TH1F("h1","h1",1,0,1);
  h1->SetLineColor(2);
  h1->SetMarkerColor(2);
  TH1F *h2=new TH1F("h2","h2",1,0,1);
  h2->SetMarkerColor(4);

  TLegend *leg1=new TLegend(0.15,0.55,0.45,0.85);
  leg1->AddEntry(h1,"PEAK");
  leg1->AddEntry(h2,"DECO");
  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);
  
  //Set directories---------------------------------------
  char* basedir = getenv("basedir");
  char* dir     = getenv("dir");

  if(!basedir){
    cout<<"Can't find basedir. Exiting"<<endl;
    exit(0);
  }
  if(!dir){
    cout<<"Can't find dir. Exiting"<<endl;
    exit(0);
  }

  cout<<"Reading files from: "<<basedir<<"/"<<dir<<endl;

  //Declare histos-----------------------------------------
  TH1F* dupeak = new TH1F("dupeak","dupeak",500,-500,500);
  TH1F* dudeco = new TH1F("dudeco","dudeco",500,-500,500);
  TH1F* dwpeak = new TH1F("dwpeak","dwpeak",1000,-1000,1000);
  TH1F* dwdeco = new TH1F("dwdeco","dwdeco",1000,-1000,1000);
  
  TH2F* duthetapeakp = new TH2F("duthetapeakp","duthetapeakp",100,-2,2,200,-500,500);
  TH2F* duthetapeakm = new TH2F("duthetapeakm","duthetapeakm",100,-2,2,200,-500,500);
  TH2F* duthetadecop = new TH2F("duthetadecop","duthetadecop",100,-2,2,200,-500,500);
  TH2F* duthetadecom = new TH2F("duthetadecom","duthetadecom",100,-2,2,200,-500,500);

  TH2F* duthetapeakp_layer[nlayers];
  TH2F* duthetapeakm_layer[nlayers];
  TH2F* duthetadecop_layer[nlayers];
  TH2F* duthetadecom_layer[nlayers];


  TH1F* tanlapeak_layer[nlayers];
  TH1F* tanladeco_layer[nlayers];

  for(int ih = 0 ; ih < nlayers ; ih++){
    duthetapeakp_layer[ih] = new TH2F(Form("duthetapeakp_%i",ih),Form("duthetapeakp_%i",ih),100,-2,2,200,-500,500);
    duthetapeakm_layer[ih] = new TH2F(Form("duthetapeakm_%i",ih),Form("duthetapeakm_%i",ih),100,-2,2,200,-500,500);
    duthetadecop_layer[ih] = new TH2F(Form("duthetadecop_%i",ih),Form("duthetadecop_%i",ih),100,-2,2,200,-500,500);
    duthetadecom_layer[ih] = new TH2F(Form("duthetadecom_%i",ih),Form("duthetadecom_%i",ih),100,-2,2,200,-500,500);

    tanlapeak_layer[ih]    = new TH1F(Form("tanlapeak_%i",ih), Form("tanlapeak_%i",ih), 10000,0.0,0.2);
    tanladeco_layer[ih]    = new TH1F(Form("tanladeco_%i",ih), Form("tanladeco_%i",ih), 10000,0.0,0.2);
  }
  
  TH2F* nstripstantrkpeak = new TH2F("nstripstantrkpeak","nstripstantrkpeak",400,-2,2,21,-0.5,20.5);
  TH2F* nstripstantrkdeco = new TH2F("nstripstantrkdeco","nstripstantrkdeco",400,-2,2,21,-0.5,20.5);
  TH2F* nstripstantrkpeak_dt[8];
  TH2F* nstripstantrkdeco_dt[8];
  TH1F* dttimepeak_dt[8];
  TH1F* dttimedeco_dt[8];

  for(int ih = 0 ; ih < 8 ; ih++){
    nstripstantrkpeak_dt[ih]= new TH2F(Form("nstripstantrkpeak_dt_%i",ih),Form("nstripstantrkpeak_dt_%i",ih),400,-2,2,21,-0.5,20.5);
    nstripstantrkdeco_dt[ih]= new TH2F(Form("nstripstantrkdeco_dt_%i",ih),Form("nstripstantrkdeco_dt_%i",ih),400,-2,2,21,-0.5,20.5);
    dttimepeak_dt[ih]       = new TH1F(Form("dttimepeak_dt_%i",ih),Form("dttimepeak_dt_%i",ih),100,-50,50);
    dttimedeco_dt[ih]       = new TH1F(Form("dttimedeco_dt_%i",ih),Form("dttimedeco_dt_%i",ih),100,-50,50);
  }

  //Make TChains------------------------------------------
  TChain *chpeak = new TChain("PeakDecoResiduals/t","Tree");
  chpeak->Add(Form("%s/%s/root/peak.root",basedir,dir));
  int npeaktot = chpeak->GetEntries();
  int npeak = 0;

  TObjArray *listOfPeakFiles = chpeak->GetListOfFiles();
  TIter peakFileIter(listOfPeakFiles);
  TChainElement* currentPeakFile = 0;

  Float_t duP;
  Float_t dwP;
  Float_t dtanthP;
  Float_t tantrkP;
  Int_t   subdetP;
  Int_t   vP;
  Int_t   nstripsP;
  Float_t dttimeP;
  Int_t   layerP;
  Float_t tanlaP;

  while((currentPeakFile = (TChainElement*)peakFileIter.Next())) {
    TFile fpeak(currentPeakFile->GetTitle());
    TTree *tpeak = (TTree*)fpeak.Get("PeakDecoResiduals/t");
    tpeak->SetBranchAddress("du",       &duP);
    tpeak->SetBranchAddress("dw",       &dwP);
    tpeak->SetBranchAddress("dtanth",   &dtanthP);
    tpeak->SetBranchAddress("subdet",   &subdetP);
    tpeak->SetBranchAddress("v",        &vP);
    tpeak->SetBranchAddress("nstrips",  &nstripsP);
    tpeak->SetBranchAddress("tantrk" ,  &tantrkP);
    tpeak->SetBranchAddress("dttime" ,  &dttimeP);
    tpeak->SetBranchAddress("layer"  ,  &layerP);
    tpeak->SetBranchAddress("tanla"  ,  &tanlaP);
    
    unsigned int npeakentries = tpeak->GetEntries()/ndiv;
    
    
    for(unsigned int i = 0; i < npeakentries; ++i) {
      
      if(npeak % 100000 == 0) cout<<"Peak event "<<npeak<<"/"<<npeaktot<<endl;
      npeak++;

      tpeak->GetEntry(i);
      
      if(subdetP == mysubdet){
    
	dupeak->Fill(dtanthP > 0 ? duP : -duP);
	dwpeak->Fill(dwP);

	if(getDTSlice(dttimeP)>-1)
	  dttimepeak_dt[getDTSlice(dttimeP)]->Fill(dttimeP);

	if(mysubdet == TIB || mysubdet == TOB)
	  tanlapeak_layer[layerP-1]->Fill(vP > 0 ? -tanlaP : tanlaP);

	if(vP>0){
	  duthetapeakp->Fill(dtanthP,duP);
	  if(mysubdet == TIB || mysubdet == TOB)
	    duthetapeakp_layer[layerP-1]->Fill(dtanthP,duP);
	  nstripstantrkpeak->Fill(tantrkP,nstripsP);
	  if(getDTSlice(dttimeP)>-1){
            //cout<<"DTSlice "<<getDTSlice(dttimeP)<<endl;
            nstripstantrkpeak_dt[getDTSlice(dttimeP)]->Fill(tantrkP,nstripsP);
          }
	}else{
	  duthetapeakm->Fill(dtanthP,duP);
	  if(mysubdet == TIB || mysubdet == TOB)
	    duthetapeakm_layer[layerP-1]->Fill(dtanthP,duP);
	  nstripstantrkpeak->Fill(-tantrkP,nstripsP);
	  if(getDTSlice(dttimeP)>-1)
	    nstripstantrkpeak_dt[getDTSlice(dttimeP)]->Fill(-tantrkP,nstripsP);
	}


	
      }
    }
  }

  TChain *chdeco = new TChain("PeakDecoResiduals/t","Tree");
  chdeco->Add(Form("%s/%s/root/deco.root",basedir,dir));
  int ndecotot = chdeco->GetEntries();
  int ndeco = 0;

  TObjArray *listOfDecoFiles = chdeco->GetListOfFiles();
  TIter decoFileIter(listOfDecoFiles);
  TChainElement* currentDecoFile = 0;

  Float_t duD;
  Float_t dwD;
  Float_t dtanthD;
  Float_t tantrkD;
  Int_t   subdetD;
  Int_t   vD;
  Int_t   nstripsD;
  Float_t dttimeD;
  Int_t   layerD;
  Float_t tanlaD;

  while((currentDecoFile = (TChainElement*)decoFileIter.Next())) {
    TFile fdeco(currentDecoFile->GetTitle());
    TTree *tdeco = (TTree*)fdeco.Get("PeakDecoResiduals/t");
    tdeco->SetBranchAddress("du",       &duD);
    tdeco->SetBranchAddress("dw",       &dwD);
    tdeco->SetBranchAddress("dtanth",   &dtanthD);
    tdeco->SetBranchAddress("subdet",   &subdetD);
    tdeco->SetBranchAddress("v",        &vD);
    tdeco->SetBranchAddress("nstrips",  &nstripsD);
    tdeco->SetBranchAddress("tantrk" ,  &tantrkD);
    tdeco->SetBranchAddress("dttime" ,  &dttimeD);
    tdeco->SetBranchAddress("layer"  ,  &layerD);
    tdeco->SetBranchAddress("tanla"  ,  &tanlaD);

    unsigned int ndecoentries = tdeco->GetEntries()/ndiv;
    
    for(unsigned int i = 0; i < ndecoentries; ++i) {

      if(ndeco % 100000 == 0) cout<<"Deco event "<<ndeco<<"/"<<ndecotot<<endl;
      ndeco++;
      
      tdeco->GetEntry(i);
      
      if(subdetD == mysubdet){
     
	dudeco->Fill(dtanthD > 0 ? duD : -duD);
	dwdeco->Fill(dwD);

	if(getDTSlice(dttimeD)>-1)
	  dttimedeco_dt[getDTSlice(dttimeD)]->Fill(dttimeD);

	
	if(mysubdet == TIB || mysubdet == TOB)
	  tanladeco_layer[layerD-1]->Fill(vD > 0 ? -tanlaD : tanlaD);
	
	if(vD>0){
	  duthetadecop->Fill(dtanthD,duD);
	  if(mysubdet == TIB || mysubdet == TOB)
	    nstripstantrkdeco->Fill(tantrkD, nstripsD);
	  if(getDTSlice(dttimeD)>-1)
	    nstripstantrkdeco_dt[getDTSlice(dttimeD)]->Fill(tantrkD, nstripsD);
	}else{
	  duthetadecom->Fill(dtanthD,duD);
	  if(mysubdet == TIB || mysubdet == TOB)
	    duthetadecom_layer[layerD-1]->Fill(dtanthD,duD);
	  nstripstantrkdeco->Fill(-tantrkD,nstripsD);
	  if(getDTSlice(dttimeD)>-1)
	    nstripstantrkdeco_dt[getDTSlice(dttimeD)]->Fill(-tantrkD,nstripsD);
	}
	
	
	
	
      }
    }
  }
  

  //delta u, delta w TH1s
  if(draw("du_dw",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
    
    can[idx]->cd(1);
    setStats(dupeak,dudeco,0.5,0.65,0.15,false);  
    plotHists(dupeak,dudeco,"TOB","#Delta u [#mum]",-500,500,1,3);
    leg1->Draw();
    can[idx]->cd(2);
    setStats(dwpeak,dwdeco,0.5,0.65,0.15,false);  
    plotHists(dwpeak,dwdeco,"TOB","#Delta w [#mum]",-2000,2000,1,3);
    leg1->Draw();

    idx++;
  }

  if(draw("nstripsvstantrk",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
    
    plotHistsTH2(nstripstantrkpeak,
 		 "TOB PEAK","tan(#theta_{trk})",
 		 "nstrips",-2,2,0,20,0,can[idx],1);
    plotHistsTH2(nstripstantrkdeco,
 		 "TOB DECO","tan(#theta_{trk})",
 		 "nstrips",-2,2,0,20,0,can[idx],2);
    
    idx++;
  }

  if(draw("nstripsvstantrktgraph",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),800,600);
    can[idx]->cd();

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)",-1,1);
    //fnstripstantrkpeak->SetParameters(0.,1.,1.,1.,1.);
    //fnstripstantrkdeco->SetParameters(0.,1.,1.,1.,1.);

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+(x<[0])*((x-[0])*[2])+(x>0)*((x-[0])*[3])",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+(x<[0])*((x-[0])*[2])+(x>0)*((x-[0])*[3])",-1,1);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.,2.);

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+[2]*abs((x-[0]))",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+[2]*abs((x-[0]))",-1,1);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.,2.);

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","[1]+[2]*abs((x-[0]))",-1,1);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","[1]+[2]*abs((x-[0]))",-1,1);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.,2.);
    //char* func="V-fit";

    //TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak","([1]+[2]*pow(x-[0],2))",-0.5,0.3);
    //TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco","([1]+[2]*pow(x-[0],2))",-0.5,0.3);
    //fnstripstantrkpeak->SetParameters(-0.1,2.,2.);
    //fnstripstantrkdeco->SetParameters(-0.1,2.,2.);

    TF1* fnstripstantrkpeak=new TF1("fnstripstantrkpeak",
    				    "[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)+[5]*pow(x-[0],8)+[6]*pow(x-[0],10)+[7]*pow(x-[0],12)+[8]*pow(x-[0],14)",-1,1);
    TF1* fnstripstantrkdeco=new TF1("fnstripstantrkdeco",
    				    "[1]+[2]*pow(x-[0],2)+[3]*pow(x-[0],4)+[4]*pow(x-[0],6)+[5]*pow(x-[0],8)+[6]*pow(x-[0],10)+[7]*pow(x-[0],12)+[8]*pow(x-[0],14)",-1,1);
    fnstripstantrkpeak->SetParameters(0.,1.,1.,1.,1.,1.,1.,1.,1.);
    fnstripstantrkdeco->SetParameters(0.,1.,1.,1.,1.,1.,1.,1.,1.);
    char* func = "f(x) = #Sigma_{i=0}^{7} [ c_{i} (x-tan(#theta_{trk}^{MIN}))^{2i} ]";

    vector<float> bins;
    for(int ibin=0;ibin<201;ibin++) bins.push_back(ibin*0.01-1.);
    
    TGraphErrors *gnstripstantrkpeak = getTGraphFromTH2(nstripstantrkpeak,bins,0);
    fnstripstantrkpeak->SetLineColor(2);
    fnstripstantrkpeak->SetLineWidth(1);
    gnstripstantrkpeak->Fit(fnstripstantrkpeak,"R");
    gnstripstantrkpeak->SetMarkerColor(2);
    gnstripstantrkpeak->SetLineColor(2);
    gnstripstantrkpeak->SetMarkerStyle(8);
    gnstripstantrkpeak->SetMarkerSize(0.1);
    gnstripstantrkpeak->Draw("AP");
    gnstripstantrkpeak->SetTitle("TOB");
    gnstripstantrkpeak->GetXaxis()->SetTitle("<tan(#theta_{trk})>");
    gnstripstantrkpeak->GetYaxis()->SetTitle("<nstrips>");

    stringstream snstripstantrkpeak;
    snstripstantrkpeak<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkpeak->GetParameter(0),5)
		      <<" #pm "<<fround(fnstripstantrkpeak->GetParError(0),5)<<endl;

    TGraphErrors *gnstripstantrkdeco = getTGraphFromTH2(nstripstantrkdeco,bins,0);
    fnstripstantrkdeco->SetLineColor(4);
    fnstripstantrkdeco->SetLineWidth(1);
    gnstripstantrkdeco->Fit(fnstripstantrkdeco,"R");
    gnstripstantrkdeco->SetMarkerColor(4);
    gnstripstantrkdeco->SetLineColor(4);
    gnstripstantrkdeco->SetMarkerStyle(8);
    gnstripstantrkdeco->SetMarkerSize(0.1);
    gnstripstantrkdeco->Draw("sameP");

    stringstream snstripstantrkdeco;
    snstripstantrkdeco<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkdeco->GetParameter(0),5)
		      <<" #pm "<<fround(fnstripstantrkdeco->GetParError(0),5)<<endl;


    TLatex *t=new TLatex();
    t->SetNDC();
    t->SetTextSize(0.04);
    t->SetTextColor(2);
    t->DrawLatex(0.3,0.82,snstripstantrkpeak.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.3,0.75,snstripstantrkdeco.str().c_str());
    t->SetTextColor(1);
    t->DrawLatex(0.3,0.68,func);

    idx++;
  }

  float tanlapeak[8];
  float tanlaerrpeak[8];
  float tanladeco[8];
  float tanlaerrdeco[8];
  float dtpeak[8];
  float dterrpeak[8];
  float dtdeco[8];
  float dterrdeco[8];

  if(draw("nstripsvstantrktgraphdt",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,600);
    can[idx]->Divide(4,2);
    //can[idx]->cd();

    //TCanvas *c1=new TCanvas("c1","",800,600);
    //TCanvas *c2=new TCanvas("c2","",800,600);
    //c1->Divide(2,2);
    //c2->Divide(2,2);

    TF1* fnstripstantrkpeak_dt[8];
    TF1* fnstripstantrkdeco_dt[8];
    TGraphErrors *gnstripstantrkpeak_dt[8];
    TGraphErrors *gnstripstantrkdeco_dt[8];
    stringstream snstripstantrkpeak_dt[8];
    stringstream snstripstantrkdeco_dt[8];
    char* func = "V-fit";
    char* titles[]={
      "TOB (-20 ns < DT Time < -12 ns)",
      "TOB (-12 ns < DT Time < -8 ns)",
      "TOB (-8 ns < DT Time < -4 ns)",
      "TOB (-4 ns < DT Time < 0 ns)",
      "TOB (0 ns < DT Time < 4 ns)",
      "TOB (4 ns < DT Time < 8 ns)",
      "TOB (8 ns < DT Time < 12 ns)",
      "TOB (12 ns < DT Time < 20 ns)"
    };

    vector<float> bins2;
    bins2.clear();
    //for(int ibin2=0;ibin2<101;ibin2++) bins2.push_back(ibin2*0.02-1.);
    for(int ibin2=0;ibin2<201;ibin2++) bins2.push_back(ibin2*0.01-1.);

    for(int ih = 0 ; ih < 8 ; ih++ ){

      cout<<"DT Slice "<<ih<<endl;
      can[idx]->cd(ih+1);
      //if(ih<4) c1->cd(ih+1);
      //if(ih>3) c2->cd(ih-3);

      fnstripstantrkpeak_dt[ih] = new TF1(Form("fnstripstantrkpeak_dt_%i",ih),"[1]+[2]*abs((x-[0]))",-1,1);
      fnstripstantrkdeco_dt[ih] = new TF1(Form("fnstripstantrkdeco_dt_%i",ih),"[1]+[2]*abs((x-[0]))",-1,1);
      
      fnstripstantrkpeak_dt[ih]->SetParameters(-0.1,2.,2.,2.);
      fnstripstantrkdeco_dt[ih]->SetParameters(-0.1,2.,2.,2.);
      
      gnstripstantrkpeak_dt[ih] = getTGraphFromTH2(nstripstantrkpeak_dt[ih],bins2,0);
      fnstripstantrkpeak_dt[ih]->SetLineColor(2);
      fnstripstantrkpeak_dt[ih]->SetLineWidth(1);
      gnstripstantrkpeak_dt[ih]->Fit(fnstripstantrkpeak_dt[ih],"R");
      gnstripstantrkpeak_dt[ih]->SetMarkerColor(2);
      gnstripstantrkpeak_dt[ih]->SetLineColor(2);
      gnstripstantrkpeak_dt[ih]->SetMarkerStyle(8);
      gnstripstantrkpeak_dt[ih]->SetMarkerSize(0.1);
      gnstripstantrkpeak_dt[ih]->Draw("AP");
      gnstripstantrkpeak_dt[ih]->SetTitle(titles[ih]);
      gnstripstantrkpeak_dt[ih]->GetXaxis()->SetTitle("<tan(#theta_{trk})>");
      gnstripstantrkpeak_dt[ih]->GetYaxis()->SetTitle("<nstrips>");
      gnstripstantrkpeak_dt[ih]->SetMinimum(1);
      gnstripstantrkpeak_dt[ih]->SetMaximum(5);

      tanlapeak[ih]    = fnstripstantrkpeak_dt[ih]->GetParameter(0);
      tanlaerrpeak[ih] = fnstripstantrkpeak_dt[ih]->GetParError(0);
      dtpeak[ih]       = dttimepeak_dt[ih]->GetMean(1);
      dterrpeak[ih]    = dttimepeak_dt[ih]->GetRMS(1)/sqrt(dttimepeak_dt[ih]->GetEntries());

      snstripstantrkpeak_dt[ih]<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkpeak_dt[ih]->GetParameter(0),5)
			       <<" #pm "<<fround(fnstripstantrkpeak_dt[ih]->GetParError(0),5)<<endl;

      gnstripstantrkdeco_dt[ih] = getTGraphFromTH2(nstripstantrkdeco_dt[ih],bins2,0);
      fnstripstantrkdeco_dt[ih]->SetLineColor(4);
      fnstripstantrkdeco_dt[ih]->SetLineWidth(1);
      gnstripstantrkdeco_dt[ih]->Fit(fnstripstantrkdeco_dt[ih],"R");
      gnstripstantrkdeco_dt[ih]->SetMarkerColor(4);
      gnstripstantrkdeco_dt[ih]->SetLineColor(4);
      gnstripstantrkdeco_dt[ih]->SetMarkerStyle(8);
      gnstripstantrkdeco_dt[ih]->SetMarkerSize(0.1);
      gnstripstantrkdeco_dt[ih]->Draw("sameP");

      tanladeco[ih]    = fnstripstantrkdeco_dt[ih]->GetParameter(0);
      tanlaerrdeco[ih] = fnstripstantrkdeco_dt[ih]->GetParError(0);
      dtdeco[ih]       = dttimedeco_dt[ih]->GetMean(1);
      dterrdeco[ih]    = dttimedeco_dt[ih]->GetRMS(1)/sqrt(dttimedeco_dt[ih]->GetEntries());
      
      snstripstantrkdeco_dt[ih]<<"tan(#theta_{trk}^{MIN}) = "<<fround(fnstripstantrkdeco_dt[ih]->GetParameter(0),5)
			       <<" #pm "<<fround(fnstripstantrkdeco_dt[ih]->GetParError(0),5)<<endl;
      
      
      TLatex *t=new TLatex();
      t->SetNDC();
      t->SetTextSize(0.04);
      t->SetTextColor(2);
      t->DrawLatex(0.3,0.82,snstripstantrkpeak_dt[ih].str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(0.3,0.75,snstripstantrkdeco_dt[ih].str().c_str());
      t->SetTextColor(1);
      t->DrawLatex(0.3,0.68,func);
    }
    
    idx++;
  }

  if(draw("tanladt",cantitles)){
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),800,600);
    can[idx]->cd();

    TGraphErrors *glapeak=new TGraphErrors(8,dtpeak,tanlapeak,dterrpeak,tanlaerrpeak);
    TGraphErrors *gladeco=new TGraphErrors(8,dtdeco,tanladeco,dterrdeco,tanlaerrdeco);
    
    glapeak->SetMarkerColor(2);
    glapeak->SetLineColor(2);
    glapeak->SetMarkerStyle(8);
    glapeak->SetMarkerSize(1);
    glapeak->Draw("AP");
    glapeak->SetTitle("tan(#theta_{LA}) from Cluster Width");
    glapeak->GetYaxis()->SetTitle("tan(#theta_{LA})");
    glapeak->GetXaxis()->SetTitle("DT time (ns)");
    glapeak->SetMinimum(-0.11);
    gladeco->SetMaximum(-0.05);

    gladeco->SetMarkerColor(4);
    gladeco->SetLineColor(4);
    gladeco->SetMarkerStyle(8);
    gladeco->SetMarkerSize(1);
    gladeco->Draw("sameP");
   
    idx++;
  }


  if(draw("duvsdtantheta_tgraph",cantitles)){
    bool addbpcorr = false; //add tgraph for BP-corrected deco
    float size = 0.1;

    vector<float> thetabins;
//         thetabins.push_back(-1.);
//         thetabins.push_back(-0.75);
//         thetabins.push_back(-0.5);
//         thetabins.push_back(-0.25);
//         thetabins.push_back(0.);
//         thetabins.push_back(0.25);
//         thetabins.push_back(0.5);
//         thetabins.push_back(0.75);
//         thetabins.push_back(1.);

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

    
    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),1200,450);
    can[idx]->Divide(2,1);
  
    can[idx]->cd(1);
    
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
    //mgdtp->GetXaxis()->SetLimits(thetabins.at(0),thetabins.at(thetabins.size())-1);
    //mgdtp->GetXaxis()->SetLimits(-0.25,0.25);
    //mgdtp->GetYaxis()->SetLimits(-30,10);
    mgdtp->SetMinimum(-20);
    mgdtp->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    //t->DrawLatex(0.5,0.8,"y = mx + b");
    t->SetTextColor(2);
    t->DrawLatex(0.25,0.85,sduthetapeakp.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.25,0.75,sduthetadecop.str().c_str());
    if(addbpcorr){
      t->SetTextColor(6);
      t->DrawLatex(0.25,0.65,sduthetadecobp_vp.str().c_str());
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
    //mgdtm->GetXaxis()->SetLimits(thetabins.at(0),thetabins.at(thetabins.size())-1);
    //mgdtm->GetXaxis()->SetLimits(-0.25,0.25);
    //mgdtm->GetYaxis()->SetLimits(-10,30);
    mgdtm->SetMinimum(-20);
    mgdtm->SetMaximum(20);

    TLatex *t=new TLatex();
    t->SetNDC();
    //t->DrawLatex(0.5,0.8,"y = mx + b");
    t->SetTextColor(2);
    t->DrawLatex(0.25,0.85,sduthetapeakm.str().c_str());
    t->SetTextColor(4);
    t->DrawLatex(0.25,0.75,sduthetadecom.str().c_str());
    if(addbpcorr){
      t->SetTextColor(6);
      t->DrawLatex(0.25,0.65,sduthetadecobp_vm.str().c_str());
    }


    float bpeakp = fduthetapeakp->GetParameter(0);
    float bpeakm = fduthetapeakm->GetParameter(0);
    float mpeakp = fduthetapeakp->GetParameter(1);
    float mpeakm = fduthetapeakm->GetParameter(1);
    float bdecop = fduthetadecop->GetParameter(0);
    float bdecom = fduthetadecom->GetParameter(0);
    float mdecop = fduthetadecop->GetParameter(1);
    float mdecom = fduthetadecom->GetParameter(1);

    float dw     = 0.5*( (mdecop - mpeakp) + (mdecom - mpeakm) );
    float db     = 0.5*( (bdecop - bpeakp) - (bdecom - bpeakm) );
    float dtanla = db/(235-dw);
    
    cout<<"dw "<<dw<<" dtanla "<<dtanla<<endl;

    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);
    idx++;
  }


  if(draw("duvsdtantheta_layers_tgraph",cantitles)){
    
    float size        = 0.1;
   

    vector<float> thetabins;
    thetabins.clear();
    thetabins.push_back(-1.);
    thetabins.push_back(-0.75);
    thetabins.push_back(-0.5);
    thetabins.push_back(-0.25);
    thetabins.push_back(0.);
    thetabins.push_back(0.25);
    thetabins.push_back(0.5);
    thetabins.push_back(0.75);
    thetabins.push_back(1.);
    

    
    TGraphErrors *gduthetapeakp_layer[nlayers];
    TGraphErrors *gduthetapeakm_layer[nlayers];
    TGraphErrors *gduthetadecop_layer[nlayers];
    TGraphErrors *gduthetadecom_layer[nlayers];

    TF1 *fduthetapeakp_layer[nlayers];
    TF1 *fduthetapeakm_layer[nlayers];
    TF1 *fduthetadecop_layer[nlayers];
    TF1 *fduthetadecom_layer[nlayers];

    stringstream sduthetapeakp_layer[nlayers];
    stringstream sduthetapeakm_layer[nlayers];
    stringstream sduthetadecop_layer[nlayers];
    stringstream sduthetadecom_layer[nlayers];

    TLatex *t=new TLatex();
    t->SetNDC();

    TCanvas *tanlacan[nlayers];
    float tanlapeak[nlayers];
    float tanlaerrpeak[nlayers];
    float tanladeco[nlayers];
    float tanlaerrdeco[nlayers];

    float bpeakp;
    float bpeakm;
    float mpeakp;
    float mpeakm;
    float bdecop;
    float bdecom;
    float mdecop;
    float mdecom;
    float bpeakerrp;
    float bpeakerrm;
    float bdecoerrp;
    float bdecoerrm;
    float bpeak; 
    float mpeak; 
    float bdeco; 
    float mdeco;    
    
    for( int ih = 0 ; ih < nlayers ; ih++ ){

      tanlapeak[ih]    = tanlapeak_layer[ih]->GetMean(1);
      tanlaerrpeak[ih] = pow(tanlapeak_layer[ih]->GetRMS(1),2);

      tanladeco[ih]    = tanladeco_layer[ih]->GetMean(1);
      tanlaerrdeco[ih] = pow(tanladeco_layer[ih]->GetRMS(1),2);

      cout<<"Layer "<<ih<<"------------------------------------"<<endl;
      cout<<"Peak LA "<<tanlapeak[ih]<<" +/- "<<tanlaerrpeak[ih]<<endl;
      cout<<"Deco LA "<<tanladeco[ih]<<" +/- "<<tanlaerrdeco[ih]<<endl;
      
      tanlacan[ih]=new TCanvas(Form("tanlacan_%i",ih),Form("tanlacan_%i",ih),1200,450);
      tanlacan[ih]->Divide(2,1);
      tanlacan[ih]->cd(1);
      
      duthetapeakp_layer[ih]  -> SetName(Form("duthetapeakp_layer_%i",ih));
      gduthetapeakp_layer[ih] =  getTGraphFromTH2(duthetapeakp_layer[ih],thetabins,0);
      gduthetapeakp_layer[ih] -> SetMarkerColor(2);
      gduthetapeakp_layer[ih] -> SetLineColor(2);
      gduthetapeakp_layer[ih] -> SetMarkerStyle(8);
      gduthetapeakp_layer[ih] -> SetMarkerSize(size);
      fduthetapeakp_layer[ih] =  new TF1(Form("fduthetapeakp_layer_%i",i),"pol1");
      fduthetapeakp_layer[ih] -> SetLineColor(2);
      gduthetapeakp_layer[ih] -> Fit(fduthetapeakp_layer[ih]);
      
      sduthetapeakp_layer[ih]<<"m = "<<fround(fduthetapeakp_layer[ih]->GetParameter(1),1)
			    <<" #pm " <<fround(fduthetapeakp_layer[ih]->GetParError(1),1)
			    <<"   b = "<<fround(fduthetapeakp_layer[ih]->GetParameter(0),1)
			    <<" #pm "<<fround(fduthetapeakp_layer[ih]->GetParError(0),2)
			    <<endl;
    
      

      duthetadecop_layer[ih]  -> SetName(Form("duthetadecop_layer_%i",ih));
      gduthetadecop_layer[ih] =  getTGraphFromTH2(duthetadecop_layer[ih],thetabins,0);
      gduthetadecop_layer[ih] -> SetMarkerColor(4);
      gduthetadecop_layer[ih] -> SetLineColor(4);
      gduthetadecop_layer[ih] -> SetMarkerStyle(8);
      gduthetadecop_layer[ih] -> SetMarkerSize(size);
      fduthetadecop_layer[ih] =  new TF1(Form("fduthetadecop_layer_%i",ih),"pol1");
      fduthetadecop_layer[ih] -> SetLineColor(4);
      gduthetadecop_layer[ih] -> Fit(fduthetadecop_layer[ih]);

      sduthetadecop_layer[ih]<<"m = "<<fround(fduthetadecop_layer[ih]->GetParameter(1),1)
			    <<" #pm " <<fround(fduthetadecop_layer[ih]->GetParError(1),1)
			    <<"   b = "<<fround(fduthetadecop_layer[ih]->GetParameter(0),1)
			    <<" #pm "<<fround(fduthetadecop_layer[ih]->GetParError(0),2)
			    <<endl;

      


      TMultiGraph *mgdtp=new TMultiGraph();
      mgdtp->Add(gduthetapeakp_layer[ih]);
      mgdtp->Add(gduthetadecop_layer[ih]);
      
      mgdtp->SetTitle("TOB (v+)");
      mgdtp->Draw("AP");
      mgdtp->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      mgdtp->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      mgdtp->SetTitle("TOB");
      mgdtp->GetXaxis()->SetLimits(-1.,1.);
      mgdtp->GetYaxis()->SetLimits(-20,20);
      
      
      t->SetTextColor(2);
      t->DrawLatex(0.25,0.85,sduthetapeakp_layer[ih].str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(0.25,0.75,sduthetadecop_layer[ih].str().c_str());
      
    
      tanlacan[ih]->cd(2);

    

      duthetapeakm_layer[ih]  -> SetName(Form("duthetapeakm_layer_%i",ih));
      gduthetapeakm_layer[ih] =  getTGraphFromTH2(duthetapeakm_layer[ih],thetabins,0);
      gduthetapeakm_layer[ih] -> SetMarkerColor(2);
      gduthetapeakm_layer[ih] -> SetLineColor(2);
      gduthetapeakm_layer[ih] -> SetMarkerStyle(8);
      gduthetapeakm_layer[ih] -> SetMarkerSize(size);
      fduthetapeakm_layer[ih] =  new TF1(Form("fduthetapeakm_layer_%i",ih),"pol1");
      fduthetapeakm_layer[ih] -> SetLineColor(2);
      gduthetapeakm_layer[ih] -> Fit(fduthetapeakm_layer[ih]);

      sduthetapeakm_layer[ih]<<"m = "<<fround(fduthetapeakm_layer[ih]->GetParameter(1),1)
			    <<" #pm " <<fround(fduthetapeakm_layer[ih]->GetParError(1),1)
			    <<"   b = "<<fround(fduthetapeakm_layer[ih]->GetParameter(0),1)
			    <<" #pm "<<fround(fduthetapeakm_layer[ih]->GetParError(0),2)
			    <<endl;

      duthetadecom_layer[ih]  -> SetName(Form("duthetadecom_layer_%i",ih));
      gduthetadecom_layer[ih] =  getTGraphFromTH2(duthetadecom_layer[ih],thetabins,0);
      gduthetadecom_layer[ih] -> SetMarkerColor(4);
      gduthetadecom_layer[ih] -> SetLineColor(4);
      gduthetadecom_layer[ih] -> SetMarkerStyle(8);
      gduthetadecom_layer[ih] -> SetMarkerSize(size);
      fduthetadecom_layer[ih] =  new TF1(Form("fduthetadecom_layer_%i",ih),"pol1");
      fduthetadecom_layer[ih] -> SetLineColor(4);
      gduthetadecom_layer[ih] -> Fit(fduthetadecom_layer[ih]);
      
      sduthetadecom_layer[ih]<<"m = "<<fround(fduthetadecom_layer[ih]->GetParameter(1),1)
			    <<" #pm " <<fround(fduthetadecom_layer[ih]->GetParError(1),1)
			    <<"   b = "<<fround(fduthetadecom_layer[ih]->GetParameter(0),1)
			    <<" #pm "<<fround(fduthetadecom_layer[ih]->GetParError(0),2)
			    <<endl;
      
      TMultiGraph *mgdtm=new TMultiGraph();
      mgdtm->Add(gduthetapeakm_layer[ih]);
      mgdtm->Add(gduthetadecom_layer[ih]);
      mgdtm->SetTitle("TOB (v-)");
      mgdtm->Draw("AP");
      mgdtm->GetXaxis()->SetTitle("<tan(#theta_{trk})-tan(#theta_{LA})>");
      mgdtm->GetYaxis()->SetTitle("<#Deltau> [#mum]");
      mgdtm->SetTitle("TOB");
      mgdtm->GetXaxis()->SetLimits(-1.,1.);
      mgdtm->GetYaxis()->SetLimits(-20,20);
      
      t->SetTextColor(2);
      t->DrawLatex(0.25,0.85,sduthetapeakm_layer[ih].str().c_str());
      t->SetTextColor(4);
      t->DrawLatex(0.25,0.75,sduthetadecom_layer[ih].str().c_str());

      //calculate LA
      bpeakp = fduthetapeakp_layer[ih]->GetParameter(0);
      bpeakm = fduthetapeakm_layer[ih]->GetParameter(0);
      mpeakp = fduthetapeakp_layer[ih]->GetParameter(1);
      mpeakm = fduthetapeakm_layer[ih]->GetParameter(1);
      bdecop = fduthetadecop_layer[ih]->GetParameter(0);
      bdecom = fduthetadecom_layer[ih]->GetParameter(0);
      mdecop = fduthetadecop_layer[ih]->GetParameter(1);
      mdecom = fduthetadecom_layer[ih]->GetParameter(1);

      bpeakerrp = fduthetapeakp_layer[ih]->GetParError(0);
      bpeakerrm = fduthetapeakm_layer[ih]->GetParError(0);
      bdecoerrp = fduthetadecop_layer[ih]->GetParError(0);
      bdecoerrm = fduthetadecom_layer[ih]->GetParError(0);
      
      bpeak    = 0.5*(bpeakp - bpeakm);
      mpeak    = 0.5*(mpeakp + mpeakm);
      bdeco    = 0.5*(bdecop - bdecom);
      mdeco    = 0.5*(mpeakp + mpeakm);
      bpeakerr = 0.5*(pow(bpeakerrp,2) + pow(bpeakerrm,2));
      bdecoerr = 0.5*(pow(bdecoerrp,2) + pow(bdecoerrm,2));
      
      tanlapeak[ih]    -= bpeak/(250. - mpeak);
      tanlaerrpeak[ih] += pow(bpeakerr/(250.-mpeak),2);
      tanladeco[ih]    -= bdeco/(250. - mdeco);
      tanlaerrdeco[ih] += pow(bdecoerr/(250.-mdeco),2);
    

      tanlaerrpeak[ih]  = sqrt(tanlaerrpeak[ih]);
      tanlaerrdeco[ih]  = sqrt(tanlaerrdeco[ih]);
    }
    
    float layers[nlayers];
    float layerserr[nlayers];

    for(int ih = 0 ; ih < nlayers ; ih++){
      layers[ih]    = ih+1;
      layerserr[ih] = 0.;
    }

    can[idx]=new TCanvas(cantitles.at(idx),cantitles.at(idx),800,600);
    can[idx]->cd();
    
    TGraphErrors *gtanlapeak=new TGraphErrors(nlayers,layers,tanlapeak,layerserr,tanlaerrpeak);
    gtanlapeak->SetTitle("TOB");
    gtanlapeak->GetYaxis()->SetTitle("tan(#theta_{LA})");
    gtanlapeak->GetXaxis()->SetTitle("Layer");
    gtanlapeak -> SetMarkerColor(2);
    gtanlapeak -> SetLineColor(2);
    gtanlapeak -> SetMarkerStyle(8);
    gtanlapeak -> SetMarkerSize(0.5);
    gtanlapeak -> SetMinimum(0.06);
    gtanlapeak -> SetMaximum(0.1);

    TGraphErrors *gtanladeco=new TGraphErrors(nlayers,layers,tanladeco,layerserr,tanlaerrdeco);
    gtanladeco -> SetMarkerColor(4);
    gtanladeco -> SetLineColor(4);
    gtanladeco -> SetMarkerStyle(8);
    gtanladeco -> SetMarkerSize(0.5);

    gtanlapeak->Draw("AP");
    gtanladeco->Draw("sameP");
    //line->DrawLine(-1,0,1,0);
    //line->DrawLine(0,-20,0,20);
    idx++;


  


}
  

  
  if(printgif){
    for(int ican=0;ican<cantitles.size();ican++){
      can[ican]->Modified();
      can[ican]->Update();
      can[ican]->Print(Form("plots/%s.gif",cantitles.at(ican)));
    }
  }
  
  if(writeTFile){
    ofile->cd();
    ofile->Write();
    ofile->Close();
  }
}
