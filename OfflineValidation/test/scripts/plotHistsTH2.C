void plotHistsTH2(TH2* h,char* title, char* xtitle, char* ytitle, float xmin, float xmax, float ymin, float ymax,int fit,TCanvas* c, int ipad){

  c->cd(ipad);
  TPad *plotpad=new TPad("plotpad","",0.1,0.1,1.,0.9);
  plotpad->Draw();
  plotpad->cd();

  if(fit>0){
    TF1 *f=new TF1("f","[0]+[1]*x",xmin,xmax);
    f->SetParNames("y-int","slope");
    f->SetLineWidth(2);
    f->SetLineColor(2);
    h->Fit(f);
  }
  
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->SetTitle("");
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("");
  h->Draw("colz");



  c->cd(ipad);
  TPad *xpad=new TPad("xpad","",0.,0.,1.,0.1);
  xpad->Draw();
  xpad->cd();
  TLatex *x=new TLatex();
  x->SetNDC();
  x->SetTextSize(0.8);
  x->DrawLatex(0.5,0.3,xtitle);

  c->cd(ipad);
  TPad *ypad=new TPad("ypad","",0.,0.,0.1,1.);
  ypad->Draw();
  ypad->cd();
  TLatex *y=new TLatex();
  y->SetNDC();
  y->SetTextSize(0.8);
  y->SetTextAngle(90);
  y->DrawLatex(0.6,0.5,ytitle);

  c->cd(ipad);
  TPad *tpad=new TPad("tpad","",0.,0.9,1.,1.);
  tpad->Draw();
  tpad->cd();
  TLatex *t=new TLatex();
  t->SetNDC();
  t->SetTextSize(0.6);
  t->DrawLatex(0.2,0.1,title);

  TLatex *l=new TLatex();
  l->SetNDC();
  l->SetTextSize(0.6);
  l->SetTextColor(2);
  float cor=h->GetCorrelationFactor(1,2);
  stringstream s;
  s<<"C="<<(int)(1000*cor)<<"#times10^{-3}"<<endl;
  l->DrawLatex(0.7,0.1,s.str().c_str());
}
