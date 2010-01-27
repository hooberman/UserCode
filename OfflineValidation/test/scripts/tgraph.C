{
gROOT->SetStyle("Plain");

int x[]={1,2,3,4,5,6,7,8,9,10,11,12};
int y[]={23,25,-26,-24,-31,-24,25,24,25,21,-26,-20};

TGraph *gr=new TGraph(12,x,y);
gr->SetTitle("");
gr->Draw("AP");
gr->SetMarkerStyle(8);
gr->SetMarkerSize(2);
gr->SetMarkerColor(2);
gr->GetYaxis()->SetTitle("<#Delta_{W}> (#mum)");
gr->GetXaxis()->SetTitle("Ring");
gr->GetXaxis()->SetTitleSize(0.06);
gr->GetXaxis()->SetTitleOffset(0.7);
gr->GetYaxis()->SetTitleSize(0.06);
gr->GetYaxis()->SetTitleOffset(0.7);
}
