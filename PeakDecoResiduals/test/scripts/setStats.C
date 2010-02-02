//void setStats(TH1* s,TH1* r, double startingY, double startingX = .1,bool fit){
void setStats(TH1* s,TH1* r, double startingY, double startingX = .1, double height = 0.15, bool fit){
  if (startingY<0){
    s->SetStats(0);
    r->SetStats(0);
  } else {
    gStyle->SetOptStat("mr");

    if (fit){
      s->Fit("gaus");
      TF1* f1 = (TF1*) s->GetListOfFunctions()->FindObject("gaus");
      f1->SetLineColor(2);
      f1->SetLineWidth(1);
    }
    s->Draw();
    gPad->Update(); 
    TPaveStats* st1 = (TPaveStats*) s->GetListOfFunctions()->FindObject("stats");
    //if (fit) {st1->SetOptFit(0010);    st1->SetOptStat(1001);}
    //st1->SetOptFit(0010);    st1->SetOptStat(1001);
    
    
    st1->SetX1NDC(startingX);
    st1->SetX2NDC(startingX+0.30);
    st1->SetY1NDC(startingY+height+0.05);
    st1->SetY2NDC(startingY+2*height+0.05);
    
    //st1->SetY1NDC(startingY+0.20);
    //st1->SetY2NDC(startingY+0.35);
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
    //st2->SetY2NDC(startingY+0.15);
    st2->SetY2NDC(startingY+height);
    st2->SetTextColor(4);
  }
}
