TH1* cloneHist(TH1* hist){

  stringstream name;
  name<<hist->GetName()<<"_clone";

  TH1F* h=new TH1F(name.str().c_str(),hist->GetTitle(),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  for(int ibin=1;ibin<=hist->GetNbinsX();ibin++){
    h->SetBinContent(ibin,hist->GetBinContent(ibin));
  }
  return h;
  
}

TH2* cloneHist(TH2* hist){

  stringstream name;
  name<<hist->GetName()<<"_clone";

  TH2F* h=new TH2F(name.str().c_str(),hist->GetTitle(),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),
		   hist->GetNbinsY(),hist->GetYaxis()->GetYmin(),hist->GetYaxis()->GetXmax());
  for(int ibinx=1;ibinx<=hist->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=hist->GetNbinsY();ibiny++){
      h->SetBinContent(ibinx,ibiny,hist->GetBinContent(ibinx,ibiny));
    }
  }
  return h;
  
}
