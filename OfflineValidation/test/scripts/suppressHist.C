TH2* suppressHist(TH2* hist,float xmin,float xmax){

  stringstream name;
  name<<hist->GetName()<<"_"<<xmin<<"_clone";

  TH2F* h=new TH2F(name.str().c_str(),hist->GetTitle(),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),
		   hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());

  //cout<<"("<<h->GetXaxis()->GetXmin()<<","<<h->GetXaxis()->GetXmax()<<")"<<endl;
  //cout<<"("<<h->GetYaxis()->GetXmin()<<","<<h->GetYaxis()->GetXmax()<<")"<<endl;
 
  for(int ibinx=1;ibinx<=hist->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=hist->GetNbinsY();ibiny++){
      if(hist->GetBinCenter(ibinx)>xmin && hist->GetBinCenter(ibinx)<xmax){
	h->SetBinContent(ibinx,ibiny,hist->GetBinContent(ibinx,ibiny));
      }else{
	h->SetBinContent(ibinx,ibiny,0);
      }
    }
  }
  return h;






}
