{
  
  char* peakdir = getenv("peakdir");
  char* decodir = getenv("decodir");
  
  if(!peakdir){
    cout<<"Can't find peakdir. Exiting"<<endl;
    exit(0);
  }
  if(!decodir){
    cout<<"Can't find decodir. Exiting"<<endl;
    exit(0);
  }
  
  cout<<"Peak directory "<<peakdir<<endl;
  cout<<"Deco directory "<<decodir<<endl;
  
  TChain *p = new TChain("PeakDecoResiduals/t","Tree");
  p->Add(Form("%s/res/*.root",peakdir));
  
  TChain *d = new TChain("PeakDecoResiduals/t","Tree");
  d->Add(Form("%s/res/*.root",decodir));
  
}
