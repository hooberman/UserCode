{

char* basedir = getenv("basedir");
char* dir     = getenv("dir");
//char* peakdir = getenv("peakdir");
//char* decodir = getenv("decodir");

if(!basedir){
  cout<<"Can't find basedir. Exiting"<<endl;
  exit(0);
}
if(!dir){
  cout<<"Can't find dir. Exiting"<<endl;
  exit(0);
}
// if(!peakdir){
//   cout<<"Can't find peakdir. Exiting"<<endl;
//   exit(0);
// }
// if(!decodir){
//   cout<<"Can't find decodir. Exiting"<<endl;
//   exit(0);
// }

// cout<<"Peak directory "<<peakdir<<endl;
// cout<<"Deco directory "<<decodir<<endl;
cout<<"Directory "<<basedir<<"/"<<dir<<endl;

TChain *p = new TChain("PeakDecoResiduals/t","Tree");
//p->Add(Form("%s/res/*.root",peakdir));
p->Add(Form("%s/%s/root/peak.root",basedir,dir));

TChain *d = new TChain("PeakDecoResiduals/t","Tree");
//d->Add(Form("%s/res/*.root",decodir));
d->Add(Form("%s/%s/root/deco.root",basedir,dir));

}
