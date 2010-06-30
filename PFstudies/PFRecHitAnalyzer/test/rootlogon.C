{
cout<<"Ben's rootlogon.C"<<endl;
gSystem->Load("libGraf");
gSystem->Load("libGpad");
gSystem->Load("libTree");

//gROOT->ProcessLine(".L /tas03/home/benhoob/benstyle.C");
//setBenStyle();

gROOT->ProcessLine(".L /tas03/home/benhoob/.root/tdrstyle.C");
setTDRStyle();

}
