{
cout<<"Ben's rootlogon.C"<<endl;
gSystem->Load("libGraf");
gSystem->Load("libGpad");
gSystem->Load("libTree");

gROOT->ProcessLine(".L /tas03/home/benhoob/.root/benstyle.C");
setBenStyle();
//gROOT->SetStyle("Plain");
}
