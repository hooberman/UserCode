void plotHists(TH1* h1, TH1* h2, 
	       char* title, char* xtitle,
	       float xmin, float xmax,int rebin,int fit){
   
  gROOT->LoadMacro("scripts/fround.C");
  if(rebin>1){
    h1->Rebin(rebin);
    h2->Rebin(rebin);
  }
  
  //h1->Sumw2();
  //h2->Sumw2();
  if(h1->Integral()>0) h1->Scale(1./h1->Integral());
  if(h2->Integral()>0) h2->Scale(1./h2->Integral());

  if(fit>0){
    if(fit==1){
      TF1 *f1=new TF1("f1","[0]*exp(-0.5*pow((x-[1])/[2],2))",xmin,xmax);
      TF1 *f2=new TF1("f2","[0]*exp(-0.5*pow((x-[1])/[2],2))",xmin,xmax);
      f1->SetParameters(h1->GetMaximum(),0,h1->GetRMS(1));
      f2->SetParameters(h2->GetMaximum(),0,h1->GetRMS(1));
    }
    if(fit==2){
      TF1 *f1=new TF1("f1",
		      "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))",
		      xmin,xmax);
      TF1 *f2=new TF1("f2",
		      "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))",
		      xmin,xmax);
      f1->SetParameters(h1->GetMaximum()/2.,0,200,h1->GetMaximum()/2.,500);
      f2->SetParameters(h2->GetMaximum()/2.,0,200,h2->GetMaximum()/2.,500);
      f1->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}");
      f2->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}");
    }
    if(fit==3){
      TF1 *f1=new TF1("f1","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",xmin,xmax);
      TF1 *f2=new TF1("f2","[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[1])/[4],2))+[5]*exp(-0.5*pow((x-[1])/[6],2))",xmin,xmax);
      f1->SetParameters(h1->GetMaximum()/3.,0,200,h1->GetMaximum()/3.,400,h1->GetMaximum()/3.,600);
      f2->SetParameters(h2->GetMaximum()/3.,0,200,h2->GetMaximum()/3.,400,h2->GetMaximum()/3.,600);
      f1->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
      f2->SetParNames("c_{1}","#bar{x}","#sigma_{1}","c_{2}","#sigma_{2}","c_{3}","#sigma_{3}");
    }

    //landau
    if(fit==4){
      TF1 *f1=new TF1("f1","landau");
      TF1 *f2=new TF1("f2","landau");
		      


    }


    f1->SetLineWidth(1);
    f2->SetLineWidth(1);
    f1->SetLineColor(2);
    f2->SetLineColor(4);
    h1->Fit(f1);
    h2->Fit(f2);
  
    //stringstream s1;
    //stringstream s2;
    //s1<<"#bar{x} "<<f1->GetParameter(1)<<"#pm"<<f1->GetParError(1)<<" #mum"<<endl;
    //s2<<"#bar{x} "<<f2->GetParameter(1)<<"#pm"<<f2->GetParError(1)<<" #mum"<<endl;
    //cout<<s1.str()<<endl;
    //cout<<s2.str()<<endl;
    //l1->DrawLatex(0.6,0.3,s1.str().c_str());
    //l1->DrawLatex(0.5,0.5,"test");
    //l1->DrawLatex(0,50,"test");


}

  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(xtitle);
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  //if(fit==0){
  h1->Draw();
  h2->Draw("sames");
  //} 
  //else{
  //  h1->Draw("E1");
  //  h2->Draw("samesE1");
  //}
  
  h1->GetXaxis()->SetRangeUser(xmin,xmax);  
  float max=h1->GetMaximum();
  if(h2->GetMaximum()>max)max=h2->GetMaximum();
  h1->SetMaximum(1.05*max);
  
  if(fit>0&&fit<4){
    stringstream s1;
    stringstream s2;
    //s1<<"#bar{x} = "<<fround(f1->GetParameter(1),1)<<"#pm"<<fround(f1->GetParError(1),1)<<endl;
    //s2<<"#bar{x} = "<<fround(f2->GetParameter(1),1)<<"#pm"<<fround(f2->GetParError(1),1)<<endl;
    s1<<"#bar{x} = "<<fround(f1->GetParameter(1),1)<<endl;
    s2<<"#bar{x} = "<<fround(f2->GetParameter(1),1)<<endl;
    TLatex *l1=new TLatex();
    TLatex *l2=new TLatex();
    l1->SetTextSize(0.06);
    l2->SetTextSize(0.06);
    l1->SetNDC();
    l2->SetNDC();
    l1->SetTextColor(2);
    l2->SetTextColor(4);
    l1->DrawLatex(0.15,0.4,s1.str().c_str());
    l2->DrawLatex(0.15,0.3,s2.str().c_str());
  }

}

