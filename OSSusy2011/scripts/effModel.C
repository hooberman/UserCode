{

  const unsigned int n = 3;
  char* samples[n]={"LM1","LM3","LM6"};


  TH1F* h = new TH1F("h","h",1,0,1);
  h->Sumw2();

  for(unsigned int i = 0 ; i < n ; ++i ){


    TChain *chreco = new TChain("t");
    chreco->Add(Form("../output/V00-02-15/highpt/%sv2_smallTree.root",samples[i]));

    TChain *chgen = new TChain("t");
    chgen->Add(Form("../output/V00-02-18/highpt/%sv2_smallTree_gen.root",samples[i]));


    TCut sel("!passz");
    TCut highmet("pfmet>275 && ht>300");
    TCut highht ("pfmet>200 && ht>600");
    TCut tight  ("pfmet>275 && ht>600");

    TCut recoweight("trgeff * weight * 4.7");
    //TCut recoweight("weight");

    cout << endl << endl;


    // dilepton only
    chreco->Draw("0.5>>h",sel*recoweight);
    //cout << "Dilepton only : reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;

    float dilreco    = h->GetBinContent(1);
    float dilrecoerr = h->GetBinError(1);

    chgen->Draw("0.5>>h","geff * weight * 4.7");
    //cout << "Dilepton only : gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl << endl;

    float dilgen    = h->GetBinContent(1);
    float dilgenerr = h->GetBinError(1);

    // high MET
    chreco->Draw("0.5>>h",(sel+highmet)*recoweight);
    //cout << "HighMET       : reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;

    float metreco    = h->GetBinContent(1);
    float metrecoerr = h->GetBinError(1);

    chgen->Draw("0.5>>h","geffmet * weight * 4.7");
    //cout << "HighMET       : gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl << endl;

    float metgen    = h->GetBinContent(1);
    float metgenerr = h->GetBinError(1);

    // high HT
    chreco->Draw("0.5>>h",(sel+highht)*recoweight);
    //cout << "HighHT        : reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;

    float htreco    = h->GetBinContent(1);
    float htrecoerr = h->GetBinError(1);

    chgen->Draw("0.5>>h","geffht * weight * 4.7");
    //cout << "HighHT        : gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl << endl;

    float htgen    = h->GetBinContent(1);
    float htgenerr = h->GetBinError(1);

    // tight
    chreco->Draw("0.5>>h",(sel+tight)*recoweight);
    //cout << "Tight         : reco  " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl;

    float tightreco    = h->GetBinContent(1);
    float tightrecoerr = h->GetBinError(1);

    chgen->Draw("0.5>>h","gefftight * weight * 4.7");
    //cout << "Tight         : gen   " << Form("$%.0f \\pm %.1f$ ",h->Integral(),h->GetBinError(1)) << endl << endl;

    float tightgen    = h->GetBinContent(1);
    float tightgenerr = h->GetBinError(1);

    cout << endl;
    cout << samples[i] << "(reco)     & " 
	 << Form("$%.0f \\pm %.1f$",dilreco,dilrecoerr) << " & "  
	 << Form("$%.0f \\pm %.1f$",metreco,metrecoerr) << " & " 
	 << Form("$%.0f \\pm %.1f$",htreco,htrecoerr) << " & " 
	 << Form("$%.0f \\pm %.1f$",tightreco,tightrecoerr) << " \\\\ " << endl;

    cout << samples[i] << "(gen)      & " 
	 << Form("$%.0f \\pm %.1f$",dilgen,dilgenerr) << " & "  
	 << Form("$%.0f \\pm %.1f$",metgen,metgenerr) << " & " 
	 << Form("$%.0f \\pm %.1f$",htgen,htgenerr) << " & " 
	 << Form("$%.0f \\pm %.1f$",tightgen,tightgenerr) << " \\\\ " << endl;

    float dilratio    = dilgen/dilreco;
    float dilratioerr = dilratio * sqrt( pow(dilgenerr/dilgen,2) + pow(dilrecoerr/dilreco,2) );

    float metratio    = metgen/metreco;
    float metratioerr = metratio * sqrt( pow(metgenerr/metgen,2) + pow(metrecoerr/metreco,2) );

    float htratio    = htgen/htreco;
    float htratioerr = htratio * sqrt( pow(htgenerr/htgen,2) + pow(htrecoerr/htreco,2) );

    float tightratio    = tightgen/tightreco;
    float tightratioerr = tightratio * sqrt( pow(tightgenerr/tightgen,2) + pow(tightrecoerr/tightreco,2) );

    cout << samples[i] << "(gen/reco) & " 
	 << Form("$%.2f \\pm %.2f$",dilratio,dilratioerr) << " & "  
	 << Form("$%.2f \\pm %.2f$",metratio,metratioerr) << " & " 
	 << Form("$%.2f \\pm %.2f$",htratio,htratioerr) << " & " 
	 << Form("$%.2f \\pm %.2f$",tightratio,tightratioerr) << " \\\\ " << endl;


  }
}
