{
  
  //float xmin = 76;
  //float xmax = 106;

  float xmin = 81;
  float xmax = 101;

  TChain *ch = new TChain("T1");
  ch->Add("output/v5/TTbar_baby.root");

  ch->Draw("TMath::Min(dilmass,299.9)>>h(300,0,300)","njets>1");
  
  int bin81  = h->FindBin( xmin );
  int bin101 = h->FindBin( xmax ) -1;
  
  float nin  = h->Integral( bin81 , bin101 );
  float ntot = h->Integral();

  cout << "n(in)  " << nin << endl;
  cout << "n(tot) " << ntot << endl;
  cout << "K      " << nin / ntot << endl;

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  h->SetLineColor(2);
  h->GetXaxis()->SetTitle("M(ll) (GeV)");
  h->Draw();

//   TLegend *leg = new TLegend(0.65,0.6,0.85,0.8);
//   leg->AddEntry(hee,"Z#rightarrowee","l");
//   leg->AddEntry(hmm,"Z#rightarrow#mu#mu","l");
//   leg->SetFillColor(0);
//   leg->SetBorderSize(1);
//   leg->Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line->DrawLine( xmin , 0 , xmin , 1.05 * h->GetMaximum() );
  line->DrawLine( xmax , 0 , xmax , 1.05 * h->GetMaximum() );

   TLatex t;
   t.SetNDC();
   t.SetTextColor(2);
   t.DrawLatex(0.5,0.8,Form("N(in)  = %.0f",nin)); 
   t.DrawLatex(0.5,0.7,Form("N(tot) = %.0f",ntot)); 
   t.DrawLatex(0.5,0.6,Form("K      = %.3f",nin/ntot)); 






}
