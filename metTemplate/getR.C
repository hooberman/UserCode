{
  
  float xmin = 76;
  float xmax = 106;

  //float xmin = 81;
  //float xmax = 101;


  TChain *ch = new TChain("T1");
  ch->Add("output/v5/lepdata_skim_baby.root");

  //ch->Draw("dilmass>>hee(100,0,200)","leptype == 0 && jetptll-ptll > -5 && jetptlt-ptlt > -5 && dilmass > 81. && dilmass < 101.");
  //ch->Draw("dilmass>>hmm(100,0,200)","leptype == 1 && dilmass > 81. && dilmass < 101.");
  
  //float nee = hee->Integral();
  //float nmm = hmm->Integral();

  ch->Draw("dilmass>>hee(200,0,200)","leptype == 0 && jetptll-ptll > -5 && jetptlt-ptlt > -5");
  ch->Draw("dilmass>>hmm(200,0,200)","leptype == 1");
  
  int bin81  = hee->FindBin( xmin );
  int bin101 = hee->FindBin( xmax ) -1;
  
  float nee = hee->Integral( bin81 , bin101 );
  float nmm = hmm->Integral( bin81 , bin101 );
  
  cout << "nee " << nee << endl;
  cout << "nmm " << nmm << endl;
  cout << "R   " << sqrt( nmm / nee ) << endl;

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  gPad->SetLogy(1);
  hee->SetLineColor(2);
  hmm->SetLineColor(4);
  hmm->Draw();
  hmm->GetXaxis()->SetTitle("M(ll) (GeV)");
  hee->Draw("same");

  TLegend *leg = new TLegend(0.65,0.6,0.85,0.8);
  leg->AddEntry(hee,"Z#rightarrowee","l");
  leg->AddEntry(hmm,"Z#rightarrow#mu#mu","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line->DrawLine( xmin , 0 , xmin , 2 * hmm->GetMaximum() );
  line->DrawLine( xmax , 0 , xmax , 2 * hmm->GetMaximum() );

  TLatex t;
  t.SetNDC();
  t.SetTextColor(2);
  t.DrawLatex(0.2,0.8,Form("N(ee) = %.0f",nee)); 
  t.SetTextColor(4);
  t.DrawLatex(0.2,0.7,Form("N(#mu#mu) = %.0f",nmm)); 
  t.SetTextColor(1);
  t.DrawLatex(0.2,0.6,Form("R = %.2f",sqrt(nmm/nee))); 






}
