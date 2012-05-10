{

  TChain *photon = new TChain("T1");
  photon->Add("../templates/V00-00-01/Photon_baby.root");

  TCut clean("abs(etag)<2 && jetpt-etg>-5 && pfjetid==1");
  

}
