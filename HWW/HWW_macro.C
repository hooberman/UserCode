{

  TChain* ww = new TChain("Events");
  ww->Add("WWTo2L2Nu_PU_testFinal_baby.root");

  TChain* hww = new TChain("Events");
  hww->Add("HToWWTo2L2NuM130_PU_testFinal_baby.root");

  TCut weight = "event_scale1fb_";
  TCut met_projpt = "(event_type!=2&met_projpt>35)|(event_type==2&met_projpt>20)";
  //all cuts in 2010 analysis h->ww m=130
  TCut cutsh130 = "lepsoft_pt>20&lephard_pt>25&dil_dphi<1.05&dil_mass<45.&jets_num==0&extralep_num==0&lowptbtags_num==0"+met_projpt;

  TCut cutsh130_presel = "jets_num==0&extralep_num==0&lowptbtags_num==0"+met_projpt;
  TCut wwcuts = "lepsoft_pt>20&lephard_pt>25&dil_dphi<1.05&dil_mass<45.";

  ww->Draw("0.5>>htotbkg(1,0,1)",cutsh130_presel*weight);
  hww->Draw("0.5>>htotsig(1,0,1)",cutsh130_presel*weight);

  ww->Draw("0.5>>hpassbkg(1,0,1)",(cutsh130_presel+wwcuts)*weight);
  hww->Draw("0.5>>hpasssig(1,0,1)",(cutsh130_presel+wwcuts)*weight);

  cout << "total WW      " << htotbkg->Integral() << endl;
  cout << "total H->WW   " << htotsig->Integral() << endl;
  cout << "pass WW       " << hpassbkg->Integral() << endl;
  cout << "pass H->WW    " << hpasssig->Integral() << endl;
  cout << "H->WW eff     " << hpasssig->Integral() / htotsig->Integral() << endl;
  cout << "WW eff        " << hpassbkg->Integral() / htotbkg->Integral() << endl;

}
