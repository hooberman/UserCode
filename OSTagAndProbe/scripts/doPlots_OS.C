{

  gROOT->ProcessLine(".L scripts/plotSimple_OS.C+");



  //------------------------
  // electrons
  //------------------------

  // //------------------------
  // // muons
  // //------------------------

  //ID eff vs nvtx
  plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("pt1", "pt1"), 
		TCut("pt2", "pt2"), 20, 0, 100,1,false,false);

  //iso eff vs. nvtx
  plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("pt1", "pt1"), 
		TCut("pt2", "pt2"), 20, 0, 100,1,true,false);

  // //ID eff vs nvtx
  // plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("nvtx", "gooddavtx"), 
  // 		TCut("nvtx", "gooddavtx"), 20, 0, 20,1,false,false);
  // //iso eff vs. nvtx
  // plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("nvtx", "gooddavtx"), 
  // 		TCut("nvtx", "gooddavtx"), 20, 0, 20,1,true,false);




}
