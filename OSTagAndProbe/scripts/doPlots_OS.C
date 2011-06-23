{

  gROOT->ProcessLine(".L scripts/plotSimple_OS.C+");



  //------------------------
  // electrons
  //------------------------

  //ID eff vs pt
  plotSimple_OS("output/baby_data_skim.root", "output/baby_dyee.root", TCut("pt1", "pt1"), 
   		TCut("pt2", "pt2"), 20, 0, 100,0,false,false);
  
  //iso eff vs. pt
  plotSimple_OS("output/baby_data_skim.root", "output/baby_dyee.root", TCut("pt1", "pt1"), 
   		TCut("pt2", "pt2"), 20, 0, 100,0,true,false);
  
  //------------------------
  // muons
  //------------------------

  //ID eff vs pt
  plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("pt1", "pt1"), 
   		TCut("pt2", "pt2"), 20, 0, 100,1,false,false);
  
  //iso eff vs. pt
  plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("pt1", "pt1"), 
   		TCut("pt2", "pt2"), 20, 0, 100,1,true,false);
  

  //electron iso eff vs. nvtx
  //plotSimple_OS("output/baby_data_skim.root", "output/baby_dyee.root", TCut("nvtx", "gooddavtx"), 
  //  		TCut("nvtx", "gooddavtx"), 20, 0, 20,0,true,true);

  //muon iso eff vs. nvtx
  //plotSimple_OS("output/baby_data_skim.root", "output/baby_dymm.root", TCut("nvtx", "gooddavtx"), 
  // 		TCut("nvtx", "gooddavtx"), 20, 0, 20,1,true,true);
  



}
