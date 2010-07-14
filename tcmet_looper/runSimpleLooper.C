#include "TChain.h"
#include "simpleLooper.C"

void runTCMETLooperTemplate(char* prefix){

  TChain* ch = new TChain("Events");
  bool isData = false;
  
  if( strcmp( prefix , "data" )  == 0 ){  

    ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_1_1_58p.root");
    // ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_127_1_JWK.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_131_1_rU2.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_142_1_fXZ.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_143_1_esF.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_146_1_nj5.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_148_1_hk8.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_151_1_0Ym.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_153_1_WjM.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_154_1_vRB.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_155_1_XG8.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_156_1_Yl3.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_157_1_aCr.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_158_1_8eR.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_159_1_f9n.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_160_1_7ng.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_165_1_sUo.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_167_1_qJf.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_168_1_CDZ.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_172_1_VZ6.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_173_1_baX.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_174_1_Eoe.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_175_1_3wL.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_176_1_NSQ.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_179_1_zlb.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_189_1_kz6.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_49_1_Ysb.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_50_1_pUX.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_5_1_2Wp.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_53_1_bfa.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_54_1_4lm.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_6_1_ryl.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_69_1_BSl.root");
//     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_70_1_bO3.root");

    //ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_71_1_hoQ.root");
    //ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_7_1_9RG.root");
    //ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_84_1_oLb.root");
    //ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_85_1_jTE.root");
    //ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_86_1_00o.root");
    //ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_87_1_57k.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_88_1_unD.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_89_1_CXM.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_90_1_3ev.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_91_1_t76.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_93_1_EOC.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_94_1_IQL.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_95_1_HEq.root");
    //     ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_96_1_B7B.root");
    //    ch->Add("/hadoop/cms/store/user/benhoob/CMS2_V03-05-02/MinimumBias_Commissioning10-GOODCOLL-Jun9thSkim_v1/ntuple_98_1_aIV.root");
  }   														       
			
  else{
    cout << "UNRECOGNIZED SAMPLE " << prefix << ", QUITTING" << endl;
    exit(0);
  }
  
  simpleLooper* looper = new simpleLooper();
  
  cout << "Running on sample " << prefix << endl;
  looper->ScanChain(ch, prefix, isData);
  
}





