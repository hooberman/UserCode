#!/bin/bash

# root -b -q "doSample_Z.C (\"data_53X_2012A\")"         > data_53X_2012A.log 2>&1 &  
## need to sleep for ~15 seconds to allow libs to compile after submitting first job
## (only if libs aren't already compiled)
# sleep 15
# root -b -q "doSample_Z.C (\"data_53X_2012B\")"         > data_53X_2012B.log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X_2012C\")"         > data_53X_2012C.log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X_2012D\")"         > data_53X_2012D.log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X_2012D_ee\")"      > data_53X_2012D_ee.log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X_2012D_mm\")"      > data_53X_2012D_mm.log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X_2012D_em\")"      > data_53X_2012D_em.log 2>&1 &  

# root -b -q "doSample_Z.C (\"ttbar_53X_slim\")"         > ttbar_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"zz2l2q_53X_slim\")"        > zz2l2q_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"zz2l2nu_53X_slim\")"       > zz2l2nu_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"zz4l_53X_slim\")"          > zz4l_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"wz2l2q_53X_slim\")"        > wz2l2q_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"wz3lnu_53X_slim\")"        > wz3lnu_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"ttw_53X_slim\")"           > ttw_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"ttz_53X_slim\")"           > ttz_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"tbz_53X_slim\")"           > tbz_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"vvv_53X_slim\")"           > vvv_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"t_53X_slim\")"             > t_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"ww_53X_slim\")"            > ww_53X_slim.log 2>&1 &  
# root -b -q "doSample_Z.C (\"zjets_small_53X_slim\")"   > zjets_small_53X_slim.log 2>&1 &  
root -b -q "doSample_Z.C (\"zjets_53X_slim\")"         > zjets_53X_slim.log 2>&1 &  

#------------------
# !!! syntax needs to be fixed for these
#------------------

# root -b -q "doSample_Z.C (\"zz2l2q_53X"              \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X"                \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_2012C_53X"          \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"data_53X_edgeSync"       \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"wz_53X"                  \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zz_53X"                  \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"gmsb"                    \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"gmsb_526"                \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"gmsb_526_v2"             \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"wzsms"                   \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"wjets_53X"               \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"ww_53X"                  \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"t_53X"                   \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"wz2l2q_53X"              \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zz2l2q_53X"              \")" > .log 2>&1 &  

# root -b -q "doSample_Z.C (\"zz4l_53X"                \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"ttW_53X"                 \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"ttZ_53X"                 \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"VVV_53X"                 \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"zjets_full_53X"          \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"ttbar_53X"               \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zjets_53X"               \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zjets_MET50_53X"         \")" > .log 2>&1 &  

# root -b -q "doSample_Z.C (\"dataskim2010"            \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"data2012cv2\")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"ttbar_massiveb"          \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"ww"                      \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"t"                       \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"wz"                      \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zz"                      \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"ttbar"                   \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zjets"                   \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"zjets_10to50"            \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"data"                    \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"dataskim"                \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"data2012c"               \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"wzsms"                   \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"gmsb"                    \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"testfilter_newJEC"       \")" > .log 2>&1 &  

# #   //--------------------------------------------------
# #   // samples for sync exercise
# #   //--------------------------------------------------

# root -b -q "doSample_Z.C (\"RelValZEE"                  \")" > .log 2>&1 &  
# root -b -q "doSample_Z.C (\"RelValZMM"                  \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"MuEG_199752"                \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"DoubleElectron_199752"      \")" > .log 2>&1 &
# root -b -q "doSample_Z.C (\"DoubleMu_199752"            \")" > .log 2>&1 &







