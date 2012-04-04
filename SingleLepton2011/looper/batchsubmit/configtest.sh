#!/bin/bash

TAG="V00-04-00"
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_SingleMu2011A_PrompRecoV4V33

echo "#!/bin/bash" > submitAll.sh
echo "source /code/osgcode/ucsdt2/gLite32/etc/profile.d/grid_env.sh " >> submitAll.sh
echo "voms-proxy-init -voms cms" >> submitAll.sh
echo "condor_submit condor_${TAG}_SingleMu2011A_PrompRecoV4V33.cmd" >> submitAll.sh
chmod +x submitAll.sh
echo "[writeAllConfig] wrote submit script submitAll.sh"