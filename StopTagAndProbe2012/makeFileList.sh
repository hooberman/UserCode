#!/bin/bash

find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleMu_Run2012A-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleMu_2012A_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleMu_Run2012B-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleMu_2012B_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleMu_Run2012C-PromptReco-v2_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleMu_2012C_list.txt

find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleElectron_Run2012A-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleEl_2012A_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleElectron_Run2012B-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleEl_2012B_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleElectron_Run2012C-PromptReco-v2_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleEl_2012C_list.txt

