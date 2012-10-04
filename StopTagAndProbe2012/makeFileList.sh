#!/bin/bash

find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleMu_Run2012A-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleMu_2012A_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleMu_Run2012B-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleMu_2012B_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleMu_Run2012C-PromptReco-v2_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleMu_2012C_list.txt

find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleElectron_Run2012A-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleEl_2012A_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleElectron_Run2012B-13Jul2012-v1_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleEl_2012B_list.txt
find /hadoop/cms/store/user/benhoob/StopTNPSkim2012/SingleElectron_Run2012C-PromptReco-v2_AOD/V05-03-13/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > SingleEl_2012C_list.txt

cat SingleMu_2012A_list.txt SingleMu_2012B_list.txt SingleMu_2012C_list.txt SingleEl_2012A_list.txt SingleEl_2012B_list.txt SingleEl_2012C_list.txt > SingleLepton_all_list.txt

find /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/SingleOrDiLepton/  -type f -print | sed 's/\/hadoop\/cms/root:\/\/xrootd.unl.edu\//' > ZJets_list.txt