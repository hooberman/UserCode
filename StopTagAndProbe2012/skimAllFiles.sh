#!/bin/bash

DIR=$1

root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleMu2012AFull_${DIR}/res\",\"merged_json.root\"\)
root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleMu2012BFull_${DIR}/res\",\"merged_json.root\"\)
root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleMu2012CFull_${DIR}/res\",\"merged_json.root\"\)
root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleMu2012DFull_${DIR}/res\",\"merged_json.root\"\)

root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleEl2012AFull_${DIR}/res\",\"merged_json.root\"\)
root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleEl2012BFull_${DIR}/res\",\"merged_json.root\"\)
root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleEl2012CFull_${DIR}/res\",\"merged_json.root\"\)
root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"SingleEl2012DFull_${DIR}/res\",\"merged_json.root\"\)

#root -b -q skimLeptonTree.C+\(\"Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt\",\"ZJetsFull_${DIR}/res\",\"merged_json.root\"\)
