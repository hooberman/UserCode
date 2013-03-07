#!/bin/bash

DIR=$1

echo "Merging SingleMu2012AFull_${DIR}"
rm SingleMu2012AFull_${DIR}/res/merged.root
./mergeFiles.sh SingleMu2012AFull_${DIR}/res/

echo "Merging SingleMu2012BFull_${DIR}"
rm SingleMu2012BFull_${DIR}/res/merged.root
./mergeFiles.sh SingleMu2012BFull_${DIR}/res/

echo "Merging SingleMu2012CFull_${DIR}"
rm SingleMu2012CFull_${DIR}/res/merged.root
./mergeFiles.sh SingleMu2012CFull_${DIR}/res/

echo "Merging SingleMu2012DFull_${DIR}"
rm SingleMu2012DFull_${DIR}/res/merged.root
./mergeFiles.sh SingleMu2012DFull_${DIR}/res/

echo "Merging SingleEl2012AFull_${DIR}"
rm SingleEl2012AFull_${DIR}/res/merged.root
./mergeFiles.sh SingleEl2012AFull_${DIR}/res/

echo "Merging SingleEl2012BFull_${DIR}"
rm SingleEl2012BFull_${DIR}/res/merged.root
./mergeFiles.sh SingleEl2012BFull_${DIR}/res/

echo "Merging SingleEl2012CFull_${DIR}"
rm SingleEl2012CFull_${DIR}/res/merged.root
./mergeFiles.sh SingleEl2012CFull_${DIR}/res/

echo "Merging SingleEl2012DFull_${DIR}"
rm SingleEl2012DFull_${DIR}/res/merged.root
./mergeFiles.sh SingleEl2012DFull_${DIR}/res/

echo "Merging ZJetsFull_${DIR}"
rm ZJetsFull_${DIR}/res/merged.root
./mergeFiles.sh ZJetsFull_${DIR}/res/

