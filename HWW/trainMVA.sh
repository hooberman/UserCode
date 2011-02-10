#!/bin/bash -f

export dir=$1
echo "Training directory" $dir

mkdir Trainings/$dir
cp TMVA_HWW.C Trainings/$dir/.

rm weights/*xml
rm weights/*C

#train subset of MVA's
#root -l TMVA_HWW.C\(\"Fisher\"\) > Trainings/$dir/TMVA_HWW.log 2>&1

#train all MVA's
root -l TMVA_HWW.C > Trainings/$dir/TMVA_HWW.log 2>&1
#root -l TMVA_HWW.C

cp -r weights/ Trainings/$dir/.
cp TMVA_HWW.root Trainings/$dir/.