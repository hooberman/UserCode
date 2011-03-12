#!/bin/bash

if [[ $1 == "" ]];  then
    echo "usage: ./makeLooper.sh RELEASE"
    exit
fi

export dir=$1

echo "Making directory" $dir
cvs co -d $dir UserCode/benhoob/tcmet_looper/simpleLooper

cd $dir

echo "Checking out CORE"
cvs co -d CORE  UserCode/JRibnik/CMS2/NtupleMacros/CORE
mv CMS2.cc CORE/.
mv CMS2.h CORE/.

echo "Checking out Tools"
cvs co -d Tools UserCode/JRibnik/CMS2/NtupleMacros/Tools

mkdir output

cd Tools/MiniFWLite
make clean
make
cd -
