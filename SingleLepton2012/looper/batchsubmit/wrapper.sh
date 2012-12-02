#!/bin/bash

#
# args
#

FILEID=$1
FILE=$2
COPYDIR=$3

echo "[wrapper] FILEID    = " ${FILEID}
echo "[wrapper] FILE      = " ${FILE}
echo "[wrapper] COPYDIR   = " ${COPYDIR}

#
# set up environment
#

echo "[wrapper] setting env"
export CMS_PATH=/code/osgcode/cmssoft/cms
export SCRAM_ARCH=slc5_amd64_gcc434
source /code/osgcode/cmssoft/cmsset_default.sh
export ROOTSYS=/code/osgcode/cmssoft/cms/slc5_amd64_gcc434/lcg/root/5.27.06b-cms21/
export LD_LIBRARY_PATH=$ROOTSYS/lib:/code/osgcode/cmssoft/cms/slc5_amd64_gcc434/external/gcc/4.3.4-cms/lib64:$LD_LIBRARY_PATH
export PATH=$HOME/bin:$ROOTSYS/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

#echo "[wrapper] printing env"
#printenv

#
# untar input sandbox
#

echo "[wrapper] extracting input sandbox"
tar -zxf input.tar.gz

#source job_input/setupenv.sh
#printenv

cd job_input
echo "[wrapper] input contents are"
ls -a
mkdir output

echo "[wrapper] directory contents are"
ls

#
# run it
#
echo "
{
gROOT->ProcessLine(\".L processBaby.C\");
processBaby(\"${FILEID}\", \"${FILE}\");
}
" > runme.C

echo "[wrapper] running"
cat runme.C

root -b -q runme.C

#
# do something with output
#

echo "[wrapper] output is"
ls
ls output

#
# clean up
#

echo "[wrapper] copying file"
OUTPUT=`ls output/ | grep ${FILEID}`
echo "[wrapper] OUTPUT = " ${OUTPUT}

if [ ! -d "${COPYDIR}" ]; then
    echo "creating output directory " ${COPYDIR}
    mkdir ${COPYDIR}
fi

lcg-cp -b -D srmv2 --vo cms -t 2400 --verbose file:`pwd`/output/${OUTPUT} srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=${COPYDIR}/${OUTPUT}

echo "[wrapper] cleaning up"
for FILE in `find . -not -name "*stderr" -not -name "*stdout"`; do rm -rf $FILE; done
echo "[wrapper] cleaned up"
ls


