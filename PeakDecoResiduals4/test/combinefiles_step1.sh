#!/bin/csh -f

export mypath=$1
echo "Combining files at" $mypath


hadd -f ${mypath}/res/merged_1.root ${mypath}/res/temp_1*.root > ${mypath}/res/log_1.txt 2>&1 &
hadd -f ${mypath}/res/merged_2.root ${mypath}/res/temp_2*.root > ${mypath}/res/log_2.txt 2>&1 &
hadd -f ${mypath}/res/merged_3.root ${mypath}/res/temp_3*.root > ${mypath}/res/log_3.txt 2>&1 &
hadd -f ${mypath}/res/merged_4.root ${mypath}/res/temp_4*.root > ${mypath}/res/log_4.txt 2>&1 &
hadd -f ${mypath}/res/merged_5.root ${mypath}/res/temp_5*.root > ${mypath}/res/log_5.txt 2>&1 &
hadd -f ${mypath}/res/merged_6.root ${mypath}/res/temp_6*.root > ${mypath}/res/log_6.txt 2>&1 &
hadd -f ${mypath}/res/merged_7.root ${mypath}/res/temp_7*.root > ${mypath}/res/log_7.txt 2>&1 &
hadd -f ${mypath}/res/merged_8.root ${mypath}/res/temp_8*.root > ${mypath}/res/log_8.txt 2>&1 &
hadd -f ${mypath}/res/merged_9.root ${mypath}/res/temp_9*.root > ${mypath}/res/log_9.txt 2>&1 &

