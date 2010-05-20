#!/usr/bin/python

#-----------------------------------------------------------------------------------
#usage: python getGlobalTag.py DATASETNAME
#This file returns the global tag corresponding to given dataset
#-----------------------------------------------------------------------------------

import readline
import os
import commands
import re
import optparse
import sys

parser = optparse.OptionParser ("Usage: %prog [options] dataset")
(options, args) = parser.parse_args()
neededArgs = 1

if len(args) < neededArgs:
    parser.error("Must provide 1 argument")
    sys.exit(0)
                
dataset         = args[0]

print 'Find global tag for dataset  :  ' + dataset

#get an LFN corresponding to given dataset

slashRE = re.compile(r'^/')
first = True

file = ''

for line in commands.getoutput('dbsql find file where dataset like %s' % dataset).split('\n'):
    if slashRE.search(line):
        line = line.strip()
        if first:
            first = False
            file = line

print 'Use LFN                      :  ' + file

#get local file corresponding to LFN
rfioRE = re.compile(r'^rfio:')
rfio = ''

for line in commands.getoutput('edmFileUtil %s' % file).split('\n'):
    if rfioRE.search(line):
        line = line.strip()
        rfio = line.split(' ')[0]

print 'Using file                   :  ' + rfio

#get global tag from edmProvDump output
globaltagRE1 = re.compile(r'GlobalTag RECO')
globaltagRE2 = re.compile(r'globaltag')
foundGlobalTag = False
globaltag = ''

for line in commands.getoutput('edmProvDump %s' % rfio).split('\n'):

    if globaltagRE1.search(line):
        line = line.strip()
        foundGlobalTag = True

    if foundGlobalTag:
        if globaltagRE2.search(line):
            line = line.strip()
            globaltag = line.split("'")[1]

    
print 'GlobalTag                    :  ' + globaltag
print globaltag
